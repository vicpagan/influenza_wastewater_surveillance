#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <zlib.h>
#include <pthread.h>

#include "global.h"
#include "options.h"
#include "msa.h"
#include "sam.h"
#include "clean_reads.h"
#include "bowtie_alignment.h"
#include "file_utils.h"
#include "align_reference.h"
#include "calculate_allele_freq.h"
#include "write_mismatch_matrix.h"

// TODO: reimplement with a new variant-sites-like file format for problematic sites
// ProblematicSites *read_in_problematic_sites(char **problematic_sites_filepaths, int num_refs){
// 	FILE* file;
// 	if (( file = fopen("problematic_sites_sarsCov2.vcf","r")) == (FILE *) NULL ) fprintf(stderr, "Problematic Sites File could not be opened.\n");
// 	char buffer[1000];
// 	char name[30];
// 	int position;
// 	char ch1[1];
// 	char ch2[1];
// 	char ch3[1];
// 	char ch4[1];
// 	char s1[10];
// 	char s2[30];
// 	int i=0;
// 	while( fgets(buffer,1000,file) != NULL){
// 		if ( buffer[0] != '#' ){
// 			sscanf(buffer,"%s\t%d\t%c\t%c\t%c\t%c\t%s\t%s",&name,&position,&ch1,&ch2,&ch3,&ch4,&s1,&s2);
// 			problematic_sites[i]=position;
// 			i++;
// 		}
// 	}
// 	fclose(file);
// 	return i;
// }


int main(int argc, char **argv)
{
	struct timespec tstart = {0, 0}, tend = {0, 0};
	Options opt;
	opt.remove_identical_sequences = 0;
	opt.paired = 0;
	opt.em_error = 0.005;
	opt.coverage = 50;
	opt.clean_reads = 0;
	opt.fasta_format = 0;
	opt.freq = 0.01;
	opt.llr = 0;
	opt.min_strains = 500;
	opt.max_strains = 10000;
	opt.num_threads = 1;
	opt.no_read_bam = 0;
	opt.deletion_threshold = 0.002;
	opt.num_references = 1;
	opt.sequence_length_threshold = 95;
	opt.trim_length = 15;
	opt.fastq_trimmer_threshold = 35;
	opt.end_region_length = 10;
	opt.end_region_error_mult = 3.0;
	opt.num_threads = 1;
	strcpy(opt.working_dir, ".");
	memset(opt.print_counts_filepath, '\0', 1000);
	memset(opt.print_deletions_filepath, '\0', 1000);
	parse_options(argc, argv, &opt);

	int i;

	if (opt.num_references <= 0)
	{
		fprintf(stderr, "Error: -N/--num-references must be set to a positive number.\n");
		exit(1);
	}

	char **msa_reference_filepaths = list_sorted_dir_files(opt.msa_reference_dir, opt.num_references, "MSA reference");
	char **bowtie2_reference_filepaths = list_sorted_dir_files(opt.bowtie2_reference_dir, opt.num_references, "Bowtie2 reference");

	if (opt.clean_reads == 1)
	{
		printf("You've selected -d to clean your FASTA/FASTQ reads. If this is not correct, please quit the program and remove the -d option. Cleaning reads...\n");
		clean_reads(opt.single_end_filepath, opt.forward_end_filepath, opt.reverse_end_filepath, opt.working_dir, opt.paired, opt.fasta_format, opt.sequence_length_threshold, opt.trim_length, opt.fastq_trimmer_threshold);
	}

	ReferenceData *reference_data_strs = (ReferenceData *)malloc(opt.num_references * sizeof(ReferenceData));

	align_all_references(reference_data_strs, opt.num_references, msa_reference_filepaths, bowtie2_reference_filepaths);

	int ref_idx;
	for (ref_idx = 0; ref_idx < opt.num_references; ref_idx++)
	{
		char sam_filename[1000];
		sprintf(sam_filename, "%s.%d", opt.sam_filepath, ref_idx);

		char *sam_path = get_filepath_in_working_dir(sam_filename, opt.working_dir);

		perform_bowtie_alignment_xeq(bowtie2_reference_filepaths[ref_idx], opt.single_end_filepath, opt.forward_end_filepath, opt.reverse_end_filepath, sam_path, opt.working_dir, opt.paired, opt.fasta_format);
		int invoke_cleaning = calculate_error_rates(sam_path, opt.end_region_length, opt.end_region_error_mult);
		if (invoke_cleaning == 1 && opt.clean_reads == 0)
		{
			// FIXME: As of now, this changes where the opt.xxx_filepath points to once clean_reads() is called. 
			// Since we now have multiple reference strains, we should determine how to handle this (either all are cleaned, or only the ones that need cleaning are cleaned)
			printf("Error rates of reads are too high! Cleaning reads...\n");
			clean_reads(opt.single_end_filepath, opt.forward_end_filepath, opt.reverse_end_filepath, opt.working_dir, opt.paired, opt.fasta_format, opt.sequence_length_threshold, opt.trim_length, opt.fastq_trimmer_threshold);
		}
		perform_bowtie_alignment(bowtie2_reference_filepaths[ref_idx], opt.single_end_filepath, opt.forward_end_filepath, opt.reverse_end_filepath, sam_path, opt.working_dir, opt.paired, opt.fasta_format);

		printf("Reading in SAM results for reference %d...\n", ref_idx);
		reference_data_strs[ref_idx].sam_results_str = read_in_sam_results(sam_path);
	}

	printf("Reading in MSA...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	char *reference_name = read_fasta_header_name(msa_reference_filepaths[0]);
	MSA msa_str = read_in_msa(opt.msa_filepath, reference_name);
	free(reference_name);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	printf("Number of strains: %d\n", msa_str.num_sequences);
	printf("MSA length: %d\n", msa_str.sequence_length);

	if (opt.remove_identical_sequences == 1)
	{
		printf("Finding identical sequences...\n");
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		remove_identical_sequences(&msa_str);
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	}

	double **allele_frequency = (double **)malloc(msa_str.sequence_length * sizeof(double *));
	for (i = 0; i < msa_str.sequence_length; i++)
	{
		allele_frequency[i] = (double *)calloc(4, sizeof(double));
	}

	int output_allele_counts = (opt.print_counts_filepath[0] != '\0');
	int output_deletions = (opt.print_deletions_filepath[0] != '\0');

	printf("Calculating allele frequencies...\n");
	if (opt.paired == 1)
	{
		calculate_allele_freq_paired(allele_frequency, &msa_str, opt.freq, tstart, tend, opt.coverage, opt.min_strains, opt.max_strains, output_allele_counts, opt.print_counts_filepath, output_deletions, opt.print_deletions_filepath, opt.deletion_threshold, reference_data_strs, opt.num_references);
	}
	else
	{
		calculate_allele_freq_single(allele_frequency, &msa_str, opt.freq, tstart, tend, opt.coverage, opt.min_strains, opt.max_strains, output_allele_counts, opt.print_counts_filepath, output_deletions, opt.print_deletions_filepath, opt.deletion_threshold, reference_data_strs, opt.num_references);
	}

	for (i = 0; i < msa_str.sequence_length; i++)
	{
		free(allele_frequency[i]);
	}
	free(allele_frequency);

	if (msa_str.num_sequences == 0)
	{
		printf("No strains remaining, exiting...\n");
		exit(1);
	}

	// NOTE: opt.no_read_bam has no effect currently
	printf("Creating mismatch matrix...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	write_mismatch_matrix(opt.output_filepath, reference_data_strs, opt.num_references, &msa_str, opt.paired, opt.num_threads);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));

	for (ref_idx = 0; ref_idx < opt.num_references; ref_idx++)
	{
		for (i = 0; i < reference_data_strs[ref_idx].sam_results_str.num_sam_lines; i++)
		{
			free(reference_data_strs[ref_idx].sam_results_str.sam_results[i]);
		}
		free(reference_data_strs[ref_idx].sam_results_str.sam_results);
		free(reference_data_strs[ref_idx].reference_index);
	}
	free(reference_data_strs);

	for (i = 0; i < msa_str.num_sequences; i++)
	{
		free(msa_str.sequences[i]);
		free(msa_str.sequence_names[i]);
	}
	free(msa_str.sequences);
	free(msa_str.sequence_names);

	for (i = 0; i < opt.num_references; i++)
	{
		free(msa_reference_filepaths[i]);
		free(bowtie2_reference_filepaths[i]);
	}
	free(msa_reference_filepaths);
	free(bowtie2_reference_filepaths);

	char *buffer = (char *)calloc(FASTA_MAXLINE, sizeof(char));
	if (opt.llr == 1)
	{
		sprintf(buffer, "Rscript ../EM_C_LLR.R -i %s -f %lf -e %lf -l -r %s -b %s", opt.output_filepath, opt.freq, opt.em_error, opt.msa_filepath, opt.print_counts_filepath);
	}
	else
	{
		sprintf(buffer, "Rscript ../EM_C_LLR.R -i %s -f %lf -e %lf -r %s -b %s", opt.output_filepath, opt.freq, opt.em_error, opt.msa_filepath, opt.print_counts_filepath);
	}
	system(buffer);
	free(buffer);

	return 0;
}