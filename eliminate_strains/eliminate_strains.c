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
#include "bowtie_alignment.h"
#include "file_utils.h"
#include "align_references.h"
#include "calculate_allele_freq.h"
#include "build_mismatch_matrix.h"
#include "calculate_proportions.h"

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
	struct timespec program_start = {0, 0}, program_end = {0, 0};
	clock_gettime(CLOCK_MONOTONIC, &program_start);

	struct timespec tstart = {0, 0}, tend = {0, 0};
	
	Options opt;
	opt.paired = 0;
	opt.em_error = 0.005;
	opt.fasta_format = 0;
	opt.llr = 0;
	opt.num_top_strains_llr = 10;
	opt.num_threads = 1;
	opt.max_num_reads = 100000;
	opt.no_read_bam = 0;
	opt.num_references = 1;
	opt.verbose = 0;
	opt.num_threads = 1;
	strcpy(opt.working_dir, ".");
	strcpy(opt.output_dir, ".");
	parse_options(argc, argv, &opt);

	int i, ref_idx;

	if (opt.em_error < 0.001 || opt.em_error > 1)
    {
        fprintf(stderr, "Error: EM error rate should be in the range (0.001, 1)!\n");
        exit(1);
    }

	printf("Reading in MSA...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	MSA msa_str = read_in_msa(opt.msa_filepath);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	printf("Number of strains in MSA: %d\n", msa_str.num_sequences);
	printf("MSA sequence length: %d\n\n", msa_str.sequence_length);

	printf("Aligning reference strains...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	ReferencesData references_data_str = align_references(reference_sequences_filepaths, opt.non_imputed_positions_filepath, &msa_str, opt.num_references);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));

	printf("Performing Bowtie2 alignments...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	perform_bowtie_alignments(opt.reference_sequences_dir, opt.bowtie2_indexes_dir, opt.reference_strain_names, opt.num_references, opt.single_end_filepath, opt.forward_end_filepath, opt.reverse_end_filepath, opt.working_dir, opt.paired, opt.fasta_format, opt.verbose)
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));

	printf("Reading in SAM results...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	references_data_str.sam_results_str = read_in_sam_results(opt.working_dir, opt.reference_strain_names, opt.num_references);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	printf("Number of lines in SAM files: %d\n", references_data_str.sam_results_str.num_sam_lines);
	printf("SAM max line length: %d\n\n", references_data_str.sam_results_str.max_sam_line_length);

	// NOTE: opt.no_read_bam has no effect currently
	printf("Building mismatch matrix...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	MismatchData mismatch_data_str = build_mismatch_matrix(&references_data_str, &msa_str, opt.paired, opt.num_threads);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));

	for (ref_idx = 0; ref_idx < opt.num_references; ref_idx++)
	{
		for (i = 0; i < references_data_str.sam_results_str.num_sam_lines; i++)
		{
			free(references_data_str.sam_results_str.sam_results[ref_idx][i]);
		}
		free(references_data_str.sam_results_str.sam_results[ref_idx]);
	}
	free(references_data_str.sam_results_str.sam_results);

	for (ref_idx = 0; ref_idx < opt.num_references; ref_idx++)
	{
		free(references_data_str.references_to_msa_positions[ref_idx]);
		free(references_data_str.reference_names[ref_idx]);
	}
	free(references_data_str.references_to_msa_positions);
	free(references_data_str.reference_names);
	free(references_data_str.reference_sequence_msa_indexes);
	free(references_data_str.reference_sequence_lengths);

	for (i = 0; i < msa_str.num_sequences; i++)
	{
		free(msa_str.sequences[i]);
		free(msa_str.sequence_names[i]);
	}
	free(msa_str.sequences);
	free(msa_str.sequence_names);

	for (ref_idx = 0; ref_idx < opt.num_references; ref_idx++)
	{
		free(reference_sequences_filepaths[ref_idx]);
	}
	free(reference_sequences_filepaths);

	printf("Calculating proportions of each strain...\n");
	calculate_proportions(&mismatch_data_str, opt.output_dir, opt.em_error, opt.llr, opt.num_top_strains_llr, 1, opt.num_threads);
	
	for (i = 0; i < mismatch_data_str.num_msa_sequences; i++)
	{
		free(mismatch_data_str.msa_sequence_names[i]);
	}
	free(mismatch_data_str.msa_sequence_names);
	for (i = 0; i < mismatch_data_str.num_reads; i++)
	{
		free(mismatch_data_str.read_names[i]);
	}
	free(mismatch_data_str.read_names);
	free(mismatch_data_str.block_sizes);

	clock_gettime(CLOCK_MONOTONIC, &program_end);
	printf("Entire program took %.5fsec\n", ((double)program_end.tv_sec + 1.0e-9 * program_end.tv_nsec) - ((double)program_start.tv_sec + 1.0e-9 * program_start.tv_nsec));

	return 0;
}