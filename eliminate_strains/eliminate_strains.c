#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <assert.h>
#include <unistd.h>
#include <pthread.h>
#include <dirent.h>

#include "global.h"
#include "options.h"
#include "align_reference.h"
#include "calculate_allele_freq.h"
#include "write_mismatch_matrix.h"

/**
 * @brief Classify a SAM FLAG value's read-pair role by inspecting specific bits.
 * @param n The SAM FLAG field value.
 * @return -1 if unmapped (bit 2 set); 2 if "not paired" (bit 3 set); otherwise
 * bit 6 of the flag: 1 = first read in pair, 0 = second read in pair.
 * @note Despite the name, this does not convert to a binary string. It extracts
 * specific SAM flag bits to answer "which mate is this?".
 */
int dec2bin(int n)
{
	int binaryNum[32];
	int i = 0;
	while (n > 0)
	{
		binaryNum[i] = n % 2;
		n = n / 2;
		i++;
	}
	if (binaryNum[2] == 1)
	{
		return -1;
	}
	else if (binaryNum[3] == 1)
	{
		// unpaired
		return 2;
	}
	else
	{
		// if 1 this is first pair, if 0 this is second pair
		return binaryNum[6];
	}
}

/**
 * @brief Determine the MSA's alignment width (number of columns) from the first entry.
 * @param MSA_file Open gzipped MSA FASTA file, positioned at the start.
 * @return Number of bases/columns in one aligned sequence.
 */
int setMSALength(gzFile MSA_file)
{
	char buffer[FASTA_MAXLINE];
	int length = 0;
	int i = 0;
	int iter = 0;
	while (gzgets(MSA_file, buffer, FASTA_MAXLINE) != NULL)
	{
		if (buffer[0] != '>')
		{
			for (i = 0; buffer[i] != '\n'; i++)
			{
				length++;
			}
		}
		else if (iter == 0)
		{
			iter++;
		}
		else
		{
			break;
		}
	}
	return length;
}

/**
 * @brief Count strains in the MSA and find the longest strain name.
 * @param MSA_file Open gzipped MSA FASTA file, positioned at the start.
 * @param strain_info Output: [0] = number of strains, [1] = max name length.
 */
void setNumStrains(gzFile MSA_file, int *strain_info)
{
	char buffer[FASTA_MAXLINE];
	int i = 0;
	int num_strains = 0;
	int max_name = 0;
	while (gzgets(MSA_file, buffer, FASTA_MAXLINE) != NULL)
	{
		if (buffer[0] == '>')
		{
			int length = strlen(buffer) - 1;
			if (length > max_name)
			{
				max_name = length;
			}
			num_strains++;
		}
	}
	strain_info[0] = num_strains;
	strain_info[1] = max_name;
}

/**
 * @brief Load the MSA into memory and locate the reference strain's row.
 *
 * Reads every strain's name and sequence into `names`/`MSA` (uppercasing
 * A/C/G/T, keeping '-' as-is, and zeroing out anything else). Separately
 * checks each name against a hardcoded list of known reference-strain names
 * to find which row is "the" reference (see @note).
 *
 * @param MSA_file Open gzipped MSA FASTA file, positioned at the start.
 * @param MSA Output: MSA[i] filled with strain i's sequence.
 * @param names Output: names[i] filled with strain i's name.
 * @param length_of_MSA Alignment width (from setMSALength()); unused here but kept for signature symmetry.
 * @return Row index of the strain matching one of the hardcoded reference names, or 0 if none matched.
 * @note The hardcoded name list is SARS-CoV-2-era; needs updating (or a config-driven
 * replacement) for influenza reference strains.
 */
int readInMSA(gzFile MSA_file, char **MSA, char **names, int length_of_MSA, char *reference_name)
{
	char buffer[FASTA_MAXLINE];
	int i = 0;
	int index = -1;
	int found_ref = 0;
	int ref_index = 0;
	while (gzgets(MSA_file, buffer, FASTA_MAXLINE) != NULL)
	{
		if (buffer[0] == '>')
		{
			index++;
			for (i = 1; buffer[i] != '\n'; i++)
			{
				names[index][i - 1] = buffer[i];
			}
			names[index][i - 1] = '\0';
			if (strcmp(names[index], reference_name) == 0)
			{
				found_ref = 1;
			}
		}
		else
		{
			int size = strlen(buffer);
			for (i = 0; i < size; i++)
			{
				if (buffer[i] == 'A' || buffer[i] == 'a')
				{
					MSA[index][i] = 'A';
					// allele_frequency[i][0]++;
				}
				else if (buffer[i] == 'G' || buffer[i] == 'g')
				{
					MSA[index][i] = 'G';
					// allele_frequency[i][1]++;
				}
				else if (buffer[i] == 'C' || buffer[i] == 'c')
				{
					MSA[index][i] = 'C';
					// allele_frequency[i][2]++;
				}
				else if (buffer[i] == 'T' || buffer[i] == 't')
				{
					MSA[index][i] = 'T';
					// allele_frequency[i][3]++;
				}
				else if (buffer[i] == '-')
				{
					MSA[index][i] = '-';
				}
				else
				{
					MSA[index][i] = '\0';
				}
			}
			if (found_ref == 1)
			{
				ref_index = index;
				/*int placement = 0;
				for(i=0; i<length_of_MSA; i++){
					if ( MSA[index][i] != '-' ){
						reference[placement] = i;
						placement++;
					}
				}*/
				found_ref = 0;
			}
		}
	}
	return ref_index;
}

/**
 * @brief Find strains with identical sequences and blank out the duplicates' names.
 * 
 * @param number_of_strains Number of rows in MSA.
 * @param length_of_MSA Number of columns in MSA.
 * @param MSA Numerically-coded alignment.
 * @param identical Scratch/output array of duplicate row indices (-1-terminated on input).
 * @param names_of_strains Strain names; duplicates get their name blanked ('\0').
 * @param maxname Length of each name buffer in names_of_strains.
 * 
 * @return Number of strains remaining after removing duplicates.
 * 
 * @note Not currently invoked -- no CLI flag sets opt.remove_identical, so this path is never reached from main().
 */
int removeIdenticalStrains(int number_of_strains, int length_of_MSA, int **MSA, int *identical, char **names_of_strains, int maxname)
{
	int i, j, k;
	int mismatch = 0;
	int index = 0;
	int *identical2 = (int *)malloc(number_of_strains * sizeof(int));
	for (i = 0; i < number_of_strains; i++)
	{
		identical2[i] = 0;
	}
	for (i = 0; i < number_of_strains; i++)
	{
		// printf("working on %d\n",i);
		for (j = i + 1; j < number_of_strains; j++)
		{
			mismatch = 0;
			for (k = 0; k < length_of_MSA; k++)
			{
				if (MSA[i][k] != MSA[j][k])
				{
					mismatch++;
				}
			}
			if (mismatch == 0)
			{
				printf("Identical seq found!\n");
				int placement = 0;
				for (k = 0; k < number_of_strains; k++)
				{
					if (identical[k] == -1)
					{
						placement = k;
						break;
					}
				}
				int found = 0;
				for (k = 0; k < number_of_strains; k++)
				{
					if (identical[k] == j)
					{
						found = 1;
					}
				}
				if (found == 0)
				{
					identical[placement] = j;
				}
			}
		}
	}
	int number_of_identical_strains = 0;
	for (i = 0; i < number_of_strains; i++)
	{
		if (identical[i] == -1)
		{
			break;
		}
		number_of_identical_strains++;
		// printf("%d\n",identical[i]);
	}
	printf("Number of identical strains: %d\n", number_of_identical_strains);
	for (i = 0; i < number_of_strains; i++)
	{
		int ident = 0;
		for (j = 0; j < number_of_identical_strains; j++)
		{
			if (i == identical[j])
			{
				ident = 1;
			}
		}
		if (ident == 0)
		{
			identical2[i] = i;
		}
		else
		{
			identical2[i] = -1;
		}
	}
	for (i = 0; i < number_of_strains; i++)
	{
		identical[i] = -1;
	}
	int placement = 0;
	for (i = 0; i < number_of_strains; i++)
	{
		if (identical2[i] != -1)
		{
			identical[placement] = identical2[i];
			placement++;
		}
		else
		{
			memset(names_of_strains[i], '\0', maxname);
		}
	}
	free(identical2);
	return number_of_strains - number_of_identical_strains;
}

/**
 * @brief Read a variant-sites file: first line is a count, remaining lines are site indices.
 * @param variant_sites_file Open file.
 * @param number_of_sites Output: [0] = count read from the first line.
 * @return Newly allocated array of site indices (caller must free()).
 */
int *readVariantSitesFile(FILE *variant_sites_file, int *number_of_sites)
{
	char buffer[20];
	int iter = 0;
	int *variant_sites;
	int placement = 0;
	while (fgets(buffer, 20, variant_sites_file) != NULL)
	{
		if (iter == 0)
		{
			number_of_sites[0] = atoi(buffer);
			variant_sites = (int *)malloc(number_of_sites[0] * sizeof(int));
			memset(variant_sites, 0, number_of_sites[0]);
			iter++;
		}
		else
		{
			variant_sites[placement] = atoi(buffer);
			placement++;
		}
	}
	return variant_sites;
}

// TODO: reimplement with a new variant-sites-like file format for problematic sites
/*int process_problematic_sites(int* problematic_sites){
	FILE* file;
	if (( file = fopen("problematic_sites_sarsCov2.vcf","r")) == (FILE *) NULL ) fprintf(stderr, "Problematic Sites File could not be opened.\n");
	char buffer[1000];
	char name[30];
	int position;
	char ch1[1];
	char ch2[1];
	char ch3[1];
	char ch4[1];
	char s1[10];
	char s2[30];
	int i=0;
	while( fgets(buffer,1000,file) != NULL){
		if ( buffer[0] != '#' ){
			sscanf(buffer,"%s\t%d\t%c\t%c\t%c\t%c\t%s\t%s",&name,&position,&ch1,&ch2,&ch3,&ch4,&s1,&s2);
			problematic_sites[i]=position;
			i++;
		}
	}
	fclose(file);
	return i;
}*/

/**
 * @brief Copy the strains that survived elimination into the pre-sized
 * global resize_MSA/resize_names_of_strains buffers, then free the originals.
 * @param number_of_strains_remaining Number of strains to copy (size of strains_kept).
 * @param strains_kept Row indices into MSA/names_of_strains of surviving strains.
 * @param MSA Full (pre-elimination) strain panel; freed at the end.
 * @param names_of_strains Full (pre-elimination) strain names; freed at the end.
 * @param maxname Length of each name buffer.
 * @param number_of_total_strains Total rows in MSA/names_of_strains before elimination.
 */
void reallocate_memory(int number_of_strains_remaining, int *strains_kept, char **MSA, char **names_of_strains, int maxname, int number_of_total_strains)
{
	int i, j;
	int counter = 0;
	for (i = 0; i < number_of_strains_remaining; i++)
	{
		strcpy(resize_names_of_strains[counter], names_of_strains[strains_kept[i]]);
		strcpy(resize_MSA[counter], MSA[strains_kept[i]]);
		counter++;
	}
	for (i = 0; i < number_of_total_strains; i++)
	{
		free(names_of_strains[i]);
		free(MSA[i]);
	}
	free(names_of_strains);
	free(MSA);
}

/**
 * @brief Load every non-header line of a SAM file into the global `sam_results` array.
 * @param sam_file Open SAM file.
 * @note Used for threaded paired-mode (opt.no_read_bam == 0); sam_results must
 * already be allocated to hold at least as many lines as are in the file.
 */
void readInSamFile(FILE *sam_file)
{
	char buffer[FASTA_MAXLINE];
	int i, j, k;
	i = 0;
	while (fgets(buffer, FASTA_MAXLINE, sam_file) != NULL)
	{
		if (buffer[0] != '@')
		{
			strcpy(sam_results[i], buffer);
			i++;
		}
	}
}
/**
 * @brief Nudge a thread's [start, end) line range so it doesn't split a read pair.
 *
 * If `start` lands on a second-in-pair SAM record, backs up one line so the
 * pair starts together; if `end` lands on a first-in-pair record, extends by
 * one line so the pair finishes together.
 *
 * @param start Proposed first line index (into sam_results) for this thread.
 * @param end Proposed last line index (into sam_results) for this thread.
 * @param return_arr Output: [0] = adjusted start, [1] = adjusted end.
 * @param max_sam_length Length of the scratch line buffer to allocate.
 */
void adjust_start_end(int start, int end, int *return_arr, int max_sam_length)
{
	char *s;
	char *buffer_copy = (char *)malloc((max_sam_length + 1) * sizeof(char));
	strcpy(buffer_copy, sam_results[start]);
	s = strtok(buffer_copy, "\t");
	s = strtok(NULL, "\t");
	int decimal = 0;
	sscanf(s, "%d", &decimal);
	decimal = dec2bin(decimal);
	if (decimal == 0)
	{ // second in pair
		return_arr[0] = start - 1;
	}
	else
	{
		return_arr[0] = start;
	}
	strcpy(buffer_copy, sam_results[end]);
	s = strtok(buffer_copy, "\t");
	s = strtok(NULL, "\t");
	decimal = 0;
	sscanf(s, "%d", &decimal);
	decimal = dec2bin(decimal);
	if (decimal == 1)
	{
		return_arr[1] = end + 1;
	}
	else
	{
		return_arr[1] = end;
	}
	free(buffer_copy);
}

/**
 * @brief Trim 15bp off each end of every read in a FASTQ file.
 * @param filename Base filename; reads "<filename>_trimmed1.fastq", writes
 * "<filename>_trimmed2.fastq".
 */
void trim_ends_fastq(char filename[])
{
	FILE *filename_ptr;
	FILE *outfile_ptr;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	char outfilename[1000];
	char infile[1000];
	strcpy(infile, filename);
	strcpy(outfilename, filename);
	strcat(infile, "_trimmed1.fastq");
	strcat(outfilename, "_trimmed2.fastq");
	filename_ptr = fopen(infile, "r");
	if (filename_ptr == NULL)
	{
		printf("Cannot open %s\n", infile);
		exit(1);
	}
	outfile_ptr = fopen(outfilename, "w");
	if (outfile_ptr == NULL)
	{
		printf("Cannot open %s\n", outfilename);
		exit(1);
	}
	int counter = 0;
	char *line0 = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *line1 = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *line2 = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *line3 = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	int i;
	for (i = 0; i < FASTA_MAXLINE; i++)
	{
		line0[i] = '\0';
		line1[i] = '\0';
		line2[i] = '\0';
		line3[i] = '\0';
	}
	while ((read = getline(&line, &len, filename_ptr)) != -1)
	{
		if (counter == 0)
		{
			strcpy(line0, line);
		}
		else if (counter == 1)
		{
			strcpy(line1, line);
		}
		else if (counter == 2)
		{
			strcpy(line2, line);
		}
		else if (counter == 3)
		{
			strcpy(line3, line);
		}
		else
		{
			printf("Error in trimming sequences!\n");
			exit(1);
		}
		if (strlen(line1) > 95 && counter == 3)
		{
			fprintf(outfile_ptr, "%s", line0);
			for (i = 15; i < strlen(line1) - 15; i++)
			{
				fprintf(outfile_ptr, "%c", line1[i]);
			}
			fprintf(outfile_ptr, "\n");
			fprintf(outfile_ptr, "%s", line2);
			for (i = 15; i < strlen(line3) - 15; i++)
			{
				fprintf(outfile_ptr, "%c", line3[i]);
			}
			fprintf(outfile_ptr, "\n");
		}
		counter = counter + 1;
		if (counter == 4)
		{
			counter = 0;
		}
	}
	fclose(filename_ptr);
	fclose(outfile_ptr);
	free(line0);
	free(line1);
	free(line2);
	free(line3);
}
/**
 * @brief Trim 15bp off each end of every sequence in a FASTA file.
 * @param filename Base filename; reads "<filename>.fasta", writes "<filename>_trimmed2.fasta".
 */
void trim_ends_fasta(char filename[])
{
	FILE *filename_ptr;
	FILE *outfile_ptr;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	char outfilename[1000];
	char infile[1000];
	strcpy(infile, filename);
	strcpy(outfilename, filename);
	strcat(infile, ".fasta");
	strcat(outfilename, "_trimmed2.fasta");
	filename_ptr = fopen(infile, "r");
	if (filename_ptr == NULL)
	{
		printf("Cannot open %s\n", infile);
		exit(1);
	}
	outfile_ptr = fopen(outfilename, "w");
	if (outfile_ptr == NULL)
	{
		printf("Cannot open %s\n", outfilename);
		exit(1);
	}
	int counter = 0;
	char *line0 = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *line1 = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	int i;
	for (i = 0; i < FASTA_MAXLINE; i++)
	{
		line0[i] = '\0';
		line1[i] = '\0';
	}
	while ((read = getline(&line, &len, filename_ptr)) != -1)
	{
		if (counter == 0)
		{
			strcpy(line0, line);
		}
		else if (counter == 1)
		{
			strcpy(line1, line);
		}
		else
		{
			printf("Error in trimming sequences!\n");
			exit(1);
		}
		if (strlen(line1) > 95 && counter == 1)
		{
			fprintf(outfile_ptr, "%s", line0);
			for (i = 15; i < strlen(line1) - 15; i++)
			{
				fprintf(outfile_ptr, "%c", line1[i]);
			}
			fprintf(outfile_ptr, "\n");
		}
		counter = counter + 1;
		if (counter == 2)
		{
			counter = 0;
		}
	}
	fclose(filename_ptr);
	fclose(outfile_ptr);
	free(line0);
	free(line1);
}
/**
 * @brief Quality-trim and end-trim the input reads (single-end or paired), in place.
 *
 * For FASTQ input, runs `fastq_quality_trimmer` then trim_ends_fastq(); for
 * FASTA input, just runs trim_ends_fasta(). Rewrites single_end_file (or
 * forward_end_file/reverse_end_file) to point at the trimmed output --
 * since these are caller-owned buffers passed by pointer, the caller sees
 * the update (unlike the previous Options-by-value version, where this
 * rewrite was silently lost on return).
 * 
 * @param paired 1 if using paired-end reads.
 * @param fasta_format 1 if reads are FASTA (vs. FASTQ).
 * @param single_end_file Single-end reads path; rewritten in place if paired == 0.
 * @param forward_end_file Forward reads path; rewritten in place if paired == 1.
 * @param reverse_end_file Reverse reads path; rewritten in place if paired == 1.
 */
void clean_reads(int paired, int fasta_format, char *single_end_file, char *forward_end_file, char *reverse_end_file)
{
	int i;
	char *buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(buffer, '\0', FASTA_MAXLINE);
	if (paired == 0)
	{
		char *prefix;
		if (single_end_file[strlen(single_end_file) - 6] == '.' && single_end_file[strlen(single_end_file) - 5] == 'f' && single_end_file[strlen(single_end_file) - 4] == 'a' && single_end_file[strlen(single_end_file) - 3] == 's' && single_end_file[strlen(single_end_file) - 2] == 't' && (single_end_file[strlen(single_end_file) - 1] == 'q' || single_end_file[strlen(single_end_file) - 1] == 'a'))
		{
			prefix = (char *)malloc((strlen(single_end_file) - 6 + 15) * sizeof(char));
			for (i = 0; i < strlen(single_end_file) - 6; i++)
			{
				prefix[i] = single_end_file[i];
			}
			for (i = strlen(single_end_file) - 6; i < strlen(single_end_file) - 6 + 15; i++)
			{
				prefix[i] = '\0';
			}
		}
		else if (single_end_file[strlen(single_end_file) - 3] == '.' && single_end_file[strlen(single_end_file) - 2] == 'f' && (single_end_file[strlen(single_end_file) - 1] == 'q' || single_end_file[strlen(single_end_file) - 1] == 'a'))
		{
			prefix = (char *)malloc((strlen(single_end_file) - 3 + 15) * sizeof(char));
			for (i = 0; i < strlen(single_end_file) - 3; i++)
			{
				prefix[i] = single_end_file[i];
			}
			for (i = strlen(single_end_file) - 3; i < strlen(single_end_file) - 3 + 15; i++)
			{
				prefix[i] = '\0';
			}
		}
		else
		{
			printf("Your reads don't end with .fastq, .fasta, .fa, or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		if (fasta_format == 0)
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t 35 -i %s -o %s_trimmed1.fastq -Q33", single_end_file, prefix);
			system(buffer);
			trim_ends_fastq(prefix);
		}
		else
		{
			trim_ends_fasta(prefix);
		}
		// sprintf(buffer, "fastx_trimmer -m 65 -t 15 -i %s_trimmed1.fastq -o %s_trimmed2.fastq",single_end_file,prefix);
		// system(buffer);
		// sprintf(buffer, "fastx_trimmer -f 15 -i %s_trimmed2.fastq -o %s_trimmed3.fastq",single_end_file,prefix);
		// system(buffer);
		char suffix[1000];
		if (fasta_format == 0)
		{
			strcpy(suffix, "_trimmed2.fastq");
		}
		else
		{
			strcpy(suffix, "_trimmed2.fasta");
		}
		strcat(prefix, suffix);
		strcpy(single_end_file, prefix);
		free(prefix);
	}
	else
	{
		char *prefix_forward;
		if (forward_end_file[strlen(forward_end_file) - 6] == '.' && forward_end_file[strlen(forward_end_file) - 5] == 'f' && forward_end_file[strlen(forward_end_file) - 4] == 'a' && forward_end_file[strlen(forward_end_file) - 3] == 's' && forward_end_file[strlen(forward_end_file) - 2] == 't' && (forward_end_file[strlen(forward_end_file) - 1] == 'q' || forward_end_file[strlen(forward_end_file) - 1] == 'a'))
		{
			prefix_forward = (char *)malloc((strlen(forward_end_file) - 6 + 15) * sizeof(char));
			for (i = 0; i < strlen(forward_end_file) - 6; i++)
			{
				prefix_forward[i] = forward_end_file[i];
			}
			for (i = strlen(forward_end_file) - 6; i < strlen(forward_end_file) - 6 + 15; i++)
			{
				prefix_forward[i] = '\0';
			}
		}
		else if (forward_end_file[strlen(forward_end_file) - 3] == '.' && forward_end_file[strlen(forward_end_file) - 2] == 'f' && (forward_end_file[strlen(forward_end_file) - 1] == 'q' || forward_end_file[strlen(forward_end_file) - 1] == 'a'))
		{
			prefix_forward = (char *)malloc((strlen(forward_end_file) - 3 + 15) * sizeof(char));
			for (i = 0; i < strlen(forward_end_file) - 3; i++)
			{
				prefix_forward[i] = forward_end_file[i];
			}
			for (i = strlen(forward_end_file) - 3; i < strlen(forward_end_file) - 3 + 15; i++)
			{
				prefix_forward[i] = '\0';
			}
		}
		else
		{
			printf("Your reads don't end with .fastq, .fasta, .fa or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		if (fasta_format == 0)
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t 35 -i %s -o %s_trimmed1.fastq -Q33", forward_end_file, prefix_forward);
			system(buffer);
			trim_ends_fastq(prefix_forward);
		}
		else
		{
			trim_ends_fasta(prefix_forward);
		}
		char suffix[1000];
		if (fasta_format == 0)
		{
			strcpy(suffix, "_trimmed2.fastq");
		}
		else
		{
			strcpy(suffix, "_trimmed2.fasta");
		}
		strcat(prefix_forward, suffix);
		strcpy(forward_end_file, prefix_forward);
		free(prefix_forward);
		char *prefix_reverse;
		if (reverse_end_file[strlen(reverse_end_file) - 6] == '.' && reverse_end_file[strlen(reverse_end_file) - 5] == 'f' && reverse_end_file[strlen(reverse_end_file) - 4] == 'a' && reverse_end_file[strlen(reverse_end_file) - 3] == 's' && reverse_end_file[strlen(reverse_end_file) - 2] == 't' && (reverse_end_file[strlen(reverse_end_file) - 1] == 'q' || reverse_end_file[strlen(reverse_end_file) - 1] == 'a'))
		{
			prefix_reverse = (char *)malloc((strlen(reverse_end_file) - 6 + 15) * sizeof(char));
			for (i = 0; i < strlen(reverse_end_file) - 6; i++)
			{
				prefix_reverse[i] = reverse_end_file[i];
			}
			for (i = strlen(reverse_end_file) - 6; i < strlen(reverse_end_file) - 6 + 15; i++)
			{
				prefix_reverse[i] = '\0';
			}
		}
		else if (reverse_end_file[strlen(reverse_end_file) - 3] == '.' && reverse_end_file[strlen(reverse_end_file) - 2] == 'f' && (reverse_end_file[strlen(reverse_end_file) - 1] == 'q' || reverse_end_file[strlen(reverse_end_file) - 1] == 'a'))
		{
			prefix_reverse = (char *)malloc((strlen(reverse_end_file) - 3 + 15) * sizeof(char));
			for (i = 0; i < strlen(reverse_end_file) - 3; i++)
			{
				prefix_reverse[i] = reverse_end_file[i];
			}
			for (i = strlen(reverse_end_file) - 3; i < strlen(reverse_end_file) - 3 + 15; i++)
			{
				prefix_reverse[i] = '\0';
			}
		}
		else
		{
			printf("Your reads don't end with .fastq, .fasta, .fa, or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		if (fasta_format == 0)
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t 35 -i %s -o %s_trimmed1.fastq -Q33", reverse_end_file, prefix_reverse);
			system(buffer);
			trim_ends_fastq(prefix_reverse);
		}
		else
		{
			trim_ends_fasta(prefix_reverse);
		}
		strcat(prefix_reverse, suffix);
		strcpy(reverse_end_file, prefix_reverse);
		free(prefix_reverse);
	}
}

/**
 * @brief Build the Bowtie2 index (if missing) and run Bowtie2 to produce sam_path.
 * @param paired 1 if using paired-end reads.
 * @param fasta_format 1 if reads are FASTA (vs. FASTQ).
 * @param bowtie_reference_path Path to the reference genome to build/align against.
 * @param forward_end_file Forward reads path (paired mode).
 * @param reverse_end_file Reverse reads path (paired mode).
 * @param single_end_file Single-end reads path.
 * @param sam_path Output SAM file path.
 */
void perform_bowtie_alignment(int paired, int fasta_format, char *bowtie_reference_path, char *forward_end_file, char *reverse_end_file, char *single_end_file, char *sam_path)
{
	char *buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(buffer, '\0', FASTA_MAXLINE);
	// FIXME: Not sure if we should keep this check but make it unhardcoded or not
	// if (access("/space/lenore/influenza/references/reference_sequences/cat_ref_19088566.fasta.1.bt2", F_OK) == 0)
	// {
	// 	printf("Wuhan reference file MN908947.3.fasta.1.bt2 exists. Not rebuilding bowtie2-build database.\n");
	// }
	// else
	// {
		sprintf(buffer, "bowtie2-build -f %s %s", bowtie_reference_path, bowtie_reference_path);
		system(buffer);
	// }
	if (paired == 1 && fasta_format == 1)
	{
		sprintf(buffer, "bowtie2 --all --no-unal -f -x %s -1 %s -2 %s -S %s", bowtie_reference_path, forward_end_file, reverse_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 0 && fasta_format == 1)
	{
		sprintf(buffer, "bowtie2 --all --no-unal -f -x %s -U %s -S %s", bowtie_reference_path, single_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 1 && fasta_format == 0)
	{
		sprintf(buffer, "bowtie2 --all --no-unal -x %s -1 %s -2 %s -S %s", bowtie_reference_path, forward_end_file, reverse_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 0 && fasta_format == 0)
	{
		sprintf(buffer, "bowtie2 --all --no-unal -x %s -U %s -S %s", bowtie_reference_path, single_end_file, sam_path);
		system(buffer);
	}
	free(buffer);
}

/**
 * @brief Same as perform_bowtie_alignment(), but with --xeq (extended CIGAR,
 * distinguishing '=' match from 'X' mismatch), used only for error-rate estimation.
 * @param paired 1 if using paired-end reads.
 * @param fasta_format 1 if reads are FASTA (vs. FASTQ).
 * @param bowtie_reference_path Path to the reference genome to build/align against.
 * @param forward_end_file Forward reads path (paired mode).
 * @param reverse_end_file Reverse reads path (paired mode).
 * @param single_end_file Single-end reads path.
 * @param sam_path Output SAM file path.
 */
void perform_bowtie_alignment_xeq(int paired, int fasta_format, char *bowtie_reference_path, char *forward_end_file, char *reverse_end_file, char *single_end_file, char *sam_path)
{
	char *buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(buffer, '\0', FASTA_MAXLINE);
	// FIXME: Not sure if we should keep this check but make it unhardcoded or not
	// if (access("/space/lenore/influenza/references/reference_sequences/cat_ref_19088566.fasta.1.bt2", F_OK) == 0)
	// {
	// 	printf("Wuhan reference file MN908947.3.fasta.1.bt2 exists. Not rebuilding bowtie2-build database.\n");
	// }
	// else
	// {
		sprintf(buffer, "bowtie2-build -f %s %s", bowtie_reference_path, bowtie_reference_path);
		system(buffer);
	// }
	if (paired == 1 && fasta_format == 1)
	{
		sprintf(buffer, "bowtie2 --all --xeq --no-unal -f -x %s -1 %s -2 %s -S %s", bowtie_reference_path, forward_end_file, reverse_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 0 && fasta_format == 1)
	{
		sprintf(buffer, "bowtie2 --all --xeq --no-unal -f -x %s -U %s -S %s", bowtie_reference_path, single_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 1 && fasta_format == 0)
	{
		sprintf(buffer, "bowtie2 --all --xeq --no-unal -x %s -1 %s -2 %s -S %s", bowtie_reference_path, forward_end_file, reverse_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 0 && fasta_format == 0)
	{
		sprintf(buffer, "bowtie2 --all --xeq --no-unal -x %s -U %s -S %s", bowtie_reference_path, single_end_file, sam_path);
		system(buffer);
	}
	free(buffer);
}

/**
 * @brief Estimate whether reads need quality/end trimming, from an --xeq SAM file.
 *
 * Counts mismatches ('X' CIGAR ops) in the first/last 10bp of each read
 * ("ends") versus the middle, and flags cleaning as needed if the per-base
 * error rate at the ends is more than 3x the middle's rate.
 *
 * @param sam_file Path to a SAM file produced with perform_bowtie_alignment_xeq().
 * @return 1 if reads should be cleaned (see clean_reads()), 0 otherwise.
 */
int calculate_error_rates(char sam_file[])
{
	FILE *sam_ptr;
	if ((sam_ptr = fopen(sam_file, "r")) == (FILE *)NULL)
		fprintf(stderr, "Sam file could not be opened.\n");
	int i, j;
	char buffer[FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	int mismatches_ends = 0;
	int mismatches_middle = 0;
	int total_aligned_reads = 0;
	int invoke_cleaning = 0;
	int total_length = 0;
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar[i] = 0;
		cigar_chars[i] = '\0';
	}
	while (fgets(buffer, FASTA_MAXLINE, sam_ptr) != NULL)
	{
		if (buffer[0] != '@')
		{
			total_aligned_reads++;
			char *buffer_copy = strdup(buffer);
			int length_of_sam = strlen(buffer);
			s = strtok(buffer, "\t");
			for (i = 0; i < 5; i++)
			{
				s = strtok(NULL, "\t");
			}
			char *cigar_string;
			cigar_string = strdup(s);
			char *copy = strdup(cigar_string);
			char *res = strtok(cigar_string, "=XID");
			int index = 0;
			while (res)
			{
				int from = res - cigar_string + strlen(res);
				int cigar_count = 0;
				sscanf(res, "%d", &cigar_count);
				res = strtok(NULL, "=XID");
				int to = res != NULL ? res - cigar_string : strlen(copy);
				char cigar_char = '\0';
				sscanf(copy + from, "%c", &cigar_char);
				cigar[index] = cigar_count;
				cigar_chars[index] = cigar_char;
				index++;
			}
			free(copy);
			free(cigar_string);
			s = strtok(buffer_copy, "\t");
			for (i = 0; i < 9; i++)
			{
				s = strtok(NULL, "\t");
			}
			int cigar_char_count = index;
			int position = 0;
			char *sequence = s;
			int length_of_seq = strlen(sequence);
			total_length = total_length + length_of_seq;
			for (i = 0; i < cigar_char_count; i++)
			{
				for (j = 0; j < cigar[i]; j++)
				{
					if (cigar_chars[i] == '=' || cigar_chars[i] == 'X' || cigar_chars[i] == 'I')
					{
						position++;
					}
					if (cigar_chars[i] == 'X' && position < 10)
					{
						mismatches_ends++;
					}
					else if (cigar_chars[i] == 'X' && position > length_of_seq - 11)
					{
						mismatches_ends++;
					}
					else if (cigar_chars[i] == 'X')
					{
						mismatches_middle++;
					}
				}
			}
			free(buffer_copy);
		}
	}
	fclose(sam_ptr);
	double average_size;
	average_size = total_length / total_aligned_reads;
	double error_rate_ends;
	double error_rate_middle;
	error_rate_ends = mismatches_ends / total_aligned_reads;
	error_rate_ends = error_rate_ends / 20;
	error_rate_middle = mismatches_middle / total_aligned_reads;
	error_rate_middle = error_rate_middle / (average_size - 20);
	printf("Error rate ends: %lf\n", error_rate_ends);
	printf("Error rate middle: %lf\n", error_rate_middle);
	if (error_rate_ends > 3 * error_rate_middle)
	{
		printf("The error rate in your 5' and 3' ends of your reads is 3x larger than the error rate in the middle of your reads. Your reads need cleaning... quality filtering and trimming ends\n");
		invoke_cleaning = 1;
	}
	return invoke_cleaning;
}

/** @brief qsort comparator for sorting an array of C strings alphabetically. */
int cmp_str(const void *a, const void *b)
{
	return strcmp(*(const char **)a, *(const char **)b);
}

/**
 * @brief Read the strain name from the header line of a single-sequence FASTA file.
 * @param fasta_path Path to a FASTA file whose first line is ">strain_name".
 * @return Newly allocated string with the header's name (no '>', no newline).
 */
char *read_fasta_header_name(char *fasta_path)
{
	FILE *file = fopen(fasta_path, "r");
	if (file == NULL)
	{
		fprintf(stderr, "Error: could not open reference file '%s' to read its name\n", fasta_path);
		exit(1);
	}
	char buffer[FASTA_MAXLINE];
	if (fgets(buffer, FASTA_MAXLINE, file) == NULL || buffer[0] != '>')
	{
		fprintf(stderr, "Error: '%s' does not start with a FASTA header line\n", fasta_path);
		exit(1);
	}
	fclose(file);
	buffer[strcspn(buffer, "\r\n")] = '\0';
	return strdup(buffer + 1); // skip the leading '>'
}

/**
 * @brief List the first `num_references` regular files in a directory, sorted
 * alphabetically, as full paths.
 *
 * Used to line up the Nth file across msa/reference/variant directories as
 * belonging to the same subtype (see main()). Exits the program with an
 * error if the directory has fewer than num_references usable files.
 *
 * @param dir_path Directory to scan.
 * @param num_references Number of files required (and returned).
 * @param dir_label Human-readable label for this directory, used in error messages.
 * @return Newly allocated array of num_references newly allocated path strings
 * (caller owns both the array and each string).
 */
char **list_sorted_dir_files(char *dir_path, int num_references, char *dir_label)
{
	DIR *dir = opendir(dir_path);
	if (dir == NULL)
	{
		fprintf(stderr, "Error: could not open %s directory '%s'\n", dir_label, dir_path);
		exit(1);
	}
 
	char **filenames = NULL;
	int count = 0;
	struct dirent *entry;
	while ((entry = readdir(dir)) != NULL)
	{
		if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
		{
			continue;
		}
		filenames = (char **)realloc(filenames, (count + 1) * sizeof(char *));
		filenames[count] = strdup(entry->d_name);
		count++;
	}
	closedir(dir);
 
	if (count < num_references)
	{
		fprintf(stderr, "Error: %s directory '%s' has %d file(s), but %d reference(s) were requested (-N).\n", dir_label, dir_path, count, num_references);
		exit(1);
	}
 
	// sort alphabetically so position i means the same subtype across all
	// four directories (msa/reference/bowtie2-reference/variant)
	qsort(filenames, count, sizeof(char *), cmp_str);
 
	if (count > num_references)
	{
		fprintf(stderr, "Warning: %s directory '%s' has %d files; only using the first %d (sorted alphabetically).\n", dir_label, dir_path, count, num_references);
	}
 
	char **paths = (char **)malloc(num_references * sizeof(char *));
	int i;
	for (i = 0; i < num_references; i++)
	{
		paths[i] = (char *)malloc((strlen(dir_path) + strlen(filenames[i]) + 2) * sizeof(char));
		sprintf(paths[i], "%s/%s", dir_path, filenames[i]);
	}
 
	for (i = 0; i < count; i++)
	{
		free(filenames[i]);
	}
	free(filenames);
 
	return paths;
}


/**
 * @brief 
 * 
 * @param sam_results 
 * @param num_sam_lines 
 * @param length_of_MSA 
 * @param number_of_strains_remaining 
 * @param map 
 * @param resize_MSA 
 * @param reference_index 
 */
void compute_min_mismatches_for_subtype(char **sam_results, int num_sam_lines, int length_of_MSA, int number_of_strains_remaining, char **resize_MSA, int *reference_index)
{
	int i, j, k;
	char buffer[FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar[i] = 0;
		cigar_chars[i] = '\0';
	}
	int *number_of_mismatches = (int *)malloc(number_of_strains_remaining * sizeof(int));
	int alignment_size;
	int first_in_pair = 0;
	int second_in_pair = 0;
	int first_seq_cigar[MAX_CIGAR];
	char first_seq_cigar_chars[MAX_CIGAR];
	int first_seq_cigar_char_count = 0;
	char *first_seq = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	memset(first_seq, '\0', MAX_READ_LENGTH);
	int first_seq_start_pos = 0;
	int visited[MAX_READ_LENGTH];
	int visited_place = 0;
	// memset(visited,-1,MAX_READ_LENGTH);
	for (i = 0; i < MAX_READ_LENGTH; i++)
	{
		visited[i] = -1;
	}
	int line_number = 0;
	int line_count;
	char *current_readname = NULL;
	char *context = NULL;
	for (line_count = 0; line_count < num_sam_lines; line_count++)
	{
		line_number++;
		char *buffer_copy = strdup(sam_results[line_count]);
		s = strtok_r(sam_results[line_count], "\t", &context);
		char *name = strdup(s);
		s = strtok_r(NULL, "\t", &context);
		int decimal = 0;
		sscanf(s, "%d", &decimal);
		decimal = dec2bin(decimal);
		if (decimal == 1)
		{
			if (current_readname != NULL)
			{
				free(current_readname);
			}
			current_readname = name;
			first_in_pair = 1;
		}
		else if (decimal == 0)
		{
			second_in_pair = 1;
			free(name);
		}
		else if (decimal == 2)
		{
			if (current_readname != NULL)
			{
				free(current_readname);
			}
			current_readname = name;
		}
		else
		{
			free(name);
		}
		for (i = 0; i < 2; i++)
		{
			s = strtok_r(NULL, "\t", &context);
		}
		int position = 0;
		sscanf(s, "%d", &position);
		position--;
		if (first_in_pair == 1)
		{
			first_seq_start_pos = position;
		}
		s = strtok_r(NULL, "\t", &context);
		s = strtok_r(NULL, "\t", &context);
		char *cigar_string;
		cigar_string = strdup(s);
		char *copy = strdup(cigar_string);
		char *res = strtok_r(cigar_string, "MID", &context);
		int index = 0;
		while (res)
		{
			int from = res - cigar_string + strlen(res);
			int cigar_count = 0;
			sscanf(res, "%d", &cigar_count);
			res = strtok_r(NULL, "MID", &context);
			int to = res != NULL ? res - cigar_string : strlen(copy);
			char cigar_char = '\0';
			sscanf(copy + from, "%c", &cigar_char);
			cigar[index] = cigar_count;
			if (first_in_pair == 1)
			{
				first_seq_cigar[index] = cigar_count;
			}
			cigar_chars[index] = cigar_char;
			if (first_in_pair == 1)
			{
				first_seq_cigar_chars[index] = cigar_char;
			}
			index++;
		}
		free(copy);
		free(cigar_string);
		s = strtok_r(buffer_copy, "\t", &context);
		for (i = 0; i < 9; i++)
		{
			s = strtok_r(NULL, "\t", &context);
		}
		char *sequence = s;
		if (first_in_pair == 1)
		{
			strcpy(first_seq, sequence);
			first_seq_cigar_char_count = index;
		}
		int cigar_char_count = index;
		index = 0;
		int start = 0;
		int start_ref = 0;
		if (decimal == 1 || decimal == 2)
		{
			alignment_size = 0;
		}
		if (decimal == 1 || decimal == 2)
		{
			for (i = 0; i < number_of_strains_remaining; i++)
			{
				number_of_mismatches[i] = 0;
			}
		}
		int l = 0;
		int m = 0;
		if (decimal != -1)
		{
			visited_place = 0;
			if (second_in_pair == 1)
			{
				int start_ref1 = 0;
				int start1 = 0;
				for (i = 0; i < number_of_strains_remaining; i++)
				{
					visited_place = 0;
					start1 = 0;
					start_ref1 = 0;
					for (j = 0; j < first_seq_cigar_char_count; j++)
					{
						for (k = 0; k < first_seq_cigar[j]; k++)
						{
							int start2 = 0;
							int start_ref2 = 0;
							if (start_ref1 + first_seq_start_pos + k >= position)
							{
								for (l = 0; l < cigar_char_count; l++)
								{
									int position_in_MSA = reference_index[start_ref1 + first_seq_start_pos + k];
									for (m = 0; m < cigar[l]; m++)
									{
										if (start_ref1 + first_seq_start_pos + k == start_ref2 + position + m && first_seq_cigar_chars[j] == 'M')
										{
											if (cigar_chars[l] == 'M')
											{
												if (first_seq[start1 + k] != sequence[start2 + m])
												{
													if (first_seq[start1 + k] == 'A' || first_seq[start1 + k] == 'a')
													{
														if (resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
														{
															if (position_in_MSA < length_of_MSA)
															{
																number_of_mismatches[i]--;
																assert(number_of_mismatches[i] >= 0);
															}
														}
													}
													else if (first_seq[start1 + k] == 'G' || first_seq[start1 + k] == 'g')
													{
														if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
														{
															if (position_in_MSA < length_of_MSA)
															{
																number_of_mismatches[i]--;
																assert(number_of_mismatches[i] >= 0);
															}
														}
													}
													else if (first_seq[start1 + k] == 'C' || first_seq[start1 + k] == 'c')
													{
														if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
														{
															if (position_in_MSA < length_of_MSA)
															{
																number_of_mismatches[i]--;
																assert(number_of_mismatches[i] >= 0);
															}
														}
													}
													else if (first_seq[start1 + k] == 'T' || first_seq[start1 + k] == 't')
													{
														if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
														{
															if (position_in_MSA < length_of_MSA)
															{
																number_of_mismatches[i]--;
																assert(number_of_mismatches[i] >= 0);
															}
														}
													}
												}
											}
											/*if (cigar_chars[l]== 'D'){
												if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0'){
													if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
														number_of_mismatches[i]--;
													}
												}
											}*/
											/*if (cigar_chars[l] == 'I'){
												if (first_seq[k+start1] == 'A' || first_seq[k+start1] == 'a'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'A' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}else if ( first_seq[start1+k]=='G' || first_seq[start1+k]=='g'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'G' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}else if ( first_seq[start1+k]=='C' || first_seq[start1+k]=='c'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'C' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}else if ( first_seq[start1+k]=='T' || first_seq[start1+k]=='t'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'T' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}
											}*/
											if (cigar_chars[l] == 'M')
											{
												visited[visited_place] = start_ref1 + first_seq_start_pos + k;
												visited_place++;
											}
										}
									}
									//}
									if (cigar_chars[l] == 'M')
									{
										start2 = start2 + cigar[l];
										start_ref2 = start_ref2 + cigar[l];
									}
									if (cigar_chars[l] == 'I')
									{
										start2 = start2 + cigar[l];
									}
									if (cigar_chars[l] == 'D')
									{
										start_ref2 = start_ref2 + cigar[l];
									}
								}
							}
							// if ( first_seq_cigar_chars[j] != 'I' ){
							//	start_ref1 = start_ref1+first_seq_cigar[j];
							//	start1 = start1 + first_seq_cigar[j];
							// }
						}
						if (first_seq_cigar_chars[j] == 'M')
						{
							start1 = start1 + first_seq_cigar[j];
							start_ref1 = start_ref1 + first_seq_cigar[j];
						}
						if (first_seq_cigar_chars[j] == 'I')
						{
							start1 = start1 + first_seq_cigar[j];
						}
						if (first_seq_cigar_chars[j] == 'D')
						{
							start_ref1 = first_seq_cigar[j] + start_ref1;
						}
					}
				}
			}
			for (i = 0; i < number_of_strains_remaining; i++)
			{
				start = 0;
				start_ref = 0;
				for (j = 0; j < cigar_char_count; j++)
				{
					for (k = 0; k < cigar[j]; k++)
					{
						int skip = 0;
						int position_in_MSA = reference_index[k + position + start_ref];
						for (l = 0; l < visited_place; l++)
						{
							if (position_in_MSA == visited[l])
							{
								skip = 1;
							}
						}
						if (cigar_chars[j] == 'M')
						{
							if (sequence[k + start] == 'A' || sequence[k + start] == 'a')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'G' || sequence[k + start] == 'g')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'C' || sequence[k + start] == 'c')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
								//		number_of_mismatches[i]++;
								//
								// }}
								if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'T' || sequence[k + start] == 't')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										number_of_mismatches[i]++;
									}
								}
							}
						}
						if (cigar_chars[j] == 'D')
						{
							// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
							//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
							//		number_of_mismatches[i]++;
							//	}
							// }
							if (resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
							{
								if (position_in_MSA < length_of_MSA && skip == 0)
								{
									//	number_of_mismatches[i]++;
								}
							}
						}
						if (cigar_chars[j] == 'I')
						{
							if (sequence[k + start] == 'A' || sequence[k + start] == 'a')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
							if (sequence[k + start] == 'G' || sequence[k + start] == 'g')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
							if (sequence[k + start] == 'C' || sequence[k + start] == 'c')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
							if (sequence[k + start] == 'T' || sequence[k + start] == 't')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
						}
					}
					if (cigar_chars[j] == 'M')
					{
						start = cigar[j] + start;
						start_ref = cigar[j] + start_ref;
						if (i == 0)
						{
							alignment_size = alignment_size + cigar[j];
						}
					}
					if (cigar_chars[j] == 'I')
					{
						start = cigar[j] + start;
					}
					if (cigar_chars[j] == 'D')
					{
						start_ref = cigar[j] + start_ref;
					}
					// if (i==0){
					//	alignment_size = alignment_size + cigar[j];
					// }
				}
			}
		}
		if (decimal == 0 || decimal == 2)
		{
			int min_mismatch = number_of_mismatches[0];
			for (i = 1; i < number_of_strains_remaining; i++)
			{
				if (number_of_mismatches[i] < min_mismatch)
				{
					min_mismatch = number_of_mismatches[i];
				}
			}
			update_best_match(map, current_readname, min_mismatch);
			free(current_readname);
			current_readname = NULL;
		}
		free(buffer_copy);
		first_in_pair = 0;
		second_in_pair = 0;
	}
	free(number_of_mismatches);
}

/**
 * @brief Entry point: runs the full strain-elimination + mismatch-matrix pipeline.
 *
 * Stages: parse CLI options -> (optionally) clean reads -> align reference to
 * build reference_index -> run Bowtie2 (twice: once --xeq to check error rates,
 * once normally) -> load the MSA panel -> compute allele frequencies and
 * eliminate incompatible strains -> build the read x strain mismatch matrix
 * (threaded if paired+read-sam mode) -> hand off to EM_C_LLR.R for the final
 * relative-abundance estimation.
 */
int main(int argc, char **argv)
{
	struct timespec tstart = {0, 0}, tend = {0.0};
	Options opt;
	opt.remove_identical = 0;
	opt.paired = 0;
	opt.error = 0.005;
	opt.coverage = 50;
	opt.clean_reads = 0;
	opt.fasta_format = 0;
	opt.freq = 0.01;
	opt.llr = 0;
	opt.min_strains = 500;
	opt.max_strains = 10000;
	opt.number_of_cores = 1;
	opt.no_read_bam = 0;
	opt.deletion_threshold = 0.002;
	opt.num_references = 0;
	memset(opt.print_counts, '\0', 1000);
	memset(opt.print_deletions, '\0', 1000);
	parse_options(argc, argv, &opt);
	int i, j, k, l;

	// check that the user specified a positive number of references
	if (opt.num_references <= 0)
	{
		fprintf(stderr, "Error: -N/--num-references must be set to a positive number.\n");
		exit(1);
	}
 
	char **MSA_paths = list_sorted_dir_files(opt.MSA_dir, opt.num_references, "MSA");
	char **MSA_reference_paths = list_sorted_dir_files(opt.MSA_reference_dir, opt.num_references, "MSA reference");
	char **bowtie2_reference_paths = list_sorted_dir_files(opt.bowtie2_reference_dir, opt.num_references, "Bowtie2 reference");
	char **variant_paths = list_sorted_dir_files(opt.variant_dir, opt.num_references, "variant sites");
 
	if (opt.clean_reads == 1)
	{
		printf("You've selected -d to clean your FASTA/FASTQ reads. If this is not correct, please quit the program and removed the -d option. Cleaning reads...\n");
		clean_reads(opt.paired, opt.fasta_format, opt.single_end_file, opt.forward_end_file, opt.reverse_end_file);
	}
	
	// no problematic sites for now
	// TODO: implement the process_problematic_sites() function to read in problematic sites from a file (similar to variant sites)
	int problematic_sites[] = {};
	int number_of_problematic_sites = 0;

	int **reference_index_per_ref = (int **)malloc(opt.num_references * sizeof(int *)); // array of reference_index arrays, one per reference
	char ***resize_MSA_per_ref = (char ***)malloc(opt.num_references * sizeof(char **)); // array of MSA matricies (as 2d arrays), one per reference
	char ***resize_names_per_ref = (char ***)malloc(opt.num_references * sizeof(char **)); // array of strain name arrays, one per reference

	int *number_of_strains_remaining_per_ref = (int *)malloc(opt.num_references * sizeof(int)); // array of number of strains remaining after elimination, one per reference
	int *length_of_MSA_per_ref = (int *)malloc(opt.num_references * sizeof(int)); // array of length of MSA, one per reference

	char ***sam_results_per_ref = (char ***)malloc(opt.num_references * sizeof(char **)); // array of sam results, one per reference
	int *num_sam_lines_per_ref = (int *)malloc(opt.num_references * sizeof(int)); // array of number of sam lines in their respective results, one per reference

	int ref_idx;
	for (ref_idx = 0; ref_idx < opt.num_references; ref_idx++)
	{
		int *reference_index = (int *)malloc(FASTA_MAXLINE * sizeof(int));
		for (i = 0; i < FASTA_MAXLINE; i++)
		{
			reference_index[i] = 0;
		}	

		align_references(number_of_problematic_sites, problematic_sites, MSA_reference_paths[ref_idx], bowtie2_reference_paths[ref_idx], reference_index);

		char sam_path[1100];
		sprintf(sam_path, "%s.%d", opt.sam, ref_idx);

		perform_bowtie_alignment_xeq(opt.paired, opt.fasta_format, bowtie2_reference_paths[ref_idx], opt.forward_end_file, opt.reverse_end_file, opt.single_end_file, sam_path);
		int invoke_cleaning = calculate_error_rates(sam_path);
		if (invoke_cleaning == 1 && opt.clean_reads == 0)
		{
			clean_reads(opt.paired, opt.fasta_format, opt.single_end_file, opt.forward_end_file, opt.reverse_end_file);
		}
		perform_bowtie_alignment(opt.paired, opt.fasta_format, bowtie2_reference_paths[ref_idx], opt.forward_end_file, opt.reverse_end_file, opt.single_end_file, sam_path);

		gzFile MSA_file = Z_NULL;
		if ((MSA_file = gzopen(MSA_paths[ref_idx], "r")) == Z_NULL)
			fprintf(stderr, "MSA index %d File could not be opened.\n", ref_idx);
		int length_of_MSA = setMSALength(MSA_file);
		gzclose(MSA_file);
		printf("Length of MSA index %d: %d\n", ref_idx, length_of_MSA);

		// TODO: implement the process_problematic_sites() function to read in problematic sites from a file (similar to variant sites)
		/*int* problematic_sites = (int*)malloc(1000*sizeof(int));
		for(i=0; i<1000; i++){
			problematic_sites[i]=-1;
		}
		int number_of_problematic_sites=process_problematic_sites(problematic_sites);*/
		
		if ((MSA_file = gzopen(MSA_paths[ref_idx], "r")) == Z_NULL)
			fprintf(stderr, "MSA index %d File could not be opened.\n", ref_idx);
		int *strain_info = (int *)malloc(2 * sizeof(int));
		setNumStrains(MSA_file, strain_info);
		gzclose(MSA_file);
		int number_of_strains = strain_info[0];
		int max_name_length = strain_info[1];
		free(strain_info);
		printf("Number of strains for reference %d: %d\n", ref_idx, number_of_strains);
		printf("Maxname length for reference %d: %d\n", ref_idx, max_name_length);

		char **MSA = (char **)malloc(number_of_strains * sizeof(char *));
		for (i = 0; i < number_of_strains; i++)
		{
			MSA[i] = (char *)malloc((length_of_MSA + 1) * sizeof(char));
		}
		char **names_of_strains = (char **)malloc(number_of_strains * sizeof(char *));
		for (i = 0; i < number_of_strains; i++)
		{
			names_of_strains[i] = (char *)malloc((max_name_length + 1) * sizeof(char));
		}

		printf("Reading in MSA index %d...\n", ref_idx);
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		if ((MSA_file = gzopen(MSA_paths[ref_idx], "r")) == Z_NULL)
			fprintf(stderr, "MSA index %d File could not be opened.\n", ref_idx);
		char *reference_name = read_fasta_header_name(MSA_reference_paths[ref_idx]);
		readInMSA(MSA_file, MSA, names_of_strains, length_of_MSA, reference_name);
		free(reference_name);
		gzclose(MSA_file);
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));

		// TODO: eventually implement the removal of identical strains for speedup/optimization
		/*int* identical = (int *)malloc(number_of_strains*sizeof(int));
		int number_of_identical_strains=number_of_strains;
		if (opt.remove_identical==0){
			for(i=0; i<number_of_strains; i++){
				identical[i]=i;
			}
		}*/
		// if (opt.remove_identical==1){
		//	printf("Finding identical sequences\n");
		//	clock_gettime(CLOCK_MONOTONIC, &tstart);
		//	for(i=0; i<number_of_strains; i++){
		//		identical[i]=-1;
		//	}
		//	number_of_identical_strains=removeIdenticalStrains(number_of_strains,length_of_MSA,MSA,identical,names_of_strains,max_name_length);
		//	clock_gettime(CLOCK_MONOTONIC, &tend);
		//	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		// }

		FILE *variant_sites_file;
		int *variant_sites;
		int *number_of_variant_sites_p = (int *)malloc(sizeof(int));
		if ((variant_sites_file = fopen(variant_paths[ref_idx], "r")) == (FILE *)NULL)
			fprintf(stderr, "variants index %d File could not be opened.\n", ref_idx);
		variant_sites = readVariantSitesFile(variant_sites_file, number_of_variant_sites_p);
		fclose(variant_sites_file);
		int this_number_of_variant_sites = number_of_variant_sites_p[0];
		free(number_of_variant_sites_p);
		printf("Number of variant sites for reference %d: %d\n", ref_idx, this_number_of_variant_sites);

		printf("Calculating allele frequencies for reference %d...\n", ref_idx);
		FILE *sam_file;
		if ((sam_file = fopen(sam_path, "r")) == (FILE *)NULL)
			fprintf(stderr, "sam index %d File could not be opened.\n", ref_idx);
		double **allele_frequency = (double **)malloc(length_of_MSA * sizeof(double *));
		for (i = 0; i < length_of_MSA; i++)
		{
			allele_frequency[i] = (double *)malloc(4 * sizeof(double));
			for (j = 0; j < 4; j++)
			{
				allele_frequency[i][j] = 0;
			}
		}
		int number_of_strains_remaining = 0;
		int *max_sam_length = (int *)malloc(2 * sizeof(int));
		max_sam_length[0] = 0;
		max_sam_length[1] = 0;

		char **sam_results = NULL;
		int sam_lines = 0;

		if (opt.paired == 1)
		{
			number_of_strains_remaining = calculateAlleleFreq_paired(sam_file, allele_frequency, length_of_MSA, MSA, number_of_strains, names_of_strains, opt.freq, max_name_length, tstart, tend, this_number_of_variant_sites, variant_sites, opt.coverage, reference_index, opt.min_strains, opt.max_strains, opt.print_counts, max_sam_length, opt.print_deletions, opt.deletion_threshold, &sam_results, &sam_lines);
		}
		else
		{
			number_of_strains_remaining = calculateAlleleFreq(sam_file, allele_frequency, length_of_MSA, MSA, number_of_strains, names_of_strains, opt.freq, max_name_length, tstart, tend, this_number_of_variant_sites, variant_sites, opt.coverage, reference_index, opt.min_strains, opt.max_strains, opt.print_counts, max_sam_length, opt.print_deletions, opt.deletion_threshold, &sam_results, &sam_lines);
		}

		fclose(sam_file);
		for (i = 0; i < length_of_MSA; i++)
		{
			free(allele_frequency[i]);
		}
		free(allele_frequency);
		// free(identical);

		int *strains_kept = (int *)malloc(number_of_strains_remaining * sizeof(int));
		int placement = 0;
		for (i = 0; i < number_of_strains; i++)
		{
			if (names_of_strains[i][0] != '\0')
			{
				strains_kept[placement] = i;
				placement++;
			}
		}

		// reallocate_memory() writes into the GLOBAL resize_MSA/resize_names_of_strains
		// (existing convention -- see writeMismatchMatrix_paired). Allocate those
		// globals for this subtype, call it, then stash the globals into this
		// subtype's slot before the next loop iteration reuses them.
		resize_names_of_strains = (char **)malloc(number_of_strains_remaining * sizeof(char *));
		resize_MSA = (char **)malloc(number_of_strains_remaining * sizeof(char *));
		for (i = 0; i < number_of_strains_remaining; i++)
		{
			resize_names_of_strains[i] = (char *)malloc((max_name_length + 1) * sizeof(char));
			resize_MSA[i] = (char *)malloc((length_of_MSA + 1) * sizeof(char));
			for (j = 0; j < length_of_MSA + 1; j++)
			{
				resize_MSA[i][j] = '\0';
			}
			for (j = 0; j < max_name_length + 1; j++)
			{
				resize_names_of_strains[i][j] = '\0';
			}
		}
		reallocate_memory(number_of_strains_remaining, strains_kept, MSA, names_of_strains, max_name_length, number_of_strains);
		free(strains_kept);
		free(max_sam_length);
 
		if (number_of_strains_remaining == 0)
		{
			printf("No strains remaining for reference %d, exiting...\n", ref_idx);
			exit(1);
		}
 
		reference_index_per_ref[ref_idx] = reference_index;
		resize_MSA_per_ref[ref_idx] = resize_MSA;
		resize_names_per_ref[ref_idx] = resize_names_of_strains;
		number_of_strains_remaining_per_ref[ref_idx] = number_of_strains_remaining;
		length_of_MSA_per_ref[ref_idx] = length_of_MSA;
		sam_results_per_ref[ref_idx] = sam_results;
		num_sam_lines_per_ref[ref_idx] = sam_lines;
	}



	printf("Creating mismatch matrix...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	FILE *outfile;
	if ((outfile = fopen(opt.outfile, "w")) == (FILE *)NULL)
		fprintf(stderr, "File could not be opened.\n");
	if ((sam_file = fopen(opt.sam, "r")) == (FILE *)NULL)
		fprintf(stderr, "File could not be opened.\n");
	// --- Threaded paired-end mode: load all SAM lines into memory, split them
	// across opt.number_of_cores threads (each running writeMismatchMatrix_paired),
	// then stitch each thread's output rows back together in order. ---
	if (opt.paired == 1 && opt.no_read_bam == 0)
	{
		sam_results = (char **)malloc(max_sam_length[1] * sizeof(char *));
		for (i = 0; i < max_sam_length[1]; i++)
		{
			sam_results[i] = (char *)malloc((max_sam_length[0] + 1) * sizeof(char));
			for (j = 0; j < max_sam_length[0] + 1; j++)
			{
				sam_results[i][j] = '\0';
			}
		}
		readInSamFile(sam_file);
		fclose(sam_file);
		fprintf(outfile, "qName\tblockSizes");
		for (i = 0; i < number_of_strains_remaining; i++)
		{
			fprintf(outfile, "\t%s", resize_names_of_strains[i]);
		}
		fprintf(outfile, "\n");
		pthread_t threads[opt.number_of_cores];	 // array of our threads
		ThreadStruct thread_str[opt.number_of_cores]; // array of stuct that contains input and output for each thread
		for (i = 0; i < opt.number_of_cores; i++)
		{
			thread_str[i].sam_line_start = 0;
			thread_str[i].sam_line_end = 0;
			thread_str[i].results_str = malloc(sizeof(struct ResultsStruct));
			thread_str[i].number_of_strains_remaining = number_of_strains_remaining;
			thread_str[i].number_of_strains = number_of_strains;
			thread_str[i].length_of_MSA = length_of_MSA;
		}
		// Split sam_results into opt.number_of_cores roughly-equal chunks, then
		// nudge each chunk's boundary (via adjust_start_end) so no thread's
		// range splits a read pair across two threads.
		int divideFile, start, end;
		divideFile = max_sam_length[1] / opt.number_of_cores;
		j = 0;
		for (i = 0; i < opt.number_of_cores; i++)
		{
			start = j;
			end = j + divideFile;
			if (i == opt.number_of_cores - 1)
			{
				end = max_sam_length[1] - 1;
			}
			int *adjust_ends = (int *)malloc(2 * sizeof(int));
			if (thread_str[i].sam_line_start == 1)
			{
				start++;
			}
			adjust_ends[0] = 0;
			adjust_ends[1] = 0;
			adjust_start_end(start, end, adjust_ends, max_sam_length[0]);
			if (start != adjust_ends[0])
			{
				thread_str[i - 1].sam_line_end = thread_str[i - 1].sam_line_end - 1;
			}
			if (end != adjust_ends[1])
			{
				thread_str[i + 1].sam_line_start = thread_str[i + 1].sam_line_start + 1;
			}
			if (i == opt.number_of_cores - 1)
			{
				adjust_ends[1] = max_sam_length[1];
			}
			thread_str[i].sam_line_start = adjust_ends[0];
			thread_str[i].sam_line_end = adjust_ends[1];
			thread_str[i].thread_index = i;
			thread_str[i].max_sam_line_length = max_sam_length[0];
			j = j + divideFile;
			free(adjust_ends);
		}
		for (i = 0; i < opt.number_of_cores; i++)
		{
			printf("thread: %d start: %d end: %d\n", thread_str[i].thread_index, thread_str[i].sam_line_start, thread_str[i].sam_line_end);
			thread_str[i].results_str->mismatch = (char **)malloc((thread_str[i].sam_line_end - thread_str[i].sam_line_start) * (sizeof(char *)));
			for (k = 0; k < thread_str[i].sam_line_end - thread_str[i].sam_line_start; k++)
			{
				thread_str[i].results_str->mismatch[k] = (char *)malloc((max_sam_length[0] + 300000) * sizeof(char));
				for (l = 0; l < max_sam_length[0] + 300000; l++)
				{
					thread_str[i].results_str->mismatch[k][l] = '\0';
				}
			}
		}
		for (i = 0; i < opt.number_of_cores; i++)
		{
			pthread_create(&threads[i], NULL, writeMismatchMatrix_paired, &thread_str[i]);
		}
		for (i = 0; i < opt.number_of_cores; i++)
		{
			pthread_join(threads[i], NULL);
		}
		for (i = 0; i < opt.number_of_cores; i++)
		{
			for (j = 0; j < (thread_str[i].sam_line_end - thread_str[i].sam_line_start); j++)
			{
				if (thread_str[i].results_str->mismatch[j][0] == '\0')
				{
					break;
				}
				fprintf(outfile, "%s\n", thread_str[i].results_str->mismatch[j]);
			}
		}
		for (i = 0; i < max_sam_length[1]; i++)
		{
			free(sam_results[i]);
		}
		free(sam_results);
		for (i = 0; i < opt.number_of_cores; i++)
		{
			for (j = 0; j < thread_str[i].sam_line_end - thread_str[i].sam_line_start; j++)
			{
				free(thread_str[i].results_str->mismatch[j]);
			}
			free(thread_str[i].results_str->mismatch);
			free(thread_str[i].results_str);
		}
	}
	else if (opt.paired == 1 && opt.no_read_bam == 1)
	{
		// Paired-end, low-memory mode: single-threaded, streams the SAM file directly.
		writeMismatchMatrix_paired_no_read_bam(outfile, sam_file, resize_MSA, length_of_MSA, number_of_strains, number_of_strains_remaining, resize_names_of_strains, reference_index);
	}
	else
	{
		// Single-end mode: no threading, no mate-pair overlap handling needed.
		writeMismatchMatrix(outfile, sam_file, MSA, strains_kept, length_of_MSA, number_of_strains, number_of_strains_remaining);
	}
	fclose(outfile);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fseconds\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	// free(reference);
	// free(variant_sites);
	// for(i=0; i<number_of_strains; i++){
	// free(names_of_strains[i]);
	//	free(MSA[i]);
	//}
	// free(names_of_strains);
	// free(MSA);
	for (i = 0; i < number_of_strains_remaining; i++)
	{
		free(resize_names_of_strains[i]);
		free(resize_MSA[i]);
	}
	free(resize_names_of_strains);
	free(resize_MSA);
	if (number_of_strains_remaining > 0)
	{
		free(strains_kept);
	}
	else
	{
		printf("No strains remaining exiting...\n");
		exit(1);
	}
	char *buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(buffer, '\0', FASTA_MAXLINE);
	// FIXME: Remove the hardcoding here 
	if (opt.llr == 1)
	{
		sprintf(buffer, "Rscript ../EM_C_LLR.R -i %s -f %lf -e %lf -l -s -v %s -r %s -b %s", opt.outfile, opt.freq, opt.error, variant_paths[0], MSA_paths[0], opt.print_counts);
	}
	else
	{
		sprintf(buffer, "Rscript ../EM_C_LLR.R -i %s -f %lf -e %lf -s -v %s -r %s -b %s", opt.outfile, opt.freq, opt.error, variant_paths[0], MSA_paths[0], opt.print_counts);
	}
	system(buffer);
	free(buffer);
}