#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "bowtie_alignment.h"
#include "file_utils.h"
#include "global.h"

/**
 * @brief 
 * 
 * @param sam_filepath 
 * @param end_region_length 
 * @param end_region_error_mult 
 * @return int 
 */
int calculate_error_rates(char *sam_filepath, int end_region_length, double end_region_error_mult)
{
	FILE *sam_file = fopen(sam_filepath, "r");
	if (sam_file == NULL)
	{
		fprintf(stderr, "Error: could not open '%s' to calculate error rates.\n", sam_filepath);
		exit(1);
	}

	char buffer[FASTA_MAXLINE];
	char *token;
	int cigar_vals[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	int mismatches_ends = 0;
	int mismatches_middle = 0;
	int total_aligned_reads = 0;
	int total_length = 0;
	int invoke_cleaning = 0;

	int i, j;
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar_vals[i] = 0;
		cigar_chars[i] = '\0';
	}

	while (fgets(buffer, FASTA_MAXLINE, sam_file) != NULL)
	{
		if (buffer[0] != '@')
		{
			total_aligned_reads++;
			char *buffer_copy = strdup(buffer);

			token = strtok(buffer, "\t");
			for (i = 0; i < 5; i++)
			{
				token = strtok(NULL, "\t");
			}

			char *cigar_string = strdup(token);
			char *cigar_string_copy = strdup(cigar_string);
			char *res = strtok(cigar_string, "=XID");

			int num_cigar_chars = 0;
			while (res)
			{
				int from = res - cigar_string + strlen(res);
				int cigar_count = 0;
				sscanf(res, "%d", &cigar_count);
				res = strtok(NULL, "=XID");
				char cigar_char = '\0';
				sscanf(cigar_string_copy + from, "%c", &cigar_char);
				cigar_vals[num_cigar_chars] = cigar_count;
				cigar_chars[num_cigar_chars] = cigar_char;
				num_cigar_chars++;
			}
			free(cigar_string_copy);
			free(cigar_string);

			token = strtok(buffer_copy, "\t");
			for (i = 0; i < 9; i++)
			{
				token = strtok(NULL, "\t");
			}
			char *sequence = token;
			int sequence_length = strlen(sequence);
			total_length = total_length + sequence_length;

			int position = 0;
			for (i = 0; i < num_cigar_chars; i++)
			{
				for (j = 0; j < cigar_vals[i]; j++)
				{
					if (cigar_chars[i] == '=' || cigar_chars[i] == 'X' || cigar_chars[i] == 'I')
					{
						position++;
					}
					if (cigar_chars[i] == 'X' && position < end_region_length)
					{
						mismatches_ends++;
					}
					else if (cigar_chars[i] == 'X' && position > sequence_length - end_region_length - 1)
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
	fclose(sam_file);

	if (total_aligned_reads == 0)
	{
		fprintf(stderr, "Warning: no aligned reads found in '%s' -- skipping error-rate check.\n", sam_filepath);
		return 0;
	}

	double average_size = (double)total_length / total_aligned_reads;
	double error_rate_ends = (double)mismatches_ends / total_aligned_reads;
	error_rate_ends = error_rate_ends / (2 * end_region_length);
	double error_rate_middle = (double)mismatches_middle / total_aligned_reads;
	error_rate_middle = error_rate_middle / (average_size - (2 * end_region_length));

	printf("Error rate ends: %lf\n", error_rate_ends);
	printf("Error rate middle: %lf\n", error_rate_middle);
	if (error_rate_ends > end_region_error_mult * error_rate_middle)
	{
		printf("The error rate in your 5' and 3' ends of your reads is %.1fx larger than the error rate in the middle of your reads. Your reads need cleaning... quality filtering and trimming ends\n", end_region_error_mult);
		invoke_cleaning = 1;
	}
	return invoke_cleaning;
}

/**
 * @brief 
 * 
 * @param bowtie2_reference_path 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param sam_results_filepath 
 * @param working_dir 
 * @param using_paired_end_reads 
 * @param using_fasta_format 
 */
void perform_bowtie_alignment(char *bowtie2_reference_path, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *sam_results_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format)
{
	char *buffer = (char *)calloc(FASTA_MAXLINE, sizeof(char));

	char *index_prefix = get_filepath_in_working_dir(bowtie2_reference_path, working_dir);

	char index_check_path[1150];
	sprintf(index_check_path, "%s.1.bt2", index_prefix);
	if (access(index_check_path, F_OK) == 0)
	{
		printf("Bowtie2 index '%s' already exists. Not rebuilding.\n", index_prefix);
	}
	else
	{
		sprintf(buffer, "bowtie2-build -f %s %s", bowtie2_reference_path, index_prefix);
		system(buffer);
	}

	if (using_paired_end_reads && using_fasta_format)
	{
		sprintf(buffer, "bowtie2 --all -f -x %s -1 %s -2 %s -S %s", index_prefix, forward_end_filepath, reverse_end_filepath, sam_results_filepath);
	}
	else if (!using_paired_end_reads && using_fasta_format)
	{
		sprintf(buffer, "bowtie2 --all -f -x %s -U %s -S %s", index_prefix, single_end_filepath, sam_results_filepath);
	}
	else if (using_paired_end_reads && !using_fasta_format)
	{
		sprintf(buffer, "bowtie2 --all -x %s -1 %s -2 %s -S %s", index_prefix, forward_end_filepath, reverse_end_filepath, sam_results_filepath);
	}
	else
	{
		sprintf(buffer, "bowtie2 --all -x %s -U %s -S %s", index_prefix, single_end_filepath, sam_results_filepath);
		
	}
	system(buffer);
	free(buffer);
	free(index_prefix);
}

/**
 * @brief 
 * 
 * @param bowtie2_reference_path 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param sam_results_filepath 
 * @param working_dir 
 * @param using_paired_end_reads 
 * @param using_fasta_format 
 */
void perform_bowtie_alignment_xeq(char *bowtie2_reference_path, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *sam_results_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format)
{
	char *buffer = (char *)calloc(FASTA_MAXLINE, sizeof(char));

	char *index_prefix = get_filepath_in_working_dir(bowtie2_reference_path, working_dir);

	char index_check_path[1150];
	sprintf(index_check_path, "%s.1.bt2", index_prefix);
	if (access(index_check_path, F_OK) == 0)
	{
		printf("Bowtie2 index '%s' already exists. Not rebuilding.\n", index_prefix);
	}
	else
	{
		sprintf(buffer, "bowtie2-build -f %s %s", bowtie2_reference_path, index_prefix);
		system(buffer);
	}

	if (using_paired_end_reads && using_fasta_format)
	{
		sprintf(buffer, "bowtie2 --all --xeq -f -x %s -1 %s -2 %s -S %s", index_prefix, forward_end_filepath, reverse_end_filepath, sam_results_filepath);
	}
	else if (!using_paired_end_reads && using_fasta_format)
	{
		sprintf(buffer, "bowtie2 --all --xeq -f -x %s -U %s -S %s", index_prefix, single_end_filepath, sam_results_filepath);
	}
	else if (using_paired_end_reads && !using_fasta_format)
	{
		sprintf(buffer, "bowtie2 --all --xeq -x %s -1 %s -2 %s -S %s", index_prefix, forward_end_filepath, reverse_end_filepath, sam_results_filepath);
	}
	else
	{
		sprintf(buffer, "bowtie2 --all --xeq -x %s -U %s -S %s", index_prefix, single_end_filepath, sam_results_filepath);
		
	}
	system(buffer);
	free(buffer);
	free(index_prefix);
}
