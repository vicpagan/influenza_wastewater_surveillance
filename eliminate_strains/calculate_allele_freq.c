#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "msa.h"
#include "calculate_allele_freq.h"

// TODO: Make it so that we calculate allele frequency based on whichever reference strain the read aligns with the best

void calculate_allele_freq_paired(double **allele, MSA *msa_str, double freq_threshold, struct timespec tstart, struct timespec tend, int coverage, int min_strains_remaining, int max_strains_remaining, int output_allele_counts, char *allele_counts_output_filepath, int output_deletions, char *deletions_output_filepath, double deletion_threshold, ReferenceData *reference_data_strs, int num_references)
{
	int i, j, k, sam_line_idx, msa_seq_idx;
	int msa_sequence_length = msa_str->sequence_length;

	char first_copy[FASTA_MAXLINE];
	char second_copy[FASTA_MAXLINE];
	char first_sam_line_cigar[FASTA_MAXLINE];
	char second_sam_line_cigar[FASTA_MAXLINE];
	char *first_sam_fields[11];
	char *second_sam_fields[11];

	char *first_sequence = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
	int first_cigar_vals[MAX_CIGAR] = {0};
	char first_cigar_chars[MAX_CIGAR] = {'\0'};

	char *second_sequence = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
	int second_cigar_vals[MAX_CIGAR] = {0};
	char second_cigar_chars[MAX_CIGAR] = {'\0'};

	int *deletions = (int *)calloc((msa_sequence_length), sizeof(int));

	int first_msa_positions[MAX_READ_LENGTH];
	char first_bases[MAX_READ_LENGTH];
	int second_msa_positions[MAX_READ_LENGTH];
	char second_bases[MAX_READ_LENGTH];
	int merged_msa_positions[2 * MAX_READ_LENGTH];
	char merged_bases[2 * MAX_READ_LENGTH];

	SAMResults sam_results_strs[num_references];
	int *reference_indexes[num_references];
	for (i = 0; i < num_references; i++)
	{
		sam_results_strs[i] = reference_data_strs[i].sam_results_str;
		reference_indexes[i] = reference_data_strs[i].reference_index;
	}
	int num_sam_lines = sam_results_strs[0].num_sam_lines;

	clock_gettime(CLOCK_MONOTONIC, &tstart);
	// --- Stage 1: tally per-site A/C/G/T counts (and deletions) from every read pair ---
	for (sam_line_idx = 0; sam_line_idx < num_sam_lines; sam_line_idx += 2)
	{
		int skip_reads = 0;

		int lowest_total_nm = INT32_MAX;
		int best_sam_index = -1;
		for (j = 0; j < num_references; j++)
		{
			strcpy(first_copy, sam_results_strs[j].sam_results[sam_line_idx]);
			strcpy(second_copy, sam_results_strs[j].sam_results[sam_line_idx + 1]);

			char *first_token = strtok(first_copy, "\t");
			char *second_token = strtok(second_copy, "\t");

			int first_nm = -1;
			int second_nm = -1;
			int total_nm = -1;

			while (first_token != NULL && first_nm == -1)
			{
				if (strncmp(first_token, "NM:i:", 5) == 0)
				{
					first_nm = atoi(first_token + 5);
				}
				first_token = strtok(NULL, "\t");
			}

			while (second_token != NULL && second_nm == -1)
			{
				if (strncmp(second_token, "NM:i:", 5) == 0)
				{
					second_nm = atoi(second_token + 5);
				}
				second_token = strtok(NULL, "\t");
			}

			total_nm = first_nm + second_nm;
			if (total_nm < lowest_total_nm && first_nm != -1 && second_nm != -1)
			{
				lowest_total_nm = total_nm;
				best_sam_index = j;
			}
		}

		if (best_sam_index != -1)
		{
			strcpy(first_copy, sam_results_strs[best_sam_index].sam_results[sam_line_idx]);
			strcpy(second_copy, sam_results_strs[best_sam_index].sam_results[sam_line_idx + 1]);

			char *first_token = strtok(first_copy, "\t");
			char *second_token = strtok(second_copy, "\t");

			j = 0;
			while (first_token != NULL && j < 11)
			{
				first_sam_fields[j] = first_token;
				j++;
				first_token = strtok(NULL, "\t");
			}

			j = 0;
			while (second_token != NULL && j < 11)
			{
				second_sam_fields[j] = second_token;
				j++;
				second_token = strtok(NULL, "\t");
			}

			int first_sequence_start_pos = atoi(first_sam_fields[3]) - 1;
			int second_sequence_start_pos = atoi(second_sam_fields[3]) - 1;

			strcpy(first_sam_line_cigar, first_sam_fields[5]);

			int first_num_cigar_chars = 0;
			int cigar_counter = 0;
			int index = 0;
			for (j = 0; first_sam_line_cigar[j] != '\0'; j++)
			{
				if (isdigit(first_sam_line_cigar[j]))
				{
					int digit = first_sam_line_cigar[j] - '0';
					cigar_counter = cigar_counter * 10 + digit;
				}
				else
				{
					first_cigar_vals[index] = cigar_counter;
					first_cigar_chars[index] = first_sam_line_cigar[j];
					first_num_cigar_chars++;

					index++;
					cigar_counter = 0;
				}
			}
			strcpy(first_sequence, first_sam_fields[9]);

			strcpy(second_sam_line_cigar, second_sam_fields[5]);

			int second_num_cigar_chars = 0;
			cigar_counter = 0;
			index = 0;
			for (j = 0; second_sam_line_cigar[j] != '\0'; j++)
			{
				if (isdigit(second_sam_line_cigar[j]))
				{
					int digit = second_sam_line_cigar[j] - '0';
					cigar_counter = cigar_counter * 10 + digit;
				}
				else
				{
					second_cigar_vals[index] = cigar_counter;
					second_cigar_chars[index] = second_sam_line_cigar[j];
					second_num_cigar_chars++;

					index++;
					cigar_counter = 0;
				}
			}
			strcpy(second_sequence, second_sam_fields[9]);

			int first_index = 0;
			int first_sequence_offset = 0;
			int first_reference_offset = 0;
			for (j = 0; j < first_num_cigar_chars; j++)
			{
				for (k = 0; k < first_cigar_vals[j]; k++)
				{
					int first_pos_in_msa = reference_indexes[best_sam_index][first_sequence_start_pos + first_reference_offset];

					if (first_cigar_chars[j] == 'M')
					{
						if (first_pos_in_msa != -1)
						{
							first_msa_positions[first_index] = first_pos_in_msa;
							first_bases[first_index] = toupper(first_sequence[first_sequence_offset]);
							first_index++;
						}
						first_sequence_offset++;
						first_reference_offset++;
					}
					else if (first_cigar_chars[j] == 'I' || first_cigar_chars[j] == 'S')
					{
						first_sequence_offset++;
					}
					else if (first_cigar_chars[j] == 'D')
					{
						if (first_pos_in_msa != -1)
						{
							deletions[first_pos_in_msa]++;
						}
						first_reference_offset++;
					}
				}
			}

			int second_index = 0;
			int second_sequence_offset = 0;
			int second_reference_offset = 0;
			for (j = 0; j < second_num_cigar_chars; j++)
			{
				for (k = 0; k < second_cigar_vals[j]; k++)
				{
					int second_pos_in_msa = reference_indexes[best_sam_index][second_sequence_start_pos + second_reference_offset];

					if (second_cigar_chars[j] == 'M')
					{
						if (second_pos_in_msa != -1)
						{
							second_msa_positions[second_index] = second_pos_in_msa;
							second_bases[second_index] = toupper(second_sequence[second_sequence_offset]);
							second_index++;
						}
						second_sequence_offset++;
						second_reference_offset++;
					}
					else if (second_cigar_chars[j] == 'I' || second_cigar_chars[j] == 'S')
					{
						second_sequence_offset++;
					}
					else if (second_cigar_chars[j] == 'D')
					{
						if (second_pos_in_msa != -1)
						{
							deletions[second_pos_in_msa]++;
						}
						second_reference_offset++;
					}
				}
			}

			int merged_index = 0;
			j = 0;
			k = 0;
			while (j < first_index && k < second_index)
			{
				if (first_msa_positions[j] < second_msa_positions[k])
				{
					merged_msa_positions[merged_index] = first_msa_positions[j];
					merged_bases[merged_index] = first_bases[j];
					merged_index++;
					j++;
				}
				else if (first_msa_positions[j] > second_msa_positions[k])
				{
					merged_msa_positions[merged_index] = second_msa_positions[k];
					merged_bases[merged_index] = second_bases[k];
					merged_index++;
					k++;
				}
				else
				{
					if (first_bases[j] == second_bases[k])
					{
						merged_msa_positions[merged_index] = first_msa_positions[j];
						merged_bases[merged_index] = first_bases[j];
						merged_index++;
					}
					j++;
					k++;
				}
			}
			while (j < first_index)
			{
				merged_msa_positions[merged_index] = first_msa_positions[j];
				merged_bases[merged_index] = first_bases[j];
				merged_index++;
				j++;
			}
			while (k < second_index)
			{
				merged_msa_positions[merged_index] = second_msa_positions[k];
				merged_bases[merged_index] = second_bases[k];
				merged_index++;
				k++;
			}

			for (j = 0; j < merged_index; j++)
			{
				if (merged_msa_positions[j] < msa_sequence_length)
				{
					switch (merged_bases[j])
					{
						case 'A':
							allele[merged_msa_positions[j]][0]++;
							break;
						case 'G':
							allele[merged_msa_positions[j]][1]++;
							break;
						case 'C':
							allele[merged_msa_positions[j]][2]++;
							break;
						case 'T':
							allele[merged_msa_positions[j]][3]++;
							break;
					}
				}
			}
		}
	}
	free(first_sequence);
	free(second_sequence);
	
	if (output_deletions)
	{
		FILE *deletion_sites_file;
		if ((deletion_sites_file = fopen(deletions_output_filepath, "w")) == (FILE *)NULL)
		{
			fprintf(stderr, "Deletion sites output file could not be opened.\n");
		}
		fprintf(deletion_sites_file, "Site\tFrequency\n");
		for (i = 0; i < msa_sequence_length; i++)
		{
			double deletion_ratio = (double)(deletions[i]) / num_sam_lines;
			if (deletion_ratio > deletion_threshold)
			{
				fprintf(deletion_sites_file, "%d\t%lf\n", i, deletion_ratio);
			}
		}
		fclose(deletion_sites_file);
	}
	free(deletions);

	if (output_allele_counts)
	{
		FILE *allele_counts_file;
		if ((allele_counts_file = fopen(allele_counts_output_filepath, "w")) == (FILE *)NULL)
		{
			fprintf(stderr, "Allele Counts output file could not be opened.\n");
		}	
		fprintf(allele_counts_file, "position\tA\tG\tC\tT\n");
		for (i = 0; i < msa_sequence_length; i++)
		{
			fprintf(allele_counts_file, "%d\t%lf\t%lf\t%lf\t%lf\n", i, allele[i][0], allele[i][1], allele[i][2], allele[i][3]);
		}
		fclose(allele_counts_file);
	}

	// --- Stage 2: convert counts to frequencies, drop low-coverage sites, mark "bad" bases ---
	int num_covered_sites = 0;
	for (i = 0; i < msa_sequence_length; i++)
	{
		double total = 0;
		for (j = 0; j < 4; j++)
		{
			total = total + allele[i][j];
		}
		if (total > 0)
		{
			num_covered_sites++;
		}
	}

	int *covered_sites = (int *)malloc(num_covered_sites * sizeof(int));
	for (i = 0; i < num_covered_sites; i++)
	{
		covered_sites[i] = -1;
	}

	k = 0;
	for (i = 0; i < msa_sequence_length; i++)
	{
		double total = 0;
		for (j = 0; j < 4; j++)
		{
			total = total + allele[i][j];
		}
		if (total >= coverage)
		{
			covered_sites[k] = i;
			k++;
		}
		for (j = 0; j < 4; j++)
		{
			allele[i][j] = allele[i][j] / total;
		}
	}
	printf("Number of sites not covered: %d\n", msa_sequence_length - k);

	int *bad_bases_count = (int *)calloc(msa_sequence_length, sizeof(int));
	char **bad_base_char = (char **)malloc(msa_sequence_length * sizeof(char *));
	for (i = 0; i < msa_sequence_length; i++)
	{
		bad_base_char[i] = (char *)malloc(4 * sizeof(char));

		j = 0;
		for (k = 0; k < 4; k++)
		{
			if (allele[i][k] < freq_threshold)
			{
				bad_bases_count[i]++;
			}
			else
			{
				switch (k)
				{
					case 0: 
						bad_base_char[i][j] = 'A'; 
						break;
					case 1: 
						bad_base_char[i][j] = 'G'; 
						break;
					case 2: 
						bad_base_char[i][j] = 'C'; 
						break;
					case 3: 
						bad_base_char[i][j] = 'T'; 
						break;
				}
				j++;
			}
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tstart);

	// --- Stage 3: iteratively eliminate strains incompatible with the "bad" bases above ---
	int num_sequences = msa_str->num_sequences;
	
	int *incompat_counter = (int *)calloc(num_sequences, sizeof(int));
	int *sequences_to_remove = (int *)calloc(num_sequences, sizeof(int));

	int num_sequences_remaining;
	int count;
	int run_loop = 1;
	int num_loop_iterations = 1;
	while (run_loop)
	{
		printf("iteration %d\n", num_loop_iterations);
		for (i = 0; i < num_sequences; i++)
		{
			incompat_counter[i] = 0;
		}
		num_sequences_remaining = num_sequences;
		for (i = 0; i < num_covered_sites; i++)
		{
			for (j = 0; j < num_sequences; j++)
			{
				if (incompat_counter[j] < num_loop_iterations)
				{
					count = 0;
					for (k = 0; k < 4 - bad_bases_count[covered_sites[i]]; k++)
					{
						if (msa_str->sequences[j][covered_sites[i]] != bad_base_char[covered_sites[i]][k])
						{
							count++;
						}
					}
					if (count == (4 - bad_bases_count[covered_sites[i]]))
					{
						incompat_counter[j]++;
					}
				}
			}
		}
		for (i = 0; i < num_sequences; i++)
		{
			if (incompat_counter[i] == num_loop_iterations)
			{
				num_sequences_remaining--;
			}
		}
		if (num_sequences_remaining >= min_strains_remaining && num_sequences_remaining < max_strains_remaining)
		{
			printf("exiting loop. %d remaining\n", num_sequences_remaining);
			run_loop = 0;
		}
		else if (num_sequences_remaining >= max_strains_remaining)
		{
			printf("%d strains remaining. exiting...\n", num_sequences_remaining);
			exit(1);
		}
		else
		{
			printf("there are %d remaining... \n", num_sequences_remaining);
			num_loop_iterations++;
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));

	num_sequences_remaining = num_sequences;
	for (i = 0; i < num_sequences; i++)
	{
		if (incompat_counter[i] == num_loop_iterations)
		{
			sequences_to_remove[i] = 1;
			num_sequences_remaining--;
		}
	}
	printf("Number remaining: %d\n", num_sequences_remaining);

	if (num_sequences != num_sequences_remaining)
	{
		prune_msa_sequences(msa_str, sequences_to_remove);
	}

	for (i = 0; i < msa_sequence_length; i++)
	{
		free(bad_base_char[i]);
	}
	free(bad_base_char);
	free(bad_bases_count);
	free(covered_sites);
	free(incompat_counter);
}


void calculate_allele_freq_single(double **allele, MSA *msa_str, double freq_threshold, struct timespec tstart, struct timespec tend, int coverage, int min_strains_remaining, int max_strains_remaining, int output_allele_counts, char *allele_counts_output_filepath, int output_deletions, char *deletions_output_filepath, double deletion_threshold, ReferenceData *reference_data_strs, int num_references)
{
	int i, j, k, sam_line_idx, msa_seq_idx;
	int msa_sequence_length = msa_str->sequence_length;

	char copy[FASTA_MAXLINE];
	char sam_line_cigar[FASTA_MAXLINE];
	char *sam_fields[11];

	char *sequence = (char *)calloc(MAX_READ_LENGTH, sizeof(char));
	int cigar_vals[MAX_CIGAR] = {0};
	char cigar_chars[MAX_CIGAR] = {'\0'};

	clock_gettime(CLOCK_MONOTONIC, &tstart);

	int *deletions = (int *)calloc((msa_sequence_length), sizeof(int));

	int msa_positions[MAX_READ_LENGTH];
	char bases[MAX_READ_LENGTH];

	SAMResults sam_results_strs[num_references];
	int *reference_indexes[num_references];
	for (int i = 0; i < num_references; i++)
	{
		sam_results_strs[i] = reference_data_strs[i].sam_results_str;
		reference_indexes[i] = reference_data_strs[i].reference_index;
	}
	int num_sam_lines = sam_results_strs[0].num_sam_lines;

	// --- Stage 1: tally per-site A/C/G/T counts (and deletions) from every read ---
	for (sam_line_idx = 0; sam_line_idx < num_sam_lines; sam_line_idx++)
	{
		int skip_reads = 0;

		int lowest_nm = INT32_MAX;
		int best_sam_index = -1;
		for (j = 0; j < num_references; j++)
		{
			strcpy(copy, sam_results_strs[j].sam_results[sam_line_idx]);

			char *token = strtok(copy, "\t");

			int nm = -1;

			while (token != NULL && nm == -1)
			{
				if (strncmp(token, "NM:i:", 5) == 0)
				{
					nm = atoi(token + 5);
				}
				token = strtok(NULL, "\t");
			}

			if (nm < lowest_nm && nm != -1)
			{
				lowest_nm = nm;
				best_sam_index = j;
			}
		}

		if (best_sam_index != -1)
		{
			strcpy(copy, sam_results_strs[best_sam_index].sam_results[sam_line_idx]);

			char *token = strtok(copy, "\t");

			j = 0;
			while (token != NULL && j < 11)
			{
				sam_fields[j] = token;
				j++;
				token = strtok(NULL, "\t");
			}

			int sequence_start_pos = atoi(sam_fields[3]) - 1;

			strcpy(sam_line_cigar, sam_fields[5]);

			int num_cigar_chars = 0;
			int cigar_counter = 0;
			int cigar_index = 0;
			for (j = 0; sam_line_cigar[j] != '\0'; j++)
			{
				if (isdigit(sam_line_cigar[j]))
				{
					int digit = sam_line_cigar[j] - '0';
					cigar_counter = cigar_counter * 10 + digit;
				}
				else
				{
					cigar_vals[cigar_index] = cigar_counter;
					cigar_chars[cigar_index] = sam_line_cigar[j];
					num_cigar_chars++;

					cigar_index++;
					cigar_counter = 0;
				}
			}
			strcpy(sequence, sam_fields[9]);

			int msa_index = 0;
			int sequence_offset = 0;
			int reference_offset = 0;
			for (j = 0; j < num_cigar_chars; j++)
			{
				for (k = 0; k < cigar_vals[j]; k++)
				{
					int pos_in_msa = reference_indexes[best_sam_index][sequence_start_pos + reference_offset];

					if (cigar_chars[j] == 'M')
					{
						if (pos_in_msa != -1)
						{
							msa_positions[msa_index] = pos_in_msa;
							bases[msa_index] = toupper(sequence[sequence_offset]);
							msa_index++;
						}
						sequence_offset++;
						reference_offset++;
					}
					else if (cigar_chars[j] == 'I' || cigar_chars[j] == 'S')
					{
						sequence_offset++;
					}
					else if (cigar_chars[j] == 'D')
					{
						if (pos_in_msa != -1)
						{
							deletions[pos_in_msa]++;
						}
						reference_offset++;
					}
				}
			}

			for (j = 0; j < msa_index; j++)
			{
				if (msa_positions[j] < msa_sequence_length)
				{
					switch (bases[j])
					{
						case 'A':
							allele[msa_positions[j]][0]++;
							break;
						case 'G':
							allele[msa_positions[j]][1]++;
							break;
						case 'C':
							allele[msa_positions[j]][2]++;
							break;
						case 'T':
							allele[msa_positions[j]][3]++;
							break;
					}
				}
			}
		}
	}
	free(sequence);
	
	if (output_deletions)
	{
		FILE *deletion_sites_file;
		if ((deletion_sites_file = fopen(deletions_output_filepath, "w")) == (FILE *)NULL)
		{
			fprintf(stderr, "Deletion sites output file could not be opened.\n");
		}
		fprintf(deletion_sites_file, "Site\tFrequency\n");
		for (i = 0; i < msa_sequence_length; i++)
		{
			double deletion_ratio = (double)(deletions[i]) / num_sam_lines;
			if (deletion_ratio > deletion_threshold)
			{
				fprintf(deletion_sites_file, "%d\t%lf\n", i, deletion_ratio);
			}
		}
		fclose(deletion_sites_file);
	}
	free(deletions);

	if (output_allele_counts)
	{
		FILE *allele_counts_file;
		if ((allele_counts_file = fopen(allele_counts_output_filepath, "w")) == (FILE *)NULL)
		{
			fprintf(stderr, "Allele Counts output file could not be opened.\n");
		}	
		fprintf(allele_counts_file, "position\tA\tG\tC\tT\n");
		for (i = 0; i < msa_sequence_length; i++)
		{
			fprintf(allele_counts_file, "%d\t%lf\t%lf\t%lf\t%lf\n", i, allele[i][0], allele[i][1], allele[i][2], allele[i][3]);
		}
		fclose(allele_counts_file);
	}

	// --- Stage 2: convert counts to frequencies, drop low-coverage sites, mark "bad" bases ---
	int num_covered_sites = 0;
	for (i = 0; i < msa_sequence_length; i++)
	{
		double total = 0;
		for (j = 0; j < 4; j++)
		{
			total = total + allele[i][j];
		}
		if (total > 0)
		{
			num_covered_sites++;
		}
	}

	int *covered_sites = (int *)malloc(num_covered_sites * sizeof(int));
	for (i = 0; i < num_covered_sites; i++)
	{
		covered_sites[i] = -1;
	}

	k = 0;
	for (i = 0; i < msa_sequence_length; i++)
	{
		double total = 0;
		for (j = 0; j < 4; j++)
		{
			total = total + allele[i][j];
		}
		if (total >= coverage)
		{
			covered_sites[k] = i;
			k++;
		}
		for (j = 0; j < 4; j++)
		{
			allele[i][j] = allele[i][j] / total;
		}
	}
	printf("Number of sites not covered: %d\n", msa_sequence_length - k);

	int *bad_bases_count = (int *)calloc(msa_sequence_length, sizeof(int));
	char **bad_base_char = (char **)malloc(msa_sequence_length * sizeof(char *));
	for (msa_seq_idx = 0; msa_seq_idx < msa_sequence_length; msa_seq_idx++)
	{
		bad_base_char[msa_seq_idx] = (char *)malloc(4 * sizeof(char));

		j = 0;
		for (k = 0; k < 4; k++)
		{
			if (allele[msa_seq_idx][k] < freq_threshold)
			{
				bad_bases_count[msa_seq_idx]++;
			}
			else
			{
				switch (k)
				{
					case 0: 
						bad_base_char[msa_seq_idx][j] = 'A'; 
						break;
					case 1: 
						bad_base_char[msa_seq_idx][j] = 'G'; 
						break;
					case 2: 
						bad_base_char[msa_seq_idx][j] = 'C'; 
						break;
					case 3: 
						bad_base_char[msa_seq_idx][j] = 'T'; 
						break;
				}
				j++;
			}
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tstart);

	// --- Stage 3: iteratively eliminate strains incompatible with the "bad" bases above ---
	int num_sequences = msa_str->num_sequences;
	
	int *incompat_counter = (int *)calloc(num_sequences, sizeof(int));
	int *sequences_to_remove = (int *)calloc(num_sequences, sizeof(int));

	int num_sequences_remaining;
	int count;
	int run_loop = 1;
	int num_loop_iterations = 1;
	while (run_loop)
	{
		printf("iteration %d\n", num_loop_iterations);
		for (i = 0; i < num_sequences; i++)
		{
			incompat_counter[i] = 0;
		}
		num_sequences_remaining = num_sequences;
		for (i = 0; i < num_covered_sites; i++)
		{
			for (j = 0; j < num_sequences; j++)
			{
				if (incompat_counter[j] < num_loop_iterations)
				{
					count = 0;
					for (k = 0; k < 4 - bad_bases_count[covered_sites[i]]; k++)
					{
						if (msa_str->sequences[j][covered_sites[i]] != bad_base_char[covered_sites[i]][k])
						{
							count++;
						}
					}
					if (count == (4 - bad_bases_count[covered_sites[i]]))
					{
						incompat_counter[j]++;
					}
				}
			}
		}
		for (i = 0; i < num_sequences; i++)
		{
			if (incompat_counter[i] == num_loop_iterations)
			{
				num_sequences_remaining--;
			}
		}
		if (num_sequences_remaining >= min_strains_remaining && num_sequences_remaining < max_strains_remaining)
		{
			printf("exiting loop. %d remaining\n", num_sequences_remaining);
			run_loop = 0;
		}
		else if (num_sequences_remaining >= max_strains_remaining)
		{
			printf("%d strains remaining. exiting...\n", num_sequences_remaining);
			exit(1);
		}
		else
		{
			printf("there are %d remaining... \n", num_sequences_remaining);
			num_loop_iterations++;
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));

	num_sequences_remaining = num_sequences;
	for (i = 0; i < num_sequences; i++)
	{
		if (incompat_counter[i] == num_loop_iterations)
		{
			sequences_to_remove[i] = 1;
			num_sequences_remaining--;
		}
	}
	printf("Number remaining: %d\n", num_sequences_remaining);

	if (num_sequences != num_sequences_remaining)
	{
		prune_msa_sequences(msa_str, sequences_to_remove);
	}

	for (i = 0; i < msa_sequence_length; i++)
	{
		free(bad_base_char[i]);
	}
	free(bad_base_char);
	free(bad_bases_count);
	free(covered_sites);
	free(incompat_counter);
}
