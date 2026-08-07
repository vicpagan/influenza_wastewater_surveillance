#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#include <pthread.h>

#include "write_mismatch_matrix.h"
#include "sam.h"
#include "msa.h"

void *write_mismatch_matrix_paired(void *ptr)
{
	int i, j, k, ref_idx, sam_line_idx, msa_seq_idx;

	MismatchMatrixThreadStruct *thread_str = (MismatchMatrixThreadStruct *)ptr;
    int thread_index = thread_str->thread_index;
	int sam_partition_start = thread_str->sam_partition_start;
	int sam_partition_end = thread_str->sam_partition_end;
	int num_references = thread_str->num_references;

	ReferenceData *reference_data_strs = thread_str->reference_data_strs;
	SAMResults sam_results_strs[num_references];
	int *reference_indexes[num_references];

	int num_sam_lines = reference_data_strs[0].sam_results_str.num_sam_lines;
	int max_sam_line_length = reference_data_strs[0].sam_results_str.max_sam_line_length;
	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		sam_results_strs[ref_idx] = reference_data_strs[ref_idx].sam_results_str;
		reference_indexes[ref_idx] = reference_data_strs[ref_idx].reference_index;
		assert(num_sam_lines == sam_results_strs[ref_idx].num_sam_lines);
		if (sam_results_strs[ref_idx].max_sam_line_length > max_sam_line_length)
		{
			max_sam_line_length = sam_results_strs[ref_idx].max_sam_line_length;
		}
	}

	MSA *msa_str = thread_str->msa_str;
	char **msa_sequences = msa_str->sequences;
	int msa_sequence_length = msa_str->sequence_length;
	int num_msa_sequences = msa_str->num_sequences;

	char *row_buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *readname = (char *)malloc(FASTA_MAXLINE * sizeof(char));

	char *first_copy = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *second_copy = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *first_sam_line_cigar = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *second_sam_line_cigar = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *first_sam_fields[11];
	char *second_sam_fields[11];

	char *first_sequence = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	char *second_sequence = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	int *first_cigar_vals = (int *)malloc(MAX_CIGAR * sizeof(int));
	char *first_cigar_chars = (char *)malloc(MAX_CIGAR * sizeof(char));
	int *second_cigar_vals = (int *)malloc(MAX_CIGAR * sizeof(int));
	char *second_cigar_chars = (char *)malloc(MAX_CIGAR * sizeof(char));

	int *first_msa_positions = (int *)malloc(MAX_READ_LENGTH * sizeof(int));
	char *first_bases = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	int *second_msa_positions = (int *)malloc(MAX_READ_LENGTH * sizeof(int));
	char *second_bases = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	int *merged_msa_positions = (int *)malloc(2 * MAX_READ_LENGTH * sizeof(int));
	char *merged_bases = (char *)malloc(2 * MAX_READ_LENGTH * sizeof(char));

	int *current_mismatch_matrix_row = (int *)malloc(num_msa_sequences * sizeof(int));

	for (sam_line_idx = sam_partition_start; sam_line_idx < sam_partition_end; sam_line_idx = sam_line_idx + 2)
	{
		for (i = 0; i < num_msa_sequences; i++)
		{
			current_mismatch_matrix_row[i] = INT32_MAX;
		}
		int best_alignment_size = -1;

		for (ref_idx = 0; ref_idx < num_references; ref_idx++)
		{
			strcpy(first_copy, sam_results_strs[ref_idx].sam_results[sam_line_idx]);
			strcpy(second_copy, sam_results_strs[ref_idx].sam_results[sam_line_idx + 1]);

			char *first_token = strtok(first_copy, "\t");
			j = 0;
			while (first_token != NULL && j < 11)
			{
				first_sam_fields[j] = first_token;
				j++;
				first_token = strtok(NULL, "\t");
			}

			char *second_token = strtok(second_copy, "\t");
			j = 0;
			while (second_token != NULL && j < 11)
			{
				second_sam_fields[j] = second_token;
				j++;
				second_token = strtok(NULL, "\t");
			}

			if (ref_idx == 0)
			{
				strcpy(readname, first_sam_fields[0]);
			}

			int first_status = parse_sam_flags(atoi(first_sam_fields[1]));
			if (first_status == 0 || first_status == 1)
			{
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

				int current_alignment_size = 0;

				int first_index = 0;
				int first_sequence_offset = 0;
				int first_reference_offset = 0;
				for (j = 0; j < first_num_cigar_chars; j++)
				{
					for (k = 0; k < first_cigar_vals[j]; k++)
					{
						int first_pos_in_msa = reference_indexes[ref_idx][first_sequence_start_pos + first_reference_offset];

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
							current_alignment_size++;
						}
						else if (first_cigar_chars[j] == 'I' || first_cigar_chars[j] == 'S')
						{
							first_sequence_offset++;
						}
						else if (first_cigar_chars[j] == 'D')
						{
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
						int second_pos_in_msa = reference_indexes[ref_idx][second_sequence_start_pos + second_reference_offset];

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
							current_alignment_size++;
						}
						else if (second_cigar_chars[j] == 'I' || second_cigar_chars[j] == 'S')
						{
							second_sequence_offset++;
						}
						else if (second_cigar_chars[j] == 'D')
						{
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
						current_alignment_size--;
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

				int nm;
				for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
				{
					nm = 0;
					for (j = 0; j < merged_index; j++)
					{
						if (merged_msa_positions[j] < msa_sequence_length)
						{
							char msa_seq_base = msa_sequences[msa_seq_idx][merged_msa_positions[j]];
							if (msa_seq_base != merged_bases[j] && msa_seq_base != '-' && msa_seq_base != '\0')
							{
								nm++;
							}
						}
					}

					if (nm < current_mismatch_matrix_row[msa_seq_idx])
					{
						current_mismatch_matrix_row[msa_seq_idx] = nm;
					}
				}

				best_alignment_size = current_alignment_size;
			}
		}

		if (best_alignment_size != -1)
		{
			sprintf(row_buffer, "%s\t%d", readname, best_alignment_size);
			for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
			{
				char num_buffer[16];
				sprintf(num_buffer, "\t%d", current_mismatch_matrix_row[msa_seq_idx]);
				strcat(row_buffer, num_buffer);
			}

			pthread_mutex_lock(thread_str->write_mutex);
			fprintf(thread_str->outfile, "%s\n", row_buffer);
			pthread_mutex_unlock(thread_str->write_mutex);
		}
	}

	free(row_buffer);
	free(readname);
	free(first_copy);
	free(second_copy);
	free(first_sam_line_cigar);
	free(second_sam_line_cigar);
	free(first_sequence);
	free(second_sequence);
	free(first_cigar_vals);
	free(first_cigar_chars);
	free(second_cigar_vals);
	free(second_cigar_chars);
	free(first_msa_positions);
	free(first_bases);
	free(second_msa_positions);
	free(second_bases);
	free(merged_msa_positions);
	free(merged_bases);
	free(current_mismatch_matrix_row);
	return NULL;
}

void *write_mismatch_matrix_single(void *ptr)
{
	int i, j, k, ref_idx, sam_line_idx, msa_seq_idx;

	MismatchMatrixThreadStruct *thread_str = (MismatchMatrixThreadStruct *)ptr;
    int thread_index = thread_str->thread_index;
	int sam_partition_start = thread_str->sam_partition_start;
	int sam_partition_end = thread_str->sam_partition_end;
	int num_references = thread_str->num_references;

	ReferenceData *reference_data_strs = thread_str->reference_data_strs;
	SAMResults sam_results_strs[num_references];
	int *reference_indexes[num_references];

	int num_sam_lines = reference_data_strs[0].sam_results_str.num_sam_lines;
	int max_sam_line_length = reference_data_strs[0].sam_results_str.max_sam_line_length;
	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		sam_results_strs[ref_idx] = reference_data_strs[ref_idx].sam_results_str;
		reference_indexes[ref_idx] = reference_data_strs[ref_idx].reference_index;
		assert(num_sam_lines == sam_results_strs[ref_idx].num_sam_lines);
		if (sam_results_strs[ref_idx].max_sam_line_length > max_sam_line_length)
		{
			max_sam_line_length = sam_results_strs[ref_idx].max_sam_line_length;
		}
	}

	MSA *msa_str = thread_str->msa_str;
	char **msa_sequences = msa_str->sequences;
	int msa_sequence_length = msa_str->sequence_length;
	int num_msa_sequences = msa_str->num_sequences;

	char *row_buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *readname = (char *)malloc(FASTA_MAXLINE * sizeof(char));

	char *copy = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *sam_line_cigar = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	char *sam_fields[11];

	char *sequence = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	int *cigar_vals = (int *)malloc(MAX_CIGAR * sizeof(int));
	char *cigar_chars = (char *)malloc(MAX_CIGAR * sizeof(char));

	int *msa_positions = (int *)malloc(MAX_READ_LENGTH * sizeof(int));
	char *bases = (char *)malloc(MAX_READ_LENGTH * sizeof(char));

	int *current_mismatch_matrix_row = (int *)malloc(num_msa_sequences * sizeof(int));

	for (sam_line_idx = sam_partition_start; sam_line_idx < sam_partition_end; sam_line_idx++)
	{
		for (i = 0; i < num_msa_sequences; i++)
		{
			current_mismatch_matrix_row[i] = INT32_MAX;
		}
		int best_alignment_size = -1;

		for (ref_idx = 0; ref_idx < num_references; ref_idx++)
		{
			strcpy(copy, sam_results_strs[ref_idx].sam_results[sam_line_idx]);

			char *token = strtok(copy, "\t");

			j = 0;
			while (token != NULL && j < 11)
			{
				sam_fields[j] = token;
				j++;
				token = strtok(NULL, "\t");
			}

			if (ref_idx == 0)
			{
				strcpy(readname, sam_fields[0]);
			}

			int status = parse_sam_flags(atoi(sam_fields[1]));
			if (status != -1)
			{
				int sequence_start_pos = atoi(sam_fields[3]) - 1;

				strcpy(sam_line_cigar, sam_fields[5]);

				int num_cigar_chars = 0;
				int cigar_counter = 0;
				int index = 0;
				for (j = 0; sam_line_cigar[j] != '\0'; j++)
				{
					if (isdigit(sam_line_cigar[j]))
					{
						int digit = sam_line_cigar[j] - '0';
						cigar_counter = cigar_counter * 10 + digit;
					}
					else
					{
						cigar_vals[index] = cigar_counter;
						cigar_chars[index] = sam_line_cigar[j];
						num_cigar_chars++;

						index++;
						cigar_counter = 0;
					}
				}
				strcpy(sequence, sam_fields[9]);

				int current_alignment_size = 0;

				int msa_index = 0;
				int sequence_offset = 0;
				int reference_offset = 0;
				for (j = 0; j < num_cigar_chars; j++)
				{
					for (k = 0; k < cigar_vals[j]; k++)
					{
						int pos_in_msa = reference_indexes[ref_idx][sequence_start_pos + reference_offset];

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
							current_alignment_size++;
						}
						else if (cigar_chars[j] == 'I' || cigar_chars[j] == 'S')
						{
							sequence_offset++;
						}
						else if (cigar_chars[j] == 'D')
						{
							reference_offset++;
						}
					}
				}

				int nm;
				for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
				{
					nm = 0;
					for (j = 0; j < msa_index; j++)
					{
						if (msa_positions[j] < msa_sequence_length)
						{
							char msa_seq_base = msa_sequences[msa_seq_idx][msa_positions[j]];
							if (msa_seq_base != bases[j] && msa_seq_base != '-' && msa_seq_base != '\0')
							{
								nm++;
							}
						}
					}

					if (nm < current_mismatch_matrix_row[msa_seq_idx])
					{
						current_mismatch_matrix_row[msa_seq_idx] = nm;
					}
				}

				best_alignment_size = current_alignment_size;
			}
		}

		if (best_alignment_size != -1)
		{
			sprintf(row_buffer, "%s\t%d", readname, best_alignment_size);
			for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
			{
				char num_buffer[16];
				sprintf(num_buffer, "\t%d", current_mismatch_matrix_row[msa_seq_idx]);
				strcat(row_buffer, num_buffer);
			}

			pthread_mutex_lock(thread_str->write_mutex);
			fprintf(thread_str->outfile, "%s\n", row_buffer);
			pthread_mutex_unlock(thread_str->write_mutex);
		}
	}

	free(row_buffer);
	free(readname);
	free(copy);
	free(sam_line_cigar);
	free(sequence);
	free(cigar_vals);
	free(cigar_chars);
	free(msa_positions);
	free(bases);
	free(current_mismatch_matrix_row);
	return NULL;
}

void write_mismatch_matrix(char *outfile_path, ReferenceData *reference_data_strs, int num_references, MSA *msa_str, int paired, int num_threads)
{
	int i;

	FILE *outfile = fopen(outfile_path, "w");
	if (outfile == NULL)
	{
		fprintf(stderr, "Error: could not open '%s' for writing the mismatch matrix\n", outfile_path);
		exit(1);
	}

	fprintf(outfile, "qName\tblockSizes");
	for (i = 0; i < msa_str->num_sequences; i++)
	{
		fprintf(outfile, "\t%s", msa_str->sequence_names[i]);
	}
	fprintf(outfile, "\n");

	int num_sam_lines = reference_data_strs[0].sam_results_str.num_sam_lines;
	int num_read_strains;
	int lines_per_read_strain;
	void *(*function)(void *);
	if (paired)
	{
		num_read_strains = num_sam_lines / 2;
		lines_per_read_strain = 2;
		function = write_mismatch_matrix_paired;
	}
	else
	{
		num_read_strains = num_sam_lines;
		lines_per_read_strain = 1;
		function = write_mismatch_matrix_single;
	}

	pthread_mutex_t write_mutex;
	pthread_mutex_init(&write_mutex, NULL);

	pthread_t *threads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
	MismatchMatrixThreadStruct *thread_strs = (MismatchMatrixThreadStruct *)malloc(num_threads * sizeof(MismatchMatrixThreadStruct));

	int read_strains_per_thread = num_read_strains / num_threads;
	int remainder = num_read_strains % num_threads;

	int read_strain_cursor = 0;
	for (i = 0; i < num_threads; i++)
	{
		int this_thread_read_strains = read_strains_per_thread;
		if (i < remainder)
		{
			this_thread_read_strains++;
		}

		int start_read_strain = read_strain_cursor;
		int end_read_strain = read_strain_cursor + this_thread_read_strains;
		read_strain_cursor = end_read_strain;

		thread_strs[i].thread_index = i;
		thread_strs[i].sam_partition_start = start_read_strain * lines_per_read_strain;
		thread_strs[i].sam_partition_end = end_read_strain * lines_per_read_strain;
		thread_strs[i].num_references = num_references;
		thread_strs[i].reference_data_strs = reference_data_strs;
		thread_strs[i].msa_str = msa_str;
		thread_strs[i].outfile = outfile;
		thread_strs[i].write_mutex = &write_mutex;

		pthread_create(&threads[i], NULL, function, &thread_strs[i]);
	}

	for (i = 0; i < num_threads; i++)
	{
		pthread_join(threads[i], NULL);
	}

	pthread_mutex_destroy(&write_mutex);
	fclose(outfile);
	free(threads);
	free(thread_strs);
}