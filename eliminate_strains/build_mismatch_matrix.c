#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <limits.h>
#include <assert.h>
#include <pthread.h>

#include "build_mismatch_matrix.h"
#include "sam.h"
#include "msa.h"

void *build_mismatch_matrix_paired(void *ptr)
{
	int i, j, k, ref_idx, sam_line_idx, msa_seq_idx;

	BuildMismatchMatrixThread *bmm_thread_str = (BuildMismatchMatrixThread *)ptr;
	int sam_partition_start = bmm_thread_str->sam_partition_start;
	int sam_partition_end = bmm_thread_str->sam_partition_end;
	int thread_id = bmm_thread_str->thread_index;

	ReferencesData *references_data_str = bmm_thread_str->references_data_str;
	int num_references = references_data_str->num_references;

	SAMResults sam_results_str = references_data_str->sam_results_str;
	int max_sam_line_length = sam_results_str.max_sam_line_length;

	MSA *msa_str = bmm_thread_str->msa_str;
	char **msa_sequences = msa_str->sequences;
	int msa_sequence_length = msa_str->sequence_length;
	int num_msa_sequences = msa_str->num_sequences;

	MismatchData *mismatch_data_str = bmm_thread_str->mismatch_data_str;

	char *readname = (char *)malloc((max_sam_line_length + 1) * sizeof(char));

	char *first_copy = (char *)malloc((max_sam_line_length + 1) * sizeof(char));
	char *second_copy = (char *)malloc((max_sam_line_length + 1) * sizeof(char));
	char *first_sam_line_cigar = (char *)malloc(MAX_CIGAR * sizeof(char));
	char *second_sam_line_cigar = (char *)malloc(MAX_CIGAR * sizeof(char));
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
	int *best_mismatch_matrix_row = (int *)malloc(num_msa_sequences * sizeof(int));

	for (sam_line_idx = sam_partition_start; sam_line_idx < sam_partition_end; sam_line_idx = sam_line_idx + 2)
	{
		int row_idx = bmm_thread_str->read_strain_offset + ((sam_line_idx - sam_partition_start) / 2);
		int best_reference_mismatch = INT32_MAX;
		int best_alignment_size = -1;

		printf("Thread %d is working on row %d\n", thread_id, row_idx);

		for (ref_idx = 0; ref_idx < num_references; ref_idx++)
		{
			printf("Thread %d is working on reference %d\n", thread_id, ref_idx);

			strcpy(first_copy, sam_results_str.sam_results[ref_idx][sam_line_idx]);
			strcpy(second_copy, sam_results_str.sam_results[ref_idx][sam_line_idx + 1]);

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
						int first_pos_in_msa = references_data_str->reference_indexes[ref_idx][first_sequence_start_pos + first_reference_offset];

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
						int second_pos_in_msa = references_data_str->reference_indexes[ref_idx][second_sequence_start_pos + second_reference_offset];

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
					current_mismatch_matrix_row[msa_seq_idx] = nm;
				}

				int current_reference_mismatch = current_mismatch_matrix_row[references_data_str->reference_sequence_msa_indexes[ref_idx]];
				if (current_reference_mismatch < best_reference_mismatch)
				{
					printf("Thread %d is replacing the best reference stuff. Current mismatch = %d      best mismatch = %d", thread_id, current_reference_mismatch, best_reference_mismatch);
					best_reference_mismatch = current_reference_mismatch;
					best_alignment_size = current_alignment_size;
					for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
					{
						best_mismatch_matrix_row[msa_seq_idx] = current_mismatch_matrix_row[msa_seq_idx];
					}
				}
			}
		}

		strcpy(mismatch_data_str->read_names[row_idx], readname);
		mismatch_data_str->block_sizes[row_idx] = best_alignment_size;
		if (best_alignment_size != -1) 
		{
			for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
			{
				mismatch_data_str->mismatch_matrix[row_idx][msa_seq_idx] = best_mismatch_matrix_row[msa_seq_idx];
			}
		}
	}

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
	free(best_mismatch_matrix_row);

	return NULL;
}

void *build_mismatch_matrix_single(void *ptr)
{
	int i, j, k, ref_idx, sam_line_idx, msa_seq_idx;

	BuildMismatchMatrixThread *bmm_thread_str = (BuildMismatchMatrixThread *)ptr;
	int sam_partition_start = bmm_thread_str->sam_partition_start;
	int sam_partition_end = bmm_thread_str->sam_partition_end;

	ReferencesData *references_data_str = bmm_thread_str->references_data_str;
	int num_references = references_data_str->num_references;

	SAMResults sam_results_str = references_data_str->sam_results_str;
	int max_sam_line_length = sam_results_str.max_sam_line_length;

	MSA *msa_str = bmm_thread_str->msa_str;
	char **msa_sequences = msa_str->sequences;
	int msa_sequence_length = msa_str->sequence_length;
	int num_msa_sequences = msa_str->num_sequences;

	MismatchData *mismatch_data_str = bmm_thread_str->mismatch_data_str;

	char *readname = (char *)malloc((max_sam_line_length + 1) * sizeof(char));

	char *copy = (char *)malloc((max_sam_line_length + 1) * sizeof(char));
	char *sam_line_cigar = (char *)malloc(MAX_CIGAR * sizeof(char));
	char *sam_fields[11];

	char *sequence = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	int *cigar_vals = (int *)malloc(MAX_CIGAR * sizeof(int));
	char *cigar_chars = (char *)malloc(MAX_CIGAR * sizeof(char));

	int *msa_positions = (int *)malloc(MAX_READ_LENGTH * sizeof(int));
	char *bases = (char *)malloc(MAX_READ_LENGTH * sizeof(char));

	int *current_mismatch_matrix_row = (int *)malloc(num_msa_sequences * sizeof(int));
	int *best_mismatch_matrix_row = (int *)malloc(num_msa_sequences * sizeof(int));

	for (sam_line_idx = sam_partition_start; sam_line_idx < sam_partition_end; sam_line_idx++)
	{
		int row_idx = bmm_thread_str->read_strain_offset + (sam_line_idx - sam_partition_start);
		int best_reference_mismatch = INT32_MAX;
		int best_alignment_size = -1;

		for (ref_idx = 0; ref_idx < num_references; ref_idx++)
		{
			strcpy(copy, sam_results_str.sam_results[ref_idx][sam_line_idx]);

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
						int pos_in_msa = references_data_str->reference_indexes[ref_idx][sequence_start_pos + reference_offset];

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
					current_mismatch_matrix_row[msa_seq_idx] = nm;
				}

				int current_reference_mismatch = current_mismatch_matrix_row[references_data_str->reference_sequence_msa_indexes[ref_idx]];
				if (current_reference_mismatch < best_reference_mismatch)
				{
					best_reference_mismatch = current_reference_mismatch;
					best_alignment_size = current_alignment_size;
					for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
					{
						best_mismatch_matrix_row[msa_seq_idx] = current_mismatch_matrix_row[msa_seq_idx];
					}
				}
			}
		}

		strcpy(mismatch_data_str->read_names[row_idx], readname);
		mismatch_data_str->block_sizes[row_idx] = best_alignment_size;
		if (best_alignment_size != -1)
		{
			for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
			{
				mismatch_data_str->mismatch_matrix[row_idx][msa_seq_idx] = best_mismatch_matrix_row[msa_seq_idx];
			}
		}
	}

	free(readname);
	free(copy);
	free(sam_line_cigar);
	free(sequence);
	free(cigar_vals);
	free(cigar_chars);
	free(msa_positions);
	free(bases);
	free(current_mismatch_matrix_row);
	free(best_mismatch_matrix_row);

	return NULL;
}

MismatchData build_mismatch_matrix(ReferencesData *references_data_str, MSA *msa_str, int paired, int num_threads)
{
	int i, ref_idx, msa_seq_idx;

	int num_references = references_data_str->num_references;
	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		references_data_str->reference_sequence_msa_indexes[ref_idx] = -1;
		msa_seq_idx = 0;
		while (msa_seq_idx < msa_str->num_sequences && references_data_str->reference_sequence_msa_indexes[ref_idx] == -1)
		{
			if (strcmp(references_data_str->reference_names[ref_idx], msa_str->sequence_names[msa_seq_idx]) == 0)
			{
				references_data_str->reference_sequence_msa_indexes[ref_idx] = msa_seq_idx;
			}
			msa_seq_idx++;
		}

		if (references_data_str->reference_sequence_msa_indexes[ref_idx] == -1)
		{
			fprintf(stderr, "Error: reference sequence '%s' was not found in the MSA.\n", references_data_str->reference_names[ref_idx]);
			exit(1);
		}
	}

	int num_sam_lines = references_data_str->sam_results_str.num_sam_lines;
	int max_sam_line_length = references_data_str->sam_results_str.max_sam_line_length;

	int num_reads;
	int lines_per_read;
	void *(*function)(void *);
	if (paired)
	{
		num_reads = num_sam_lines / 2;
		lines_per_read = 2;
		function = build_mismatch_matrix_paired;
	}
	else
	{
		num_reads = num_sam_lines;
		lines_per_read = 1;
		function = build_mismatch_matrix_single;
	}

	printf("DEBUG: num_sam_lines = %d, num_reads = %d\n", num_sam_lines, num_reads);

	MismatchData mismatch_data;
	mismatch_data.num_reads = num_reads;
	mismatch_data.num_msa_sequences = msa_str->num_sequences;
	mismatch_data.msa_sequence_names = msa_str->sequence_names;

	mismatch_data.block_sizes = (int *)malloc(num_reads * sizeof(int));

	mismatch_data.read_names = (char **)malloc(num_reads * sizeof(char *));
	mismatch_data.mismatch_matrix = (int **)malloc(num_reads * sizeof(int *));
	for (i = 0; i < num_reads; i++)
	{
		mismatch_data.read_names[i] = (char *)malloc(max_sam_line_length * sizeof(char));
		mismatch_data.mismatch_matrix[i] = (int *)malloc(msa_str->num_sequences * sizeof(int));
	}

	pthread_t *pthreads = (pthread_t *)malloc(num_threads * sizeof(pthread_t));
	BuildMismatchMatrixThread *bmm_thread_strs = (BuildMismatchMatrixThread *)malloc(num_threads * sizeof(BuildMismatchMatrixThread));

	int reads_per_thread = num_reads / num_threads;
	int remainder = num_reads % num_threads;

	int read_cursor = 0;
	for (i = 0; i < num_threads; i++)
	{
		printf("Looping for thread number %d\n", i);
		int this_thread_num_reads = reads_per_thread;
		if (i < remainder)
		{
			this_thread_num_reads++;
		}

		bmm_thread_strs[i].thread_index = i;
		bmm_thread_strs[i].sam_partition_start = read_cursor * lines_per_read;
		bmm_thread_strs[i].sam_partition_end = (read_cursor + this_thread_num_reads) * lines_per_read;
		bmm_thread_strs[i].read_strain_offset = read_cursor;
		bmm_thread_strs[i].references_data_str = references_data_str;
		bmm_thread_strs[i].msa_str = msa_str;
		bmm_thread_strs[i].mismatch_data_str = &mismatch_data;

		read_cursor += this_thread_num_reads;

		pthread_create(&pthreads[i], NULL, function, &bmm_thread_strs[i]);
	}

	for (i = 0; i < num_threads; i++)
	{
		pthread_join(pthreads[i], NULL);
	}

	free(pthreads);
	free(bmm_thread_strs);

	return mismatch_data;
}