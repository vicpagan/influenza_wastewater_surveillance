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

/* ======================= paired-end mismatch counting =======================
 * Counts mismatches exactly like the original implementation's writeMismatchMatrix_paired()
 * (github.com/lpipes/SARS_CoV_2_wastewater_surveillance), including its quirks:
 *   - only A/C/G/T read bases are compared (N and other symbols are ignored);
 *   - overlapping mates: positions covered by M in both mates are recorded as REFERENCE
 *     coordinates; where the mates' bases differ there, mate 1's mismatch is taken back;
 *     mate 2 then skips the positions whose MSA COLUMN is among those reference coordinates.
 *     So overlapping mismatches are counted once per mate except where a reference coordinate
 *     happens to equal its MSA column, and a mate-2 position outside the overlap is skipped if
 *     its MSA column equals an overlap coordinate;
 *   - a pair with one mate unmapped is counted from the mapped mate alone;
 *   - aligned size = M bases of both mates minus the overlapping positions.
 *
 * For each read, the counters below record how often each of these rules changed the
 * mismatch counts of the row that was kept (the best reference). "Mismatches" are summed over
 * all strains in the MSA: one read position that mismatches 500 strains counts 500. */

enum
{
	STAT_SINGLE_MATE_READS,     // reads counted from one mate (the other was unmapped)
	STAT_SINGLE_MATE_MM,        //   mismatches those reads contributed
	STAT_OVERLAP_TWICE_POS,     // overlapping positions counted for both mates (bases agree)
	STAT_OVERLAP_TWICE_MM,      //   extra mismatches from counting them twice
	STAT_DISAGREE_POS,          // overlapping positions where the mates' bases differ
	STAT_DISAGREE_TAKEBACK_MM,  //   mismatches taken back from mate 1 there
	STAT_DISAGREE_MATE2_MM,     //   mismatches counted from mate 2 there
	STAT_WRONG_SKIP_POS,        // mate-2 positions outside the overlap skipped anyway
	STAT_WRONG_SKIP_MM,         //   mismatches not counted because of that
	STAT_NON_ACGT_POS,          // non-A/C/G/T read bases ignored
	STAT_NON_ACGT_MM,           //   mismatches they would have added if compared
	STAT_TOTAL_MM               // all mismatches in the kept rows (for scale)
};
#if NUM_COUNTING_STATS != 12
#error "NUM_COUNTING_STATS in global.h must match the counters in build_mismatch_matrix.c"
#endif

/* what an entry does to a row, and which counter its mismatches go to (-1: none) */
enum { ENTRY_COUNT, ENTRY_SINGLE, ENTRY_TWICE, ENTRY_TAKEBACK, ENTRY_DISAGREE_MATE2, ENTRY_WRONG_SKIP, ENTRY_NON_ACGT };
static const int entry_weight[] = {1, 1, 1, -1, 1, 0, 0};
static const int entry_stat[] = {-1, STAT_SINGLE_MATE_MM, STAT_OVERLAP_TWICE_MM, STAT_DISAGREE_TAKEBACK_MM,
								 STAT_DISAGREE_MATE2_MM, STAT_WRONG_SKIP_MM, STAT_NON_ACGT_MM};

typedef struct PairWork
{
	int refs1[MAX_READ_LENGTH], cols1[MAX_READ_LENGTH];
	int refs2[MAX_READ_LENGTH], cols2[MAX_READ_LENGTH];
	char bases1[MAX_READ_LENGTH], bases2[MAX_READ_LENGTH];
	int in_overlap2[MAX_READ_LENGTH]; // for mate-2 positions: 0 = not overlapping, 1 = bases agree, 2 = differ
	int visited[MAX_READ_LENGTH];
	int entry_cols[3 * MAX_READ_LENGTH];
	char entry_bases[3 * MAX_READ_LENGTH];
	int entry_types[3 * MAX_READ_LENGTH];
} PairWork;

static int is_acgt(char base)
{
	return base == 'A' || base == 'C' || base == 'G' || base == 'T' ||
		   base == 'a' || base == 'c' || base == 'g' || base == 't';
}

static int compare_ints(const void *a, const void *b)
{
	int x = *(const int *)a, y = *(const int *)b;
	return (x > y) - (x < y);
}

/* One mate's M positions: reference coordinate, MSA column (-1 if none) and read base. */
static int parse_mate(char **fields, const int *ref_to_msa, int ref_len, int *refs, int *cols, char *bases)
{
	int start = atoi(fields[3]) - 1;
	const char *cigar = fields[5];
	const char *sequence = fields[9];
	int n = 0, count = 0, seq_off = 0, ref_off = 0;
	for (const char *c = cigar; *c != '\0'; c++)
	{
		if (isdigit((unsigned char)*c))
		{
			count = count * 10 + (*c - '0');
			continue;
		}
		for (int k = 0; k < count; k++)
		{
			if (*c == 'M')
			{
				if (n >= MAX_READ_LENGTH)
				{
					fprintf(stderr, "Error: read '%s' has more than %d aligned bases (MAX_READ_LENGTH).\n", fields[0], MAX_READ_LENGTH);
					exit(1);
				}
				int r = start + ref_off;
				refs[n] = r;
				cols[n] = (r >= 0 && r < ref_len) ? ref_to_msa[r] : -1;
				bases[n] = sequence[seq_off];
				n++;
				seq_off++;
				ref_off++;
			}
			else if (*c == 'I' || *c == 'S')
			{
				seq_off++;
			}
			else if (*c == 'D')
			{
				ref_off++;
			}
		}
		count = 0;
	}
	return n;
}

/* One reference's mismatch row for a read pair, and the counters for it.
 * Returns the aligned size, or -1 if neither mate is mapped to this reference. */
static int count_pair_mismatches(char **first_fields, char **second_fields, const int *ref_to_msa, int ref_len,
								 char **msa_sequences, int msa_sequence_length, int num_msa_sequences,
								 int *row, long long *stats, PairWork *w)
{
	int first_status = parse_sam_flags(atoi(first_fields[1]));
	int second_status = parse_sam_flags(atoi(second_fields[1]));
	int n_entries = 0, alignment_size;
	memset(stats, 0, NUM_COUNTING_STATS * sizeof(long long));

#define ADD_ENTRY(col, base, type)                                   \
	do                                                               \
	{                                                                \
		w->entry_cols[n_entries] = (col);                            \
		w->entry_bases[n_entries] = (char)toupper((unsigned char)(base)); \
		w->entry_types[n_entries++] = (type);                        \
	} while (0)
#define VALID_COL(col) ((col) >= 0 && (col) < msa_sequence_length)

	if (first_status == 0 || first_status == 1)
	{
		/* both mates mapped */
		int n1 = parse_mate(first_fields, ref_to_msa, ref_len, w->refs1, w->cols1, w->bases1);
		int n2 = parse_mate(second_fields, ref_to_msa, ref_len, w->refs2, w->cols2, w->bases2);

		/* mate 1: every A/C/G/T base with an MSA column */
		for (int i = 0; i < n1; i++)
		{
			if (!VALID_COL(w->cols1[i]))
			{
				continue;
			}
			if (is_acgt(w->bases1[i]))
			{
				ADD_ENTRY(w->cols1[i], w->bases1[i], ENTRY_COUNT);
			}
			else
			{
				stats[STAT_NON_ACGT_POS]++;
				ADD_ENTRY(w->cols1[i], w->bases1[i], ENTRY_NON_ACGT);
			}
		}

		/* overlap: same reference coordinate in both mates */
		int n_visited = 0, i = 0, k = 0;
		memset(w->in_overlap2, 0, n2 * sizeof(int));
		while (i < n1 && k < n2)
		{
			if (w->refs1[i] < w->refs2[k])
			{
				i++;
			}
			else if (w->refs1[i] > w->refs2[k])
			{
				k++;
			}
			else
			{
				w->visited[n_visited++] = w->refs1[i]; /* a REFERENCE coordinate */
				if (w->bases1[i] != w->bases2[k])
				{
					w->in_overlap2[k] = 2;
					stats[STAT_DISAGREE_POS]++;
					if (is_acgt(w->bases1[i]) && VALID_COL(w->cols1[i]))
					{
						/* mates disagree: take back mate 1's mismatch here */
						ADD_ENTRY(w->cols1[i], w->bases1[i], ENTRY_TAKEBACK);
					}
				}
				else
				{
					w->in_overlap2[k] = 1;
				}
				i++;
				k++;
			}
		}
		qsort(w->visited, n_visited, sizeof(int), compare_ints);

		/* mate 2: skip positions whose MSA COLUMN is among the visited REFERENCE coordinates */
		for (k = 0; k < n2; k++)
		{
			int col = w->cols2[k];
			if (!VALID_COL(col))
			{
				continue;
			}
			if (!is_acgt(w->bases2[k]))
			{
				stats[STAT_NON_ACGT_POS]++;
				ADD_ENTRY(col, w->bases2[k], ENTRY_NON_ACGT);
				continue;
			}
			int skip = bsearch(&col, w->visited, n_visited, sizeof(int), compare_ints) != NULL;
			if (skip)
			{
				if (w->in_overlap2[k] == 0)
				{
					/* skipped although mate 1 doesn't cover this position */
					stats[STAT_WRONG_SKIP_POS]++;
					ADD_ENTRY(col, w->bases2[k], ENTRY_WRONG_SKIP);
				}
				/* else: an overlapping position, counted once as intended */
			}
			else if (w->in_overlap2[k] == 1)
			{
				stats[STAT_OVERLAP_TWICE_POS]++;
				ADD_ENTRY(col, w->bases2[k], ENTRY_TWICE);
			}
			else if (w->in_overlap2[k] == 2)
			{
				ADD_ENTRY(col, w->bases2[k], ENTRY_DISAGREE_MATE2);
			}
			else
			{
				ADD_ENTRY(col, w->bases2[k], ENTRY_COUNT);
			}
		}
		alignment_size = n1 + n2 - n_visited;
	}
	else if (first_status == 2 || (first_status == -1 && second_status == 2))
	{
		/* only one mate mapped: count it alone */
		char **fields = (first_status == 2) ? first_fields : second_fields;
		int n = parse_mate(fields, ref_to_msa, ref_len, w->refs1, w->cols1, w->bases1);
		for (int i = 0; i < n; i++)
		{
			if (!VALID_COL(w->cols1[i]))
			{
				continue;
			}
			if (is_acgt(w->bases1[i]))
			{
				ADD_ENTRY(w->cols1[i], w->bases1[i], ENTRY_SINGLE);
			}
			else
			{
				stats[STAT_NON_ACGT_POS]++;
				ADD_ENTRY(w->cols1[i], w->bases1[i], ENTRY_NON_ACGT);
			}
		}
		stats[STAT_SINGLE_MATE_READS] = 1;
		alignment_size = n;
	}
	else
	{
		return -1;
	}
#undef ADD_ENTRY
#undef VALID_COL

	for (int s = 0; s < num_msa_sequences; s++)
	{
		int nm = 0;
		for (int e = 0; e < n_entries; e++)
		{
			char msa_base = msa_sequences[s][w->entry_cols[e]];
			if (msa_base != w->entry_bases[e] && msa_base != '-' && msa_base != '\0')
			{
				int type = w->entry_types[e];
				nm += entry_weight[type];
				if (entry_stat[type] >= 0)
				{
					stats[entry_stat[type]]++;
				}
			}
		}
		row[s] = nm;
		stats[STAT_TOTAL_MM] += nm;
	}
	return alignment_size;
}

void *build_mismatch_matrix_paired(void *ptr)
{
	int j, ref_idx, sam_line_idx;

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
	char *first_copy = (char *)malloc((max_sam_line_length + 1) * sizeof(char));
	char *second_copy = (char *)malloc((max_sam_line_length + 1) * sizeof(char));
	char *first_sam_fields[11];
	char *second_sam_fields[11];

	int *current_mismatch_matrix_row = (int *)malloc(num_msa_sequences * sizeof(int));
	int *best_mismatch_matrix_row = (int *)malloc(num_msa_sequences * sizeof(int));
	PairWork *work = (PairWork *)malloc(sizeof(PairWork));
	long long current_stats[NUM_COUNTING_STATS] = {0}, best_stats[NUM_COUNTING_STATS] = {0};
	memset(bmm_thread_str->counting_stats, 0, sizeof(bmm_thread_str->counting_stats));

	for (sam_line_idx = sam_partition_start; sam_line_idx < sam_partition_end; sam_line_idx = sam_line_idx + 2)
	{
		int row_idx = bmm_thread_str->read_strain_offset + ((sam_line_idx - sam_partition_start) / 2);
		int best_reference_mismatch = INT32_MAX;
		int best_alignment_size = -1;

		for (ref_idx = 0; ref_idx < num_references; ref_idx++)
		{
			strcpy(first_copy, sam_results_str.sam_results[ref_idx][sam_line_idx]);
			strcpy(second_copy, sam_results_str.sam_results[ref_idx][sam_line_idx + 1]);

			char *first_local_save = NULL;
			char *first_token = strtok_r(first_copy, "\t", &first_local_save);
			j = 0;
			while (first_token != NULL && j < 11)
			{
				first_sam_fields[j] = first_token;
				j++;
				first_token = strtok_r(NULL, "\t", &first_local_save);
			}

			char *second_local_save = NULL;
			char *second_token = strtok_r(second_copy, "\t", &second_local_save);
			j = 0;
			while (second_token != NULL && j < 11)
			{
				second_sam_fields[j] = second_token;
				j++;
				second_token = strtok_r(NULL, "\t", &second_local_save);
			}

			if (ref_idx == 0)
			{
				strcpy(readname, first_sam_fields[0]);
			}

			int current_alignment_size = count_pair_mismatches(first_sam_fields, second_sam_fields,
				references_data_str->references_to_msa_positions[ref_idx], references_data_str->reference_sequence_lengths[ref_idx],
				msa_sequences, msa_sequence_length, num_msa_sequences, current_mismatch_matrix_row, current_stats, work);
			if (current_alignment_size < 0)
			{
				continue;
			}

			/* keep the reference whose own strain the read matches best (ties: earlier reference) */
			int current_reference_mismatch = current_mismatch_matrix_row[references_data_str->reference_sequence_msa_indexes[ref_idx]];
			if (current_reference_mismatch < best_reference_mismatch)
			{
				best_reference_mismatch = current_reference_mismatch;
				best_alignment_size = current_alignment_size;
				memcpy(best_mismatch_matrix_row, current_mismatch_matrix_row, num_msa_sequences * sizeof(int));
				memcpy(best_stats, current_stats, sizeof(best_stats));
			}
		}

		strcpy(mismatch_data_str->read_names[row_idx], readname);
		mismatch_data_str->alignment_sizes[row_idx] = best_alignment_size;
		if (best_alignment_size != -1)
		{
			memcpy(mismatch_data_str->mismatch_matrix[row_idx], best_mismatch_matrix_row, num_msa_sequences * sizeof(int));
			for (int st = 0; st < NUM_COUNTING_STATS; st++)
			{
				bmm_thread_str->counting_stats[st] += best_stats[st];
			}
		}
	}

	free(readname);
	free(first_copy);
	free(second_copy);
	free(current_mismatch_matrix_row);
	free(best_mismatch_matrix_row);
	free(work);

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

			char *local_save = NULL;
			char *token = strtok_r(copy, "\t", &local_save);

			j = 0;
			while (token != NULL && j < 11)
			{
				sam_fields[j] = token;
				j++;
				token = strtok_r(NULL, "\t", &local_save);
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
						int pos_in_msa = references_data_str->references_to_msa_positions[ref_idx][sequence_start_pos + reference_offset];

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
		mismatch_data_str->alignment_sizes[row_idx] = best_alignment_size;
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

	MismatchData mismatch_data;
	mismatch_data.num_reads = num_reads;
	mismatch_data.num_msa_sequences = msa_str->num_sequences;
	
	mismatch_data.msa_sequence_names = (char **)malloc(msa_str->num_sequences * sizeof(char *));
	for (msa_seq_idx = 0; msa_seq_idx < msa_str->num_sequences; msa_seq_idx++)
	{
		mismatch_data.msa_sequence_names[msa_seq_idx] = strdup(msa_str->sequence_names[msa_seq_idx]);
	}

	mismatch_data.alignment_sizes = (int *)malloc(num_reads * sizeof(int));

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

	if (paired)
	{
		long long totals[NUM_COUNTING_STATS] = {0};
		for (i = 0; i < num_threads; i++)
		{
			for (int st = 0; st < NUM_COUNTING_STATS; st++)
			{
				totals[st] += bmm_thread_strs[i].counting_stats[st];
			}
		}
		printf("Mismatch counting (the original implementation's rules), in the kept rows;\n");
		printf("mismatches are summed over all %d strains (%lld in the whole matrix):\n", msa_str->num_sequences, totals[STAT_TOTAL_MM]);
		printf("  reads counted from one mate (other unmapped):  %lld reads, %lld mismatches\n",
			   totals[STAT_SINGLE_MATE_READS], totals[STAT_SINGLE_MATE_MM]);
		printf("  overlapping positions counted for both mates:  %lld positions, %lld extra mismatches\n",
			   totals[STAT_OVERLAP_TWICE_POS], totals[STAT_OVERLAP_TWICE_MM]);
		printf("  overlapping positions where the mates differ:  %lld positions, %lld mismatches taken back from mate 1, %lld counted from mate 2\n",
			   totals[STAT_DISAGREE_POS], totals[STAT_DISAGREE_TAKEBACK_MM], totals[STAT_DISAGREE_MATE2_MM]);
		printf("  mate-2 positions outside the overlap skipped:  %lld positions, %lld mismatches not counted\n",
			   totals[STAT_WRONG_SKIP_POS], totals[STAT_WRONG_SKIP_MM]);
		printf("  non-A/C/G/T read bases ignored:                %lld bases, %lld mismatches not counted\n\n",
			   totals[STAT_NON_ACGT_POS], totals[STAT_NON_ACGT_MM]);
	}

	free(pthreads);
	free(bmm_thread_strs);

	// remove reads that did not align to any reference
	int write_idx = 0;
	for (int read_idx = 0; read_idx < num_reads; read_idx++)
	{
		if (mismatch_data.alignment_sizes[read_idx] != -1)
		{
			if (write_idx != read_idx)
			{
				mismatch_data.read_names[write_idx] = mismatch_data.read_names[read_idx];
				mismatch_data.alignment_sizes[write_idx] = mismatch_data.alignment_sizes[read_idx];
				mismatch_data.mismatch_matrix[write_idx] = mismatch_data.mismatch_matrix[read_idx];
			}
			write_idx++;
		}
		else
		{
			free(mismatch_data.read_names[read_idx]);
			free(mismatch_data.mismatch_matrix[read_idx]);
		}
	}
	mismatch_data.num_reads = write_idx;
	mismatch_data.read_names = (char **)realloc(mismatch_data.read_names, write_idx * sizeof(char *));
	mismatch_data.alignment_sizes = (int *)realloc(mismatch_data.alignment_sizes, write_idx * sizeof(int));
	mismatch_data.mismatch_matrix = (int **)realloc(mismatch_data.mismatch_matrix, write_idx * sizeof(int *));

	return mismatch_data;
}