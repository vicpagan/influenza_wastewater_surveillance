#ifndef WRITE_MISMATCH_MATRIX_H
#define WRITE_MISMATCH_MATRIX_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "global.h"


/**
 * @brief Thread worker: build mismatch-matrix rows for this thread's slice of the
 * preloaded SAM lines (global `sam_results`).
 *
 * For each read pair, counts mismatches (per strain, in the global `resize_MSA`
 * panel) at every aligned position, avoiding double-counting bases covered by
 * both mates of a pair (see the `visited` array). Formats one output row per
 * read pair as "readname\\talignment_size\\tmismatch_count...", written into
 * this thread's ThreadStruct.results_str.
 *
 * @param ptr Pointer to this thread's ThreadStruct (cast internally).
 * @return NULL always; real output is written into ptr->results_str, not returned.
 */
void *writeMismatchMatrix_paired(void *ptr);

/**
 * @brief Single-threaded, streaming counterpart to writeMismatchMatrix_paired().
 *
 * Same per-pair mismatch-counting logic, but reads SAM lines directly from
 * `samfile` one at a time (instead of the preloaded `sam_results` array) and
 * writes rows straight to `outfile` -- used for -n/--no-read-sam mode, which
 * trades speed for lower memory use.
 *
 * @param outfile Output mismatch-matrix file.
 * @param samfile Open SAM file (paired-end alignments).
 * @param MSA Strain panel (post-elimination), one row per remaining strain.
 * @param length_of_MSA Number of MSA columns.
 * @param number_of_strains Total strains before elimination (unused here; kept for signature symmetry).
 * @param number_of_strains_remaining Strains left after elimination (columns in the matrix).
 * @param names_of_strains Names of the remaining strains (written as the header row).
 * @param reference_index Maps SAM alignment position -> MSA column.
 */
void writeMismatchMatrix_paired_no_read_bam(FILE *outfile, FILE *samfile, char **MSA, int length_of_MSA, int number_of_strains, int number_of_strains_remaining, char **names_of_strains, int *reference_index);

/**
 * @brief Single-end counterpart to writeMismatchMatrix_paired_no_read_bam(): builds
 * the mismatch matrix one read at a time, with no mate-pair overlap handling needed.
 * @param outfile Output mismatch-matrix file.
 * @param samfile Open SAM file (single-end alignments).
 * @param MSA Full (pre-elimination) strain panel; indexed via strains_kept.
 * @param strains_kept Row indices into MSA of the strains that survived elimination.
 * @param length_of_MSA Number of MSA columns.
 * @param number_of_strains Total strains before elimination.
 * @param number_of_strains_remaining Strains left after elimination (columns in the matrix).
 */
void writeMismatchMatrix(FILE *outfile, FILE *samfile, char **MSA, int *strains_kept, int length_of_MSA, int number_of_strains, int number_of_strains_remaining);


#endif // WRITE_MISMATCH_MATRIX_H