#ifndef CALCULATE_ALLELE_FREQ_H
#define CALCULATE_ALLELE_FREQ_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "global.h"

/**
 * @brief Compute per-site allele frequencies from paired-end reads, then iteratively
 * eliminate reference strains incompatible with those frequencies.
 *
 * Three stages, in order:
 *  1. Walk the SAM file read-pair by read-pair, tallying A/C/G/T counts (and deletions)
 *     at each MSA site the reads cover, translating SAM positions to MSA columns via
 *     reference_index[]. Overlapping mate pairs are only counted once (see `visited`).
 *  2. Convert counts to frequencies, drop sites below `coverage`, and mark any base
 *     with frequency < freq_threshold as "bad" at that site.
 *  3. Iteratively remove (blank the name of) any strain carrying too many "bad" bases
 *     at the variant sites, growing the incompatibility tolerance each pass until the
 *     number of strains remaining falls within [min_strains_remaining, max_strains_remaining).
 *
 * @param sam Open SAM file (paired-end alignments).
 * @param allele Output/scratch: per-site A/C/G/T counts, then converted to frequencies.
 * @param length_of_MSA Number of MSA columns.
 * @param MSA Numerically-coded... (actually char-coded) alignment; see names_of_strains.
 * @param number_of_strains Total strains before elimination.
 * @param names_of_strains Strain names; eliminated strains get their name blanked.
 * @param freq_threshold Below this frequency, a base is "bad" at a site.
 * @param maxname Length of each name buffer.
 * @param tstart Scratch timing variable (reused internally).
 * @param tend Scratch timing variable (reused internally).
 * @param number_of_variant_sites Length of variant_sites.
 * @param variant_sites Sites to check for strain elimination; entries get set to -1
 * if not adequately covered.
 * @param coverage Minimum read depth to trust a site's allele frequency.
 * @param reference_index Maps SAM alignment position -> MSA column (see align_reference.c).
 * @param min_strains_remaining Lower bound on strains remaining; loop stops once reached.
 * @param max_strains_remaining Upper bound; exceeding this aborts the whole program.
 * @param print_counts If non-empty, path to dump per-site allele counts.
 * @param max_sam_length Output: [0] = longest raw SAM line seen, [1] = number of SAM records.
 * @param print_deletions If non-empty, path to dump sites with deletion frequency above threshold.
 * @param deletion_threshold Frequency threshold for reporting a deletion site.
 * @param sam_results_out Output: newly allocated array of every raw SAM
 * line read (same format readInSamFile() would produce), so callers can reuse
 * it later without a second disk read. Pass NULL to skip caching.
 * @param num_sam_lines_out Output: number of lines in cached_sam_lines_out.
 * 
 * @return Number of strains remaining after elimination.
 */
int calculateAlleleFreq_paired(FILE *sam, double **allele, int length_of_MSA, char **MSA, int number_of_strains, char **names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int number_of_variant_sites, int *variant_sites, int coverage, int *reference_index, int min_strains_remaining, int max_strains_remaining, char print_counts[], char print_deletions[], double deletion_threshold, char ***sam_results, int *num_sam_lines, int *max_sam_line_length);


/**
 * @brief Single-end counterpart of calculateAlleleFreq_paired(): same three stages
 * (tally allele counts -> threshold "bad" bases -> iteratively eliminate strains),
 * without the mate-pair overlap handling paired-end reads need.
 * 
 * @see calculateAlleleFreq_paired for full parameter documentation (identical here),
 * including sam_results_out/num_sam_lines_out.
 */
int calculateAlleleFreq(FILE *sam, double **allele, int length_of_MSA, char **MSA, int number_of_strains, char **names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int number_of_variant_sites, int *variant_sites, int coverage, int *reference_index, int min_strains_remaining, int max_strains_remaining, char print_counts[], char print_deletions[], double deletion_threshold, char ***sam_results, int *num_sam_lines, int *max_sam_line_length);


#endif // CALCULATE_ALLELE_FREQ_H