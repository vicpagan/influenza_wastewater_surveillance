#ifndef CALCULATE_ALLELE_FREQ_H
#define CALCULATE_ALLELE_FREQ_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>

#include "global.h"

/**
 * @brief 
 * 
 * @param allele 
 * @param msa_str 
 * @param freq_threshold 
 * @param tstart 
 * @param tend 
 * @param coverage 
 * @param min_strains_remaining 
 * @param max_strains_remaining 
 * @param output_allele_counts 
 * @param allele_counts_output_filepath 
 * @param output_deletions 
 * @param deletions_output_filepath 
 * @param deletion_threshold 
 * @param reference_data_str 
 * @return int 
 */
int calculate_allele_freq_paired(double **allele, MSA *msa_str, double freq_threshold, struct timespec tstart, struct timespec tend, int coverage, int min_strains_remaining, int max_strains_remaining, int output_allele_counts, char *allele_counts_output_filepath, int output_deletions, char *deletions_output_filepath, double deletion_threshold, ReferenceData *reference_data_str);


/**
 * @brief 
 * 
 * @param allele 
 * @param msa_str 
 * @param freq_threshold 
 * @param tstart 
 * @param tend 
 * @param coverage 
 * @param min_strains_remaining 
 * @param max_strains_remaining 
 * @param output_allele_counts 
 * @param allele_counts_output_filepath 
 * @param output_deletions 
 * @param deletions_output_filepath 
 * @param deletion_threshold 
 * @param reference_data_str 
 * @return int 
 */
int calculate_allele_freq(double **allele, MSA *msa_str, double freq_threshold, struct timespec tstart, struct timespec tend, int coverage, int min_strains_remaining, int max_strains_remaining, int output_allele_counts, char *allele_counts_output_filepath, int output_deletions, char *deletions_output_filepath, double deletion_threshold, ReferenceData *reference_data_str);


#endif // CALCULATE_ALLELE_FREQ_H