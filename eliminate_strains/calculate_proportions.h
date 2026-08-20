#ifndef CALCULATE_PROPORTIONS_H
#define CALCULATE_PROPORTIONS_H

#include "global.h"

/**
 * @brief 
 * 
 * @param mismatch_data_str 
 * @param output_csv_filepath 
 * @param error_rate 
 * @param filter 
 * @param compute_strain_llr 
 * @param compute_site_llr 
 * @param num_plot 
 * @param num_threads 
 */
void calculate_proportions(MismatchData *mismatch_data_str, char *output_csv_filepath, double error_rate, double filter, int compute_strain_llr, int compute_site_llr, int num_plot, int num_threads);



#endif