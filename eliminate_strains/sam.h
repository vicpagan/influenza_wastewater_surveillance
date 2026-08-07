#ifndef SAM_H
#define SAM_H

#include <zlib.h>
#include "global.h"

/**
 * @brief 
 * 
 * @param flag_value 
 * @return int 
 */
int parse_sam_flags(int flag_value);

/**
 * @brief 
 * 
 * @param sam_results_file 
 * @param sam_results_str 
 */
void parse_sam_info(gzFile sam_results_file, SAMResults *sam_results_str);

/**
 * @brief 
 * 
 * @param sam_results_file 
 * @param sam_results_str 
 */
void read_sam_lines(gzFile sam_results_file, SAMResults *sam_results_str);

/**
 * @brief 
 * 
 * @param sam_results_filepath 
 * @return SAMResults 
 */
SAMResults read_in_sam_results(char *sam_results_filepath);

#endif // SAM_H