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
 * @param sam_results_filepaths 
 * @param num_references 
 * @return SAMResults 
 */
SAMResults read_in_sam_results(char **sam_results_filepaths, int num_references);


#endif // SAM_H