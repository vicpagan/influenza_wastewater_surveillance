#ifndef BUILD_MISMATCH_MATRIX_H
#define BUILD_MISMATCH_MATRIX_H

#include "global.h"

/**
 * @brief 
 * 
 * @param ptr 
 * @return void* 
 */
void *build_mismatch_matrix_single(void *ptr);

/**
 * @brief 
 * 
 * @param ptr 
 * @return void* 
 */
void *build_mismatch_matrix_paired(void *ptr);

/**
 * @brief 
 * 
 */
MismatchData build_mismatch_matrix(ReferencesData *references_data_str, MSA *msa_str, int paired, int num_threads);


#endif // BUILD_MISMATCH_MATRIX_H