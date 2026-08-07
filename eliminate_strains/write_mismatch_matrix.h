#ifndef WRITE_MISMATCH_MATRIX_H
#define WRITE_MISMATCH_MATRIX_H

#include "global.h"

/**
 * @brief 
 * 
 * @param ptr 
 * @return void* 
 */
void *write_mismatch_matrix_single(void *ptr);

/**
 * @brief 
 * 
 * @param ptr 
 * @return void* 
 */
void *write_mismatch_matrix_paired(void *ptr);

/**
 * @brief 
 * 
 */
void write_mismatch_matrix(char *outfile_path, ReferenceData *reference_data_strs, int num_references, MSA *msa_str, int paired, int num_threads);


#endif // WRITE_MISMATCH_MATRIX_H