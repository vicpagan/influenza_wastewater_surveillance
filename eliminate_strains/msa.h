#ifndef MSA_H
#define MSA_H

#include <zlib.h>
#include "global.h"

/**
 * @brief 
 * 
 * @param msa_file 
 * @param msa 
 */
void parse_msa_info(gzFile msa_file, MSA *msa_str);

/**
 * @brief 
 * 
 * @param msa_file 
 * @param msa_str 
 * @param reference_sequence_name 
 */
void read_msa_sequences(gzFile msa_file, MSA *msa_str);

/**
 * @brief 
 * 
 * @param msa_filepath 
 * @param reference_sequence_name 
 * @return MSA 
 */
MSA read_in_msa(char *msa_filepath);

/**
 * @brief 
 * 
 * @param msa_str 
 * @param sequences_to_remove 
 */
void prune_msa_sequences(MSA *msa_str, int *sequences_to_remove);

/**
 * @brief 
 * 
 * @param msa_str 
 */
void remove_identical_sequences(MSA *msa_str);

#endif // MSA_H