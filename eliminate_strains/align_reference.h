#ifndef ALIGN_REFERENCE_H
#define ALIGN_REFERENCE_H

#include "global.h"

/**
 * @brief 
 * 
 * @param reference_data_str 
 * @param msa_reference_filepath 
 * @param bowtie_reference_filepath 
 */
void align_reference(ReferenceData *reference_data_str, char *msa_reference_filepath, char *bowtie_reference_filepath);

/**
 * @brief 
 * 
 * @param reference_data_strs 
 * @param num_references 
 * @param msa_reference_filepaths 
 * @param bowtie2_reference_filepaths 
 */
void align_all_references(ReferenceData *reference_data_strs, int num_references, char **msa_reference_filepaths, char **bowtie2_reference_filepaths);

#endif // ALIGN_REFERENCE_H