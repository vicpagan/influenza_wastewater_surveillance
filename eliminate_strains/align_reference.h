#ifndef ALIGN_REFERENCE_H
#define ALIGN_REFERENCE_H

#include "global.h"

/**
 * @brief 
 * 
 * @param num_references 
 * @param msa_reference_filepaths 
 * @param bowtie2_reference_filepaths 
 * @return ReferencesData 
 */
ReferencesData align_references(char **reference_sequences_filepaths, MSA *msa_str, int num_references);

#endif // ALIGN_REFERENCE_H