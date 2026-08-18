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
ReferencesData align_references(char **msa_reference_filepaths, char **bowtie2_reference_filepaths, int num_references);

#endif // ALIGN_REFERENCE_H