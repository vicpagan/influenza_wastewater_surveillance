#ifndef ALIGN_REFERENCE_H
#define ALIGN_REFERENCE_H

#include "needleman_wunsch.h"
#include "global.h"

/**
 * @brief Aligns the MSA reference sequence to the Bowtie2 reference sequence and builds reference_index
 *
 * reference_index[i] = the column in the MSA reference that 
 * corresponds to position i in the Bowtie2 reference
 *
 * @param number_of_problematic_sites number of known problematic genome sites
 * @param problematic_sites list of known problematic genome sites
 * @param msa_reference_path path to the MSA reference file
 * @param bowtie_reference_path path to the Bowtie2 reference file
 * @param reference_index Output: array of indices mapping Bowtie2 reference positions to MSA reference columns
 */
void align_references(int number_of_problematic_sites, int problematic_sites[], char *msa_reference_path, char *bowtie_reference_path, int *reference_index);

#endif