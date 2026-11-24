#ifndef ANCHOR_ALIGNMENT_H
#define ANCHOR_ALIGNMENT_H

#include "global.h"

// Anchor modes: the reads are aligned once, with Bowtie2, to one index of all strains, and only
// the best alignment of each read is kept (no --all). That alignment only says which MSA
// columns the read covers; the mismatches are then counted against every MSA row on those
// columns, as in the other modes. No reference strains are chosen.
//
//   ALIGNMENT_MODE_ANCHOR_IMPUTED  the index holds every MSA row without its gaps, so imputed
//                                  bases can place a read. Position i of a row's sequence is
//                                  the row's i-th column that is not a gap.
//   ALIGNMENT_MODE_ANCHOR_REAL     the index holds the non-imputed sequence of every strain in
//                                  the reference sequences directory. Position i of a strain's
//                                  sequence is entry i of its line in the positions file.
//
// The index is kept in the Bowtie2 indexes directory under the names below and is built if
// missing. It belongs to one MSA (or one reference sequences directory): a small .info file
// next to it records which, and a run with a different one stops with an error.
#define ANCHOR_IMPUTED_INDEX_NAME "all_imputed_strains"
#define ANCHOR_REAL_INDEX_NAME "all_real_strains"

typedef struct AnchorTargets AnchorTargets;

/**
 * @brief Name of the index of the anchor mode, also the name of its SAM file (<working_dir>/<name>.sam)
 */
const char *anchor_index_name(int alignment_mode);

/**
 * @brief Builds the index of all strains if it is missing, then aligns the reads to it,
 * writing <working_dir>/<anchor_index_name>.sam
 */
void perform_anchor_alignment(Options *opt, MSA *msa_str);

/**
 * @brief Prepares the MSA positions of the targets the reads aligned to. In the real mode,
 * reads the positions lines of those strains only (the positions file is large).
 *
 * @param sam_results_str the SAM lines of the anchor alignment (one file)
 * @return AnchorTargets* to be stored in ReferencesData.anchor_targets
 */
AnchorTargets *prepare_anchor_targets(Options *opt, MSA *msa_str, SAMResults *sam_results_str);

/**
 * @brief MSA positions of a target: element i is the MSA column of base i of its sequence
 *
 * @param target_name RNAME of a SAM line
 * @param buffer room for msa_str->sequence_length positions, used in the imputed mode
 * @param length output: number of positions
 * @return const int* the positions (buffer, or memory owned by targets); NULL if the target is unknown
 */
const int *anchor_target_positions(const AnchorTargets *targets, const char *target_name, int *buffer, int *length);

void free_anchor_targets(AnchorTargets *targets);

#endif // ANCHOR_ALIGNMENT_H
