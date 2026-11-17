#ifndef SELECT_REFERENCES_H
#define SELECT_REFERENCES_H

#include "global.h"

// minimap2 options of the first pass (all strains) unless --minimap2-args replaces them
// (-t, the index and the read files are always added by the program). They are the options
// EsViritu aligns with, without its base-level alignment:
//   -x sr              short-read preset
//   --secondary=yes    report more than the single best strain of a read
//   -N 100             at most 100 secondary alignments per read: its best-scoring strains
//   -p 0.9             only strains scoring close to the read's best one
//   -f 10000           skip seeds that occur more than 10,000 times in the database
#define DEFAULT_MINIMAP2_ARGS "-x sr --secondary=yes -N 100 -p 0.9 -f 10000"

// added after the first-pass options for the second pass (provisional strains only), where a
// later option replaces an earlier one: no limit on the alignments per read (-N), no seeds
// skipped (-f), and small batches (-K), since minimap2 holds the output of a whole batch in memory
#define SECOND_PASS_MINIMAP2_ARGS "-N 1000000 -f 1000000 -K 500k"

/**
 * @brief Chooses the reference strains for a sample, in two passes of minimap2 over the
 * non-imputed sequences of the strains (every <STRAIN>.fasta in the reference sequences
 * directory):
 *
 * 1. All reads against all strains, keeping only each read's best-scoring strains. A greedy
 *    set cover of these lists (repeatedly the strain that covers the most reads not covered
 *    yet) gives a provisional choice. Reads that aligned to fewer than
 *    opt->min_strains_per_read strains are left out.
 * 2. All reads against the provisional strains only, keeping every strain a read aligns to.
 *    The greedy set cover of these complete lists gives the final choice.
 *
 * A strain can be chosen only if it has an MSA row: its own, or that of its group of redundant
 * strains. The minimap2 index of all sequences is kept in the Bowtie2 indexes directory
 * (all_reference_sequences.mmi) and built if missing. Writes selected_references.tsv to the
 * output directory.
 *
 * @param opt program options (directories, reads, threads, ...)
 * @param msa_str the MSA (for its row names)
 * @param redundant_groups groups of redundant strains (may have no members)
 * @param selected_names output: newly allocated list of the chosen strain names
 * @return int number of chosen strains
 */
int select_reference_strains(Options *opt, MSA *msa_str, RedundantGroups *redundant_groups, char ***selected_names);

#endif // SELECT_REFERENCES_H
