#ifndef SELECT_REFERENCES_H
#define SELECT_REFERENCES_H

#include "global.h"

// minimap2 options for aligning all reads to all strains, unless --minimap2-args replaces
// them (-t, the index and the read files are always added by the program). They are the
// options EsViritu aligns with, writing PAF with base-level alignment (-c):
//   -x sr              short-read preset
//   -c                 base-level alignment, needed for the identity of an alignment
//   --secondary=yes    report more than the single best strain of a read
//   -N 100             at most 100 secondary alignments per read: its best-scoring strains
//   -p 0.9             only strains scoring close to the read's best one
//   -f 10000           skip seeds that occur more than 10,000 times in the database
#define DEFAULT_MINIMAP2_ARGS "-x sr -c --secondary=yes -N 100 -p 0.9 -f 10000"

// added after those options when aligning to a few strains only (the given reference strains,
// the provisional strains), where a later option replaces an earlier one: no limit on the
// alignments per read (-N), no seeds skipped (-f), and small batches (-K), since minimap2
// holds the output of a whole batch in memory
#define UNCAPPED_MINIMAP2_ARGS "-N 1000000 -f 1000000 -K 500k"

// An alignment counts only if at least this fraction of the read is aligned, with at least
// this identity (matching bases / alignment length); these are EsViritu's limits. A read with
// such an alignment is an influenza read.
#define MIN_ALIGNED_FRACTION 0.9
#define MIN_INFLUENZA_IDENTITY 0.8

/**
 * @brief Keeps the influenza reads and completes the list of reference strains, using
 * minimap2 on the non-imputed sequences of the strains (every <STRAIN>.fasta in the
 * reference sequences directory):
 *
 * 1. All reads are aligned to all strains (each read's best-scoring strains are reported) and,
 *    without limits, to the reference strains given with -r. A read is an influenza read if it
 *    aligns to at least opt->min_strains_per_read strains or to a given reference strain. The
 *    influenza reads are written to the working directory, and opt's read file paths are
 *    changed to those files, so that the rest of the program uses only them.
 * 2. A strain covers a read if the alignment has at least opt->min_reference_identity. For the
 *    influenza reads that no given reference strain covers, strains are added by greedy set
 *    cover (repeatedly the strain that covers the most reads not covered yet): first a
 *    provisional choice from each read's best strains, then, after aligning those reads to the
 *    provisional strains without limits, the final choice. A strain is added only if it covers
 *    at least opt->min_reads_per_added_reference reads not covered yet.
 *
 * A strain can be a reference only if it has an MSA row: its own, or that of its group of
 * redundant strains. The minimap2 index of all sequences is kept in the Bowtie2 indexes
 * directory (all_reference_sequences.mmi) and built if missing. Writes selected_references.tsv
 * to the output directory.
 *
 * @param opt program options; its reference strain names (-r) are freed, its read file paths changed
 * @param msa_str the MSA (for its row names)
 * @param redundant_groups groups of redundant strains (may have no members)
 * @param selected_names output: newly allocated list of the reference strain names, given ones first
 * @return int number of reference strains
 */
int select_reference_strains(Options *opt, MSA *msa_str, RedundantGroups *redundant_groups, char ***selected_names);

#endif // SELECT_REFERENCES_H
