#ifndef BOWTIE_ALIGNMENT_H
#define BOWTIE_ALIGNMENT_H

/**
 * @brief Runs Bowtie2 on the reads against one index, writing a SAM file. Stops the program if
 * Bowtie2 does not finish successfully.
 *
 * @param report_all 1: every alignment of a read is reported (--all); 0: only its best one
 * @param bowtie2_args extra Bowtie2 options (-B/--bowtie2-args), placed after the program's own
 * options so that they take precedence; may be empty
 */
void run_bowtie2(char *index_prefix, char *sam_filepath, int report_all, char *bowtie2_args, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, int using_paired_end_reads, int using_fasta_format, int num_threads, int verbose);

void perform_bowtie_alignments(char *reference_sequences_dir, char* bowtie2_indexes_dir, char **reference_strain_names, int num_references, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format, int num_threads, int verbose, char *bowtie2_args);


#endif // BOWTIE_ALIGNMENT_H
