#ifndef BOWTIE_ALIGNMENT_H
#define BOWTIE_ALIGNMENT_H

void perform_bowtie_alignments(char *reference_sequences_dir, char* bowtie2_indexes_dir, char **reference_strain_names, int num_references, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format, int verbose);


#endif // BOWTIE_ALIGNMENT_H