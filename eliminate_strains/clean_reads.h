#ifndef CLEAN_READS_H
#define CLEAN_READS_H

/**
 * @brief 
 * 
 * @param filename 
 * @param sequence_length_threshold 
 * @param trim_length 
 */
void trim_ends_and_filter_fastq(const char *filename, int sequence_length_threshold, int trim_length);

/**
 * @brief 
 * 
 * @param filename 
 * @param sequence_length_threshold 
 * @param trim_length 
 */
void trim_ends_and_filter_fasta(const char *filename, int sequence_length_threshold, int trim_length);

/**
 * @brief 
 * 
 * @param filepath 
 * @return char* 
 */
char *get_fastx_prefix(const char *filepath);

/**
 * @brief 
 * 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param working_dir 
 * @param using_paired_end_reads 
 * @param using_fasta_format 
 * @param sequence_length_threshold 
 * @param trim_length 
 * @param fastq_trimmer_threshold 
 */
void clean_reads(char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format, int sequence_length_threshold, int trim_length, int fastq_trimmer_threshold);

#endif // CLEAN_READS_H