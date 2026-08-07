#ifndef BOWTIE_ALIGNMENT_H
#define BOWTIE_ALIGNMENT_H

/**
 * @brief 
 * 
 * @param sam_filepath 
 * @param end_region_length 
 * @param end_region_error_mult 
 * @return int 
 */
int calculate_error_rates(char *sam_filepath, int end_region_length, double end_region_error_mult);

/**
 * @brief 
 * 
 * @param bowtie2_reference_path 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param sam_results_filepath 
 * @param working_dir 
 * @param using_paired_end_reads 
 * @param using_fasta_format 
 */
void perform_bowtie_alignment(char *bowtie2_reference_path, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *sam_results_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format, int verbose);

/**
 * @brief 
 * 
 * @param bowtie2_reference_path 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param sam_results_filepath 
 * @param working_dir 
 * @param using_paired_end_reads 
 * @param using_fasta_format 
 */
void perform_bowtie_alignment_xeq(char *bowtie2_reference_path, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *sam_results_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_forma, int verbose);

#endif // BOWTIE_ALIGNMENT_H