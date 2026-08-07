#ifndef BOWTIE_ALIGNMENT_H
#define BOWTIE_ALIGNMENT_H

/**
 * @brief 
 * 
 * @param sam_filepath 
 * @return int 
 */
int calculate_error_rates(char *sam_filepath);

/**
 * @brief 
 * 
 * @param bowtie2_reference_path 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param sam_results_filepath 
 * @param using_paired_end_reads 
 * @param using_fasta_format 
 */
void perform_bowtie_alignment(char *bowtie2_reference_path, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *sam_results_filepath, int using_paired_end_reads, int using_fasta_format);

/**
 * @brief 
 * 
 * @param bowtie2_reference_path 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param sam_results_filepath 
 * @param using_paired_end_reads 
 * @param using_fasta_format 
 */
void perform_bowtie_alignment_xeq(char *bowtie2_reference_path, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *sam_results_filepath, int using_paired_end_reads, int using_fasta_format);

#endif // BOWTIE_ALIGNMENT_H