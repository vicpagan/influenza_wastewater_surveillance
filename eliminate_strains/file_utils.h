#ifndef FILE_UTILS_H
#define FILE_UTILS_H

/**
 * @brief 
 * 
 * @param fasta_path 
 * @return char* 
 */
char *read_fastx_header_name(char *fasta_filepath);

/**
 * @brief 
 * 
 * @param dir_path 
 * @param num_references 
 * @param dir_label 
 * @return char** 
 */
char **list_sorted_dir_files(char *dir_path, int num_references, char *dir_label);

/**
 * @brief
 * 
 * @param filepath 
 * @param output_filepath 
 */
char *get_filepath_in_working_dir(char *filepath, char *working_dir);


#endif // FILE_UTILS_H