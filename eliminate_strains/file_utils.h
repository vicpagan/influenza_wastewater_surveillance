#ifndef FILE_UTILS_H
#define FILE_UTILS_H

#include <stddef.h>
#include <zlib.h>

/**
 * @brief Reads one whole line of any length, without its line ending
 *
 * @param line buffer, allocated and enlarged as needed (start with *line == NULL)
 * @param capacity size of the buffer
 * @return int 1 if a line was read, 0 at the end of the file
 */
int gz_read_line(gzFile file, char **line, size_t *capacity);

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