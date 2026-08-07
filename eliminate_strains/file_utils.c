#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

#include "file_utils.h"
#include "global.h"

/**
 * @brief 
 * 
 * @param fasta_path 
 * @return char* 
 */
char *read_fasta_header_name(char *fasta_filepath)
{
	FILE *file = fopen(fasta_filepath, "r");
	if (file == NULL)
	{
		fprintf(stderr, "Error: could not open reference file '%s' to read its name\n", fasta_filepath);
		exit(1);
	}
	char buffer[FASTA_MAXLINE];
	if (fgets(buffer, FASTA_MAXLINE, file) == NULL || buffer[0] != '>')
	{
		fprintf(stderr, "Error: '%s' does not start with a FASTA header line\n", fasta_filepath);
		exit(1);
	}
	fclose(file);
	buffer[strcspn(buffer, "\r\n")] = '\0';
	return strdup(buffer + 1);
}

/**
 * @brief 
 * 
 * @param a 
 * @param b 
 * @return int 
 */
static int cmp_str(const void *a, const void *b)
{
	return strcmp(*(const char **)a, *(const char **)b);
}

/**
 * @brief 
 * 
 * @param dir_path 
 * @param num_references 
 * @param dir_label 
 * @return char** 
 */
char **list_sorted_dir_files(char *dir_path, int num_references, char *dir_label)
{
	DIR *dir = opendir(dir_path);
	if (dir == NULL)
	{
		fprintf(stderr, "Error: could not open %s directory '%s'\n", dir_label, dir_path);
		exit(1);
	}
 
	char **filenames = NULL;
	int count = 0;
	struct dirent *entry;
	while ((entry = readdir(dir)) != NULL)
	{
		if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
		{
			continue;
		}
		filenames = (char **)realloc(filenames, (count + 1) * sizeof(char *));
		filenames[count] = strdup(entry->d_name);
		count++;
	}
	closedir(dir);
 
	if (count < num_references)
	{
		fprintf(stderr, "Error: %s directory '%s' has %d file(s), but %d reference(s) were requested (-N).\n", dir_label, dir_path, count, num_references);
		exit(1);
	}
 
	// sort alphabetically so position i means the same subtype across all
	// three directories (msa/reference/bowtie2-reference)
	qsort(filenames, count, sizeof(char *), cmp_str);
 
	if (count > num_references)
	{
		fprintf(stderr, "Warning: %s directory '%s' has %d files; only using the first %d (sorted alphabetically).\n", dir_label, dir_path, count, num_references);
	}
 
	char **paths = (char **)malloc(num_references * sizeof(char *));
	int i;
	for (i = 0; i < num_references; i++)
	{
		paths[i] = (char *)malloc((strlen(dir_path) + strlen(filenames[i]) + 2) * sizeof(char));
		sprintf(paths[i], "%s/%s", dir_path, filenames[i]);
	}
 
	for (i = 0; i < count; i++)
	{
		free(filenames[i]);
	}
	free(filenames);
 
	return paths;
}
