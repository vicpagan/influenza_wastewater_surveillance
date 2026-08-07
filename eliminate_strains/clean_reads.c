#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "clean_reads.h"
#include "file_utils.h"
#include "global.h"

/**
 * @brief 
 * 
 * @param prefix 
 * @param sequence_length_threshold 
 * @param trim_length 
 */
void trim_ends_and_filter_fastq(const char *prefix, int sequence_length_threshold, int trim_length)
{
	char input_filepath[1000];
	char output_filepath[1000];

	FILE *input_file;
	FILE *output_file;
	
	strcpy(input_filepath, prefix);
	strcpy(output_filepath, prefix);
	strcat(input_filepath, "_trimmed1.fastq");
	strcat(output_filepath, "_trimmed2.fastq");
	
	if ((input_file = fopen(input_filepath, "r")) == (FILE *)NULL)
	{
		printf("Cannot open %s\n", input_filepath);
		exit(1);
	}
	if ((output_file = fopen(output_filepath, "w")) == (FILE *)NULL)
	{
		printf("Cannot open %s\n", output_filepath);
		exit(1);
	}
	
	char *strain_name = NULL;
	char *sequence = NULL;
	char *plus = NULL;
	char *quality = NULL;

	size_t strain_name_len = 0;
	size_t sequence_len = 0;
	size_t plus_len = 0;
	size_t quality_len = 0;
	
	int i;
	while (getline(&strain_name, &strain_name_len, input_file) != -1
		&& getline(&sequence, &sequence_len, input_file) != -1
		&& getline(&plus, &plus_len, input_file) != -1
		&& getline(&quality, &quality_len, input_file) != -1)
    {
		int seq_length = strlen(sequence);

		if (seq_length >= sequence_length_threshold)
		{
			fprintf(output_file, "%s", strain_name);

			for (i = trim_length; i < (seq_length - trim_length); i++)
			{
				fprintf(output_file, "%c", sequence[i]);
			}
			fprintf(output_file, "\n");

			fprintf(output_file, "%s", plus);

			for (i = trim_length; i < (seq_length - trim_length); i++)
			{
				fprintf(output_file, "%c", quality[i]);
			}
			fprintf(output_file, "\n");
		}
	}
	fclose(input_file);
	fclose(output_file);
	
	free(strain_name);
	free(sequence);
	free(plus);
	free(quality);
}

/**
 * @brief 
 * 
 * @param prefix 
 * @param sequence_length_threshold 
 * @param trim_length 
 */
void trim_ends_and_filter_fasta(const char *prefix, int sequence_length_threshold, int trim_length)
{
	char input_filepath[1000];
	char output_filepath[1000];

	FILE *input_file;
	FILE *output_file;
	
	strcpy(input_filepath, prefix);
	strcpy(output_filepath, prefix);
	strcat(input_filepath, "_trimmed1.fasta");
	strcat(output_filepath, "_trimmed2.fasta");
	
	if ((input_file = fopen(input_filepath, "r")) == (FILE *)NULL)
	{
		printf("Cannot open %s\n", input_filepath);
		exit(1);
	}
	if ((output_file = fopen(output_filepath, "w")) == (FILE *)NULL)
	{
		printf("Cannot open %s\n", output_filepath);
		exit(1);
	}
	
	char *strain_name = NULL;
	char *sequence = NULL;

	size_t strain_name_len = 0;
	size_t sequence_len = 0;

	int i;
	while (getline(&strain_name, &strain_name_len, input_file) != -1
		&& getline(&sequence, &sequence_len, input_file) != -1)
    {
		int seq_length = strlen(sequence);

		if (seq_length >= sequence_length_threshold)
		{
			fprintf(output_file, "%s", strain_name);

			for (i = trim_length; i < (seq_length - trim_length); i++)
			{
				fprintf(output_file, "%c", sequence[i]);
			}
			fprintf(output_file, "\n");
		}
	}
	fclose(input_file);
	fclose(output_file);
	
	free(strain_name);
	free(sequence);
}

char *get_fastx_prefix(const char *filepath)
{
    const char *file_ext = strrchr(filepath, '.');

    if (strcmp(file_ext, ".fastq") != 0 && strcmp(file_ext, ".fasta") != 0 && strcmp(file_ext, ".fq") != 0 && strcmp(file_ext, ".fa") != 0)
    {
        return NULL;
    }

    size_t file_ext_len = strlen(file_ext);
	size_t filepath_len = strlen(filepath);
	size_t prefix_len = filepath_len - file_ext_len;

	// allocate extra space for "_trimmed<1/2>.fast<a/q>" suffix and file extension
    char *prefix = malloc(prefix_len + 16);

    memcpy(prefix, filepath, prefix_len);
    prefix[prefix_len] = '\0';

    return prefix;
}

/**
 * @brief 
 * 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param working_dir 
 * @param using_paired 
 * @param using_fasta 
 * @param sequence_length_threshold 
 * @param trim_length 
 * @param fastq_trimmer_threshold 
 */
void clean_reads(char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format, int sequence_length_threshold, int trim_length, int fastq_trimmer_threshold)
{
	char *buffer = (char *)calloc(FASTA_MAXLINE, sizeof(char));

	if (using_paired_end_reads)
	{
		char *forward_prefix_orig = get_fastx_prefix(forward_end_filepath);
		if (forward_prefix_orig == NULL)
		{
			printf("Your reads don't end with .fastq, .fasta, .fa, or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		char *forward_prefix = get_filepath_in_working_dir(forward_prefix_orig, working_dir);
		free(forward_prefix_orig);

		if (using_fasta_format)
		{
			trim_ends_and_filter_fasta(forward_prefix, sequence_length_threshold, trim_length);
		}
		else
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t %d -i %s -o %s_trimmed1.fastq -Q33", fastq_trimmer_threshold, forward_end_filepath, forward_prefix);
			system(buffer);
			trim_ends_and_filter_fastq(forward_prefix, sequence_length_threshold, trim_length);
		}

		char forward_suffix[16];
		if (using_fasta_format)
		{
			strcpy(forward_suffix, "_trimmed2.fasta");
		}
		else
		{
			strcpy(forward_suffix, "_trimmed2.fastq");
		}

		char *forward_final = (char *)malloc(strlen(forward_prefix) + strlen(forward_suffix) + 1);
		sprintf(forward_final, "%s%s", forward_prefix, forward_suffix);
		strcpy(forward_end_filepath, forward_final);
		free(forward_final);
		free(forward_prefix);

		char *reverse_prefix_orig = get_fastx_prefix(reverse_end_filepath);
		if (reverse_prefix_orig == NULL)
		{
			printf("Your reads don't end with .fastq, .fasta, .fa, or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		char *reverse_prefix = get_filepath_in_working_dir(reverse_prefix_orig, working_dir);
		free(reverse_prefix_orig);

		if (using_fasta_format)
		{
			trim_ends_and_filter_fasta(reverse_prefix, sequence_length_threshold, trim_length);
		}
		else
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t %d -i %s -o %s_trimmed1.fastq -Q33", fastq_trimmer_threshold, reverse_end_filepath, reverse_prefix);
			system(buffer);
			trim_ends_and_filter_fastq(reverse_prefix, sequence_length_threshold, trim_length);
		}

		char reverse_suffix[16];
		if (using_fasta_format)
		{
			strcpy(reverse_suffix, "_trimmed2.fasta");
		}
		else
		{
			strcpy(reverse_suffix, "_trimmed2.fastq");
		}

		char *reverse_final = (char *)malloc(strlen(reverse_prefix) + strlen(reverse_suffix) + 1);
		sprintf(reverse_final, "%s%s", reverse_prefix, reverse_suffix);
		strcpy(reverse_end_filepath, reverse_final);
		free(reverse_final);
		free(reverse_prefix);
	}
	else
	{
		char *prefix_orig = get_fastx_prefix(single_end_filepath);
		if (prefix_orig == NULL)
		{
			printf("Your reads don't end with .fastq, .fasta, .fa, or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		char *prefix = get_filepath_in_working_dir(prefix_orig, working_dir);
		free(prefix_orig);

		if (using_fasta_format)
		{
			trim_ends_and_filter_fasta(prefix, sequence_length_threshold, trim_length);
		}
		else
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t %d -i %s -o %s_trimmed1.fastq -Q33", fastq_trimmer_threshold, single_end_filepath, prefix);
			system(buffer);
			trim_ends_and_filter_fastq(prefix, sequence_length_threshold, trim_length);
		}

		char suffix[16];
		if (using_fasta_format)
		{
			strcpy(suffix, "_trimmed2.fasta");
		}
		else
		{
			strcpy(suffix, "_trimmed2.fastq");
		}

		char *final_path = (char *)malloc(strlen(prefix) + strlen(suffix) + 1);
		sprintf(final_path, "%s%s", prefix, suffix);
		strcpy(single_end_filepath, final_path);
		free(final_path);
		free(prefix);
	}
	free(buffer);
}