#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sam.h"

/**
 * @brief 
 * 
 * @param flag_value 
 * @return int 
 */
int parse_sam_flags(int flag_value)
{
	if (flag_value & (1 << 2))
    {
        // read was not mapped
        return -1;
    }
    else if (flag_value & (1 << 3))
    {
        // read was mapped, but mate was not mapped
        return 2;
    }
    else
    {
        // 1 if this is first read in pair, 0 if this is second read in pair
        return (flag_value & (1 << 6)) != 0;
    }
}

/**
 * @brief 
 * 
 * @param sam_results_filepath 
 * @return SAMResults 
 */
SAMResults read_in_sam_results(char **sam_results_filepaths, int num_references)
{
	int i, ref_idx;
	char buffer[FASTA_MAXLINE];

	SAMResults sam_results_str;
	sam_results_str.max_sam_line_length = 0;
	sam_results_str.sam_results = (char ***)malloc(num_references * sizeof(char **));

	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		gzFile sam_results_file;
		if ((sam_results_file = gzopen(sam_results_filepaths[ref_idx], "r")) == (gzFile)NULL)
		{
			fprintf(stderr, "SAM results file '%s' could not be opened.\n", sam_results_filepaths[ref_idx]);
			exit(1);
		}
		else
		{
			int num_sam_lines = 0;
			int max_sam_line_length = 0;

			while (gzgets(sam_results_file, buffer, FASTA_MAXLINE) != NULL)
			{
				if (buffer[0] != '@')
				{
					int sam_line_length = 0;
					for (i = 0; buffer[i] != '\n'; i++)
					{
						sam_line_length++;
					}

					if (sam_line_length > max_sam_line_length)
					{
						max_sam_line_length = sam_line_length;
					}

					num_sam_lines++;
				}
			}

			if (ref_idx == 0)
			{
				sam_results_str.num_sam_lines = num_sam_lines;
			}

			if (sam_results_str.num_sam_lines != num_sam_lines)
			{
				fprintf(stderr, "Error: Number of reads in SAM file '%s' is different than previous SAM file '%s'", sam_results_filepaths[ref_idx], sam_results_filepaths[ref_idx - 1]);
				exit(1);
			}
			if (sam_results_str.max_sam_line_length < max_sam_line_length)
			{
				sam_results_str.max_sam_line_length = max_sam_line_length;
			}

			sam_results_str.sam_results[ref_idx] = (char **)malloc(sam_results_str.num_sam_lines * sizeof(char *));
			char **current_ref_sam_results = sam_results_str.sam_results[ref_idx];
			for (int i = 0; i < sam_results_str.num_sam_lines; i++)
			{
				current_ref_sam_results[i] = (char *)calloc((max_sam_line_length + 1), sizeof(char));
			}

			gzrewind(sam_results_file);

			i = 0;
			while (gzgets(sam_results_file, buffer, FASTA_MAXLINE) != NULL)
			{
				if (buffer[0] != '@')
				{
					buffer[strcspn(buffer, "\r\n")] = '\0';
					strcpy(current_ref_sam_results[i], buffer);
					i++;
				}
			}

			gzclose(sam_results_file);
		}
	}
	return sam_results_str;
}
