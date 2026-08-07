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
 * @param sam_results_file 
 * @param sam_results 
 */
void parse_sam_info(gzFile sam_results_file, SAMResults *sam_results_str)
{
	char buffer[FASTA_MAXLINE];

	int num_sam_lines = 0;
	int max_sam_line_length = 0;

	int i;
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

	sam_results_str->num_sam_lines = num_sam_lines;
	sam_results_str->max_sam_line_length = max_sam_line_length;
}

/**
 * @brief 
 * 
 * @param sam_results_file 
 * @param sam_results_str 
 */
void read_sam_lines(gzFile sam_results_file, SAMResults *sam_results_str)
{
	char buffer[FASTA_MAXLINE];

	int i = 0;
	while (gzgets(sam_results_file, buffer, FASTA_MAXLINE) != NULL)	
	{
		if (buffer[0] != '@')
		{
			strcpy(sam_results_str->sam_results[i], buffer);
			i++;
		}
	}
}

/**
 * @brief 
 * 
 * @param sam_results_filepath 
 * @return SAMResults 
 */
SAMResults read_in_sam_results(char *sam_results_filepath)
{
	SAMResults sam_results_str;

	gzFile sam_results_file;
	if ((sam_results_file = gzopen(sam_results_filepath, "r")) == (gzFile)NULL)
	{
		fprintf(stderr, "SAM results File could not be opened.\n");

		sam_results_str.num_sam_lines = -1;
		sam_results_str.max_sam_line_length = -1;
		sam_results_str.sam_results = NULL;
	}
	else
	{
		parse_sam_info(sam_results_file, &sam_results_str);

		sam_results_str.sam_results = malloc(sam_results_str.num_sam_lines * sizeof(char *));
		for (int i = 0; i < sam_results_str.num_sam_lines; i++)
		{
			sam_results_str.sam_results[i] = malloc((sam_results_str.max_sam_line_length + 1) * sizeof(char));
		}

		gzrewind(sam_results_file);

		read_sam_lines(sam_results_file, &sam_results_str);

		gzclose(sam_results_file);
	}
	return sam_results_str;
}
