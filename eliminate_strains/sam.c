#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "sam.h"
#include "file_utils.h"

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
 * @brief Reads the alignment lines of one SAM file, keeping one line per read (single-end) or
 * one line per mate (paired-end, mate 1 first).
 *
 * Bowtie2 writes the lines of a read together (--reorder), its best alignment first. With
 * --all a read can have more lines than that; only the first line of each mate is kept, so
 * every SAM file of the same reads gives the same number of lines in the same order.
 *
 * @param num_lines output: number of lines kept
 * @param max_line_length output: length of the longest line kept
 * @param num_lines_dropped output: further lines of reads, not kept
 * @return char** the lines kept
 */
static char **read_sam_file(char *sam_filepath, char *strain_name, int using_paired_end_reads, int *num_lines, int *max_line_length, long *num_lines_dropped)
{
	gzFile sam_results_file;
	if ((sam_results_file = gzopen(sam_filepath, "r")) == (gzFile)NULL)
	{
		fprintf(stderr, "SAM results file for reference strain '%s' could not be opened.\n", strain_name);
		exit(1);
	}
	gzbuffer(sam_results_file, 1 << 20);

	int capacity = 1 << 16;
	int count = 0;
	char **lines = (char **)malloc((size_t)capacity * sizeof(char *));
	int longest = 0;
	long dropped = 0;

	// the lines kept of the read being read: [0] mate 1 (or the read), [1] mate 2
	char *pending[2] = {NULL, NULL};
	char *current_name = NULL;
	size_t current_name_capacity = 0;
	size_t current_name_length = 0;
	int have_current = 0;

	char *line = NULL;
	size_t line_capacity = 0;
	int at_end = 0;
	while (!at_end)
	{
		at_end = !gz_read_line(sam_results_file, &line, &line_capacity);
		if (!at_end && (line[0] == '@' || line[0] == '\0'))
		{
			continue;
		}

		const char *first_tab = at_end ? NULL : strchr(line, '\t');
		if (!at_end && first_tab == NULL)
		{
			fprintf(stderr, "Error: SAM file '%s' has a line without fields.\n", sam_filepath);
			exit(1);
		}
		size_t name_length = at_end ? 0 : (size_t)(first_tab - line);

		int same_read = !at_end && have_current && name_length == current_name_length && strncmp(line, current_name, name_length) == 0;
		if (!same_read && have_current)
		{
			// the previous read is complete
			if (count + 2 > capacity)
			{
				capacity *= 2;
				lines = (char **)realloc(lines, (size_t)capacity * sizeof(char *));
			}
			if (using_paired_end_reads && (pending[0] == NULL || pending[1] == NULL))
			{
				fprintf(stderr, "Error: SAM file '%s' has no line for mate %d of read '%s'.\n", sam_filepath, pending[0] == NULL ? 1 : 2, current_name);
				exit(1);
			}
			lines[count++] = pending[0];
			if (using_paired_end_reads)
			{
				lines[count++] = pending[1];
			}
			pending[0] = NULL;
			pending[1] = NULL;
			have_current = 0;
		}
		if (at_end)
		{
			break;
		}

		if (!same_read)
		{
			if (name_length + 1 > current_name_capacity)
			{
				current_name_capacity = 2 * (name_length + 1);
				current_name = (char *)realloc(current_name, current_name_capacity);
			}
			memcpy(current_name, line, name_length);
			current_name[name_length] = '\0';
			current_name_length = name_length;
			have_current = 1;
		}

		int mate = 0;
		if (using_paired_end_reads && (atoi(first_tab + 1) & (1 << 7)))
		{
			mate = 1;
		}
		if (pending[mate] != NULL)
		{
			dropped++;
			continue;
		}

		int line_length = (int)strlen(line);
		if (line_length > longest)
		{
			longest = line_length;
		}
		pending[mate] = (char *)malloc((size_t)line_length + 1);
		memcpy(pending[mate], line, (size_t)line_length + 1);
	}

	free(line);
	free(current_name);
	gzclose(sam_results_file);

	*num_lines = count;
	*max_line_length = longest;
	*num_lines_dropped = dropped;
	return lines;
}

/**
 * @brief Reads the SAM file of each reference strain (<working_dir>/<strain>.sam) and keeps
 * the reads that aligned to at least one of them
 * 
 * @return SAMResults 
 */
SAMResults read_in_sam_results(char *working_dir, char **reference_strain_names, int num_references, int using_paired_end_reads)
{
	int i, ref_idx, sam_line_idx;

	SAMResults sam_results_str;
	sam_results_str.num_sam_lines = 0;
	sam_results_str.max_sam_line_length = 0;
	sam_results_str.sam_results = (char ***)malloc(num_references * sizeof(char **));

	char *current_sam_filepath = (char *)malloc(strlen(working_dir) + FASTA_MAXLINE);
	long total_lines_dropped = 0;

	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		char *current_strain_name = reference_strain_names[ref_idx];

		sprintf(current_sam_filepath, "%s/%s.sam", working_dir, current_strain_name);

		int num_sam_lines = 0;
		int max_sam_line_length = 0;
		long num_lines_dropped = 0;
		sam_results_str.sam_results[ref_idx] = read_sam_file(current_sam_filepath, current_strain_name, using_paired_end_reads, &num_sam_lines, &max_sam_line_length, &num_lines_dropped);
		total_lines_dropped += num_lines_dropped;

		if (ref_idx == 0)
		{
			sam_results_str.num_sam_lines = num_sam_lines;
		}

		if (sam_results_str.num_sam_lines != num_sam_lines)
		{
			fprintf(stderr, "Error: Number of reads in SAM file for reference strain '%s' is different than previous reference strain '%s'\n", current_strain_name, reference_strain_names[ref_idx - 1]);
			exit(1);
		}
		if (sam_results_str.max_sam_line_length < max_sam_line_length)
		{
			sam_results_str.max_sam_line_length = max_sam_line_length;
		}
	}
	free(current_sam_filepath);

	if (total_lines_dropped > 0)
	{
		printf("Further alignments of reads that were not used (only the first, best one of each read is): %ld SAM lines\n", total_lines_dropped);
	}

	int lines_per_read = 1;
	if (using_paired_end_reads)
	{
		lines_per_read = 2;
	}

	int *read_idx_aligned = (int *)calloc((sam_results_str.num_sam_lines / lines_per_read), sizeof(int));

	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		for (sam_line_idx = 0; sam_line_idx < sam_results_str.num_sam_lines; sam_line_idx += lines_per_read)
		{
			const char *first_tab_read_1 = strchr(sam_results_str.sam_results[ref_idx][sam_line_idx], '\t');
			int parsed_sam_flags_read_1 = parse_sam_flags(atoi(first_tab_read_1 + 1));

			int read_1_aligned = 0;
			if (parsed_sam_flags_read_1 != -1)
			{
				read_1_aligned = 1;
			}

			int read_2_aligned = 0;
			if (using_paired_end_reads)
			{
				const char *first_tab_read_2 = strchr(sam_results_str.sam_results[ref_idx][sam_line_idx + 1], '\t');
				int parsed_sam_flags_read_2 = parse_sam_flags(atoi(first_tab_read_2 + 1));
				if (parsed_sam_flags_read_2 != -1)
				{
					read_2_aligned = 1;
				}
			}

			if (read_1_aligned || read_2_aligned)
			{
				read_idx_aligned[sam_line_idx / lines_per_read] = 1;
			}
		}
	}

	int num_kept_lines = 0;
	for (sam_line_idx = 0; sam_line_idx < sam_results_str.num_sam_lines; sam_line_idx += lines_per_read)
	{
		if (read_idx_aligned[sam_line_idx / lines_per_read] == 0)
		{
			for (ref_idx = 0; ref_idx < num_references; ref_idx++)
			{
				for (i = 0; i < lines_per_read; i++)
				{
					free(sam_results_str.sam_results[ref_idx][sam_line_idx + i]);
				}
			}
		}
		else
		{
			if (num_kept_lines != sam_line_idx)
			{
				for (ref_idx = 0; ref_idx < num_references; ref_idx++)
				{
					for (i = 0; i < lines_per_read; i++)
					{
						sam_results_str.sam_results[ref_idx][num_kept_lines + i] = sam_results_str.sam_results[ref_idx][sam_line_idx + i];
					}
				}
			}
			num_kept_lines += lines_per_read;
		}
	}

	printf("Input number of reads = %d\n", sam_results_str.num_sam_lines);
	printf("Number of reads kept = %d\n", num_kept_lines);

	if (num_kept_lines == 0)
	{
		fprintf(stderr, "Error: no reads aligned to any reference strain.\n");
		exit(1);
	}

	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		sam_results_str.sam_results[ref_idx] = realloc(sam_results_str.sam_results[ref_idx], num_kept_lines * sizeof(char *));
	}
	sam_results_str.num_sam_lines = num_kept_lines;

	free(read_idx_aligned);

	return sam_results_str;
}
