#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <zlib.h>

#include "debug.h"

void write_mismatch_matrix(const MismatchData *mismatch_data_str, const char *output_dir)
{
    char mismatch_matrix_filepath[1100];
	snprintf(mismatch_matrix_filepath, sizeof(mismatch_matrix_filepath), "%s/mismatch_matrix.tsv", output_dir);

	FILE *mismatch_matrix_file = fopen(mismatch_matrix_filepath, "w");
	if (mismatch_matrix_file == NULL)
	{
		fprintf(stderr, "Error: could not open '%s' to write the mismatch matrix.\n", mismatch_matrix_filepath);
		exit(1);
	}

	fprintf(mismatch_matrix_file, "read\talignment_size");
	for (int j = 0; j < mismatch_data_str->num_msa_sequences; j++)
	{
		fprintf(mismatch_matrix_file, "\t%s", mismatch_data_str->msa_sequence_names[j]);
	}
	fputc('\n', mismatch_matrix_file);

	for (int i = 0; i < mismatch_data_str->num_reads; i++)
	{
		fprintf(mismatch_matrix_file, "%s\t%d", mismatch_data_str->read_names[i], mismatch_data_str->alignment_sizes[i]);
		for (int j = 0; j < mismatch_data_str->num_msa_sequences; j++)
		{
			fprintf(mismatch_matrix_file, "\t%d", mismatch_data_str->mismatch_matrix[i][j]);
		}
		fputc('\n', mismatch_matrix_file);
	}

	if (fclose(mismatch_matrix_file) != 0)
	{
		fprintf(stderr, "Error: could not finish writing '%s'.\n", mismatch_matrix_filepath);
		exit(1);
	}
}

MismatchData read_mismatch_matrix(const char *input_mismatch_matrix_filepath, int max_num_reads)
{
	MismatchData m;
	const char *delims = " \t";
	size_t cap = 1 << 20;
	char *line = (char *)malloc(cap);
	long max_capacity = (long)max_num_reads + 1;
	if (max_capacity < 1)
	{
		max_capacity = 1;
	}
	int capacity = (int)(max_capacity < 1024 ? max_capacity : 1024);
	int expected_fields = 0;
 
	m.num_reads = 0;
	m.num_msa_sequences = 0;
	m.msa_sequence_names = NULL;
	m.read_names = (char **)malloc(capacity * sizeof(char *));
	m.alignment_sizes = (int *)malloc(capacity * sizeof(int));
	m.mismatch_matrix = (int **)malloc(capacity * sizeof(int *));
 
	gzFile input_mismatch_matrix_file = gzopen(input_mismatch_matrix_filepath, "r");
	if (input_mismatch_matrix_file == NULL)
	{
		fprintf(stderr, "Error: cannot open mismatch matrix '%s'.\n", input_mismatch_matrix_filepath);
		exit(1);
	}
 
	while (1)
	{
		/* read one whole line, however long (rows can be hundreds of kilobytes) */
		size_t len = 0;
		line[0] = '\0';
		while (gzgets(input_mismatch_matrix_file, line + len, (int)(cap - len)) != NULL)
		{
			len += strlen(line + len);
			if (len > 0 && line[len - 1] == '\n')
			{
				break;
			}
			if (len + 1 >= cap)
			{
				cap *= 2;
				line = (char *)realloc(line, cap);
			}
		}
		if (len == 0)
		{
			break; /* end of file */
		}
		while (len > 0 && (line[len - 1] == '\n' || line[len - 1] == '\r'))
		{
			line[--len] = '\0';
		}
 
		/* count its fields (separated by spaces or tabs); skip blank lines */
		int num_fields = 0, in_field = 0;
		for (char *c = line; *c != '\0'; c++)
		{
			if (*c == ' ' || *c == '\t')
			{
				in_field = 0;
			}
			else if (!in_field)
			{
				in_field = 1;
				num_fields++;
			}
		}
		if (num_fields == 0)
		{
			continue;
		}
 
		/* header: read-name column, aligned-size column, then one column per strain */
		if (m.msa_sequence_names == NULL)
		{
			char *read_col = strtok(line, delims);
			char *size_col = strtok(NULL, delims);
			if (num_fields < 3 ||
				(strcmp(read_col, "read") != 0 && strcmp(read_col, "qName") != 0) ||
				size_col == NULL ||
				(strcmp(size_col, "alignment_size") != 0 && strcmp(size_col, "blockSizes") != 0))
			{
				fprintf(stderr, "Error: '%s' does not look like a mismatch matrix: its header should start with "
						"'read alignment_size' or 'qName blockSizes'.\n", input_mismatch_matrix_filepath);
				exit(1);
			}
			expected_fields = num_fields;
			m.num_msa_sequences = num_fields - 2;
			m.msa_sequence_names = (char **)malloc(m.num_msa_sequences * sizeof(char *));
			for (int j = 0; j < m.num_msa_sequences; j++)
			{
				m.msa_sequence_names[j] = strdup(strtok(NULL, delims));
			}
			if (strcmp(m.msa_sequence_names[m.num_msa_sequences - 1], "source") == 0)
			{
				/* merge_matrices.py's "source" column is not a strain; its values are skipped below */
				m.num_msa_sequences--;
				free(m.msa_sequence_names[m.num_msa_sequences]);
			}
			if (m.num_msa_sequences == 0)
			{
				fprintf(stderr, "Error: mismatch matrix '%s' has no strain columns.\n", input_mismatch_matrix_filepath);
				exit(1);
			}
			continue;
		}
 
		/* data row: read name, aligned size, then one mismatch count per strain */
		int row = m.num_reads + 1;
		if (num_fields != expected_fields)
		{
			fprintf(stderr, "Error: '%s', data row %d has %d columns but the header has %d.\n",
					input_mismatch_matrix_filepath, row, num_fields, expected_fields);
			exit(1);
		}
		if (m.num_reads == capacity)
		{
			if (capacity >= max_capacity)
			{
				fprintf(stderr, "Error: Too many reads in the mismatch matrix '%s'!\n", input_mismatch_matrix_filepath);
				fprintf(stderr, "Max number of reads: %d (raise it with -m)\n", max_num_reads);
				fprintf(stderr, "Stopped reading at data row %d.\n", row);
				exit(1);
			}
			capacity *= 2;
			m.read_names = (char **)realloc(m.read_names, capacity * sizeof(char *));
			m.alignment_sizes = (int *)realloc(m.alignment_sizes, capacity * sizeof(int));
			m.mismatch_matrix = (int **)realloc(m.mismatch_matrix, capacity * sizeof(int *));
		}
		int r = m.num_reads;
		m.read_names[r] = strdup(strtok(line, delims));
		m.mismatch_matrix[r] = (int *)malloc(m.num_msa_sequences * sizeof(int));
		for (int j = -1; j < m.num_msa_sequences; j++) /* j = -1 is the aligned-size column */
		{
			char *token = strtok(NULL, delims);
			char *end;
			double value = strtod(token, &end);
			/* counts may be written as decimals ("10.0") but must be whole and non-negative */
			if (end == token || *end != '\0' || value != floor(value) || value < 0 || value > 2147483647.0)
			{
				fprintf(stderr, "Error: '%s', data row %d: '%s' is not a whole, non-negative number.\n",
						input_mismatch_matrix_filepath, row, token);
				exit(1);
			}
			if (j < 0)
			{
				m.alignment_sizes[r] = (int)value;
			}
			else
			{
				m.mismatch_matrix[r][j] = (int)value;
			}
		}
		m.num_reads++;
	}
	gzclose(input_mismatch_matrix_file);
	free(line);
 
	if (m.msa_sequence_names == NULL)
	{
		fprintf(stderr, "Error: mismatch matrix '%s' is empty.\n", input_mismatch_matrix_filepath);
		exit(1);
	}
	if (m.num_reads == 0)
	{
		fprintf(stderr, "Error: mismatch matrix '%s' has no data rows.\n", input_mismatch_matrix_filepath);
		exit(1);
	}
	return m;
}