#include <stdio.h>
#include <stdlib.h>

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