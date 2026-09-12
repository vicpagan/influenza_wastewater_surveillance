#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <zlib.h>

#include "align_references.h"
#include "file_utils.h"

/**
 * @brief 
 * 
 * @param num_references 
 * @param msa_reference_filepaths 
 * @param bowtie2_reference_filepaths 
 * @return ReferencesData 
 */
ReferencesData align_references(char **reference_sequences_filepaths, char *non_imputed_positions_filepath, MSA *msa_str, int num_references)
{
	int i, ref_idx, msa_seq_idx, site_idx;

	char buffer[FASTA_MAXLINE];

	ReferencesData references_data_str;
	references_data_str.reference_names = (char **)malloc(num_references * sizeof(char *));
	references_data_str.reference_sequence_msa_indexes = (int *)malloc(num_references * sizeof(int));
	references_data_str.references_to_msa_positions = (int **)malloc(num_references * sizeof(int *));
	references_data_str.reference_sequence_lengths = (int *)malloc(num_references * sizeof(int));
	references_data_str.num_references = num_references;

	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		gzFile reference_sequence_file;
		if ((reference_sequence_file = gzopen(reference_sequences_filepaths[ref_idx], "r")) == (gzFile)NULL)
		{
			fprintf(stderr, "Error! Cannot open reference file '%s'.\n", reference_sequences_filepaths[ref_idx]);
			exit(1);
		}

		references_data_str.reference_names[ref_idx] = (char *)calloc(FASTA_MAXLINE, sizeof(char));
		while (gzgets(reference_sequence_file, buffer, FASTA_MAXLINE) != NULL)
		{
			buffer[strcspn(buffer, "\r\n")] = '\0';
			if (buffer[0] == '>')
			{
				for (i = 1; buffer[i] != '\0'; i++)
				{
					references_data_str.reference_names[ref_idx][i - 1] = buffer[i];
				}
				references_data_str.reference_names[ref_idx][i - 1] = '\0';
				references_data_str.reference_names[ref_idx] = realloc(references_data_str.reference_names[ref_idx], i);
			}
			else
			{
				references_data_str.reference_sequence_lengths[ref_idx] = strlen(buffer);
			}
		}
		gzclose(reference_sequence_file);

		references_data_str.reference_sequence_msa_indexes[ref_idx] = -1;
		msa_seq_idx = 0;

		while (msa_seq_idx < msa_str->num_sequences && references_data_str.reference_sequence_msa_indexes[ref_idx] == -1)
		{
			if (strcmp(references_data_str.reference_names[ref_idx], msa_str->sequence_names[msa_seq_idx]) == 0)
			{
				references_data_str.reference_sequence_msa_indexes[ref_idx] = msa_seq_idx;
			}
			msa_seq_idx++;
		}

		if (references_data_str.reference_sequence_msa_indexes[ref_idx] == -1)
		{
			fprintf(stderr, "Error: Reference '%s' not found in MSA.\n", references_data_str.reference_names[ref_idx]);
			exit(1);
		}

		// read in reference sequence
		gzFile non_imputed_positions_file;
		if ((non_imputed_positions_file = gzopen(non_imputed_positions_filepath, "r")) == (gzFile)NULL)
		{
			fprintf(stderr, "Error! Cannot open non-imputed positions file.\n");
			exit(1);
		}

		int found_strain = 0;
		int processed_strain = 0;
		while (gzgets(non_imputed_positions_file, buffer, FASTA_MAXLINE) != NULL && !processed_strain)
		{
			buffer[strcspn(buffer, "\r\n")] = '\0';
			if (buffer[0] == '>')
			{
				if (strcmp(buffer + 1, references_data_str.reference_names[ref_idx]) == 0)
				{
					found_strain = 1;
				}
			}
			else
			{
				if (found_strain)
				{
					int num_positions = references_data_str.reference_sequence_lengths[ref_idx];
					references_data_str.references_to_msa_positions[ref_idx] = (int *)malloc(num_positions * sizeof(int));
					memset(references_data_str.references_to_msa_positions[ref_idx], -1, num_positions * sizeof(int));

					char *delim = strtok(buffer, ",");
					site_idx = 0;
					while (delim != NULL)
					{
						references_data_str.references_to_msa_positions[ref_idx][site_idx] = atoi(delim); // NOTE: if file is one-indexed, subtract 1
						delim = strtok(NULL, ",");
						site_idx++;
					}

					if (site_idx != num_positions)
					{
						fprintf(stderr, "Error: Not enough positions in non-imputed positions file for reference '%s'.\n", references_data_str.reference_names[ref_idx]);
						fprintf(stderr, "num_positions = %d, site_idx = %d\n", num_positions, site_idx);
						exit(1);
					}
					processed_strain = 1;
				}
			}
		}
		gzclose(non_imputed_positions_file);

		if (!processed_strain)
		{
			fprintf(stderr, "Error: Could not find reference '%s' in non-imputed positions file '%s'.\n", references_data_str.reference_names[ref_idx], non_imputed_positions_filepath);
			exit(1);
		}
	}
	
	return references_data_str;
}
