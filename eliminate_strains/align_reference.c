#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>

#include "align_reference.h"
#include "file_utils.h"
#include "external/needleman_wunsch.h"

/**
 * @brief 
 * 
 * @param num_references 
 * @param msa_reference_filepaths 
 * @param bowtie2_reference_filepaths 
 * @return ReferencesData 
 */
ReferencesData align_references(char **reference_sequences_filepaths, MSA *msa_str, int num_references)
{
	int ref_idx, msa_seq_idx;

	ReferencesData references_data_str;
	references_data_str.reference_names = (char **)malloc(num_references * sizeof(char *));
	references_data_str.reference_sequence_msa_indexes = (int *)malloc(num_references * sizeof(int));
	references_data_str.reference_indexes = (int **)malloc(num_references * sizeof(int *));
	references_data_str.num_references = num_references;

	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		references_data_str.reference_names[ref_idx] = read_fastx_header_name(reference_sequences_filepaths[ref_idx]);
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
			fprintf(stderr, "Error: reference '%s' not found in MSA.\n", references_data_str.reference_names[ref_idx]);
			exit(1);
		}
		
		references_data_str.reference_indexes[ref_idx] = (int *)malloc(FASTA_MAXLINE * sizeof(int));
		if (references_data_str.reference_indexes[ref_idx] == NULL)
		{
			fprintf(stderr, "Memory allocation for reference index array failed.\n");
			exit(1);
		}
		memset(references_data_str.reference_indexes[ref_idx], -1, FASTA_MAXLINE * sizeof(int));

		char buffer[FASTA_MAXLINE];

		// read in reference sequence
		FILE *reference_sequence_file;
		if ((reference_sequence_file = fopen(reference_sequences_filepaths[ref_idx], "r")) == (FILE *)NULL)
		{
			fprintf(stderr, "Error! Cannot open MSA reference file.");
			exit(1);
		}

		char *reference_sequence = (char *)calloc(FASTA_MAXLINE, sizeof(char));
		if (!reference_sequence)
		{
			fprintf(stderr, "Memory allocation for MSA reference sequence failed\n");
			exit(1);
		}
		while (fgets(buffer, FASTA_MAXLINE, reference_sequence_file) != NULL)
		{
			if (buffer[0] != '>')
			{
				buffer[strcspn(buffer, "\r\n")] = '\0';
                strcpy(reference_sequence, buffer);
			}
		}
		fclose(reference_sequence_file);

		// read in imputed reference sequence from msa
		char *msa_imputed_reference_sequence = strdup(msa_str->sequences[references_data_str.reference_sequence_msa_indexes[ref_idx]]);

		// NOTE: needleman-wunsch may not be necessary here, can just do a linear scan/walk?
		// use needleman-wunsch alignment
		nw_aligner_t *nw = needleman_wunsch_new();
		alignment_t *result = alignment_create(256);
		int match = 1;
		int mismatch = -2;
		int gap_open = -4;
		int gap_extend = -1;
		char no_start_gap_penalty = 1;
		char no_end_gap_penalty = 1;
		char no_gaps_in_a = 0, no_gaps_in_b = 0;
		char no_mismatches = 0;
		char case_sensitive = 0;
		scoring_t scoring;
		scoring_init(&scoring, match, mismatch, gap_open, gap_extend, no_start_gap_penalty, no_end_gap_penalty, no_gaps_in_a, no_gaps_in_b, no_mismatches, case_sensitive);
		needleman_wunsch_align(reference_sequence, msa_imputed_reference_sequence, &scoring, nw, result);
		// printf("seqA: %s\n", result->result_a);
		// printf("seqB: %s\n", result->result_b);
		printf("alignment score: %i\n", result->score);

		// fill reference indicies
		// ProblematicSites *problematic_sites_str = &reference_data_str->problematic_sites_str;
		int length_alignment = strlen(result->result_b);
		int j = 0;

		int *current_reference_index = references_data_str.reference_indexes[ref_idx];

		int i, k;
		for (i = 0; i < length_alignment; i++)
		{
			if (result->result_a[i] == '-')
			{
				current_reference_index[i] = -1;
			}
			else
			{
				if (toupper(result->result_a[i]) != toupper(result->result_b[i]))
				{
					current_reference_index[i] = -1;
				}
				else
				{
					current_reference_index[i] = j;

					// TODO: implement problematic sites 
					// for (k = 0; k < problematic_sites_str->num_problematic_sites; k++)
					// {
					// 	if (i == problematic_sites_str->problematic_sites[k] - 1)
					// 	{
					// 		current_reference_index[i] = -1;
					// 	}
					// }
				}
				j++;
			}
		}

		free(reference_sequence);
		free(msa_imputed_reference_sequence);
		needleman_wunsch_free(nw);
		alignment_free(result);
	}
	
	return references_data_str;
}
