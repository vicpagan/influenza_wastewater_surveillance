#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "align_reference.h"

/**
 * @brief Aligns the MSA reference sequence to the Bowtie2 reference sequence and builds a reference index array
 *
 * reference_index[i] = the column in the MSA reference that 
 * corresponds to position i in the Bowtie2 reference
 * 
 * @param reference_data_str 
 * @param msa_reference_filepath 
 * @param bowtie2_reference_filepath 
 */
void align_reference(ReferenceData *reference_data_str, char *msa_reference_filepath, char *bowtie2_reference_filepath)
{
	reference_data_str->reference_index = (int *)calloc(FASTA_MAXLINE, sizeof(int));
	if (reference_data_str->reference_index == NULL)
	{
		fprintf(stderr, "Memory allocation for reference index array failed.\n");
		exit(1);
	}

	char buffer[FASTA_MAXLINE];

	// read in msa reference sequence
	FILE *msa_reference_file;
	if ((msa_reference_file = fopen(msa_reference_filepath, "r")) == (FILE *)NULL)
	{
		fprintf(stderr, "Error! Cannot open MSA reference file.");
		exit(1);
	}

	char *msa_reference_sequence = (char *)calloc(FASTA_MAXLINE, sizeof(char));
	if (!msa_reference_sequence)
	{
		fprintf(stderr, "Memory allocation for MSA reference sequence failed\n");
		exit(1);
	}
	while (fgets(buffer, FASTA_MAXLINE, msa_reference_file) != NULL)
	{
		if (buffer[0] != '>')
		{
			strcpy(msa_reference_sequence, buffer);
		}
	}
	fclose(msa_reference_file);

	// read in bowtie2 reference sequence
	FILE *bowtie2_reference_file;
	if ((bowtie2_reference_file = fopen(bowtie2_reference_filepath, "r")) == (FILE *)NULL)
	{
		fprintf(stderr, "Error! Cannot open Bowtie2 reference file.");
		exit(1);
	}

	char *bowtie2_reference_sequence = (char *)calloc(FASTA_MAXLINE, sizeof(char));
	if (!bowtie2_reference_sequence)
	{
		fprintf(stderr, "Memory allocation for Bowtie2 reference sequence failed\n");
		exit(1);
	}
	while (fgets(buffer, FASTA_MAXLINE, bowtie2_reference_file) != NULL)
	{
		if (buffer[0] != '>')
		{
			strcpy(bowtie2_reference_sequence, buffer);
		}
	}
	fclose(bowtie2_reference_file);

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
	needleman_wunsch_align(msa_reference_sequence, bowtie2_reference_sequence, &scoring, nw, result);
	printf("seqA: %s\n", result->result_a);
	printf("seqB: %s\n", result->result_b);
	printf("alignment score: %i\n", result->score);

	// fill reference indicies
	// ProblematicSites *problematic_sites_str = &reference_data_str->problematic_sites_str;
	int length_alignment = strlen(result->result_b);
	int j = 0;

	int i, k;
	for (i = 0; i < length_alignment; i++)
	{
		if (result->result_a[i] == '-')
		{
			reference_data_str->reference_index[i] = -1;
		}
		else
		{
			if (result->result_a[i] != result->result_b[i])
			{
				reference_data_str->reference_index[i] = -1;
			}
			else
			{
				reference_data_str->reference_index[i] = j;

				// TODO: implement problematic sites 
				// for (k = 0; k < problematic_sites_str->num_problematic_sites; k++)
				// {
				// 	if (i == problematic_sites_str->problematic_sites[k] - 1)
				// 	{
				// 		reference_data_str->reference_index[i] = -1;
				// 	}
				// }
			}
			j++;
		}
	}

	free(msa_reference_sequence);
	free(bowtie2_reference_sequence);
	needleman_wunsch_free(nw);
	alignment_free(result);
}

/**
 * @brief
 * 
 * @param reference_data_strs
 * @param num_references
 * @param msa_reference_filepaths
 * @param bowtie2_reference_filepaths
 */
void align_all_references(ReferenceData **reference_data_strs, int num_references, char **msa_reference_filepaths, char **bowtie2_reference_filepaths)
{
	int ref_idx;
	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		align_reference(reference_data_strs[ref_idx], msa_reference_filepaths[ref_idx], bowtie2_reference_filepaths[ref_idx]);
	}
}
