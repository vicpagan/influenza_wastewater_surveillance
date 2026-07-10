#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "align_reference.h"

/**
 * @brief Aligns the MSA reference sequence to the Bowtie2 reference sequence and builds reference_index
 *
 * reference_index[i] = the column in the MSA reference that 
 * corresponds to position i in the Bowtie2 reference
 *
 * @param number_of_problematic_sites number of known problematic genome sites
 * @param problematic_sites list of known problematic genome sites
 * @param MSA_reference_path path to the MSA reference file
 * @param bowtie2_reference_path path to the Bowtie2 reference file
 */
void align_references(int number_of_problematic_sites, int problematic_sites[], char *MSA_reference_path, char *bowtie2_reference_path)
{
	// read in MSA reference sequence
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	
	FILE *MSA_reference_file = fopen(MSA_reference_path, "r");
	if (MSA_reference_file == NULL)
	{
		printf("Error! Cannot open MSA reference file.");
		exit(1);
	}

	char *MSA_reference_seq = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	int i, j;
	for (i = 0; i < FASTA_MAXLINE; i++)
	{
		MSA_reference_seq[i] = '\0';
	}
	while ((read = getline(&line, &len, MSA_reference_file)) != -1)
	{
		j = strlen(line);
		for (i = 0; i < j; i++)
		{
			MSA_reference_seq[i] = line[i];
		}
	}
	fclose(MSA_reference_file);

	// read in bowtie2 reference sequence
	FILE *bowtie2_reference_file = fopen(bowtie2_reference_path, "r");
	if (bowtie2_reference_file == NULL)
	{
		printf("Error! Cannot open bowtie2 reference file. Please make sure this file is in the current directory.");
		exit(1);
	}
	char *bowtie2_reference_seq = malloc(FASTA_MAXLINE * sizeof(char));
	if (!bowtie2_reference_seq)
	{
		fprintf(stderr, "Memory allocation failed\n");
	}

	bowtie2_reference_seq[0] = '\0'; // initialize as empty string
	char line2[FASTA_MAXLINE];

	while (fgets(line2, sizeof(line2), bowtie2_reference_file))
	{
		if (line2[0] == '>')
		{
			// skip FASTA header line
			continue;
		}
		// strip newline
		line2[strcspn(line2, "\r\n")] = '\0';
		strcat(bowtie2_reference_seq, line2);
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
	needleman_wunsch_align(MSA_reference_seq, bowtie2_reference_seq, &scoring, nw, result);
	printf("seqA: %s\n", result->result_a);
	printf("seqB: %s\n", result->result_b);
	printf("alignment score: %i\n", result->score);

	// fill reference indicies
	int length_alignment = strlen(result->result_b);
	j = 0;
	int k = 0;
	for (i = 0; i < length_alignment; i++)
	{
		if (result->result_a[i] == '-')
		{
			reference_index[i] = -1;
		}
		else
		{
			if (result->result_a[i] != result->result_b[i])
			{
				reference_index[i] = -1;
			}
			else
			{
				reference_index[i] = j;
				for (k = 0; k < number_of_problematic_sites; k++)
				{
					if (i == problematic_sites[k] - 1)
					{
						reference_index[i] = -1;
					}
				}
			}
			j++;
		}
	}

	free(MSA_reference_seq);
	free(bowtie2_reference_seq);
	free(line);
	needleman_wunsch_free(nw);
	alignment_free(result);
}
