#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include "needleman_wunsch.h"
#include "global.h"
void align_references(int number_of_problematic_sites, int problematic_sites[], char MSA_reference[], Options opt)
{
	FILE *Alignment_ref;
	char *line = NULL;
	size_t len = 0;
	ssize_t read;
	Alignment_ref = fopen(MSA_reference, "r");
	if (Alignment_ref == NULL)
	{
		printf("Error! Cannot open MSA reference file.");
		exit(1);
	}
	char *EPI = (char *)malloc(30000 * sizeof(char));
	int i, j;
	for (i = 0; i < 30000; i++)
	{
		EPI[i] = '\0';
	}
	while ((read = getline(&line, &len, Alignment_ref)) != -1)
	{
		j = strlen(line);
		for (i = 0; i < j; i++)
		{
			EPI[i] = line[i];
		}
	}
	fclose(Alignment_ref);
	FILE *Wuhan_file;
	// Wuhan_file = fopen("EPI_ISL_19407925.fasta","r");
	// Wuhan_file = fopen("EPI_ISL_18115435.fasta","r");
	// Wuhan_file = fopen("EPI_ISL_19407925.fasta","r");
	// Wuhan_file = fopen("/space/lenore/influenza/references/reference_sequences/cat_ref_19088566.fasta","r");
	// Wuhan_file = fopen("/space/lenore/influenza/references/reference_sequences/cat_ref_19407907.fasta","r");
	// Wuhan_file = fopen("/space/lenore/influenza/no_HA_db/reference/EPI_ISL_19407925.fasta","r");

	char *bowtieref_file = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(bowtieref_file, '\0', FASTA_MAXLINE);
	sprintf(bowtieref_file, "%s", opt.bowtie2_reference);
	Wuhan_file = fopen(bowtieref_file, "r");
	// Wuhan_file = fopen("EPI_ISL_18147004.fasta","r");
	if (Wuhan_file == NULL)
	{
		printf("Error! Cannot open MN908947.3.fasta. Please make sure this file is in the current directory.");
		exit(1);
	}
	char *sequence = malloc(30000 * sizeof(char));
	if (!sequence)
	{
		fprintf(stderr, "Memory allocation failed\n");
	}

	sequence[0] = '\0'; // initialize as empty string
	char line2[10024];

	while (fgets(line2, sizeof(line2), Wuhan_file))
	{
		if (line2[0] == '>')
		{
			// skip FASTA header line
			continue;
		}
		// strip newline
		line2[strcspn(line2, "\r\n")] = '\0';
		strcat(sequence, line2);
	}

	fclose(Wuhan_file);
	char *Wuhan = (char *)malloc(29903 * sizeof(char));
	for (i = 0; i < 29903; i++)
	{
		Wuhan[i] = '\0';
	}
	/*while((read=getline(&line, &len, Wuhan_file)) != -1){
		j=strlen(line);
		for(i=0; i<j; i++){
			Wuhan[i]=line[i];
		}
	}
	fclose(Wuhan_file);*/
	if (access("reference.fasta", F_OK) == 0)
	{
		printf("reference file exists not writing a new file...\n");
	}
	else
	{
		FILE *Wuhan_file;
		if ((Wuhan_file = fopen("MN908947.fasta", "w")) == (FILE *)NULL)
			fprintf(stderr, "Wuhan Reference File could not be opened for writing.\n");
		fprintf(Wuhan_file, ">MN908947.3\n");
		for (i = 0; i < 29903; i++)
		{
			fprintf(Wuhan_file, "%c", Wuhan[i]);
		}
		fprintf(Wuhan_file, "\n");
		fclose(Wuhan_file);
	}
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
	needleman_wunsch_align(EPI, sequence, &scoring, nw, result);
	printf("seqA: %s\n", result->result_a);
	printf("seqB: %s\n", result->result_b);
	printf("alignment score: %i\n", result->score);
	int length_alignment = strlen(result->result_b);
	// int *reference_index = (int*)malloc(length_alignment*sizeof(int));
	// for(i=0; i<length_alignment; i++){
	//	reference_index[i]=0;
	// }
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
	free(EPI);
	free(Wuhan);
	needleman_wunsch_free(nw);
	alignment_free(result);
}
