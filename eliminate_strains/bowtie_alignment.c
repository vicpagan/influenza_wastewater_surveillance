#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#include "bowtie_alignment.h"
#include "file_utils.h"
#include "global.h"

/**
 * @brief 
 * 
 * @param reference_sequences_dir 
 * @param bowtie2_indexes_dir 
 * @param reference_strain_names 
 * @param num_references 
 * @param single_end_filepath 
 * @param forward_end_filepath 
 * @param reverse_end_filepath 
 * @param working_dir 
 * @param using_paired_end_reads 
 * @param using_fasta_format 
 * @param verbose 
 */
void perform_bowtie_alignments(char *reference_sequences_dir, char* bowtie2_indexes_dir, char **reference_strain_names, int num_references, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format, int verbose)
{
	int ref_idx;
	char *buffer = (char *)calloc(FASTA_MAXLINE, sizeof(char));

	char quiet_flag[16];
	if (verbose)
	{
		strcpy(quiet_flag, "");
	}
	else
	{
		strcpy(quiet_flag, " --quiet");
	}

	char current_reference_filepath[1024];
	char current_index_prefix[1024];
	char current_index_check_path[1024];
	char current_sam_filepath[1024];

	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		char *current_strain_name = reference_strain_names[ref_idx];

		// reference sequence filepath: <reference_sequences_directory>/<strain>.fasta
		sprintf(current_reference_filepath, "%s/%s.fasta", reference_sequences_dir, current_strain_name);

		// bowtie2 index prefix: <bowtie2_indexes_directory>/<strain>
		sprintf(current_index_prefix, "%s/%s", bowtie2_indexes_dir, current_strain_name);

		// bowtie2 index check filepath: <bowtie2_indexes_directory>/<strain>.1.bt2
		sprintf(current_index_check_path, "%s.1.bt2", current_index_prefix);

		// SAM file output filepath: <working_dir>/<strain>.sam
		sprintf(current_sam_filepath, "%s/%s.sam", working_dir, current_strain_name);

		if (access(current_index_check_path, F_OK) != 0)
		{
			printf("Bowtie2 index not built for strain '%s'. Building now...\n", current_strain_name);

			sprintf(buffer, "bowtie2-build%s -f %s %s", quiet_flag, current_reference_filepath, current_index_prefix);
			system(buffer);
		}
		
		if (using_paired_end_reads && using_fasta_format)
		{
			sprintf(buffer, "bowtie2 --all%s -f -x %s -1 %s -2 %s -S %s", quiet_flag, current_index_prefix, forward_end_filepath, reverse_end_filepath, current_sam_filepath);
		}
		else if (!using_paired_end_reads && using_fasta_format)
		{
			sprintf(buffer, "bowtie2 --all%s -f -x %s -U %s -S %s", quiet_flag, current_index_prefix, single_end_filepath, current_sam_filepath);
		}
		else if (using_paired_end_reads && !using_fasta_format)
		{
			sprintf(buffer, "bowtie2 --all%s -x %s -1 %s -2 %s -S %s", quiet_flag, current_index_prefix, forward_end_filepath, reverse_end_filepath, current_sam_filepath);
		}
		else
		{
			sprintf(buffer, "bowtie2 --all%s -x %s -U %s -S %s", quiet_flag, current_index_prefix, single_end_filepath, current_sam_filepath);
		}
		system(buffer);
	}

	free(buffer);
}