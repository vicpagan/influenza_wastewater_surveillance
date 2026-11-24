#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "bowtie_alignment.h"
#include "file_utils.h"
#include "global.h"

void run_bowtie2(char *index_prefix, char *sam_filepath, int report_all, char *bowtie2_args, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, int using_paired_end_reads, int using_fasta_format, int num_threads, int verbose)
{
	size_t command_size = strlen(index_prefix) + strlen(sam_filepath) + strlen(bowtie2_args) + strlen(single_end_filepath) + strlen(forward_end_filepath) + strlen(reverse_end_filepath) + 256;
	char *command = (char *)malloc(command_size);

	int length = snprintf(command, command_size, "bowtie2%s%s%s -p %d --reorder", report_all ? " --all" : "", verbose ? "" : " --quiet", using_fasta_format ? " -f" : "", num_threads);
	if (bowtie2_args[0] != '\0')
	{
		length += snprintf(command + length, command_size - length, " %s", bowtie2_args);
	}
	if (using_paired_end_reads)
	{
		snprintf(command + length, command_size - length, " -x %s -1 %s -2 %s -S %s", index_prefix, forward_end_filepath, reverse_end_filepath, sam_filepath);
	}
	else
	{
		snprintf(command + length, command_size - length, " -x %s -U %s -S %s", index_prefix, single_end_filepath, sam_filepath);
	}

	// a SAM file left by an earlier run must not be read if this run fails
	remove(sam_filepath);

	fflush(stdout);
	int status = system(command);
	if (status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
	{
		fprintf(stderr, "Error: Bowtie2 did not finish successfully%s. The command was:\n%s\n", verbose ? "" : " (run again with -v to see its messages)", command);
		exit(1);
	}

	free(command);
}

/**
 * @brief Aligns the reads to each reference strain, writing <working_dir>/<strain>.sam. The
 * Bowtie2 index of a strain is built if it is missing.
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
 * @param bowtie2_args extra Bowtie2 options, may be empty
 */
void perform_bowtie_alignments(char *reference_sequences_dir, char* bowtie2_indexes_dir, char **reference_strain_names, int num_references, char *single_end_filepath, char *forward_end_filepath, char *reverse_end_filepath, char *working_dir, int using_paired_end_reads, int using_fasta_format, int num_threads, int verbose, char *bowtie2_args)
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

	char current_reference_filepath[4200];
	char current_index_prefix[4200];
	char current_index_check_path[4300];
	char current_sam_filepath[4200];

	for (ref_idx = 0; ref_idx < num_references; ref_idx++)
	{
		char *current_strain_name = reference_strain_names[ref_idx];

		// reference sequence filepath: <reference_sequences_directory>/<strain>.fasta
		snprintf(current_reference_filepath, sizeof(current_reference_filepath), "%s/%s.fasta", reference_sequences_dir, current_strain_name);

		// bowtie2 index prefix: <bowtie2_indexes_directory>/<strain>
		snprintf(current_index_prefix, sizeof(current_index_prefix), "%s/%s", bowtie2_indexes_dir, current_strain_name);

		// bowtie2 index check filepath: <bowtie2_indexes_directory>/<strain>.1.bt2
		snprintf(current_index_check_path, sizeof(current_index_check_path), "%s.1.bt2", current_index_prefix);

		// SAM file output filepath: <working_dir>/<strain>.sam
		snprintf(current_sam_filepath, sizeof(current_sam_filepath), "%s/%s.sam", working_dir, current_strain_name);

		if (access(current_index_check_path, F_OK) != 0)
		{
			printf("Bowtie2 index not built for strain '%s'. Building now...\n", current_strain_name);
			fflush(stdout);

			snprintf(buffer, FASTA_MAXLINE, "bowtie2-build%s -f %s %s", quiet_flag, current_reference_filepath, current_index_prefix);
			system(buffer);
		}

		// every alignment of a read is reported (--all); sam.c keeps the first (best) one
		run_bowtie2(current_index_prefix, current_sam_filepath, 1, bowtie2_args, single_end_filepath, forward_end_filepath, reverse_end_filepath, using_paired_end_reads, using_fasta_format, num_threads, verbose);
	}

	free(buffer);
}
