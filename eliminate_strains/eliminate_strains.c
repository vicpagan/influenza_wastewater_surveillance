#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <zlib.h>
#include <assert.h>
#include <unistd.h>
#include <pthread.h>
#include <dirent.h>

#include "global.h"
#include "options.h"
#include "align_reference.h"
#include "calculate_allele_freq.h"
#include "write_mismatch_matrix.h"


/**
 * @brief Parse the MSA file for metadata information
 * 
 * @param msa_file The MSA file to read from
 * @param msa Reference to the stored MSA instance
 */
void parse_msa_info(gzFile msa_file, MSA *msa)
{
	char buffer[FASTA_MAXLINE];

	int num_sequences = 0;
	int sequence_length = 0;
	int max_sequence_name_length = 0;

	int i;
	while (gzgets(msa_file, buffer, FASTA_MAXLINE) != NULL)
    {
        if (buffer[0] == '>')
        {
			// read in sequence name line
            int sequence_name_length = 0;
			for (i = 1; buffer[i] != '\n'; i++)
			{
				sequence_name_length++;
			}

            if (sequence_name_length > max_sequence_name_length)
            {
                max_sequence_name_length = sequence_name_length;
            }

            num_sequences++;
        }
        else if (num_sequences == 1)
        {
			// read in sequence line
            for (i = 0; buffer[i] != '\n'; i++)
            {
                sequence_length++;
            }
        }
    }

	msa->num_sequences = num_sequences;
	msa->sequence_length = sequence_length;
	msa->max_sequence_name_length = max_sequence_name_length;
}

/**
 * @brief Reads in the MSA sequences and sequence names from the input file to the given MSA struct
 * 
 * @param msa_file The MSA file to read from
 * @param msa The reference to the stored MSA instance
 * @param ref_seq_name The name of a reference sequence
 */
void read_msa_sequences(gzFile msa_file, MSA *msa_str, const char *reference_sequence_name)
{
	char buffer[FASTA_MAXLINE];

	int seq_idx = -1;
	// int found_reference = 0;

	int i;
	while (gzgets(msa_file, buffer, FASTA_MAXLINE) != NULL)
	{
		if (buffer[0] == '>')
		{
			// read in sequence name line
			seq_idx++;

			for (i = 1; buffer[i] != '\n'; i++)
			{
				msa_str->sequence_names[seq_idx][i - 1] = buffer[i];
			}
			msa_str->sequence_names[seq_idx][i - 1] = '\0';

			// if (strcmp(msa_str->sequence_names[seq_idx], reference_sequence_name) == 0)
			// {
			// 	found_reference = 1;
			// }
		}
		else
		{
			// read in sequence line
			int length = strlen(buffer);

			for (i = 0; i < length; i++)
			{
				switch (buffer[i])
                {
                    case 'A':
                    case 'a':
                        msa_str->sequences[seq_idx][i] = 'A';
						// allele_frequency[i][0]++;
                        break;

                    case 'C':
                    case 'c':
                        msa_str->sequences[seq_idx][i] = 'C';
						// allele_frequency[i][1]++;
                        break;

                    case 'G':
                    case 'g':
                        msa_str->sequences[seq_idx][i] = 'G';
						// allele_frequency[i][2]++;
                        break;

                    case 'T':
                    case 't':
                        msa_str->sequences[seq_idx][i] = 'T';
						// allele_frequency[i][3]++;
                        break;

                    case '-':
                        msa_str->sequences[seq_idx][i] = '-';
                        break;

                    default:
                        msa_str->sequences[seq_idx][i] = '\0';
                        break;
                }
			}

			// if (found_reference)
			// {
			// 	msa_str->reference_sequence_index = index;
			// 	found_reference = 0;
			// }
		}
	}
}

/**
 * @brief Reads in MSA and its important information into a MSA struct instance created by the function
 * 
 * @param msa_file The MSA file to read from
 * @param ref_seq_name The name of a reference sequence
 * @return MSA Instance of an MSA struct that stores the MSA data
 */
MSA read_in_msa(char *msa_filepath, const char *reference_sequence_name)
{
	MSA msa_str;
	
	gzFile msa_file;
	if ((msa_file = gzopen(msa_filepath, "r")) == (gzFile)NULL)
	{
		fprintf(stderr, "MSA File could not be opened.\n");

		msa_str.num_sequences = -1;
		msa_str.sequence_length = -1;
		msa_str.max_sequence_name_length = -1;
		// msa_str.reference_sequence_index = -1;

		msa_str.sequence_names = NULL;
		msa_str.sequences = NULL;
	}
	else
	{
		// parse and store MSA metadata info
		parse_msa_info(msa_file, &msa_str);

		// allocate storage for sequences and their names
		msa_str.sequence_names = malloc(msa_str.num_sequences * sizeof(char *));
		msa_str.sequences = malloc(msa_str.num_sequences * sizeof(char *));
		for (int i = 0; i < msa_str.num_sequences; i++)
		{
			msa_str.sequence_names[i] = malloc((msa_str.max_sequence_name_length + 1) * sizeof(char));
			msa_str.sequences[i] = malloc((msa_str.sequence_length + 1) * sizeof(char));
		}

		gzrewind(msa_file);

		// parse and store MSA sequences and their names
		read_msa_sequences(msa_file, &msa_str, reference_sequence_name);

		gzclose(msa_file);
	}
    return msa_str;
}

/**
 * @brief Finds identical sequences in the MSA and marks them for removal
 * 
 * @param msa_str The reference to the stored MSA instance
 * @param sequences_to_remove Direct-indexed array of sequences to remove ((sequences_to_remove[i] == 1) = remove sequence i)
 * 
 * @return The number of sequences marked for removal
 */
int remove_identical_sequences(MSA *msa_str, int *sequences_to_remove)
{
	int num_sequences = msa_str->num_sequences;
	int sequence_length = msa_str->sequence_length;

	int num_sequences_removed = 0;
	
	// compare every sequence to each other and mark rightmost duplicates for removal
	int i, j, k, mismatch;
	for (i = 0; i < num_sequences; i++)
	{
		if (sequences_to_remove[i] == 0)
		{
			for (j = i + 1; j < num_sequences; j++)
			{
				if (sequences_to_remove[j] == 0)
				{
					mismatch = 0;
					for (k = 0; k < sequence_length; k++)
					{
						if (msa_str->sequences[i][k] != msa_str->sequences[j][k])
						{
							mismatch++;
						}
					}

					if (mismatch == 0)
					{
						sequences_to_remove[j] = 1;
						num_sequences_removed++;
					}
				}
			}
		}	
	}

	return num_sequences_removed;
}

/**
 * @brief Prunes the list of sequences and sequence names in the stored MSA based on a provided list of sequences to remove
 * 
 * @param msa_str The reference to the stored MSA instance
 * @param sequences_to_remove Direct-indexed array of sequences to remove ((sequences_to_remove[i] == 1) = remove sequence i)
 */
void prune_msa_sequences(MSA *msa_str, int *sequences_to_remove)
{
	int num_sequences = msa_str->num_sequences;

	int l = 0; // left pointer (writing)
	int r; // right pointer (reading)
	for (r = 0; r < num_sequences; r++)
	{
		// keep the sequence at position r
		if (sequences_to_remove[r] == 0)
		{
			if (l != r)
			{
				msa_str->sequences[l] = msa_str->sequences[r];
				msa_str->sequence_names[l] = msa_str->sequence_names[r];
			}
			l++;
		}
	}
	// at this point, l represents the number of sequences remaining
	
	// free up remaining space past what we are keeping
	for (r = l; r < num_sequences; r++)
	{
		free(msa_str->sequences[r]);
		free(msa_str->sequence_names[r]);
	}
	msa_str->sequences = realloc(msa_str->sequences, l * sizeof(char *));
	msa_str->sequence_names = realloc(msa_str->sequence_names, l * sizeof(char *));

	msa_str->num_sequences = l;
}

/**
 * @brief 
 * 
 * @param variant_sites_filepath
 * @param num_references 
 * @return VariantSites
 */
VariantSites read_in_variant_sites(char *variant_sites_filepath)
{
	char buffer[16];

	VariantSites variant_sites_str;

	FILE *variant_sites_file;
	if ((variant_sites_file = fopen(variant_sites_filepath, "r")) == (FILE *)NULL)
	{
		fprintf(stderr, "variant sites File could not be opened.\n");

		variant_sites_str.num_variant_sites = -1;
		variant_sites_str.variant_sites = NULL;
	}
	else
	{
		int first_line = 1;
		int i = 0;
		while (fgets(buffer, 16, variant_sites_file) != NULL)
		{
			if (first_line)
			{
				variant_sites_str.num_variant_sites = atoi(buffer);
				variant_sites_str.variant_sites = (int *)malloc(variant_sites_str.num_variant_sites * sizeof(int));

				first_line = 0;
			}
			else
			{
				variant_sites_str.variant_sites[i] = atoi(buffer);
				i++;
			}
		}
	}
	fclose(variant_sites_file);

	return variant_sites_str;
}

// TODO: reimplement with a new variant-sites-like file format for problematic sites
// ProblematicSites *read_in_problematic_sites(char **problematic_sites_filepaths, int num_refs){
// 	FILE* file;
// 	if (( file = fopen("problematic_sites_sarsCov2.vcf","r")) == (FILE *) NULL ) fprintf(stderr, "Problematic Sites File could not be opened.\n");
// 	char buffer[1000];
// 	char name[30];
// 	int position;
// 	char ch1[1];
// 	char ch2[1];
// 	char ch3[1];
// 	char ch4[1];
// 	char s1[10];
// 	char s2[30];
// 	int i=0;
// 	while( fgets(buffer,1000,file) != NULL){
// 		if ( buffer[0] != '#' ){
// 			sscanf(buffer,"%s\t%d\t%c\t%c\t%c\t%c\t%s\t%s",&name,&position,&ch1,&ch2,&ch3,&ch4,&s1,&s2);
// 			problematic_sites[i]=position;
// 			i++;
// 		}
// 	}
// 	fclose(file);
// 	return i;
// }

/**
 * @brief Parses the SAM flags for the state of the current read
 * 
 * @param flag_value Base 10 representation of the flag bits
 * @return int Returns 1 if first read in pair, 0 if second read in pair, -1 if this read was not mapped, 2 if this read was mapped but this read's mate was not
 */
int parse_sam_flags(int flag_value)
{
	int bit_string[32] = {0};

	// transform flag into binary string
	int i = 0;
	while (flag_value > 0)
	{
		bit_string[i] = flag_value % 2;
		flag_value >>= 1;
		i++;
	}

	if (bit_string[2] == 1)
	{
		// read was not mapped
		return -1;
	}
	else if (bit_string[3] == 1)
	{
		// read was mapped, but mate was not mapped
		return 2;
	}
	else
	{
		// 1 if this is first read in pair, 0 if this is second read in pair
		return bit_string[6];
	}
}

/**
 * @brief 
 * 
 * @param sam_results_file 
 * @param sam_results 
 */
void parse_sam_info(gzFile sam_results_file, SAMResults *sam_results)
{
	char buffer[FASTA_MAXLINE];

	int num_sam_lines = 0;
	int max_sam_line_length = 0;

	int i;
	while (gzgets(sam_results_file, buffer, FASTA_MAXLINE) != NULL)
    {
        if (buffer[0] != '@')
        {
            int sam_line_length = 0;
			for (i = 0; buffer[i] != '\n'; i++)
			{
				sam_line_length++;
			}

            if (sam_line_length > max_sam_line_length)
            {
                max_sam_line_length = sam_line_length;
            }

            num_sam_lines++;
        }
    }

	sam_results->num_sam_lines = num_sam_lines;
	sam_results->max_sam_line_length = max_sam_line_length;
}

/**
 * @brief 
 * 
 * @param sam_results_file 
 * @param sam_results_str 
 */
void read_sam_lines(gzFile sam_results_file, SAMResults *sam_results_str)
{
	char buffer[FASTA_MAXLINE];

	int i = 0;
	while (gzgets(buffer, FASTA_MAXLINE, sam_results_file) != NULL)
	{
		if (buffer[0] != '@')
		{
			strcpy(sam_results_str->sam_results[i], buffer);
			i++;
		}
	}
}

/**
 * @brief 
 * 
 * @param sam_results_filepath 
 * @return SAMResults 
 */
SAMResults read_in_sam_results(char *sam_results_filepath)
{
	SAMResults sam_results_str;

	gzFile sam_results_file;
	if ((sam_results_file = gzopen(sam_results_filepath, "r")) == (gzFile)NULL)
	{
		fprintf(stderr, "SAM results File could not be opened.\n");

		sam_results_str.num_sam_lines = -1;
		sam_results_str.max_sam_line_length = -1;
		sam_results_str.sam_results = NULL;
	}
	else
	{
		parse_sam_info(sam_results_file, &sam_results_str);

		sam_results_str.sam_results = malloc(sam_results_str.num_sam_lines * sizeof(char *));
		for (int i = 0; i < sam_results_str.num_sam_lines; i++)
		{
			sam_results_str.sam_results[i] = malloc((sam_results_str.max_sam_line_length + 1) * sizeof(char));
		}

		gzrewind(sam_results_file);

		read_sam_lines(sam_results_file, &sam_results_str);

		gzclose(sam_results_file);
	}
	return sam_results_str;
}


/**
 * @brief Nudge a thread's [start, end) line range so it doesn't split a read pair.
 *
 * If `start` lands on a second-in-pair SAM record, backs up one line so the
 * pair starts together; if `end` lands on a first-in-pair record, extends by
 * one line so the pair finishes together.
 *
 * @param start Proposed first line index (into sam_results) for this thread.
 * @param end Proposed last line index (into sam_results) for this thread.
 * @param return_arr Output: [0] = adjusted start, [1] = adjusted end.
 * @param max_sam_length Length of the scratch line buffer to allocate.
 */
void adjust_start_end(int start, int end, int *return_arr, int max_sam_length)
{
	char *s;
	char *buffer_copy = (char *)malloc((max_sam_length + 1) * sizeof(char));
	strcpy(buffer_copy, sam_results[start]);
	s = strtok(buffer_copy, "\t");
	s = strtok(NULL, "\t");
	int decimal = 0;
	sscanf(s, "%d", &decimal);
	decimal = parse_sam_flags(decimal);
	if (decimal == 0)
	{ // second in pair
		return_arr[0] = start - 1;
	}
	else
	{
		return_arr[0] = start;
	}
	strcpy(buffer_copy, sam_results[end]);
	s = strtok(buffer_copy, "\t");
	s = strtok(NULL, "\t");
	decimal = 0;
	sscanf(s, "%d", &decimal);
	decimal = parse_sam_flags(decimal);
	if (decimal == 1)
	{
		return_arr[1] = end + 1;
	}
	else
	{
		return_arr[1] = end;
	}
	free(buffer_copy);
}

/**
 * @brief 
 * 
 * @param filename 
 */
void trim_ends_fastq(const char *filename)
{
	char input_filepath[1000];
	char output_filepath[1000];

	FILE *input_file;
	FILE *output_file;
	
	strcpy(input_filepath, filename);
	strcpy(output_filepath, filename);
	strcat(input_filepath, "_trimmed1.fastq");
	strcat(output_filepath, "_trimmed2.fastq");
	
	if ((input_file = fopen(input_filepath, "r")) == (FILE *)NULL)
	{
		printf("Cannot open %s\n", input_filepath);
		exit(1);
	}
	if ((output_file = fopen(output_filepath, "w")) == (FILE *)NULL)
	{
		printf("Cannot open %s\n", output_filepath);
		exit(1);
	}
	
	char *strain_name = (char *)calloc(FASTA_MAXLINE, sizeof(char));
	char *sequence = (char *)calloc(FASTA_MAXLINE, sizeof(char));
	char *plus = (char *)calloc(2, sizeof(char));
	char *quality = (char *)calloc(FASTA_MAXLINE,  sizeof(char));

	size_t len = 0;
	int i;
	while (getline(&strain_name, &len, input_file) != -1 && getline(&sequence, &len, input_file) != -1 && getline(&plus, &len, input_file) != -1 && getline(&quality, &len, input_file) != -1)
    {
		int seq_length = strlen(sequence);

		fprintf(output_file, "%s", strain_name);

		if (seq_length > 95)
		{
			for (i = 15; i < (seq_length - 15); i++)
			{
				fprintf(output_file, "%c", sequence[i]);
			}
			fprintf(output_file, "\n");
		}
		else 
		{
			fprintf(output_file, "%s", sequence);
		}

		fprintf(output_file, "%s", plus);

		if (seq_length > 95)
		{
			for (i = 15; i < (seq_length - 15); i++)
			{
				fprintf(output_file, "%c", quality[i]);
			}
			fprintf(output_file, "\n");
		}
		else 
		{
			fprintf(output_file, "%s", quality);
		}
	}
	fclose(input_file);
	fclose(output_file);
	
	free(strain_name);
	free(sequence);
	free(plus);
	free(quality);
}

/**
 * @brief 
 * 
 * @param filename 
 */
void trim_ends_fasta(const char *filename)
{
	char input_filepath[1000];
	char output_filepath[1000];

	FILE *input_file;
	FILE *output_file;
	
	strcpy(input_filepath, filename);
	strcpy(output_filepath, filename);
	strcat(input_filepath, "_trimmed1.fasta");
	strcat(output_filepath, "_trimmed2.fasta");
	
	if ((input_file = fopen(input_filepath, "r")) == (FILE *)NULL)
	{
		printf("Cannot open %s\n", input_filepath);
		exit(1);
	}
	if ((output_file = fopen(output_filepath, "w")) == (FILE *)NULL)
	{
		printf("Cannot open %s\n", output_filepath);
		exit(1);
	}
	
	char *strain_name = (char *)calloc(FASTA_MAXLINE, sizeof(char));
	char *sequence = (char *)calloc(FASTA_MAXLINE, sizeof(char));

	size_t len = 0;
	int i;
	while (getline(&strain_name, &len, input_file) != -1 && getline(&sequence, &len, input_file) != -1)
    {
		int seq_length = strlen(sequence);

		fprintf(output_file, "%s", strain_name);

		if (seq_length > 95)
		{
			for (i = 15; i < (seq_length - 15); i++)
			{
				fprintf(output_file, "%c", sequence[i]);
			}
			fprintf(output_file, "\n");
		}
		else 
		{
			fprintf(output_file, "%s", sequence);
		}
	}
	fclose(input_file);
	fclose(output_file);
	
	free(strain_name);
	free(sequence);
}

/**
 * @brief Quality-trim and end-trim the input reads (single-end or paired), in place.
 *
 * For FASTQ input, runs `fastq_quality_trimmer` then trim_ends_fastq(); for
 * FASTA input, just runs trim_ends_fasta(). Rewrites single_end_file (or
 * forward_end_file/reverse_end_file) to point at the trimmed output --
 * since these are caller-owned buffers passed by pointer, the caller sees
 * the update (unlike the previous Options-by-value version, where this
 * rewrite was silently lost on return).
 * 
 * @param paired 1 if using paired-end reads.
 * @param fasta_format 1 if reads are FASTA (vs. FASTQ).
 * @param single_end_file Single-end reads path; rewritten in place if paired == 0.
 * @param forward_end_file Forward reads path; rewritten in place if paired == 1.
 * @param reverse_end_file Reverse reads path; rewritten in place if paired == 1.
 */
void clean_reads(int paired, int fasta_format, char *single_end_file, char *forward_end_file, char *reverse_end_file)
{
	int i;
	char *buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(buffer, '\0', FASTA_MAXLINE);
	if (paired == 0)
	{
		char *prefix;
		if (single_end_file[strlen(single_end_file) - 6] == '.' && single_end_file[strlen(single_end_file) - 5] == 'f' && single_end_file[strlen(single_end_file) - 4] == 'a' && single_end_file[strlen(single_end_file) - 3] == 's' && single_end_file[strlen(single_end_file) - 2] == 't' && (single_end_file[strlen(single_end_file) - 1] == 'q' || single_end_file[strlen(single_end_file) - 1] == 'a'))
		{
			prefix = (char *)malloc((strlen(single_end_file) - 6 + 15) * sizeof(char));
			for (i = 0; i < strlen(single_end_file) - 6; i++)
			{
				prefix[i] = single_end_file[i];
			}
			for (i = strlen(single_end_file) - 6; i < strlen(single_end_file) - 6 + 15; i++)
			{
				prefix[i] = '\0';
			}
		}
		else if (single_end_file[strlen(single_end_file) - 3] == '.' && single_end_file[strlen(single_end_file) - 2] == 'f' && (single_end_file[strlen(single_end_file) - 1] == 'q' || single_end_file[strlen(single_end_file) - 1] == 'a'))
		{
			prefix = (char *)malloc((strlen(single_end_file) - 3 + 15) * sizeof(char));
			for (i = 0; i < strlen(single_end_file) - 3; i++)
			{
				prefix[i] = single_end_file[i];
			}
			for (i = strlen(single_end_file) - 3; i < strlen(single_end_file) - 3 + 15; i++)
			{
				prefix[i] = '\0';
			}
		}
		else
		{
			printf("Your reads don't end with .fastq, .fasta, .fa, or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		if (fasta_format == 0)
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t 35 -i %s -o %s_trimmed1.fastq -Q33", single_end_file, prefix);
			system(buffer);
			trim_ends_fastq(prefix);
		}
		else
		{
			trim_ends_fasta(prefix);
		}
		// sprintf(buffer, "fastx_trimmer -m 65 -t 15 -i %s_trimmed1.fastq -o %s_trimmed2.fastq",single_end_file,prefix);
		// system(buffer);
		// sprintf(buffer, "fastx_trimmer -f 15 -i %s_trimmed2.fastq -o %s_trimmed3.fastq",single_end_file,prefix);
		// system(buffer);
		char suffix[1000];
		if (fasta_format == 0)
		{
			strcpy(suffix, "_trimmed2.fastq");
		}
		else
		{
			strcpy(suffix, "_trimmed2.fasta");
		}
		strcat(prefix, suffix);
		strcpy(single_end_file, prefix);
		free(prefix);
	}
	else
	{
		char *prefix_forward;
		if (forward_end_file[strlen(forward_end_file) - 6] == '.' && forward_end_file[strlen(forward_end_file) - 5] == 'f' && forward_end_file[strlen(forward_end_file) - 4] == 'a' && forward_end_file[strlen(forward_end_file) - 3] == 's' && forward_end_file[strlen(forward_end_file) - 2] == 't' && (forward_end_file[strlen(forward_end_file) - 1] == 'q' || forward_end_file[strlen(forward_end_file) - 1] == 'a'))
		{
			prefix_forward = (char *)malloc((strlen(forward_end_file) - 6 + 15) * sizeof(char));
			for (i = 0; i < strlen(forward_end_file) - 6; i++)
			{
				prefix_forward[i] = forward_end_file[i];
			}
			for (i = strlen(forward_end_file) - 6; i < strlen(forward_end_file) - 6 + 15; i++)
			{
				prefix_forward[i] = '\0';
			}
		}
		else if (forward_end_file[strlen(forward_end_file) - 3] == '.' && forward_end_file[strlen(forward_end_file) - 2] == 'f' && (forward_end_file[strlen(forward_end_file) - 1] == 'q' || forward_end_file[strlen(forward_end_file) - 1] == 'a'))
		{
			prefix_forward = (char *)malloc((strlen(forward_end_file) - 3 + 15) * sizeof(char));
			for (i = 0; i < strlen(forward_end_file) - 3; i++)
			{
				prefix_forward[i] = forward_end_file[i];
			}
			for (i = strlen(forward_end_file) - 3; i < strlen(forward_end_file) - 3 + 15; i++)
			{
				prefix_forward[i] = '\0';
			}
		}
		else
		{
			printf("Your reads don't end with .fastq, .fasta, .fa or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		if (fasta_format == 0)
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t 35 -i %s -o %s_trimmed1.fastq -Q33", forward_end_file, prefix_forward);
			system(buffer);
			trim_ends_fastq(prefix_forward);
		}
		else
		{
			trim_ends_fasta(prefix_forward);
		}
		char suffix[1000];
		if (fasta_format == 0)
		{
			strcpy(suffix, "_trimmed2.fastq");
		}
		else
		{
			strcpy(suffix, "_trimmed2.fasta");
		}
		strcat(prefix_forward, suffix);
		strcpy(forward_end_file, prefix_forward);
		free(prefix_forward);
		char *prefix_reverse;
		if (reverse_end_file[strlen(reverse_end_file) - 6] == '.' && reverse_end_file[strlen(reverse_end_file) - 5] == 'f' && reverse_end_file[strlen(reverse_end_file) - 4] == 'a' && reverse_end_file[strlen(reverse_end_file) - 3] == 's' && reverse_end_file[strlen(reverse_end_file) - 2] == 't' && (reverse_end_file[strlen(reverse_end_file) - 1] == 'q' || reverse_end_file[strlen(reverse_end_file) - 1] == 'a'))
		{
			prefix_reverse = (char *)malloc((strlen(reverse_end_file) - 6 + 15) * sizeof(char));
			for (i = 0; i < strlen(reverse_end_file) - 6; i++)
			{
				prefix_reverse[i] = reverse_end_file[i];
			}
			for (i = strlen(reverse_end_file) - 6; i < strlen(reverse_end_file) - 6 + 15; i++)
			{
				prefix_reverse[i] = '\0';
			}
		}
		else if (reverse_end_file[strlen(reverse_end_file) - 3] == '.' && reverse_end_file[strlen(reverse_end_file) - 2] == 'f' && (reverse_end_file[strlen(reverse_end_file) - 1] == 'q' || reverse_end_file[strlen(reverse_end_file) - 1] == 'a'))
		{
			prefix_reverse = (char *)malloc((strlen(reverse_end_file) - 3 + 15) * sizeof(char));
			for (i = 0; i < strlen(reverse_end_file) - 3; i++)
			{
				prefix_reverse[i] = reverse_end_file[i];
			}
			for (i = strlen(reverse_end_file) - 3; i < strlen(reverse_end_file) - 3 + 15; i++)
			{
				prefix_reverse[i] = '\0';
			}
		}
		else
		{
			printf("Your reads don't end with .fastq, .fasta, .fa, or .fq. Please decompress your files if they are gzipped.\n");
			exit(1);
		}
		if (fasta_format == 0)
		{
			sprintf(buffer, "fastq_quality_trimmer -v -t 35 -i %s -o %s_trimmed1.fastq -Q33", reverse_end_file, prefix_reverse);
			system(buffer);
			trim_ends_fastq(prefix_reverse);
		}
		else
		{
			trim_ends_fasta(prefix_reverse);
		}
		strcat(prefix_reverse, suffix);
		strcpy(reverse_end_file, prefix_reverse);
		free(prefix_reverse);
	}
}

/**
 * @brief Build the Bowtie2 index (if missing) and run Bowtie2 to produce sam_path.
 * @param paired 1 if using paired-end reads.
 * @param fasta_format 1 if reads are FASTA (vs. FASTQ).
 * @param bowtie_reference_path Path to the reference genome to build/align against.
 * @param forward_end_file Forward reads path (paired mode).
 * @param reverse_end_file Reverse reads path (paired mode).
 * @param single_end_file Single-end reads path.
 * @param sam_path Output SAM file path.
 */
void perform_bowtie_alignment(int paired, int fasta_format, char *bowtie_reference_path, char *forward_end_file, char *reverse_end_file, char *single_end_file, char *sam_path)
{
	char *buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(buffer, '\0', FASTA_MAXLINE);

	sprintf(buffer, "bowtie2-build -f %s %s", bowtie_reference_path, bowtie_reference_path);
	system(buffer);

	if (paired == 1 && fasta_format == 1)
	{
		sprintf(buffer, "bowtie2 --all -f -x %s -1 %s -2 %s -S %s", bowtie_reference_path, forward_end_file, reverse_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 0 && fasta_format == 1)
	{
		sprintf(buffer, "bowtie2 --all -f -x %s -U %s -S %s", bowtie_reference_path, single_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 1 && fasta_format == 0)
	{
		sprintf(buffer, "bowtie2 --all -x %s -1 %s -2 %s -S %s", bowtie_reference_path, forward_end_file, reverse_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 0 && fasta_format == 0)
	{
		sprintf(buffer, "bowtie2 --all -x %s -U %s -S %s", bowtie_reference_path, single_end_file, sam_path);
		system(buffer);
	}
	free(buffer);
}

/**
 * @brief Same as perform_bowtie_alignment(), but with --xeq (extended CIGAR,
 * distinguishing '=' match from 'X' mismatch), used only for error-rate estimation.
 * @param paired 1 if using paired-end reads.
 * @param fasta_format 1 if reads are FASTA (vs. FASTQ).
 * @param bowtie_reference_path Path to the reference genome to build/align against.
 * @param forward_end_file Forward reads path (paired mode).
 * @param reverse_end_file Reverse reads path (paired mode).
 * @param single_end_file Single-end reads path.
 * @param sam_path Output SAM file path.
 */
void perform_bowtie_alignment_xeq(int paired, int fasta_format, char *bowtie_reference_path, char *forward_end_file, char *reverse_end_file, char *single_end_file, char *sam_path)
{
	char *buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(buffer, '\0', FASTA_MAXLINE);

	sprintf(buffer, "bowtie2-build -f %s %s", bowtie_reference_path, bowtie_reference_path);
	system(buffer);

	if (paired == 1 && fasta_format == 1)
	{
		sprintf(buffer, "bowtie2 --all --xeq -f -x %s -1 %s -2 %s -S %s", bowtie_reference_path, forward_end_file, reverse_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 0 && fasta_format == 1)
	{
		sprintf(buffer, "bowtie2 --all --xeq -f -x %s -U %s -S %s", bowtie_reference_path, single_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 1 && fasta_format == 0)
	{
		sprintf(buffer, "bowtie2 --all --xeq -x %s -1 %s -2 %s -S %s", bowtie_reference_path, forward_end_file, reverse_end_file, sam_path);
		system(buffer);
	}
	else if (paired == 0 && fasta_format == 0)
	{
		sprintf(buffer, "bowtie2 --all --xeq -x %s -U %s -S %s", bowtie_reference_path, single_end_file, sam_path);
		system(buffer);
	}
	free(buffer);
}

/**
 * @brief Estimate whether reads need quality/end trimming, from an --xeq SAM file.
 *
 * Counts mismatches ('X' CIGAR ops) in the first/last 10bp of each read
 * ("ends") versus the middle, and flags cleaning as needed if the per-base
 * error rate at the ends is more than 3x the middle's rate.
 *
 * @param sam_file Path to a SAM file produced with perform_bowtie_alignment_xeq().
 * @return 1 if reads should be cleaned (see clean_reads()), 0 otherwise.
 */
int calculate_error_rates(char sam_file[])
{
	FILE *sam_ptr;
	if ((sam_ptr = fopen(sam_file, "r")) == (FILE *)NULL)
		fprintf(stderr, "Sam file could not be opened.\n");
	int i, j;
	char buffer[FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	int mismatches_ends = 0;
	int mismatches_middle = 0;
	int total_aligned_reads = 0;
	int invoke_cleaning = 0;
	int total_length = 0;
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar[i] = 0;
		cigar_chars[i] = '\0';
	}
	while (fgets(buffer, FASTA_MAXLINE, sam_ptr) != NULL)
	{
		if (buffer[0] != '@')
		{
			total_aligned_reads++;
			char *buffer_copy = strdup(buffer);
			int length_of_sam = strlen(buffer);
			s = strtok(buffer, "\t");
			for (i = 0; i < 5; i++)
			{
				s = strtok(NULL, "\t");
			}
			char *cigar_string;
			cigar_string = strdup(s);
			char *copy = strdup(cigar_string);
			char *res = strtok(cigar_string, "=XID");
			int index = 0;
			while (res)
			{
				int from = res - cigar_string + strlen(res);
				int cigar_count = 0;
				sscanf(res, "%d", &cigar_count);
				res = strtok(NULL, "=XID");
				int to = res != NULL ? res - cigar_string : strlen(copy);
				char cigar_char = '\0';
				sscanf(copy + from, "%c", &cigar_char);
				cigar[index] = cigar_count;
				cigar_chars[index] = cigar_char;
				index++;
			}
			free(copy);
			free(cigar_string);
			s = strtok(buffer_copy, "\t");
			for (i = 0; i < 9; i++)
			{
				s = strtok(NULL, "\t");
			}
			int cigar_char_count = index;
			int position = 0;
			char *sequence = s;
			int length_of_seq = strlen(sequence);
			total_length = total_length + length_of_seq;
			for (i = 0; i < cigar_char_count; i++)
			{
				for (j = 0; j < cigar[i]; j++)
				{
					if (cigar_chars[i] == '=' || cigar_chars[i] == 'X' || cigar_chars[i] == 'I')
					{
						position++;
					}
					if (cigar_chars[i] == 'X' && position < 10)
					{
						mismatches_ends++;
					}
					else if (cigar_chars[i] == 'X' && position > length_of_seq - 11)
					{
						mismatches_ends++;
					}
					else if (cigar_chars[i] == 'X')
					{
						mismatches_middle++;
					}
				}
			}
			free(buffer_copy);
		}
	}
	fclose(sam_ptr);
	double average_size;
	average_size = total_length / total_aligned_reads;
	double error_rate_ends;
	double error_rate_middle;
	error_rate_ends = mismatches_ends / total_aligned_reads;
	error_rate_ends = error_rate_ends / 20;
	error_rate_middle = mismatches_middle / total_aligned_reads;
	error_rate_middle = error_rate_middle / (average_size - 20);
	printf("Error rate ends: %lf\n", error_rate_ends);
	printf("Error rate middle: %lf\n", error_rate_middle);
	if (error_rate_ends > 3 * error_rate_middle)
	{
		printf("The error rate in your 5' and 3' ends of your reads is 3x larger than the error rate in the middle of your reads. Your reads need cleaning... quality filtering and trimming ends\n");
		invoke_cleaning = 1;
	}
	return invoke_cleaning;
}

/** @brief qsort comparator for sorting an array of C strings alphabetically. */
int cmp_str(const void *a, const void *b)
{
	return strcmp(*(const char **)a, *(const char **)b);
}

/**
 * @brief Read the strain name from the header line of a single-sequence FASTA file.
 * @param fasta_path Path to a FASTA file whose first line is ">strain_name".
 * @return Newly allocated string with the header's name (no '>', no newline).
 */
char *read_fasta_header_name(char *fasta_path)
{
	FILE *file = fopen(fasta_path, "r");
	if (file == NULL)
	{
		fprintf(stderr, "Error: could not open reference file '%s' to read its name\n", fasta_path);
		exit(1);
	}
	char buffer[FASTA_MAXLINE];
	if (fgets(buffer, FASTA_MAXLINE, file) == NULL || buffer[0] != '>')
	{
		fprintf(stderr, "Error: '%s' does not start with a FASTA header line\n", fasta_path);
		exit(1);
	}
	fclose(file);
	buffer[strcspn(buffer, "\r\n")] = '\0';
	return strdup(buffer + 1); // skip the leading '>'
}

/**
 * @brief List the first `num_references` regular files in a directory, sorted
 * alphabetically, as full paths.
 *
 * Used to line up the Nth file across msa/reference/variant directories as
 * belonging to the same subtype (see main()). Exits the program with an
 * error if the directory has fewer than num_references usable files.
 *
 * @param dir_path Directory to scan.
 * @param num_references Number of files required (and returned).
 * @param dir_label Human-readable label for this directory, used in error messages.
 * @return Newly allocated array of num_references newly allocated path strings
 * (caller owns both the array and each string).
 */
char **list_sorted_dir_files(char *dir_path, int num_references, char *dir_label)
{
	DIR *dir = opendir(dir_path);
	if (dir == NULL)
	{
		fprintf(stderr, "Error: could not open %s directory '%s'\n", dir_label, dir_path);
		exit(1);
	}
 
	char **filenames = NULL;
	int count = 0;
	struct dirent *entry;
	while ((entry = readdir(dir)) != NULL)
	{
		if (strcmp(entry->d_name, ".") == 0 || strcmp(entry->d_name, "..") == 0)
		{
			continue;
		}
		filenames = (char **)realloc(filenames, (count + 1) * sizeof(char *));
		filenames[count] = strdup(entry->d_name);
		count++;
	}
	closedir(dir);
 
	if (count < num_references)
	{
		fprintf(stderr, "Error: %s directory '%s' has %d file(s), but %d reference(s) were requested (-N).\n", dir_label, dir_path, count, num_references);
		exit(1);
	}
 
	// sort alphabetically so position i means the same subtype across all
	// four directories (msa/reference/bowtie2-reference/variant)
	qsort(filenames, count, sizeof(char *), cmp_str);
 
	if (count > num_references)
	{
		fprintf(stderr, "Warning: %s directory '%s' has %d files; only using the first %d (sorted alphabetically).\n", dir_label, dir_path, count, num_references);
	}
 
	char **paths = (char **)malloc(num_references * sizeof(char *));
	int i;
	for (i = 0; i < num_references; i++)
	{
		paths[i] = (char *)malloc((strlen(dir_path) + strlen(filenames[i]) + 2) * sizeof(char));
		sprintf(paths[i], "%s/%s", dir_path, filenames[i]);
	}
 
	for (i = 0; i < count; i++)
	{
		free(filenames[i]);
	}
	free(filenames);
 
	return paths;
}

/**
 * @brief Entry point: runs the full strain-elimination + mismatch-matrix pipeline.
 *
 * Stages: parse CLI options -> (optionally) clean reads -> align reference to
 * build reference_index -> run Bowtie2 (twice: once --xeq to check error rates,
 * once normally) -> load the MSA panel -> compute allele frequencies and
 * eliminate incompatible strains -> build the read x strain mismatch matrix
 * (threaded if paired+read-sam mode) -> hand off to EM_C_LLR.R for the final
 * relative-abundance estimation.
 */
int main(int argc, char **argv)
{
	struct timespec tstart = {0, 0}, tend = {0.0};
	Options opt;
	opt.remove_identical = 0;
	opt.paired = 0;
	opt.error = 0.005;
	opt.coverage = 50;
	opt.clean_reads = 0;
	opt.fasta_format = 0;
	opt.freq = 0.01;
	opt.llr = 0;
	opt.min_strains = 500;
	opt.max_strains = 10000;
	opt.number_of_cores = 1;
	opt.no_read_bam = 0;
	opt.deletion_threshold = 0.002;
	opt.num_references = 0;
	memset(opt.print_counts, '\0', 1000);
	memset(opt.print_deletions, '\0', 1000);
	parse_options(argc, argv, &opt);
	int i, j, k, l;

	// check that the user specified a positive number of references
	if (opt.num_references <= 0)
	{
		fprintf(stderr, "Error: -N/--num-references must be set to a positive number.\n");
		exit(1);
	}
 
	char *msa_filepath = opt.msa_filepath;
	char **msa_reference_filepaths = list_sorted_dir_files(opt.msa_reference_dir, opt.num_references, "MSA reference");
	char **bowtie2_reference_filepaths = list_sorted_dir_files(opt.bowtie2_reference_dir, opt.num_references, "Bowtie2 reference");
	char **variant_sites_filepaths = list_sorted_dir_files(opt.variant_sites_dir, opt.num_references, "variant sites");
 
	if (opt.clean_reads == 1)
	{
		printf("You've selected -d to clean your FASTA/FASTQ reads. If this is not correct, please quit the program and removed the -d option. Cleaning reads...\n");
		clean_reads(opt.paired, opt.fasta_format, opt.single_end_filepath, opt.forward_end_filepath, opt.reverse_end_filepath);
	}
	
	// no problematic sites for now
	// TODO: implement the process_problematic_sites() function to read in problematic sites from a file (similar to variant sites)
	int problematic_sites[] = {};
	int number_of_problematic_sites = 0;

	int num_sam_lines = 0;

	ReferenceData *ref_data_str = (ReferenceData *)malloc(opt.num_references * sizeof(ReferenceData));

	int ref_idx;
	for (ref_idx = 0; ref_idx < opt.num_references; ref_idx++)
	{
		align_references(&ref_data_str[ref_idx], msa_reference_filepaths[ref_idx], bowtie2_reference_filepaths[ref_idx]);

		char sam_path[1100];
		sprintf(sam_path, "%s.%d", opt.sam, ref_idx);

		perform_bowtie_alignment_xeq(opt.paired, opt.fasta_format, bowtie2_reference_paths[ref_idx], opt.forward_end_file, opt.reverse_end_file, opt.single_end_file, sam_path);
		int invoke_cleaning = calculate_error_rates(sam_path);
		if (invoke_cleaning == 1 && opt.clean_reads == 0)
		{
			clean_reads(opt.paired, opt.fasta_format, opt.single_end_file, opt.forward_end_file, opt.reverse_end_file);
		}
		perform_bowtie_alignment(opt.paired, opt.fasta_format, bowtie2_reference_paths[ref_idx], opt.forward_end_file, opt.reverse_end_file, opt.single_end_file, sam_path);

		gzFile MSA_file = Z_NULL;
		if ((MSA_file = gzopen(MSA_paths[ref_idx], "r")) == Z_NULL)
			fprintf(stderr, "MSA index %d File could not be opened.\n", ref_idx);
		int msa_seq_length = parse_msa_seq_length(MSA_file);
		gzclose(MSA_file);
		printf("Length of MSA index %d: %d\n", ref_idx, length_of_MSA);

		// TODO: implement the process_problematic_sites() function to read in problematic sites from a file (similar to variant sites)
		/*int* problematic_sites = (int*)malloc(1000*sizeof(int));
		for(i=0; i<1000; i++){
			problematic_sites[i]=-1;
		}
		int number_of_problematic_sites=process_problematic_sites(problematic_sites);*/
		
		if ((MSA_file = gzopen(MSA_paths[ref_idx], "r")) == Z_NULL)
			fprintf(stderr, "MSA index %d File could not be opened.\n", ref_idx);
		int num_strains = 0;
		int max_strain_name_length = 0;
		setNumStrains(MSA_file, &num_strains, &max_strain_name_length);
		gzclose(MSA_file);
		printf("Number of strains for reference %d: %d\n", ref_idx, num_strains);
		printf("Maxname length for reference %d: %d\n", ref_idx, max_strain_name_length);

		char **MSA = (char **)malloc(num_strains * sizeof(char *));
		for (i = 0; i < num_strains; i++)
		{
			MSA[i] = (char *)malloc((length_of_MSA + 1) * sizeof(char));
		}
		char **names_of_strains = (char **)malloc(num_strains * sizeof(char *));
		for (i = 0; i < num_strains; i++)
		{
			names_of_strains[i] = (char *)malloc((max_strain_name_length + 1) * sizeof(char));
		}

		printf("Reading in MSA index %d...\n", ref_idx);
		clock_gettime(CLOCK_MONOTONIC, &tstart);
		if ((MSA_file = gzopen(MSA_paths[ref_idx], "r")) == Z_NULL)
			fprintf(stderr, "MSA index %d File could not be opened.\n", ref_idx);
		char *reference_name = read_fasta_header_name(MSA_reference_paths[ref_idx]);
		readInMSA(MSA_file, MSA, names_of_strains, length_of_MSA, reference_name);
		free(reference_name);
		gzclose(MSA_file);
		clock_gettime(CLOCK_MONOTONIC, &tend);
		printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));

		// TODO: eventually implement the removal of identical strains for speedup/optimization
		/*int* identical = (int *)malloc(number_of_strains*sizeof(int));
		int number_of_identical_strains=number_of_strains;
		if (opt.remove_identical==0){
			for(i=0; i<number_of_strains; i++){
				identical[i]=i;
			}
		}*/
		// if (opt.remove_identical==1){
		//	printf("Finding identical sequences\n");
		//	clock_gettime(CLOCK_MONOTONIC, &tstart);
		//	for(i=0; i<number_of_strains; i++){
		//		identical[i]=-1;
		//	}
		//	number_of_identical_strains=removeIdenticalStrains(number_of_strains,length_of_MSA,MSA,identical,names_of_strains,max_name_length);
		//	clock_gettime(CLOCK_MONOTONIC, &tend);
		//	printf("Took %.5fsec\n",((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec));
		// }

		FILE *variant_sites_file;
		int *variant_sites;
		int *number_of_variant_sites_p = (int *)malloc(sizeof(int));
		if ((variant_sites_file = fopen(variant_paths[ref_idx], "r")) == (FILE *)NULL)
			fprintf(stderr, "variants index %d File could not be opened.\n", ref_idx);
		variant_sites = readVariantSitesFile(variant_sites_file, number_of_variant_sites_p);
		fclose(variant_sites_file);
		int this_number_of_variant_sites = number_of_variant_sites_p[0];
		free(number_of_variant_sites_p);
		printf("Number of variant sites for reference %d: %d\n", ref_idx, this_number_of_variant_sites);

		printf("Calculating allele frequencies for reference %d...\n", ref_idx);
		FILE *sam_file;
		if ((sam_file = fopen(sam_path, "r")) == (FILE *)NULL)
			fprintf(stderr, "sam index %d File could not be opened.\n", ref_idx);
		double **allele_frequency = (double **)malloc(length_of_MSA * sizeof(double *));
		for (i = 0; i < length_of_MSA; i++)
		{
			allele_frequency[i] = (double *)malloc(4 * sizeof(double));
			for (j = 0; j < 4; j++)
			{
				allele_frequency[i][j] = 0;
			}
		}
		int number_of_strains_remaining = 0;

		char **sam_results = NULL;
		int cached_num_sam_lines = 0;
		int cached_max_sam_line_length = -1;

		if (opt.paired == 1)
		{
			number_of_strains_remaining = calculateAlleleFreq_paired(sam_file, allele_frequency, length_of_MSA, MSA, num_strains, names_of_strains, opt.freq, max_strain_name_length, tstart, tend, this_number_of_variant_sites, variant_sites, opt.coverage, reference_index, opt.min_strains, opt.max_strains, opt.print_counts, opt.print_deletions, opt.deletion_threshold, &sam_results, &num_sam_lines, &max_sam_line_length);
		}
		else
		{
			number_of_strains_remaining = calculateAlleleFreq(sam_file, allele_frequency, length_of_MSA, MSA, num_strains, names_of_strains, opt.freq, max_strain_name_length, tstart, tend, this_number_of_variant_sites, variant_sites, opt.coverage, reference_index, opt.min_strains, opt.max_strains, opt.print_counts, opt.print_deletions, opt.deletion_threshold, &sam_results, &num_sam_lines, &max_sam_line_length);
		}

		fclose(sam_file);
		for (i = 0; i < length_of_MSA; i++)
		{
			free(allele_frequency[i]);
		}
		free(allele_frequency);
		// free(identical);

		int *strains_kept = (int *)malloc(number_of_strains_remaining * sizeof(int));
		int placement = 0;
		for (i = 0; i < num_strains; i++)
		{
			if (names_of_strains[i][0] != '\0')
			{
				strains_kept[placement] = i;
				placement++;
			}
		}

		// FIXME: Uses global variables still, fix this
		resize_names_of_strains = (char **)malloc(number_of_strains_remaining * sizeof(char *));
		resize_MSA = (char **)malloc(number_of_strains_remaining * sizeof(char *));
		for (i = 0; i < number_of_strains_remaining; i++)
		{
			resize_names_of_strains[i] = (char *)malloc((max_strain_name_length + 1) * sizeof(char));
			resize_MSA[i] = (char *)malloc((length_of_MSA + 1) * sizeof(char));
			for (j = 0; j < length_of_MSA + 1; j++)
			{
				resize_MSA[i][j] = '\0';
			}
			for (j = 0; j < max_strain_name_length + 1; j++)
			{
				resize_names_of_strains[i][j] = '\0';
			}
		}
		reallocate_memory(number_of_strains_remaining, strains_kept, MSA, names_of_strains, max_strain_name_length, num_strains);
		free(strains_kept);
 
		if (number_of_strains_remaining == 0)
		{
			printf("No strains remaining for reference %d, exiting...\n", ref_idx);
			exit(1);
		}
 
		reference_data[ref_idx].reference_index = reference_index;
		reference_data[ref_idx].resize_MSA = resize_MSA;
		reference_data[ref_idx].resize_names_of_strains = resize_names_of_strains;
		reference_data[ref_idx].number_of_strains_remaining = number_of_strains_remaining;
		reference_data[ref_idx].length_of_MSA = length_of_MSA;
		reference_data[ref_idx].sam_results = sam_results;
		
	}

	printf("Creating mismatch matrix...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	FILE *outfile;
	if ((outfile = fopen(opt.outfile, "w")) == (FILE *)NULL)
		fprintf(stderr, "File could not be opened.\n");
	if (opt.paired == 1 && opt.no_read_bam == 0)
	{
		pthread_t threads[opt.number_of_cores];	 // array of our threads
		MismatchMatrixThreadStruct thread_str[opt.number_of_cores]; // array of stuct that contains input and output for each thread
		for (i = 0; i < opt.number_of_cores; i++)
		{
			thread_str[i].sam_line_start = 0;
			thread_str[i].sam_line_end = 0;
			thread_str[i].thread_index = i;
			thread_str[i].mismatch_matrix_row_partition = malloc((opt.number_of_cores + 1) * sizeof(int));
		}
		// Split sam_results into opt.number_of_cores roughly-equal chunks, then
		// nudge each chunk's boundary (via adjust_start_end) so no thread's
		// range splits a read pair across two threads.
		int divideFile, start, end;
		divideFile = max_sam_length[1] / opt.number_of_cores;
		j = 0;
		for (i = 0; i < opt.number_of_cores; i++)
		{
			start = j;
			end = j + divideFile;
			if (i == opt.number_of_cores - 1)
			{
				end = max_sam_length[1] - 1;
			}
			int *adjust_ends = (int *)malloc(2 * sizeof(int));
			if (thread_str[i].sam_line_start == 1)
			{
				start++;
			}
			adjust_ends[0] = 0;
			adjust_ends[1] = 0;
			adjust_start_end(start, end, adjust_ends, max_sam_length[0]);
			if (start != adjust_ends[0])
			{
				thread_str[i - 1].sam_line_end = thread_str[i - 1].sam_line_end - 1;
			}
			if (end != adjust_ends[1])
			{
				thread_str[i + 1].sam_line_start = thread_str[i + 1].sam_line_start + 1;
			}
			if (i == opt.number_of_cores - 1)
			{
				adjust_ends[1] = max_sam_length[1];
			}
			thread_str[i].sam_line_start = adjust_ends[0];
			thread_str[i].sam_line_end = adjust_ends[1];
			thread_str[i].thread_index = i;
			thread_str[i].max_sam_line_length = max_sam_length[0];
			j = j + divideFile;
			free(adjust_ends);
		}
		for (i = 0; i < opt.number_of_cores; i++)
		{
			printf("thread: %d start: %d end: %d\n", thread_str[i].thread_index, thread_str[i].sam_line_start, thread_str[i].sam_line_end);
			thread_str[i].results_str->mismatch = (char **)malloc((thread_str[i].sam_line_end - thread_str[i].sam_line_start) * (sizeof(char *)));
			for (k = 0; k < thread_str[i].sam_line_end - thread_str[i].sam_line_start; k++)
			{
				thread_str[i].results_str->mismatch[k] = (char *)malloc((max_sam_length[0] + 300000) * sizeof(char));
				for (l = 0; l < max_sam_length[0] + 300000; l++)
				{
					thread_str[i].results_str->mismatch[k][l] = '\0';
				}
			}
		}
		for (i = 0; i < opt.number_of_cores; i++)
		{
			pthread_create(&threads[i], NULL, writeMismatchMatrix_paired, &thread_str[i]);
		}
		for (i = 0; i < opt.number_of_cores; i++)
		{
			pthread_join(threads[i], NULL);
		}
		for (i = 0; i < opt.number_of_cores; i++)
		{
			for (j = 0; j < (thread_str[i].sam_line_end - thread_str[i].sam_line_start); j++)
			{
				if (thread_str[i].results_str->mismatch[j][0] == '\0')
				{
					break;
				}
				fprintf(outfile, "%s\n", thread_str[i].results_str->mismatch[j]);
			}
		}
		for (i = 0; i < max_sam_length[1]; i++)
		{
			free(sam_results[i]);
		}
		free(sam_results);
		for (i = 0; i < opt.number_of_cores; i++)
		{
			for (j = 0; j < thread_str[i].sam_line_end - thread_str[i].sam_line_start; j++)
			{
				free(thread_str[i].results_str->mismatch[j]);
			}
			free(thread_str[i].results_str->mismatch);
			free(thread_str[i].results_str);
		}
	}
	else if (opt.paired == 1 && opt.no_read_bam == 1)
	{
		// Paired-end, low-memory mode: single-threaded, streams the SAM file directly.
		writeMismatchMatrix_paired_no_read_bam(outfile, sam_file, resize_MSA, length_of_MSA, number_of_strains, number_of_strains_remaining, resize_names_of_strains, reference_index);
	}
	else
	{
		// Single-end mode: no threading, no mate-pair overlap handling needed.
		writeMismatchMatrix(outfile, sam_file, MSA, strains_kept, length_of_MSA, number_of_strains, number_of_strains_remaining);
	}
	fclose(outfile);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fseconds\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	// free(reference);
	// free(variant_sites);
	// for(i=0; i<number_of_strains; i++){
	// free(names_of_strains[i]);
	//	free(MSA[i]);
	//}
	// free(names_of_strains);
	// free(MSA);
	for (i = 0; i < number_of_strains_remaining; i++)
	{
		free(resize_names_of_strains[i]);
		free(resize_MSA[i]);
	}
	free(resize_names_of_strains);
	free(resize_MSA);
	if (number_of_strains_remaining > 0)
	{
		free(strains_kept);
	}
	else
	{
		printf("No strains remaining exiting...\n");
		exit(1);
	}
	char *buffer = (char *)malloc(FASTA_MAXLINE * sizeof(char));
	memset(buffer, '\0', FASTA_MAXLINE);
	// FIXME: Remove the hardcoding here 
	if (opt.llr == 1)
	{
		sprintf(buffer, "Rscript ../EM_C_LLR.R -i %s -f %lf -e %lf -l -s -v %s -r %s -b %s", opt.outfile, opt.freq, opt.error, variant_paths[0], MSA_paths[0], opt.print_counts);
	}
	else
	{
		sprintf(buffer, "Rscript ../EM_C_LLR.R -i %s -f %lf -e %lf -s -v %s -r %s -b %s", opt.outfile, opt.freq, opt.error, variant_paths[0], MSA_paths[0], opt.print_counts);
	}
	system(buffer);
	free(buffer);
}