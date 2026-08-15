#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "msa.h"

/**
 * @brief Parse the MSA file for metadata information
 * 
 * @param msa_file The MSA file to read from
 * @param msa Reference to the stored MSA instance
 */
void parse_msa_info(gzFile msa_file, MSA *msa_str)
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

	msa_str->num_sequences = num_sequences;
	msa_str->sequence_length = sequence_length;
	msa_str->max_sequence_name_length = max_sequence_name_length;
}

/**
 * @brief Reads in the MSA sequences and sequence names from the input file to the given MSA struct
 * 
 * @param msa_file The MSA file to read from
 * @param msa The reference to the stored MSA instance
 * @param ref_seq_name The name of a reference sequence
 */
void read_msa_sequences(gzFile msa_file, MSA *msa_str)
{
	char buffer[FASTA_MAXLINE];

	int seq_idx = -1;

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
		}
		else
		{
			// read in sequence line
			buffer[strcspn(buffer, "\r\n")] = '\0';

			int length = strlen(buffer);
			for (i = 0; i < length; i++)
			{
				switch (buffer[i])
                {
                    case 'A':
                    case 'a':
                        msa_str->sequences[seq_idx][i] = 'A';
                        break;

                    case 'C':
                    case 'c':
                        msa_str->sequences[seq_idx][i] = 'C';
                        break;

                    case 'G':
                    case 'g':
                        msa_str->sequences[seq_idx][i] = 'G';
                        break;

                    case 'T':
                    case 't':
                        msa_str->sequences[seq_idx][i] = 'T';
                        break;

                    case '-':
                        msa_str->sequences[seq_idx][i] = '-';
                        break;

                    default:
                        msa_str->sequences[seq_idx][i] = '\0';
                        break;
                }
			}
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
MSA read_in_msa(char *msa_filepath)
{
	MSA msa_str;
	
	gzFile msa_file;
	if ((msa_file = gzopen(msa_filepath, "r")) == (gzFile)NULL)
	{
		fprintf(stderr, "MSA File could not be opened.\n");

		msa_str.num_sequences = -1;
		msa_str.sequence_length = -1;
		msa_str.max_sequence_name_length = -1;

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
			msa_str.sequence_names[i] = calloc((msa_str.max_sequence_name_length + 1), sizeof(char));
			msa_str.sequences[i] = calloc((msa_str.sequence_length + 1), sizeof(char));
		}

		gzrewind(msa_file);

		// parse and store MSA sequences and their names
		read_msa_sequences(msa_file, &msa_str);

		gzclose(msa_file);
	}
    return msa_str;
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

	int l = 0;
	int r;
	for (r = 0; r < num_sequences; r++)
	{
		if (sequences_to_remove[r] == 1)
		{
			free(msa_str->sequences[r]);
			free(msa_str->sequence_names[r]);
		}
		else
		{
			if (l != r)
			{
				msa_str->sequences[l] = msa_str->sequences[r];
				msa_str->sequence_names[l] = msa_str->sequence_names[r];
			}
			l++;
		}
	}

	msa_str->sequences = realloc(msa_str->sequences, l * sizeof(char *));
	msa_str->sequence_names = realloc(msa_str->sequence_names, l * sizeof(char *));
	msa_str->num_sequences = l;
}

/**
 * @brief Finds identical sequences in the MSA and marks them for removal
 * 
 * @param msa_str The reference to the stored MSA instance
 */
void remove_identical_sequences(MSA *msa_str)
{
	int num_sequences = msa_str->num_sequences;
	int sequence_length = msa_str->sequence_length;
	int *sequences_to_remove = (int *)calloc(num_sequences, sizeof(int));

	int num_sequences_removed = 0;
	
	// compare every sequence to each other and mark rightmost duplicates for removal
	int i, j, k, nm;
	for (i = 0; i < num_sequences; i++)
	{
		if (sequences_to_remove[i] == 0)
		{
			for (j = i + 1; j < num_sequences; j++)
			{
				if (sequences_to_remove[j] == 0)
				{
					nm = 0;
					for (k = 0; k < sequence_length; k++)
					{
						if (msa_str->sequences[i][k] != msa_str->sequences[j][k])
						{
							nm++;
						}
					}

					if (nm == 0)
					{
						printf("Removing sequence %d\n", j);
						sequences_to_remove[j] = 1;
						num_sequences_removed++;
					}
				}
			}
		}	
	}

	if (num_sequences_removed != 0)
	{
		prune_msa_sequences(msa_str, sequences_to_remove);
	}
	free(sequences_to_remove);
}
