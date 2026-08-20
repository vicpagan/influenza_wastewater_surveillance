#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "calculate_proportions.h"
#include "external/covid_em.h"
#include "external/hashmap.h"


void calculate_proportions(MismatchData *mismatch_data_str, char *output_csv_filepath, double error_rate, double filter, int compute_strain_llr, int compute_site_llr, int num_plot, int num_threads)
{
    int i, read_idx, msa_seq_idx;
	int num_reads = mismatch_data_str->num_reads;
    if (num_reads <= 0)
    {
        printf("Error: No reads in mismatch matrix!\n");
        exit(1);
    }

    int num_msa_sequences = mismatch_data_str->num_msa_sequences;

    HASHMAP(char, HashmapEntry) hashmap;
    hashmap_init(&hashmap, hashmap_hash_string, strcmp);
    
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        char *current_column = (char *)calloc((num_reads * 12) + 1, sizeof(char));

        int offset = 0;
        int length;
        for (read_idx = 0; read_idx < num_reads; read_idx++)
        {
            length = sprintf(current_column + offset, "%d,", mismatch_data_str->mismatch_matrix[read_idx][msa_seq_idx]);
		    offset += length;
        }

        HashmapEntry *entry = hashmap_get(&hashmap, current_column);
        if (entry == NULL)
        {
            HashmapEntry *new_entry = (HashmapEntry *)malloc(sizeof(HashmapEntry));

            int new_entry_mismatch_column_length = strlen(current_column);
            new_entry->mismatch_column = (char *)malloc(new_entry_mismatch_column_length + 1);
            strcpy(new_entry->mismatch_column, current_column);

            int new_entry_strain_name_length = strlen(mismatch_data_str->msa_sequence_names[msa_seq_idx]);
            new_entry->msa_strain_names_length = new_entry_strain_name_length;
            new_entry->msa_strain_names = (char *)malloc(new_entry_strain_name_length + 1);
            strcpy(new_entry->msa_strain_names, mismatch_data_str->msa_sequence_names[msa_seq_idx]);

            hashmap_put(&hashmap, new_entry->mismatch_column, new_entry);
        }
        else
        {
            int duplicate_strain_name_length = strlen(mismatch_data_str->msa_sequence_names[msa_seq_idx]);
            entry->msa_strain_names = (char *)realloc(entry->msa_strain_names, entry->msa_strain_names_length + duplicate_strain_name_length + 2);
            strcat(entry->msa_strain_names, ":");
            strcat(entry->msa_strain_names, mismatch_data_str->msa_sequence_names[msa_seq_idx]);
            entry->msa_strain_names_length = entry->msa_strain_names_length + duplicate_strain_name_length + 1;
        }
        free(current_column);
    }

    char **strain_names_to_group = (char **)malloc(num_msa_sequences * sizeof(char *));
    char **strain_names_to_remove = (char **)malloc(num_msa_sequences * sizeof(char *));
    int strain_names_to_remove_idx = 0;
    int strain_names_to_group_idx = 0;

    char *mismatch_column;
    HashmapEntry *entry;
    void *pos;
    hashmap_foreach_safe(mismatch_column, entry, &hashmap, pos)
    {
        char *delimiter = strchr(entry->msa_strain_names, ':');
        if (delimiter != NULL)
        {
            int strain_name_length = delimiter - entry->msa_strain_names;
            strain_names_to_group[strain_names_to_group_idx] = (char *)malloc(strain_name_length + 1);
            memcpy(strain_names_to_group[strain_names_to_group_idx], entry->msa_strain_names, strain_name_length);
            strain_names_to_group[strain_names_to_group_idx][strain_name_length] = '\0';
            strain_names_to_group_idx++;

            do
            {
                char *next_delimiter = strchr(delimiter + 1, ':');
                if (next_delimiter == NULL)
                {
                    strain_name_length = strlen(delimiter + 1);
                }
                else
                {
                    strain_name_length = next_delimiter - (delimiter + 1);
                }
                strain_names_to_remove[strain_names_to_remove_idx] = (char *)malloc(strain_name_length + 1);
                memcpy(strain_names_to_remove[strain_names_to_remove_idx], delimiter + 1, strain_name_length);
                strain_names_to_remove[strain_names_to_remove_idx][strain_name_length] = '\0';
                strain_names_to_remove_idx++;

                delimiter = next_delimiter;
            }
            while(delimiter != NULL);
        }

        hashmap_remove(&hashmap, entry->mismatch_column);
        free(entry->mismatch_column);
        free(entry->msa_strain_names);
        free(entry);
    }
    hashmap_cleanup(&hashmap);

    int *columns_to_remove = (int *)calloc(num_msa_sequences, sizeof(int));
    int current_group_idx = 1;
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        char *current_msa_seq_name = mismatch_data_str->msa_sequence_names[msa_seq_idx];
        int found_sequence = 0;
        int strain_names_to_remove_counter = 0;
        while (found_sequence == 0 && strain_names_to_remove_counter < strain_names_to_remove_idx)
        {
            if (strcmp(current_msa_seq_name, strain_names_to_remove[strain_names_to_remove_counter]) == 0)
            {
                columns_to_remove[msa_seq_idx] = 1;
                found_sequence = 1;
            }
            strain_names_to_remove_counter++;
        }

        found_sequence = 0;
        int strain_names_to_group_counter = 0;
        while (found_sequence == 0 && strain_names_to_group_counter < strain_names_to_group_idx)
        {
            if (strcmp(current_msa_seq_name, strain_names_to_group[strain_names_to_group_counter]) == 0)
            {
                int num_digits = (int)log10(abs(current_group_idx)) + 1;
                current_msa_seq_name = (char *)realloc(current_msa_seq_name, num_digits + 7);

                sprintf(current_msa_seq_name, "Group %d", current_group_idx);
                mismatch_data_str->msa_sequence_names[msa_seq_idx] = current_msa_seq_name;

                current_group_idx++;
                found_sequence = 1;
            }
            strain_names_to_group_counter++;
        }
    }

    for (i = 0; i < strain_names_to_group_idx; i++) 
    {
        free(strain_names_to_group[i]);
    }
    free(strain_names_to_group);
    for (i = 0; i < strain_names_to_remove_idx; i++) 
    {
        free(strain_names_to_remove[i]);
    }   
    free(strain_names_to_remove);

    int new_num_msa_sequences = 0;
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        if (columns_to_remove[msa_seq_idx] == 0)
        {
            mismatch_data_str->msa_sequence_names[new_num_msa_sequences] = mismatch_data_str->msa_sequence_names[msa_seq_idx];
            new_num_msa_sequences++;
        }
        else
        {
            free(mismatch_data_str->msa_sequence_names[msa_seq_idx]);
        }
    }
    mismatch_data_str->msa_sequence_names = (char **)realloc(mismatch_data_str->msa_sequence_names, new_num_msa_sequences * sizeof(char *));

    for (read_idx = 0; read_idx < num_reads; read_idx++)
    {
        int new_msa_seq_idx = 0;

        for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
        {
            if (columns_to_remove[msa_seq_idx] == 0)
            {
                mismatch_data_str->mismatch_matrix[read_idx][new_msa_seq_idx] = mismatch_data_str->mismatch_matrix[read_idx][msa_seq_idx];
                new_msa_seq_idx++;
            }
        }

        if (new_msa_seq_idx != new_num_msa_sequences)
        {
            fprintf(stderr, "Error: Row length for row %d is not as expected!", read_idx);
            exit(1);
        }

        mismatch_data_str->mismatch_matrix[read_idx] = (int *)realloc(mismatch_data_str->mismatch_matrix[read_idx], new_msa_seq_idx * sizeof(int));
    }

    mismatch_data_str->num_msa_sequences = new_num_msa_sequences;
    num_msa_sequences = new_num_msa_sequences;
    free(columns_to_remove);

}