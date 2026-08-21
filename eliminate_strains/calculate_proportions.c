#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include "calculate_proportions.h"
#include "external/covid_em.h"
#include "external/hashmap.h"

static int parameters_valid(const double *p, int n)
{
    double sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        if (!isfinite(p[i]))
        {
            return 0;
        }

        if (p[i] < 0.0)
        {
            return 0;
        }

        sum += p[i];
    }

    if (!isfinite(sum))
    {
        return 0;
    }

    if (fabs(sum - 1.0) > 1.0e-9)
    {
        return 0;
    }

    return 1;
}

static int vector_is_finite(const double *x, int n)
{
    for (int i = 0; i < n; ++i)
    {
        if (!isfinite(x[i]))
        {
            return 0;
        }
    }

    return 1;
}

double *run_squarem(const double *p0, const double **l_matrix, int num_reads, int num_msa_sequences)
{
    int msa_seq_idx;

    const double step_min0 = 1.0;
    const double step_max0 = 1.0;
    const double mstep = 4.0;
    const double objfn_inc = 1.0;
    const double tolerance = 1.0e-7;
    const double max_time = 14400;

    double step_min = step_min0;
    double step_max = step_max0;

    double *p = (double *)malloc(num_msa_sequences * sizeof(double));
    double *p1 = (double *)malloc(num_msa_sequences * sizeof(double));
    double *p2 = (double *)malloc(num_msa_sequences * sizeof(double));

    double *q1 = (double *)malloc(num_msa_sequences * sizeof(double));
    double *q2 = (double *)malloc(num_msa_sequences * sizeof(double));
    double *v = (double *)malloc(num_msa_sequences * sizeof(double));

    double *p_new = (double *)malloc(num_msa_sequences * sizeof(double));
    double *p_temp = (double *)malloc(num_msa_sequences * sizeof(double));
    double *result = (double *)malloc(num_msa_sequences * sizeof(double));

    memcpy(p, p0, num_msa_sequences * sizeof(double));

    double prev_obj = negative_log_likelihood(p, l_matrix, num_reads, num_msa_sequences);

    struct timespec start_time;
    struct timespec current_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    int converged = 0;
    double seconds_elapsed;

    do
    {
        em_update(p, l_matrix, num_reads, num_msa_sequences, p1);
        em_update(p1, l_matrix, num_reads, num_msa_sequences, p2);

        double q1_norm_sq = 0.0;
        double v_norm_sq = 0.0;
        for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
        {
            q1[msa_seq_idx] = p1[msa_seq_idx] - p[msa_seq_idx];
            q2[msa_seq_idx] = p2[msa_seq_idx] - p1[msa_seq_idx];
            v[msa_seq_idx] = q2[msa_seq_idx] - q1[msa_seq_idx];

            q1_norm_sq += q1[msa_seq_idx] * q1[msa_seq_idx];
            v_norm_sq += v[msa_seq_idx] * v[msa_seq_idx];
        }

        double alpha = 1.0;
        if (v_norm_sq <= 0.0)
        {
            memcpy(p_new, p2, num_msa_sequences * sizeof(double));
        }
        else
        {
            alpha = sqrt(q1_norm_sq / v_norm_sq);

            if (alpha < step_min)
            {
                alpha = step_min;
            }
            if (alpha > step_max)
            {
                alpha = step_max;
            }

            for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
            {
                p_new[msa_seq_idx] = p[msa_seq_idx] + 2.0 * alpha * q1[msa_seq_idx] + alpha * alpha * v[msa_seq_idx];
            }

            if (fabs(alpha - 1.0) > 0.01)
            {
                em_update(p_new, l_matrix, num_reads, num_msa_sequences, p_temp);

                memcpy(p_new, p_temp, num_msa_sequences * sizeof(double));
            }
        }

        double curr_obj = negative_log_likelihood(p_new, l_matrix, num_reads, num_msa_sequences);
        if (curr_obj > prev_obj + objfn_inc)
        {
            memcpy(p_new, p2, num_msa_sequences * sizeof(double));

            curr_obj = negative_log_likelihood(p_new, l_matrix, num_reads, num_msa_sequences);

            if (fabs(alpha - step_max) < 1.0e-12)
            {
                step_max = fmax(step_max0, step_max / mstep);
            }

            alpha = 1.0;
        }

        if (fabs(alpha - step_max) < 1.0e-12)
        {
            step_max *= mstep;
        }

        double parameter_difference_sq = 0.0;

        for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
        {
            double difference = p_new[msa_seq_idx] - p[msa_seq_idx];
            parameter_difference_sq += difference * difference;
        }

        double parameter_difference = sqrt(parameter_difference_sq);

        memcpy(p, p_new, num_msa_sequences * sizeof(double));

        prev_obj = curr_obj;

        if (parameter_difference < tolerance)
        {
            converged = 1;
        }

        clock_gettime(CLOCK_MONOTONIC, &current_time);
        seconds_elapsed = (double)(current_time.tv_sec - start_time.tv_sec) + (double)(current_time.tv_nsec - start_time.tv_nsec) * 1.0e-9;
    }
    while (seconds_elapsed < max_time && converged == 0);

    memcpy(result, p, num_msa_sequences * sizeof(double));

    free(p);
    free(p1);
    free(p2);
    free(q1);
    free(q2);
    free(v);
    free(p_new);
    free(p_temp);

    return result;
}

void calculate_proportions(MismatchData *mismatch_data_str, char *output_csv_filepath, double error_rate, double filter, int compute_strain_llr, int compute_site_llr, int num_plot, int num_threads)
{
    srand((unsigned int)time(NULL));

    int i, read_idx, msa_seq_idx;
	int num_reads = mismatch_data_str->num_reads;
    if (num_reads <= 0)
    {
        printf("Error: No reads in mismatch matrix!\n");
        exit(1);
    }

    // FIXME: Should this check be made in the beginning? Yknow so people dont waste all their time with the rest of the code just for it to fail here.
    if (error_rate < 0.001 || error_rate > 1)
    {
        printf("Error: EM error rate should be in the range (0.001, 1)!\n");
        exit(1);
    }

    int num_msa_sequences = mismatch_data_str->num_msa_sequences;

    char **group_strains = (char **)malloc(num_msa_sequences * sizeof(char *)); 

    /////////////////////////////////////////////////////////////////////
    ////////////////// REMOVING UNIDENTIFIABLE STRAINS //////////////////
    /////////////////////////////////////////////////////////////////////

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
    int strain_names_to_group_idx = 0;
    int strain_names_to_remove_idx = 0;

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

            group_strains[strain_names_to_group_idx] = strdup(entry->msa_strain_names);

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

    strain_names_to_group = (char **)realloc(strain_names_to_group, strain_names_to_group_idx * sizeof(char *));
    strain_names_to_remove = (char **)realloc(strain_names_to_remove, strain_names_to_remove_idx * sizeof(char *));
    group_strains = (char **)realloc(group_strains, strain_names_to_group_idx * sizeof(char *));

    int *columns_to_remove = (int *)calloc(num_msa_sequences, sizeof(int));
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
                int num_digits = (int)log10(strain_names_to_group_counter + 1) + 1;
                current_msa_seq_name = (char *)realloc(current_msa_seq_name, num_digits + 7);

                sprintf(current_msa_seq_name, "Group %d", strain_names_to_group_counter);
                mismatch_data_str->msa_sequence_names[msa_seq_idx] = current_msa_seq_name;

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

    /////////////////////////////////////////////////////////////////////
    /////////////////// CALCULATING LIKELIHOOD MATRIX ///////////////////
    /////////////////////////////////////////////////////////////////////

    double **likelihood_matrix = (double **)malloc(num_reads * sizeof(double *));
    
    for (read_idx = 0; read_idx < num_reads; read_idx++)
    {
        likelihood_matrix[read_idx] = (double *)malloc(num_msa_sequences * sizeof(double));

        for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
        {
            int mismatches = mismatch_data_str->mismatch_matrix[read_idx][msa_seq_idx];
            int matches = mismatch_data_str->block_sizes[read_idx] - mismatches;

            likelihood_matrix[read_idx][msa_seq_idx] = pow(error_rate, mismatches) * pow(1.0 - error_rate, matches);
            if (likelihood_matrix[read_idx][msa_seq_idx] == 0.0)
            {
                likelihood_matrix[read_idx][msa_seq_idx] = DBL_MIN;
            }
        }
        free(mismatch_data_str->mismatch_matrix[read_idx]);
    }
    free(mismatch_data_str->mismatch_matrix);

    ////////////////////////////////////////////////////////////////////
    /////////// INITIALIZING PROPORTIONS FOR EACH MSA STRAIN ///////////
    ////////////////////////////////////////////////////////////////////

    double *p0 = (double *)malloc(num_msa_sequences * sizeof(double));
    double *p0_rand = (double *)malloc(num_msa_sequences * sizeof(double));

    double p0_sum = 0.0;
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        double column_sum = 0.0;
        for (read_idx = 0; read_idx < num_reads; read_idx++)
        {
            column_sum += likelihood_matrix[read_idx][msa_seq_idx];
        }
        p0[msa_seq_idx] = column_sum;
        p0_sum += column_sum;
    }
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        p0[msa_seq_idx] /= p0_sum;
    }

    double p0_rand_sum = 0.0;
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        p0_rand[msa_seq_idx] = (double)rand() / ((double)RAND_MAX + 1.0);
        p0_rand_sum += p0_rand[msa_seq_idx];
    }
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        p0_rand[msa_seq_idx] /= p0_rand_sum;
    }

    ////////////////////////////////////////////////////////////////////
    /////////////////////////// SQUAREM RUNS ///////////////////////////
    ////////////////////////////////////////////////////////////////////

    double *p = (double *)malloc(num_msa_sequences * sizeof(double));
    double *p1 = (double *)malloc(num_msa_sequences * sizeof(double));
    double *p2 = (double *)malloc(num_msa_sequences * sizeof(double));

    double *q1 = (double *)malloc(num_msa_sequences * sizeof(double));
    double *q2 = (double *)malloc(num_msa_sequences * sizeof(double));
    double *v = (double *)malloc(num_msa_sequences * sizeof(double));

    double *p_new = (double *)malloc(num_msa_sequences * sizeof(double));
    double *p_temp = (double *)malloc(num_msa_sequences * sizeof(double));
    double *result = (double *)malloc(num_msa_sequences * sizeof(double));

    prinf("Sleeping for 30 seconds to check peak memory...\n");
    sleep(30);

    free(p0);
    free(p0_rand);

    free(p);
    free(p1);
    free(p2);
    free(q1);
    free(q2);
    free(v);
    free(p_new);
    free(p_temp);

    for (read_idx = 0; read_idx < num_reads; read_idx++)
    {
        free(likelihood_matrix[read_idx]);
    }
    free(likelihood_matrix);
}