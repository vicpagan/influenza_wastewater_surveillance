#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

#include "calculate_proportions.h"
#include "external/covid_em.h"
#include "external/hashmap.h"

static int compare_proportions_desc(const void *a, const void *b)
{
    const ProportionData *pa = a;
    const ProportionData *pb = b;

    if (pa->proportion < pb->proportion)
        return 1;

    if (pa->proportion > pb->proportion)
        return -1;

    return 0;
}

static double squarem_objfn(const double *theta_1, const double **l_matrix, int num_reads, int num_msa_sequences)
{
    for (int i = 0; i < num_msa_sequences; ++i)
    {
        if (theta_1[i] <= 0.0 || theta_1[i] >= 1.0)
        {
            return 1.0e10;
        }
    }

    return negative_log_likelihood(theta_1, l_matrix, num_reads, num_msa_sequences);
}

double *run_squarem(const double *theta_0, const double **l_matrix, int num_reads, int num_msa_sequences)
{
    const int max_iter = 1500;
    const double convergence_tol = 1.0e-7;
    const double equality_tol = 1.0e-12;

    const double step_min0 = 1.0;
    const double step_max0 = 1.0;
    const double mstep = 4.0;

    const double objfn_inc = 1.0e-10;

    double step_min = step_min0;
    double step_max = step_max0;

    double *theta = malloc(num_msa_sequences * sizeof(double));
    double *theta_1 = malloc(num_msa_sequences * sizeof(double));
    double *theta_2 = malloc(num_msa_sequences * sizeof(double));

    double *r = malloc(num_msa_sequences * sizeof(double));
    double *v = malloc(num_msa_sequences * sizeof(double));

    double *theta_temp = malloc(num_msa_sequences * sizeof(double));
    double *theta_new = malloc(num_msa_sequences * sizeof(double));

    if (!theta || !theta_1 || !theta_2 || !r || !v || !theta_temp || !theta_new)
    {
        free(theta);
        free(theta_1);
        free(theta_2);
        free(r);
        free(v);
        free(theta_temp);
        free(theta_new);

        return NULL;
    }
    memcpy(theta, theta_0,  num_msa_sequences * sizeof(double));

    double obj_old = squarem_objfn(theta, l_matrix, num_reads, num_msa_sequences);

    for (int iter = 0; iter < max_iter; ++iter)
    {
        em_update(theta, l_matrix, num_reads, num_msa_sequences, theta_1);
        em_update(theta_1, l_matrix, num_reads, num_msa_sequences, theta_2);

        double q1_norm_sq = 0.0;
        double q2_norm_sq = 0.0;
        for (int i = 0; i < num_msa_sequences; ++i)
        {
            r[i] = theta_1[i] - theta[i];
            v[i] = theta_2[i] - theta_1[i];

            q1_norm_sq += r[i] * r[i];
            q2_norm_sq += v[i] * v[i];
        }

        double q1_norm = sqrt(q1_norm_sq);
        if (q1_norm < convergence_tol)
        {
            memcpy(theta, theta_1, num_msa_sequences * sizeof(double));
            break;
        }

        double q2_norm = sqrt(q2_norm_sq);
        if (q2_norm < convergence_tol)
        {
            memcpy(theta, theta_2, num_msa_sequences * sizeof(double));
            break;
        }

        double sv2 = 0.0;
        double srv = 0.0;
        for (int i = 0; i < num_msa_sequences; ++i)
        {
            v[i] -= r[i];

            sv2 += v[i] * v[i];
            srv += r[i] * v[i];
        }

        double alpha;
        if (sv2 > 0.0)
        {
            alpha = sqrt(q1_norm_sq / sv2);
        }
        else
        {
            alpha = 1.0;
        }
        alpha = fmax(step_min, fmin(step_max, alpha));

        memcpy(theta_new, theta, num_msa_sequences * sizeof(double));

        for (int i = 0; i < num_msa_sequences; ++i)
        {
            theta_new[i] += 2.0 * alpha * r[i];
            theta_new[i] += alpha * alpha * v[i];
        }

        if (fabs(alpha - 1.0) > 0.01)
        {
            em_update(theta_new, l_matrix, num_reads, num_msa_sequences, theta_temp);
            memcpy(theta_new, theta_temp, num_msa_sequences * sizeof(double));
        }

        int invalid_point = 0;
        for (int i = 0; i < num_msa_sequences; ++i)
        {
            if (isnan(theta_new[i]) || isinf(theta_new[i]))
            {
                invalid_point = 1;
                break;
            }
        }

        if (invalid_point)
        {
            memcpy(theta_new, theta_2, num_msa_sequences * sizeof(double));

            step_max = fmax(step_max0, step_max / mstep);
        }

        double obj_new = squarem_objfn(theta_new, l_matrix, num_reads, num_msa_sequences);;

        int extrapolation_succeeded = 1;
        if (obj_new > obj_old + objfn_inc)
        {
            memcpy(theta, theta_2, num_msa_sequences * sizeof(double));

            extrapolation_succeeded = 0;

            step_max = fmax(step_max0, step_max / mstep);
        }
        else
        {
            memcpy(theta, theta_new, num_msa_sequences * sizeof(double));

            obj_old = obj_new;
        }

        if (extrapolation_succeeded)
        {
            if (fabs(alpha - step_max) < equality_tol)
            {
                step_max = step_max * mstep;
            }
            
            if (fabs(alpha - step_min) < equality_tol && step_min < 0.0)
            {
                step_min = step_min * mstep;
            }
        }
    }

    free(theta_1);
    free(theta_2);
    free(r);
    free(v);
    free(theta_new);
    free(theta_temp);

    return theta;
}

void calculate_proportions(MismatchData *mismatch_data_str, char *output_dir, double error_rate, double filter, int compute_strain_llr, int num_top_strains_llr, int num_plot, int num_threads)
{
    srand((unsigned int)time(NULL));

    struct timespec start = {0, 0}, end = {0, 0};

    int i, read_idx, msa_seq_idx;
	int num_reads = mismatch_data_str->num_reads;
    if (num_reads <= 0)
    {
        printf("Error: No reads in mismatch matrix!\n");
        exit(1);
    }

    int num_msa_sequences = mismatch_data_str->num_msa_sequences;

    char **group_strains = (char **)malloc(num_msa_sequences * sizeof(char *)); 

    /////////////////////////////////////////////////////////////////////
    ////////////////// REMOVING UNIDENTIFIABLE STRAINS //////////////////
    /////////////////////////////////////////////////////////////////////

    printf("Grouping unidentifiable strains...\n");
	clock_gettime(CLOCK_MONOTONIC, &start);

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


    char unidentifiable_strains_path[4096];
    snprintf(unidentifiable_strains_path, sizeof(unidentifiable_strains_path), "%s/unidentifiable_strains.txt", output_dir);

    FILE *unidentifiable_strains_file = fopen(unidentifiable_strains_path, "w");
    if (unidentifiable_strains_file == NULL)
    {
        fprintf(stderr, "Error: Could not open unidentifiable strains file at '%s'\n", unidentifiable_strains_path);
        exit(1);
    }

    for (int group_idx = 0; group_idx < strain_names_to_group_idx; group_idx++)
    {
        fprintf(unidentifiable_strains_file, "Group %d:\n", group_idx);

        const char *current_char = group_strains[group_idx];
        while (*current_char != '\0')
        {
            if (*current_char == ':')
            {
                fputc('\n', unidentifiable_strains_file);
            }
            else
            {
                fputc(*current_char, unidentifiable_strains_file);
            }
            current_char++;
        }

        fputc('\n\n', unidentifiable_strains_file);

        free(group_strains[group_idx]);
    }

    fclose(unidentifiable_strains_file);
    free(group_strains);

    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("Took %.5fsec\n", ((double)end.tv_sec + 1.0e-9 * end.tv_nsec) - ((double)start.tv_sec + 1.0e-9 * start.tv_nsec));

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

    double *theta_0 = (double *)malloc(num_msa_sequences * sizeof(double));
    double *theta_0_rand = (double *)malloc(num_msa_sequences * sizeof(double));

    double theta_0_sum = 0.0;
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        double column_sum = 0.0;
        for (read_idx = 0; read_idx < num_reads; read_idx++)
        {
            column_sum += likelihood_matrix[read_idx][msa_seq_idx];
        }
        theta_0[msa_seq_idx] = column_sum;
        theta_0_sum += column_sum;
    }
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        theta_0[msa_seq_idx] /= theta_0_sum;
    }

    double theta_0_rand_sum = 0.0;
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        theta_0_rand[msa_seq_idx] = (double)rand() / ((double)RAND_MAX + 1.0);
        theta_0_rand_sum += theta_0_rand[msa_seq_idx];
    }
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        theta_0_rand[msa_seq_idx] /= theta_0_rand_sum;
    }

    ////////////////////////////////////////////////////////////////////
    /////////////////////////// SQUAREM RUNS ///////////////////////////
    ////////////////////////////////////////////////////////////////////

    printf("Starting SQUAREM step...\n");
	clock_gettime(CLOCK_MONOTONIC, &start);
    double *proportions = run_squarem(theta_0, likelihood_matrix, num_reads, num_msa_sequences);
    double *proportions_rand = run_squarem(theta_0_rand, likelihood_matrix, num_reads, num_msa_sequences);
    clock_gettime(CLOCK_MONOTONIC, &end);
    printf("Took %.5fsec\n", ((double)end.tv_sec + 1.0e-9 * end.tv_nsec) - ((double)start.tv_sec + 1.0e-9 * start.tv_nsec));

    if (!proportions || !proportions_rand)
    {
        fprintf(stderr, "Error: SQUAREM algorithm failed!\n");
        exit(1);
    }

    double difference = 0.0;
    for (int i = 0; i < num_msa_sequences; ++i)
    {
        difference += fabs(proportions[i] - proportions_rand[i]);
    }
    printf("\nDifference between two optimizations from different starting points %.2g\n", difference);\

    free(theta_0_rand);
    free(proportions_rand);

    ////////////////////////////////////////////////////////////////////
    ///////////////////////////// LLR STEP /////////////////////////////
    ////////////////////////////////////////////////////////////////////

    double log_likelihood_with = log_likelihood(proportions, likelihood_matrix, num_reads, num_msa_sequences);

    // TODO: Implement per-strain LLR step

    ////////////////////////////////////////////////////////////////////
    //////////////////////// OUTPUT PROPORTIONS ////////////////////////
    ////////////////////////////////////////////////////////////////////

    ProportionData *proportions_data = (ProportionData *)malloc(num_msa_sequences * sizeof(ProportionData));
    for (msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        proportions_data[msa_seq_idx].msa_strain_name = mismatch_data_str->msa_sequence_names[msa_seq_idx];
        proportions_data[msa_seq_idx].proportion = proportions[msa_seq_idx];
    }
    free(theta_0);
    free(proportions);

    qsort(proportions_data, (size_t)num_msa_sequences, sizeof(ProportionData), compare_proportions_desc);

    char output_csv_filepath[1050];
    snprintf(output_csv_filepath, sizeof(output_csv_filepath), "%s/proportions.csv", output_dir);

    FILE *output_csv_file = fopen(output_csv_filepath, "w");
    if (output_csv_file == NULL)
    {
        fprintf(stderr, "Error: Could not open output csv file at '%s'!\n", output_csv_filepath);
        exit(1);
    }

    fprintf(output_csv_file, "names.arg,p\n");
    for (int msa_seq_idx = 0; msa_seq_idx < num_msa_sequences; msa_seq_idx++)
    {
        fprintf(output_csv_file, "\"%s\",%.3f\n", proportions_data[msa_seq_idx].msa_strain_name, proportions_data[msa_seq_idx].proportion);
    }

    fclose(output_csv_file);

    ///////////////////////////////////////////////////////////////////
    /////////////////////////// FREE MEMORY ///////////////////////////
    ///////////////////////////////////////////////////////////////////

    free(proportions_data);

    for (read_idx = 0; read_idx < num_reads; read_idx++)
    {
        free(likelihood_matrix[read_idx]);
    }
    free(likelihood_matrix);
}