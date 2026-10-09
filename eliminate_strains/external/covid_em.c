#include <math.h>
#include <float.h>
#include <omp.h>

#include "covid_em.h"

void em_update(const double *proportions, const double **q_matrix, int num_reads, int num_strains, double *updated_proportions, int num_threads)
{
	#pragma omp parallel for num_threads(num_threads)
    for (int j = 0; j < num_strains; j++)
    {
        updated_proportions[j] = 0.0;
    }

    #pragma omp parallel for num_threads(num_threads) reduction(+:updated_proportions[:num_strains])
    for (int i = 0; i < num_reads; i++)
    {
        double row_sum = 0.0;
        for (int j = 0; j < num_strains; j++)
        {
            row_sum += q_matrix[i][j] * proportions[j];
        }

        if (row_sum <= 0.0)
        {
            row_sum = DBL_MIN;
        }

        for (int j = 0; j < num_strains; j++)
        {
            updated_proportions[j] += q_matrix[i][j] * proportions[j] / row_sum;
        }
    }

	#pragma omp parallel for num_threads(num_threads)
    for (int j = 0; j < num_strains; j++)
    {
        updated_proportions[j] /= (double)num_reads;
    }
}

double log_likelihood(const double *proportions, const double **q_matrix, int num_reads, int num_strains, int num_threads)
{
    double log_likelihood = 0.0;

    #pragma omp parallel for num_threads(num_threads) reduction(+:log_likelihood)
    for (int i = 0; i < num_reads; i++)
    {
        double row_sum = 0.0;
        for (int j = 0; j < num_strains; j++)
        {
            row_sum += q_matrix[i][j] * proportions[j];
        }

        if (row_sum <= 0.0)
        {
            row_sum = DBL_MIN;
        }

        log_likelihood += log(row_sum);
    }
    return log_likelihood;
}

double negative_log_likelihood(const double *proportions, const double **q_matrix, int num_reads, int num_strains, int num_threads)
{
    return -log_likelihood(proportions, q_matrix, num_reads, num_strains, num_threads);
}