#include <math.h>

#include "covid_em.h"

void calculate_q(const int *block_sizes, int **mismatch_matrix, int num_reads, int num_strains, double error_rate, double **q_matrix)
{
	int i,j;

	for (i = 0; i < num_reads; i++)
	{
		for (j = 0; j < num_strains; j++)
		{
			int mismatches = mismatch_matrix[i][j];
			int matches = block_sizes[i] - mismatches;

			q_matrix[i][j] = pow(error_rate, mismatches) * pow(1.0 - error_rate, matches);
		}
	}
}

void em_update(const double *proportions, double **q_matrix, int num_reads, int num_strains, double *updated_proportions)
{
	int i,j;

	for (j = 0; j < num_strains; j++)
	{
		updated_proportions[j] = 0.0;
	}

	for (i = 0; i < num_reads; i++)
	{
		double row_sum = 0.0;

		for (j = 0; j < num_strains; j++)
		{
			row_sum += q_matrix[i][j] * proportions[j];
		}

		for (j = 0; j < num_strains; j++)
		{
			updated_proportions[j] += q_matrix[i][j] * proportions[j] / row_sum;
		}
	}

	for (j = 0; j < num_strains; j++)
	{
		updated_proportions[j] /= (double)num_reads;
	}
}

double log_likelihood(const double *proportions, double **q_matrix, int num_reads, int num_strains)
{
	int i,j;
	double log_likelihood = 0.0;

	for (i = 0; i < num_reads; i++)
	{
		double row_sum = 0.0;

		for (j = 0; j < num_strains; j++)
		{
			row_sum += q_matrix[i][j] * proportions[j];
		}

		log_likelihood += log(row_sum);
	}

	return log_likelihood;
}

double negative_log_likelihood(const double *proportions, double **q_matrix, int num_reads, int num_strains)
{
	return -log_likelihood(proportions, q_matrix, num_reads, num_strains);
}