#include <math.h>
#include <float.h>

#include "covid_em.h"

void em_update(const double *proportions, const double **q_matrix, int num_reads, int num_strains, double *updated_proportions)
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

		// Added this line
		if (row_sum <= 0.0)
		{
			row_sum = DBL_MIN;
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

double log_likelihood(const double *proportions, const double **q_matrix, int num_reads, int num_strains)
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

double negative_log_likelihood(const double *proportions, const double **q_matrix, int num_reads, int num_strains)
{
	return -log_likelihood(proportions, q_matrix, num_reads, num_strains);
}
