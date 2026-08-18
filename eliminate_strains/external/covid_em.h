#ifndef COVID_EM_H
#define COVID_EM_H


void calculate_q(const int *block_sizes, const int **mismatch_matrix, int num_reads, int num_strains, double error_rate, double **q_matrix);

void em_update(const double *proportions, const double **q_matrix, int num_reads, int num_strains, double *updated_proportions);

double log_likelihood(const double *proportions, const double **q_matrix, int num_reads, int num_strains);

double negative_log_likelihood(const double *proportions, const double **q_matrix, int num_reads, int num_strains);

#endif