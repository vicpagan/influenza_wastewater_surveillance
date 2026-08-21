#ifndef COVID_EM_H
#define COVID_EM_H

void em_update(const double *proportions, const double **q_matrix, int num_reads, int num_strains, double *updated_proportions);

double log_likelihood(const double *proportions, const double **q_matrix, int num_reads, int num_strains);

double negative_log_likelihood(const double *proportions, const double **q_matrix, int num_reads, int num_strains);

#endif