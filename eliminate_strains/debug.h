#ifndef DEBUG_H
#define DEBUG_H

#include "global.h"

void write_mismatch_matrix(const MismatchData *mismatch_data_str, const char *filepath);

MismatchData read_mismatch_matrix(const char *filepath, int max_num_reads);

#endif // DEBUG_H