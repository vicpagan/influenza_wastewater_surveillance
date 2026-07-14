#include "calculate_allele_freq.h"

/**
 * @brief Compute per-site allele frequencies from paired-end reads, then iteratively
 * eliminate reference strains incompatible with those frequencies.
 *
 * Three stages, in order:
 *  1. Walk the SAM file read-pair by read-pair, tallying A/C/G/T counts (and deletions)
 *     at each MSA site the reads cover, translating SAM positions to MSA columns via
 *     reference_index[]. Overlapping mate pairs are only counted once (see `visited`).
 *  2. Convert counts to frequencies, drop sites below `coverage`, and mark any base
 *     with frequency < freq_threshold as "bad" at that site.
 *  3. Iteratively remove (blank the name of) any strain carrying too many "bad" bases
 *     at the variant sites, growing the incompatibility tolerance each pass until the
 *     number of strains remaining falls within [min_strains_remaining, max_strains_remaining).
 *
 * @param sam Open SAM file (paired-end alignments).
 * @param allele Output/scratch: per-site A/C/G/T counts, then converted to frequencies.
 * @param length_of_MSA Number of MSA columns.
 * @param MSA Numerically-coded... (actually char-coded) alignment; see names_of_strains.
 * @param number_of_strains Total strains before elimination.
 * @param names_of_strains Strain names; eliminated strains get their name blanked.
 * @param freq_threshold Below this frequency, a base is "bad" at a site.
 * @param maxname Length of each name buffer.
 * @param tstart Scratch timing variable (reused internally).
 * @param tend Scratch timing variable (reused internally).
 * @param number_of_variant_sites Length of variant_sites.
 * @param variant_sites Sites to check for strain elimination; entries get set to -1
 * if not adequately covered.
 * @param coverage Minimum read depth to trust a site's allele frequency.
 * @param reference_index Maps SAM alignment position -> MSA column (see align_reference.c).
 * @param min_strains_remaining Lower bound on strains remaining; loop stops once reached.
 * @param max_strains_remaining Upper bound; exceeding this aborts the whole program.
 * @param print_counts If non-empty, path to dump per-site allele counts.
 * @param max_sam_length Output: [0] = longest raw SAM line seen, [1] = number of SAM records.
 * @param print_deletions If non-empty, path to dump sites with deletion frequency above threshold.
 * @param deletion_threshold Frequency threshold for reporting a deletion site.
 * @param sam_results_out Output: newly allocated array of every raw SAM
 * line read (same format readInSamFile() would produce), so callers can reuse
 * it later without a second disk read. Pass NULL to skip caching.
 * @param num_sam_lines_out Output: number of lines in sam_results_out.
 * 
 * @return Number of strains remaining after elimination.
 */
int calculateAlleleFreq_paired(FILE *sam, double **allele, int length_of_MSA, char **MSA, int number_of_strains, char **names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int number_of_variant_sites, int *variant_sites, int coverage, int *reference_index, int min_strains_remaining, int max_strains_remaining, char print_counts[], int *max_sam_length, char print_deletions[], double deletion_threshold, char ***sam_results_out, int *num_sam_lines_out)
{
	int i, j;
	char buffer[FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar[i] = 0;
		cigar_chars[i] = '\0';
	}
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	int first_in_pair = 0;
	int second_in_pair = 0;
	int first_seq_length = 0;
	int first_seq_cigar[MAX_CIGAR];
	char first_seq_cigar_chars[MAX_CIGAR];
	int first_seq_cigar_char_count = 0;
	char *first_seq = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	memset(first_seq, '\0', MAX_READ_LENGTH);
	int first_end_pos = 0;
	int first_seq_start_pos = 0;
	int visited[MAX_READ_LENGTH];
	// memset(visited,-1,MAX_READ_LENGTH);
	for (i = 0; i < MAX_READ_LENGTH; i++)
	{
		visited[i] = -1;
	}
	double *deletions = (double *)malloc((length_of_MSA + 1) * sizeof(double));
	for (i = 0; i < length_of_MSA + 1; i++)
	{
		deletions[i] = 0;
	}
	// Growable cache of every raw SAM line, filled in the same pass as the
	// tallying below -- this is what lets callers avoid a second disk read
	// (e.g. via readInSamFile()) later to build the mismatch matrix.
	char **cached_lines = NULL;
	int cached_capacity = 0;
	int cached_count_local = 0;

	// --- Stage 1: tally per-site A/C/G/T counts (and deletions) from every read pair ---
	while (fgets(buffer, FASTA_MAXLINE, sam) != NULL)
	{
		if (buffer[0] != '@')
		{
			max_sam_length[1]++;
			if (sam_results_out != NULL)
			{
				if (cached_count_local == cached_capacity)
				{
					cached_capacity = cached_capacity == 0 ? 1024 : cached_capacity * 2;
					cached_lines = (char **)realloc(cached_lines, cached_capacity * sizeof(char *));
				}
				cached_lines[cached_count_local] = strdup(buffer);
				cached_count_local++;
			}
			char *buffer_copy = strdup(buffer);
			int length_of_sam = strlen(buffer);
			if (length_of_sam > max_sam_length[0])
			{
				max_sam_length[0] = length_of_sam;
			}
			s = strtok(buffer, "\t");
			// char* name = strdup(s);
			s = strtok(NULL, "\t");
			int decimal = 0;
			sscanf(s, "%d", &decimal);
			decimal = dec2bin(decimal);
			if (decimal == 1)
			{
				first_in_pair = 1;
			}
			else if (decimal == 0)
			{
				second_in_pair = 1;
			}
			// free(name);
			for (i = 0; i < 2; i++)
			{
				s = strtok(NULL, "\t");
			}
			int position = 0;
			sscanf(s, "%d", &position);
			position--;
			if (first_in_pair == 1)
			{
				first_seq_start_pos = position;
			}
			s = strtok(NULL, "\t");
			s = strtok(NULL, "\t");
			char *cigar_string;
			cigar_string = strdup(s);
			// printf("%s\n",cigar_string);
			if (strcmp(cigar_string, "*") != 0)
			{
				char *copy = strdup(cigar_string);
				char *res = strtok(cigar_string, "MID");
				int index = 0;
				while (res)
				{
					int from = res - cigar_string + strlen(res);
					// printf("%s\n",res);
					int cigar_count = 0;
					sscanf(res, "%d", &cigar_count);
					// printf("%d\n",cigar_count);
					res = strtok(NULL, "MID");
					int to = res != NULL ? res - cigar_string : strlen(copy);
					// printf("%.*s\n", to-from, copy+from);
					char cigar_char = '\0';
					sscanf(copy + from, "%c", &cigar_char);
					// printf("cigar char: %c\n",cigar_char);
					cigar[index] = cigar_count;
					if (first_in_pair == 1)
					{
						first_seq_cigar[index] = cigar_count;
					}
					cigar_chars[index] = cigar_char;
					if (first_in_pair == 1)
					{
						first_seq_cigar_chars[index] = cigar_char;
					}
					index++;
				}
				// printf("%s\n",cigar_string);
				free(copy);
				// printf("S: %s\n",buffer_copy);
				s = strtok(buffer_copy, "\t");
				// printf("S: %s\n",s);
				for (i = 0; i < 9; i++)
				{
					s = strtok(NULL, "\t");
				}
				char *sequence = s;
				if (first_in_pair == 1)
				{
					strcpy(first_seq, sequence);
					first_seq_length = strlen(first_seq);
				}
				int cigar_char_count = index;
				if (first_in_pair == 1)
				{
					first_seq_cigar_char_count = index;
				}
				double alignment_score = 0;
				int matches = 0;
				int ns = 0;
				int is = 0;
				for (i = 0; i < cigar_char_count; i++)
				{
					if (cigar_chars[i] == 'M')
					{
						matches = matches + cigar[i];
					}
					if (cigar_chars[i] == 'I')
					{
						ns = ns + cigar[i];
					}
					if (cigar_chars[i] == 'N')
					{
						is = is + cigar[i];
					}
				}
				if (cigar_char_count > 0)
				{
					alignment_score = matches / (matches + ns + is);
				}
				index = 0;
				int start = 0;
				int start_ref = 0;
				int k = 0;
				int l = 0;
				int visited_place = 0;
				if (second_in_pair == 1 /*&& alignment_score > 0.8*/)
				{
					int start_ref1 = 0;
					int start1 = 0;
					for (i = 0; i < first_seq_cigar_char_count; i++)
					{
						for (j = 0; j < first_seq_cigar[i]; j++)
						{
							int position_in_MSA = reference_index[start_ref1 + first_seq_start_pos + j];
							if (position_in_MSA != -1 && position_in_MSA >= position)
							{
								int start2 = 0;
								int start_ref2 = 0;
								for (k = 0; k < cigar_char_count; k++)
								{
									for (l = 0; l < cigar[k]; l++)
									{
										if (position_in_MSA != -1 && position_in_MSA == reference_index[start_ref2 + position + l] && cigar_chars[k] == 'M' && first_seq_cigar_chars[i] == 'M')
										{
											if (first_seq[start1 + j] != sequence[start2 + l] && first_seq_start_pos + j <= variant_sites[number_of_variant_sites - 1])
											{
												if (first_seq[start1 + j] == 'A' || first_seq[start1 + j] == 'a')
												{
													allele[position_in_MSA][0]--;
												}
												else if (first_seq[start1 + j] == 'G' || first_seq[start1 + j] == 'g')
												{
													allele[position_in_MSA][1]--;
												}
												else if (first_seq[start1 + j] == 'C' || first_seq[start1 + j] == 'c')
												{
													allele[position_in_MSA][2]--;
												}
												else if (first_seq[start1 + j] == 'T' || first_seq[start1 + j] == 't')
												{
													allele[position_in_MSA][3]--;
												}
											}
											visited[visited_place] = position_in_MSA;
											visited_place++;
										}
									}
									if (cigar_chars[k] == 'M')
									{
										start2 = start2 + cigar[k];
										start_ref2 = start_ref2 + cigar[k];
									}
									if (cigar_chars[k] == 'I')
									{
										start2 = start2 + cigar[k];
									}
									if (cigar_chars[k] == 'D')
									{
										start_ref2 = start_ref2 + cigar[k];
									}
								}
							}
						}
						if (first_seq_cigar_chars[i] == 'M')
						{
							start1 = start1 + first_seq_cigar[i];
							start_ref1 = start_ref1 + first_seq_cigar[i];
						}
						if (first_seq_cigar_chars[i] == 'I')
						{
							start1 = start1 + first_seq_cigar[i];
						}
						if (first_seq_cigar_chars[i] == 'D')
						{
							start_ref1 = first_seq_cigar[i] + start_ref1;
						}
					}
				}
				// if ( alignment_score > 0.8){
				for (i = 0; i < cigar_char_count; i++)
				{
					// printf("cigar[%d]=%d\n",i,cigar[i]);
					for (j = 0; j < cigar[i]; j++)
					{
						// printf("cigar_chars[%d]: %c\n",i,cigar_chars[i]);
						int skip = 0;
						int position_in_MSA = reference_index[j + start_ref + position];
						for (k = 0; k < visited_place; k++)
						{
							if (position_in_MSA != -1 && visited[k] == position_in_MSA)
							{
								skip = 1;
							}
						}
						if (cigar_chars[i] == 'M' /*|| cigar_chars[i] == 'I'*/)
						{
							if (sequence[j + start] == 'A' || sequence[j + start] == 'a')
							{
								// if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								//	allele[reference[j+start_ref+position]][0]++;
								// }
								if (position_in_MSA != -1 && position_in_MSA < length_of_MSA && skip == 0)
								{
									allele[position_in_MSA][0]++;
								}
							}
							else if (sequence[j + start] == 'G' || sequence[j + start] == 'g')
							{
								// if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								//	allele[reference[j+start_ref+position]][1]++;
								// }
								if (position_in_MSA != -1 && position_in_MSA < length_of_MSA && skip == 0)
								{
									allele[position_in_MSA][1]++;
								}
							}
							else if (sequence[j + start] == 'C' || sequence[j + start] == 'c')
							{
								// if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								//	allele[reference[j+start_ref+position]][2]++;
								// }
								if (position_in_MSA != -1 && position_in_MSA < length_of_MSA && skip == 0)
								{
									allele[position_in_MSA][2]++;
								}
							}
							else if (sequence[j + start] == 'T' || sequence[j + start] == 't')
							{
								// if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
								//	allele[reference[j+start_ref+position]][3]++;
								// }
								if (position_in_MSA != -1 && position_in_MSA < length_of_MSA && skip == 0)
								{
									allele[position_in_MSA][3]++;
								}
							}
							first_end_pos = position_in_MSA;
						}
						if (cigar_chars[i] == 'I')
						{
							first_end_pos = position_in_MSA;
						}
						if (cigar_chars[i] == 'D')
						{
							printf("POS: %d\n", position_in_MSA);
							if (position_in_MSA != -1)
							{
								deletions[position_in_MSA]++;
							}
						}
					}
					if (cigar_chars[i] == 'M')
					{
						start = cigar[i] + start;
						start_ref = cigar[i] + start_ref;
					}
					if (cigar_chars[i] == 'I')
					{
						start = cigar[i] + start;
					}
					if (cigar_chars[i] == 'D')
					{
						start_ref = cigar[i] + start_ref;
					}
				}
			}
			//}
			free(buffer_copy);
			free(cigar_string);
			if (second_in_pair == 1)
			{
				memset(first_seq, '\0', MAX_READ_LENGTH);
				first_seq_cigar_char_count = 0;
			}
			first_in_pair = 0;
			second_in_pair = 0;
		}
	}
	free(first_seq);
	if (print_deletions[0] != '\0')
	{
		FILE *deletion_sites_file;
		if ((deletion_sites_file = fopen(print_deletions, "w")) == (FILE *)NULL)
			fprintf(stderr, "Deletion sites file could not be opened.\n");
		fprintf(deletion_sites_file, "Site\tFrequency\n");
		for (i = 0; i < length_of_MSA; i++)
		{
			// double deletion_freq;
			// deletion_freq = deletions[i]/max_sam_length[1];
			if (deletions[i] / max_sam_length[1] > deletion_threshold)
			{
				fprintf(deletion_sites_file, "%d\t%lf\n", i, deletions[i] / max_sam_length[1]);
			}
		}
		fclose(deletion_sites_file);
	}
	// free(deletions);
	// --- Stage 2: convert counts to frequencies, drop low-coverage sites, mark "bad" bases ---
	int covered = 0;
	for (i = 0; i < length_of_MSA; i++)
	{
		double total = 0;
		for (j = 0; j < 4; j++)
		{
			total = total + allele[i][j];
		}
		if (total > 0)
		{
			covered++;
		}
	}
	int *covered_sites = (int *)malloc(covered * sizeof(int));
	for (i = 0; i < covered; i++)
	{
		covered_sites[i] = -1;
	}
	int k = 0;
	if (print_counts[0] != '\0')
	{
		FILE *allele_counts_file;
		if ((allele_counts_file = fopen(print_counts, "w")) == (FILE *)NULL)
			fprintf(stderr, "Allele Counts file could not be opened.\n");
		fprintf(allele_counts_file, "position\tA\tG\tC\tT\n");
		for (i = 0; i < length_of_MSA; i++)
		{
			fprintf(allele_counts_file, "%d\t%lf\t%lf\t%lf\t%lf\n", i, allele[i][0], allele[i][1], allele[i][2], allele[i][3]);
		}
		fclose(allele_counts_file);
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		double total = 0;
		for (j = 0; j < 4; j++)
		{
			total = total + allele[i][j];
		}
		if (total >= coverage)
		{
			covered_sites[k] = i;
			k++;
		}
		for (j = 0; j < 4; j++)
		{
			allele[i][j] = allele[i][j] / total;
		}
	}
	printf("Number of sites not covered: %d\n", length_of_MSA - k);
	for (i = 0; i < number_of_variant_sites; i++)
	{
		int found = 0;
		for (j = 0; j < covered; j++)
		{
			if (variant_sites[i] == covered_sites[j])
			{
				found = 1;
			}
		}
		if (found == 0)
		{
			variant_sites[i] = -1;
		}
	}
	free(covered_sites);
	int temp_num_var_sites = 0;
	for (i = 0; i < number_of_variant_sites; i++)
	{
		if (variant_sites[i] >= 0 && variant_sites[i] < length_of_MSA - 1)
		{
			temp_num_var_sites++;
		}
	}
	int *variant_sites_updated = (int *)malloc(temp_num_var_sites * sizeof(int));
	k = 0;
	for (i = 0; i < number_of_variant_sites; i++)
	{
		if (variant_sites[i] != -1 && variant_sites[i] < length_of_MSA - 1)
		{
			variant_sites_updated[k] = variant_sites[i];
			k++;
		}
	}
	free(variant_sites);
	int **bad_bases = (int **)malloc(length_of_MSA * sizeof(int *));
	for (i = 0; i < length_of_MSA; i++)
	{
		bad_bases[i] = (int *)malloc(4 * sizeof(int));
		for (j = 0; j < 4; j++)
		{
			bad_bases[i][j] = 0;
		}
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		for (j = 0; j < 4; j++)
		{
			if (allele[i][j] < freq_threshold)
			{
				bad_bases[i][j] = 1;
			}
		}
	}
	int *bad_bases_count = (int *)malloc(length_of_MSA * sizeof(int));
	for (i = 0; i < length_of_MSA; i++)
	{
		bad_bases_count[i] = 0;
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		for (j = 0; j < 4; j++)
		{
			bad_bases_count[i] += bad_bases[i][j];
		}
	}
	char **bad_base_char = (char **)malloc(length_of_MSA * sizeof(char *));
	for (i = 0; i < length_of_MSA; i++)
	{
		bad_base_char[i] = (char *)malloc(4 * sizeof(char));
		for (j = 0; j < 4; j++)
		{
			bad_base_char[i][j] = '\0';
		}
	}
 
	for (i = 0; i < length_of_MSA; i++)
	{
		j = 0;
		// for(j=0; j<4-bad_bases_count[i]; j++){
		for (k = 0; k < 4; k++)
		{
			if (bad_bases[i][k] == 0)
			{
				if (k == 0)
				{
					bad_base_char[i][j] = 'A';
					j++;
				}
				else if (k == 1)
				{
					bad_base_char[i][j] = 'G';
					j++;
				}
				else if (k == 2)
				{
					bad_base_char[i][j] = 'C';
					j++;
				}
				else if (k == 3)
				{
					bad_base_char[i][j] = 'T';
					j++;
				}
				bad_bases[i][k] = 1;
			}
		}
		//}
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		free(bad_bases[i]);
	}
	free(bad_bases);
	// int length_of_reference = 0;
	// for (i=0; i<length_of_MSA; i++){
	//	if ( reference[i]==-1 ){
	//		break;
	//	}else{
	//		length_of_reference++;
	//	}
	// }
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	// --- Stage 3: iteratively eliminate strains incompatible with the "bad" bases above ---
	int number_remaining = number_of_strains;
	int base;
	int count;
	int var_count = 0;
	int number_removed = 0;
	int run_loop = 1;
	int number_of_iterations = 1;
	int *incompat_counter = (int *)malloc(number_of_strains * sizeof(int));
	for (i = 0; i < number_of_strains; i++)
	{
		incompat_counter[i] = 0;
	}
	while (run_loop == 1)
	{
		printf("iteration %d\n", number_of_iterations);
		for (i = 0; i < number_of_strains; i++)
		{
			incompat_counter[i] = 0;
		}
		number_remaining = number_of_strains;
		for (i = 0; i < temp_num_var_sites; i++)
		{
			number_removed = 0;
			for (j = 0; j < number_of_strains; j++)
			{
				/*if ( MSA[identical[j]][reference[i]] != '-' ){
					int base;
					if ( MSA[identical[j]][reference[i]] == 'A' ){
						base=0;
					}else if ( MSA[identical[j]][reference[i]] == 'G' ){
						base=1;
					}else if ( MSA[identical[j]][reference[i]] == 'C' ){
						base=2;
					}else if ( MSA[identical[j]][reference[i]] == 'T' ){
						base=3;
					}
					if ( allele[reference[i]][base] < freq_threshold ){
						//memset(names_of_strains[identical[j]],'\0',maxname);
						names_of_strains[identical[j]][0] = '\0';
					}
				}*/
				// if ( names_of_strains[j][0] != '\0'){
				if (incompat_counter[j] < number_of_iterations)
				{
					count = 0;
					for (k = 0; k < 4 - bad_bases_count[variant_sites_updated[i]]; k++)
					{
						if (MSA[j][variant_sites_updated[i]] != bad_base_char[variant_sites_updated[i]][k])
						{
							count++;
						}
					}
					if (count == (4 - bad_bases_count[variant_sites_updated[i]]))
					{
						/*if ( number_remaining ==1 ){
							printf("site %d\n",variant_sites_updated[i]);
							printf("last one: %s\n",names_of_strains[j]);
							printf("MSA[%d][%d]: %c\n",j,variant_sites_updated[i],MSA[j][variant_sites_updated[i]]);
							printf("bad_base_char[%d][0]: %c\n",variant_sites_updated[i],bad_base_char[variant_sites_updated[i]][0]);
						}*/
						// printf("removing %d: %s at site %d\n",j,names_of_strains[j],variant_sites_updated[i]);
						// if ( strcmp(names_of_strains[j],"EPI_ISL_4510987")==0){
						//	printf("removing %d: %s at site %d\n",j,names_of_strains[j],variant_sites_updated[i]);
						// }
						// names_of_strains[j][0] = '\0';
						// if ( incompat_counter[j]==number_of_iterations){
						// number_remaining--;
						// number_removed++;
						incompat_counter[j]++;
						//}
					}
				}
			}
		}
		for (i = 0; i < number_of_strains; i++)
		{
			if (incompat_counter[i] == number_of_iterations)
			{
				number_remaining--;
				number_removed++;
			}
		}
		if (number_remaining >= min_strains_remaining && number_remaining < max_strains_remaining)
		{
			printf("exiting loop. %d remaining\n", number_remaining);
			run_loop = 0;
		}
		else if (number_remaining >= max_strains_remaining)
		{
			printf("%d strains remaining. exiting...\n", number_remaining);
			exit(1);
		}
		else
		{
			printf("there are %d remaining... \n", number_remaining);
			number_of_iterations++;
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	for (i = 0; i < number_of_strains; i++)
	{
		if (incompat_counter[i] == number_of_iterations)
		{
			names_of_strains[i][0] = '\0';
		}
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		free(bad_base_char[i]);
	}
	free(bad_base_char);
	free(bad_bases_count);
	free(variant_sites_updated);
	number_remaining = 0;
	for (i = 0; i < number_of_strains; i++)
	{
		if (names_of_strains[i][0] != '\0')
		{
			// printf("Remaining strain: %s\n",names_of_strains[i]);
			number_remaining++;
		}
	}
	printf("Number remaining: %d\n", number_remaining);
	free(incompat_counter);
	if (sam_results_out != NULL)
	{
		*sam_results_out = cached_lines;
		*num_sam_lines_out = cached_count_local;
	}
	return number_remaining;
}

/**
 * @brief Single-end counterpart of calculateAlleleFreq_paired(): same three stages
 * (tally allele counts -> threshold "bad" bases -> iteratively eliminate strains),
 * without the mate-pair overlap handling paired-end reads need.
 * @see calculateAlleleFreq_paired for full parameter documentation (identical here),
 * including sam_results_out/num_sam_lines_out.
 */
int calculateAlleleFreq(FILE *sam, double **allele, int length_of_MSA, char **MSA, int number_of_strains, char **names_of_strains, double freq_threshold, int maxname, struct timespec tstart, struct timespec tend, int number_of_variant_sites, int *variant_sites, int coverage, int *reference_index, int min_strains_remaining, int max_strains_remaining, char print_counts[], int *max_sam_length, char print_deletions[], double deletion_threshold, char ***sam_results_out, int *num_sam_lines_out)
{
	int i, j;
	char buffer[FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar[i] = 0;
		cigar_chars[i] = '\0';
	}
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	double *deletions = (double *)malloc(length_of_MSA * sizeof(double));
	for (i = 0; i < length_of_MSA; i++)
	{
		deletions[i] = 0;
	}
	// Growable cache of every raw SAM line, filled in the same pass as the
	// tallying below -- this is what lets callers avoid a second disk read
	// (e.g. via readInSamFile()) later to build the mismatch matrix.
	char **cached_lines = NULL;
	int cached_capacity = 0;
	int cached_count_local = 0;
	// --- Stage 1: tally per-site A/C/G/T counts (and deletions) from every read ---
	while (fgets(buffer, FASTA_MAXLINE, sam) != NULL)
	{
		if (buffer[0] != '@')
		{
			max_sam_length[1]++;
			if (sam_results_out != NULL)
			{
				if (cached_count_local == cached_capacity)
				{
					cached_capacity = cached_capacity == 0 ? 1024 : cached_capacity * 2;
					cached_lines = (char **)realloc(cached_lines, cached_capacity * sizeof(char *));
				}
				cached_lines[cached_count_local] = strdup(buffer);
				cached_count_local++;
			}
			char *buffer_copy = strdup(buffer);
			int length_of_sam = strlen(buffer);
			if (length_of_sam > max_sam_length[0])
			{
				max_sam_length[0] = length_of_sam;
			}
			s = strtok(buffer, "\t");
			for (i = 0; i < 3; i++)
			{
				s = strtok(NULL, "\t");
			}
			int position = 0;
			sscanf(s, "%d", &position);
			position--;
			s = strtok(NULL, "\t");
			s = strtok(NULL, "\t");
			char *cigar_string;
			cigar_string = strdup(s);
			char *copy = strdup(cigar_string);
			char *res = strtok(cigar_string, "MID");
			int index = 0;
			while (res)
			{
				int from = res - cigar_string + strlen(res);
				int cigar_count = 0;
				sscanf(res, "%d", &cigar_count);
				res = strtok(NULL, "MID");
				int to = res != NULL ? res - cigar_string : strlen(copy);
				char cigar_char = '\0';
				sscanf(copy + from, "%c", &cigar_char);
				cigar[index] = cigar_count;
				cigar_chars[index] = cigar_char;
				index++;
			}
			free(copy);
			free(cigar_string);
			s = strtok(buffer_copy, "\t");
			for (i = 0; i < 9; i++)
			{
				s = strtok(NULL, "\t");
			}
			int cigar_char_count = index;
			int start = 0;
			int start_ref = 0;
			char *sequence = s;
			for (i = 0; i < cigar_char_count; i++)
			{
				for (j = 0; j < cigar[i]; j++)
				{
					if (cigar_chars[i] == 'M')
					{
						int position_in_MSA = reference_index[j + start_ref + position];
						if (sequence[j + start] == 'A' || sequence[j + start] == 'a')
						{
							// if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][0]++;
							// }
							if (position_in_MSA != -1 && position_in_MSA < length_of_MSA)
							{
								allele[position_in_MSA][0]++;
							}
						}
						else if (sequence[j + start] == 'G' || sequence[j + start] == 'g')
						{
							// if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][1]++;
							// }
							if (position_in_MSA != -1 && position_in_MSA < length_of_MSA)
							{
								allele[position_in_MSA][1]++;
							}
						}
						else if (sequence[j + start] == 'C' || sequence[j + start] == 'c')
						{
							// if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][2]++;
							// }
							if (position_in_MSA != -1 && position_in_MSA < length_of_MSA)
							{
								allele[position_in_MSA][2]++;
							}
						}
						else if (sequence[j + start] == 'T' || sequence[j + start] == 't')
						{
							// if ( reference[j+start_ref+position] < length_of_MSA && reference[j+start_ref+position] != -1){
							//	allele[reference[j+start_ref+position]][3]++;
							// }
							if (position_in_MSA != -1 && position_in_MSA < length_of_MSA)
							{
								allele[position_in_MSA][3]++;
							}
						}
					}
				}
				if (cigar_chars[i] == 'M')
				{
					start = cigar[i] + start;
					start_ref = cigar[i] + start_ref;
				}
				if (cigar_chars[i] == 'I')
				{
					start = cigar[i] + start;
				}
				if (cigar_chars[i] == 'D')
				{
					start_ref = cigar[i] + start_ref;
					for (j = 0; j < cigar[i]; j++)
					{
						deletions[reference_index[j + start_ref + position]]++;
					}
				}
			}
			free(buffer_copy);
		}
	}
	if (print_counts[0] != '\0')
	{
		FILE *allele_counts_file;
		if ((allele_counts_file = fopen(print_counts, "w")) == (FILE *)NULL)
			fprintf(stderr, "Allele Counts file could not be opened.\n");
		fprintf(allele_counts_file, "position\tA\tG\tC\tT\n");
		for (i = 0; i < length_of_MSA; i++)
		{
			fprintf(allele_counts_file, "%d\t%lf\t%lf\t%lf\t%lf\n", i, allele[i][0], allele[i][1], allele[i][2], allele[i][3]);
		}
		fclose(allele_counts_file);
	}
	if (print_deletions[0] != '\0')
	{
		FILE *deletion_sites_file;
		if ((deletion_sites_file = fopen(print_deletions, "w")) == (FILE *)NULL)
			fprintf(stderr, "Deletion sites file could not be opened.\n");
		fprintf(deletion_sites_file, "Site\tFrequency\n");
		for (i = 0; i < length_of_MSA; i++)
		{
			if (deletions[i] / max_sam_length[1] > deletion_threshold)
			{
				fprintf(deletion_sites_file, "%d\t%lf\n", i, deletions[i] / max_sam_length[1]);
			}
		}
		fclose(deletion_sites_file);
	}
	// --- Stage 2: convert counts to frequencies, drop low-coverage sites, mark "bad" bases ---
	int covered = 0;
	for (i = 0; i < length_of_MSA; i++)
	{
		double total = 0;
		for (j = 0; j < 4; j++)
		{
			total = total + allele[i][j];
		}
		if (total > 0)
		{
			covered++;
		}
	}
	int *covered_sites = (int *)malloc(covered * sizeof(int));
	for (i = 0; i < covered; i++)
	{
		covered_sites[i] = -1;
	}
	int k = 0;
	for (i = 0; i < length_of_MSA; i++)
	{
		double total = 0;
		for (j = 0; j < 4; j++)
		{
			total = total + allele[i][j];
		}
		if (total >= coverage)
		{
			covered_sites[k] = i;
			k++;
		}
		for (j = 0; j < 4; j++)
		{
			allele[i][j] = allele[i][j] / total;
		}
	}
	printf("Number of sites not covered: %d\n", length_of_MSA - k);
	for (i = 0; i < number_of_variant_sites; i++)
	{
		int found = 0;
		for (j = 0; j < covered; j++)
		{
			if (variant_sites[i] == covered_sites[j])
			{
				found = 1;
			}
		}
		if (found == 0)
		{
			variant_sites[i] = -1;
		}
	}
	free(covered_sites);
	int temp_num_var_sites = 0;
	for (i = 0; i < number_of_variant_sites; i++)
	{
		if (variant_sites[i] >= 0 && variant_sites[i] < length_of_MSA - 1)
		{
			temp_num_var_sites++;
		}
	}
	int *variant_sites_updated = (int *)malloc(temp_num_var_sites * sizeof(int));
	k = 0;
	for (i = 0; i < number_of_variant_sites; i++)
	{
		if (variant_sites[i] != -1 && variant_sites[i] < length_of_MSA - 1)
		{
			variant_sites_updated[k] = variant_sites[i];
			k++;
		}
	}
	free(variant_sites);
	
	int **bad_bases = (int **)malloc(length_of_MSA * sizeof(int *));
	for (i = 0; i < length_of_MSA; i++)
	{
		bad_bases[i] = (int *)malloc(4 * sizeof(int));
		for (j = 0; j < 4; j++)
		{
			bad_bases[i][j] = 0;
		}
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		for (j = 0; j < 4; j++)
		{
			if (allele[i][j] < freq_threshold)
			{
				bad_bases[i][j] = 1;
			}
		}
	}
	int *bad_bases_count = (int *)malloc(length_of_MSA * sizeof(int));
	for (i = 0; i < length_of_MSA; i++)
	{
		bad_bases_count[i] = 0;
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		for (j = 0; j < 4; j++)
		{
			bad_bases_count[i] += bad_bases[i][j];
		}
	}
	char **bad_base_char = (char **)malloc(length_of_MSA * sizeof(char *));
	for (i = 0; i < length_of_MSA; i++)
	{
		bad_base_char[i] = (char *)malloc(4 * sizeof(char));
		for (j = 0; j < 4; j++)
		{
			bad_base_char[i][j] = '\0';
		}
	}
 
	for (i = 0; i < length_of_MSA; i++)
	{
		j = 0;
		// for(j=0; j<4-bad_bases_count[i]; j++){
		for (k = 0; k < 4; k++)
		{
			if (bad_bases[i][k] == 0)
			{
				if (k == 0)
				{
					bad_base_char[i][j] = 'A';
					j++;
				}
				else if (k == 1)
				{
					bad_base_char[i][j] = 'G';
					j++;
				}
				else if (k == 2)
				{
					bad_base_char[i][j] = 'C';
					j++;
				}
				else if (k == 3)
				{
					bad_base_char[i][j] = 'T';
					j++;
				}
				bad_bases[i][k] = 1;
			}
		}
		//}
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		free(bad_bases[i]);
	}
	free(bad_bases);
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	printf("Eliminating strains...\n");
	clock_gettime(CLOCK_MONOTONIC, &tstart);
	// --- Stage 3: iteratively eliminate strains incompatible with the "bad" bases above ---
	int number_remaining = 0;
	int base;
	int count;
	int var_count = 0;
	/*int* nums_of_strains = (int*)malloc(number_of_different_strains*sizeof(int));
	for(i=0; i<number_of_different_strains; i++){
		nums_of_strains[i]=-1;
	}
	int place=0;*/
	int number_removed = 0;
	int run_loop = 1;
	int number_of_iterations = 1;
	int *incompat_counter = (int *)malloc(number_of_strains * sizeof(int));
	while (run_loop == 1)
	{
		printf("iteration %d\n", number_of_iterations);
		for (i = 0; i < number_of_strains; i++)
		{
			incompat_counter[i] = 0;
		}
		number_remaining = number_of_strains;
		for (i = 0; i < temp_num_var_sites; i++)
		{
			// if ( variant_sites[i] == covered_sites[next] ){
			for (j = 0; j < number_of_strains; j++)
			{
				// if ( names_of_strains[j][0] != '\0'){
				if (incompat_counter[j] < number_of_iterations)
				{
					/*if ( MSA[j][variant_sites_updated[i]] == 'A' ){
						base=0;
					}else if ( MSA[j][variant_sites_updated[i]] == 'G' ){
						base=1;
					}else if ( MSA[j][variant_sites_updated[i]] == 'C' ){
						base=2;
					}else if ( MSA[j][variant_sites_updated[i]] == 'T' ){
						base=3;
					}
					if ( allele[variant_sites_updated[i]][base] < freq_threshold ){
						//memset(names_of_strains[identical[j]],'\0',maxname);
						names_of_strains[j][0] = '\0';
					}*/
					count = 0;
					for (k = 0; k < 4 - bad_bases_count[variant_sites_updated[i]]; k++)
					{
						if (MSA[j][variant_sites_updated[i]] != bad_base_char[variant_sites_updated[i]][k])
						{
							count++;
						}
					}
					if (count == (4 - bad_bases_count[variant_sites_updated[i]]))
					{
						// names_of_strains[j][0] = '\0';
						incompat_counter[j]++;
					}
				}
			}
		}
		for (i = 0; i < number_of_strains; i++)
		{
			if (incompat_counter[i] == number_of_iterations)
			{
				number_remaining--;
				number_removed++;
			}
		}
		if (number_remaining >= min_strains_remaining && number_remaining < max_strains_remaining)
		{
			printf("exiting loop. %d remaining\n", number_remaining);
			run_loop = 0;
		}
		else if (number_remaining >= max_strains_remaining)
		{
			printf("%d strains remaining. exiting...\n", number_remaining);
			exit(1);
		}
		else
		{
			printf("there are %d remaining... \n", number_remaining);
			number_of_iterations++;
		}
	}
	clock_gettime(CLOCK_MONOTONIC, &tend);
	printf("Took %.5fsec\n", ((double)tend.tv_sec + 1.0e-9 * tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9 * tstart.tv_nsec));
	for (i = 0; i < number_of_strains; i++)
	{
		if (incompat_counter[i] == number_of_iterations)
		{
			names_of_strains[i][0] = '\0';
		}
	}
	for (i = 0; i < length_of_MSA; i++)
	{
		free(bad_base_char[i]);
	}
	free(bad_base_char);
	free(bad_bases_count);
	free(variant_sites_updated);
	number_remaining = 0;
	for (i = 0; i < number_of_strains; i++)
	{
		if (names_of_strains[i][0] != '\0')
		{
			// printf("Remaining strain: %s\n",names_of_strains[i]);
			number_remaining++;
		}
	}
	printf("Number of strains remaining is %d\n", number_remaining);
	if (sam_results_out != NULL)
	{
		*sam_results_out = cached_lines;
		*num_sam_lines_out = cached_count_local;
	}
	return number_remaining;
}
