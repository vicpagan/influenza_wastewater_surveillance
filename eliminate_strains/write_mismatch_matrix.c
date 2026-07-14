#include "write_mismatch_matrix.h"

/**
 * @brief Thread worker: build mismatch-matrix rows for this thread's slice of the
 * preloaded SAM lines (global `sam_results`).
 *
 * For each read pair, counts mismatches (per strain, in the global `resize_MSA`
 * panel) at every aligned position, avoiding double-counting bases covered by
 * both mates of a pair (see the `visited` array). Formats one output row per
 * read pair as "readname\\talignment_size\\tmismatch_count...", written into
 * this thread's ThreadStruct.results_str.
 *
 * @param ptr Pointer to this thread's ThreadStruct (cast internally).
 * @return NULL always; real output is written into thread_str->results_str, not returned.
 */
void *writeMismatchMatrix_paired(void *ptr)
{
	int i, j, k;
	struct ThreadStruct *thread_str = (ThreadStruct *)ptr;
	ResultsStruct *results_str = thread_str->results_str;
	int max_sam_line_length = thread_str->max_sam_line_length;
	int length_of_MSA = thread_str->length_of_MSA;
	int number_of_strains = thread_str->number_of_strains;
	int number_of_strains_remaining = thread_str->number_of_strains_remaining;
	int thread_index = thread_str->thread_index;
	char buffer[FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar[i] = 0;
		cigar_chars[i] = '\0';
	}
	int *number_of_mismatches = (int *)malloc(number_of_strains_remaining * sizeof(int));
	int alignment_size;
	int first_in_pair = 0;
	int second_in_pair = 0;
	int first_seq_cigar[MAX_CIGAR];
	char first_seq_cigar_chars[MAX_CIGAR];
	int first_seq_cigar_char_count = 0;
	char *first_seq = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	memset(first_seq, '\0', MAX_READ_LENGTH);
	int first_seq_start_pos = 0;
	int visited[MAX_READ_LENGTH];
	int visited_place = 0;
	// memset(visited,-1,MAX_READ_LENGTH);
	for (i = 0; i < MAX_READ_LENGTH; i++)
	{
		visited[i] = -1;
	}
	int line_number = 0;
	// while( fgets(buffer,FASTA_MAXLINE,samfile) != NULL ){
	int line_count;
	char *resultsPath = (char *)malloc((max_sam_line_length + 300000) * sizeof(char));
	int index_mismatch = 0;
	char *context = NULL;
	for (line_count = thread_str->sam_line_start; line_count < thread_str->sam_line_end; line_count++)
	{
		line_number++;
		char *buffer_copy = strdup(sam_results[line_count]);
		s = strtok_r(sam_results[line_count], "\t", &context);
		char *name = strdup(s);
		s = strtok_r(NULL, "\t", &context);
		int decimal = 0;
		sscanf(s, "%d", &decimal);
		decimal = dec2bin(decimal);
		if (decimal == 1)
		{
			for (i = 0; i < max_sam_line_length + 300000; i++)
			{
				resultsPath[i] = '\0';
			}
			strcpy(resultsPath, name);
			strcat(resultsPath, "\t");
			first_in_pair = 1;
		}
		else if (decimal == 0)
		{
			second_in_pair = 1;
		}
		if (decimal == 2)
		{
			for (i = 0; i < max_sam_line_length + 300000; i++)
			{
				resultsPath[i] = '\0';
			}
			strcpy(resultsPath, name);
			strcat(resultsPath, "\t");
		}
		free(name);
		for (i = 0; i < 2; i++)
		{
			s = strtok_r(NULL, "\t", &context);
		}
		int position = 0;
		sscanf(s, "%d", &position);
		position--;
		if (first_in_pair == 1)
		{
			first_seq_start_pos = position;
		}
		s = strtok_r(NULL, "\t", &context);
		s = strtok_r(NULL, "\t", &context);
		char *cigar_string;
		cigar_string = strdup(s);
		char *copy = strdup(cigar_string);
		char *res = strtok_r(cigar_string, "MID", &context);
		int index = 0;
		while (res)
		{
			int from = res - cigar_string + strlen(res);
			int cigar_count = 0;
			sscanf(res, "%d", &cigar_count);
			res = strtok_r(NULL, "MID", &context);
			int to = res != NULL ? res - cigar_string : strlen(copy);
			char cigar_char = '\0';
			sscanf(copy + from, "%c", &cigar_char);
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
		free(copy);
		free(cigar_string);
		s = strtok_r(buffer_copy, "\t", &context);
		for (i = 0; i < 9; i++)
		{
			s = strtok_r(NULL, "\t", &context);
		}
		char *sequence = s;
		if (first_in_pair == 1)
		{
			strcpy(first_seq, sequence);
			first_seq_cigar_char_count = index;
		}
		int cigar_char_count = index;
		index = 0;
		int start = 0;
		int start_ref = 0;
		if (decimal == 1 || decimal == 2)
		{
			alignment_size = 0;
		}
		if (decimal == 1 || decimal == 2)
		{
			for (i = 0; i < number_of_strains_remaining; i++)
			{
				number_of_mismatches[i] = 0;
			}
		}
		int l = 0;
		int m = 0;
		if (decimal != -1)
		{
			visited_place = 0;
			if (second_in_pair == 1)
			{
				int start_ref1 = 0;
				int start1 = 0;
				for (i = 0; i < number_of_strains_remaining; i++)
				{
					visited_place = 0;
					start1 = 0;
					start_ref1 = 0;
					for (j = 0; j < first_seq_cigar_char_count; j++)
					{
						for (k = 0; k < first_seq_cigar[j]; k++)
						{
							int start2 = 0;
							int start_ref2 = 0;
							if (start_ref1 + first_seq_start_pos + k >= position)
							{
								for (l = 0; l < cigar_char_count; l++)
								{
									int position_in_MSA = reference_index[start_ref1 + first_seq_start_pos + k];
									for (m = 0; m < cigar[l]; m++)
									{
										if (start_ref1 + first_seq_start_pos + k == start_ref2 + position + m && first_seq_cigar_chars[j] == 'M')
										{
											if (cigar_chars[l] == 'M')
											{
												if (first_seq[start1 + k] != sequence[start2 + m])
												{
													if (first_seq[start1 + k] == 'A' || first_seq[start1 + k] == 'a')
													{
														if (resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
														{
															if (position_in_MSA < length_of_MSA)
															{
																number_of_mismatches[i]--;
																assert(number_of_mismatches[i] >= 0);
															}
														}
													}
													else if (first_seq[start1 + k] == 'G' || first_seq[start1 + k] == 'g')
													{
														if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
														{
															if (position_in_MSA < length_of_MSA)
															{
																number_of_mismatches[i]--;
																assert(number_of_mismatches[i] >= 0);
															}
														}
													}
													else if (first_seq[start1 + k] == 'C' || first_seq[start1 + k] == 'c')
													{
														if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
														{
															if (position_in_MSA < length_of_MSA)
															{
																number_of_mismatches[i]--;
																assert(number_of_mismatches[i] >= 0);
															}
														}
													}
													else if (first_seq[start1 + k] == 'T' || first_seq[start1 + k] == 't')
													{
														if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
														{
															if (position_in_MSA < length_of_MSA)
															{
																number_of_mismatches[i]--;
																assert(number_of_mismatches[i] >= 0);
															}
														}
													}
												}
											}
											/*if (cigar_chars[l]== 'D'){
												if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0'){
													if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
														number_of_mismatches[i]--;
													}
												}
											}*/
											/*if (cigar_chars[l] == 'I'){
												if (first_seq[k+start1] == 'A' || first_seq[k+start1] == 'a'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'A' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}else if ( first_seq[start1+k]=='G' || first_seq[start1+k]=='g'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'G' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}else if ( first_seq[start1+k]=='C' || first_seq[start1+k]=='c'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'C' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}else if ( first_seq[start1+k]=='T' || first_seq[start1+k]=='t'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'T' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}
											}*/
											if (cigar_chars[l] == 'M')
											{
												visited[visited_place] = start_ref1 + first_seq_start_pos + k;
												visited_place++;
											}
										}
									}
									//}
									if (cigar_chars[l] == 'M')
									{
										start2 = start2 + cigar[l];
										start_ref2 = start_ref2 + cigar[l];
									}
									if (cigar_chars[l] == 'I')
									{
										start2 = start2 + cigar[l];
									}
									if (cigar_chars[l] == 'D')
									{
										start_ref2 = start_ref2 + cigar[l];
									}
								}
							}
							// if ( first_seq_cigar_chars[j] != 'I' ){
							//	start_ref1 = start_ref1+first_seq_cigar[j];
							//	start1 = start1 + first_seq_cigar[j];
							// }
						}
						if (first_seq_cigar_chars[j] == 'M')
						{
							start1 = start1 + first_seq_cigar[j];
							start_ref1 = start_ref1 + first_seq_cigar[j];
						}
						if (first_seq_cigar_chars[j] == 'I')
						{
							start1 = start1 + first_seq_cigar[j];
						}
						if (first_seq_cigar_chars[j] == 'D')
						{
							start_ref1 = first_seq_cigar[j] + start_ref1;
						}
					}
				}
			}
			for (i = 0; i < number_of_strains_remaining; i++)
			{
				start = 0;
				start_ref = 0;
				for (j = 0; j < cigar_char_count; j++)
				{
					for (k = 0; k < cigar[j]; k++)
					{
						int skip = 0;
						int position_in_MSA = reference_index[k + position + start_ref];
						for (l = 0; l < visited_place; l++)
						{
							if (position_in_MSA == visited[l])
							{
								skip = 1;
							}
						}
						if (cigar_chars[j] == 'M')
						{
							if (sequence[k + start] == 'A' || sequence[k + start] == 'a')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'G' || sequence[k + start] == 'g')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'C' || sequence[k + start] == 'c')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
								//		number_of_mismatches[i]++;
								//
								// }}
								if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'T' || sequence[k + start] == 't')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										number_of_mismatches[i]++;
									}
								}
							}
						}
						if (cigar_chars[j] == 'D')
						{
							// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
							//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
							//		number_of_mismatches[i]++;
							//	}
							// }
							if (resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
							{
								if (position_in_MSA < length_of_MSA && skip == 0)
								{
									//	number_of_mismatches[i]++;
								}
							}
						}
						if (cigar_chars[j] == 'I')
						{
							if (sequence[k + start] == 'A' || sequence[k + start] == 'a')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
							if (sequence[k + start] == 'G' || sequence[k + start] == 'g')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
							if (sequence[k + start] == 'C' || sequence[k + start] == 'c')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
							if (sequence[k + start] == 'T' || sequence[k + start] == 't')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
						}
					}
					if (cigar_chars[j] == 'M')
					{
						start = cigar[j] + start;
						start_ref = cigar[j] + start_ref;
						if (i == 0)
						{
							alignment_size = alignment_size + cigar[j];
						}
					}
					if (cigar_chars[j] == 'I')
					{
						start = cigar[j] + start;
					}
					if (cigar_chars[j] == 'D')
					{
						start_ref = cigar[j] + start_ref;
					}
					// if (i==0){
					//	alignment_size = alignment_size + cigar[j];
					// }
				}
			}
		}
		if (decimal == 0 || decimal == 2)
		{
			if (visited_place > 0)
			{
				alignment_size = alignment_size - visited_place;
			}
			// fprintf(outfile,"\t%d",alignment_size);
			char *num = NULL;
			asprintf(&num, "%d", alignment_size);
			strcat(resultsPath, num);
			free(num);
			for (i = 0; i < number_of_strains_remaining; i++)
			{
				// fprintf(outfile,"\t%d",number_of_mismatches[i]);
				char *num2 = NULL;
				asprintf(&num2, "\t%d", number_of_mismatches[i]);
				strcat(resultsPath, num2);
				free(num2);
			}
			// fprintf(outfile,"\n");
			strcpy(results_str->mismatch[index_mismatch], resultsPath);
			index_mismatch++;
		}
		free(buffer_copy);
		first_in_pair = 0;
		second_in_pair = 0;
	}
	free(number_of_mismatches);
	return NULL;
}

/**
 * @brief Single-threaded, streaming counterpart to writeMismatchMatrix_paired().
 *
 * Same per-pair mismatch-counting logic, but reads SAM lines directly from
 * `samfile` one at a time (instead of the preloaded `sam_results` array) and
 * writes rows straight to `outfile` -- used for -n/--no-read-sam mode, which
 * trades speed for lower memory use.
 *
 * @param outfile Output mismatch-matrix file.
 * @param samfile Open SAM file (paired-end alignments).
 * @param MSA Strain panel (post-elimination), one row per remaining strain.
 * @param length_of_MSA Number of MSA columns.
 * @param number_of_strains Total strains before elimination (unused here; kept for signature symmetry).
 * @param number_of_strains_remaining Strains left after elimination (columns in the matrix).
 * @param names_of_strains Names of the remaining strains (written as the header row).
 * @param reference_index Maps SAM alignment position -> MSA column.
 */
void writeMismatchMatrix_paired_no_read_bam(FILE *outfile, FILE *samfile, char **MSA, int length_of_MSA, int number_of_strains, int number_of_strains_remaining, char **names_of_strains, int *reference_index)
{
	int i, j, k;
	char buffer[FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar[i] = 0;
		cigar_chars[i] = '\0';
	}
	int *number_of_mismatches = (int *)malloc(number_of_strains_remaining * sizeof(int));
	fprintf(outfile, "qName\tblockSizes");
	for (i = 0; i < number_of_strains_remaining; i++)
	{
		fprintf(outfile, "\t%s", names_of_strains[i]);
	}
	fprintf(outfile, "\n");
	int alignment_size;
	int first_in_pair = 0;
	int second_in_pair = 0;
	int first_seq_cigar[MAX_CIGAR];
	char first_seq_cigar_chars[MAX_CIGAR];
	int first_seq_cigar_char_count = 0;
	char *first_seq = (char *)malloc(MAX_READ_LENGTH * sizeof(char));
	memset(first_seq, '\0', MAX_READ_LENGTH);
	int first_seq_start_pos = 0;
	int visited[MAX_READ_LENGTH];
	int visited_place = 0;
	// memset(visited,-1,MAX_READ_LENGTH);
	for (i = 0; i < MAX_READ_LENGTH; i++)
	{
		visited[i] = -1;
	}
	int line_number = 0;
	while (fgets(buffer, FASTA_MAXLINE, samfile) != NULL)
	{
		line_number++;
		if (buffer[0] != '@')
		{
			char *buffer_copy = strdup(buffer);
			s = strtok(buffer, "\t");
			char *name = strdup(s);
			s = strtok(NULL, "\t");
			int decimal = 0;
			sscanf(s, "%d", &decimal);
			decimal = dec2bin(decimal);
			if (decimal == 1)
			{
				fprintf(outfile, "%s", name);
				first_in_pair = 1;
			}
			else if (decimal == 0)
			{
				second_in_pair = 1;
			}
			if (decimal == 2)
			{
				fprintf(outfile, "%s", name);
			}
			free(name);
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
			free(copy);
			free(cigar_string);
			s = strtok(buffer_copy, "\t");
			for (i = 0; i < 9; i++)
			{
				s = strtok(NULL, "\t");
			}
			char *sequence = s;
			if (first_in_pair == 1)
			{
				strcpy(first_seq, sequence);
				first_seq_cigar_char_count = index;
			}
			int cigar_char_count = index;
			index = 0;
			int start = 0;
			int start_ref = 0;
			if (decimal == 1 || decimal == 2)
			{
				alignment_size = 0;
			}
			if (decimal == 1 || decimal == 2)
			{
				for (i = 0; i < number_of_strains_remaining; i++)
				{
					number_of_mismatches[i] = 0;
				}
			}
			int l = 0;
			int m = 0;
			if (decimal != -1)
			{
				visited_place = 0;
				if (second_in_pair == 1)
				{
					int start_ref1 = 0;
					int start1 = 0;
					for (i = 0; i < number_of_strains_remaining; i++)
					{
						visited_place = 0;
						start1 = 0;
						start_ref1 = 0;
						for (j = 0; j < first_seq_cigar_char_count; j++)
						{
							for (k = 0; k < first_seq_cigar[j]; k++)
							{
								int start2 = 0;
								int start_ref2 = 0;
								if (start_ref1 + first_seq_start_pos + k >= position)
								{
									for (l = 0; l < cigar_char_count; l++)
									{
										int position_in_MSA = reference_index[start_ref1 + first_seq_start_pos + k];
										for (m = 0; m < cigar[l]; m++)
										{
											if (start_ref1 + first_seq_start_pos + k == start_ref2 + position + m && first_seq_cigar_chars[j] == 'M')
											{
												if (cigar_chars[l] == 'M')
												{
													if (first_seq[start1 + k] != sequence[start2 + m])
													{
														if (first_seq[start1 + k] == 'A' || first_seq[start1 + k] == 'a')
														{
															if (MSA[i][position_in_MSA] != 'A' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
															{
																if (position_in_MSA < length_of_MSA)
																{
																	number_of_mismatches[i]--;
																	assert(number_of_mismatches[i] >= 0);
																}
															}
														}
														else if (first_seq[start1 + k] == 'G' || first_seq[start1 + k] == 'g')
														{
															if (MSA[i][position_in_MSA] != 'G' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
															{
																if (position_in_MSA < length_of_MSA)
																{
																	number_of_mismatches[i]--;
																	assert(number_of_mismatches[i] >= 0);
																}
															}
														}
														else if (first_seq[start1 + k] == 'C' || first_seq[start1 + k] == 'c')
														{
															if (MSA[i][position_in_MSA] != 'C' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
															{
																if (position_in_MSA < length_of_MSA)
																{
																	number_of_mismatches[i]--;
																	assert(number_of_mismatches[i] >= 0);
																}
															}
														}
														else if (first_seq[start1 + k] == 'T' || first_seq[start1 + k] == 't')
														{
															if (MSA[i][position_in_MSA] != 'T' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
															{
																if (position_in_MSA < length_of_MSA)
																{
																	number_of_mismatches[i]--;
																	assert(number_of_mismatches[i] >= 0);
																}
															}
														}
													}
												}
												/*if (cigar_chars[l]== 'D'){
													if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0'){
														if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
															number_of_mismatches[i]--;
														}
													}
												}*/
												/*if (cigar_chars[l] == 'I'){
													if (first_seq[k+start1] == 'A' || first_seq[k+start1] == 'a'){
														if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'A' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
															if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
																number_of_mismatches[i]--;
															}
														}
													}else if ( first_seq[start1+k]=='G' || first_seq[start1+k]=='g'){
														if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'G' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
															if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
																number_of_mismatches[i]--;
															}
														}
													}else if ( first_seq[start1+k]=='C' || first_seq[start1+k]=='c'){
														if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'C' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
															if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
																number_of_mismatches[i]--;
															}
														}
													}else if ( first_seq[start1+k]=='T' || first_seq[start1+k]=='t'){
														if ( MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != 'T' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '-' && MSA[strains_kept[i]][k+first_seq_start_pos+start_ref1] != '\0' ){
															if ( k+first_seq_start_pos+start_ref1 < length_of_MSA){
																number_of_mismatches[i]--;
															}
														}
													}
												}*/
												if (cigar_chars[l] == 'M')
												{
													visited[visited_place] = start_ref1 + first_seq_start_pos + k;
													visited_place++;
												}
											}
										}
										//}
										if (cigar_chars[l] == 'M')
										{
											start2 = start2 + cigar[l];
											start_ref2 = start_ref2 + cigar[l];
										}
										if (cigar_chars[l] == 'I')
										{
											start2 = start2 + cigar[l];
										}
										if (cigar_chars[l] == 'D')
										{
											start_ref2 = start_ref2 + cigar[l];
										}
									}
								}
								// if ( first_seq_cigar_chars[j] != 'I' ){
								//	start_ref1 = start_ref1+first_seq_cigar[j];
								//	start1 = start1 + first_seq_cigar[j];
								// }
							}
							if (first_seq_cigar_chars[j] == 'M')
							{
								start1 = start1 + first_seq_cigar[j];
								start_ref1 = start_ref1 + first_seq_cigar[j];
							}
							if (first_seq_cigar_chars[j] == 'I')
							{
								start1 = start1 + first_seq_cigar[j];
							}
							if (first_seq_cigar_chars[j] == 'D')
							{
								start_ref1 = first_seq_cigar[j] + start_ref1;
							}
						}
					}
				}
				for (i = 0; i < number_of_strains_remaining; i++)
				{
					start = 0;
					start_ref = 0;
					for (j = 0; j < cigar_char_count; j++)
					{
						for (k = 0; k < cigar[j]; k++)
						{
							int skip = 0;
							int position_in_MSA = reference_index[k + position + start_ref];
							for (l = 0; l < visited_place; l++)
							{
								if (position_in_MSA == visited[l])
								{
									skip = 1;
								}
							}
							if (cigar_chars[j] == 'M')
							{
								if (sequence[k + start] == 'A' || sequence[k + start] == 'a')
								{
									// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
									//		number_of_mismatches[i]++;
									//	}
									// }
									if (MSA[i][position_in_MSA] != 'A' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
									{
										if (position_in_MSA < length_of_MSA && skip == 0)
										{
											number_of_mismatches[i]++;
										}
									}
								}
								else if (sequence[k + start] == 'G' || sequence[k + start] == 'g')
								{
									// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
									//		number_of_mismatches[i]++;
									//	}
									// }
									if (MSA[i][position_in_MSA] != 'G' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
									{
										if (position_in_MSA < length_of_MSA && skip == 0)
										{
											number_of_mismatches[i]++;
										}
									}
								}
								else if (sequence[k + start] == 'C' || sequence[k + start] == 'c')
								{
									// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
									//		number_of_mismatches[i]++;
									//
									// }}
									if (MSA[i][position_in_MSA] != 'C' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
									{
										if (position_in_MSA < length_of_MSA && skip == 0)
										{
											number_of_mismatches[i]++;
										}
									}
								}
								else if (sequence[k + start] == 'T' || sequence[k + start] == 't')
								{
									// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T' /*&& MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'*/){
									//		number_of_mismatches[i]++;
									//	}
									// }
									if (MSA[i][position_in_MSA] != 'T' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
									{
										if (position_in_MSA < length_of_MSA && skip == 0)
										{
											number_of_mismatches[i]++;
										}
									}
								}
							}
							if (cigar_chars[j] == 'D')
							{
								// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
								//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
								//		number_of_mismatches[i]++;
								//	}
								// }
								if (MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA && skip == 0)
									{
										// number_of_mismatches[i]++;
									}
								}
							}
							if (cigar_chars[j] == 'I')
							{
								if (sequence[k + start] == 'A' || sequence[k + start] == 'a')
								{
									// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A'){
									//		number_of_mismatches[i]++;
									//	}
									// }
									if (MSA[i][position_in_MSA] != 'A' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
									{
										if (position_in_MSA < length_of_MSA && skip == 0)
										{
											//	number_of_mismatches[i]++;
										}
									}
								}
								if (sequence[k + start] == 'G' || sequence[k + start] == 'g')
								{
									// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G'){
									//		number_of_mismatches[i]++;
									//	}
									// }
									if (MSA[i][position_in_MSA] != 'G' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
									{
										if (position_in_MSA < length_of_MSA && skip == 0)
										{
											//	number_of_mismatches[i]++;
										}
									}
								}
								if (sequence[k + start] == 'C' || sequence[k + start] == 'c')
								{
									// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C'){
									//		number_of_mismatches[i]++;
									//	}
									// }
									if (MSA[i][position_in_MSA] != 'C' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
									{
										if (position_in_MSA < length_of_MSA && skip == 0)
										{
											//	number_of_mismatches[i]++;
										}
									}
								}
								if (sequence[k + start] == 'T' || sequence[k + start] == 't')
								{
									// if ( reference[k+position+start_ref] < length_of_MSA && reference[k+position+start_ref] != -1){
									//	if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T'){
									//		number_of_mismatches[i]++;
									//	}
									// }
									if (MSA[i][position_in_MSA] != 'T' && MSA[i][position_in_MSA] != '-' && MSA[i][position_in_MSA] != '\0')
									{
										if (position_in_MSA < length_of_MSA && skip == 0)
										{
											//	number_of_mismatches[i]++;
										}
									}
								}
							}
						}
						if (cigar_chars[j] == 'M')
						{
							start = cigar[j] + start;
							start_ref = cigar[j] + start_ref;
							if (i == 0)
							{
								alignment_size = alignment_size + cigar[j];
							}
						}
						if (cigar_chars[j] == 'I')
						{
							start = cigar[j] + start;
						}
						if (cigar_chars[j] == 'D')
						{
							start_ref = cigar[j] + start_ref;
						}
						// if (i==0){
						//	alignment_size = alignment_size + cigar[j];
						// }
					}
				}
			}
			if (decimal == 0 || decimal == 2)
			{
				if (visited_place > 0)
				{
					alignment_size = alignment_size - visited_place;
				}
				fprintf(outfile, "\t%d", alignment_size);
				for (i = 0; i < number_of_strains_remaining; i++)
				{
					fprintf(outfile, "\t%d", number_of_mismatches[i]);
				}
				fprintf(outfile, "\n");
			}
			free(buffer_copy);
			first_in_pair = 0;
			second_in_pair = 0;
		}
	}
	free(number_of_mismatches);
}

/**
 * @brief Single-end counterpart to writeMismatchMatrix_paired_no_read_bam(): builds
 * the mismatch matrix one read at a time, with no mate-pair overlap handling needed.
 * @param outfile Output mismatch-matrix file.
 * @param samfile Open SAM file (single-end alignments).
 * @param MSA Full (pre-elimination) strain panel; indexed via strains_kept.
 * @param strains_kept Row indices into MSA of the strains that survived elimination.
 * @param length_of_MSA Number of MSA columns.
 * @param number_of_strains Total strains before elimination.
 * @param number_of_strains_remaining Strains left after elimination (columns in the matrix).
 */
void writeMismatchMatrix(FILE *outfile, FILE *samfile, char **MSA, int *strains_kept, int length_of_MSA, int number_of_strains, int number_of_strains_remaining)
{
	int i, j, k;
	char buffer[FASTA_MAXLINE];
	char *s;
	int cigar[MAX_CIGAR];
	char cigar_chars[MAX_CIGAR];
	for (i = 0; i < MAX_CIGAR; i++)
	{
		cigar[i] = 0;
		cigar_chars[i] = '\0';
	}
	int *number_of_mismatches = (int *)malloc(number_of_strains_remaining * sizeof(int));
	fprintf(outfile, "qName\tblockSizes");
	for (i = 0; i < number_of_strains_remaining; i++)
	{
		fprintf(outfile, "\t%s", resize_names_of_strains[i]);
	}
	fprintf(outfile, "\n");
	int alignment_size;
	while (fgets(buffer, FASTA_MAXLINE, samfile) != NULL)
	{
		if (buffer[0] != '@')
		{
			char *buffer_copy = strdup(buffer);
			s = strtok(buffer, "\t");
			fprintf(outfile, "%s", s);
			s = strtok(NULL, "\t");
			int decimal = 0;
			sscanf(s, "%d", &decimal);
			decimal = dec2bin(decimal);
			for (i = 0; i < 2; i++)
			{
				s = strtok(NULL, "\t");
			}
			int position = 0;
			sscanf(s, "%d", &position);
			position--;
			for (i = 0; i < 2; i++)
			{
				s = strtok(NULL, "\t");
			}
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
			char *sequence = s;
			int cigar_char_count = index;
			int start = 0;
			int start_ref = 0;
			alignment_size = 0;
			for (i = 0; i < number_of_strains_remaining; i++)
			{
				number_of_mismatches[i] = 0;
			}
			free(buffer_copy);
			/*char alignment_size [MAX_CIGAR];
			sscanf(s, "%s", &(alignment_size));
			int size = strlen(alignment_size);
			char number_to_convert [MAX_CIGAR];
			memset(number_to_convert,'\0',MAX_CIGAR);
			int placement = 0;
			int alignment_size_int = 0;
			for(i=0; i<size; i++){
				if (isalpha(alignment_size[i])){
					if (alignment_size[i]=='M'){
						alignment_size_int += atoi(number_to_convert);
						memset(number_to_convert,'\0',MAX_CIGAR);
						placement=0;
					}
				}else{
					number_to_convert[placement]=alignment_size[i];
					placement++;
				}
			}
			fprintf(outfile,"\t%d",alignment_size_int);
			for(i=0; i<4; i++){
				s = strtok(NULL,"\t");
			}
			char* sequence = s;
			size = strlen(sequence);*/
			for (i = 0; i < number_of_strains_remaining; i++)
			{
				start = 0;
				start_ref = 0;
				for (j = 0; j < cigar_char_count; j++)
				{
					for (k = 0; k < cigar[j]; k++)
					{
						int position_in_MSA = reference_index[k + position + start_ref];
						if (cigar_chars[j] == 'M')
						{
							if (sequence[k + start] == 'A' || sequence[k + start] == 'a')
							{
								// if ( MSA[strains_kept[i]][reference[k+position+start_ref]] != 'A' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-' ){
								if (resize_MSA[i][position_in_MSA] != 'A' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									// if ( reference[k+start_ref+position] < length_of_MSA ){
									if (position_in_MSA < length_of_MSA)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'G' || sequence[k + start] == 'g')
							{
								// if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'G' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
								if (resize_MSA[i][position_in_MSA] != 'G' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									if (position_in_MSA < length_of_MSA)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'C' || sequence[k + start] == 'c')
							{
								// if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'C' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
								if (resize_MSA[i][position_in_MSA] != 'C' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									// if ( reference[k+start_ref+position] < length_of_MSA ){
									if (position_in_MSA < length_of_MSA)
									{
										number_of_mismatches[i]++;
									}
								}
							}
							else if (sequence[k + start] == 'T' || sequence[k + start] == 't')
							{
								// if (MSA[strains_kept[i]][reference[k+position+start_ref]] != 'T' && MSA[strains_kept[i]][reference[k+position+start_ref]] != '-'){
								if (resize_MSA[i][position_in_MSA] != 'T' && resize_MSA[i][position_in_MSA] != '-' && resize_MSA[i][position_in_MSA] != '\0')
								{
									// if ( reference[k+start_ref+position] < length_of_MSA ){
									if (position_in_MSA < length_of_MSA)
									{
										number_of_mismatches[i]++;
									}
								}
							}
						}
					}
					if (cigar_chars[j] == 'M')
					{
						start = cigar[j] + start;
						start_ref = cigar[j] + start_ref;
						if (i == 0)
						{
							alignment_size = alignment_size + cigar[j];
						}
					}
					if (cigar_chars[j] == 'I')
					{
						start = cigar[j] + start;
					}
					if (cigar_chars[j] == 'D')
					{
						start_ref = cigar[j] + start_ref;
					}
				}
			}
			fprintf(outfile, "\t%d", alignment_size);
			for (i = 0; i < number_of_strains_remaining; i++)
			{
				fprintf(outfile, "\t%d", number_of_mismatches[i]);
			}
			fprintf(outfile, "\n");
		}
	}
	free(number_of_mismatches);
}