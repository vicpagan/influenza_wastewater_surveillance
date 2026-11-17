#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <dirent.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "select_references.h"
#include "redundant_groups.h"
#include "external/hashmap.h"

// PAF columns (1-based) that are used
#define PAF_TARGET_NAME_COLUMN 6
#define PAF_NUM_MATCHES_COLUMN 10

#define REFERENCE_FILE_EXTENSION ".fasta"

typedef HASHMAP(char, int) NameToIndexMap;

/**
 * @brief Which strains each read aligned to in one minimap2 run. The strains of read r are
 * read_strains[read_offsets[r]] ... read_strains[read_offsets[r + 1] - 1], each strain once.
 */
typedef struct ReadHits
{
	int num_reads;
	char **read_names;      // first pass only
	int *first_pass_reads;  // second pass only: the index each read had in the first pass
	int64_t *read_offsets;
	int *read_strains;
	int64_t num_hits;

	int64_t num_paf_lines;
	int64_t num_lines_to_strains_without_row;
	int64_t num_lines_to_unknown_strains;
	int num_other_reads;    // second pass only: reads with alignments that are not among the reads to cover
} ReadHits;

static int compare_names(const void *a, const void *b)
{
	return strcmp(*(char *const *)a, *(char *const *)b);
}

/**
 * @brief Lists the strains that have a <STRAIN>.fasta file in the reference sequences directory
 *
 * @param num_strains output: number of strains
 * @return char** the strain names, sorted
 */
static char **list_strains_in_reference_sequences_dir(char *reference_sequences_dir, int *num_strains)
{
	DIR *dir = opendir(reference_sequences_dir);
	if (dir == NULL)
	{
		fprintf(stderr, "Error! Cannot open reference sequences directory '%s'.\n", reference_sequences_dir);
		exit(1);
	}

	int capacity = 1 << 12;
	int count = 0;
	char **names = (char **)malloc(capacity * sizeof(char *));
	size_t extension_length = strlen(REFERENCE_FILE_EXTENSION);
	struct dirent *entry;
	while ((entry = readdir(dir)) != NULL)
	{
		size_t length = strlen(entry->d_name);
		if (length <= extension_length || strcmp(entry->d_name + length - extension_length, REFERENCE_FILE_EXTENSION) != 0)
		{
			continue;
		}
		if (count == capacity)
		{
			capacity *= 2;
			names = (char **)realloc(names, capacity * sizeof(char *));
		}
		names[count] = strndup(entry->d_name, length - extension_length);
		count++;
	}
	closedir(dir);

	qsort(names, count, sizeof(char *), compare_names);
	*num_strains = count;
	return names;
}

/**
 * @brief Writes the sequences of some strains to one FASTA file, each record named after its file
 *
 * @param strains indexes into strain_names of the strains to write, or NULL for all of them
 * @param num_strains number of strains to write
 */
static void write_sequences_fasta(char *fasta_filepath, char *reference_sequences_dir, char **strain_names, int *strains, int num_strains)
{
	int i;
	FILE *fasta_file = fopen(fasta_filepath, "w");
	if (fasta_file == NULL)
	{
		fprintf(stderr, "Error: cannot write '%s'.\n", fasta_filepath);
		exit(1);
	}

	char *buffer = (char *)malloc(1 << 16);
	char *reference_filepath = (char *)malloc(strlen(reference_sequences_dir) + 4096);
	for (i = 0; i < num_strains; i++)
	{
		char *strain_name = strain_names[strains != NULL ? strains[i] : i];
		sprintf(reference_filepath, "%s/%s%s", reference_sequences_dir, strain_name, REFERENCE_FILE_EXTENSION);
		gzFile reference_file = gzopen(reference_filepath, "r");
		if (reference_file == (gzFile)NULL)
		{
			fprintf(stderr, "Error! Cannot open file for reference sequence '%s'.\n", strain_name);
			exit(1);
		}

		// the sequence lines are copied as they are; header lines are replaced by one header, the file's name
		fprintf(fasta_file, ">%s\n", strain_name);
		int at_line_start = 1;
		int in_header = 0;
		while (gzgets(reference_file, buffer, 1 << 16) != NULL)
		{
			size_t length = strlen(buffer);
			if (at_line_start)
			{
				in_header = (buffer[0] == '>');
			}
			if (!in_header)
			{
				fputs(buffer, fasta_file);
			}
			at_line_start = (length > 0 && buffer[length - 1] == '\n');
		}
		if (!at_line_start && !in_header)
		{
			fputc('\n', fasta_file);
		}
		gzclose(reference_file);
	}
	free(buffer);
	free(reference_filepath);
	if (fclose(fasta_file) != 0)
	{
		fprintf(stderr, "Error: could not finish writing '%s'.\n", fasta_filepath);
		exit(1);
	}
}

/**
 * @brief Builds the minimap2 index of all sequences in the reference sequences directory:
 * writes them to one FASTA file, indexes that, and removes the FASTA file. Works under
 * temporary names so that an interrupted or parallel run never leaves a partial index under
 * the final name.
 */
static void build_minimap2_index(char *index_filepath, char *reference_sequences_dir, char **strain_names, int num_strains, int num_threads, int verbose)
{
	char *temporary_fasta_filepath = (char *)malloc(strlen(index_filepath) + 64);
	char *temporary_index_filepath = (char *)malloc(strlen(index_filepath) + 64);
	sprintf(temporary_fasta_filepath, "%s.%d.fasta", index_filepath, (int)getpid());
	sprintf(temporary_index_filepath, "%s.%d.tmp", index_filepath, (int)getpid());

	write_sequences_fasta(temporary_fasta_filepath, reference_sequences_dir, strain_names, NULL, num_strains);

	char *command = (char *)malloc(2 * strlen(index_filepath) + 512);
	sprintf(command, "minimap2 -x sr -t %d -d %s %s%s", num_threads, temporary_index_filepath, temporary_fasta_filepath, verbose ? "" : " > /dev/null 2>&1");
	int status = system(command);
	if (status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
	{
		remove(temporary_fasta_filepath);
		remove(temporary_index_filepath);
		fprintf(stderr, "Error: building the minimap2 index did not finish successfully%s. The command was:\n%s\n", verbose ? "" : " (run again with -v to see its messages)", command);
		exit(1);
	}
	free(command);

	if (rename(temporary_index_filepath, index_filepath) != 0)
	{
		fprintf(stderr, "Error: could not write '%s'.\n", index_filepath);
		exit(1);
	}
	remove(temporary_fasta_filepath);
	free(temporary_fasta_filepath);
	free(temporary_index_filepath);
}

/**
 * @brief Runs minimap2 and records which strains each read aligned to. minimap2 writes PAF to
 * stdout, which is read here as it is produced.
 *
 * First pass (first_pass_read_indexes == NULL): every read with an alignment is recorded, with
 * its name. Second pass: only reads that are in first_pass_read_indexes and marked in
 * read_is_wanted are recorded, each with the index it had in the first pass.
 *
 * @param minimap2_args minimap2 options
 * @param target_filepath minimap2 index or FASTA file to align to
 * @param strain_name_to_index strain name -> index of the strain
 * @param msa_row_of_strain for each strain its MSA row name, or NULL if it has none (alignments to it are ignored)
 * @param num_matches_of_strain output: matching bases over all alignments of each strain (added to)
 */
static ReadHits collect_read_hits(Options *opt, const char *minimap2_args, char *target_filepath, NameToIndexMap *strain_name_to_index, const char **msa_row_of_strain, int num_strains, double *num_matches_of_strain, NameToIndexMap *first_pass_read_indexes, char *read_is_wanted, int num_first_pass_reads)
{
	int strain_idx;
	int is_first_pass = (first_pass_read_indexes == NULL);

	char *command = (char *)malloc(16384);
	if (opt->paired)
	{
		snprintf(command, 16384, "minimap2 %s -t %d %s %s %s%s", minimap2_args, opt->num_threads, target_filepath, opt->forward_end_filepath, opt->reverse_end_filepath, opt->verbose ? "" : " 2>/dev/null");
	}
	else
	{
		snprintf(command, 16384, "minimap2 %s -t %d %s %s%s", minimap2_args, opt->num_threads, target_filepath, opt->single_end_filepath, opt->verbose ? "" : " 2>/dev/null");
	}
	printf("Running: %s\n", command);
	fflush(stdout);

	FILE *paf = popen(command, "r");
	if (paf == NULL)
	{
		fprintf(stderr, "Error: could not start minimap2.\n");
		exit(1);
	}
	setvbuf(paf, NULL, _IOFBF, 1 << 22);

	ReadHits hits;
	memset(&hits, 0, sizeof(ReadHits));
	int reads_capacity = 1 << 16;
	int64_t hits_capacity = 1 << 20;
	hits.read_offsets = (int64_t *)malloc(((size_t)reads_capacity + 1) * sizeof(int64_t));
	hits.read_strains = (int *)malloc((size_t)hits_capacity * sizeof(int));
	if (is_first_pass)
	{
		hits.read_names = (char **)malloc((size_t)reads_capacity * sizeof(char *));
	}
	else
	{
		hits.first_pass_reads = (int *)malloc((size_t)reads_capacity * sizeof(int));
	}

	// to detect a read whose alignments are not on consecutive lines
	HASHMAP(char, char) seen_read_names;                                                             // first pass
	hashmap_init(&seen_read_names, hashmap_hash_string, strcmp);
	char *read_was_seen = (char *)calloc(num_first_pass_reads > 0 ? num_first_pass_reads : 1, 1);    // second pass

	int *last_read_of_strain = (int *)malloc(num_strains * sizeof(int));  // so a read counts once per strain
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		last_read_of_strain[strain_idx] = -1;
	}

	char *line = NULL;
	size_t line_capacity = 0;
	ssize_t line_length;
	char *current_read_name = NULL;
	int read_idx = -1;  // -1: the alignments of the current read are skipped

	while ((line_length = getline(&line, &line_capacity, paf)) != -1)
	{
		if (line_length > 0 && line[line_length - 1] == '\n')
		{
			line[--line_length] = '\0';
		}
		if (line_length == 0)
		{
			continue;
		}

		// split off the columns that are needed
		char *columns[PAF_NUM_MATCHES_COLUMN];
		int num_columns = 0;
		char *field = line;
		while (num_columns < PAF_NUM_MATCHES_COLUMN)
		{
			columns[num_columns++] = field;
			char *tab = strchr(field, '\t');
			if (tab == NULL)
			{
				break;
			}
			*tab = '\0';
			field = tab + 1;
		}
		if (num_columns < PAF_NUM_MATCHES_COLUMN)
		{
			fprintf(stderr, "Error: unexpected minimap2 output (not PAF): '%s'\n", line);
			exit(1);
		}
		hits.num_paf_lines++;

		// the two mates of a pair are one read here: drop a trailing /1 or /2
		char *read_name = columns[0];
		size_t read_name_length = strlen(read_name);
		if (opt->paired && read_name_length > 2 && read_name[read_name_length - 2] == '/' && (read_name[read_name_length - 1] == '1' || read_name[read_name_length - 1] == '2'))
		{
			read_name[read_name_length - 2] = '\0';
		}

		if (current_read_name == NULL || strcmp(read_name, current_read_name) != 0)
		{
			// a new read starts
			int first_pass_read = -1;
			if (is_first_pass)
			{
				if (hashmap_get(&seen_read_names, read_name) != NULL)
				{
					fprintf(stderr, "Error: the alignments of read '%s' are not on consecutive lines of the minimap2 output, or two reads share this name.\n", read_name);
					exit(1);
				}
			}
			else
			{
				free(current_read_name);
				current_read_name = strdup(read_name);
				int *known_read = hashmap_get(first_pass_read_indexes, read_name);
				if (known_read == NULL || !read_is_wanted[*known_read])
				{
					hits.num_other_reads++;
					read_idx = -1;
					continue;
				}
				first_pass_read = *known_read;
				if (read_was_seen[first_pass_read])
				{
					fprintf(stderr, "Error: the alignments of read '%s' are not on consecutive lines of the minimap2 output.\n", read_name);
					exit(1);
				}
				read_was_seen[first_pass_read] = 1;
			}

			if (hits.num_reads == reads_capacity)
			{
				reads_capacity *= 2;
				hits.read_offsets = (int64_t *)realloc(hits.read_offsets, ((size_t)reads_capacity + 1) * sizeof(int64_t));
				if (is_first_pass)
				{
					hits.read_names = (char **)realloc(hits.read_names, (size_t)reads_capacity * sizeof(char *));
				}
				else
				{
					hits.first_pass_reads = (int *)realloc(hits.first_pass_reads, (size_t)reads_capacity * sizeof(int));
				}
			}
			read_idx = hits.num_reads++;
			hits.read_offsets[read_idx] = hits.num_hits;
			if (is_first_pass)
			{
				hits.read_names[read_idx] = strdup(read_name);
				current_read_name = hits.read_names[read_idx];
				hashmap_put(&seen_read_names, current_read_name, current_read_name);
			}
			else
			{
				hits.first_pass_reads[read_idx] = first_pass_read;
			}
		}
		if (read_idx == -1)
		{
			continue;
		}

		int *strain_index = hashmap_get(strain_name_to_index, columns[PAF_TARGET_NAME_COLUMN - 1]);
		if (strain_index == NULL)
		{
			hits.num_lines_to_unknown_strains++;
			continue;
		}
		strain_idx = *strain_index;
		if (msa_row_of_strain[strain_idx] == NULL)
		{
			hits.num_lines_to_strains_without_row++;
			continue;
		}

		num_matches_of_strain[strain_idx] += atof(columns[PAF_NUM_MATCHES_COLUMN - 1]);
		if (last_read_of_strain[strain_idx] != read_idx)
		{
			last_read_of_strain[strain_idx] = read_idx;
			if (hits.num_hits == hits_capacity)
			{
				hits_capacity *= 2;
				hits.read_strains = (int *)realloc(hits.read_strains, (size_t)hits_capacity * sizeof(int));
				if (hits.read_strains == NULL)
				{
					fprintf(stderr, "Error: out of memory while storing the minimap2 alignments (%lld read-strain combinations so far).\n", (long long)hits.num_hits);
					exit(1);
				}
			}
			hits.read_strains[hits.num_hits++] = strain_idx;
		}
	}
	free(line);
	hits.read_offsets[hits.num_reads] = hits.num_hits;
	if (!is_first_pass)
	{
		free(current_read_name);
	}

	int status = pclose(paf);
	if (status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
	{
		fprintf(stderr, "Error: minimap2 did not finish successfully%s. The command was:\n%s\n", opt->verbose ? "" : " (run again with -v to see its messages)", command);
		exit(1);
	}
	free(command);
	free(last_read_of_strain);
	free(read_was_seen);
	hashmap_cleanup(&seen_read_names);

	return hits;
}

/**
 * @brief Greedy set cover: repeatedly takes the strain that the most not-yet-covered reads
 * aligned to, until every read is covered. Ties go to the strain with the most matching
 * bases over all its alignments, then to the earlier strain.
 *
 * @param read_offsets,read_strains the strains each read aligned to (reads without strains are not counted)
 * @param chosen_strains output: the chosen strains, in the order chosen
 * @param num_reads_newly_covered output: for each chosen strain, the reads it covered that no earlier one did
 * @param num_reads_of_strain output: for each strain, the number of reads that aligned to it
 * @return int number of chosen strains
 */
static int greedy_cover(int num_strains, int num_reads, int64_t *read_offsets, int *read_strains, double *num_matches_of_strain, int **chosen_strains, int **num_reads_newly_covered, int **num_reads_of_strain)
{
	int strain_idx, read_idx;
	int64_t i, j;

	// the same information by strain: the reads of strain s are
	// strain_reads[strain_offsets[s]] ... strain_reads[strain_offsets[s + 1] - 1]
	int64_t *strain_offsets = (int64_t *)calloc((size_t)num_strains + 1, sizeof(int64_t));
	for (i = 0; i < read_offsets[num_reads]; i++)
	{
		strain_offsets[read_strains[i] + 1]++;
	}
	int *num_uncovered_reads_of_strain = (int *)malloc(num_strains * sizeof(int));
	*num_reads_of_strain = (int *)malloc(num_strains * sizeof(int));
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		num_uncovered_reads_of_strain[strain_idx] = (int)strain_offsets[strain_idx + 1];
		(*num_reads_of_strain)[strain_idx] = (int)strain_offsets[strain_idx + 1];
		strain_offsets[strain_idx + 1] += strain_offsets[strain_idx];
	}
	int *strain_reads = (int *)malloc((size_t)(strain_offsets[num_strains] > 0 ? strain_offsets[num_strains] : 1) * sizeof(int));
	int64_t *next_slot = (int64_t *)malloc(num_strains * sizeof(int64_t));
	if (strain_reads == NULL || next_slot == NULL)
	{
		fprintf(stderr, "Error: out of memory while indexing the minimap2 alignments by strain.\n");
		exit(1);
	}
	memcpy(next_slot, strain_offsets, num_strains * sizeof(int64_t));
	int num_reads_to_cover = 0;
	for (read_idx = 0; read_idx < num_reads; read_idx++)
	{
		if (read_offsets[read_idx + 1] > read_offsets[read_idx])
		{
			num_reads_to_cover++;
		}
		for (i = read_offsets[read_idx]; i < read_offsets[read_idx + 1]; i++)
		{
			strain_reads[next_slot[read_strains[i]]++] = read_idx;
		}
	}
	free(next_slot);

	char *read_is_covered = (char *)calloc(num_reads > 0 ? num_reads : 1, sizeof(char));
	int num_chosen = 0;
	int chosen_capacity = 16;
	*chosen_strains = (int *)malloc(chosen_capacity * sizeof(int));
	*num_reads_newly_covered = (int *)malloc(chosen_capacity * sizeof(int));
	int num_covered_reads = 0;

	while (num_covered_reads < num_reads_to_cover)
	{
		int best_strain = -1;
		for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
		{
			if (num_uncovered_reads_of_strain[strain_idx] <= 0)
			{
				continue;
			}
			if (best_strain == -1 || num_uncovered_reads_of_strain[strain_idx] > num_uncovered_reads_of_strain[best_strain] ||
				(num_uncovered_reads_of_strain[strain_idx] == num_uncovered_reads_of_strain[best_strain] && num_matches_of_strain[strain_idx] > num_matches_of_strain[best_strain]))
			{
				best_strain = strain_idx;
			}
		}

		if (best_strain == -1)
		{
			break;
		}

		// mark the reads of the strain as covered, and lower the counts of the other strains they aligned to
		int num_newly_covered = 0;
		for (i = strain_offsets[best_strain]; i < strain_offsets[best_strain + 1]; i++)
		{
			read_idx = strain_reads[i];
			if (read_is_covered[read_idx])
			{
				continue;
			}
			read_is_covered[read_idx] = 1;
			num_newly_covered++;
			for (j = read_offsets[read_idx]; j < read_offsets[read_idx + 1]; j++)
			{
				num_uncovered_reads_of_strain[read_strains[j]]--;
			}
		}
		num_covered_reads += num_newly_covered;

		if (num_chosen == chosen_capacity)
		{
			chosen_capacity *= 2;
			*chosen_strains = (int *)realloc(*chosen_strains, chosen_capacity * sizeof(int));
			*num_reads_newly_covered = (int *)realloc(*num_reads_newly_covered, chosen_capacity * sizeof(int));
		}
		(*chosen_strains)[num_chosen] = best_strain;
		(*num_reads_newly_covered)[num_chosen] = num_newly_covered;
		num_chosen++;
	}

	free(strain_offsets);
	free(strain_reads);
	free(num_uncovered_reads_of_strain);
	free(read_is_covered);

	return num_chosen;
}

static void print_ignored_alignments(ReadHits *hits, char *index_filepath)
{
	if (hits->num_lines_to_strains_without_row > 0)
	{
		printf("Alignments ignored because the strain has no MSA row: %lld\n", (long long)hits->num_lines_to_strains_without_row);
	}
	if (hits->num_lines_to_unknown_strains > 0)
	{
		printf("Alignments ignored because the strain has no file in the reference sequences directory: %lld. The minimap2 index is out of date; delete %s to have it rebuilt.\n", (long long)hits->num_lines_to_unknown_strains, index_filepath);
	}
}

int select_reference_strains(Options *opt, MSA *msa_str, RedundantGroups *redundant_groups, char ***selected_names)
{
	int i, strain_idx, read_idx;
	int64_t hit_idx;
	const char *read_unit = opt->paired ? "read pairs" : "reads";

	// the candidates: every strain with a sequence file
	int num_strains;
	char **strain_names = list_strains_in_reference_sequences_dir(opt->reference_sequences_dir, &num_strains);
	if (num_strains == 0)
	{
		fprintf(stderr, "Error: no %s files in reference sequences directory '%s'.\n", REFERENCE_FILE_EXTENSION, opt->reference_sequences_dir);
		exit(1);
	}

	// a strain can be a reference only if it has an MSA row: its own, or that of its group of redundant strains
	HASHMAP(char, char) msa_names;
	hashmap_init(&msa_names, hashmap_hash_string, strcmp);
	for (i = 0; i < msa_str->num_sequences; i++)
	{
		hashmap_put(&msa_names, msa_str->sequence_names[i], msa_str->sequence_names[i]);
	}

	NameToIndexMap strain_name_to_index;
	hashmap_init(&strain_name_to_index, hashmap_hash_string, strcmp);
	int *strain_indexes = (int *)malloc(num_strains * sizeof(int));
	const char **msa_row_of_strain = (const char **)malloc(num_strains * sizeof(char *));  // NULL: cannot be a reference
	int num_strains_with_own_row = 0;
	int num_strains_in_groups = 0;
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		strain_indexes[strain_idx] = strain_idx;
		hashmap_put(&strain_name_to_index, strain_names[strain_idx], &strain_indexes[strain_idx]);

		msa_row_of_strain[strain_idx] = hashmap_get(&msa_names, strain_names[strain_idx]);
		if (msa_row_of_strain[strain_idx] != NULL)
		{
			num_strains_with_own_row++;
			continue;
		}
		const char *group_name = find_group_of_strain(redundant_groups, strain_names[strain_idx]);
		if (group_name != NULL)
		{
			msa_row_of_strain[strain_idx] = hashmap_get(&msa_names, group_name);
		}
		if (msa_row_of_strain[strain_idx] != NULL)
		{
			num_strains_in_groups++;
		}
	}
	hashmap_cleanup(&msa_names);

	int num_strains_without_row = num_strains - num_strains_with_own_row - num_strains_in_groups;
	printf("Strains in the reference sequences directory: %d (%d with their own MSA row, %d in a group of redundant strains)\n", num_strains, num_strains_with_own_row, num_strains_in_groups);
	if (num_strains_without_row > 0)
	{
		printf("%d strains have no MSA row and cannot be chosen%s\n", num_strains_without_row, opt->redundant_groups_filepath[0] == '\0' ? "; if they are members of groups of redundant strains, give the groups file with -G" : "");
	}

	// the minimap2 index of all the sequences is kept next to the Bowtie2 indexes
	char *index_filepath = (char *)malloc(strlen(opt->bowtie2_indexes_dir) + 64);
	sprintf(index_filepath, "%s/all_reference_sequences.mmi", opt->bowtie2_indexes_dir);
	if (access(index_filepath, F_OK) != 0)
	{
		printf("minimap2 index of all reference sequences not built. Building now (%s)...\n", index_filepath);
		fflush(stdout);
		build_minimap2_index(index_filepath, opt->reference_sequences_dir, strain_names, num_strains, opt->num_threads, opt->verbose);
	}

	// ---- first pass: all reads against all strains, keeping only each read's best strains (minimap2's limits)
	printf("\nFirst pass: all strains, each read's best strains only\n");
	double *num_matches_of_strain = (double *)calloc(num_strains, sizeof(double));
	ReadHits first = collect_read_hits(opt, opt->minimap2_args, index_filepath, &strain_name_to_index, msa_row_of_strain, num_strains, num_matches_of_strain, NULL, NULL, 0);

	// reads that aligned to fewer than the required number of strains are left out
	char *read_is_eligible = (char *)calloc(first.num_reads > 0 ? first.num_reads : 1, sizeof(char));
	int num_eligible_reads = 0;
	int num_reads_without_usable_strain = 0;
	int64_t *eligible_read_offsets = (int64_t *)malloc(((size_t)first.num_reads + 1) * sizeof(int64_t));
	int *eligible_read_strains = (int *)malloc((size_t)(first.num_hits > 0 ? first.num_hits : 1) * sizeof(int));
	int64_t num_eligible_hits = 0;
	for (read_idx = 0; read_idx < first.num_reads; read_idx++)
	{
		int64_t num_strains_of_read = first.read_offsets[read_idx + 1] - first.read_offsets[read_idx];
		eligible_read_offsets[read_idx] = num_eligible_hits;
		if (num_strains_of_read == 0)
		{
			num_reads_without_usable_strain++;
		}
		if (num_strains_of_read > 0 && num_strains_of_read >= opt->min_strains_per_read)
		{
			read_is_eligible[read_idx] = 1;
			num_eligible_reads++;
			for (hit_idx = first.read_offsets[read_idx]; hit_idx < first.read_offsets[read_idx + 1]; hit_idx++)
			{
				eligible_read_strains[num_eligible_hits++] = first.read_strains[hit_idx];
			}
		}
	}
	eligible_read_offsets[first.num_reads] = num_eligible_hits;

	printf("Alignments reported by minimap2: %lld\n", (long long)first.num_paf_lines);
	printf("Number of %s with an alignment: %d\n", read_unit, first.num_reads);
	printf("Read-strain combinations: %lld (%.1f strains per read on average)\n", (long long)first.num_hits, first.num_reads > 0 ? (double)first.num_hits / first.num_reads : 0.0);
	print_ignored_alignments(&first, index_filepath);
	if (num_reads_without_usable_strain > 0)
	{
		printf("Number of %s aligned only to strains that were ignored: %d\n", read_unit, num_reads_without_usable_strain);
	}
	printf("Number of %s aligned to at least %d strain(s), to be covered: %d (%d left out)\n", read_unit, opt->min_strains_per_read, num_eligible_reads, first.num_reads - num_eligible_reads);

	if (num_eligible_reads == 0)
	{
		fprintf(stderr, "Error: no %s to choose reference strains from.\n", read_unit);
		exit(1);
	}

	// a provisional choice: it covers every read, but with more strains than needed, because a read that
	// matches many strains equally well was given only some of them
	int *provisional_strains, *provisional_newly_covered, *provisional_num_reads_of_strain;
	int num_provisional = greedy_cover(num_strains, first.num_reads, eligible_read_offsets, eligible_read_strains, num_matches_of_strain, &provisional_strains, &provisional_newly_covered, &provisional_num_reads_of_strain);
	printf("Provisional choice: %d strain(s)\n", num_provisional);
	free(provisional_newly_covered);
	free(provisional_num_reads_of_strain);

	// ---- second pass: all reads against the provisional strains only, with no limit on the strains per read
	printf("\nSecond pass: the %d provisional strain(s), every strain a read aligns to\n", num_provisional);
	char *provisional_fasta_filepath = (char *)malloc(strlen(opt->working_dir) + 128);
	sprintf(provisional_fasta_filepath, "%s/provisional_reference_sequences.%d.fasta", opt->working_dir, (int)getpid());
	write_sequences_fasta(provisional_fasta_filepath, opt->reference_sequences_dir, strain_names, provisional_strains, num_provisional);

	NameToIndexMap first_pass_read_indexes;
	hashmap_init(&first_pass_read_indexes, hashmap_hash_string, strcmp);
	int *read_indexes = (int *)malloc((size_t)(first.num_reads > 0 ? first.num_reads : 1) * sizeof(int));
	for (read_idx = 0; read_idx < first.num_reads; read_idx++)
	{
		read_indexes[read_idx] = read_idx;
		hashmap_put(&first_pass_read_indexes, first.read_names[read_idx], &read_indexes[read_idx]);
	}

	char *second_pass_args = (char *)malloc(strlen(opt->minimap2_args) + strlen(SECOND_PASS_MINIMAP2_ARGS) + 2);
	sprintf(second_pass_args, "%s %s", opt->minimap2_args, SECOND_PASS_MINIMAP2_ARGS);
	memset(num_matches_of_strain, 0, num_strains * sizeof(double));
	ReadHits second = collect_read_hits(opt, second_pass_args, provisional_fasta_filepath, &strain_name_to_index, msa_row_of_strain, num_strains, num_matches_of_strain, &first_pass_read_indexes, read_is_eligible, first.num_reads);
	free(second_pass_args);
	remove(provisional_fasta_filepath);
	free(provisional_fasta_filepath);

	// the final lists: for each read to cover, the provisional strains it aligned to in the second pass,
	// together with those it aligned to in the first pass (so every read keeps a strain that covers it)
	char *strain_is_provisional = (char *)calloc(num_strains, sizeof(char));
	for (i = 0; i < num_provisional; i++)
	{
		strain_is_provisional[provisional_strains[i]] = 1;
	}
	int *second_pass_position = (int *)malloc((size_t)(first.num_reads > 0 ? first.num_reads : 1) * sizeof(int));
	for (read_idx = 0; read_idx < first.num_reads; read_idx++)
	{
		second_pass_position[read_idx] = -1;
	}
	for (i = 0; i < second.num_reads; i++)
	{
		second_pass_position[second.first_pass_reads[i]] = i;
	}

	int64_t *final_read_offsets = (int64_t *)malloc(((size_t)first.num_reads + 1) * sizeof(int64_t));
	int *final_read_strains = (int *)malloc((size_t)(second.num_hits + num_eligible_hits > 0 ? second.num_hits + num_eligible_hits : 1) * sizeof(int));
	int *last_read_of_strain = (int *)malloc(num_strains * sizeof(int));
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		last_read_of_strain[strain_idx] = -1;
	}
	int64_t num_final_hits = 0;
	int num_reads_missing_from_second_pass = 0;
	for (read_idx = 0; read_idx < first.num_reads; read_idx++)
	{
		final_read_offsets[read_idx] = num_final_hits;
		if (!read_is_eligible[read_idx])
		{
			continue;
		}
		int position = second_pass_position[read_idx];
		if (position == -1)
		{
			num_reads_missing_from_second_pass++;
		}
		else
		{
			for (hit_idx = second.read_offsets[position]; hit_idx < second.read_offsets[position + 1]; hit_idx++)
			{
				strain_idx = second.read_strains[hit_idx];
				if (strain_is_provisional[strain_idx])
				{
					last_read_of_strain[strain_idx] = read_idx;
					final_read_strains[num_final_hits++] = strain_idx;
				}
			}
		}
		for (hit_idx = first.read_offsets[read_idx]; hit_idx < first.read_offsets[read_idx + 1]; hit_idx++)
		{
			strain_idx = first.read_strains[hit_idx];
			if (strain_is_provisional[strain_idx] && last_read_of_strain[strain_idx] != read_idx)
			{
				last_read_of_strain[strain_idx] = read_idx;
				final_read_strains[num_final_hits++] = strain_idx;
			}
		}
	}
	final_read_offsets[first.num_reads] = num_final_hits;

	printf("Alignments reported by minimap2: %lld\n", (long long)second.num_paf_lines);
	printf("Read-strain combinations: %lld (%.1f provisional strains per read on average)\n", (long long)num_final_hits, (double)num_final_hits / num_eligible_reads);
	print_ignored_alignments(&second, index_filepath);
	if (second.num_other_reads > 0)
	{
		printf("Number of %s with an alignment here that are not among those to be covered: %d\n", read_unit, second.num_other_reads);
	}
	if (num_reads_missing_from_second_pass > 0)
	{
		printf("Number of %s to be covered without an alignment here (their first-pass strains are used): %d\n", read_unit, num_reads_missing_from_second_pass);
	}

	int *chosen_strains, *num_reads_newly_covered, *num_reads_of_strain;
	int num_selected = greedy_cover(num_strains, first.num_reads, final_read_offsets, final_read_strains, num_matches_of_strain, &chosen_strains, &num_reads_newly_covered, &num_reads_of_strain);

	char output_filepath[4096];
	snprintf(output_filepath, sizeof(output_filepath), "%s/selected_references.tsv", opt->output_dir);
	FILE *output_file = fopen(output_filepath, "w");
	if (output_file == NULL)
	{
		fprintf(stderr, "Error: cannot write '%s'.\n", output_filepath);
		exit(1);
	}
	fprintf(output_file, "rank\tstrain\tmsa_row\treads_newly_covered\treads_covered\tfraction_covered\treads_aligned_to_strain\n");
	printf("\nrank\tstrain\tmsa_row\tnewly_covered\tcovered\tfraction\n");
	char **selected = (char **)malloc((num_selected > 0 ? num_selected : 1) * sizeof(char *));
	int num_covered_reads = 0;
	for (i = 0; i < num_selected; i++)
	{
		strain_idx = chosen_strains[i];
		num_covered_reads += num_reads_newly_covered[i];
		selected[i] = strdup(strain_names[strain_idx]);
		fprintf(output_file, "%d\t%s\t%s\t%d\t%d\t%.6f\t%d\n", i + 1, strain_names[strain_idx], msa_row_of_strain[strain_idx], num_reads_newly_covered[i], num_covered_reads, (double)num_covered_reads / num_eligible_reads, num_reads_of_strain[strain_idx]);
		printf("%d\t%s\t%s\t%d\t%d\t%.4f\n", i + 1, strain_names[strain_idx], msa_row_of_strain[strain_idx], num_reads_newly_covered[i], num_covered_reads, (double)num_covered_reads / num_eligible_reads);
	}
	fclose(output_file);

	printf("\nChose %d reference strain(s), covering %d of %d %s.\n", num_selected, num_covered_reads, num_eligible_reads, read_unit);
	printf("List written to %s\n", output_filepath);
	fflush(stdout);

	hashmap_cleanup(&first_pass_read_indexes);
	hashmap_cleanup(&strain_name_to_index);
	for (read_idx = 0; read_idx < first.num_reads; read_idx++)
	{
		free(first.read_names[read_idx]);
	}
	free(first.read_names);
	free(first.read_offsets);
	free(first.read_strains);
	free(second.first_pass_reads);
	free(second.read_offsets);
	free(second.read_strains);
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		free(strain_names[strain_idx]);
	}
	free(strain_names);
	free(strain_indexes);
	free(msa_row_of_strain);
	free(index_filepath);
	free(num_matches_of_strain);
	free(read_is_eligible);
	free(read_indexes);
	free(eligible_read_offsets);
	free(eligible_read_strains);
	free(provisional_strains);
	free(strain_is_provisional);
	free(second_pass_position);
	free(final_read_offsets);
	free(final_read_strains);
	free(last_read_of_strain);
	free(chosen_strains);
	free(num_reads_newly_covered);
	free(num_reads_of_strain);

	*selected_names = selected;
	return num_selected;
}
