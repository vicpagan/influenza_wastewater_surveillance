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
#include "file_utils.h"
#include "external/hashmap.h"

// PAF columns (0-based) that are used
#define PAF_QUERY_LENGTH 1
#define PAF_QUERY_START 2
#define PAF_QUERY_END 3
#define PAF_TARGET_NAME 5
#define PAF_NUM_MATCHES 9
#define PAF_ALIGNMENT_LENGTH 10
#define PAF_NUM_COLUMNS 11

#define REFERENCE_FILE_EXTENSION ".fasta"

typedef HASHMAP(char, int) NameToIndexMap;

/**
 * @brief The reads that had an alignment in any minimap2 run, by name
 */
typedef struct ReadTable
{
	int num_reads;
	int capacity;
	char **names;
	int *last_run;               // the last minimap2 run in which the read had alignments
	HASHMAP(char, void) indexes; // name -> index + 1, kept in the pointer itself
} ReadTable;

/**
 * @brief What one minimap2 run found for each read that had an alignment in it
 */
typedef struct ReadHits
{
	int num_reads;
	int *read_ids;             // index of each read in the read table
	int *num_strains_aligned;  // strains each read aligned to well enough to count as an influenza read
	int64_t *read_offsets;     // the strains that cover read r are read_strains[read_offsets[r]] ... read_strains[read_offsets[r + 1] - 1]
	int *read_strains;
	int64_t num_hits;

	int64_t num_paf_lines;
	int64_t num_lines_rejected;
	int64_t num_lines_to_strains_without_row;
	int64_t num_lines_to_unknown_strains;
} ReadHits;

static int compare_names(const void *a, const void *b)
{
	return strcmp(*(char *const *)a, *(char *const *)b);
}

static int find_read(ReadTable *table, const char *name)
{
	void *value = hashmap_get(&table->indexes, name);
	return value == NULL ? -1 : (int)((intptr_t)value - 1);
}

static int add_read(ReadTable *table, const char *name)
{
	if (table->num_reads == table->capacity)
	{
		table->capacity *= 2;
		table->names = (char **)realloc(table->names, (size_t)table->capacity * sizeof(char *));
		table->last_run = (int *)realloc(table->last_run, (size_t)table->capacity * sizeof(int));
	}
	int read_idx = table->num_reads++;
	table->names[read_idx] = strdup(name);
	table->last_run[read_idx] = 0;
	hashmap_put(&table->indexes, table->names[read_idx], (void *)((intptr_t)read_idx + 1));
	return read_idx;
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
 * @brief Runs minimap2 on the reads (opt's read files) and records, for each read, how many
 * strains it aligned to well enough to count as an influenza read, and which strains cover it.
 * minimap2 writes PAF to stdout, which is read here as it is produced.
 *
 * An alignment counts only if at least MIN_ALIGNED_FRACTION of the read is aligned. With an
 * identity of at least MIN_INFLUENZA_IDENTITY the read counts as aligned to the strain; with at
 * least opt->min_reference_identity the strain also covers the read.
 *
 * @param run_number 1, 2, ...: a different number for each run over the same read table
 * @param minimap2_args minimap2 options
 * @param target_filepath minimap2 index or FASTA file to align to
 * @param strain_name_to_index strain name -> index of the strain
 * @param msa_row_of_strain for each strain its MSA row name, or NULL if it has none (alignments to it are ignored)
 * @param num_matches_of_strain output: matching bases over the covering alignments of each strain (added to)
 * @param read_table the reads seen so far; reads new to it are added, unless read_is_wanted is given
 * @param read_is_wanted NULL, or for each read of the table whether its alignments are to be recorded
 */
static ReadHits collect_read_hits(Options *opt, int run_number, const char *minimap2_args, char *target_filepath, NameToIndexMap *strain_name_to_index, const char **msa_row_of_strain, int num_strains, double *num_matches_of_strain, ReadTable *read_table, char *read_is_wanted)
{
	int strain_idx;

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
	hits.read_ids = (int *)malloc((size_t)reads_capacity * sizeof(int));
	hits.num_strains_aligned = (int *)malloc((size_t)reads_capacity * sizeof(int));
	hits.read_offsets = (int64_t *)malloc(((size_t)reads_capacity + 1) * sizeof(int64_t));
	hits.read_strains = (int *)malloc((size_t)hits_capacity * sizeof(int));

	// so that a read counts once per strain
	int *last_read_aligned_to_strain = (int *)malloc(num_strains * sizeof(int));
	int *last_read_covered_by_strain = (int *)malloc(num_strains * sizeof(int));
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		last_read_aligned_to_strain[strain_idx] = -1;
		last_read_covered_by_strain[strain_idx] = -1;
	}

	char *line = NULL;
	size_t line_capacity = 0;
	ssize_t line_length;
	char *current_read_name = NULL;
	int read_idx = -1;  // position of the current read in this run; -1: its alignments are skipped

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
		char *columns[PAF_NUM_COLUMNS];
		int num_columns = 0;
		char *field = line;
		while (num_columns < PAF_NUM_COLUMNS)
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
		if (num_columns < PAF_NUM_COLUMNS)
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
			free(current_read_name);
			current_read_name = strdup(read_name);

			int read_id = find_read(read_table, read_name);
			if (read_id == -1 && read_is_wanted == NULL)
			{
				read_id = add_read(read_table, read_name);
			}
			if (read_id == -1 || (read_is_wanted != NULL && !read_is_wanted[read_id]))
			{
				read_idx = -1;
				continue;
			}
			if (read_table->last_run[read_id] == run_number)
			{
				fprintf(stderr, "Error: the alignments of read '%s' are not on consecutive lines of the minimap2 output, or two reads share this name.\n", read_name);
				exit(1);
			}
			read_table->last_run[read_id] = run_number;

			if (hits.num_reads == reads_capacity)
			{
				reads_capacity *= 2;
				hits.read_ids = (int *)realloc(hits.read_ids, (size_t)reads_capacity * sizeof(int));
				hits.num_strains_aligned = (int *)realloc(hits.num_strains_aligned, (size_t)reads_capacity * sizeof(int));
				hits.read_offsets = (int64_t *)realloc(hits.read_offsets, ((size_t)reads_capacity + 1) * sizeof(int64_t));
			}
			read_idx = hits.num_reads++;
			hits.read_ids[read_idx] = read_id;
			hits.num_strains_aligned[read_idx] = 0;
			hits.read_offsets[read_idx] = hits.num_hits;
		}
		if (read_idx == -1)
		{
			continue;
		}

		// how much of the read is aligned, and how well
		double read_length = atof(columns[PAF_QUERY_LENGTH]);
		double aligned_length = atof(columns[PAF_QUERY_END]) - atof(columns[PAF_QUERY_START]);
		double num_matches = atof(columns[PAF_NUM_MATCHES]);
		double alignment_length = atof(columns[PAF_ALIGNMENT_LENGTH]);
		double identity = alignment_length > 0 ? num_matches / alignment_length : 0;
		if (read_length <= 0 || aligned_length < MIN_ALIGNED_FRACTION * read_length || identity < MIN_INFLUENZA_IDENTITY)
		{
			hits.num_lines_rejected++;
			continue;
		}

		int *strain_index = hashmap_get(strain_name_to_index, columns[PAF_TARGET_NAME]);
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

		if (last_read_aligned_to_strain[strain_idx] != read_idx)
		{
			last_read_aligned_to_strain[strain_idx] = read_idx;
			hits.num_strains_aligned[read_idx]++;
		}
		if (identity < opt->min_reference_identity)
		{
			continue;
		}

		num_matches_of_strain[strain_idx] += num_matches;
		if (last_read_covered_by_strain[strain_idx] != read_idx)
		{
			last_read_covered_by_strain[strain_idx] = read_idx;
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
	free(current_read_name);
	hits.read_offsets[hits.num_reads] = hits.num_hits;

	int status = pclose(paf);
	if (status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0)
	{
		fprintf(stderr, "Error: minimap2 did not finish successfully%s. The command was:\n%s\n", opt->verbose ? "" : " (run again with -v to see its messages)", command);
		exit(1);
	}
	free(command);
	free(last_read_aligned_to_strain);
	free(last_read_covered_by_strain);

	return hits;
}

static void free_read_hits(ReadHits *hits)
{
	free(hits->read_ids);
	free(hits->num_strains_aligned);
	free(hits->read_offsets);
	free(hits->read_strains);
}

static void print_run_summary(ReadHits *hits, const char *read_unit, char *index_filepath)
{
	printf("Alignments reported by minimap2: %lld (%lld left out: under %.0f%% of the read aligned, or identity under %.2f)\n", (long long)hits->num_paf_lines, (long long)hits->num_lines_rejected, 100 * MIN_ALIGNED_FRACTION, MIN_INFLUENZA_IDENTITY);
	if (hits->num_lines_to_strains_without_row > 0)
	{
		printf("Alignments ignored because the strain has no MSA row: %lld\n", (long long)hits->num_lines_to_strains_without_row);
	}
	if (hits->num_lines_to_unknown_strains > 0)
	{
		printf("Alignments ignored because the strain has no file in the reference sequences directory: %lld. The minimap2 index is out of date; delete %s to have it rebuilt.\n", (long long)hits->num_lines_to_unknown_strains, index_filepath);
	}
}

/**
 * @brief Greedy set cover: repeatedly takes the strain that covers the most not-yet-covered
 * reads, until every read is covered or the best strain would cover fewer than
 * min_newly_covered more. Ties go to the strain with the most matching bases over all its
 * alignments, then to the earlier strain.
 *
 * @param read_offsets,read_strains the strains that cover each read (reads without strains are not counted)
 * @param chosen_strains output: the chosen strains, in the order chosen
 * @param num_reads_newly_covered output: for each chosen strain, the reads it covers that no earlier one does
 * @param num_reads_of_strain output: for each strain, the number of reads it covers
 * @return int number of chosen strains
 */
static int greedy_cover(int num_strains, int num_reads, int64_t *read_offsets, int *read_strains, double *num_matches_of_strain, int min_newly_covered, int **chosen_strains, int **num_reads_newly_covered, int **num_reads_of_strain)
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
		if (best_strain == -1 || num_uncovered_reads_of_strain[best_strain] < min_newly_covered)
		{
			break;
		}

		// mark the reads of the strain as covered, and lower the counts of the other strains that cover them
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

/**
 * @brief Copies the reads that are marked as kept from one read file to another
 *
 * @return int number of reads written
 */
static int write_kept_reads(char *input_filepath, char *output_filepath, ReadTable *read_table, char *read_is_kept, int paired, int fasta_format)
{
	gzFile input_file = gzopen(input_filepath, "r");
	if (input_file == (gzFile)NULL)
	{
		fprintf(stderr, "Error! Cannot open read file '%s'.\n", input_filepath);
		exit(1);
	}
	FILE *output_file = fopen(output_filepath, "w");
	if (output_file == NULL)
	{
		fprintf(stderr, "Error: cannot write '%s'.\n", output_filepath);
		exit(1);
	}

	char *line = NULL;
	size_t line_capacity = 0;
	char *name = NULL;
	size_t name_capacity = 0;
	int num_written = 0;
	int keep = 0;
	int64_t line_in_record = 0;  // FASTQ: 0 header, 1 sequence, 2 plus line, 3 qualities

	while (gz_read_line(input_file, &line, &line_capacity))
	{
		int is_header = fasta_format ? (line[0] == '>') : (line_in_record == 0);
		if (is_header)
		{
			// the read's name: the header up to the first blank, without a trailing /1 or /2 for pairs
			size_t length = strcspn(line + 1, " \t");
			if (length + 1 > name_capacity)
			{
				name_capacity = 2 * (length + 1);
				name = (char *)realloc(name, name_capacity);
			}
			memcpy(name, line + 1, length);
			name[length] = '\0';
			if (paired && length > 2 && name[length - 2] == '/' && (name[length - 1] == '1' || name[length - 1] == '2'))
			{
				name[length - 2] = '\0';
			}
			int read_id = find_read(read_table, name);
			keep = (read_id != -1 && read_is_kept[read_id]);
			num_written += keep;
		}
		if (keep)
		{
			fputs(line, output_file);
			fputc('\n', output_file);
		}
		if (!fasta_format)
		{
			line_in_record = (line_in_record + 1) % 4;
		}
	}
	free(line);
	free(name);
	gzclose(input_file);
	if (fclose(output_file) != 0)
	{
		fprintf(stderr, "Error: could not finish writing '%s'.\n", output_filepath);
		exit(1);
	}

	return num_written;
}

int select_reference_strains(Options *opt, MSA *msa_str, RedundantGroups *redundant_groups, char ***selected_names)
{
	int i, strain_idx, read_id;
	int64_t hit_idx;
	const char *read_unit = opt->paired ? "read pairs" : "reads";

	// every strain with a sequence file
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
		printf("%d strains have no MSA row and are ignored%s\n", num_strains_without_row, opt->redundant_groups_filepath[0] == '\0' ? "; if they are members of groups of redundant strains, give the groups file with -G" : "");
	}

	// the reference strains given with -r
	int num_given = opt->num_references;
	int *given_strains = (int *)malloc((num_given > 0 ? num_given : 1) * sizeof(int));
	char *strain_is_given = (char *)calloc(num_strains, sizeof(char));
	for (i = 0; i < num_given; i++)
	{
		int *strain_index = hashmap_get(&strain_name_to_index, opt->reference_strain_names[i]);
		if (strain_index == NULL)
		{
			fprintf(stderr, "Error: reference strain '%s' has no file '%s/%s%s'.\n", opt->reference_strain_names[i], opt->reference_sequences_dir, opt->reference_strain_names[i], REFERENCE_FILE_EXTENSION);
			exit(1);
		}
		if (msa_row_of_strain[*strain_index] == NULL)
		{
			fprintf(stderr, "Error: reference strain '%s' has no MSA row%s.\n", opt->reference_strain_names[i], opt->redundant_groups_filepath[0] == '\0' ? " (if it is a member of a group of redundant strains, give the groups file with -G)" : "");
			exit(1);
		}
		if (strain_is_given[*strain_index])
		{
			fprintf(stderr, "Error: reference strain '%s' was given twice.\n", opt->reference_strain_names[i]);
			exit(1);
		}
		given_strains[i] = *strain_index;
		strain_is_given[*strain_index] = 1;
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

	ReadTable read_table;
	read_table.num_reads = 0;
	read_table.capacity = 1 << 16;
	read_table.names = (char **)malloc((size_t)read_table.capacity * sizeof(char *));
	read_table.last_run = (int *)malloc((size_t)read_table.capacity * sizeof(int));
	hashmap_init(&read_table.indexes, hashmap_hash_string, strcmp);

	char *uncapped_args = (char *)malloc(strlen(opt->minimap2_args) + strlen(UNCAPPED_MINIMAP2_ARGS) + 2);
	sprintf(uncapped_args, "%s %s", opt->minimap2_args, UNCAPPED_MINIMAP2_ARGS);
	char *temporary_fasta_filepath = (char *)malloc(strlen(opt->working_dir) + 128);

	// ---- all reads against all strains, keeping each read's best strains (minimap2's limits)
	printf("\nAll reads against all strains (each read's best strains only)\n");
	double *num_matches_all = (double *)calloc(num_strains, sizeof(double));
	ReadHits all = collect_read_hits(opt, 1, opt->minimap2_args, index_filepath, &strain_name_to_index, msa_row_of_strain, num_strains, num_matches_all, &read_table, NULL);
	print_run_summary(&all, read_unit, index_filepath);

	// ---- all reads against the given reference strains only, with no limits: which reads do they cover?
	ReadHits given;
	memset(&given, 0, sizeof(ReadHits));
	double *num_matches_given = (double *)calloc(num_strains, sizeof(double));
	if (num_given > 0)
	{
		printf("\nAll reads against the %d given reference strain(s)\n", num_given);
		sprintf(temporary_fasta_filepath, "%s/given_reference_sequences.%d.fasta", opt->working_dir, (int)getpid());
		write_sequences_fasta(temporary_fasta_filepath, opt->reference_sequences_dir, strain_names, given_strains, num_given);
		given = collect_read_hits(opt, 2, uncapped_args, temporary_fasta_filepath, &strain_name_to_index, msa_row_of_strain, num_strains, num_matches_given, &read_table, NULL);
		remove(temporary_fasta_filepath);
		print_run_summary(&given, read_unit, index_filepath);
	}

	// where each read is in the two runs
	int num_reads = read_table.num_reads;
	int *position_in_all = (int *)malloc((size_t)(num_reads > 0 ? num_reads : 1) * sizeof(int));
	int *position_in_given = (int *)malloc((size_t)(num_reads > 0 ? num_reads : 1) * sizeof(int));
	for (read_id = 0; read_id < num_reads; read_id++)
	{
		position_in_all[read_id] = -1;
		position_in_given[read_id] = -1;
	}
	for (i = 0; i < all.num_reads; i++)
	{
		position_in_all[all.read_ids[i]] = i;
	}
	for (i = 0; i < given.num_reads; i++)
	{
		position_in_given[given.read_ids[i]] = i;
	}

	// influenza reads: aligned to enough strains, or to a given reference strain. Each is covered by the
	// first given reference strain that covers it, if any
	char *read_is_influenza = (char *)calloc(num_reads > 0 ? num_reads : 1, sizeof(char));
	int *covering_given_strain = (int *)malloc((size_t)(num_reads > 0 ? num_reads : 1) * sizeof(int));  // position in the -r list, or -1
	int *rank_of_given_strain = (int *)malloc(num_strains * sizeof(int));
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		rank_of_given_strain[strain_idx] = -1;
	}
	for (i = 0; i < num_given; i++)
	{
		rank_of_given_strain[given_strains[i]] = i;
	}
	int num_influenza_reads = 0;
	int *given_newly_covered = (int *)calloc(num_given > 0 ? num_given : 1, sizeof(int));
	int *given_num_reads = (int *)calloc(num_given > 0 ? num_given : 1, sizeof(int));
	int *last_read_of_strain = (int *)malloc(num_strains * sizeof(int));
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		last_read_of_strain[strain_idx] = -1;
	}
	for (read_id = 0; read_id < num_reads; read_id++)
	{
		int in_all = position_in_all[read_id];
		int in_given = position_in_given[read_id];
		covering_given_strain[read_id] = -1;
		if ((in_all != -1 && all.num_strains_aligned[in_all] >= opt->min_strains_per_read) || (in_given != -1 && given.num_strains_aligned[in_given] > 0))
		{
			read_is_influenza[read_id] = 1;
			num_influenza_reads++;
		}
		if (!read_is_influenza[read_id])
		{
			continue;
		}

		// the given strains that cover the read, in either run
		for (int run = 0; run < 2; run++)
		{
			ReadHits *hits = run == 0 ? &given : &all;
			int position = run == 0 ? in_given : in_all;
			if (position == -1)
			{
				continue;
			}
			for (hit_idx = hits->read_offsets[position]; hit_idx < hits->read_offsets[position + 1]; hit_idx++)
			{
				strain_idx = hits->read_strains[hit_idx];
				int rank = rank_of_given_strain[strain_idx];
				if (rank == -1 || last_read_of_strain[strain_idx] == read_id)
				{
					continue;
				}
				last_read_of_strain[strain_idx] = read_id;
				given_num_reads[rank]++;
				if (covering_given_strain[read_id] == -1 || rank < covering_given_strain[read_id])
				{
					covering_given_strain[read_id] = rank;
				}
			}
		}
		if (covering_given_strain[read_id] != -1)
		{
			given_newly_covered[covering_given_strain[read_id]]++;
		}
	}

	if (num_influenza_reads == 0)
	{
		fprintf(stderr, "Error: no influenza %s found (none aligned well enough to a strain).\n", read_unit);
		exit(1);
	}

	// ---- the influenza reads are written to the working directory; everything from here on uses only them
	char *kept_filepaths[2];
	int num_read_files = opt->paired ? 2 : 1;
	char *input_filepaths[2] = {opt->paired ? opt->forward_end_filepath : opt->single_end_filepath, opt->reverse_end_filepath};
	int num_written[2] = {0, 0};
	for (i = 0; i < num_read_files; i++)
	{
		kept_filepaths[i] = (char *)malloc(strlen(opt->working_dir) + 64);
		if (opt->paired)
		{
			sprintf(kept_filepaths[i], "%s/influenza_reads_%d.%s", opt->working_dir, i + 1, opt->fasta_format ? "fasta" : "fastq");
		}
		else
		{
			sprintf(kept_filepaths[i], "%s/influenza_reads.%s", opt->working_dir, opt->fasta_format ? "fasta" : "fastq");
		}
		num_written[i] = write_kept_reads(input_filepaths[i], kept_filepaths[i], &read_table, read_is_influenza, opt->paired, opt->fasta_format);
	}
	if (opt->paired && num_written[0] != num_written[1])
	{
		fprintf(stderr, "Error: %d influenza reads found in '%s' but %d in '%s'; the two files must hold the same reads, in the same order.\n", num_written[0], opt->forward_end_filepath, num_written[1], opt->reverse_end_filepath);
		exit(1);
	}
	for (i = 0; i < num_read_files; i++)
	{
		if (strlen(kept_filepaths[i]) >= sizeof(opt->single_end_filepath))
		{
			fprintf(stderr, "Error: working directory path is too long.\n");
			exit(1);
		}
	}
	if (opt->paired)
	{
		strcpy(opt->forward_end_filepath, kept_filepaths[0]);
		strcpy(opt->reverse_end_filepath, kept_filepaths[1]);
	}
	else
	{
		strcpy(opt->single_end_filepath, kept_filepaths[0]);
	}

	printf("\nInfluenza %s (aligned to at least %d strain(s)%s): %d, written to %s%s%s\n", read_unit, opt->min_strains_per_read, num_given > 0 ? ", or to a given reference strain" : "", num_influenza_reads, kept_filepaths[0], opt->paired ? " and " : "", opt->paired ? kept_filepaths[1] : "");
	if (num_written[0] != num_influenza_reads)
	{
		printf("Warning: %d of them were found in the read file(s) by name.\n", num_written[0]);
	}

	// ---- reference strains to add: for the influenza reads that no given reference strain covers.
	// A provisional choice from each read's best strains; it covers every read that can be covered, but
	// with more strains than needed, because a read that matches many strains equally well was given
	// only some of them
	char *read_needs_reference = (char *)calloc(num_reads > 0 ? num_reads : 1, sizeof(char));
	int64_t *wanted_offsets = (int64_t *)malloc(((size_t)num_reads + 1) * sizeof(int64_t));
	int *wanted_strains = (int *)malloc((size_t)(all.num_hits > 0 ? all.num_hits : 1) * sizeof(int));
	int64_t num_wanted_hits = 0;
	int num_reads_covered_by_given = 0;
	int num_reads_needing_reference = 0;
	for (read_id = 0; read_id < num_reads; read_id++)
	{
		wanted_offsets[read_id] = num_wanted_hits;
		if (!read_is_influenza[read_id])
		{
			continue;
		}
		if (covering_given_strain[read_id] != -1)
		{
			num_reads_covered_by_given++;
			continue;
		}
		int position = position_in_all[read_id];
		if (position == -1 || all.read_offsets[position + 1] == all.read_offsets[position])
		{
			continue;  // no strain covers it
		}
		read_needs_reference[read_id] = 1;
		num_reads_needing_reference++;
		for (hit_idx = all.read_offsets[position]; hit_idx < all.read_offsets[position + 1]; hit_idx++)
		{
			wanted_strains[num_wanted_hits++] = all.read_strains[hit_idx];
		}
	}
	wanted_offsets[num_reads] = num_wanted_hits;
	int num_reads_no_strain_covers = num_influenza_reads - num_reads_covered_by_given - num_reads_needing_reference;

	if (num_given > 0)
	{
		printf("Covered by a given reference strain (at least %.0f%% of the read aligned, identity at least %.2f): %d\n", 100 * MIN_ALIGNED_FRACTION, opt->min_reference_identity, num_reads_covered_by_given);
	}
	printf("Not covered by a given reference strain, but by other strains: %d\n", num_reads_needing_reference);
	printf("Not covered by any strain (identity under %.2f to all): %d\n", opt->min_reference_identity, num_reads_no_strain_covers);

	int num_added = 0;
	int *added_strains = NULL, *added_newly_covered = NULL, *added_num_reads = NULL;
	if (num_reads_needing_reference > 0)
	{
		int *provisional_strains, *provisional_newly_covered, *provisional_num_reads;
		int num_provisional = greedy_cover(num_strains, num_reads, wanted_offsets, wanted_strains, num_matches_all, 1, &provisional_strains, &provisional_newly_covered, &provisional_num_reads);
		free(provisional_newly_covered);
		free(provisional_num_reads);

		// the influenza reads against the provisional strains only, with no limit on the strains per read
		printf("\nInfluenza reads against the %d provisional strain(s) (every strain a read aligns to)\n", num_provisional);
		sprintf(temporary_fasta_filepath, "%s/provisional_reference_sequences.%d.fasta", opt->working_dir, (int)getpid());
		write_sequences_fasta(temporary_fasta_filepath, opt->reference_sequences_dir, strain_names, provisional_strains, num_provisional);
		double *num_matches_provisional = (double *)calloc(num_strains, sizeof(double));
		ReadHits provisional = collect_read_hits(opt, 3, uncapped_args, temporary_fasta_filepath, &strain_name_to_index, msa_row_of_strain, num_strains, num_matches_provisional, &read_table, read_needs_reference);
		remove(temporary_fasta_filepath);
		print_run_summary(&provisional, read_unit, index_filepath);

		// the final lists: for each read, the provisional strains that cover it in this run, together
		// with those that covered it in the first run
		char *strain_is_provisional = (char *)calloc(num_strains, sizeof(char));
		for (i = 0; i < num_provisional; i++)
		{
			strain_is_provisional[provisional_strains[i]] = 1;
		}
		int *position_in_provisional = (int *)malloc((size_t)num_reads * sizeof(int));
		for (read_id = 0; read_id < num_reads; read_id++)
		{
			position_in_provisional[read_id] = -1;
		}
		for (i = 0; i < provisional.num_reads; i++)
		{
			position_in_provisional[provisional.read_ids[i]] = i;
		}
		for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
		{
			last_read_of_strain[strain_idx] = -1;
		}
		int64_t *final_offsets = (int64_t *)malloc(((size_t)num_reads + 1) * sizeof(int64_t));
		int *final_strains = (int *)malloc((size_t)(provisional.num_hits + num_wanted_hits > 0 ? provisional.num_hits + num_wanted_hits : 1) * sizeof(int));
		int64_t num_final_hits = 0;
		for (read_id = 0; read_id < num_reads; read_id++)
		{
			final_offsets[read_id] = num_final_hits;
			if (!read_needs_reference[read_id])
			{
				continue;
			}
			int position = position_in_provisional[read_id];
			if (position != -1)
			{
				for (hit_idx = provisional.read_offsets[position]; hit_idx < provisional.read_offsets[position + 1]; hit_idx++)
				{
					strain_idx = provisional.read_strains[hit_idx];
					if (strain_is_provisional[strain_idx])
					{
						last_read_of_strain[strain_idx] = read_id;
						final_strains[num_final_hits++] = strain_idx;
					}
				}
			}
			for (hit_idx = wanted_offsets[read_id]; hit_idx < wanted_offsets[read_id + 1]; hit_idx++)
			{
				strain_idx = wanted_strains[hit_idx];
				if (strain_is_provisional[strain_idx] && last_read_of_strain[strain_idx] != read_id)
				{
					last_read_of_strain[strain_idx] = read_id;
					final_strains[num_final_hits++] = strain_idx;
				}
			}
		}
		final_offsets[num_reads] = num_final_hits;

		num_added = greedy_cover(num_strains, num_reads, final_offsets, final_strains, num_matches_provisional, opt->min_reads_per_added_reference, &added_strains, &added_newly_covered, &added_num_reads);

		free_read_hits(&provisional);
		free(num_matches_provisional);
		free(provisional_strains);
		free(strain_is_provisional);
		free(position_in_provisional);
		free(final_offsets);
		free(final_strains);
	}

	// ---- the reference strains: the given ones, then the added ones
	int num_selected = num_given + num_added;
	if (num_selected == 0)
	{
		fprintf(stderr, "Error: no reference strains: none was given with -r, and no strain covers at least %d %s (-M).\n", opt->min_reads_per_added_reference, read_unit);
		exit(1);
	}

	char output_filepath[4096];
	snprintf(output_filepath, sizeof(output_filepath), "%s/selected_references.tsv", opt->output_dir);
	FILE *output_file = fopen(output_filepath, "w");
	if (output_file == NULL)
	{
		fprintf(stderr, "Error: cannot write '%s'.\n", output_filepath);
		exit(1);
	}
	fprintf(output_file, "rank\tstrain\tmsa_row\tchosen\treads_newly_covered\treads_covered\tfraction_covered\treads_covered_by_strain\n");
	printf("\nrank\tstrain\tmsa_row\tchosen\tnewly_covered\tcovered\tfraction\n");
	char **selected = (char **)malloc(num_selected * sizeof(char *));
	int num_covered_reads = 0;
	for (i = 0; i < num_selected; i++)
	{
		int is_given = i < num_given;
		strain_idx = is_given ? given_strains[i] : added_strains[i - num_given];
		int num_newly_covered = is_given ? given_newly_covered[i] : added_newly_covered[i - num_given];
		int num_reads_of_strain = is_given ? given_num_reads[i] : added_num_reads[strain_idx];
		num_covered_reads += num_newly_covered;
		selected[i] = strdup(strain_names[strain_idx]);
		fprintf(output_file, "%d\t%s\t%s\t%s\t%d\t%d\t%.6f\t%d\n", i + 1, strain_names[strain_idx], msa_row_of_strain[strain_idx], is_given ? "given" : "added", num_newly_covered, num_covered_reads, (double)num_covered_reads / num_influenza_reads, num_reads_of_strain);
		printf("%d\t%s\t%s\t%s\t%d\t%d\t%.4f\n", i + 1, strain_names[strain_idx], msa_row_of_strain[strain_idx], is_given ? "given" : "added", num_newly_covered, num_covered_reads, (double)num_covered_reads / num_influenza_reads);
	}
	fclose(output_file);

	printf("\n%d reference strain(s) (%d given, %d added), covering %d of %d influenza %s.\n", num_selected, num_given, num_added, num_covered_reads, num_influenza_reads, read_unit);
	printf("List written to %s\n", output_filepath);
	fflush(stdout);

	hashmap_cleanup(&read_table.indexes);
	hashmap_cleanup(&strain_name_to_index);
	for (read_id = 0; read_id < num_reads; read_id++)
	{
		free(read_table.names[read_id]);
	}
	free(read_table.names);
	free(read_table.last_run);
	free_read_hits(&all);
	if (num_given > 0)
	{
		free_read_hits(&given);
	}
	for (strain_idx = 0; strain_idx < num_strains; strain_idx++)
	{
		free(strain_names[strain_idx]);
	}
	free(strain_names);
	free(strain_indexes);
	free(msa_row_of_strain);
	free(given_strains);
	free(strain_is_given);
	free(index_filepath);
	free(uncapped_args);
	free(temporary_fasta_filepath);
	free(num_matches_all);
	free(num_matches_given);
	free(position_in_all);
	free(position_in_given);
	free(read_is_influenza);
	free(covering_given_strain);
	free(rank_of_given_strain);
	free(given_newly_covered);
	free(given_num_reads);
	free(last_read_of_strain);
	for (i = 0; i < num_read_files; i++)
	{
		free(kept_filepaths[i]);
	}
	free(read_needs_reference);
	free(wanted_offsets);
	free(wanted_strains);
	free(added_strains);
	free(added_newly_covered);
	free(added_num_reads);

	// from here on, the chosen strains are the reference strains
	for (i = 0; i < opt->num_references; i++)
	{
		free(opt->reference_strain_names[i]);
	}
	free(opt->reference_strain_names);

	*selected_names = selected;
	return num_selected;
}
