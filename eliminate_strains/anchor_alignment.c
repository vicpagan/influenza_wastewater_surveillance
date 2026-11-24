#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <zlib.h>
#include <sys/types.h>
#include <sys/wait.h>

#include "anchor_alignment.h"
#include "bowtie_alignment.h"
#include "select_references.h"
#include "file_utils.h"
#include "external/hashmap.h"

typedef HASHMAP(char, int) TargetNameMap;

struct AnchorTargets
{
	int mode;
	MSA *msa_str;

	// target name -> index: the MSA row (imputed mode), or the index into positions (real mode)
	TargetNameMap name_to_index;
	int num_names;
	char **names;
	int *indexes;

	// real mode: the positions line of each target the reads aligned to
	int **positions;
	int *lengths;
};

const char *anchor_index_name(int alignment_mode)
{
	return alignment_mode == ALIGNMENT_MODE_ANCHOR_IMPUTED ? ANCHOR_IMPUTED_INDEX_NAME : ANCHOR_REAL_INDEX_NAME;
}

/**
 * @brief Name of a target as Bowtie2 writes it in a SAM file: a FASTA header up to its first whitespace
 */
static char *first_word(const char *name)
{
	return strndup(name, strcspn(name, " \t\r\n"));
}

static int index_exists(const char *index_prefix)
{
	char *filepath = (char *)malloc(strlen(index_prefix) + 16);
	sprintf(filepath, "%s.1.bt2", index_prefix);
	int exists = access(filepath, F_OK) == 0;
	if (!exists)
	{
		// large index (more than about 4 billion bases)
		sprintf(filepath, "%s.1.bt2l", index_prefix);
		exists = access(filepath, F_OK) == 0;
	}
	free(filepath);
	return exists;
}

/**
 * @brief Builds a Bowtie2 index under a temporary name and renames it when it is complete, so
 * that an interrupted build never leaves an index that looks usable
 */
static void build_bowtie2_index(const char *fasta_filepath, const char *index_prefix, int num_threads, int verbose)
{
	static const char *parts[] = {".2", ".3", ".4", ".rev.1", ".rev.2", ".1"};  // .1 last: it marks the index as present
	static const char *extensions[] = {".bt2", ".bt2l"};
	int part_idx, extension_idx;

	char *temporary_prefix = (char *)malloc(strlen(index_prefix) + 64);
	sprintf(temporary_prefix, "%s.building.%d", index_prefix, (int)getpid());

	char *command = (char *)malloc(strlen(fasta_filepath) + strlen(temporary_prefix) + 256);
	sprintf(command, "bowtie2-build%s --threads %d -f %s %s%s", verbose ? "" : " --quiet", num_threads, fasta_filepath, temporary_prefix, verbose ? "" : " > /dev/null");
	fflush(stdout);
	int status = system(command);

	char *from = (char *)malloc(strlen(temporary_prefix) + 32);
	char *to = (char *)malloc(strlen(index_prefix) + 32);
	int failed = status == -1 || !WIFEXITED(status) || WEXITSTATUS(status) != 0;
	int num_renamed = 0;
	for (part_idx = 0; part_idx < 6; part_idx++)
	{
		for (extension_idx = 0; extension_idx < 2; extension_idx++)
		{
			sprintf(from, "%s%s%s", temporary_prefix, parts[part_idx], extensions[extension_idx]);
			sprintf(to, "%s%s%s", index_prefix, parts[part_idx], extensions[extension_idx]);
			if (access(from, F_OK) != 0)
			{
				continue;
			}
			if (failed)
			{
				remove(from);
			}
			else if (rename(from, to) != 0)
			{
				fprintf(stderr, "Error: could not write '%s'.\n", to);
				exit(1);
			}
			else
			{
				num_renamed++;
			}
		}
	}
	if (failed || num_renamed != 6)
	{
		fprintf(stderr, "Error: building the Bowtie2 index did not finish successfully%s. The command was:\n%s\n", verbose ? "" : " (run again with -v to see its messages)", command);
		exit(1);
	}

	free(from);
	free(to);
	free(command);
	free(temporary_prefix);
}

/**
 * @brief Writes every MSA row without its gaps to a FASTA file. Base i of a record is the
 * row's i-th column that is not a gap; a column with an ambiguous base is written as N.
 *
 * @return int number of rows written (rows with gaps only are left out)
 */
static int write_imputed_sequences_fasta(const char *fasta_filepath, MSA *msa_str)
{
	int row_idx, column_idx;
	FILE *fasta_file = fopen(fasta_filepath, "w");
	if (fasta_file == NULL)
	{
		fprintf(stderr, "Error: cannot write '%s'.\n", fasta_filepath);
		exit(1);
	}

	int num_written = 0;
	char *sequence = (char *)malloc((size_t)msa_str->sequence_length + 2);
	for (row_idx = 0; row_idx < msa_str->num_sequences; row_idx++)
	{
		const char *row = msa_str->sequences[row_idx];
		int length = 0;
		for (column_idx = 0; column_idx < msa_str->sequence_length; column_idx++)
		{
			char base = row[column_idx];
			if (base == '-')
			{
				continue;
			}
			sequence[length++] = (base == 'A' || base == 'C' || base == 'G' || base == 'T') ? base : 'N';
		}
		if (length == 0)
		{
			continue;
		}
		sequence[length++] = '\n';

		char *name = first_word(msa_str->sequence_names[row_idx]);
		fprintf(fasta_file, ">%s\n", name);
		free(name);
		fwrite(sequence, 1, (size_t)length, fasta_file);
		num_written++;
	}
	free(sequence);

	if (fclose(fasta_file) != 0)
	{
		fprintf(stderr, "Error: could not finish writing '%s' (is the disk full?).\n", fasta_filepath);
		exit(1);
	}
	return num_written;
}

/**
 * @brief Names of the strains that have a line in the positions file, as a set
 */
static void read_names_in_positions_file(char *positions_filepath, TargetNameMap *names_set, char ***names_out, int *num_names_out)
{
	gzFile file = gzopen(positions_filepath, "r");
	if (file == (gzFile)NULL)
	{
		fprintf(stderr, "Error! Cannot open non-imputed positions file.\n");
		exit(1);
	}
	gzbuffer(file, 1 << 20);

	int capacity = 1 << 16;
	int count = 0;
	char **names = (char **)malloc((size_t)capacity * sizeof(char *));
	static int present = 1;

	char *line = NULL;
	size_t line_capacity = 0;
	while (gz_read_line(file, &line, &line_capacity))
	{
		if (line[0] != '>')
		{
			continue;
		}
		if (count == capacity)
		{
			capacity *= 2;
			names = (char **)realloc(names, (size_t)capacity * sizeof(char *));
		}
		names[count] = strdup(line + 1);
		if (hashmap_get(names_set, names[count]) == NULL)
		{
			hashmap_put(names_set, names[count], &present);
		}
		count++;
	}
	free(line);
	gzclose(file);

	*names_out = names;
	*num_names_out = count;
}

/**
 * @brief Checks the .info file of an existing index against what this run would build it
 * from, or writes it for a new index
 *
 * @param description what the index was built from, one line
 */
static void check_or_write_index_info(const char *index_prefix, const char *description, int index_is_new)
{
	char *info_filepath = (char *)malloc(strlen(index_prefix) + 16);
	sprintf(info_filepath, "%s.info", index_prefix);

	if (index_is_new)
	{
		FILE *info_file = fopen(info_filepath, "w");
		if (info_file != NULL)
		{
			fprintf(info_file, "%s\n", description);
			fclose(info_file);
		}
	}
	else
	{
		FILE *info_file = fopen(info_filepath, "r");
		if (info_file == NULL)
		{
			printf("Note: '%s' is missing, so it cannot be checked that the index was built from this MSA / these reference sequences.\n", info_filepath);
		}
		else
		{
			char found[512];
			found[0] = '\0';
			if (fgets(found, sizeof(found), info_file) != NULL)
			{
				found[strcspn(found, "\r\n")] = '\0';
			}
			fclose(info_file);
			if (strcmp(found, description) != 0)
			{
				fprintf(stderr, "Error: the index '%s' was built from different sequences.\n  index: %s\n  this run: %s\nDelete the index files (%s.*) to have it built again.\n", index_prefix, found, description, index_prefix);
				exit(1);
			}
		}
	}
	free(info_filepath);
}

void perform_anchor_alignment(Options *opt, MSA *msa_str)
{
	int i;
	const char *index_name = anchor_index_name(opt->alignment_mode);

	char *index_prefix = (char *)malloc(strlen(opt->bowtie2_indexes_dir) + strlen(index_name) + 2);
	sprintf(index_prefix, "%s/%s", opt->bowtie2_indexes_dir, index_name);
	char *sam_filepath = (char *)malloc(strlen(opt->working_dir) + strlen(index_name) + 8);
	sprintf(sam_filepath, "%s/%s.sam", opt->working_dir, index_name);
	char *fasta_filepath = (char *)malloc(strlen(opt->working_dir) + strlen(index_name) + 64);
	sprintf(fasta_filepath, "%s/%s.%d.fasta", opt->working_dir, index_name, (int)getpid());

	int index_is_new = !index_exists(index_prefix);
	char description[512];

	if (opt->alignment_mode == ALIGNMENT_MODE_ANCHOR_IMPUTED)
	{
		snprintf(description, sizeof(description), "imputed MSA rows: %d rows, %d columns", msa_str->num_sequences, msa_str->sequence_length);
		if (index_is_new)
		{
			printf("Bowtie2 index of all MSA rows not built. Building now (%s)...\n", index_prefix);
			int num_written = write_imputed_sequences_fasta(fasta_filepath, msa_str);
			printf("MSA rows in the index: %d (%d rows are gaps only and left out)\n", num_written, msa_str->num_sequences - num_written);
			build_bowtie2_index(fasta_filepath, index_prefix, opt->num_threads, opt->verbose);
			remove(fasta_filepath);
		}
	}
	else
	{
		// every strain of the reference sequences directory that has a positions line
		int num_strains_in_dir = 0;
		char **strain_names = list_strains_in_reference_sequences_dir(opt->reference_sequences_dir, &num_strains_in_dir);
		snprintf(description, sizeof(description), "non-imputed sequences: %d strains in the reference sequences directory", num_strains_in_dir);
		if (index_is_new)
		{
			printf("Bowtie2 index of all non-imputed sequences not built. Building now (%s)...\n", index_prefix);
			fflush(stdout);

			TargetNameMap names_with_positions;
			hashmap_init(&names_with_positions, hashmap_hash_string, strcmp);
			char **position_names = NULL;
			int num_position_names = 0;
			read_names_in_positions_file(opt->non_imputed_positions_filepath, &names_with_positions, &position_names, &num_position_names);

			int *strains_to_write = (int *)malloc((size_t)(num_strains_in_dir > 0 ? num_strains_in_dir : 1) * sizeof(int));
			int num_to_write = 0;
			for (i = 0; i < num_strains_in_dir; i++)
			{
				if (hashmap_get(&names_with_positions, strain_names[i]) != NULL)
				{
					strains_to_write[num_to_write++] = i;
				}
			}
			printf("Strains in the index: %d (%d strains of the reference sequences directory have no line in the positions file and are left out)\n", num_to_write, num_strains_in_dir - num_to_write);
			if (num_to_write == 0)
			{
				fprintf(stderr, "Error: no strain in '%s' has a line in the positions file.\n", opt->reference_sequences_dir);
				exit(1);
			}

			write_sequences_fasta(fasta_filepath, opt->reference_sequences_dir, strain_names, strains_to_write, num_to_write);
			build_bowtie2_index(fasta_filepath, index_prefix, opt->num_threads, opt->verbose);
			remove(fasta_filepath);

			free(strains_to_write);
			hashmap_cleanup(&names_with_positions);
			for (i = 0; i < num_position_names; i++)
			{
				free(position_names[i]);
			}
			free(position_names);
		}
		for (i = 0; i < num_strains_in_dir; i++)
		{
			free(strain_names[i]);
		}
		free(strain_names);
	}
	check_or_write_index_info(index_prefix, description, index_is_new);

	// only the best alignment of each read: no --all
	run_bowtie2(index_prefix, sam_filepath, 0, opt->bowtie2_args, opt->single_end_filepath, opt->forward_end_filepath, opt->reverse_end_filepath, opt->paired, opt->fasta_format, opt->num_threads, opt->verbose);

	free(index_prefix);
	free(sam_filepath);
	free(fasta_filepath);
}

/**
 * @brief Target (RNAME) of a SAM line
 *
 * @param length output: length of the name
 * @return const char* start of the name within the line; NULL if the line has no target ("*")
 */
static const char *sam_line_target(const char *sam_line, size_t *length)
{
	const char *field = strchr(sam_line, '\t');
	if (field != NULL)
	{
		field = strchr(field + 1, '\t');
	}
	if (field == NULL)
	{
		return NULL;
	}
	field++;
	*length = strcspn(field, "\t");
	if (*length == 0 || (*length == 1 && field[0] == '*'))
	{
		return NULL;
	}
	return field;
}

/**
 * @brief Reads one line of comma-separated positions
 */
static int *parse_positions_line(const char *line, int *num_positions)
{
	int count = 1;
	const char *cursor;
	for (cursor = line; *cursor != '\0'; cursor++)
	{
		if (*cursor == ',')
		{
			count++;
		}
	}

	int *positions = (int *)malloc((size_t)count * sizeof(int));
	int position_idx = 0;
	cursor = line;
	while (*cursor != '\0' && position_idx < count)
	{
		char *end = NULL;
		positions[position_idx++] = (int)strtol(cursor, &end, 10);
		cursor = (*end == ',') ? end + 1 : end;
		if (end == cursor && *cursor != '\0')
		{
			break;  // not a number: stop instead of looping forever
		}
	}
	*num_positions = position_idx;
	return positions;
}

AnchorTargets *prepare_anchor_targets(Options *opt, MSA *msa_str, SAMResults *sam_results_str)
{
	int i, sam_line_idx;

	AnchorTargets *targets = (AnchorTargets *)calloc(1, sizeof(AnchorTargets));
	targets->mode = opt->alignment_mode;
	targets->msa_str = msa_str;
	hashmap_init(&targets->name_to_index, hashmap_hash_string, strcmp);

	char **sam_lines = sam_results_str->sam_results[0];
	int num_sam_lines = sam_results_str->num_sam_lines;
	char *name = (char *)malloc((size_t)sam_results_str->max_sam_line_length + 2);

	if (targets->mode == ALIGNMENT_MODE_ANCHOR_IMPUTED)
	{
		// every MSA row can be a target
		targets->num_names = msa_str->num_sequences;
		targets->names = (char **)malloc((size_t)msa_str->num_sequences * sizeof(char *));
		targets->indexes = (int *)malloc((size_t)msa_str->num_sequences * sizeof(int));
		for (i = 0; i < msa_str->num_sequences; i++)
		{
			targets->names[i] = first_word(msa_str->sequence_names[i]);
			targets->indexes[i] = i;
			if (hashmap_get(&targets->name_to_index, targets->names[i]) == NULL)
			{
				hashmap_put(&targets->name_to_index, targets->names[i], &targets->indexes[i]);
			}
		}

		char *row_is_target = (char *)calloc((size_t)msa_str->num_sequences, sizeof(char));
		int num_rows_hit = 0;
		for (sam_line_idx = 0; sam_line_idx < num_sam_lines; sam_line_idx++)
		{
			size_t length = 0;
			const char *target = sam_line_target(sam_lines[sam_line_idx], &length);
			if (target == NULL)
			{
				continue;
			}
			memcpy(name, target, length);
			name[length] = '\0';
			int *row = hashmap_get(&targets->name_to_index, name);
			if (row == NULL)
			{
				fprintf(stderr, "Error: a read aligned to '%s', which is not a row of the MSA. The index '%s/%s' was built from a different MSA: delete its files to have it built again.\n", name, opt->bowtie2_indexes_dir, ANCHOR_IMPUTED_INDEX_NAME);
				exit(1);
			}
			if (!row_is_target[*row])
			{
				row_is_target[*row] = 1;
				num_rows_hit++;
			}
		}
		free(row_is_target);
		printf("MSA rows that are the best hit of at least one read: %d\n", num_rows_hit);
	}
	else
	{
		// the strains the reads aligned to
		int capacity = 1 << 10;
		targets->names = (char **)malloc((size_t)capacity * sizeof(char *));
		for (sam_line_idx = 0; sam_line_idx < num_sam_lines; sam_line_idx++)
		{
			size_t length = 0;
			const char *target = sam_line_target(sam_lines[sam_line_idx], &length);
			if (target == NULL)
			{
				continue;
			}
			memcpy(name, target, length);
			name[length] = '\0';
			if (hashmap_get(&targets->name_to_index, name) != NULL)
			{
				continue;
			}
			if (targets->num_names == capacity)
			{
				capacity *= 2;
				targets->names = (char **)realloc(targets->names, (size_t)capacity * sizeof(char *));
			}
			targets->names[targets->num_names] = strdup(name);
			// the value is set below, once the array of indexes no longer moves; any non-NULL pointer marks the name as seen
			hashmap_put(&targets->name_to_index, targets->names[targets->num_names], &targets->num_names);
			targets->num_names++;
		}

		int num_targets = targets->num_names;
		targets->indexes = (int *)malloc((size_t)(num_targets > 0 ? num_targets : 1) * sizeof(int));
		targets->positions = (int **)calloc((size_t)(num_targets > 0 ? num_targets : 1), sizeof(int *));
		targets->lengths = (int *)calloc((size_t)(num_targets > 0 ? num_targets : 1), sizeof(int));
		hashmap_clear(&targets->name_to_index);
		for (i = 0; i < num_targets; i++)
		{
			targets->indexes[i] = i;
			hashmap_put(&targets->name_to_index, targets->names[i], &targets->indexes[i]);
		}
		printf("Strains that are the best hit of at least one read: %d\n", num_targets);

		// their lines of the positions file
		gzFile positions_file = gzopen(opt->non_imputed_positions_filepath, "r");
		if (positions_file == (gzFile)NULL)
		{
			fprintf(stderr, "Error! Cannot open non-imputed positions file.\n");
			exit(1);
		}
		gzbuffer(positions_file, 1 << 20);

		char *line = NULL;
		size_t line_capacity = 0;
		int num_found = 0;
		int *wanted = NULL;
		while (num_found < num_targets && gz_read_line(positions_file, &line, &line_capacity))
		{
			if (line[0] == '>')
			{
				wanted = hashmap_get(&targets->name_to_index, line + 1);
			}
			else if (wanted != NULL)
			{
				if (targets->positions[*wanted] == NULL)
				{
					targets->positions[*wanted] = parse_positions_line(line, &targets->lengths[*wanted]);
					num_found++;
				}
				wanted = NULL;
			}
		}
		free(line);
		gzclose(positions_file);

		for (i = 0; i < num_targets; i++)
		{
			if (targets->positions[i] == NULL)
			{
				fprintf(stderr, "Error: a read aligned to strain '%s', which has no line in the positions file '%s'. The index '%s/%s' was built with a different positions file: delete its files to have it built again.\n", targets->names[i], opt->non_imputed_positions_filepath, opt->bowtie2_indexes_dir, ANCHOR_REAL_INDEX_NAME);
				exit(1);
			}
		}
	}

	free(name);
	return targets;
}

const int *anchor_target_positions(const AnchorTargets *targets, const char *target_name, int *buffer, int *length)
{
	// hashmap_get does not change the map, so threads may look names up at the same time
	int *index = hashmap_get((TargetNameMap *)&targets->name_to_index, target_name);
	if (index == NULL)
	{
		*length = 0;
		return NULL;
	}

	if (targets->mode == ALIGNMENT_MODE_ANCHOR_REAL)
	{
		*length = targets->lengths[*index];
		return targets->positions[*index];
	}

	// imputed mode: base i of the row's sequence is its i-th column that is not a gap
	const char *row = targets->msa_str->sequences[*index];
	int sequence_length = targets->msa_str->sequence_length;
	int column_idx;
	int count = 0;
	for (column_idx = 0; column_idx < sequence_length; column_idx++)
	{
		if (row[column_idx] != '-')
		{
			buffer[count++] = column_idx;
		}
	}
	*length = count;
	return buffer;
}

void free_anchor_targets(AnchorTargets *targets)
{
	int i;
	if (targets == NULL)
	{
		return;
	}
	hashmap_cleanup(&targets->name_to_index);
	for (i = 0; i < targets->num_names; i++)
	{
		free(targets->names[i]);
		if (targets->positions != NULL)
		{
			free(targets->positions[i]);
		}
	}
	free(targets->names);
	free(targets->indexes);
	free(targets->positions);
	free(targets->lengths);
	free(targets);
}
