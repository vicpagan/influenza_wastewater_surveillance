#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "options.h"
#include "select_references.h"

static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"MSA-filepath", required_argument, 0, 'i'},
	{"non-imputed-positions-filepath", required_argument, 0, 'p'},
	{"reference-sequences-dir", required_argument, 0, 'g'},
	{"bowtie2-indexes-dir", required_argument, 0, 'b'},
	{"output-dir", required_argument, 0, 'o'},
	{"paired", no_argument, 0, 'P'},
	{"max-num-reads", required_argument, 0, 'm'},
	{"single-end-filepath", required_argument, 0, '0'},
	{"forward-end-filepath", required_argument, 0, '1'},
	{"reverse-end-filepath", required_argument, 0, '2'},
	{"EM-error", required_argument, 0, 'e'},
	{"fasta", no_argument, 0, 'a'},
	{"llr", no_argument, 0, 'l'},
	{"num-top-strains-llr", required_argument, 0, 'z'},
	{"cores", required_argument, 0, 't'},
	{"no-read-sam", no_argument, 0, 'n'},
	{"working-dir", required_argument, 0, 'w'},
	{"verbose", no_argument, 0, 'v'},
	{"print-mismatch-matrix", no_argument, 0, 'd'},
	{"redundant-groups-filepath", required_argument, 0, 'G'},
	{"min-strains-per-read", required_argument, 0, 'k'},
	{"minimap2-args", required_argument, 0, 'x'},
	{"min-reference-identity", required_argument, 0, 'I'},
	{"min-reads-per-added-reference", required_argument, 0, 'M'},
	{"alignment-mode", required_argument, 0, 'A'},
	{"bowtie2-args", required_argument, 0, 'B'},
	{0, 0, 0, 0}
};

// TODO: add the following line once problematic sites aspect is implemented
// -S, --problematic_sites_dir [REQUIRED,DIR]	Directory of lists of problematic sites\n
char usage[] = "\neliminate_strains [OPTIONS]\n\
	\n\
	-h, --help				\n\
	-i, --MSA-filepath [REQUIRED,FILE]		Filepath of MSA FASTA of influenza reference strains\n\
	-p, --non-imputed-positions-filepath [REQUIRED,FILE]		Filepath of the non-imputed positions for each strain within the MSA (not needed in mode 'anchor-imputed')\n\
	-g, --reference-sequences-dir [REQUIRED,DIR]		Directory of the non-imputed sequence of every strain (<STRAIN>.fasta)\n\
	-b, --bowtie2-indexes-dir [REQUIRED,DIR]		Directory of the Bowtie2 indexes (and the minimap2 index); a missing index is built\n\
	-A, --alignment-mode [MODE]		How the reads are placed on the MSA [default: references]\n\
			references	Bowtie2 to each reference strain given with -r (which is required)\n\
			minimap2	minimap2 against all strains first: only the influenza reads are kept, and reference strains are added to those given with -r (optional) for the reads they do not cover; then as 'references'\n\
			anchor-imputed	Bowtie2 to one index of all MSA rows (imputed sequences, gaps removed); each read is placed by its best hit. No reference strains: -r, -p and -g are not used\n\
			anchor-real	Bowtie2 to one index of the non-imputed sequences of all strains in -g; each read is placed by its best hit. No reference strains: -r is not used\n\
	-B, --bowtie2-args [STRING]		Extra Bowtie2 options for aligning the reads, in every mode (quote them, e.g. -B \"--score-min L,-0.6,-1.0\"). They come after the program's own options\n\
	-r, --reference-strains [NAME NAME ...]		Reference strains (e.g., -r EPI_ISL_19088566 EPI_ISL_19407907 ...). Required in mode 'references', optional in mode 'minimap2'\n\
	-G, --redundant-groups-filepath [FILE]		Redundant groups file of the MSA (\">cN\" lines, each followed by its member strains). With it, strains that are members of a group can be used as reference strains\n\
	-I, --min-reference-identity [decimal]		Mode 'minimap2': a reference strain covers a read if at least 90 percent of the read aligns to it with at least this identity [default: 0.93]\n\
	-M, --min-reads-per-added-reference [int]		Mode 'minimap2': a reference strain is added only if it covers at least this many reads that no reference strain covers yet [default: 5]\n\
	-k, --min-strains-per-read [int]		Mode 'minimap2': a read is kept as an influenza read only if it aligns to at least this many strains, or to a strain given with -r [default: 1]\n\
	-x, --minimap2-args [STRING]		Mode 'minimap2': minimap2 options for aligning the reads to all strains, replacing the defaults (quote them; keep -c) [default: \"" DEFAULT_MINIMAP2_ARGS "\"]\n\
	-o, --output-dir [DIR]		Directory to place output files [default: .]\n\
	-w, --working-dir [DIR]		Directory for intermediate/working files [default: .]\n\
	-P, --paired		Using paired-reads\n\
	-m, --max-num-reads [int]		Maximum number of reads from the read file(s) to process [default: 100,000]\n\
	-0, --single-end-filepath [FILE]		Single-end reads\n\
	-1, --forward-end-filepath [FILE]		If using paired-reads, the forward reads file\n\
	-2, --reverse-end-filepath [FILE]		If using paired-reads, the reverse reads file\n\
	-e, --EM-error [decimal]		Error rate for EM algorithm (must be in the range (0.001, 1)) [default: 0.005]\n\
	-a, --fasta		Reads are in FASTA format [default: FASTQ]\n\
	-l, --llr		Perform the LLR procedure\n\
	-z, --num-top-strains-llr [int]		Number of highest proportion strains to calculate the per-strain LLR for [default: 10]\n\
	-t, --cores [int]		Number of cores [default: 1]\n\
	-n, --no-read-sam		Don't read in sam file to memory (currently not working)\n\
	-v, --verbose		Show verbose debug output from Bowtie2 commands\n\
	-d, --print-mismatch-matrix		Print out the contents of the computed mismatch matrix to the output directory\n\
	\n";

/**
 * @brief Prints the help/usage text to CLI
 * 
 */
void print_help_statement()
{
	printf("%s", &usage[0]);
	return;
}

static void extract_reference_strains_args(int *argc_ptr, char **argv, Options *opt)
{
	int argc = *argc_ptr;

	int i;
	for (i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--reference-strains") == 0)
		{
			if (opt->reference_strain_names != NULL)
			{
				fprintf(stderr, "Error: -r/--reference-strains was given more than once; list all strain names after a single -r.\n");
				exit(1);
			}

			int flag_idx = i;
			int j = i + 1;
			int capacity = 16;
			opt->reference_strain_names = (char **)malloc(capacity * sizeof(char *));

			while (j < argc && argv[j][0] != '-')
			{
				if (opt->num_references == capacity)
				{
					capacity *= 2;
					opt->reference_strain_names = (char **)realloc(opt->reference_strain_names, capacity * sizeof(char *));
				}
				opt->reference_strain_names[opt->num_references] = strdup(argv[j]);
				opt->num_references++;
				j++;
			}

			if (opt->num_references == 0)
			{
				fprintf(stderr, "Error: -r/--reference-strains was given with no strain names after it.\n");
				exit(1);
			}

			int num_consumed = j - flag_idx;
			int k;
			for (k = flag_idx; k + num_consumed < argc; k++)
			{
				argv[k] = argv[k + num_consumed];
			}
			argc -= num_consumed;

			i--;
		}
	}

	*argc_ptr = argc;
}

/**
 * @brief Parses CLI arguments into an Options struct
 * 
 * @param argc arg count
 * @param argv arg vector
 * @param opt output options instance
 */
void parse_options(int argc, char **argv, Options *opt)
{
	int option_index, success;
	int c;
	if (argc == 1)
	{
		print_help_statement();
		exit(0);
	}

	opt->msa_filepath[0] = '\0';
	opt->non_imputed_positions_filepath[0] = '\0';
	opt->reference_sequences_dir[0] = '\0';
	opt->bowtie2_indexes_dir[0] = '\0';
	opt->single_end_filepath[0] = '\0';
	opt->forward_end_filepath[0] = '\0';
	opt->reverse_end_filepath[0] = '\0';
	opt->reference_strain_names = NULL;
	opt->num_references = 0;
	opt->redundant_groups_filepath[0] = '\0';
	opt->min_reference_identity = 0.93;
	opt->min_reads_per_added_reference = 5;
	opt->alignment_mode = ALIGNMENT_MODE_REFERENCES;
	opt->bowtie2_args[0] = '\0';

	extract_reference_strains_args(&argc, argv, opt);

	strcpy(opt->minimap2_args, DEFAULT_MINIMAP2_ARGS);
	opt->min_strains_per_read = 1;

	while (1)
	{
		c = getopt_long(argc, argv, "hi:p:g:b:o:Pm:0:1:2:e:alz:t:nw:vdG:k:x:I:M:A:B:", long_options, &option_index);
		if (c == -1)
			break;
		switch (c)
		{
		case 'h':
			print_help_statement();
			exit(0);
			break;
		case 'P':
			opt->paired = 1;
			break;
		case 'a':
			opt->fasta_format = 1;
			break;
		case 'l':
			opt->llr = 1;
			break;
		case 'g':
			success = sscanf(optarg, "%s", opt->reference_sequences_dir);
			if (!success)
				fprintf(stderr, "Invalid reference sequences directory\n");
			break;
		case 'b':
			success = sscanf(optarg, "%s", opt->bowtie2_indexes_dir);
			if (!success)
				fprintf(stderr, "Invalid bowtie2 indexes directory\n");
			break;	
		case 'n':
			opt->no_read_bam = 1;
			break;
		case 'i':
			success = sscanf(optarg, "%s", opt->msa_filepath);
			if (!success)
				fprintf(stderr, "Invalid MSA filepath\n");
			break;
		case 'p':
			success = sscanf(optarg, "%s", opt->non_imputed_positions_filepath);
			if (!success)
				fprintf(stderr, "Invalid non-imputed positions filepath\n");
			break;	
		case '0':
			success = sscanf(optarg, "%s", opt->single_end_filepath);
			if (!success)
				fprintf(stderr, "Invalid FASTA filepath\n");
			break;
		case '1':
			success = sscanf(optarg, "%s", opt->forward_end_filepath);
			if (!success)
				fprintf(stderr, "Invalid FASTA filepath\n");
			break;
		case '2':
			success = sscanf(optarg, "%s", opt->reverse_end_filepath);
			if (!success)
				fprintf(stderr, "Invalid FASTA filepath\n");
			break;
		case 't':
			success = sscanf(optarg, "%d", &(opt->num_threads));
			if (!success)
				fprintf(stderr, "Invalid number of cores\n");
			break;
		case 'm':
			success = sscanf(optarg, "%d", &(opt->max_num_reads));
			if (!success)
				fprintf(stderr, "Invalid maximum number of reads to process\n");
			break;		
		case 'e':
			success = sscanf(optarg, "%lf", &(opt->em_error));
			if (!success)
				fprintf(stderr, "Invalid error rate\n");
			break;
		case 'o':
			success = sscanf(optarg, "%s", opt->output_dir);
			if (!success)
				fprintf(stderr, "Invalid output directory\n");
			break;
		case 'w':
			success = sscanf(optarg, "%s", opt->working_dir);
			if (!success)
				fprintf(stderr, "Invalid working directory\n");
			break;
		case 'v':
			opt->verbose = 1;
			break;
		case 'd':
			opt->print_mismatch_matrix = 1;
			break;
		case 'z':
			success = sscanf(optarg, "%d", &(opt->num_top_strains_llr));
			if (!success)
				fprintf(stderr, "Invalid number of top strains to compute LLR for\n");
			break;
		case 'G':
			success = sscanf(optarg, "%2047s", opt->redundant_groups_filepath);
			if (!success)
				fprintf(stderr, "Invalid redundant groups filepath\n");
			break;
		case 'k':
			if (sscanf(optarg, "%d", &(opt->min_strains_per_read)) != 1 || opt->min_strains_per_read < 1)
			{
				fprintf(stderr, "Error: -k/--min-strains-per-read must be a whole number, 1 or more.\n");
				exit(1);
			}
			break;
		case 'I':
			if (sscanf(optarg, "%lf", &(opt->min_reference_identity)) != 1 || opt->min_reference_identity <= 0 || opt->min_reference_identity > 1)
			{
				fprintf(stderr, "Error: -I/--min-reference-identity must be a number above 0 and at most 1.\n");
				exit(1);
			}
			break;
		case 'M':
			if (sscanf(optarg, "%d", &(opt->min_reads_per_added_reference)) != 1 || opt->min_reads_per_added_reference < 1)
			{
				fprintf(stderr, "Error: -M/--min-reads-per-added-reference must be a whole number, 1 or more.\n");
				exit(1);
			}
			break;
		case 'x':
			if (strlen(optarg) >= sizeof(opt->minimap2_args))
			{
				fprintf(stderr, "Error: -x/--minimap2-args is too long.\n");
				exit(1);
			}
			strcpy(opt->minimap2_args, optarg);
			break;
		case 'A':
			if (strcmp(optarg, "references") == 0)
			{
				opt->alignment_mode = ALIGNMENT_MODE_REFERENCES;
			}
			else if (strcmp(optarg, "minimap2") == 0)
			{
				opt->alignment_mode = ALIGNMENT_MODE_MINIMAP2;
			}
			else if (strcmp(optarg, "anchor-imputed") == 0)
			{
				opt->alignment_mode = ALIGNMENT_MODE_ANCHOR_IMPUTED;
			}
			else if (strcmp(optarg, "anchor-real") == 0)
			{
				opt->alignment_mode = ALIGNMENT_MODE_ANCHOR_REAL;
			}
			else
			{
				fprintf(stderr, "Error: -A/--alignment-mode must be one of: references, minimap2, anchor-imputed, anchor-real.\n");
				exit(1);
			}
			break;
		case 'B':
			if (strlen(optarg) >= sizeof(opt->bowtie2_args))
			{
				fprintf(stderr, "Error: -B/--bowtie2-args is too long.\n");
				exit(1);
			}
			strcpy(opt->bowtie2_args, optarg);
			break;
		case '?':
			// unknown option (getopt has printed it)
			fprintf(stderr, "Run with -h for usage.\n");
			exit(1);
		}
	}

	int anchor_mode = opt->alignment_mode == ALIGNMENT_MODE_ANCHOR_IMPUTED || opt->alignment_mode == ALIGNMENT_MODE_ANCHOR_REAL;
	if (opt->msa_filepath[0] == '\0' || opt->bowtie2_indexes_dir[0] == '\0')
	{
		fprintf(stderr, "Error: -i and -b are required. Run with -h for usage.\n");
		exit(1);
	}
	if (opt->alignment_mode != ALIGNMENT_MODE_ANCHOR_IMPUTED &&
		(opt->non_imputed_positions_filepath[0] == '\0' || opt->reference_sequences_dir[0] == '\0'))
	{
		fprintf(stderr, "Error: -p and -g are required in this alignment mode. Run with -h for usage.\n");
		exit(1);
	}
	if (opt->alignment_mode == ALIGNMENT_MODE_REFERENCES && opt->num_references == 0)
	{
		fprintf(stderr, "Error: no reference strains given. Give them with -r, or choose another -A/--alignment-mode. Run with -h for usage.\n");
		exit(1);
	}
	if (anchor_mode && opt->num_references > 0)
	{
		// not an error, so that scripts that always pass -r can run every mode
		printf("Note: the reference strains given with -r are not used in this alignment mode.\n\n");
	}
	if (optind < argc)
	{
		fprintf(stderr, "Error: unexpected argument '%s'. Run with -h for usage.\n", argv[optind]);
		exit(1);
	}
	if (opt->paired && (opt->forward_end_filepath[0] == '\0' || opt->reverse_end_filepath[0] == '\0'))
	{
		fprintf(stderr, "Error: -P requires both -1 and -2.\n");
		exit(1);
	}
	if (!opt->paired && opt->single_end_filepath[0] == '\0')
	{
		fprintf(stderr, "Error: single-end mode requires -0 (or use -P with -1 and -2 for paired reads).\n");
		exit(1);
	}
}