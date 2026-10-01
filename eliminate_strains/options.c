#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>

#include "options.h"

static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"MSA-filepath", required_argument, 0, 'i'},
	{"non-imputed-positions-filepath", required_argument, 0, 'p'},
	{"reference-sequences-dir", required_argument, 0, 'g'},
	{"bowtie2-indexes-dir", required_argument, 0, 'b'},
	{"sam-prefix-filepath", required_argument, 0, 's'},
	{"output-directory", required_argument, 0, 'o'},
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
	{"num-references", required_argument, 0, 'N'},
	{"working-directory", required_argument, 0, 'w'},
	{"verbose", no_argument, 0, 'v'},
	{0, 0, 0, 0}
};

// TODO: add the following line once problematic sites aspect is implemented
// -S, --problematic_sites_dir [REQUIRED,DIR]	Directory of lists of problematic sites\n
char usage[] = "\neliminate_strains [OPTIONS]\n\
	\n\
	-h, --help				\n\
	-i, --MSA-filepath [REQUIRED,FILE]		Filepath of MSA FASTA of influenza reference strains\n\
	-p, --non-imputed-positions-filepath [REQUIRED,FILE]	Filepath of the non-imputed positions for each strain within the MSA\n\
	-g, --reference-sequences-dir [REQUIRED,DIR]	Directory of reference sequences\n\
	-b, --bowtie2-indexes-dir [REQUIRED,DIR]	Directory of precomputed bowtie2 indexes for each reference sequence\n\
	-r, --reference-strains [REQUIRED, NAME NAME ...]	List of the names of reference strains to use (e.g., -r EPI_ISL_19088566 EPI_ISL_19407907 ...)\n\
	-s, --sam-prefix-filepath [REQUIRED, FILE]		Output sam file to print alignments\n\
	-o, --output-dir [DIR]		Directory to place output files [default: .]\n\
	-N, --num-references [REQUIRED,int]	Number of reference strains to use for alignment\n\
	-P, --paired				Using paired-reads\n\
	-m, --max-num-reads [int]	Maximum number of reads from the read file(s) to process [default: 100,000]\n\
	-0, --single-end-filepath [FILE]		Single-end reads\n\
	-1, --forward-end-filepath [FILE]		If using paired-reads, the forward reads file\n\
	-2, --reverse-end-filepath [FILE]		If using paired-reads, the reverse reads file\n\
	-e, --EM-error [decimal]		Error rate for EM algorithm [default: 0.005]\n\
	-a, --fasta				Reads are in FASTA format [default: FASTQ]\n\
	-l, --llr				Perform the LLR procedure\n\
	-z, --num-top-strains-llr [int]	Number of highest proportion strains to calculate the per-strain LLR for [default: 10]\n\
	-t, --cores [decimal]			Number of cores [default: 1]\n\
	-n, --no-read-sam			Don't read in sam file to memory\n\
	-w, --working-dir [DIR]			Directory for intermediate/working files [default: .]\n\
	-v, --verbose				Show verbose debug output from Bowtie2 commands\n\
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
	opt->num_references = 0;

	int i;
	for (i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-r") == 0 || strcmp(argv[i], "--reference-strains") == 0)
		{
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
	char c;
	if (argc == 1)
	{
		print_help_statement();
		exit(0);
	}
	while (1)
	{
		c = getopt_long(argc, argv, "hPdlna:i:p:s:f:o:0:1:2:e:t:c:m:x:b:g:r:j:N:k:w:y:W:q:u:z:vR", long_options, &option_index);
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
		case 's':
			success = sscanf(optarg, "%s", opt->sam_prefix_filepath);
			if (!success)
				fprintf(stderr, "Invalid SAM filepath\n");
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
		case 'N':
			success = sscanf(optarg, "%d", &(opt->num_references));
			if (!success)
				fprintf(stderr, "Invalid number of references\n");
			break;
		case 'w':
			success = sscanf(optarg, "%s", opt->working_dir);
			if (!success)
				fprintf(stderr, "Invalid working directory\n");
			break;
		case 'v':
			opt->verbose = 1;
			break;
		case 'z':
			success = sscanf(optarg, "%d", &(opt->num_top_strains_llr));
			if (!success)
				fprintf(stderr, "Invalid number of top strains to compute LLR for\n");
			break;
		}
	}
}