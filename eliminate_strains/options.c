#include "options.h"

static struct option long_options[] =
{
	{"help", no_argument, 0, 'h'},
	{"MSA-filepath", required_argument, 0, 'i'},
	{"samfile", required_argument, 0, 's'},
	{"freq", required_argument, 0, 'f'},
	{"outfile", required_argument, 0, 'o'},
	{"variant_dir", required_argument, 0, 'v'},
	{"paired", no_argument, 0, 'p'},
	{"single_end", required_argument, 0, '0'},
	{"forward_read", required_argument, 0, '1'},
	{"reverse_read", required_argument, 0, '2'},
	{"EM-error", required_argument, 0, 'e'},
	{"coverage", required_argument, 0, 'c'},
	{"fasta", no_argument, 0, 'a'},
	{"llr", no_argument, 0, 'l'},
	{"min", required_argument, 0, 'm'},
	{"max", required_argument, 0, 'x'},
	{"print-allele-counts", required_argument, 0, 'b'},
	{"cores", required_argument, 0, 't'},
	{"MSA-reference-dir", required_argument, 0, 'g'},
	{"no-read-sam", no_argument, 0, 'n'},
	{"print-deletions", required_argument, 0, 'r'},
	{"clean-my-reads", no_argument, 0, 'd'},
	{"bowtie2-alignment_dir", required_argument, 0, 'B'},
	{"num-references", required_argument, 0, 'N'},
	{0, 0, 0, 0}
};

// TODO: add the following line once problematic sites aspect is implemented
// -p, --problematic_sites_dir [REQUIRED,DIR]	Directory of lists of problematic sites\n
char usage[] = "\neliminate_strains [OPTIONS]\n\
	\n\
	-h, --help				\n\
	-i, --MSA-filepath [REQUIRED,FILE]		Filepath of MSA FASTA of influenza reference strains\n\
	-s, --sam-filepath [REQUIRED,FILE]		Output sam file to print alignments\n\
	-f, --freq [REQUIRED,decimal]		Allele frequency to filter unlikely strains [default: 0.01]\n\
	-o, --output-filepath [REQUIRED,FILE]		Output file to print mismatch matrix for EM algorithm\n\
	-v, --variant_sites_dir [REQUIRED,DIR]	Directory of lists of variant sites\n\
	-g, --MSA-reference-dir [REQUIRED,DIR]	Directory of MSA reference sequences\n\
	-N, --num-references [REQUIRED,int]	Number of reference strains to use for alignment\n\
	-P, --paired				Using paired-reads\n\
	-0, --single_end_file [FILE]		Single-end reads\n\
	-1, --forward_file [FILE]		If using paired-reads, the forward reads file\n\
	-2, --reverse_file [FILE]		If using paired-reads, the reverse reads file\n\
	-e, --EM-error [decimal]		Error rate for EM algorithm\n\
	-d, --clean-my-reads                    Clean reads with fastq_quality_trimmer [must have FASTQ reads]\n\
	-c, --coverage [integer]		Number of reads needed to calculate allele freq [default: 50]\n\
	-a, --fasta				Reads are in FASTA format [default: FASTQ]\n\
	-l, --llr				Perform the LLR procedure\n\
	-m, --min [decimal]			Minimum strains remaining to invoke iterative procedure [default: 100]\n\
	-x, --max [decimal]			Maximum strains remaining for EM algorithm [default: 10000]\n\
	-b, --print-allele-counts [FILE]	Print allele counts to file\n\
	-t, --cores [decimal]			Number of cores [default: 1]\n\
	-n, --no-read-sam			Don't thread, don't read in sam file to memory\n\
	-r, --print-deletions [FILE]		Print sites with deletions\n\
	-j, --threshold-for-deleted-sites	Threshold to print deleted sites [default: 0.001]\n\
	-B, --bowtie2-alignment_dir [REQUIRED,DIR]		Bowtie2 reference\n\
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

// TODO: add the following lines once problematic sites aspect is implemented (additionally make sure 'p' is in the getopt_long() string)
// case 'v':
// 			success = sscanf(optarg, "%s", opt->variant_sites_dir);
// 			if (!success)
// 				fprintf(stderr, "Invalid variant sites directory\n");
// 			break;

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
		c = getopt_long(argc, argv, "hPdlnaB:i:s:f:o:v:0:1:2:e:t:c:m:x:b:g:r:j:N:", long_options, &option_index);
		if (c == -1)
			break;
		switch (c)
		{
		case 'h':
			print_help_statement();
			exit(0);
			break;
		case 'd':
			opt->clean_reads = 1;
			break;
		case 'b':
			success = sscanf(optarg, "%s", opt->print_counts);
			if (!success)
				fprintf(stderr, "Invalid counts file\n");
			break;
		case 'r':
			success = sscanf(optarg, "%s", opt->print_deletions);
			if (!success)
				fprintf(stderr, "Invalid deletions file\n");
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
			success = sscanf(optarg, "%s", opt->msa_reference_dir);
			if (!success)
				fprintf(stderr, "Invalid MSA reference directory\n");
			break;
		case 'B':
			success = sscanf(optarg, "%s", opt->bowtie2_reference_dir);
			if (!success)
				fprintf(stderr, "Invalid reference directory\n");
			break;
		case 'n':
			opt->no_read_bam = 1;
			break;
		case 'i':
			success = sscanf(optarg, "%s", opt->msa_filepath);
			if (!success)
				fprintf(stderr, "Invalid MSA filepath\n");
			break;
		case 's':
			success = sscanf(optarg, "%s", opt->sam_filepath);
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
		case 'f':
			success = sscanf(optarg, "%lf", &(opt->freq));
			if (!success)
				fprintf(stderr, "Invalid freq\n");
			break;
		case 'j':
			success = sscanf(optarg, "%lf", &(opt->deletion_threshold));
			if (!success)
				fprintf(stderr, "Invalid threshold\n");
			break;
		case 't':
			success = sscanf(optarg, "%d", &(opt->number_of_cores));
			if (!success)
				fprintf(stderr, "Invalid number of cores\n");
			break;
		case 'c':
			success = sscanf(optarg, "%d", &(opt->coverage));
			if (!success)
				fprintf(stderr, "Invalid freq\n");
			break;
		case 'm':
			success = sscanf(optarg, "%d", &(opt->min_strains));
			if (!success)
				fprintf(stderr, "Invalid min strains\n");
			break;
		case 'x':
			success = sscanf(optarg, "%d", &(opt->max_strains));
			if (!success)
				fprintf(stderr, "Invalid min strains\n");
			break;
		case 'e':
			success = sscanf(optarg, "%lf", &(opt->error));
			if (!success)
				fprintf(stderr, "Invalid error rate\n");
			break;
		case 'o':
			success = sscanf(optarg, "%s", opt->output_filepath);
			if (!success)
				fprintf(stderr, "Invalid out file\n");
			break;
		case 'v':
			success = sscanf(optarg, "%s", opt->variant_sites_dir);
			if (!success)
				fprintf(stderr, "Invalid variant sites directory\n");
			break;
		case 'N':
			success = sscanf(optarg, "%d", &(opt->num_references));
			if (!success)
				fprintf(stderr, "Invalid number of references\n");
			break;
		}
	}
}
