#ifndef _GLOBAL_
#define _GLOBAL_

#define FASTA_MAXLINE 30000 // max length of a single line sequence
#define MAX_CIGAR 1000 // max number of CIGAR operations in an alignment
#define MAX_READ_LENGTH 1000 // max length of a single read

typedef struct Options
{
	// reference and alignment files
	char fasta[1000];
	char MSA_reference[1000];
	char bowtie2_reference[1000];
	char variant[1000];
	
	// read inputs
	int paired;
	int fasta_format;
	int clean_reads;
	char single_end_file[1000];
	char forward_end_file[1000];
	char reverse_end_file[1000];

	// SAM file to write/read alignments
	char sam[1000];

	// output files
	char outfile[1000];
	char print_counts[1000];
	char print_deletions[1000];

	// algorithm parameters
	double freq;
	double error;
	int coverage;
	double deletion_threshold;
	int min_strains;
	int max_strains;
	int llr;
	
	// performance parameters
	int number_of_cores;
	int no_read_bam;
	int remove_identical;
} Options;

typedef struct ResultsStruct
{
	char **mismatch;
} ResultsStruct;

typedef struct ThreadStruct
{
	int start;
	int end;
	int thread_id;
	int max_sam_length;
	int length_of_MSA;
	int number_of_strains;
	int number_of_strains_remaining;
	ResultsStruct *results_str;
} ThreadStruct;

extern char **resize_MSA;
extern char **resize_names_of_strains;
extern int *reference_index;
extern char **sam_results;

#endif /* _GLOBAL_ */
