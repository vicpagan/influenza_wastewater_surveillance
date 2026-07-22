#ifndef _GLOBAL_
#define _GLOBAL_

#define FASTA_MAXLINE 30000 // max length of a single line sequence
#define MAX_CIGAR 1000 // max number of CIGAR operations in an alignment
#define MAX_READ_LENGTH 1000 // max length of a single read

/**
 * @brief Struct to hold command line options
 * 
 */
typedef struct Options
{
	// MSA, reference, and alignment files
	char msa_filepath[1000];
	char msa_reference_dir[1000];
	char bowtie2_reference_dir[1000];
	char variant_sites_dir[1000];
	char problematic_sites_dir[1000];

	// SAM file to write/read alignments
	char sam_filepath[1000];
	
	// read inputs
	int paired;
	int fasta_format;
	int clean_reads;
	char single_end_filepath[1000];
	char forward_end_filepath[1000];
	char reverse_end_filepath[1000];

	// output files
	char output_filepath[1000];
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
	int num_references;
	
	// performance parameters
	int number_of_cores;
	int no_read_bam;
	int remove_identical;
} Options;

/**
 * @brief 
 * 
 */
typedef struct MSA
{
	int num_sequences;
	int sequence_length;
	int max_sequence_name_length;
	// int reference_sequence_index;

	char **sequences;
	char **sequence_names;

} MSA;

// TODO: Future improvement - instead of storing the entire SAM lines and parsing them when calculating the mismatch, just store the important parts
// typedef struct SAMRecord
// {
// 	char *qname;
// 	int flag;
// 	int pos;
// 	char *cigar;
// 	char *seq;
// 	int edit_distance;
// } SAMRecord;

typedef struct SAMResults
{
	int num_sam_lines;
	int max_sam_line_length;

	char **sam_results;
} SAMResults;


typedef struct VariantSites
{
	int *variant_sites;
	int num_variant_sites;
} VariantSites;

// TODO: Implement problematic sites considerations
// typedef struct ProblematicSites
// {
// 	int *problematic_sites;
// 	int num_problematic_sites;
// } ProblematicSites;

/**
 * @brief 
 * 
 */
typedef struct ReferenceData
{
	int *reference_index;
	VariantSites variant_sites_str;
	// ProblematicSites problematic_sites_str;
	SAMResults sam_results_str;
} ReferenceData;

/**
 * @brief Struct to hold thread parameters for parallel processing the mismatch matrix
 * 
 */
typedef struct MismatchMatrixThreadStruct
{
	int sam_line_start;
	int sam_line_end;
	int thread_index;
	char **mismatch_matrix_row_partition;
} MismatchMatrixThreadStruct;

#endif /* _GLOBAL_ */
