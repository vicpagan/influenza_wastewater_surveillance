#ifndef _GLOBAL_
#define _GLOBAL_

#include <pthread.h>

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
	char problematic_sites_dir[1000];

	// SAM file to write/read alignments
	char sam_prefix_filepath[1000];
	
	// read inputs
	int paired;
	int fasta_format;
	int clean_reads;
	char single_end_filepath[1000];
	char forward_end_filepath[1000];
	char reverse_end_filepath[1000];

	// output files
	char print_counts_filepath[1000];
	char print_deletions_filepath[1000];
	char output_dir[1000];
	char working_dir[1000];

	// algorithm parameters
	double freq;
	double em_error;
	int coverage;
	double deletion_threshold;
	int min_strains;
	int max_strains;
	int llr;
	int num_references;

	// clean_reads parameters
	int end_region_length;
	double end_region_error_mult;
	int sequence_length_threshold;
	int trim_length;
	int fastq_trimmer_threshold;
	int verbose;
	
	// performance parameters
	int num_threads;
	int no_read_bam;
	int remove_identical_sequences;
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

// TODO: Add SAM results filepath for executions where we DONT want the entire SAM file written in

// TODO: Make a single SAMResults struct hold an array of SAMFile structs that hold the actual lines in the file
// This is because the num sam lines and the max sam line length should be the same value for every SAM file anyway
// Maybe just use the SAMRecord stuff from above?
typedef struct SAMResults
{
	int num_sam_lines;
	int max_sam_line_length;

	char ***sam_results;
} SAMResults;


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
typedef struct ReferencesData
{
	int num_references;
	int *reference_sequence_msa_indexes;
	char **reference_names;

	int **reference_indexes;
	SAMResults sam_results_str;
	// ProblematicSites problematic_sites_str;
} ReferencesData;

typedef struct MismatchData
{
	int num_reads;
	int num_msa_sequences;

	char **read_names;
	char **msa_sequence_names;
	int *block_sizes;

	int **mismatch_matrix;
} MismatchData;

/**
 * @brief Struct to hold thread parameters for parallel processing the mismatch matrix
 * 
 */
typedef struct BuildMismatchMatrixThread
{
	int sam_partition_start;
	int sam_partition_end;
	int read_strain_offset;
	int thread_index;
	
	ReferencesData *references_data_str;
	MSA *msa_str;

	MismatchData *mismatch_data_str;
} BuildMismatchMatrixThread;

typedef struct HashmapEntry
{
	char *mismatch_column;
	char *msa_strain_names;
	int msa_strain_names_length;
} HashmapEntry;

typedef struct ProportionData
{
	char *msa_strain_name;
	double proportion;
} ProportionData;


#endif /* _GLOBAL_ */
