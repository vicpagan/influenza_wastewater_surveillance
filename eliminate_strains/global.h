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
	char msa_filepath[2048];
	char non_imputed_positions_filepath[2048];
	char reference_sequences_dir[2048];
	char bowtie2_indexes_dir[2048];
	char problematic_sites_dir[2048];

	// SAM file to write/read alignments
	char sam_prefix_filepath[2048];
	
	// read inputs
	int paired;
	int fasta_format;
	char single_end_filepath[2048];
	char forward_end_filepath[2048];
	char reverse_end_filepath[2048];

	// output files
	char output_dir[2048];
	char working_dir[2048];

	// algorithm parameters
	double em_error;
	int llr;
	int num_top_strains_llr;
	int num_references;
	char **reference_strain_names;

	// debug parameters
	int verbose;
	
	// performance parameters
	int num_threads;
	int max_num_reads;
	int no_read_bam;
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

	int **references_to_msa_positions;
	int *reference_sequence_lengths;
	SAMResults sam_results_str;
	// ProblematicSites problematic_sites_str;
} ReferencesData;

typedef struct MismatchData
{
	int num_reads;
	int num_msa_sequences;

	char **read_names;
	char **msa_sequence_names;
	int *alignment_sizes;

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
	int column_index;
	char *msa_strain_name;
	double proportion;
} ProportionData;


#endif /* _GLOBAL_ */
