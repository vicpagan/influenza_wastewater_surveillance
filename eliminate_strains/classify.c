#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>

// A prime number for the hash map size, suitable for ~100,000 reads.
#define HASH_MAP_SIZE 200003
// A large buffer to read an entire line, which could be very long.
#define MAX_LINE_LEN (1024 * 512) // 512 KB buffer

// Data structure to store the result for each read.
// It holds the best mismatch value found in each of the three files.
// It's a node in a linked list to handle hash collisions.
typedef struct ClassificationNode
{
    char *readname;
    int min_mismatch_A;
    int min_mismatch_B;
    int min_mismatch_C;
    struct ClassificationNode *next;
} ClassificationNode;

// The hash map structure itself.
typedef struct
{
    ClassificationNode **buckets;
    int size;
} HashMap;

// Hash function (djb2 algorithm) for strings.
unsigned long hash_function(const char *str)
{
    unsigned long hash = 5381;
    int c;
    while ((c = *str++))
    {
        hash = ((hash << 5) + hash) + c; // hash * 33 + c
    }
    return hash;
}

// Creates and initializes a new hash map.
HashMap *create_hash_map(int size)
{
    HashMap *map = (HashMap *)malloc(sizeof(HashMap));
    map->size = size;
    map->buckets = (ClassificationNode **)calloc(size, sizeof(ClassificationNode *));
    if (!map->buckets)
    {
        perror("Failed to allocate memory for HashMap buckets");
        exit(EXIT_FAILURE);
    }
    return map;
}

// Finds a node for a readname or creates a new one if it doesn't exist.
ClassificationNode *get_or_create_node(HashMap *map, const char *readname)
{
    unsigned long index = hash_function(readname) % map->size;
    ClassificationNode *current = map->buckets[index];

    // Search for an existing node
    while (current != NULL)
    {
        if (strcmp(current->readname, readname) == 0)
        {
            return current;
        }
        current = current->next;
    }

    // Node not found, so create it
    ClassificationNode *newNode = (ClassificationNode *)malloc(sizeof(ClassificationNode));
    newNode->readname = strdup(readname);
    // Initialize all mismatches to a very high value
    newNode->min_mismatch_A = INT_MAX;
    newNode->min_mismatch_B = INT_MAX;
    newNode->min_mismatch_C = INT_MAX;

    // Insert the new node at the front of the bucket's list
    newNode->next = map->buckets[index];
    map->buckets[index] = newNode;

    return newNode;
}

// Processes a single subtype file.
// The 'subtype_id' tells the function which mismatch field to update.
void process_file(const char *filename, HashMap *map, char subtype_id)
{
    FILE *fp = fopen(filename, "r");
    if (!fp)
    {
        fprintf(stderr, "Warning: Could not open file %s. Skipping.\n", filename);
        return;
    }

    printf("Processing %s for Subtype '%c' classification...\n", filename, subtype_id);

    // Use a statically allocated large buffer for performance.
    static char line[MAX_LINE_LEN];

    // Skip the header line, as we no longer need the individual virus names.
    if (fgets(line, sizeof(line), fp) == NULL)
    {
        fprintf(stderr, "Warning: File %s is empty or unreadable.\n", filename);
        fclose(fp);
        return;
    }

    // Read data rows
    while (fgets(line, sizeof(line), fp) != NULL)
    {
        line[strcspn(line, "\r\n")] = 0; // Remove newline
        char *temp_line = strdup(line);  // strtok modifies the string, so use a copy

        // Get the readname from the first column
        char *readname = strtok(temp_line, "\t");
        if (readname == NULL)
        {
            free(temp_line);
            continue; // Skip empty lines
        }

        int min_mismatch_in_row = INT_MAX;
        char *token;

        // Find the minimum mismatch value across all other columns in this row
        while ((token = strtok(NULL, "\t")) != NULL)
        {
            int mismatch = atoi(token);
            if (mismatch < min_mismatch_in_row)
            {
                min_mismatch_in_row = mismatch;
            }
        }

        free(temp_line);

        // If a valid minimum was found for this read in this file
        if (min_mismatch_in_row != INT_MAX)
        {
            // Get the node for this read (or create it)
            ClassificationNode *node = get_or_create_node(map, readname);

            // Update the correct subtype field with the minimum value found
            if (subtype_id == 'A')
            {
                node->min_mismatch_A = min_mismatch_in_row;
            }
            else if (subtype_id == 'B')
            {
                node->min_mismatch_B = min_mismatch_in_row;
            }
            else if (subtype_id == 'C')
            {
                node->min_mismatch_C = min_mismatch_in_row;
            }
        }
    }

    fclose(fp);
    printf("Finished processing %s.\n", filename);
}

// Performs final classification and writes results to a file.
void classify_and_write_results(const char *filename, HashMap *map)
{
    FILE *fp = fopen(filename, "w");
    if (!fp)
    {
        perror("FATAL: Could not open output file for writing");
        exit(EXIT_FAILURE);
    }

    printf("Classifying all reads and writing to %s...\n", filename);

    fprintf(fp, "readname\tsubtype\n");

    // Iterate through every bucket in the hash map
    for (int i = 0; i < map->size; i++)
    {
        ClassificationNode *current = map->buckets[i];
        // Iterate through every node in the bucket's linked list
        while (current != NULL)
        {
            int min_overall = INT_MAX;
            const char *classified_subtype = "Undetermined";

            // Compare the best mismatch from each file to find the winner
            if (current->min_mismatch_A < min_overall)
            {
                min_overall = current->min_mismatch_A;
                classified_subtype = "A/South_Africa/K056869/2023_EPI_ISL_18147004_3C.2a1b.2a.2a.3a.1_MP_A/South_Africa/K056869/2023_EPI2718577";
            }
            if (current->min_mismatch_B < min_overall)
            {
                min_overall = current->min_mismatch_B;
                classified_subtype = "A/black-headed_gull/Leningrad_region/RII-WD319M/2023_EPI_ISL_18115435_2.3.4.4b_MP_A/black-headed_gull/Leningrad_region/RII-WD319M/2023_EPI2691318";
            }
            if (current->min_mismatch_C < min_overall)
            {
                min_overall = current->min_mismatch_C;
                classified_subtype = "A/New_York/PV153865/2024_EPI_ISL_19407925_pdm09_6B.1A.5a.2a.1_MP_68396_PV153865_MP_EPI3554131";
            }

            fprintf(fp, "%s\t%s\n", current->readname, classified_subtype);
            current = current->next;
        }
    }
    fclose(fp);
}

// Frees all memory used by the hash map.
void free_hash_map(HashMap *map)
{
    for (int i = 0; i < map->size; i++)
    {
        ClassificationNode *current = map->buckets[i];
        while (current != NULL)
        {
            ClassificationNode *temp = current;
            current = current->next;
            free(temp->readname);
            free(temp);
        }
    }
    free(map->buckets);
    free(map);
}

int main()
{
    // Create one central hash map to store results for every read
    HashMap *classification_map = create_hash_map(HASH_MAP_SIZE);

    // Process each file, telling the function which subtype it represents
    process_file("mismatch_South_Africa.txt", classification_map, 'A');
    process_file("mismatch_black-headed_gull.txt", classification_map, 'B');
    process_file("mismatch_New_York.txt", classification_map, 'C');

    // Perform the final comparison and write the output
    classify_and_write_results("classification_result.txt", classification_map);

    // Clean up all allocated memory
    free_hash_map(classification_map);

    printf("Classification complete.\n");

    return 0;
}
