#ifndef REDUNDANT_GROUPS_H
#define REDUNDANT_GROUPS_H

#include "global.h"

/**
 * @brief Reads the redundant groups file of a nonredundant MSA: a ">cN" line for each group,
 * followed by the names of its member strains (one or more per line). Members are identical
 * after imputation and share the MSA row named after their group.
 *
 * @param redundant_groups_filepath the file (may be gzip-compressed); an empty string gives no groups
 * @return RedundantGroups members sorted by name, each with the name of its group
 */
RedundantGroups read_redundant_groups(char *redundant_groups_filepath);

/**
 * @brief Name of the group (MSA row) a strain belongs to
 *
 * @return const char* the group name, or NULL if the strain is in no group
 */
const char *find_group_of_strain(RedundantGroups *redundant_groups, const char *strain_name);

#endif // REDUNDANT_GROUPS_H
