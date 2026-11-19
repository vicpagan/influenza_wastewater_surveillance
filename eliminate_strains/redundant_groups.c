#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "redundant_groups.h"
#include "file_utils.h"

typedef struct GroupMember
{
	char *member_name;
	char *group_name;
} GroupMember;

static int compare_group_members(const void *a, const void *b)
{
	return strcmp(((const GroupMember *)a)->member_name, ((const GroupMember *)b)->member_name);
}

RedundantGroups read_redundant_groups(char *redundant_groups_filepath)
{
	RedundantGroups redundant_groups;
	redundant_groups.num_members = 0;
	redundant_groups.num_groups = 0;
	redundant_groups.member_names = NULL;
	redundant_groups.group_names = NULL;

	if (redundant_groups_filepath[0] == '\0')
	{
		return redundant_groups;
	}

	gzFile file;
	if ((file = gzopen(redundant_groups_filepath, "r")) == (gzFile)NULL)
	{
		fprintf(stderr, "Error! Cannot open redundant groups file '%s'.\n", redundant_groups_filepath);
		exit(1);
	}

	int capacity = 1 << 12;
	GroupMember *members = (GroupMember *)malloc(capacity * sizeof(GroupMember));
	int num_members = 0;
	char *current_group_name = NULL;  // shared by the members of the group
	char *line = NULL;
	size_t line_capacity = 0;

	while (gz_read_line(file, &line, &line_capacity))
	{
		if (line[0] == '>')
		{
			char *name = strtok(line + 1, " \t");
			current_group_name = name != NULL ? strdup(name) : NULL;
			redundant_groups.num_groups++;
			continue;
		}

		char *name = strtok(line, " \t");
		while (name != NULL)
		{
			if (current_group_name == NULL)
			{
				fprintf(stderr, "Error: redundant groups file '%s' has a strain ('%s') before the first '>' group line.\n", redundant_groups_filepath, name);
				exit(1);
			}
			if (num_members == capacity)
			{
				capacity *= 2;
				members = (GroupMember *)realloc(members, capacity * sizeof(GroupMember));
			}
			members[num_members].member_name = strdup(name);
			members[num_members].group_name = current_group_name;
			num_members++;
			name = strtok(NULL, " \t");
		}
	}
	free(line);
	gzclose(file);

	qsort(members, num_members, sizeof(GroupMember), compare_group_members);

	redundant_groups.num_members = num_members;
	redundant_groups.member_names = (char **)malloc((num_members > 0 ? num_members : 1) * sizeof(char *));
	redundant_groups.group_names = (char **)malloc((num_members > 0 ? num_members : 1) * sizeof(char *));
	int i;
	for (i = 0; i < num_members; i++)
	{
		redundant_groups.member_names[i] = members[i].member_name;
		redundant_groups.group_names[i] = members[i].group_name;
	}
	free(members);

	return redundant_groups;
}

const char *find_group_of_strain(RedundantGroups *redundant_groups, const char *strain_name)
{
	int low = 0;
	int high = redundant_groups->num_members - 1;
	while (low <= high)
	{
		int middle = low + (high - low) / 2;
		int comparison = strcmp(strain_name, redundant_groups->member_names[middle]);
		if (comparison == 0)
		{
			return redundant_groups->group_names[middle];
		}
		if (comparison < 0)
		{
			high = middle - 1;
		}
		else
		{
			low = middle + 1;
		}
	}
	return NULL;
}
