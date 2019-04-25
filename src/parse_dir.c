#include <dirent.h>
#include "parse_dir.h"
#include "verbose.h"

static int compare_str(const void *p, const void *q)
{
	const char **l = (const char **) p;
	const char **r = (const char **) q;
	return strcmp (*l, *r);
}

int get_type(char *string, int len)
{
	int i;
	for (i = 0; i + 2 < len; ++i) {
		if (string[i] != '_')
			continue;

		if (string[i + 1] != 'R')
			continue;
		if (string[i + 2] == '1')
			return 0;
		
		if (string[i + 2] == '2')
			return 1;
		++i;
	}

	return -1;
}

void parse_dir(struct input_t *input)
{
	DIR *dirp;
	struct dirent *dp;
	int l, name_len, type, idx;
	char **input_file[2];
	int n[2] = {0, 0};

	name_len = strlen(input->name);
	dirp = opendir(input->in_dir);
	input_file[0] = calloc(1, sizeof(char*));
	input_file[1] = calloc(1, sizeof(char*));

	while ((dp = readdir(dirp)) != NULL) {
		l = strlen(dp->d_name);
		if (l < name_len + 3 || strncmp(dp->d_name, input->name, name_len))
			continue;
		
		type = get_type(dp->d_name + name_len, l - name_len);

		if (type == -1)
			continue;

		idx = n[type];
		++n[type];
		input_file[type] = realloc(input_file[type], n[type] * sizeof (char*));
		input_file[type][idx] = malloc(l + 1);
		memcpy(input_file[type][idx], dp->d_name, l + 1);
	}

	closedir(dirp);

	if (n[0] != n[1])
		__ERROR("The numbers of input files for read1 and read2 are not equal\n");

	if (!n[0])
		__ERROR("No files matched for input sample\n");
	qsort(input_file[0], n[0], sizeof (char *), compare_str);
	qsort(input_file[1], n[1], sizeof (char *), compare_str);

	input->left_file = input_file[0];
	input->right_file = input_file[1];
	input->n_files = n[0];
}
