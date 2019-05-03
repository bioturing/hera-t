#ifndef _PARSE_DIR_H_
#define _PARSE_DIR_H_

struct input_t {
	char *in_dir;
	char *name;
	int n_files;
	char **left_file;
	char **right_file;
	char *ref;
};

void parse_dir(struct input_t *input);

#endif