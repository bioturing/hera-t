#ifndef _OPT_H_
#define _OPT_H_

#include "library_type.h"

struct opt_index_t {
	char *genome;
	char *gtf;
	char *prefix;
	char *idx_dir;
	int k;
	int bwt;
};

struct opt_count_t {
	int n_files;
	char **left_file;
	char **right_file;
	int n_threads;
	char *index;
	char *out_dir;
	char *prefix;
	char *temp_dir;
	int is_dump_align;
	int count_intron;

	// Library type
	struct library_t lib;
};


void print_info();

void print_usage();

struct opt_index_t *get_opt_index(int argc, char *argv[]);

struct opt_count_t *get_opt_count(int argc, char *argv[]);

#endif
