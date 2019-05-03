#ifndef _OPT_H_
#define _OPT_H_

#include "parse_dir.h"
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
	struct library_t lib;
	int n_threads;
	char *out_dir;
	int count_intron;

	// gene expression
	struct input_t *rna;
	// cell hashing
	struct input_t *cell;
	// protein quant
	struct input_t *protein;
	// crispr
	struct input_t *crispr;
};

void print_info();

void print_usage();

struct opt_index_t *get_opt_index(int argc, char *argv[]);

struct opt_count_t *get_opt_count(int argc, char *argv[]);

#endif
