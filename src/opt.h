#ifndef _OPT_H_
#define _OPT_H_

#include "library_type.h"

#ifdef _WIN32
struct opt_index_t {
	TCHAR *genome;
	TCHAR *gtf;
	TCHAR *prefix;
	TCHAR *idx_dir;
	int k;
	int bwt;
};
#else
struct opt_index_t {
	char *genome;
	char *gtf;
	char *prefix;
	char *idx_dir;
	int k;
	int bwt;
};
#endif

#ifdef _WIN32
struct opt_count_t {
	int n_files;
	TCHAR **left_file;
	TCHAR **right_file;
	int n_threads;
	TCHAR *index;
	TCHAR *out_dir;
	TCHAR *prefix;
	TCHAR *temp_dir;
	int is_dump_align;
	int count_intron;

	// Library type
	struct library_t lib;
};
#else
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
#endif


void print_usage();

#ifdef _WIN32
struct opt_index_t *get_opt_index(int argc, TCHAR *argv[]);

struct opt_count_t *get_opt_count(int argc, TCHAR *argv[]);
#else
struct opt_index_t *get_opt_index(int argc, char *argv[]);

struct opt_count_t *get_opt_count(int argc, char *argv[]);
#endif

#endif
