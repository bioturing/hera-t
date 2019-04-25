#ifndef _READ_TSV_H_
#define _READ_TSV_H_

#include "io_utils.h"

struct content_t {
	char *s;
	int l;
};

struct tsv_t {
	FILE *f;
	char *buf;

	int *row_content;
	int l;
	int m;
	int sep;
};

struct tsv_t *init_readTSV(const char *path, int sep);

void destroy_readTSV(struct tsv_t *f);

ssize_t get_row_content(struct tsv_t *f);

struct content_t get_col_content(struct tsv_t *f, int i);

#endif