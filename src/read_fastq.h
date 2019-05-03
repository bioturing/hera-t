#ifndef _READ_FASTQ_H
#define _READ_FASTQ_H

#include "opt.h"
#include "dqueue.h"

struct pair_buffer_t {
	char *buf1;
	char *buf2;
	int input_format;
};

struct dqueue_t *init_read_fastq(int n_threads, struct input_t *input);

struct pair_buffer_t *init_pair_buffer();

void free_pair_buffer(struct pair_buffer_t *p);

void finish_read_fastq(struct dqueue_t *q);

#endif