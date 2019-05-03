#ifndef _PROCESS_TAG_H
#define _PROCESS_TAG_H

#include "khash.h"
#include "opt.h"
#include "barcode.h"

void build_tag_ref(struct input_t *input, struct ref_info_t *ref, int type);

void print_tag_count(int thread_num);

void print_tag_stat(int n_threads);

int align_tag(struct read_t *read, int thread_num);

void init_tag_threads(int n_threads);

void destroy_tag_ref();

#endif