#ifndef _PROCESS_RNA_H_
#define _PROCESS_RNA_H_

#include "verbose.h"
#include "barcode.h"

void load_rna_index(const char *prefix, int32_t count_intron);

void add_rna_ref(struct ref_info_t *ref);

void init_rna_threads(int n_threads);

int align_rna(struct read_t *read, int thread_num);

void destroy_rna_index(int n_threads);

void print_rna_count(int thread_num);

void print_rna_stat(int n_threads, int count_intron);

#endif