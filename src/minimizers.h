//
// Created by BioTuring on 2019-11-10.
//

#ifndef SKIPPING_MINIMIZERS_H
#define SKIPPING_MINIMIZERS_H

#include <stdint.h>
#include <stdint.h>
#include <limits.h>
#include <stdio.h>
#include <stddef.h>
#include "khash.h"

#define RATIO_OF_CONFIDENT 0.85
#define MIN_NUMBER_SINGLETON 2
#define EMPTY_BX UINT64_MAX

KHASH_MAP_INIT_INT64(mm_align, struct mm_align_t *);

struct mm_db_t {
    uint64_t *mm;
    uint32_t *p;
    size_t n;
    size_t size;
    int k;
};

struct mm_db_t * mm_index_char_str(char *s, int k, int w, int l);
struct mm_db_edge_t *mm_index_edges(struct asm_graph_t *g, int k, int w);
void mm_db_destroy(struct mm_db_t *db);

#endif //SKIPPING_MINIMIZERS_H
