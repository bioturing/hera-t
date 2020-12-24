#ifndef _BC_HASH_H
#define _BC_HASH_H

#include "khash.h"

KHASH_MAP_INIT_INT64(bc_umi, int);

struct umi_hash_t {
	khash_t(bc_umi) *h;
	int type;
	int count;
	int64_t idx;
};

struct bc_hash_t {
	struct mini_hash_t *h;
	int n_bc;
	struct umi_hash_t *umi;
};

struct bc_hash_t *init_bc_hash();

void add_bc_umi(struct bc_hash_t *h, int64_t bc, int64_t umi_ref, int type);

void add_umi(struct umi_hash_t *umi, int64_t umi_ref, int32_t incr, int type);

void destroy_bc_hash(struct bc_hash_t *h);

#endif