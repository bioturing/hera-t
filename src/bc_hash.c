#include "bc_hash.h"
#include "utils.h"
#include "mini_hash.h"
#include "atomic.h"

struct bc_hash_t *init_bc_hash()
{
	struct bc_hash_t *bc_hash = malloc(sizeof(struct bc_hash_t));
	init_mini_hash(&bc_hash->h, 10);
	bc_hash->n_bc = 0;
	bc_hash->umi = malloc(1);

	return bc_hash;
}

int32_t get_bc(struct bc_hash_t *bc_hash, int64_t bc)
{
	int32_t n;
	uint64_t *k;
	
	struct mini_hash_t *h = bc_hash->h;
	k = mini_get(h, bc);

	if (k != (uint64_t *)EMPTY_SLOT)
		return *k;

	k = mini_put(&bc_hash->h, bc);
	n = bc_hash->n_bc;
	*k = n;
	bc_hash->umi = realloc (bc_hash->umi, (n + 1) * sizeof(struct umi_hash_t));
	init_mini_hash(&bc_hash->umi[n].h, 7);
	bc_hash->umi[n].type = 0;
	bc_hash->umi[n].count = 0;
	bc_hash->umi[n].idx = bc;
	++bc_hash->n_bc;

	return n;
}

void add_umi(struct umi_hash_t *umi, int64_t umi_ref, int32_t incr, int type)
{
	uint64_t *slot;

	struct mini_hash_t *h = umi->h;
	slot = mini_get(h, umi_ref);

	if (slot == (uint64_t *)EMPTY_SLOT) {
		slot = mini_put(&umi->h, umi_ref);
		if (type == RNA_PRIOR)
			__sync_fetch_and_add64(&umi->count, 1);
	}

	__sync_fetch_and_add64(slot, incr);
}

void add_bc_umi(struct bc_hash_t *bc_hash, int64_t bc, int64_t umi_ref, int32_t type)
{
	int32_t idx = get_bc(bc_hash, bc);
	bc_hash->umi[idx].type |= type;

	add_umi(bc_hash->umi + idx, umi_ref, 1, type);
}

void destroy_bc_hash(struct bc_hash_t *bc_hash)
{
	int i;
	for (i = 0; i < bc_hash->n_bc; ++i)
		destroy_mini_hash(bc_hash->umi[i].h);
	destroy_mini_hash(bc_hash->h);
	free(bc_hash);
}