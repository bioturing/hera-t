#include "bc_hash.h"

struct bc_hash_t *init_bc_hash()
{
	struct bc_hash_t *bc_hash = malloc(sizeof(struct bc_hash_t));
	bc_hash->h = kh_init(bc_umi);
	bc_hash->n_bc = 0;
	bc_hash->umi = malloc(1);

	return bc_hash;
}

int32_t get_bc(struct bc_hash_t *bc_hash, int64_t bc)
{
	khiter_t k;
	int32_t ret, n;
	
	khash_t(bc_umi) *h = bc_hash->h;
	k = kh_get(bc_umi, h, bc);

	if (k != kh_end(h))
		return kh_value(h, k);

	k = kh_put(bc_umi, h, bc, &ret);
	n = bc_hash->n_bc;
	kh_value(h, k) = n;
	bc_hash->umi = realloc (bc_hash->umi, (n + 1) * sizeof(struct umi_hash_t));
	bc_hash->umi[n].h = kh_init(bc_umi);
	bc_hash->umi[n].type = 0;
	bc_hash->umi[n].count = 0;
	bc_hash->umi[n].idx = bc;
	++bc_hash->n_bc;

	return n;
}

void add_umi(struct umi_hash_t *umi, int64_t umi_ref, int32_t incr)
{
	khiter_t k;
	int32_t ret;
	
	khash_t(bc_umi) *h = umi->h;
	k = kh_get(bc_umi, h, umi_ref);

	if (k == kh_end(h)) {
		k = kh_put(bc_umi, h, umi_ref, &ret);
		kh_value(h, k) = 0;
		++umi->count;
	}
	
	kh_value(h, k) += incr;
}

void add_bc_umi(struct bc_hash_t *bc_hash, int64_t bc, int64_t umi_ref, int32_t type)
{
	int32_t idx = get_bc(bc_hash, bc);
	bc_hash->umi[idx].type |= type;

	add_umi(bc_hash->umi + idx, umi_ref, 1);
}

void destroy_bc_hash(struct bc_hash_t *bc_hash)
{
	int i;
	for (i = 0; i < bc_hash->n_bc; ++i)
		kh_destroy(bc_umi, bc_hash->umi[i].h);
	kh_destroy(bc_umi, bc_hash->h);
	free(bc_hash);
}