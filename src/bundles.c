#include "bundles.h"

/*****************************************
*            WORKING WITH BUNDLE         *
*****************************************/

struct raw_alg_t *init_raw_alg()
{
	struct raw_alg_t *ret;
	ret = malloc(sizeof(struct raw_alg_t));
	ret->n = 0;
	ret->m = 64;
	ret->max_score = 0;
	ret->cands = malloc(ret->m * sizeof(struct align_t));
	return ret;
}

void reinit_align_array(struct raw_alg_t *p)
{
	assert(p != NULL);
	p->n = 0;
	p->max_score = 0;
}

void destroy_raw_alg(struct raw_alg_t *p)
{
	if (!p) return;
	free(p->cands);
	free(p);
}

struct interval_t *init_intron_array()
{
	struct interval_t *ret;
	ret = malloc(sizeof(struct interval_t));
	ret->n = 0;
	ret->m = 64;
	ret->id = malloc(ret->m * sizeof(int));
	return ret;
}

void reinit_intron_array(struct interval_t *p)
{
	p->n = 0;
}

void destroy_intron_array(struct interval_t *p)
{
	if (!p) return;

	free(p->id);
	free(p);
}

struct array_2D_t *init_array_2D(int nrow, int ncol, int word)
{
	struct array_2D_t *ret = malloc(sizeof(struct array_2D_t));
	ret->len = nrow * ncol * word;
	ret->data = malloc(ret->len);
	ret->nrow = nrow;
	ret->rows = malloc(ret->nrow * sizeof(void *));
	int i;
	for (i = 0; i < nrow; ++i)
#ifndef _MSC_VER
		ret->rows[i] = (void *)(ret->data + i * ncol * word);
#else
		ret->rows[i] = (int *)((char*)(ret->data) + i * ncol * word);
#endif
	return ret;
}

void **resize_array_2D(struct array_2D_t *p, int nrow, int ncol, int word)
{
	if (!p)
		return NULL;
	int len = nrow * ncol * word;
	if (len > p->len) {
		p->len = len;
		p->data = realloc(p->data, p->len);
	}
	if (nrow > p->nrow) {
		p->nrow = nrow;
		p->rows = realloc(p->rows, p->nrow * sizeof(void *));
	}

	int i;
	for (i = 0; i < nrow; ++i)
		p->rows[i] = (void *)((char *)p->data + i * ncol * word);
	return p->rows;
}

void destroy_array_2D(struct array_2D_t *p)
{
	if (!p) return;
	free(p->data);
	free(p->rows);
	free(p);
}

struct recycle_bin_t *init_recycle_bin()
{
	struct recycle_bin_t *ret = malloc(sizeof(struct recycle_bin_t));
	ret->n = 0;
	ret->m = 64;
	ret->cands = malloc(ret->m * sizeof(struct recycle_alg_t));
	return ret;
}

void reinit_recycle_bin(struct recycle_bin_t *p)
{
	assert(p != NULL);
	p->n = 0;
}

void destroy_recycle_bin(struct recycle_bin_t *p)
{
	if (!p) return;
	free(p->cands);
	free(p);
}

struct seed_t *init_seed()
{
	struct seed_t *ret;
	ret = calloc(1, sizeof(struct seed_t));
	ret->m = 0x100;
	ret->beg = malloc(ret->m * sizeof(int));
	ret->ancs = malloc(ret->m * sizeof(struct anchor_t));
	ret->m_seed = 0x8;
	ret->n_hit = malloc(ret->m_seed * sizeof(int));
	ret->offset = malloc(ret->m_seed * sizeof(int));
	ret->hits = malloc(ret->m_seed * sizeof(int *));
	return ret;
}

void reinit_seed(struct seed_t *p)
{
	assert(p != NULL);
	p->n = p->n_seed = 0;
}

void destroy_seed(struct seed_t *p)
{
	if (!p) return;
	free(p->beg);
	free(p->ancs);
	free(p->n_hit);
	free(p->offset);
	free(p->hits);
	free(p);
}

void init_bundle(struct bundle_data_t *bundle)
{
	bundle->alg_array = init_raw_alg();
	bundle->intron_array = init_intron_array();
	bundle->tmp_array = init_array_2D(100, 100, 4);
	bundle->recycle_bin = init_recycle_bin();
	bundle->seed_cons = init_seed();
}

void reinit_bundle(struct bundle_data_t *bundle)
{
	reinit_align_array(bundle->alg_array);
	reinit_intron_array(bundle->intron_array);
	reinit_recycle_bin(bundle->recycle_bin);
	reinit_seed(bundle->seed_cons);
}

void destroy_bundle(struct bundle_data_t *bundle)
{
	destroy_raw_alg(bundle->alg_array);
	destroy_intron_array(bundle->intron_array);
	destroy_array_2D(bundle->tmp_array);
	destroy_recycle_bin(bundle->recycle_bin);
	destroy_seed(bundle->seed_cons);
}