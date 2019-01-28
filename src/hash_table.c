#include <string.h>
#include <stdlib.h>
#include <stdint.h>

#include "hash_table.h"
#include "io_utils.h"
#include "verbose.h"

static struct cons_build_t *bcons;

static struct cons_table_t *hcons;

/* CONS HASH */

void init_cons_hash(int l2_size)
{
	// assert(size > 0);
	// assert((size & (size - 1)) == 0);
	extern struct cons_build_t *bcons;
	bcons = calloc(1, sizeof(struct cons_build_t));
	bcons->l2_size = l2_size;
	int size = 1 << l2_size;
	bcons->mask = size - 1;
	bcons->buckets = calloc(size, sizeof(struct cons_bucket_t *));
	bcons->bsize = calloc(size + 1, sizeof(int));
}

int insert_cons_hash(uint64_t id)
{
	extern struct cons_build_t *bcons;
	uint32_t p = id & bcons->mask;
	uint32_t kid = id >> bcons->l2_size;
	struct cons_bucket_t *buck = bcons->buckets[p];
	int n = bcons->bsize[p], l = 0, r = n, mid, ret;

	while (l < r) {
		mid = (l + r) >> 1;
		if (buck[mid].id == kid) {
			l = mid;
			break;
		}
		if (buck[mid].id > kid)
			r = mid;
		else
			l = mid + 1;
	}

	if (l == n || buck[l].id > kid) {
		bcons->buckets[p] = realloc(bcons->buckets[p], (n + 1) *
					    sizeof(struct cons_bucket_t));
		if (l < n)
			memmove(bcons->buckets[p] + l + 1, bcons->buckets[p] + l,
				(n - l) * sizeof(struct cons_bucket_t));
		bcons->buckets[p][l].id = kid;
		bcons->buckets[p][l].head = 1;
		++bcons->bsize[p];
		ret = 1;
	} else {
		++bcons->buckets[p][l].head;
		ret = 0;
	}

	return ret;
}

void recount_cons_hash()
{
	extern struct cons_build_t *bcons;
	int i, k, tmp, size = bcons->mask + 1;
	bcons->npos = 0;

	for (i = 0; i < size; ++i) {
		for (k = 0; k < bcons->bsize[i]; ++k) {
			tmp = bcons->npos;
			bcons->npos += bcons->buckets[i][k].head;
			bcons->buckets[i][k].head = tmp;
		}
	}

	bcons->pos = malloc(bcons->npos * sizeof(int));
}

void addpos_cons_hash(uint64_t id, int pos)
{
	extern struct cons_build_t *bcons;
	uint32_t p = (uint32_t)(id & bcons->mask);
	uint32_t kid = id >> bcons->l2_size;
	struct cons_bucket_t *buck = bcons->buckets[p];
	int n = bcons->bsize[p], l = 0, r = n, mid;

	while (l < r) {
		mid = (l + r) >> 1;
		if (buck[mid].id == kid) {
			l = mid;
			break;
		}
		if (buck[mid].id > kid)
			r = mid;
		else
			l = mid + 1;
	}

	bcons->pos[bcons->buckets[p][l].head++] = pos;
}

int query_cons_hash(uint64_t id, int **pos)
{
	extern struct cons_table_t *hcons;
	uint32_t p = id & hcons->mask;
	uint32_t kid = id >> hcons->l2_size;
	uint32_t *buck = hcons->id;
	int *head = hcons->head, l = hcons->bpos[p], r = hcons->bpos[p + 1], mid;

	while (l < r) {
		mid = (l + r) >> 1;
		if (buck[mid] == kid) {
			*pos = hcons->pos + head[mid];
			return head[mid + 1] - head[mid];
		}
		if (buck[mid] > kid)
			r = mid;
		else
			l = mid + 1;
	}

	*pos = NULL;
	return 0;
}

void store_cons_hash(const char *file_path, int kcons)
{
	extern struct cons_build_t *bcons;
	FILE *fi = xfopen(file_path, "wb");

	xfwrite(&kcons, sizeof(int), 1, fi);
	int size = bcons->mask + 1, tmp, sum, max, k, i;
	for (i = sum = 0; i < size; ++i) {
		for (k = 0; k < bcons->bsize[i]; ++k) {
			tmp = bcons->buckets[i][k].head;
			bcons->buckets[i][k].head = sum;
			sum = tmp;
		}
	}

	for (i = sum = max = 0; i <= size; ++i) {
		if (bcons->bsize[i] > max)
			max = bcons->bsize[i];
		tmp = bcons->bsize[i];
		bcons->bsize[i] = sum;
		sum += tmp;
	}

	uint32_t *id = malloc(max * sizeof(uint32_t));
	int *head = malloc(max * sizeof(int));
	int bsize;
	xfwrite(&(bcons->l2_size), sizeof(int), 1, fi);
	xfwrite(bcons->bsize, sizeof(int), size + 1, fi);
	for (i = 0; i < size; ++i) {
		bsize = bcons->bsize[i + 1] - bcons->bsize[i];
		if (bsize > 0) {
			for (k = 0; k < bsize; ++k)
				id[k] = bcons->buckets[i][k].id;
			xfwrite(id, sizeof(uint32_t), bsize, fi);
		}
	}
	for (i = 0; i < size; ++i) {
		bsize = bcons->bsize[i + 1] - bcons->bsize[i];
		if (bsize > 0) {
			for (k = 0; k < bsize; ++k)
				head[k] = bcons->buckets[i][k].head;
			xfwrite(head, sizeof(int), bsize, fi);
		}
	}

	xfwrite(&(bcons->npos), sizeof(int), 1, fi);
	xfwrite(bcons->pos, sizeof(int), bcons->npos, fi);

	xfclose(fi);
	free(id);
	free(head);
}

void load_cons_hash(const char *file_path, int *kcons)
{
	extern struct cons_table_t *hcons;
	FILE *fi = xfopen(file_path, "rb");

	hcons = calloc(1, sizeof(struct cons_table_t));
	xfread(kcons, sizeof(int), 1, fi);
	
	int size;
	xfread(&(hcons->l2_size), sizeof(int), 1, fi);
	size = 1 << hcons->l2_size;
	hcons->mask = size - 1;

	hcons->bpos = malloc((size + 1) * sizeof(int));
	xfread(hcons->bpos, sizeof(int), size + 1, fi);

	hcons->id = malloc(hcons->bpos[size] * sizeof(uint32_t));
	xfread(hcons->id, sizeof(uint32_t), hcons->bpos[size], fi);

	hcons->head = malloc((hcons->bpos[size] + 1) * sizeof(int));
	xfread(hcons->head, sizeof(int), hcons->bpos[size] + 1, fi);

	hcons->pos = malloc(hcons->head[hcons->bpos[size]] * sizeof(int));
	xfread(hcons->pos, sizeof(int), hcons->head[hcons->bpos[size]], fi);

	xfclose(fi);
}

void free_cons_hash_index()
{
	extern struct cons_build_t *bcons;
	if (!bcons)
		return;

	int size = bcons->mask + 1, i;
	for (i = 0; i < size; ++i)
		free(bcons->buckets[i]);
	free(bcons->buckets);
	free(bcons->pos);
	free(bcons->bsize);
	free(bcons);
	bcons = NULL;
}

void free_cons_hash()
{
	extern struct cons_table_t *hcons;
	if (!hcons)
		return;

	free(hcons->id);
	free(hcons->head);
	free(hcons->bpos);
	free(hcons->pos);
	free(hcons);
	hcons = NULL;
}
