#include "utils.h"

#define BUF_SIZE 		1048576

int8_t nt4_table[256] = {
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4
};

char *nt4_char = "ACGTN";
char *rev_nt4_char = "TGCAN";

char *str_concate(const char *str1, const char *str2)
{
	/* TODO: put this to depricated */
	size_t len1 = strlen(str1);
	size_t len2 = strlen(str2);
	char *str3 = malloc(len1 + len2 + 1);
	strcpy(str3, str1);
	strcpy(str3 + len1, str2);
	str3[len1 + len2] = '\0';
	return str3;
}

char *get_rev(const char *seq, int len)
{
	if (seq == NULL)
		return NULL;

	int i, k;
	char *ret = malloc(len + 1);
	for (i = 0, k = len - 1; i < len; ++i, --k)
		ret[i] = seq[k];
	ret[len] = '\0';
	return ret;
}

char *get_rev_complement(const char *seq, int len)
{
	if (seq == NULL)
		return NULL;

	char *ret = malloc(len + 1);
	int i, k;
	for (i = 0, k = len - 1; i < len; ++i, --k) 
		ret[i] = rev_nt4_char[(int)nt4_table[seq[k]]];
	ret[len] = '\0';
	return ret;
}

double realtime()
{
	struct timeval tp;
#if defined(_MSC_VER)
	_gettimeofday(&tp,  NULL);
#else
	gettimeofday(&tp,  NULL);
#endif
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

int check_valid_nu(const char *seq, int len)
{
	int i;
	for (i = 0; i < len; ++i)
		if (nt4_table[(int)seq[i]] >= NNU)
			return 0;
	return 1;
}

int64_t seq2num(const char *seq, int len)
{
	int64_t ret = 0;
	int i;
	for (i = 0; i < len; ++i)
		ret = ret * 5 + nt4_table[(int)seq[i]];
	return ret;
}

char *num2seq(int64_t num, int len)
{
	char *ret = malloc(len + 1);
	int i;
	for (i = 0; i < len; ++i) {
		ret[len - i - 1] = nt4_char[num % 5];
		num /= 5;
	}
	ret[len] = '\0';
	return ret;
}

struct reader_pos_t {
	int pos;
	int cnt;
};

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

void free_pair_buffer(struct pair_buffer_t *p)
{
	if (!p) return;
	free(p->buf1);
	free(p->buf2);
	free(p);
}

struct pair_buffer_t *init_pair_buffer()
{
	struct pair_buffer_t *ret = malloc(sizeof(struct pair_buffer_t));
	ret->buf1 = malloc(BUF_SIZE + 1);
	ret->buf2 = malloc(BUF_SIZE + 1);
	return ret;
}

struct dqueue_t *init_dqueue_PE(int cap)
{
	struct dqueue_t *ret = init_dqueue(cap);
	struct pair_buffer_t *p;
	int i;
	for (i = 0; i < cap; ++i) {
		p = init_pair_buffer();
		d_enqueue_out(ret, p);
	}
	return ret;
}

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

struct interval_t *init_intron_array()
{
	struct interval_t *ret;
	ret = malloc(sizeof(struct interval_t));
	ret->n = 0;
	ret->m = 64;
	ret->id = malloc(ret->m * sizeof(int));
	return ret;
}

void destroy_raw_alg(struct raw_alg_t *p)
{
	if (!p) return;
	free(p->cands);
	free(p);
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
	if (!p) return NULL;
	int len = nrow * ncol * word;
	if (len > p->len) {
		p->len = len;
		p->data = realloc(p->data, p->len);
	}
	if (nrow > p->nrow) {
		p->nrow = nrow;
		p->rows = realloc(p->rows, p->nrow * sizeof(void *));
	}
//	memset(p->data, 0, len);
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

void reinit_align_array(struct raw_alg_t *p)
{
	assert(p != NULL);
	p->n = 0;
	p->max_score = 0;
}

void reinit_intron_array(struct interval_t *p)
{
	p->n = 0;
}

void reinit_recycle_bin(struct recycle_bin_t *p)
{
	assert(p != NULL);
	p->n = 0;
}

void reinit_seed(struct seed_t *p)
{
	assert(p != NULL);
	p->n = p->n_seed = 0;
}

void destroy_recycle_bin(struct recycle_bin_t *p)
{
	if (!p) return;
	free(p->cands);
	free(p);
}

void init_bundle(struct worker_bundle_t *bundle)
{
	bundle->alg_array = init_raw_alg();
	bundle->intron_array = init_intron_array();
	bundle->tmp_array = init_array_2D(100, 100, 4);
	bundle->recycle_bin = init_recycle_bin();
	bundle->seed_cons = init_seed();
}

void reinit_bundle(struct worker_bundle_t *bundle)
{
	reinit_align_array(bundle->alg_array);
	reinit_intron_array(bundle->intron_array);
	reinit_recycle_bin(bundle->recycle_bin);
	reinit_seed(bundle->seed_cons);
}

void destroy_bundle(struct worker_bundle_t *bundle)
{
	destroy_raw_alg(bundle->alg_array);
	destroy_intron_array(bundle->intron_array);
	destroy_array_2D(bundle->tmp_array);
	destroy_recycle_bin(bundle->recycle_bin);
	destroy_seed(bundle->seed_cons);
}
