#include <assert.h>
#include <malloc.h>

#include "bwt.h"
#include "interval_tree.h"
#include "genome.h"
#include "io_utils.h"
#include "radix_sort.h"
#include "utils.h"

struct gn_anchor_t {
	int offset;
	bioint_t pos;
	int len;
};

#define ga_get_block(p, s, mask) ((p).pos >> (s) & (mask))
#define ga_less_than(x, y) ((x).pos < (y).pos)

RS_IMPL(anchor2, struct gn_anchor_t, 32, 8, ga_less_than, ga_get_block)

struct gn_seed_t {
	int offset;
	int len;
	bioint_t l, r;
};

static int k_spl = 13;
static int k_gn = 19;
static int k_ext = 39;
static int max_occ = 5000;
static int max_check_len = 100000;
static int intron = 0;

static int err_sub[5][5] = {
				{0, 1, 1, 1, 1},
				{1, 0, 1, 1, 1},
				{1, 1, 0, 1, 1},
				{1, 1, 1, 0, 1},
				{1, 1, 1, 1, 1}
			};

static struct bwt_t bwt;

void genome_init_bwt(struct bwt_t *b, int32_t count_intron)
{
	extern struct bwt_t bwt;
	memcpy(&bwt, b, sizeof(struct bwt_t));
	intron = count_intron;
}

int genome_linear(const uint8_t *pac, const char *seq, int len,
			struct gn_anchor_t *s, int n, int max_err)
{
	int pbeg, pend, err, sbeg, i, j, k, cs, cr;
	bioint_t rbeg;
	rbeg = s->pos;

	pbeg = len;
	pend = 0;
	for (i = 0; i < n; ++i) {
		pbeg = __min(pbeg, s[i].offset);
		pend = __max(pend, s[i].offset + s[i].len);
	}
	assert(pbeg < pend);

	err = 0;

	for (i = 0, k = pbeg; i < n; ++i) {
		sbeg = s[i].offset;
		if (k < sbeg) {
			for (j = k; j < sbeg; ++j) {
				cs = nt4_table[(int)seq[j]];
				cr = __get_pac(pac, rbeg + j);
				err += err_sub[cs][cr];
				if (err >= max_err) return max_err;
			}
			k = sbeg;
		}
		k = __max(sbeg + s[i].len, k);
	}

	if (pbeg) {
		for (i = 0; i < pbeg; ++i) {
			cs = nt4_table[(int)seq[i]];
			cr = __get_pac(pac, rbeg + i);
			err += err_sub[cs][cr];
			if (err >= max_err) return max_err;
		}
	}

	if (len - pend) {
		for (i = pend; i < len; ++i) {
			cs = nt4_table[(int)seq[i]];
			cr = __get_pac(pac, rbeg + i);
			err += err_sub[cs][cr];
			if (err >= max_err) return max_err;
		}
	}

	return err;
}

struct gn_anchor_t *get_anchor(struct gn_seed_t *se, int m, int *n)
{
#define __anchor_lt(x, y) ((x).pos < (y).pos || ((x).pos == (y).pos && (x).offset < (y).offset))
	extern struct bwt_t bwt;
	struct gn_anchor_t **tmp, *ret;
	bioint_t pos, x, j;
	bioint_t *len, *iter;
	int i, k, n_occ, b_i;
	tmp = alloca(m * sizeof(struct gn_anchor_t *));
	ret = NULL;
	len = alloca(m * sizeof(bioint_t));
	*n = 0;
	for (i = 0; i < m; ++i) {
		k = m - i - 1;
		n_occ = se[k].r - se[k].l + 1;
		if (n_occ < max_occ) {
			tmp[i] = malloc(n_occ * sizeof(struct gn_anchor_t));
			x = 0;
			for (j = se[k].l; j <= se[k].r; ++j) {
				pos = bwt_sa(&bwt, j);
				if ((int)pos < se[k].offset) continue;
				tmp[i][x].pos = pos - se[k].offset;
				tmp[i][x].offset = se[k].offset;
				tmp[i][x].len = se[k].len;
				++x;
			}
			if (x)	rs_sort(anchor2, tmp[i], tmp[i] + x);
			*n += x;
			len[i] = x;
		} else {
			len[i] = 0;
			tmp[i] = NULL;
		}
	}

	ret = malloc(*n * sizeof(struct gn_seed_t));
	k = 0;
	iter = alloca(m * sizeof(bioint_t));
	memset(iter, 0, sizeof(bioint_t) * m);
	while (k < *n) {
		b_i = -1;
		for (i = 0; i < m; ++i) {
			if (iter[i] < len[i] && (b_i == -1 ||
				__anchor_lt(tmp[b_i][iter[b_i]],
						tmp[i][iter[i]])))
				b_i = i;
		}
		assert(b_i != -1);
		ret[k++] = tmp[b_i][iter[b_i]++];
	}
	for (i = 0; i < m; ++i)
		free(tmp[i]);
	return ret;
#undef __anchor_lt
}

/*int count_intron(struct gn_anchor_t *s, int n, const char *seq, int len,
		int max_err, struct worker_bundle_t *bundle, char r_str)
{
	int i, k, min_len, err, s_len, l, max;
#ifdef _MSC_VER
	int *start = calloc
#else
	int start[n];
#endif

	min_len = 40;
	max = max_err;
	for (i = l = 0; i < n;) {
		s_len = 0;
		for (k = i; k + 1 < n && s[k].pos == s[k + 1].pos; ++k)
			s_len += s[k].len;
		s_len += s[k].len;
		if (s_len >= min_len) {
			err = genome_linear(bwt.pac, seq, len, s + i, k - i + 1, max);
			if (err < max) {
				max = err;
				l = 1;
				start[0] = s[i].pos;
			} else if (err < max_err && err == max){
				start[l++] = s[i].pos;
			}
		}
		i = k + 1;
	}

	if (!l)
		return -1;

	for (i = 0; i < l && bundle->intron_array->n < 2; ++i)
		query_interval(start[i], start[i] + len - 1, 
				0, bundle->intron_array, r_str);

	return bundle->intron_array->n? 2: 3;
}*/

int check_genome(struct gn_anchor_t *s, int n, const char *seq,
					 int len, int max_err)
{
	int i, k, min_len, err, s_len;

	min_len = 40;
	for (i = 0; i < n;) {
		s_len = 0;
		for (k = i; k + 1 < n && s[k].pos == s[k + 1].pos; ++k)
			s_len += s[k].len;
		s_len += s[k].len;
		if (s_len >= min_len) {
			err = genome_linear(bwt.pac, seq, len, s + i, k - i + 1, max_err);
			if (err < max_err)
				return 3;
		}
		i = k + 1;
	}

	return -1;
}

int get_align_genome(const char *seq, int len, int max_err,
		     struct worker_bundle_t *bundle, char r_str)
{
	extern struct bwt_t bwt;
	bioint_t l, r, o_l, o_r;
	uint8_t c;
	int i, k, m, n, s_len, ret;
	struct gn_anchor_t *s;
	struct gn_seed_t *se;

	ret = -1;
	m = len / k_spl + 1;
	se = malloc(m * sizeof(struct gn_seed_t));
	m = 0;
	s = NULL;
	
	// special case for the 'longest' kmer
	l = 0; r = bwt.seq_len;
	for (i = len - 1; i >= 0; --i) {
		c = nt4_table[(int)seq[i]];
		if (c > 3) break;
		bwt_2occ(&bwt, l - 1, r, c, &o_l, &o_r);
		l = bwt.CC[c] + o_l + 1;
		r = bwt.CC[c] + o_r;
		if (l > r) break;
		s_len = len - i;
		if (s_len >= k_gn) {
			m = 1;
			se[0].offset = i;
			se[0].l = l;
			se[0].r = r;
			se[0].len = s_len;
		}
	}

	if (i < 0 && l <= r) {
		ret = 3;
		goto done;
	}

	for (i = len - k_spl; i - k_spl >= 0; i -= k_spl) {
		l = 0; r = bwt.seq_len;
		for (k = i - 1; k >= 0; --k) {
			c = nt4_table[(int)seq[k]];
			if (c > 3) break;
			bwt_2occ(&bwt, l - 1, r, c, &o_l, &o_r);
			o_l = bwt.CC[c] + o_l + 1;
			o_r = bwt.CC[c] + o_r;
			if (o_l > o_r) break;
			l = o_l;
			r = o_r;
		}
		assert(l <= r);
		s_len = i - k - 1;
		++k;
		if (s_len >= k_gn) {
			se[m].offset = k;
			se[m].len = s_len;
			se[m].l = l;
			se[m].r = r;
			++m;
		}
	}

	s = get_anchor(se, m, &n);
	if (!n)
		goto done;

	/*if (intron)
		ret = count_intron(s, n, seq, len, max_err, bundle, r_str);
	else*/
		ret = check_genome(s, n, seq, len, max_err);
done:
	free(s);
	free(se);
	return ret;
}

int genome_map_err(struct read_t *read, int max_err,
		   struct worker_bundle_t *bundle)
{
	if (!intron && max_err < 0)
		return 0;

	if (max_err == 0)
		return 1;

	int ret;
	ret = max_err < 0? 0: 1;
	max_err = __abs(max_err);
	ret = __max(ret,
		get_align_genome(read->seq, read->len, max_err, bundle, 0));
	if (ret <= 1){
		char *tmp = get_rev_complement(read->seq, read->len);
		ret = __max(ret,
			get_align_genome(tmp, read->len, max_err, bundle, 1));
		free(tmp);
	}

	return ret;
}

// int check_linear_genome(const char *seq, int len, int max_err)
// {
	
// }

// /*
//  * Check whether read has better mapping position on genome
//  * return 1 if found any
//  * return 0 if not
//  */
// int check_on_genome(struct read_t *read, int max_err)
// {
// 	int ret;
// 	if (check_linear_genome(read->seq, read->len, max_err))
// 		return 1;

// 	char *tmp = get_rev_complement(read->seq, read->len);
// 	ret = check_linear_genome(tmp, read->len, max_err);
// 	free(tmp);
// 	if (ret)
// 		return 1;
// 	return 0;
// }
