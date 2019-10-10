#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>

#include "alignment.h"
#include "dynamic_alignment.h"
#include "genome.h"
#include "hash_table.h"
#include "io_utils.h"
#include "barcode.h"
#include "radix_sort.h"
#include "verbose.h"

#define recycle_get_block(p, s, mask) ((p).ref_pos >> (s) & (mask))
#define recycle_less_than(x, y) ((x).ref_pos < (y).ref_pos)

RS_IMPL(recycle, struct recycle_alg_t, 32, 8, recycle_less_than, recycle_get_block)

static struct gene_info_t genes;
static struct transcript_info_t trans;

static int kcons;
static uint64_t kcons_mask;

#define ERROR_RATIO		0.055
#define ACCEPT_RATIO		0.65
#define PARTIAL_RATIO		0.4
#define DROP_RATIO		0.8
#define MAX_DIST		1000

void alignment_init_ref_info(struct gene_info_t *g, struct transcript_info_t *t)
{
	extern struct gene_info_t genes;
	extern struct transcript_info_t trans;
	memcpy(&genes, g, sizeof(struct gene_info_t));
	memcpy(&trans, t, sizeof(struct transcript_info_t));
}

void alignment_init_hash(const char *path)
{
	extern int kcons;
	extern uint64_t kcons_mask;
	load_cons_hash(path, &kcons);
	kcons_mask = (1ull << (kcons << 1)) - 1;
}

static inline uint64_t get_index_cons(const char *seq)
{
	uint64_t ret = 0;
	int i;
	uint8_t c;
	for (i = 0; i < kcons; ++i) {
		c = nt4_table[(int)seq[i]];
		// if (c >= 4)
		// 	c = lrand48() & 3;
		if (c > 3)
			return (uint64_t)(-1);
		ret = (ret << 2) | c;
	}
	return ret;
}

void get_cons_seed(struct read_t *read, struct seed_t *rs, int step)
{
	extern struct transcript_info_t trans;
	uint64_t idx;
	int i, k, s, n, m;
	int *ret;
	n = rs->n_seed;
	if (n + 1 > rs->m_seed) {
		rs->m_seed = n + 1;
		rs->n_hit = realloc(rs->n_hit, rs->m_seed * sizeof(int));
		rs->offset = realloc(rs->offset, rs->m_seed * sizeof(int));
		rs->hits = realloc(rs->hits, rs->m_seed * sizeof(int *));
	}

	for (i = 0; i < n; ++i) {
		// s = i * seed_info->step;
		s = i * step;
		rs->offset[i] = s;

		idx = get_index_cons(read->seq + s);
		if (idx == (uint64_t)(-1))
			m = 0;
		else
			m = query_cons_hash(idx, &ret);

		rs->n_hit[i] = m;

		if (rs->n + m > rs->m) {
			rs->m = rs->n + m;
			__round_up_32(rs->m);
			rs->beg = realloc(rs->beg, rs->m * sizeof(int));
			rs->ancs = realloc(rs->ancs,
					rs->m * sizeof(struct anchor_t));
		}

		for (k = 0; k < m; ++k)
			rs->beg[rs->n + k] = ret[k] - s;
		rs->n += m;
	}

	rs->hits[0] = rs->beg;
	for (i = 0; i < n; ++i)
		rs->hits[i + 1] = rs->hits[i] + rs->n_hit[i];
}

void find_cons_seeds(struct read_t *read, struct seed_t *ret)
{
	extern int kcons;
	int step;
	step = (kcons + 1) / 2;
	ret->n_seed = (read->len - kcons) / step + 1;
	if (ret->n_seed == 1) {
		ret->n_seed = 2;
		step = read->len - kcons;
	}
	// ret->scons = malloc(ret->ncons * sizeof(struct hit_t));
	get_cons_seed(read, ret, step);
}

void merge_seed(struct seed_t *s)
{
	/* make sure s->ancs have size at least s->n */
	int n, ns, k, best_i, i;
	int *n_hit, *c;
	int **hits;
	struct anchor_t *a;

	n_hit = s->n_hit;
	hits = s->hits;
	a = s->ancs;
	n = s->n;
	ns = s->n_seed;
	c = alloca(ns * sizeof(int));
	memset(c, 0, ns * sizeof(int));
	k = 0;
	while (k < n) {
		best_i = -1;
		for (i = 0; i < ns; ++i) {
			if (c[i] < n_hit[i] && (best_i == -1 ||
			     hits[i][c[i]] < hits[best_i][c[best_i]]))
				best_i = i;
		}
		assert(best_i != -1);
		a[k].beg = hits[best_i][c[best_i]];
		a[k].offset = s->offset[best_i];
		++k;
		++c[best_i];
	}
}

void linear_cons_anchor(char *seq, int len, struct anchor_t *a, int n,
			struct raw_alg_t *ret, struct recycle_bin_t *bin,
			int max_err, int partial_score, int clip)
{
	extern struct transcript_info_t trans;
	int ref_id, ref_len, ref_pos;
	int expected_score, aligned_base, clipped_base, err, drop_thres, error_quota;
	int score, pbeg, pend, seq_beg, seq_end;
	int chain_beg, chain_end, chain_score;
	int i, k, sbeg;
	char *ref;
	struct extend_align_t ext_alg;

	ref_id = trans.idx[a->beg + a->offset];
	ref_len = trans.tran_len[ref_id];
	ref_pos = a->beg - trans.tran_beg[ref_id];

	/* literally, there is some kind of clip 30S70M at postion 10...
	 * Ignore ref_pos < 0
	 */

	/* only check substring [seq_beg, seq_end) */
	if (ref_pos + len > ref_len)
		seq_end = ref_len - ref_pos;
	else
		seq_end = len;
	if (ref_pos < 0)
		seq_beg = -ref_pos;
	else
		seq_beg = 0;
	/*  ref have offset [seq_beg] compare to seq */
	ref = trans.seq + trans.tran_beg[ref_id] + (ref_pos + seq_beg);

	error_quota = max_err * SUB_GAP;
	drop_thres = error_quota * DROP_RATIO;
	score = 0;

	pbeg = seq_end;
	pend = seq_beg;
	for (i = 0; i < n; ++i) {
		pbeg = __min(pbeg, a[i].offset);
		pend = __max(pend, a[i].offset + kcons);
	}
	assert(pbeg < pend);

	for (i = 0, k = pbeg; i < n; ++i) {
		sbeg = a[i].offset;
		if (k < sbeg) { /* linear check segment that is lack of info */
			score += b2b_check_nocigar(seq + k, ref + (k - seq_beg),
						   sbeg - k, &error_quota);
			if (error_quota < 0) return;
			k = sbeg;
		}
		assert(k <= sbeg + kcons);
		score += (sbeg + kcons - k) * SUB_MAX;
		k = sbeg + kcons;
	}

	chain_beg = pbeg;
	chain_end = pend;
	chain_score = score;

	// linear to begin
	if (pbeg - seq_beg) {
		k = pbeg - seq_beg;
		error_quota -= align_linear_bw(ref + k, seq + pbeg, k,
					drop_thres, error_quota, &ext_alg);
		score += ext_alg.score;
		pbeg -= ext_alg.seq_len;
	}

	if (seq_end - pend) {
		k = pend - seq_beg;
		error_quota -= align_linear_fw(ref + k, seq + pend,
				seq_end - pend, drop_thres, error_quota, &ext_alg);
		score += ext_alg.score;
		pend += ext_alg.seq_len;
	}

	aligned_base = pend - pbeg;
	clipped_base = len - aligned_base;

	if (aligned_base == len) {
		expected_score = len * SUB_MAX - max_err * SUB_GAP;
	} else {
		err = max_err / 2 - (len - aligned_base) / 5;
		err = __max(0, err);
		expected_score = aligned_base * SUB_MAX - err * SUB_GAP;
	}

	if (score >= expected_score && clipped_base <= clip &&
	    (pbeg == 0 || pend == len)) {
		if (ret->m == ret->n) {
			ret->m <<= 1;
			ret->cands = realloc(ret->cands,
					ret->m * sizeof(struct align_t));
		}

		ret->cands[ret->n].pos = trans.tran_beg[ref_id] + ref_pos + pbeg;
		ret->cands[ret->n].score = score;
		ret->max_score = __max(ret->max_score, score);
		++ret->n;
	} else if (score >= partial_score || chain_score >= partial_score) {
		if (score < expected_score || score < partial_score) {
			pbeg = chain_beg;
			pend = chain_end;
			score = chain_score;
		}
		if (bin->m == bin->n) {
			bin->m <<= 1;
			bin->cands = realloc(bin->cands,
					bin->m * sizeof(struct recycle_alg_t));
		}
		// bin->cands[bin->n].ref_pos = ref_pos + trans.tran_beg[ref_id];
		bin->cands[bin->n].ref_pos = a->beg;
		bin->cands[bin->n].read_pos = pbeg;
		bin->cands[bin->n].len = pend - pbeg;
		bin->cands[bin->n].score = score;
		++bin->n;
	}
}

void rescue_perfect(struct read_t *read, struct raw_alg_t *ret,
			struct recycle_bin_t *bin)
{
	extern struct transcript_info_t trans;
	struct recycle_alg_t *cands;
	cands = bin->cands;
	int i, k, n, len, score, pbeg, pend, seq_beg, seq_end;
	int ref_pos, ref_len, ref_id, ref_beg, ref_end;
	int error_quota, err, max_err, drop_thres, clip;
	int max_score, expected_score, aligned_base, clipped_base;
	char *ref, *seq;
	struct extend_align_t ext_alg;
	len = read->len;
	n = bin->n;

	max_err = len * ERROR_RATIO;
	max_score = len * SUB_MAX;
	error_quota = max_err * SUB_GAP;
	drop_thres = error_quota * DROP_RATIO;
	clip = max_err;
	seq = read->seq;

	int ibin = n;
	for (i = 0; i < n; ++i) {
		pbeg = cands[i].read_pos;
		pend = pbeg + cands[i].len;
		score = cands[i].score;
		ref_pos = cands[i].ref_pos;
		ref_id = trans.idx[ref_pos + pbeg];
		ref_len = trans.tran_len[ref_id];
		ref_pos -= trans.tran_beg[ref_id];

		ref_beg = ref_pos + pbeg;
		ref_end = ref_beg + cands[i].len;
		if (ref_pos + read->len > ref_len)
			seq_end = ref_len - ref_pos;
		else
			seq_end = read->len;
		if (ref_pos < 0)
			seq_beg = -ref_pos;
		else
			seq_beg = 0;
		ref = trans.seq + trans.tran_beg[ref_id] + (ref_pos + seq_beg);

		if (pbeg - seq_beg) {
			k = pbeg - seq_beg;
			// linear_bw_nocigar(ref + k, read->seq + pbeg, k, error_quota, &ext_alg);
			error_quota -= align_linear_bw(ref + k, seq + pbeg, k,
						drop_thres, error_quota, &ext_alg);
			score += ext_alg.score;
			pbeg -= ext_alg.seq_len;
			ref_beg -= ext_alg.seq_len;
		}
		if (seq_end - pend) {
			k = pend - seq_beg;
			// linear_fw_nocigar(ref + k, read->seq + pend, seq_end - pend, error_quota, &ext_alg);
			error_quota -= align_linear_fw(ref + k, seq + pend,
				seq_end - pend, drop_thres, error_quota, &ext_alg);
			score += ext_alg.score;
			pend += ext_alg.seq_len;
			ref_end += ext_alg.seq_len;
		}
		aligned_base = pend - pbeg;
		clipped_base = len - aligned_base;

		if (aligned_base == len) {
			expected_score = max_score - max_err * SUB_GAP;
		} else {
			err = max_err / 2 - clipped_base / 5;
			err = __max(0, max_err);
			expected_score = aligned_base * SUB_MAX - err * SUB_GAP;
		}

		if (score >= expected_score
		    && clipped_base <= clip && (pbeg == 0 || pend == len)) {
			if (ret->m == ret->n) {
				ret->m <<= 1;
				ret->cands = realloc(ret->cands, ret->m * sizeof(struct align_t));
			}
			// ret->cands[ret->n].pos = ref_beg;
			// ret->cands[ret->n].tid = ref_id;
			ret->cands[ret->n].pos = ref_beg + trans.tran_beg[ref_id];
			ret->cands[ret->n].score = score;
			ret->max_score = __max(ret->max_score, score);
			++ret->n;
			cands[i] = cands[--ibin];
		}
	}
	bin->n = ibin;
}

void semi_indel_map(char *seq, int len, struct recycle_alg_t *seed,
			struct raw_alg_t *ret, struct array_2D_t *tmp_array)
{
	extern struct transcript_info_t trans;
	int ref_id, ref_len, ref_pos, ref_sent;
	int max_score, accept_score, error_quota;
	int score, pbeg, pend, ref_beg, ref_end, ext_len, ext_err;
	int aligned_base, clipped_base, expected_score, max_err, err, drop_thres;
	char *ref;
	struct extend_align_t ext_alg;

	pbeg = seed->read_pos;
	pend = pbeg + seed->len;

	ref_pos = seed->ref_pos;
	ref_id = trans.idx[ref_pos + pbeg];
	ref_len = trans.tran_len[ref_id];
	ref_pos -= trans.tran_beg[ref_id];
	ref_beg = ref_pos + pbeg;
	ref_end = ref_beg + seed->len;

	ref_sent = ref_beg;

	/* ref have offset [pbeg] compare to read->seq
	 */
	ref = trans.seq + trans.tran_beg[ref_id] + ref_sent;

	max_score = len * SUB_MAX;
	max_err = ERROR_RATIO * len;
	accept_score = max_score * ACCEPT_RATIO;
	error_quota = max_err * SUB_GAP;
	drop_thres = error_quota * DROP_RATIO;

	ext_len = max_err * SUB_GAP_RATIO;

	score = seed->score;
	if (pbeg) {
		ext_err = align_banded_bw(ref, seq + pbeg, __min(ref_sent,
								pbeg + ext_len),
					  pbeg, ext_len, drop_thres,
					  error_quota, tmp_array, &ext_alg);
		error_quota -= ext_err;
		score += ext_alg.score;
		ref_beg -= ext_alg.ref_len;
		pbeg -= ext_alg.seq_len;
	}

	if (len - pend) {
		ext_err = align_banded_fw(ref + seed->len, seq + pend,
					  __min(len - pend + ext_len,
						ref_len - ref_sent - seed->len),
					  len - pend, ext_len, drop_thres,
					  error_quota, tmp_array, &ext_alg);
		error_quota -= ext_err;
		score += ext_alg.score;
		ref_end += ext_alg.ref_len;
		pend += ext_alg.seq_len;
	}
	aligned_base = pend - pbeg;
	clipped_base = len - aligned_base;

	if (aligned_base == len) {
		expected_score = max_score - max_err * SUB_GAP;
	} else {
		err = max_err / 2 - clipped_base / 5;
		err = __max(0, err);
		expected_score = aligned_base * SUB_MAX - err * SUB_GAP;
	}

	if (score >= accept_score && score >= expected_score            /* Score is good */
	    && (pbeg == 0 || pend == len)) {                      /* Clip at one end */
		if (ret->m == ret->n) {
			ret->m <<= 1;
			ret->cands = realloc(ret->cands, ret->m * sizeof(struct align_t));
		}

		ret->cands[ret->n].pos = ref_beg + trans.tran_beg[ref_id];
		ret->cands[ret->n].score = score;
		ret->max_score = __max(ret->max_score, score);
		++ret->n;
	}
}

void get_perfect_map(struct read_t *read, struct seed_t *s,
		     struct worker_bundle_t *bundle)
{
	int n, i, k, cs, n_seed;
	int max_score, partial_score;
	struct anchor_t *a;
	struct raw_alg_t *ret;
	struct recycle_bin_t *bin;
	bin = bundle->recycle_bin;
	ret = bundle->alg_array;
	a = s->ancs;
	n = s->n;
	n_seed = s->n_seed;

	max_score = read->len * SUB_MAX;
	partial_score = max_score * PARTIAL_RATIO;
	for (i = 0; i < n;) {
		for (k = i; k + 1 < n && a[k].beg == a[k + 1].beg &&
			trans.idx[a[k].beg + a[k].offset] == trans.idx[a[k + 1].beg + a[k + 1].offset]; ++k);
		cs = k - i + 1;
		if (cs == n_seed)
			linear_cons_anchor(read->seq, read->len, a + i, cs,
					ret, bin, 0, partial_score, 0);
		i = k + 1;
	}
}

void get_perfect_map_alt(struct read_t *read, struct seed_t *s,
			 struct worker_bundle_t *bundle)
{
	int n, i, k, cs, n_seed;
	int max_score, max_err, partial_score, clip;
	struct anchor_t *a;
	struct raw_alg_t *ret;
	struct recycle_bin_t *bin;
	bin = bundle->recycle_bin;
	ret = bundle->alg_array;
	a = s->ancs;
	n = s->n;
	n_seed = s->n_seed;

	max_score = read->len * SUB_MAX;
	max_err = read->len * ERROR_RATIO;
	partial_score = max_score * PARTIAL_RATIO;
	clip = max_err;
	for (i = 0; i < n;) {
		for (k = i; k + 1 < n && a[k].beg == a[k + 1].beg &&
			trans.idx[a[k].beg + a[k].offset] == trans.idx[a[k + 1].beg + a[k + 1].offset]; ++k);
		cs = k - i + 1;
		if (cs >= n_seed - 2)
			linear_cons_anchor(read->seq, read->len, a + i, cs,
					ret, bin, max_err, partial_score, clip);
		i = k + 1;
	}
}

void get_linear_map(struct read_t *read, struct seed_t *s, int min_s, int max_s,
			struct worker_bundle_t *bundle)
{
	int n, i, k, cs;
	int max_score, max_err, partial_score, clip;
	struct anchor_t *a;
	struct raw_alg_t *ret;
	struct recycle_bin_t *bin;
	bin = bundle->recycle_bin;
	ret = bundle->alg_array;
	a = s->ancs;
	n = s->n;

	max_score = read->len * SUB_MAX;
	max_err = read->len * ERROR_RATIO;
	partial_score = max_score * PARTIAL_RATIO;
	clip = max_err;
	for (i = 0; i < n;) {
		for (k = i; k + 1 < n && a[k].beg == a[k + 1].beg &&
			trans.idx[a[k].beg + a[k].offset] == trans.idx[a[k + 1].beg + a[k + 1].offset]; ++k);
		cs = k - i + 1;
		if (cs >= min_s && cs <= max_s)
			linear_cons_anchor(read->seq, read->len, a + i, cs,
					ret, bin, max_err, partial_score, clip);
		i = k + 1;
	}
}

int check_linear_map(struct read_t *read, struct worker_bundle_t *bundle)
{
	int err, max_err, max_score;
	struct raw_alg_t *algs;
	struct seed_t *s_cons;

	s_cons = bundle->seed_cons;
	algs = bundle->alg_array;
	max_score = read->len * SUB_MAX;
	max_err = read->len * ERROR_RATIO;

	get_perfect_map(read, s_cons, bundle);
	if (algs->n)
		goto genome_check;

	rescue_perfect(read, bundle->alg_array, bundle->recycle_bin);
	get_perfect_map_alt(read, s_cons, bundle);
	if (algs->n && algs->max_score >= max_score - 2 * SUB_GAP)
		goto genome_check;

	get_linear_map(read, s_cons, 1, s_cons->n_seed - 3, bundle);
	if (algs->n)
		goto genome_check;

	return 0;

genome_check:
	err = (max_score - algs->max_score + SUB_GAP - 1) / SUB_GAP;
	return genome_map_err(read, __min(max_err, err), bundle);
}

int check_indel_map(struct read_t *read, struct worker_bundle_t *bundle)
{
	int max_err = read->len * ERROR_RATIO;

	if (!bundle->recycle_bin->n)
		return genome_map_err(read, -max_err, bundle);

	struct recycle_bin_t *bin;
	struct array_2D_t *tmp_array;
	struct raw_alg_t *ret;
	struct recycle_alg_t *cands;
	int n, i, max_score, err;

	bin = bundle->recycle_bin;
	tmp_array = bundle->tmp_array;
	ret = bundle->alg_array;
	cands = bin->cands;
	n = bin->n;

	for (i = 0; i < n; ++i)
		semi_indel_map(read->seq, read->len, cands + i, ret, tmp_array);

	if (ret->n) {
		max_score = read->len * SUB_MAX;
		err = (max_score - ret->max_score + SUB_GAP - 1) / SUB_GAP;
		err = __min(max_err, err);
	} else {
		err = -max_err;
	}

	return genome_map_err(read, err, bundle);
}

void store_read_chromium(struct read_t *r, struct raw_alg_t *alg,
			struct kmhash_t *bc_table, pthread_mutex_t *lock_hash,
			struct library_t lib)
{
	int i, g, gene;
	uint64_t bc_idx, umi_gene_idx;

	gene = -1;
	for (i = 0; i < alg->n; ++i) {
		if (alg->cands[i].score < alg->max_score)
			continue;
		assert(alg->cands[i].score == alg->max_score);
		g = trans.gene_idx[trans.idx[alg->cands[i].pos]];
		if (gene == -1)
			gene = g;
		else if (gene != g) {
			gene = -1;
			break;
		}
	}
	if (gene == -1)
		return;
	bc_idx = umi_gene_idx = 0;
	for (i = 0; i < lib.bc_len; ++i)
		bc_idx = bc_idx * 5 + nt4_table[(int)r->seq[i]];
	for (i = 0; i < lib.umi_len; ++i)
		umi_gene_idx = umi_gene_idx * 5 + nt4_table[(int)r->seq[lib.bc_len + i]];
	umi_gene_idx = umi_gene_idx << GENE_BIT_LEN | gene;
	kmhash_put_bc_umi(bc_table, lock_hash, bc_idx, umi_gene_idx);
}

/*
void store_exon(struct read_t *read1, struct raw_alg_t *alg)
{
	int32_t i, gene, t, stat;
	int64_t bc_id;

	idx = seq2num(read1->seq, read1->len);
	gene = -1;
	for (i = 0; i < alg->n; ++i) {
		t = alg->cands[i].pos;
		if (gene == -1)
			gene = trans.gene_idx[trans.idx[t]];
		else if (gene != trans.gene_idx[trans.idx[t]])
			break;
	}
	
	if (i != alg->n)
		return;
	
	// add_bc_umi(idx, gene);
}

void store_intron(struct read_t *read1, struct interval_t *intron)
{
	if (intron->n != 1)
		return;
	
	// add_bc_umi(seq2num(read1->seq, read1->len), intron->id[0]);
}
*/

void write_aligns(struct shared_fstream_t *fstream, struct read_t *r1,
		  struct read_t *r2, struct raw_alg_t *a)
{
	if (!fstream) return;
	extern struct transcript_info_t trans;
	extern struct gene_info_t genes;
	int l, alg_len, rem_buf, i, tid, gid;
	l = fstream->buf_len;
	rem_buf = SFS_BUF_SZ - l;
		 /* r1->name       + "\t" */
	alg_len = strlen(r1->name) + 1    + a->n * (genes.l_id + trans.l_id + 10 + 3);
	if (alg_len > rem_buf) { // not enough expected buffer
		sfs_flush(fstream);
		rem_buf = SFS_BUF_SZ;
		l = 0;
	}
	if (rem_buf < alg_len)
		__ERROR("Wrtting alignments: insufficient amount of buffer, please report to us.");
	l += sprintf(fstream->buf + l, "%s", r1->name);
	for (i = 0; i < a->n; ++i) {
		tid = trans.idx[a->cands[i].pos];
		gid = trans.gene_idx[tid];
		l += sprintf(fstream->buf + l, "\t%s\t%s\t%d\t%d",
				trans.tran_id + tid * trans.l_id,
				genes.gene_id + gid * genes.l_id,
				a->cands[i].pos - trans.tran_beg[tid],
				a->cands[i].score);
	}
	fstream->buf[l++] = '\n';
	fstream->buf_len = l;
}

void align_chromium_read(struct read_t *read1, struct read_t *read2,
			 struct worker_bundle_t *bundle)
{
	if (!read1->name || !read1->seq || !read2->name || !read2->seq)
		return;

	int r1_len = bundle->lib.bc_len + bundle->lib.umi_len;
	if (read1->len < r1_len)
		__ERROR("Read lenght of %s is not consistent with library type.\n Expect >= %u.\n Receive %u.\n", read1->name, r1_len, read1->len);

	++bundle->result->nread;
	reinit_bundle(bundle);

	int ret;
	struct seed_t *s_cons;

	s_cons = bundle->seed_cons;
	find_cons_seeds(read2, s_cons);
	merge_seed(s_cons);

	ret = check_linear_map(read2, bundle);
	if (ret == 0)
		ret = check_indel_map(read2, bundle);

	if (ret == 1){
		store_read_chromium(read1, bundle->alg_array, bundle->bc_table,
					bundle->lock_hash, bundle->lib);
		// store_exon(read1, bundle->alg_array);
		++bundle->result->exon;
	} else if (ret == 2) {
		// store_intron(read1, bundle->intron_array);
		++bundle->result->intron;
	} else if (ret == 3) 
		++bundle->result->intergenic;
	else
		++bundle->result->unmap;
}
