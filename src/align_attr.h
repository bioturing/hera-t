#ifndef _ALIGN_ATTR_H_
#define _ALIGN_ATTR_H_

#include <stdint.h>

struct anchor_t {
	int beg;			// where should read begin?
	int offset;			// where this anchor begin on read?
};

struct seed_t {
	struct anchor_t *ancs;
	int *beg;
	int n;
	int m;

	int *n_hit;
	int *offset;
	int **hits;
	int n_seed;
	int m_seed;
};

struct align_stat_t {
	int64_t nread;
	int64_t exon;
	int64_t unmap;
	int64_t intron;
	int64_t intergenic;
	int s;
};

struct align_t {
	int pos;			// Position on concatenate sequences
	int score;			// Align score
};

struct raw_alg_t {
	struct align_t *cands;
	int n;
	int m;
	int max_score;
};

struct recycle_alg_t {
	int ref_pos;
	int read_pos;
	int len;
	int score;
};

struct recycle_bin_t {
	struct recycle_alg_t *cands;
	int m;
	int n;
};

#endif
