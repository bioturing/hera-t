#ifndef _DYNAMIC_ALIGNMENT_H_
#define _DYNAMIC_ALIGNMENT_H_

#include "attribute.h"
#include "bundles.h"

#define SUB_MAX		1
#define SUB_MIN		-2
#define SUB_GAP		3
#define GAP_E		-3
#define GAP_O		-1
#define SUB_GAP_RATIO	0.75

struct extend_align_t {
	int score;
	int ref_len;
	int seq_len;
};

int b2b_check_nocigar(const char *ref, const char *seq, int len, int *err_quota);

int align_linear_fw(const char *ref, const char *seq, int len, int drop,
			int error_quota, struct extend_align_t *ret);

int align_linear_bw(const char *ref, const char *seq, int len, int drop,
			int error_quota, struct extend_align_t *ret);

int align_banded_fw(const char *ref, const char *seq, int lref, int lseq,
		int width, int drop, int err_quota, struct array_2D_t *array,
		struct extend_align_t *ret);

int align_banded_bw(const char *ref, const char *seq, int lref, int lseq,
		int width, int drop, int err_quota, struct array_2D_t *array,
		struct extend_align_t *ret);

#endif
