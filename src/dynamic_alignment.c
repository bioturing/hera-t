#include <assert.h>
#include "dynamic_alignment.h"
#include "utils.h"

int sub_mat[5][5] = {
			{  1, -2, -2, -2, -1 },
			{ -2,  1, -2, -2, -1 },
			{ -2, -2,  1, -2, -1 },
			{ -2, -2, -2,  1, -1 },
			{ -1, -1, -1, -1, -1 }
		    };

int b2b_check_nocigar(const char *ref, const char *seq, int len, int *err_quota)
{
	int i, ret, s;
	ret = 0;
	for (i = 0; i < len; ++i) {
		s = sub_mat[nt4_table[(int)ref[i]]][nt4_table[(int)seq[i]]];
		ret += s;
		*err_quota -= (SUB_MAX - s);
		// Out of error quota, should end soon
		if (*err_quota < 0)
			return ret;
	}
	return ret;
}

int align_linear_fw(const char *ref, const char *seq, int len, int drop,
			int error_quota, struct extend_align_t *ret)
{
	int bscore, score, seq_len, i;
	score = bscore = 0; // best score
	seq_len = 0;   // best pos
	for (i = 0; i < len; ++i) {
		score += sub_mat[nt4_table[(int)ref[i]]][nt4_table[(int)seq[i]]];
		if (score >= bscore) {
			bscore = score;
			seq_len = i + 1;
		}
		if (bscore - score > drop) {
			++i;
			break;
		}
	}

	if (score >= len * SUB_MAX - error_quota && i == len) { // can extend to the break point
		ret->seq_len = len;
		ret->score = score;
		score = len * SUB_MAX - score;
	} else if (bscore >= seq_len * SUB_MAX - error_quota) { // only keep the peak
		ret->seq_len = seq_len;
		ret->score = bscore;
		score = seq_len * SUB_MAX - bscore;
	} else { // discard
		ret->score = ret->seq_len = 0;
		score = 0;
	}
	return score;
}

int align_linear_bw(const char *ref, const char *seq, int len, int drop,
			int error_quota, struct extend_align_t *ret)
{
	int bscore, score, seq_len, i;
	score = bscore = 0; // best score
	seq_len = 0;   // best pos
	for (i = 0; i < len; ++i) {
		score += sub_mat[nt4_table[(int)ref[-i - 1]]][nt4_table[(int)seq[-i - 1]]];
		if (score >= bscore) {
			bscore = score;
			seq_len = i + 1;
		}
		if (bscore - score > drop) {
			++i;
			break;
		}
	}

	if (score >= len * SUB_MAX - error_quota && i == len) { // can extend to the break point
		ret->seq_len = len;
		ret->score = score;
		score = len * SUB_MAX - score;
	} else if (bscore >= seq_len * SUB_MAX - error_quota) { // only keep the peak
		ret->seq_len = seq_len;
		ret->score = bscore;
		score = seq_len * SUB_MAX - bscore;
	} else { // discard
		ret->score = ret->seq_len = 0;
		score = 0;
	}
	return score;
}

int **get_2D(int m, int n)
{
	int **ret;
	int i;
	ret = malloc(m * sizeof(int *));
	for (i = 0; i < m; ++i)
		ret[i] = calloc(n, sizeof(int));
	return ret;
}

void destroy_2D(int **a, int m)
{
	if (!a) return;
	int i;
	for (i = 0; i < m; ++i)
		free(a[i]);
	free(a);
}

int align_banded_fw(const char *ref, const char *seq, int lref, int lseq,
			int width, int drop, int max_err, struct array_2D_t *tmp_array,
			struct extend_align_t *ret)
{
	assert(lref >= 0 && lseq >= 0);
	if (!lref || !lseq) {
		ret->score = 0;
		ret->ref_len = ret->seq_len = 0;
		return 0;
	}
	int i, k, l, r, bscore, ref_len, seq_len, escore, cur_score;
	int **f = (int **)resize_array_2D(tmp_array, lref + 1, lseq + 1, 4);
	for (i = 0; i < lref + 1; ++i)
		f[i][0] = i * GAP_E;
	for (i = 0; i < lseq + 1; ++i)
		f[0][i] = i * GAP_E;
	bscore = 0;
	ref_len = seq_len = 0;
	for (i = 0; i < lref; ++i) {
		l = __max(0, i - width);
		r = __min(lseq, i + width + 1);
		for (k = l; k < r; ++k)
			f[i + 1][k + 1] = f[i][k] + sub_mat[nt4_table[(int)ref[i]]][nt4_table[(int)seq[k]]];
		for (k = l; k + 1 < r; ++k)
			f[i + 1][k + 1] = __max(f[i + 1][k + 1],
				 		f[i][k + 1] + GAP_E);
		for (k = l + 1; k < r; ++k)
			f[i + 1][k + 1] = __max(f[i + 1][k + 1],
						f[i + 1][k] + GAP_E);
		cur_score = 0;
		for (k = l; k < r; ++k) {
			if (f[i + 1][k + 1] > bscore) {
				bscore = f[i + 1][k + 1];
				seq_len = k + 1;
				ref_len = i + 1;
			}
			if (f[i + 1][k + 1] > cur_score)
				cur_score = f[i + 1][k + 1];
		}
		if (bscore - cur_score > drop) {
			if (bscore >= seq_len * SUB_MAX - max_err) {
				ret->score = bscore;
				ret->seq_len = seq_len;
				ret->ref_len = ref_len;
			} else {
				ret->score = ret->seq_len = ret->ref_len = 0;
			}
			return ret->seq_len * SUB_MAX - ret->score;
		}
	}
	escore = f[0][lseq];
	ret->ref_len = 0;
	for (i = 0; i < lref; ++i) {
		if (i + width + 1 < lseq) continue;
		if (f[i + 1][lseq] > escore) {
			escore = f[i + 1][lseq];
			ret->ref_len = i + 1;
		}
	}
	if (escore >= lseq * SUB_MAX - max_err) {
		ret->score = escore;
		ret->seq_len = lseq;
	} else if (bscore >= seq_len * SUB_MAX - max_err) {
		ret->score = bscore;
		ret->seq_len = seq_len;
		ret->ref_len = ref_len;
	} else {
		ret->score = ret->seq_len = ret->ref_len = 0;
	}
	return ret->seq_len * SUB_MAX - ret->score;
}

int align_banded_bw(const char *ref, const char *seq, int lref, int lseq,
			int width, int drop, int max_err, struct array_2D_t *tmp_array,
			struct extend_align_t *ret)
{
	assert(lref >= 0 && lseq >= 0);
	if (!lref || !lseq) {
		ret->score = 0;
		ret->ref_len = ret->seq_len = 0;
		return 0;
	}
	int i, k, l, r, bscore, ref_len, seq_len, escore, cur_score;
	int **f = (int **)resize_array_2D(tmp_array, lref + 1, lseq + 1, 4);
	for (i = 0; i < lref + 1; ++i)
		f[i][0] = i * GAP_E;
	for (i = 0; i < lseq + 1; ++i)
		f[0][i] = i * GAP_E;
	bscore = 0;
	ref_len = seq_len = 0;
	for (i = 0; i < lref; ++i) {
		l = __max(0, i - width);
		r = __min(lseq, i + width + 1);
		for (k = l; k < r; ++k)
			f[i + 1][k + 1] = f[i][k] +
				sub_mat[nt4_table[(int)ref[-i - 1]]][nt4_table[(int)seq[-k - 1]]];
		for (k = l; k + 1 < r; ++k)
			f[i + 1][k + 1] = __max(f[i + 1][k + 1],
				 		f[i][k + 1] + GAP_E);
		for (k = l + 1; k < r; ++k)
			f[i + 1][k + 1] = __max(f[i + 1][k + 1],
						f[i + 1][k] + GAP_E);
		cur_score = 0;
		for (k = l; k < r; ++k) {
			if (f[i + 1][k + 1] > bscore) {
				bscore = f[i + 1][k + 1];
				seq_len = k + 1;
				ref_len = i + 1;
			}
			if (f[i + 1][k + 1] > cur_score)
				cur_score = f[i + 1][k + 1];
		}
		if (bscore - cur_score > drop) {
			if (bscore >= seq_len * SUB_MAX - max_err) {
				ret->score = bscore;
				ret->seq_len = seq_len;
				ret->ref_len = ref_len;
			} else {
				ret->score = ret->seq_len = ret->ref_len = 0;
			}
			return ret->seq_len * SUB_MAX - ret->score;
		}
	}
	escore = f[0][lseq];
	ret->ref_len = 0;
	for (i = 0; i < lref; ++i) {
		if (i + width + 1 < lseq) continue;
		if (f[i + 1][lseq] > escore) {
			escore = f[i + 1][lseq];
			ret->ref_len = i + 1;
		}
	}
	if (escore >= lseq * SUB_MAX - max_err) {
		ret->score = escore;
		ret->seq_len = lseq;
	} else if (bscore >= seq_len * SUB_MAX - max_err) {
		ret->score = bscore;
		ret->seq_len = seq_len;
		ret->ref_len = ref_len;
	} else {
		ret->score = ret->seq_len = ret->ref_len = 0;
	}
	return ret->seq_len * SUB_MAX - ret->score;
}

