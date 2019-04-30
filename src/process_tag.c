#include <pthread.h>

#include "process_tag.h"
#include "readTSV.h"
#include "attribute.h"
#include "verbose.h"
#include "utils.h"

#define MAX_ERROR 	1

KHASH_MAP_INIT_INT(tag, int);

struct tag_ref_t {
	khash_t(tag) *h;

	int trim;
	char *left_pat;
	char *right_pat;
	int left_len;
	int right_len;
	int ref_len;
};

struct tag_stat_t {
	int64_t nread;
	int64_t map;
	int64_t unmap;
};

const uint64_t _pow5_r[] = {1ull, 5ull, 25ull, 125ull, 625ull, 3125ull, 15625ull, 78125ull, 390625ull, 1953125ull, 9765625ull, 48828125ull, 244140625ull, 1220703125ull, 6103515625ull, 30517578125ull, 152587890625ull, 762939453125ull, 3814697265625ull, 19073486328125ull, 95367431640625ull, 476837158203125ull, 2384185791015625ull, 11920928955078125ull, 59604644775390625ull, 298023223876953125ull, 1490116119384765625ull};

static struct tag_ref_t *tag_ref;

int check_left_pattern(struct read_t *read, int start, char *pattern, int l)
{
	int i, k, error, max_iter;
	
	max_iter = read->len - tag_ref->ref_len - l + 1;

	if (max_iter < start)
		return -1;

	for (i = start; i < max_iter; ++i) {
		for (k = error = 0; k < l && error <= MAX_ERROR; ++k)
			if (nt4_table[(int)read->seq[i + k]] != nt4_table[(int)pattern[k]])
				++error;
		if (k == l)
			return start + l;
	}

	return -1;
}

int check_right_pattern(struct read_t *read, int start, char *pattern, int l)
{
	int i, k, error, max_iter;
	
	max_iter  = read->len - l + 1;
	start += tag_ref->ref_len;

	if (max_iter < start)
		return -1;

	for (i = start; i < max_iter; ++i) {
		for (k = error = 0; k < l && error <= MAX_ERROR; ++k)
			if (nt4_table[(int)read->seq[i + k]] != nt4_table[(int)pattern[k]])
				++error;
		if (k == l)
			return i - tag_ref->ref_len;
	}

	return -1;
}

/*****************************************
*            BUILD REF INDEX             *
*****************************************/

struct tag_ref_t *init_reference()
{
	struct tag_ref_t *ref = malloc(sizeof(struct tag_ref_t));
	ref->h = kh_init(tag);
	ref->trim = ref->ref_len = ref->left_len = ref->right_len = 0;
	ref->left_pat = ref->right_pat = NULL;

	return ref;
}

void get_col_idx(int *col_idx, struct tsv_t *f)
{
	ssize_t ret = get_row_content(f);
	if (!ret)
		__ERROR("Reference file is empty.\n");
	int i;
	struct content_t value;

	memset(col_idx, -1, 4 * sizeof(int));
	for (i = 0; i < f->l; ++i) {
		value = get_col_content(f, i);
		if (!strncmp(value.s, "id", value.l)) {
			col_idx[0] = i;
		} else if(!strncmp(value.s, "name", value.l)) {
			col_idx[1] = i;
		} else if(!strncmp(value.s, "pattern", value.l)) {
			col_idx[2] = i;
		} else if(!strncmp(value.s, "sequence", value.l)) {
			col_idx[3] = i;
		}
	}

	if (col_idx[2] == -1)
		__ERROR("Can not find column 'pattern'. Please re-check the input file header.\n");
	if (col_idx[3] == -1)
		__ERROR("Can not find column 'sequence'. Please re-check the input file header.\n");

	if (col_idx[0] == -1)
		col_idx[0] = col_idx[1];
	if (col_idx[1] == -1)
		col_idx[1] = col_idx[0];

	if (col_idx[0] == -1)
		__ERROR("Can not found column 'id' or 'name'. Please re-check the input file header.\n");
}

void add_hash(khash_t(tag) *h, int32_t idx, int len, int order)
{
	khiter_t k;
	int32_t ret, i, tmp_idx, new_idx, ch, c;

	k = kh_get(tag, h, idx);

	if (k != kh_end(h))
		__ERROR("There is pair of reference tags that is 1-hamming-distance away from each other");

	k = kh_put(tag, h, idx, &ret);
	kh_value(h, k) = order;

	for (i = 0; i < len; ++i) {
		ch = (idx / _pow5_r[i]) % 5;
		tmp_idx = idx - ch * _pow5_r[i];
		for (c = 0; c < NNU; ++c) {
			if (ch == c)
				continue;
			new_idx = tmp_idx + _pow5_r[i] * c;

			k = kh_put(tag, h, new_idx, &ret);
			kh_value(h, k) = order;
		}
	}
}

void parse_pattern(struct tag_ref_t *ref, struct content_t pattern)
{
	int i, l, len;
	int trim = 0;

	len = pattern.l;
	i = 0;
	l = 0;
	if (!strncmp(pattern.s, "5P", 2)) {
		i = 2;
		while(i < len && pattern.s[i] == 'N') {
			++l;
			++i;
		}

		if (i == len)
			__ERROR("Invalid reference pattern %s\n", pattern.s);

		trim = l;
		len = len;
	}

	if (!strncmp(pattern.s + (len - 2), "3P", 2)) {
		if (ref->trim)
			__ERROR("Should only 3P or 5P appear in reference pattern: %s\n", pattern.s);
		i = len - 1;
		while (i >= 0 && pattern.s[i] == 'N') {
			++l;
			--i;
		}

		if (i < 0)
			__ERROR("Invalid reference pattern %s\n", pattern.s);

		trim = -l;
		len -= i;
	}

	if (trim) {
		if (!ref->trim)
			ref->trim = trim;
		else if (ref->trim != trim)
			__ERROR("There are difference types of pattern.\n");
	}

	l = 0;
	while (i + l + 4 < len && pattern.s[i + l] != '(')
		++l;
	if (i + l == len)
		__ERROR("Invalid reference pattern %s\n", pattern.s);

	if (l) {
		if (!ref->left_len) {
			ref->left_pat = malloc(l);
			memcpy(ref->left_pat, pattern.s + i, l);
			ref->left_len = l;
		} else if ((ref->left_len != l ||
			strncmp(ref->left_pat, pattern.s + i, l))) {
			__ERROR("There are difference types of pattern.\n");
		}
	}

	i += l;
	if (strncmp(pattern.s + i, "(BC)", 4))
		__ERROR("Should only (BC) appears in pattern: %s\n", pattern.s);
	i += 4;
	l = 0;
	while (i + l < len && pattern.s[i + l] != 'N')
		++l;
	if (l) {
		if (!ref->right_len) {
			ref->right_pat = malloc(l);
			memcpy(ref->right_pat, pattern.s + i, l);
			ref->right_len = l;
		} else if ((ref->right_len != l ||
			strncmp(ref->right_pat, pattern.s + i, l))) {
			__ERROR("There are difference types of pattern.\n");
		}
	}
}

void add_ref(struct ref_info_t *ref, struct content_t id, struct content_t name)
{
	int n = ref->n_refs;
	++ref->n_refs;

	//id
	++id.l;
	ref->id_iter = realloc(ref->id_iter, ref->n_refs * sizeof(int));
	ref->id_iter[n] = ref->id_len;
	ref->id_len += id.l;
	ref->ref_id = realloc(ref->ref_id, ref->id_len);
	memcpy(ref->ref_id + ref->id_iter[n], id.s, id.l);

	//name
	++name.l;
	ref->text_iter = realloc(ref->text_iter, ref->n_refs * sizeof(int));
	ref->text_iter[n] = ref->text_len;
	ref->text_len += name.l;
	ref->ref_text = realloc(ref->ref_text, ref->text_len);
	memcpy(ref->ref_text + ref->text_iter[n], name.s, name.l);
}

void build_tag_ref(struct input_t *input, struct ref_info_t *ref, int type)
{
	tag_ref = init_reference();
	struct tsv_t *f = init_readTSV(input->ref, ',');
	struct content_t value;
	int col_idx[4]; //[id, name, pattern, sequence];

	get_col_idx(col_idx, f);

	while(get_row_content(f) > 0) {
		value = get_col_content(f, col_idx[3]);
		if (tag_ref->ref_len && tag_ref->ref_len != value.l)
			__ERROR("Reference sequence lengths are not equal");

		if (!check_valid_nu(value.s, value.l))
			__ERROR("Reference tags sequence is not valid %s", value.s);

		tag_ref->ref_len = value.l;
		add_hash(tag_ref->h, seq2num(value.s, value.l), value.l, ref->n_refs);
		parse_pattern(tag_ref, get_col_content(f, col_idx[2]));
		add_ref(ref, get_col_content(f, col_idx[0]),
			get_col_content(f, col_idx[1]));
	}

	ref->type[type] = ref->n_refs;
	destroy_readTSV(f);
};

/*****************************************
*                 MAP TAGS               *
*****************************************/

struct tag_stat_t *tag_count;

void init_tag_threads(int n_threads)
{
	tag_count = calloc(n_threads + 1, sizeof(struct tag_stat_t));
}

int32_t get_tag_idx(struct read_t *read)
{
	int start;
	if (tag_ref->trim >=0)
		start = tag_ref->trim;
	else
		start = read->len + tag_ref->trim + 1 -
			(tag_ref->left_len + tag_ref->right_len + tag_ref->ref_len);
	
	if (tag_ref->left_len)
		start = check_left_pattern(read, start, tag_ref->left_pat, tag_ref->left_len);
	else if (tag_ref->right_len)	
		start = check_right_pattern(read, start, tag_ref->right_pat, tag_ref->right_len);

	if (start == -1)
		return -1;
	return seq2num(read->seq + start, tag_ref->ref_len);
}

int align_tag(struct read_t *read, int thread_num)
{
	struct tag_stat_t *c = tag_count + (thread_num + 1);

	int32_t tag_idx;
	khiter_t k;

	++c->nread;

	tag_idx = get_tag_idx(read);
	if (tag_idx != -1) {
		k = kh_get(tag, tag_ref->h, tag_idx);
		if (k != kh_end(tag_ref->h)) {
			++c->map;
			return kh_value(tag_ref->h, k);
		}
	}

	++c->unmap;
	return -1;
}

void update_tag_result(struct tag_stat_t *res, struct tag_stat_t *add)
{
	res->nread += add->nread;
	res->map += add->map;
	res->unmap += add->unmap;
}

void print_tag_count(int thread_num)
{
	struct tag_stat_t *c = tag_count + (thread_num + 1);
	update_tag_result(tag_count, c);
	memset(c, 0, sizeof(struct tag_stat_t));
	__VERBOSE("\r Mapped reads: %ld / %ld",
			 tag_count[0].map, tag_count[0].nread);
}

void print_tag_stat(int n_threads)
{
	int i;
	for (i = 1; i <= n_threads; ++i)
		update_tag_result(tag_count, tag_count + i);
	__VERBOSE("\n");
	__VERBOSE_LOG("INFO", "Total number of reads        : %ld\n", tag_count[0].nread);
	__VERBOSE_LOG("INFO", "Number of mapped reads       : %ld\n", tag_count[0].map);
	__VERBOSE_LOG("INFO", "Number of unmapped reads     : %ld\n", tag_count[0].unmap);
}

void destroy_tag_ref()
{
	kh_destroy(tag, tag_ref->h);
	free(tag_ref);
	free(tag_count);
}