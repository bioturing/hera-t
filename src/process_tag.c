#include <pthread.h>

#include "atomic.h"

#include "process_tag.h"
#include "readTSV.h"
#include "attribute.h"
#include "verbose.h"
#include "utils.h"
#include "mini_hash.h"
#include "minimizers.h"

#define MAX_ERROR 	1

KHASH_MAP_INIT_INT(tag, int);
#define MINIMIZER_KMER 4
#define MINIMIZER_WINDOW 3
#define MIN_MM_MAP 1
#define MIN_REF_BASE 10

static int is_10x_data;

struct tag_ref_t {
	khash_t(tag) *h;

	int trim;
	char *left_pat;
	char *right_pat;
	int left_len;
	int right_len;
	int ref_len[2];
	struct mini_hash_t *mh;
};

struct tag_stat_t {
	uint64_t nread;
        uint64_t map;
        uint64_t unmap;
};

const uint64_t _pow5_r[] = {1ull, 5ull, 25ull, 125ull, 625ull, 3125ull, 15625ull, 78125ull, 390625ull, 1953125ull, 9765625ull, 48828125ull, 244140625ull, 1220703125ull, 6103515625ull, 30517578125ull, 152587890625ull, 762939453125ull, 3814697265625ull, 19073486328125ull, 95367431640625ull, 476837158203125ull, 2384185791015625ull, 11920928955078125ull, 59604644775390625ull, 298023223876953125ull, 1490116119384765625ull};

static struct tag_ref_t *tag_ref;

void check_left_pattern(struct read_t *read, int *range)
{
	int i, k, len, error, max_iter;
	char *pattern = tag_ref->left_pat;
	range[0] = range[1] = -1;

	if (tag_ref->trim >= 0) {
		i = tag_ref->trim;
		len = i + tag_ref->left_len + tag_ref->right_len + tag_ref->ref_len[1];
	} else {
		len = read->len + tag_ref->trim;
		i = len - tag_ref->left_len - tag_ref->right_len - tag_ref->ref_len[1];
	}

	if (len < tag_ref->left_len + tag_ref->ref_len[0])
		return;

	max_iter = len - tag_ref->ref_len[0] - tag_ref->left_len;
	for (; i <= max_iter; ++i) {
		for (k = error = 0; k < tag_ref->left_len && error <= MAX_ERROR; ++k)
			if (nt4_table[(int)read->seq[i + k]] != nt4_table[(int)pattern[k]])
				++error;
		if (k == tag_ref->left_len) {
			range[0] = i + tag_ref->left_len;
			range[1] = len;
			return;
		}
	}
}

void check_right_pattern(struct read_t *read, int *range)
{
	int i, k, len, error, max_iter;
	char *pattern = tag_ref->right_pat;
	range[0] = range[1] = -1;

	if (tag_ref->trim >= 0) {
		i = tag_ref->trim;
		len = i + tag_ref->left_len + tag_ref->right_len + tag_ref->ref_len[1];
	} else {
		len = read->len + tag_ref->trim;
		i = len - tag_ref->left_len - tag_ref->right_len - tag_ref->ref_len[1];
	}

	if (len < tag_ref->right_len + tag_ref->ref_len[0])
		return ;

	max_iter  = read->len - tag_ref->right_len;
	i += tag_ref->ref_len[0];

	for (; i <= max_iter; ++i) {
		for (k = error = 0; k < tag_ref->right_len && error <= MAX_ERROR; ++k)
			if (nt4_table[(int)read->seq[i + k]] != nt4_table[(int)pattern[k]])
				++error;
		if (k == tag_ref->right_len) {
			range[0] = -(i - tag_ref->ref_len[1]);
			range[1] = i;
			return;
		}
	}

}

/*****************************************
*            BUILD REF INDEX             *
*****************************************/

struct tag_ref_t *init_reference()
{
	struct tag_ref_t *ref = malloc(sizeof(struct tag_ref_t));
	ref->h = kh_init(tag);
	ref->trim = ref->left_len = ref->right_len = 0;
	ref->ref_len[0] = ref->ref_len[1] = 0;
	ref->left_pat = ref->right_pat = NULL;
	init_mini_hash(&ref->mh, 3);

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

	/* Add 1-hamming distance tag */
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
	if (!strncmp(pattern.s, "5P", 2) || !strncmp(pattern.s, "^", 1)) {
		i = (!strncmp(pattern.s, "5P", 2)) ? 2 : 1;
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

void add_tag_len(int len)
{
	if (!tag_ref->ref_len[0]) {
		tag_ref->ref_len[0] = tag_ref->ref_len[1] = len;
		return;
	}

	tag_ref->ref_len[0] = __min(tag_ref->ref_len[0], len);
	tag_ref->ref_len[1] = __max(tag_ref->ref_len[1], len);
}

int add_minimizer(struct mini_hash_t **h, char *s, int l, int tag)
{
	struct mm_db_t *db = mm_index_char_str(s, MINIMIZER_KMER, MINIMIZER_WINDOW, l);
	int i;
	uint64_t *slot;
	for (i = 0; i < db->n; ++i) {
		slot = mini_put(h, db->mm[i]);
		if (*slot != 0)
			*slot = TOME_STONE;
		else
			*slot = tag + 1; //Avoid confusion of zero: empty slot or the first tag
	}
	return db->n;
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
		add_tag_len(value.l);

		if (!check_valid_nu(value.s, value.l))
			__ERROR("Reference tags sequence is not valid %s", value.s);

		add_hash(tag_ref->h, seq2num(value.s, value.l), value.l, ref->n_refs);
		if (value.l > MIN_REF_BASE) //tag ref must be long enough to be index by minimizers
			add_minimizer(&tag_ref->mh, value.s, value.l, ref->n_refs);
		parse_pattern(tag_ref, get_col_content(f, col_idx[2]));
		add_ref(ref, get_col_content(f, col_idx[0]),
			get_col_content(f, col_idx[1]));
	}

	mm_stats(tag_ref->mh);
	ref->type[type] = ref->n_refs;
	destroy_readTSV(f);
};

void mm_stats(struct mini_hash_t *h)
{
	int i;
	int n = 0, single = 0;
	for (i = 0; i < h->size; ++i) {
		if (h->key[i] != EMPTY_SLOT) ++n;
		if (h->h[i] != 0 && h->h[i] != TOME_STONE) ++single;
	}
	__VERBOSE("Total of minimizers: %d\n", n);
	__VERBOSE("Total of single-ton minimizers: %d\n", single);
}

/*****************************************
*                 MAP TAGS               *
*****************************************/

struct tag_stat_t tag_count;

void get_tag_range(struct read_t *read, int *range)
{
	if (tag_ref->left_len) {
		check_left_pattern(read, range);
	} else if (tag_ref->right_len) {
		check_right_pattern(read, range);
	} else if (tag_ref->trim >= 0) {
		range[0] = tag_ref->trim;
		range[1] = read->len;
	} else {
		range[1] = read->len + tag_ref->trim;
		range[0] = -(range[1] - tag_ref->ref_len[1]);
	}
}

int mm_map(struct mm_db_t *db, struct mini_hash_t *h)
{
	int i, n_map = 0;
	uint64_t ref = 0;
	uint64_t *slot;

	i = 0;
	while ((ref == 0 || ref == TOME_STONE) && i < db->n)  {
		slot =  mini_get(h, db->mm[i]);
		if (slot != (uint64_t *)EMPTY_SLOT && *slot != TOME_STONE)
			ref = *slot;
		++i;
	}

	if (ref == 0 || ref == TOME_STONE)
		return 0;
	else
		n_map = 1;

	for (; i < db->n; ++i) {
		slot = mini_get(h, db->mm[i]);
		if (slot != (uint64_t *)EMPTY_SLOT) {
			if (*slot == TOME_STONE)
				continue;
			if (*slot != ref) {
				return 0;
			}
			++n_map;
		}
	}

	if (n_map > MIN_MM_MAP)
		return ref;
	return 0;
}

int align_tag(struct read_t *read, int thread_num)
{

	int32_t tag_idx, i;
	int32_t range[2];
	khiter_t k;

	atomic_add_and_fetch64(&(tag_count.nread), 1);

	get_tag_range(read, range);

	if (range[0] == -1 || range[1] == -1) {
		atomic_add_and_fetch64(&(tag_count.unmap), 1);
		return -1;
	}

	if (range[0] >= 0) {
		for (i = tag_ref->ref_len[0]; i <= tag_ref->ref_len[1]; ++i) {
			if (range[0] + i >= range[1])
				break;
			tag_idx = seq2num(read->seq + range[0], i);
			k = kh_get(tag, tag_ref->h, tag_idx);
			if (k != kh_end(tag_ref->h)) {
				atomic_add_and_fetch64(&(tag_count.map), 1);
				return kh_value(tag_ref->h, k);
			}
		}
	} else {
		for (i = tag_ref->ref_len[1]; i >= tag_ref->ref_len[0]; --i) {
			if (range[1] - i < range[0])
				break;
			tag_idx = seq2num(read->seq - range[0], i);
			k = kh_get(tag, tag_ref->h, tag_idx);
			if (k != kh_end(tag_ref->h)) {
				atomic_add_and_fetch64(&(tag_count.map), 1);
				return kh_value(tag_ref->h, k);
			}
			--range[0];
		}
	}
	// For unmap read in 10x data
	if (i > MIN_REF_BASE) {
		struct mm_db_t *db = mm_index_char_str(read->seq + range[0], MINIMIZER_KMER, MINIMIZER_WINDOW, i);
		int ref = mm_map(db, tag_ref->mh);
		mm_db_destroy(db);
		if (ref != 0) {
			atomic_add_and_fetch64(&(tag_count.map), 1);
			return ref - 1; // tag index for minimizer map is 1-based
		}
	}

	atomic_add_and_fetch64(&(tag_count.unmap), 1);
	return -1;
}

void init_tag_threads(int n_threads)
{
	tag_count.nread = 0;
	tag_count.map = 0;
	tag_count.unmap = 0;
}

void print_tag_count(int thread_num)
{
	__VERBOSE("Mapped %llu / %llu\n", tag_count.map, tag_count.nread);
}

void print_tag_stat(int n_threads)
{
	__VERBOSE("\n");
	__VERBOSE_LOG("INFO", "Total number of reads        : %llu\n", tag_count.nread);
	__VERBOSE_LOG("INFO", "Number of mapped reads       : %llu\n", tag_count.map);
	__VERBOSE_LOG("INFO", "Number of unmapped reads     : %llu\n", tag_count.unmap);
}

void destroy_tag_ref()
{
	kh_destroy(tag, tag_ref->h);
	free(tag_ref);
}
