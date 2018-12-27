#ifndef _RADIX_SORT_H_
#define _RADIX_SORT_H_

#define RS_MIN_SIZE		64

#define RS_PROTO(name, rs_type)						       \
void rs_sort_##name(rs_type *beg, rs_type *end);

#define RS_IMPL(name, rs_type, width, block, less_than, get_block)	       \
struct rs_bucket_##name {						       \
	rs_type *b, *e;							       \
}; \
static inline void is_sort_##name(rs_type *beg, rs_type *end) \
{ \
	rs_type *i, *j, tmp; \
	for (i = beg + 1; i < end; ++i) { \
		if (less_than(*i, *(i - 1))) { \
			tmp = *i; \
			for (j = i; j > beg && less_than(tmp, *(j - 1)); --j) \
				*j = *(j - 1); \
			*j = tmp; \
		} \
	} \
} \
static void recursive_sort_##name(rs_type *beg, rs_type *end, int n_bits, int s) \
{ \
	rs_type *i; \
	int size = 1 << n_bits, m = size - 1; \
	struct rs_bucket_##name *k; \
	struct rs_bucket_##name *b = \
		(struct rs_bucket_##name*)alloca(size * sizeof(struct rs_bucket_##name)); \
	struct rs_bucket_##name *be = b + size; \
	for (k = b; k != be; ++k) k->b = k->e = beg; \
	for (i = beg; i != end; ++i) ++b[get_block(*i, s, m)].e; \
	for (k = b + 1; k != be; ++k) \
		k->e += (k - 1)->e - beg, k->b = (k - 1)->e; \
	for (k = b; k != be;) { \
		if (k->b != k->e) { \
			struct rs_bucket_##name *l; \
			if ((l = b + get_block(*k->b, s, m)) != k) { \
				rs_type tmp = *k->b, swap; \
				do { \
					swap = tmp; \
					tmp = *l->b; \
					*l->b++ = swap; \
					l = b + get_block(tmp, s, m); \
				} while (l != k); \
				*k->b++ = tmp; \
			} else ++k->b; \
		} else ++k; \
	} \
	for (b->b = beg, k = b + 1; k != be; ++k) k->b = (k - 1)->e; \
	if (s) { \
		s -= n_bits; \
		for (k = b; k != be; ++k) { \
			if (k->e - k->b > RS_MIN_SIZE) \
				recursive_sort_##name(k->b, k->e, n_bits, s); \
			else if (k->e - k->b > 1) \
				is_sort_##name(k->b, k->e); \
		} \
	} \
} \
void rs_sort_##name(rs_type *beg, rs_type *end) \
{ \
	if (end - beg > RS_MIN_SIZE) \
		recursive_sort_##name(beg, end, block, width - block); \
	else \
		is_sort_##name(beg, end); \
}

#define rs_sort(name, beg, end) rs_sort_##name(beg, end)
#endif
