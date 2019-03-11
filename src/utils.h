#ifndef _UTILS_H_
#define _UTILS_H_

#if defined(_MSC_VER)
#pragma warning(disable:4996)
#endif

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <pthread.h>

#include "align_attr.h"
#include "attribute.h"
#include "dqueue.h"
#include "get_buffer.h"
#include "interval_tree.h"
#include "kmhash.h"
#include "khash.h"
#include "pthread_barrier.h"
#include "library_type.h"
#include "antibody.h"

#if defined(_MSC_VER)
#include <time.h>
#include <windows.h>
#include <getopt.h>
#include <BaseTsd.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#endif /* _MSC_VER */

#define MAX_INT32		2147483647
#define MIN_INT32		-2147483648

#define MASK32			4294967295ULL

#define BUFSZ			4096

#define THREAD_STACK_SIZE	16777216

#define FORWARD			0
#define REVERSE			1
#define LEFT			0
#define RIGHT			1

//#define WRITE_TEXT
//#define RUN_AFTER_ANALYSIS

/*
 * Built in macros
 */

#define __abs(x) 		((x) < 0 ? -(x) : (x))

#define __min(a, b) 		((a) < (b) ? (a) : (b))

#define __max(a, b) 		((a) > (b) ? (a) : (b))

#define __min3(a, b, c)		__min(__min((a), (b)), (c))

#define __max3(a, b, c)		__max(__max((a), (b)), (c))

#define __round_up_32(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16, ++(x))

#define __is_sep(c)		((c) == ' ' || (c) == '\t')

#define normalize_mapq(s)	do {					       \
	if ((s) < 0) (s) = 0;						       \
	if ((s) > 60) (s) = 60;						       \
} while (0)

/*
 * Built-in macros function
 */

#define __ALLOC(ptr, sz)	(ptr) = xmalloc(sizeof(*(ptr)) * (sz))

#define __REALLOC(ptr, sz)	(ptr) = xrealloc((ptr), sizeof(*(ptr)) * (sz))

/* push back val to ptr, ptr has sz element, realloc + 1 */
#define __PUSH_BACK(ptr, sz, val) do {					       \
	assert((sz) >= 0);						       \
	__REALLOC((ptr), (sz) + 1);					       \
	(ptr)[(sz)++] = (val);						       \
} while(0)

#define __FREE_AND_NULL(ptr) do {					       \
	free(p);							       \
	(p) = NULL;							       \
} while (0)

#define __SWAP(x, y) do {						       \
	assert(sizeof(x) == sizeof(y));					       \
	int8_t temp[sizeof(x)];						       \
	memcpy(temp, &(y), sizeof(x));					       \
	memcpy(&(y), &(x), sizeof(x));					       \
	memcpy(&(x), temp, sizeof(x));					       \
} while (0)

/*
 * Built-in function
 */

/* get time */
double realtime();

/* reverse compelemnt */
char *get_rev_complement(const char *seq, int len);

/* reverse string */
char *get_rev(const char *seq, int len);

/* return new char* concate s1 and s2 */
char *str_concate(const char *s1, const char *s2);

/* convert from [ACGTN]+ seq to number */
int64_t seq2num(const char *seq, int len);

/* convert from number to [ACGTN]+ to number */
char *num2seq(int64_t num, int len);

int check_valid_nu(const char *seq, int len);
/*
 * Global variable
 */

extern int8_t nt4_table[256];
extern char *nt4_char, *rev_nt4_char;

struct array_2D_t {
	void *data;
	int len;
	int nrow;
	void **rows;
};

struct producer_bundle_t {
	int *n_consumer;
	void *stream;
	pthread_barrier_t *barrier;
	pthread_mutex_t *lock;
	struct dqueue_t *q;
};

struct worker_bundle_t {
	struct dqueue_t *q;
	struct kmhash_t *bc_table;
	pthread_mutex_t *lock_count;
	pthread_mutex_t *lock_hash;
	struct align_stat_t *result;
	struct raw_alg_t *alg_array;
	struct interval_t *intron_array;
	struct recycle_bin_t *recycle_bin;
	struct seed_t *seed_cons;
	struct array_2D_t *tmp_array;
	struct shared_fstream_t *align_fstream;
	struct library_t lib;
	struct antibody_lib_t *antibody_lib;
	// Function pointer
	void (*init_bundle)(struct worker_bundle_t*);
	void (*destroy_bundle)(struct worker_bundle_t*);
	void (*map_read)(struct read_t*, struct read_t*, struct worker_bundle_t*);
};

struct worker_data_t {
	struct kmhash_t *bc_table;
	struct opt_count_t *opt;
	struct dqueue_t *q;
	struct worker_bundle_t *worker_bundles;
	pthread_t *worker_threads;
	pthread_t *producer_threads;
	pthread_mutex_t *lock_count;
	pthread_attr_t *attr;
	struct align_stat_t result;
	struct antibody_lib_t *lib;
};

struct quant_data_t {
	int n_files;
	char **left_file;
	char **right_file;
	struct antibody_lib_t *lib;
	// Function pointer
	void (*action)(struct worker_data_t);
};

struct pair_buffer_t {
	char *buf1;
	char *buf2;
	int input_format;
};

struct dqueue_t *init_dqueue_PE(int cap);

struct pair_buffer_t *init_pair_buffer();

void free_pair_buffer(struct pair_buffer_t *p);

void **resize_array_2D(struct array_2D_t *p, int nrow, int ncol, int word);

void init_bundle(struct worker_bundle_t *bundle);

void destroy_bundle(struct worker_bundle_t *bundle);

void reinit_bundle(struct worker_bundle_t *bundle);

#endif /* _UTILS_H_ */
