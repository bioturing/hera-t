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
#include "khash.h"
#include "pthread_barrier.h"
#include "library_type.h"

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

#define THREAD_STACK_SIZE	16777216

#define RNA_PRIOR	4
#define TAG_PRIOR	2

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

/* concat two strings */
void concat_str(char *s1, int l1, char *s2, int l2);

/* check if string contain N */
int check_valid_nu(const char *seq, int len);

/* check if string contain polyA */
int check_polyA(char *str, int len);

extern int8_t nt4_table[256];
extern char *nt4_char, *rev_nt4_char;

#endif /* _UTILS_H_ */
