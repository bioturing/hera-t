//
// Created by BioTuring on 2019-11-08.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <limits.h>
#include <inttypes.h>

#include "atomic.h"
#include "utils.h"
#include "minimizers.h"
#include "mini_hash.h"

#ifdef DEBUG
#define DEBUG_MM printf
#else
#define DEBUG_MM
#endif

void mm_db_insert(struct mm_db_t *db, uint64_t km, uint32_t p)
{
	if (db->n == db->size) {
		db->mm =  realloc(db->mm, (db->size << 1) * sizeof(uint64_t));
		db->p =  realloc(db->p, (db->size << 1) * sizeof(uint64_t));
		db->size <<= 1;
	}
	db->mm[db->n] = km;
	db->p[db->n++] = p;
}

struct mm_db_t * mm_db_init()
{
	struct mm_db_t *db = calloc(1, sizeof(struct mm_db_t));
	db->n = 0;
	db->size = 8;
	db->mm = calloc(db->size, sizeof(struct mm_db_t));
	db->p = calloc(db->size, sizeof(struct mm_db_t));
	return db;
}

void mm_db_destroy(struct mm_db_t *db)
{
	free(db->mm);
	free(db->p);
	free(db);
}

/**
 * Get the i kmer of the DNA sequence s
 * @param s DNA sequence s
 * @param i pos
 * @param k k size
 * @return encoded kmer
 */
static inline uint64_t get_km_i_str(char *s, int i, int k)
{
	uint64_t km, c;
	int j;
	int pad = (32 - k - 1)*2;
	km = 0;
	for (j = 0; j < k; ++j) {
		c = (uint64_t)nt4_table[s[i + j]];
		km |= c;
		km <<= 2;
	}
	km <<= pad;
	return km;
}

/**
 * Helper function that dump the true DNA sequence of an encoded sequence s
 * @param s binary encoded sequence s
 * @param l length of the encoded sequence s
 * @return true DNA sequence
 */
static inline char *mm_dump_seq(uint64_t s, uint32_t l)
{
	char *seq = calloc(l, sizeof(char));
	for (int i = 0; i < l; ++i) {
		uint8_t c = (uint8_t)((s >> (64 - i*2 - 2)) & (uint64_t)3);
		seq[i] = nt4_char[c];
	}
	return seq;
}

#define HASH64(k) twang_mix64(k)
#define DEBUG_PRINT

/**
 * The core function that index minimizers for a given string s
 * @param s DNA sequence s
 * @param k k size for minimizer indexing
 * @param w window size for minimizer indexing
 * @param l length of the DNA sequence s
 * @return minimizer collection of the string s
 */
struct mm_db_t * mm_index_char_str(char *s, int k, int w, int l)
{
	struct mm_db_t *db = mm_db_init();
	db->k = k;

	int i, j, p = -1;
	uint64_t km, mm, c;
	uint64_t km_h, mm_h;
	int pad = (32 - k - 1)*2;

	uint64_t tmp = (uint64_t)k;
	mm_h = km_h = HASH64(tmp);
	for (i = 0; i < l - w + 1; ++i) {
		DEBUG_PRINT("[i = %d]\n", i);
		if (i + w + k - 1 >= l)
			break;
		if (p < i) {
			km = mm = get_km_i_str(s, i, k);
			mm_h = km_h = HASH64(mm);
			p = i;
			for (j = 0; j < w; ++j) {
				c = (uint64_t) nt4_table[s[i + j + k - 1]];
				km |= ((uint64_t) c << (pad + 2));
				km_h = HASH64(km);
				if (km_h < mm_h) {
					mm = km;
					mm_h = km_h;
					p = i + j;
					mm_db_insert(db, km, p);
				}
				km <<= 2;
			}
			DEBUG_PRINT("[1]minimizers at window %d: %d\n", i, p);

			continue;
		} else {
			c = (uint64_t) nt4_table[s[i + w + k - 2]];
			km |= ((uint64_t) c << (pad + 2));
			km_h = HASH64(km);
			if (km_h < mm_h){
				p = i + w - 1;
				mm = km;
				mm_h = km_h;
				mm_db_insert(db, km, p);
				DEBUG_PRINT("[2]minimizers at window %d: %d\n", i, p);
			}
			km <<= 2;
		}
	}
	//mm_print(db);
	return db;
}
