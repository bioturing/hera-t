//
// Created by BioTuring on 2020-02-01.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>

#include "atomic.h"
#include "verbose.h"

#include "mini_hash.h"

/*
 Credit for primes table: Aaron Krowne
 http://br.endernet.org/~akrowne/
 http://planetmath.org/encyclopedia/GoodHashTablePrimes.html
 */
static const uint64_t primes[] = { 53, 97, 193, 389, 769, 1543, 3079, 6151,
                                   12289, 24593, 49157, 98317, 196613, 393241, 786433, 1572869, 3145739,
                                   6291469, 12582917, 25165843, 50331653, 100663319, 201326611, 402653189,
                                   805306457, 1610612741 };
#define N_PRIMES_NUMBER 26

#define MAX_LOAD_FACTOR 0.8
#define FATAL_LOAD_FACTOR 0.9

pthread_mutex_t h_table_mut;

/**
 * @brief Init a hash table with a pre-define size in the prime number table
 * @param p_index index of the prime number table
 * @return
 */
void init_mini_hash(struct mini_hash_t **h_table, uint32_t p_index)
{
	struct mini_hash_t *table = calloc(1, sizeof(struct mini_hash_t));
	table->prime_index = p_index;
	uint32_t h_size = primes[table->prime_index];
	table->h = calloc(h_size, sizeof(uint64_t));
	table->key = calloc(h_size, sizeof(uint64_t));
	memset(table->key, 255, sizeof(uint64_t) * h_size);
	table->size = h_size;
	table->max_cnt = (uint64_t) (table->size * MAX_LOAD_FACTOR);
	table->count = 0;
	*h_table = table;
}

void destroy_mini_hash(struct mini_hash_t *h_table)
{
	/** for (int i = 0; i < h_table->size; ++i){
		if (h_table->h[i] != EMPTY_BX)
			free((struct mini_hash_t *) h_table->h[i]);
	}
	 **/
	free(h_table->h);
	free(h_table->key);
}

/*
 * Thomas Wang 64 bit mix hash function
 */

uint64_t twang_mix64(uint64_t key) {
	key = (~key) + (key << 21); // key *= (1 << 21) - 1; key -= 1;
	key = key ^ (key >> 24);
	key = key + (key << 3) + (key << 8); // key *= 1 + (1 << 3) + (1 << 8)
	key = key ^ (key >> 14);
	key = key + (key << 2) + (key << 4); // key *= 1 + (1 << 2) + (1 << 4)
	key = key ^ (key >> 28);
	key = key + (key << 31); // key *= 1 + (1 << 31)
	return key;
}

/*
 * Inverse of twang_mix64
 *
 * Note that twang_unmix64 is significantly slower than twang_mix64.
 */

uint64_t twang_unmix64(uint64_t key) {
	// See the comments in jenkins_rev_unmix32 for an explanation as to how this
	// was generated
	key *= 4611686016279904257U;
	key ^= (key >> 28) ^ (key >> 56);
	key *= 14933078535860113213U;
	key ^= (key >> 14) ^ (key >> 28) ^ (key >> 42) ^ (key >> 56);
	key *= 15244667743933553977U;
	key ^= (key >> 24) ^ (key >> 48);
	key = (key + 1) * 9223367638806167551U;
	return key;
}

/**
 * @brief Expand the hash table by re-hash and doubling size
 * when the load factor reach MAX_LOAD_FACTOR (0.65 by default)
 * key must be hashed by MurMurHash64
 */
void mini_expand(struct mini_hash_t **new_table_ptr)
{
	uint32_t i;
	uint64_t *slot;
	uint64_t key, data, val;
	uint32_t log_size = 1<<16;
	struct mini_hash_t *new_table;
	struct mini_hash_t *h_table = *new_table_ptr;
	init_mini_hash(&new_table, h_table->prime_index + 1);
	assert(h_table->h != NULL);
	for (i = 0; i < h_table->size; ++i) {
		if (!((i + 1) % log_size)) {
			log_size <<= 1;
		}
		if (__mini_empty(h_table, i))
			continue;
		val = h_table->h[i];
		data = h_table->key[i];
		key = twang_mix64(data);
		slot = mini_put_by_key(new_table, data, key);
		*slot = val;
	}
	destroy_mini_hash(h_table);
	*new_table_ptr = new_table;
}

/**
 * @brief Try expanding the hash table by doubling in size
 */
inline void try_expanding(struct mini_hash_t **h_table)
{
	struct mini_hash_t *table = *h_table;
	if (table->count == table->max_cnt) {
		if (table->prime_index < N_PRIMES_NUMBER - 1) {
			mini_expand(h_table);
		} else {
			if ((double)table->count > FATAL_LOAD_FACTOR * (table->size)) {
				__VERBOSE("Hash table size reached limit!\n");
			}
		}
	}
}

/**
 * Increase the the count of one barcode by 1
 * @param data  barcode encoded as an uint64_t number
 * @param key   hash(data)
 */
uint64_t *mini_put_by_key(struct mini_hash_t *h_table, uint64_t data, uint64_t key)
{
	uint64_t i;
	uint64_t mask = h_table->size;
	uint64_t slot = key % mask;
	uint64_t is_empty = atomic_bool_CAS64(h_table->key + slot, EMPTY_SLOT, data);
	if (is_empty) { // slot is empty -> fill in
		atomic_add_and_fetch64(&(h_table->count), 1);
	} else if (!atomic_bool_CAS64(h_table->key + slot, data, data)) { // slot is reserved
		//linear probing
		for (i = slot + 1; i < h_table->size && !atomic_bool_CAS64(h_table->key + i, data, data); ++i) {
			is_empty = atomic_bool_CAS64(h_table->key + i, EMPTY_SLOT, data);
			if (is_empty)
				break;
		}
		if (i == h_table->size) {
			for (i = 0; i < slot && !atomic_bool_CAS64(h_table->key + i, data, data); ++i) {
				is_empty = atomic_bool_CAS64(h_table->key + i, EMPTY_SLOT, data);
				if (is_empty)
					break;
			}
		}
		assert(!atomic_bool_CAS64(&i, slot, slot));
		if (is_empty) //room at probe is empty -> fill in
			atomic_add_and_fetch64(&(h_table->count), 1);
		slot = i;
	}
	return h_table->h + slot;
}

/**
 * @param data  barcode encoded as an uint64_t number
 * @param key   hash(data)
 */
uint64_t *mini_get_by_key(struct mini_hash_t *h_table, uint64_t data, uint64_t key)
{
	uint64_t i;
	uint64_t mask = h_table->size;
	uint64_t slot = key % mask;
	uint64_t is_empty = atomic_bool_CAS64(h_table->key + slot, EMPTY_SLOT, EMPTY_SLOT);
	if (is_empty) { // slot is empty -> fill in
		return (uint64_t *)EMPTY_SLOT;
	} else if (!atomic_bool_CAS64(h_table->key + slot, data, data)) { // slot is reserved
		//linear probing
		for (i = slot + 1; i < h_table->size && !atomic_bool_CAS64(h_table->key + i, data, data); ++i) {
			is_empty = atomic_bool_CAS64(h_table->key + i, EMPTY_SLOT, EMPTY_SLOT);
			if (is_empty)
				break;
		}
		if (i == h_table->size) {
			for (i = 0; i < slot && !atomic_bool_CAS64(h_table->key + i, data, data); ++i) {
				is_empty = atomic_bool_CAS64(h_table->key + i, EMPTY_SLOT, EMPTY_SLOT);
				if (is_empty)
					break;
			}
		}
		assert(!atomic_bool_CAS64(&i, slot, slot));
		if (is_empty) //room at probe is empty -> fill in
			return (uint64_t *)EMPTY_SLOT;
		slot = i;
	}
	return h_table->h + slot;
}

/**
 * @brief
 * @param data
 * @param key
 * @return count of data
 */

uint64_t *mini_get(struct mini_hash_t *h_table, uint64_t data)
{
	uint64_t key = twang_mix64(data);
	uint64_t *slot = mini_get_by_key(h_table, data, key);
	return slot;
}

/**
 * @brief Increase the count of data to 1
 * @param data byte array of data
 * @param len length in byte of data
 */
uint64_t *mini_put(struct mini_hash_t **h_table, uint64_t data)
{
	struct mini_hash_t *table = *h_table;
	uint64_t key = twang_mix64(data);
	if (atomic_bool_CAS64(&table->count, table->max_cnt, table->max_cnt)){
		pthread_mutex_lock(&h_table_mut);
		try_expanding(h_table);
		pthread_mutex_unlock(&h_table_mut);
	}
	return mini_put_by_key(*h_table, data, key);
}

