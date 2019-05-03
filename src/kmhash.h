#ifndef _KMHASH_H_
#define _KMHASH_H_

#include <pthread.h>
#include <stdint.h>

#include "atomic.h"
#include "pthread_barrier.h"
#include "semaphore_wrapper.h"

typedef uint32_t kmint_t;
#define __sync_fetch_and_add_kmint __sync_fetch_and_add32
#define __sync_val_compare_and_swap_kmint __sync_val_compare_and_swap32
#define __sync_bool_compare_and_swap_kmint __sync_bool_compare_and_swap32

#define KMHASH_MAX_SIZE			UINT32_C(0x80000000)
#define KMHASH_UMIHASH_SIZE		UINT32_C(0x10)
#define KMHASH_SIZE		UINT32_C(0x10000)
#define KMHASH_SINGLE_RESIZE		UINT32_C(0x100000)

#define KMHASH_UPPER			0.77

#define TOMB_STONE			((kmkey_t)-1)

typedef uint64_t kmkey_t;
#define __sync_fetch_and_add_kmkey __sync_fetch_and_add64
#define __sync_val_compare_and_swap_kmkey __sync_val_compare_and_swap64
#define __sync_bool_compare_and_swap_kmkey __sync_bool_compare_and_swap64
typedef uint64_t kmval_t;

struct umi_hash_t {
	kmint_t size;
	kmint_t n_items;
	kmkey_t *bucks;
};

struct kmbucket_t {
	kmkey_t idx;
	pthread_mutex_t lock;
	struct umi_hash_t *umis;
};

struct kmhash_t {
	kmint_t size;
	kmint_t old_size;
	kmint_t n_items;
	kmint_t n_probe;
	struct kmbucket_t *bucks;
	// struct kmbucket_t *old_bucks;
	uint8_t *flag;
	int status;
	int n_workers;
	// struct sem_wrap_t gsem;
	pthread_mutex_t *locks;
	int *pos;
};

struct kmresize_bundle_t {
	struct kmhash_t *h;
	int n_threads;
	int thread_no;
	pthread_barrier_t *barrier;
};

struct kmhash_t *init_kmhash(kmint_t size, int n_threads);

void kmhash_destroy(struct kmhash_t *h);

void kmhash_put_bc_umi(struct kmhash_t *h, pthread_mutex_t *lock,
						kmkey_t bc, kmkey_t umi);

void umihash_put_umi_single(struct umi_hash_t *h, kmkey_t key);

kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key);

kmint_t umihash_get(struct umi_hash_t *h, kmkey_t key);

void umihash_destroy(struct umi_hash_t *h);

#endif
