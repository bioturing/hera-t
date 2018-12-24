#include <pthread.h>
#include <stdlib.h>

#include "kmhash.h"
#include "pthread_barrier.h"
#include "semaphore_wrapper.h"
#include "utils.h"
#include "verbose.h"
#include "atomic.h"

#define __hash_int(x) (kmint_t)((x) >> 33 ^ (x) ^ (x) << 11)

#define __round_up_kmint(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 ++(x))

#define HM_MAGIC_1			UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MAGIC_2			UINT64_C(0x94d049bb133111eb)

#define KMHASH_IDLE			0
#define KMHASH_BUSY			1

static kmkey_t __hash_int2(kmkey_t k)
{
	kmkey_t x = k;
	x = (x ^ (x >> 30)) * HM_MAGIC_1;
	x = (x ^ (x >> 27)) * HM_MAGIC_2;
	x ^= (x > 31);
	return x;
}

static void init_bc_bucks(struct kmbucket_t *b, int n_threads)
{
	int i;
	kmint_t k;
	struct umi_hash_t *umis;
	b->umis = calloc(1, sizeof(struct umi_hash_t));
	umis = b->umis;
	umis->size = KMHASH_UMIHASH_SIZE;
	umis->bucks = malloc(umis->size * sizeof(kmkey_t));
	umis->status = KMHASH_IDLE;
	for (k = 0; k < umis->size; ++k)
		umis->bucks[k] = TOMB_STONE;
	for (i = 0; i < n_threads; ++i)
		sem_wrap_post(&(b->bsem));
}

static kmint_t umihash_put_umi(struct umi_hash_t *h, kmkey_t key)
{
	if ((kmint_t)(h->size * KMHASH_UPPER) < h->n_items)
		return KMHASH_MAX_SIZE;
	kmint_t mask, step, i, last;
	kmkey_t cur_key, k;

	mask = h->size - 1;
	k = __hash_int(key);
	last = i = k & mask;
	{
		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i]), TOMB_STONE, key);
		if (cur_key == TOMB_STONE || cur_key == key) {
			if (cur_key == TOMB_STONE)
				__sync_fetch_and_add_kmint(&(h->n_items), 1);
			return i;
		}
	}
	step = 0;
	do {
		i = (i + (++step)) & mask;
		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i]), TOMB_STONE, key);
	} while (cur_key != TOMB_STONE && cur_key != key && i != last);
	if (cur_key == TOMB_STONE || cur_key == key) {
		if (cur_key == TOMB_STONE)
			__sync_fetch_and_add_kmint(&(h->n_items), 1);
		return i;
	}
	return KMHASH_MAX_SIZE;
}

static kmint_t kmhash_put_bc_no_init(struct kmhash_t *h, kmkey_t key)
{
	if ((kmint_t)(h->size * KMHASH_UPPER) < h->n_items)
		return KMHASH_MAX_SIZE;
	kmint_t mask, step, i, last;
	kmkey_t cur_key, k;

	mask = h->size - 1;
	k = __hash_int(key);
	last = i = k & mask;
	{
		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
		if (cur_key == TOMB_STONE || cur_key == key) {
			if (cur_key == TOMB_STONE) {
				// init_bc_bucks(h->bucks + i, h->n_workers);
				__sync_fetch_and_add_kmint(&(h->n_items), 1);
			}
			return i;
		}
	}
	step = 0;
	do {
		i = (i + (++step)) & mask;
		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
	} while (cur_key != TOMB_STONE && cur_key != key && i != last);
	if (cur_key == TOMB_STONE || cur_key == key) {
		if (cur_key == TOMB_STONE) {
			// init_bc_bucks(h->bucks + i, h->n_workers);
			__sync_fetch_and_add_kmint(&(h->n_items), 1);
		}
		return i;
	}
	return KMHASH_MAX_SIZE;
}

kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key)
{
	kmint_t mask, step, i, last;
	kmkey_t k;

	mask = h->size - 1;
	k = __hash_int(key);
	last = i = k & mask;
	if (h->bucks[i].idx == key)
		return i;
	if (h->bucks[i].idx == TOMB_STONE)
		return KMHASH_MAX_SIZE;
	step = 0;
	do {
		i = (i + (++step)) & mask;
		if (h->bucks[i].idx == key)
			return i;
	} while (i != last && h->bucks[i].idx != TOMB_STONE);
	return KMHASH_MAX_SIZE;
}

kmint_t umihash_get(struct umi_hash_t *h, kmkey_t key)
{
	kmint_t mask, step, i, last;
	kmkey_t k;

	mask = h->size - 1;
	k = __hash_int(key);
	last = i = k & mask;
	if (h->bucks[i] == key)
		return i;
	if (h->bucks[i] == TOMB_STONE)
		return KMHASH_MAX_SIZE;
	step = 0;
	do {
		i = (i + (++step)) & mask;
		if (h->bucks[i] == key)
			return i;
	} while (i != last && h->bucks[i] != TOMB_STONE);
	return KMHASH_MAX_SIZE;
}

static kmint_t kmhash_put_bc(struct kmhash_t *h, kmkey_t key)
{
	if ((kmint_t)(h->size * KMHASH_UPPER) < h->n_items)
		return KMHASH_MAX_SIZE;
	kmint_t mask, step, i, last;
	kmkey_t cur_key, k;

	mask = h->size - 1;
	k = __hash_int(key);
	last = i = k & mask;
	{
		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
		if (cur_key == TOMB_STONE || cur_key == key) {
			if (cur_key == TOMB_STONE) {
				init_bc_bucks(h->bucks + i, h->n_workers);
				__sync_fetch_and_add_kmint(&(h->n_items), 1);
			}
			return i;
		}
	}
	step = 0;
	do {
		i = (i + (++step)) & mask;
		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
	} while (cur_key != TOMB_STONE && cur_key != key && i != last);
	if (cur_key == TOMB_STONE || cur_key == key) {
		if (cur_key == TOMB_STONE) {
			init_bc_bucks(h->bucks + i, h->n_workers);
			__sync_fetch_and_add_kmint(&(h->n_items), 1);
		}
		return i;
	}
	return KMHASH_MAX_SIZE;
}

void *kmresize_worker(void *data)
{
	struct kmresize_bundle_t *bundle = (struct kmresize_bundle_t *)data;
	struct kmhash_t *h;
	int j;
	kmint_t i, k, l, r, cap;

	h = bundle->h;
	// Init all bucket
	cap = h->size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->size);
	for (i = l; i < r; ++i) {
		h->bucks[i].idx = TOMB_STONE;
		sem_wrap_init(&(h->bucks[i].bsem), 0);
	}

	pthread_barrier_wait(bundle->barrier);

	// Fill buckets
	cap = h->old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->old_size);
	for (i = l; i < r; ++i) {
		if (h->old_bucks[i].idx == TOMB_STONE)
			continue;
		k = kmhash_put_bc_no_init(h, h->old_bucks[i].idx);
		if (k == KMHASH_MAX_SIZE)
			__ERROR("Resizing barcodes hash table fail");
		h->bucks[k].umis = h->old_bucks[i].umis;
		for (j = 0; j < bundle->n_threads; ++j)
			sem_wrap_post(&(h->bucks[k].bsem));
		sem_wrap_destroy(&(h->old_bucks[i].bsem));
	}

	pthread_exit(NULL);
}

void kmhash_resize_multi(struct kmhash_t *h)
{
	int n_threads, i;
	n_threads = h->n_workers;

	h->old_size = h->size;
	h->old_bucks = h->bucks;

	h->size <<= 1;
	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));

	h->n_items = 0;

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *t;
	t = calloc(h->n_workers, sizeof(pthread_t));

	pthread_barrier_t barrier;
	pthread_barrier_init(&barrier, NULL, n_threads);

	struct kmresize_bundle_t *bundles;
	bundles = calloc(n_threads, sizeof(struct kmresize_bundle_t));

	for (i = 0; i < n_threads; ++i) {
		bundles[i].n_threads = n_threads;
		bundles[i].thread_no = i;
		bundles[i].h = h;
		bundles[i].barrier = &barrier;
		pthread_create(t + i, &attr, kmresize_worker, bundles + i);
	}

	for (i = 0; i < n_threads; ++i)
		pthread_join(t[i], NULL);

	pthread_attr_destroy(&attr);
	pthread_barrier_destroy(&barrier);

	free(t);
	free(bundles);
	free(h->old_bucks);
}

void kmhash_resize_single(struct kmhash_t *h)
{
	kmint_t i, k;
	int j;

	h->old_size = h->size;
	h->old_bucks = h->bucks;

	h->size <<= 1;
	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));

	// Initilize new buckets
	h->n_items = 0;
	for (i = 0; i < h->size; ++i) {
		h->bucks[i].idx = TOMB_STONE;
		sem_wrap_init(&(h->bucks[i].bsem), 0);
	}

	for (i = 0; i < h->old_size; ++i) {
		if (h->old_bucks[i].idx == TOMB_STONE)
			continue;
		k = kmhash_put_bc_no_init(h, h->old_bucks[i].idx);
		if (k == KMHASH_MAX_SIZE)
			__ERROR("Resizing barcodes hash table fail");
		h->bucks[k].umis = h->old_bucks[i].umis;
		for (j = 0; j < h->n_workers; ++j)
			sem_wrap_post(&(h->bucks[k].bsem));
		sem_wrap_destroy(&(h->old_bucks[i].bsem));
	}
	free(h->old_bucks);
}

void kmhash_resize(struct kmhash_t *h)
{
	int i;
	for (i = 0; i < h->n_workers; ++i)
		sem_wrap_wait(&(h->gsem));

	if (h->size == KMHASH_MAX_SIZE)
		__ERROR("The barcodes hash table is too big (exceeded %llu)",
			(unsigned long long)KMHASH_MAX_SIZE);

	if (h->size <= KMHASH_SINGLE_RESIZE)
		kmhash_resize_single(h);
	else
		kmhash_resize_multi(h);

	for (i = 0; i < h->n_workers; ++i)
		sem_wrap_post(&(h->gsem));
}

void *umiresize_worker(void *data)
{
	struct umiresize_bundle_t *bundle = (struct umiresize_bundle_t *)data;
	struct umi_hash_t *h;
	kmint_t i, k, l, r, cap;

	h = bundle->h;

	// Init new buckets
	cap = h->size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->size);
	for (i = l; i < r; ++i)
		h->bucks[i] = TOMB_STONE;

	pthread_barrier_wait(bundle->barrier);

	// Fill new buckets
	cap = h->old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->old_size);
	for (i = l; i < r; ++i) {
		if (h->old_bucks[i] == TOMB_STONE)
			continue;
		k = umihash_put_umi(h, h->old_bucks[i]);
		if (k == KMHASH_MAX_SIZE)
			__ERROR("[Multi] Resizing UMIs hash table fail");
	}

	pthread_exit(NULL);
}

void umihash_resize_multi(struct umi_hash_t *h, int n_threads)
{
	int i;
	h->old_size = h->size;
	h->old_bucks = h->bucks;

	h->size <<= 1;
	h->bucks = malloc(h->size * sizeof(kmkey_t));

	h->n_items = 0;

	pthread_attr_t attr;
	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	pthread_t *t;
	t = calloc(n_threads, sizeof(pthread_t));

	pthread_barrier_t barrier;
	pthread_barrier_init(&barrier, NULL, n_threads);

	struct umiresize_bundle_t *bundles;
	bundles = calloc(n_threads, sizeof(struct umiresize_bundle_t));

	for (i = 0; i < n_threads; ++i) {
		bundles[i].n_threads = n_threads;
		bundles[i].thread_no = i;
		bundles[i].h = h;
		bundles[i].barrier = &barrier;
		pthread_create(t + i, &attr, umiresize_worker, bundles + i);
	}

	for (i = 0; i < n_threads; ++i)
		pthread_join(t[i], NULL);

	pthread_attr_destroy(&attr);
	pthread_barrier_destroy(&barrier);

	free(t);
	free(bundles);
	free(h->old_bucks);
}

void umihash_resize_single(struct umi_hash_t *h)
{
	kmint_t i, k;

	h->old_size = h->size;
	h->old_bucks = h->bucks;

	h->size <<= 1;
	h->bucks = malloc(h->size * sizeof(kmkey_t));

	// Initilize new buckets
	h->n_items = 0;
	for (i = 0; i < h->size; ++i)
		h->bucks[i] = TOMB_STONE;

	for (i = 0; i < h->old_size; ++i) {
		if (h->old_bucks[i] == TOMB_STONE)
			continue;
		k = umihash_put_umi(h, h->old_bucks[i]);
		if (k == KMHASH_MAX_SIZE)
			__ERROR("[Single] Resizing UMIs hash table fail");
	}
	free(h->old_bucks);
}

void umihash_resize(struct kmbucket_t *b, int n_threads)
{
	struct umi_hash_t *h;
	int i;

	h = b->umis;
	for (i = 0; i < n_threads; ++i)
		sem_wrap_wait(&(b->bsem));

	if (h->size == KMHASH_MAX_SIZE)
		__ERROR("The UMIs hash table is too big (exceeded %llu)",
			(unsigned long long)KMHASH_MAX_SIZE);

	if (h->size <= KMHASH_SINGLE_RESIZE)
		umihash_resize_single(h);
	else
		umihash_resize_multi(h, n_threads);

	for (i = 0; i < n_threads; ++i)
		sem_wrap_post(&(b->bsem));
}

static void kmhash_put_umi(struct kmbucket_t *b, kmkey_t umi, int n_threads)
{
	struct umi_hash_t *umis;
	kmint_t k;
	umis = b->umis;

	sem_wrap_wait(&(b->bsem));
	k = umihash_put_umi(umis, umi);
	sem_wrap_post(&(b->bsem));

	if (k == KMHASH_MAX_SIZE) {
		do {
			if (__sync_bool_compare_and_swap32(&(b->umis->status), KMHASH_IDLE, KMHASH_BUSY)) {
				umihash_resize(b, n_threads);
				__sync_val_compare_and_swap32(&(b->umis->status), KMHASH_BUSY, KMHASH_IDLE);
			}

			sem_wrap_wait(&(b->bsem));
			k = umihash_put_umi(umis, umi);
			sem_wrap_post(&(b->bsem));
		} while (k == KMHASH_MAX_SIZE);
	}
}

void umihash_put_umi_single(struct umi_hash_t *h, kmkey_t key)
{
	kmint_t k;
	k = umihash_put_umi(h, key);
	if (k == KMHASH_MAX_SIZE) {
		umihash_resize_single(h);
		k = umihash_put_umi(h, key);
	}
	if (k == KMHASH_MAX_SIZE)
		__ERROR("Fatal error: unable to resize umi hash");
}

void kmhash_put_bc_umi(struct kmhash_t *h, kmkey_t bc, kmkey_t umi)
{
	kmint_t k;

	sem_wrap_wait(&(h->gsem));
	k = kmhash_put_bc(h, bc);
	if (k < KMHASH_MAX_SIZE)
		kmhash_put_umi(h->bucks + k, umi, h->n_workers);
	sem_wrap_post(&(h->gsem));

	if (k == KMHASH_MAX_SIZE) {
		do {
			if (__sync_bool_compare_and_swap32(&(h->status), KMHASH_IDLE, KMHASH_BUSY)) {
				kmhash_resize(h);
				__sync_val_compare_and_swap32(&(h->status), KMHASH_BUSY, KMHASH_IDLE);
			}

			sem_wrap_wait(&(h->gsem));
			k = kmhash_put_bc(h, bc);
			if (k < KMHASH_MAX_SIZE)
				kmhash_put_umi(h->bucks + k, umi, h->n_workers);
			sem_wrap_post(&(h->gsem));
		} while (k == KMHASH_MAX_SIZE);
	}
}

// kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key)
// {
// 	kmint_t mask, step, i, n_probe;
// 	kmkey_t k;
// 	mask = h->size - 1;
// 	k = __hash_int2(key);
// 	i = k & mask;
// 	n_probe = h->n_probe;
// 	if (h->bucks[i].idx == key)
// 		return i;
// 	step = 1;
// 	do {
// 		i = (i + (step * (step + 1)) / 2) & mask;
// 		if (h->bucks[i].idx == key)
// 			return i;
// 		++step;
// 	} while (step < n_probe);
// 	return h->size;
// }

// struct kmhash_t *init_kmhash(kmint_t size, int n_threads)
// {
//         __VERBOSE("Initilizing hash table\n");
// 	struct kmhash_t *h;
// 	kmint_t i;
// 	kmkey_t tombstone;
// 	h = calloc(1, sizeof(struct kmhash_t));
// 	h->size = size;
// 	__round_up_kmint(h->size);
// 	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));
// 	h->n_probe = estimate_probe(h->size);
//         __VERBOSE("Probe number: %d\n", (int)h->n_probe);

// 	tombstone = (kmkey_t)-1;
// 	for (i = 0; i < h->size; ++i)
// 		h->bucks[i].idx = tombstone;

// 	h->n_workers = n_threads;
// 	h->locks = calloc(n_threads, sizeof(pthread_mutex_t));
// 	int k;
// 	for (k = 0; k < n_threads; ++k)
// 		pthread_mutex_init(h->locks + k, NULL);
// 	h->status = KMHASH_IDLE;

// 	return h;
// }

struct kmhash_t *init_kmhash(kmint_t size, int n_threads)
{
	struct kmhash_t *h;
	kmint_t i;

	h = calloc(1, sizeof(struct kmhash_t));
	h->size = size;
	__round_up_kmint(h->size);
	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));
	for (i = 0; i < h->size; ++i) {
		h->bucks[i].idx = TOMB_STONE;
		sem_wrap_init(&(h->bucks[i].bsem), 0);
	}

	h->n_workers = n_threads;
	sem_wrap_init(&(h->gsem), n_threads);
	h->status = KMHASH_IDLE;

	return h;
}

void umihash_destroy(struct umi_hash_t *h)
{
	if (!h) return;

	free(h->bucks);
	free(h);
}

void kmhash_destroy(struct kmhash_t *h)
{
	if (!h) return;
	kmint_t i;

	sem_wrap_destroy(&(h->gsem));
	for (i = 0; i < h->size; ++i) {
		if (h->bucks[i].idx == TOMB_STONE)
			continue;
		sem_wrap_destroy(&(h->bucks[i].bsem));
		umihash_destroy(h->bucks[i].umis);
	}
	free(h->bucks);
	free(h->pos);
	free(h);
}

// void kmhash_destroy(struct kmhash_t *h)
// {
// 	if (!h) return;
// 	free(h->bucks);
// 	int i;
// 	for (i = 0; i < h->n_workers; ++i)
// 		pthread_mutex_destroy(h->locks + i);
// 	free(h->locks);
// 	free(h);
// }
