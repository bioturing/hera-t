#include <pthread.h>
#include <stdlib.h>

#include "kmhash.h"
#include "pthread_barrier.h"
#include "semaphore_wrapper.h"
#include "utils.h"
#include "verbose.h"
#include "atomic.h"

#define __hash_int(x) (uint64_t)((x) >> 33 ^ (x) ^ (x) << 11)

#define KMFLAG_EMPTY			0
#define KMFLAG_OLD			1
#define KMFLAG_NEW			2
#define KMFLAG_LOADING			3

#define __round_up_kmint(x) 	(--(x), (x) |= (x) >> 1,		       \
				 (x) |= (x) >> 2, (x) |= (x) >> 4,	       \
				 (x) |= (x) >> 8, (x) |= (x) >> 16,	       \
				 ++(x))

#define rs_set_old(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMFLAG_OLD << (((i) & 15) << 1)))

#define rs_set_new(x, i) ((x)[(i) >> 4] = ((x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1)))) |       \
				((uint32_t)KMFLAG_NEW << (((i) & 15) << 1)))

#define rs_set_empty(x, i) ((x)[(i) >> 4] = (x)[(i) >> 4] &		       \
				(~((uint32_t)3 << (((i) & 15) << 1))))

#define rs_get_flag(x, i) (((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)

#define rs_is_old(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)  \
							== (uint32_t)KMFLAG_OLD)

#define rs_is_empty(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)  \
							== (uint32_t)KMFLAG_EMPTY)

#define rs_is_new(x, i) ((((x)[(i) >> 4] >> (((i) & 15) << 1)) & (uint32_t)3)  \
							== (uint32_t)KMFLAG_NEW)

#define HM_MAGIC_1			UINT64_C(0xbf58476d1ce4e5b9)
#define HM_MAGIC_2			UINT64_C(0x94d049bb133111eb)

#define KMHASH_IDLE			0
#define KMHASH_BUSY			1

static inline uint64_t __hash_int2(kmkey_t k)
{
	kmkey_t x = k;
	x = (x ^ (x >> 30)) * HM_MAGIC_1;
	x = (x ^ (x >> 27)) * HM_MAGIC_2;
	x ^= (x > 31);
	return x;
}

static inline kmint_t estimate_probe_3(kmint_t size)
{
	kmint_t s, i;
	i = s = 0;
	while (s < size) {
		++i;
		s += i * i * i * 64;
	}
	return i;
}

static struct umi_hash_t *init_umi_hash()
{
	struct umi_hash_t *umis;
	umis = calloc(1, sizeof(struct umi_hash_t));
	umis->size = KMHASH_UMIHASH_SIZE;
	umis->bucks = malloc(umis->size * sizeof(kmkey_t));
	umis->n_items = 0;
	memset(umis->bucks, 255, sizeof(kmkey_t) * umis->size);
	return umis;
}

static kmint_t internal_umihash_put(struct umi_hash_t *h, kmkey_t key)
{
	if (h->n_items >= (kmint_t)(h->size * KMHASH_UPPER))
		return KMHASH_MAX_SIZE;

	kmint_t mask, i, last, step = 0;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_int(key);
	last = i = k & mask;
	if (h->bucks[i] == TOMB_STONE) {
		h->bucks[i] = key;
		++h->n_items;
		return i;
	}
	while (h->bucks[i] != TOMB_STONE && h->bucks[i] != key) {
		i = (i + (++step)) & mask;
		if (i == last)
			break;
	}
	if (h->bucks[i] == TOMB_STONE) {
		h->bucks[i] = key;
		++h->n_items;
		return i;
	}
	return h->bucks[i] == key ? i : KMHASH_MAX_SIZE;
}

// static kmint_t internal_umihash_put(struct umi_hash_t *h, kmkey_t key)
// {
// 	if (h->n_items >= (kmint_t)(h->size * KMHASH_UPPER))
// 		return KMHASH_MAX_SIZE;
// 	kmint_t mask, i, last, step = 0;
// 	kmkey_t cur_key;
// 	uint64_t k;

// 	mask = h->size - 1;
// 	k = __hash_int2(key);
// 	last = i = k & mask;
// 	{
// 		cur_key = __sync_val_compare_and_swap_kmkey(h->bucks + i, TOMB_STONE, key);
// 		if (cur_key == TOMB_STONE || cur_key == key) {
// 			if (cur_key == TOMB_STONE)
// 				__sync_fetch_and_add_kmint(&(h->n_items), 1);
// 			return i;
// 		}
// 	}
// 	do {
// 		i = (i + (++step)) & mask;
// 		cur_key = __sync_val_compare_and_swap_kmkey(h->bucks + i, TOMB_STONE, key);
// 	} while (cur_key != TOMB_STONE && cur_key != key && i != last);
// 	if (cur_key == TOMB_STONE || cur_key == key) {
// 		if (cur_key == TOMB_STONE)
// 			__sync_fetch_and_add_kmint(&(h->n_items), 1);
// 		return i;
// 	}
// 	return KMHASH_MAX_SIZE;
// }

// static kmint_t umihash_put_umi(struct umi_hash_t *h, kmkey_t key)
// {
// 	if ((kmint_t)(h->size * KMHASH_UPPER) < h->n_items)
// 		return KMHASH_MAX_SIZE;
// 	kmint_t mask, step, i, last;
// 	kmkey_t cur_key, k;

// 	mask = h->size - 1;
// 	k = __hash_int(key);
// 	last = i = k & mask;
// 	{
// 		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i]), TOMB_STONE, key);
// 		if (cur_key == TOMB_STONE || cur_key == key) {
// 			if (cur_key == TOMB_STONE)
// 				__sync_fetch_and_add_kmint(&(h->n_items), 1);
// 			return i;
// 		}
// 	}
// 	step = 0;
// 	do {
// 		i = (i + (++step)) & mask;
// 		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i]), TOMB_STONE, key);
// 	} while (cur_key != TOMB_STONE && cur_key != key && i != last);
// 	if (cur_key == TOMB_STONE || cur_key == key) {
// 		if (cur_key == TOMB_STONE)
// 			__sync_fetch_and_add_kmint(&(h->n_items), 1);
// 		return i;
// 	}
// 	return KMHASH_MAX_SIZE;
// }

// static kmint_t kmhash_put_bc_no_init(struct kmhash_t *h, kmkey_t key)
// {
// 	if ((kmint_t)(h->size * KMHASH_UPPER) < h->n_items)
// 		return KMHASH_MAX_SIZE;
// 	kmint_t mask, step, i, last;
// 	kmkey_t cur_key, k;

// 	mask = h->size - 1;
// 	k = __hash_int(key);
// 	last = i = k & mask;
// 	{
// 		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
// 		if (cur_key == TOMB_STONE || cur_key == key) {
// 			if (cur_key == TOMB_STONE) {
// 				// init_bc_bucks(h->bucks + i, h->n_workers);
// 				__sync_fetch_and_add_kmint(&(h->n_items), 1);
// 			}
// 			return i;
// 		}
// 	}
// 	step = 0;
// 	do {
// 		i = (i + (++step)) & mask;
// 		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
// 	} while (cur_key != TOMB_STONE && cur_key != key && i != last);
// 	if (cur_key == TOMB_STONE || cur_key == key) {
// 		if (cur_key == TOMB_STONE) {
// 			// init_bc_bucks(h->bucks + i, h->n_workers);
// 			__sync_fetch_and_add_kmint(&(h->n_items), 1);
// 		}
// 		return i;
// 	}
// 	return KMHASH_MAX_SIZE;
// }

// kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key)
// {
// 	kmint_t mask, step, i, last;
// 	kmkey_t k;

// 	mask = h->size - 1;
// 	k = __hash_int(key);
// 	last = i = k & mask;
// 	if (h->bucks[i].idx == key)
// 		return i;
// 	if (h->bucks[i].idx == TOMB_STONE)
// 		return KMHASH_MAX_SIZE;
// 	step = 0;
// 	do {
// 		i = (i + (++step)) & mask;
// 		if (h->bucks[i].idx == key)
// 			return i;
// 	} while (i != last && h->bucks[i].idx != TOMB_STONE);
// 	return KMHASH_MAX_SIZE;
// }

kmint_t kmhash_get(struct kmhash_t *h, kmkey_t key)
{
	kmint_t mask, i, step = 0;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_int2(key);
	i = k & mask;
	do {
		i = (i + step * (step + 1) / 2) & mask;
		if (h->bucks[i].idx == key)
			return i;
		++step;
	} while (step <= h->n_probe && h->bucks[i].idx != TOMB_STONE);
	return KMHASH_MAX_SIZE;
}

kmint_t umihash_get(struct umi_hash_t *h, kmkey_t key)
{
	kmint_t mask, i, last, step = 0;
	uint64_t k;
	mask = h->size - 1;
	k = __hash_int(key);
	last = i = k & mask;
	while (h->bucks[i] != TOMB_STONE && h->bucks[i] != key) {
		i = (i + (++step)) & mask;
		if (i == last)
			return KMHASH_MAX_SIZE;
	}
	return h->bucks[i] == key ? i : KMHASH_MAX_SIZE;
}

// kmint_t umihash_get(struct umi_hash_t *h, kmkey_t key)
// {
// 	kmint_t mask, step, i, last;
// 	kmkey_t k;

// 	mask = h->size - 1;
// 	k = __hash_int(key);
// 	last = i = k & mask;
// 	if (h->bucks[i] == key)
// 		return i;
// 	if (h->bucks[i] == TOMB_STONE)
// 		return KMHASH_MAX_SIZE;
// 	step = 0;
// 	do {
// 		i = (i + (++step)) & mask;
// 		if (h->bucks[i] == key)
// 			return i;
// 	} while (i != last && h->bucks[i] != TOMB_STONE);
// 	return KMHASH_MAX_SIZE;
// }

static kmint_t internal_kmhash_put(struct kmhash_t *h, kmkey_t key)
{
	kmint_t mask, i, step = 0;
	kmkey_t cur_key;
	uint64_t k;

	mask = h->size - 1;
	k = __hash_int2(key);
	i = k & mask;
	do {
		i = (i + step * (step + 1) / 2) & mask;
		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
		++step;
	} while (step <= h->n_probe && cur_key != key && cur_key != TOMB_STONE);
	if (cur_key == TOMB_STONE || cur_key == key) {
		if (cur_key == TOMB_STONE) {
			// init_bc_bucks(h->bucks + i, h->n_workers);
			__sync_fetch_and_add_kmint(&(h->n_items), 1);
		}
		return  i;
	}
	return KMHASH_MAX_SIZE;
}

// static kmint_t kmhash_put_bc(struct kmhash_t *h, kmkey_t key)
// {
// 	if ((kmint_t)(h->size * KMHASH_UPPER) < h->n_items)
// 		return KMHASH_MAX_SIZE;
// 	kmint_t mask, step, i, last;
// 	kmkey_t cur_key, k;

// 	mask = h->size - 1;
// 	k = __hash_int(key);
// 	last = i = k & mask;
// 	{
// 		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
// 		if (cur_key == TOMB_STONE || cur_key == key) {
// 			if (cur_key == TOMB_STONE) {
// 				init_bc_bucks(h->bucks + i, h->n_workers);
// 				__sync_fetch_and_add_kmint(&(h->n_items), 1);
// 			}
// 			return i;
// 		}
// 	}
// 	step = 0;
// 	do {
// 		i = (i + (++step)) & mask;
// 		cur_key = __sync_val_compare_and_swap_kmkey(&(h->bucks[i].idx), TOMB_STONE, key);
// 	} while (cur_key != TOMB_STONE && cur_key != key && i != last);
// 	if (cur_key == TOMB_STONE || cur_key == key) {
// 		if (cur_key == TOMB_STONE) {
// 			init_bc_bucks(h->bucks + i, h->n_workers);
// 			__sync_fetch_and_add_kmint(&(h->n_items), 1);
// 		}
// 		return i;
// 	}
// 	return KMHASH_MAX_SIZE;
// }

// void *kmresize_worker(void *data)
// {
// 	struct kmresize_bundle_t *bundle = (struct kmresize_bundle_t *)data;
// 	struct kmhash_t *h;
// 	int j;
// 	kmint_t i, k, l, r, cap;

// 	h = bundle->h;
// 	// Init all bucket
// 	cap = h->size / bundle->n_threads + 1;
// 	l = cap * bundle->thread_no;
// 	r = __min(cap * (bundle->thread_no + 1), h->size);
// 	for (i = l; i < r; ++i) {
// 		h->bucks[i].idx = TOMB_STONE;
// 		sem_wrap_init(&(h->bucks[i].bsem), 0);
// 	}

// 	pthread_barrier_wait(bundle->barrier);

// 	// Fill buckets
// 	cap = h->old_size / bundle->n_threads + 1;
// 	l = cap * bundle->thread_no;
// 	r = __min(cap * (bundle->thread_no + 1), h->old_size);
// 	for (i = l; i < r; ++i) {
// 		if (h->old_bucks[i].idx == TOMB_STONE)
// 			continue;
// 		k = kmhash_put_bc_no_init(h, h->old_bucks[i].idx);
// 		if (k == KMHASH_MAX_SIZE)
// 			__ERROR("Resizing barcodes hash table fail");
// 		h->bucks[k].umis = h->old_bucks[i].umis;
// 		for (j = 0; j < bundle->n_threads; ++j)
// 			sem_wrap_post(&(h->bucks[k].bsem));
// 		sem_wrap_destroy(&(h->old_bucks[i].bsem));
// 	}

// 	pthread_exit(NULL);
// }

void *kmresize_worker(void *data)
{
	struct kmresize_bundle_t *bundle = (struct kmresize_bundle_t *)data;
	struct kmhash_t *h;
	kmint_t i, l, r, mask, cap;

	h = bundle->h;
	mask = h->size - 1;
	/* Init all buckets */
	cap = (h->size - h->old_size) / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->size - h->old_size);
	for (i = l; i < r; ++i) {
		h->bucks[h->old_size + i].idx = TOMB_STONE;
		h->bucks[h->old_size + i].umis = NULL;
		pthread_mutex_init(&(h->bucks[h->old_size + i].lock), NULL);
	}

	cap = h->old_size / bundle->n_threads + 1;
	l = cap * bundle->thread_no;
	r = __min(cap * (bundle->thread_no + 1), h->old_size);
	for (i = l; i < r; ++i) {
		if (h->bucks[i].idx != TOMB_STONE)
			h->flag[i] = KMFLAG_OLD;
	}

	pthread_barrier_wait(bundle->barrier);

	for (i = l; i < r; ++i) {
		if (__sync_val_compare_and_swap8(h->flag + i, KMFLAG_OLD, KMFLAG_LOADING)
								== KMFLAG_OLD) {
			kmkey_t x = h->bucks[i].idx, xt;
			struct umi_hash_t *y = h->bucks[i].umis, *yt;
			h->bucks[i].idx = TOMB_STONE;
			h->bucks[i].umis = NULL;
			h->flag[i] = KMFLAG_EMPTY;
			while (1) {
				uint64_t k = __hash_int2(x);
				kmint_t j = k & mask, step = 0;
				uint8_t current_flag = KMFLAG_NEW;
				while (step <= h->n_probe) {
					if ((current_flag = __sync_val_compare_and_swap8(h->flag + j,
							KMFLAG_OLD, KMFLAG_NEW))
								== KMFLAG_OLD) {
						h->bucks[j].idx = x;
						h->bucks[j].umis = y;
						break;
					} else if ((current_flag == __sync_val_compare_and_swap8(h->flag +j,
							KMFLAG_EMPTY, KMFLAG_NEW))
								== KMFLAG_EMPTY) {
						xt = h->bucks[j].idx;
						yt = h->bucks[j].umis;
						h->bucks[j].idx = x;
						h->bucks[j].umis = y;
						x = xt;
						y = yt;
						break;
					} else if (current_flag == KMFLAG_LOADING) {
						/* not sure */
						continue;
					}
					++step;
					j = (j + step * (step + 1) / 2) & mask;
				}
				if (current_flag == KMFLAG_EMPTY)
					break;
				else if (current_flag == KMFLAG_NEW)
					__ERROR("[Multi-thread] resize kmhash failed");
			}
		}
	}
	pthread_exit(NULL);
}

void kmhash_resize_multi(struct kmhash_t *h)
{
	int n_threads, i;
	n_threads = h->n_workers;

	h->old_size = h->size;
	h->size <<= 1;
	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));
	h->flag = calloc(h->size, sizeof(uint8_t));

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
	free(h->flag);
	h->flag = NULL;
}

// void kmhash_resize_multi(struct kmhash_t *h)
// {
// 	int n_threads, i;
// 	n_threads = h->n_workers;

// 	h->old_size = h->size;
// 	h->old_bucks = h->bucks;

// 	h->size <<= 1;
// 	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));

// 	h->n_items = 0;

// 	pthread_attr_t attr;
// 	pthread_attr_init(&attr);
// 	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
// 	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

// 	pthread_t *t;
// 	t = calloc(h->n_workers, sizeof(pthread_t));

// 	pthread_barrier_t barrier;
// 	pthread_barrier_init(&barrier, NULL, n_threads);

// 	struct kmresize_bundle_t *bundles;
// 	bundles = calloc(n_threads, sizeof(struct kmresize_bundle_t));

// 	for (i = 0; i < n_threads; ++i) {
// 		bundles[i].n_threads = n_threads;
// 		bundles[i].thread_no = i;
// 		bundles[i].h = h;
// 		bundles[i].barrier = &barrier;
// 		pthread_create(t + i, &attr, kmresize_worker, bundles + i);
// 	}

// 	for (i = 0; i < n_threads; ++i)
// 		pthread_join(t[i], NULL);

// 	pthread_attr_destroy(&attr);
// 	pthread_barrier_destroy(&barrier);

// 	free(t);
// 	free(bundles);
// 	free(h->old_bucks);
// }

// void kmhash_resize_single(struct kmhash_t *h)
// {
// 	kmint_t i, k;
// 	int j;

// 	h->old_size = h->size;
// 	h->old_bucks = h->bucks;

// 	h->size <<= 1;
// 	h->bucks = malloc(h->size * sizeof(struct kmbucket_t));

// 	// Initilize new buckets
// 	h->n_items = 0;
// 	for (i = 0; i < h->size; ++i) {
// 		h->bucks[i].idx = TOMB_STONE;
// 		sem_wrap_init(&(h->bucks[i].bsem), 0);
// 	}

// 	for (i = 0; i < h->old_size; ++i) {
// 		if (h->old_bucks[i].idx == TOMB_STONE)
// 			continue;
// 		k = kmhash_put_bc_no_init(h, h->old_bucks[i].idx);
// 		if (k == KMHASH_MAX_SIZE)
// 			__ERROR("Resizing barcodes hash table fail");
// 		h->bucks[k].umis = h->old_bucks[i].umis;
// 		for (j = 0; j < h->n_workers; ++j)
// 			sem_wrap_post(&(h->bucks[k].bsem));
// 		sem_wrap_destroy(&(h->old_bucks[i].bsem));
// 	}
// 	free(h->old_bucks);
// }

void kmhash_resize_single(struct kmhash_t *h)
{
	kmint_t old_size, mask, i;
	old_size = h->size;
	h->size <<= 1;
	mask = h->size - 1;
	h->n_probe = estimate_probe_3(h->size);

	h->bucks = realloc(h->bucks, h->size * sizeof(struct kmbucket_t));
	uint8_t *flag = calloc(h->size, sizeof(uint8_t));

	for (i = old_size; i < h->size; ++i) {
		h->bucks[i].idx = TOMB_STONE;
		h->bucks[i].umis = NULL;
		pthread_mutex_init(&(h->bucks[i].lock), NULL);
	}

	for (i = 0; i < old_size; ++i)
		if (h->bucks[i].idx != TOMB_STONE)
			flag[i] = KMFLAG_OLD;

	for (i = 0; i < old_size; ++i) {
		if (flag[i] == KMFLAG_OLD) {
			kmkey_t x = h->bucks[i].idx, xt;
			struct umi_hash_t *y = h->bucks[i].umis, *yt;
			h->bucks[i].idx = TOMB_STONE;
			h->bucks[i].umis = NULL;
			flag[i] = KMFLAG_EMPTY;
			while (1) { /* kick-out process */
				uint64_t k = __hash_int2(x);
				kmint_t j = k & mask, step = 0;
				uint8_t current_flag = KMFLAG_NEW;
				while (step <= h->n_probe) {
					j = (j + step * (step + 1) / 2) & mask;
					current_flag = flag[j];
					if (current_flag == KMFLAG_EMPTY) {
						flag[j] = KMFLAG_NEW;
						h->bucks[j].idx = x;
						h->bucks[j].umis = y;
						break;
					} else if (current_flag == KMFLAG_OLD) {
						flag[j] = KMFLAG_NEW;
						xt = h->bucks[j].idx;
						yt = h->bucks[j].umis;
						h->bucks[j].idx = x;
						h->bucks[j].umis = y;
						x = xt;
						y = yt;
						break;
					}
					++step;
				}
				if (current_flag == KMFLAG_EMPTY)
					break;
				else if (current_flag == KMFLAG_NEW)
					__ERROR("[Single-thread] Resizing kmhash failed");
			}
		}
	}
	free(flag);
}

void kmhash_resize(struct kmhash_t *h)
{
	int i;
	for (i = 0; i < h->n_workers; ++i)
		pthread_mutex_lock(h->locks + i);

	if (h->size == KMHASH_MAX_SIZE)
		__ERROR("The barcodes hash table is too big (exceeded %llu)",
			(unsigned long long)KMHASH_MAX_SIZE);

	if (h->size <= KMHASH_SINGLE_RESIZE)
		kmhash_resize_single(h);
	else
		kmhash_resize_multi(h);

	for (i = 0; i < h->n_workers; ++i)
		pthread_mutex_unlock(h->locks + i);
}

// void kmhash_resize(struct kmhash_t *h)
// {
// 	int i;
// 	for (i = 0; i < h->n_workers; ++i)
// 		sem_wrap_wait(&(h->gsem));

// 	if (h->size == KMHASH_MAX_SIZE)
// 		__ERROR("The barcodes hash table is too big (exceeded %llu)",
// 			(unsigned long long)KMHASH_MAX_SIZE);

// 	if (h->size <= KMHASH_SINGLE_RESIZE)
// 		kmhash_resize_single(h);
// 	else
// 		kmhash_resize_multi(h);

// 	for (i = 0; i < h->n_workers; ++i)
// 		sem_wrap_post(&(h->gsem));
// }

// void *umiresize_worker(void *data)
// {
// 	struct umiresize_bundle_t *bundle = (struct umiresize_bundle_t *)data;
// 	struct umi_hash_t *h;
// 	kmint_t i, k, l, r, cap;

// 	h = bundle->h;

// 	// Init new buckets
// 	cap = h->size / bundle->n_threads + 1;
// 	l = cap * bundle->thread_no;
// 	r = __min(cap * (bundle->thread_no + 1), h->size);
// 	for (i = l; i < r; ++i)
// 		h->bucks[i] = TOMB_STONE;

// 	pthread_barrier_wait(bundle->barrier);

// 	// Fill new buckets
// 	cap = h->old_size / bundle->n_threads + 1;
// 	l = cap * bundle->thread_no;
// 	r = __min(cap * (bundle->thread_no + 1), h->old_size);
// 	for (i = l; i < r; ++i) {
// 		if (h->old_bucks[i] == TOMB_STONE)
// 			continue;
// 		k = umihash_put_umi(h, h->old_bucks[i]);
// 		if (k == KMHASH_MAX_SIZE)
// 			__ERROR("[Multi] Resizing UMIs hash table fail");
// 	}

// 	pthread_exit(NULL);
// }

// void umihash_resize_multi(struct umi_hash_t *h, int n_threads)
// {
// 	int i;
// 	h->old_size = h->size;
// 	h->old_bucks = h->bucks;

// 	h->size <<= 1;
// 	h->bucks = malloc(h->size * sizeof(kmkey_t));

// 	h->n_items = 0;

// 	pthread_attr_t attr;
// 	pthread_attr_init(&attr);
// 	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
// 	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

// 	pthread_t *t;
// 	t = calloc(n_threads, sizeof(pthread_t));

// 	pthread_barrier_t barrier;
// 	pthread_barrier_init(&barrier, NULL, n_threads);

// 	struct umiresize_bundle_t *bundles;
// 	bundles = calloc(n_threads, sizeof(struct umiresize_bundle_t));

// 	for (i = 0; i < n_threads; ++i) {
// 		bundles[i].n_threads = n_threads;
// 		bundles[i].thread_no = i;
// 		bundles[i].h = h;
// 		bundles[i].barrier = &barrier;
// 		pthread_create(t + i, &attr, umiresize_worker, bundles + i);
// 	}

// 	for (i = 0; i < n_threads; ++i)
// 		pthread_join(t[i], NULL);

// 	pthread_attr_destroy(&attr);
// 	pthread_barrier_destroy(&barrier);

// 	free(t);
// 	free(bundles);
// 	free(h->old_bucks);
// }

// void umihash_resize_single(struct umi_hash_t *h)
// {
// 	kmint_t i, k;

// 	h->old_size = h->size;
// 	h->old_bucks = h->bucks;

// 	h->size <<= 1;
// 	h->bucks = malloc(h->size * sizeof(kmkey_t));

// 	// Initilize new buckets
// 	h->n_items = 0;
// 	for (i = 0; i < h->size; ++i)
// 		h->bucks[i] = TOMB_STONE;

// 	for (i = 0; i < h->old_size; ++i) {
// 		if (h->old_bucks[i] == TOMB_STONE)
// 			continue;
// 		k = umihash_put_umi(h, h->old_bucks[i]);
// 		if (k == KMHASH_MAX_SIZE)
// 			__ERROR("[Single] Resizing UMIs hash table fail");
// 	}
// 	free(h->old_bucks);
// }

// void umihash_resize(struct kmbucket_t *b, int n_threads)
// {
// 	struct umi_hash_t *h;
// 	int i;

// 	h = b->umis;
// 	for (i = 0; i < n_threads; ++i)
// 		sem_wrap_wait(&(b->bsem));

// 	if (h->size == KMHASH_MAX_SIZE)
// 		__ERROR("The UMIs hash table is too big (exceeded %llu)",
// 			(unsigned long long)KMHASH_MAX_SIZE);

// 	if (h->size <= KMHASH_SINGLE_RESIZE)
// 		umihash_resize_single(h);
// 	else
// 		umihash_resize_multi(h, n_threads);

// 	for (i = 0; i < n_threads; ++i)
// 		sem_wrap_post(&(b->bsem));
// }

static void umihash_resize(struct umi_hash_t *h)
{
	kmint_t old_size, mask, i;
	old_size = h->size;
	h->size <<= 1;
	mask = h->size - 1;
	h->bucks = realloc(h->bucks, sizeof(kmkey_t) * h->size);
	uint32_t *flag = calloc(h->size >> 4, sizeof(uint32_t));

	for (i = old_size; i < h->size; ++i)
		h->bucks[i] = TOMB_STONE;

	for (i = 0; i < old_size; ++i) {
		if (h->bucks[i] != TOMB_STONE)
			rs_set_old(flag, i);
	}

	for (i = 0; i < old_size; ++i) {
		if (rs_is_old(flag, i)) {
			kmkey_t x = h->bucks[i], xt;
			h->bucks[i] = TOMB_STONE;
			rs_set_empty(flag, i);
			while (1) {
				uint64_t k = __hash_int(x);
				kmint_t j, last, step = 0;
				j = last = k & mask;
				while (!rs_is_empty(flag, j) && !rs_is_old(flag, j)) {
					j = (j + (++step)) & mask;
					if (j == last)
						break;
				}
				if (rs_is_empty(flag, j)) {
					rs_set_new(flag, j);
					h->bucks[j] = x;
					break;
				} else if (rs_is_old(flag, j)) {
					rs_set_new(flag, j);
					xt = h->bucks[j];
					h->bucks[j] = x;
					x = xt;
				} else {
					__ERROR("Resize UMI hash failed");
				}
			}
		}
	}
	free(flag);
}

static void internal_kmhash_put_umi(struct kmbucket_t *b, kmkey_t umi)
{
	kmint_t k;
	pthread_mutex_lock(&(b->lock));
	if (b->umis == NULL)
		b->umis = init_umi_hash();
	k = internal_umihash_put(b->umis, umi);
	while (k == KMHASH_MAX_SIZE) {
		umihash_resize(b->umis);
		k = internal_umihash_put(b->umis, umi);
	}
	pthread_mutex_unlock(&(b->lock));
}

// static void kmhash_put_umi(struct kmbucket_t *b, kmkey_t umi, int n_threads)
// {
// 	struct umi_hash_t *umis;
// 	kmint_t k;
// 	umis = b->umis;

// 	sem_wrap_wait(&(b->bsem));
// 	k = umihash_put_umi(umis, umi);
// 	sem_wrap_post(&(b->bsem));

// 	if (k == KMHASH_MAX_SIZE) {
// 		do {
// 			if (__sync_bool_compare_and_swap32(&(b->umis->status), KMHASH_IDLE, KMHASH_BUSY)) {
// 				umihash_resize(b, n_threads);
// 				__sync_val_compare_and_swap32(&(b->umis->status), KMHASH_BUSY, KMHASH_IDLE);
// 			}

// 			sem_wrap_wait(&(b->bsem));
// 			k = umihash_put_umi(umis, umi);
// 			sem_wrap_post(&(b->bsem));
// 		} while (k == KMHASH_MAX_SIZE);
// 	}
// }

void umihash_put_umi_single(struct umi_hash_t *h, kmkey_t key)
{
	kmint_t k;
	k = internal_umihash_put(h, key);
	while (k == KMHASH_MAX_SIZE) {
		umihash_resize(h);
		k = internal_umihash_put(h, key);
	}
}

// void umihash_put_umi_single(struct umi_hash_t *h, kmkey_t key)
// {
// 	kmint_t k;
// 	k = umihash_put_umi(h, key);
// 	if (k == KMHASH_MAX_SIZE) {
// 		umihash_resize_single(h);
// 		k = umihash_put_umi(h, key);
// 	}
// 	if (k == KMHASH_MAX_SIZE)
// 		__ERROR("Fatal error: unable to resize umi hash");
// }

void kmhash_put_bc_umi(struct kmhash_t *h, pthread_mutex_t *lock,
						kmkey_t bc, kmkey_t umi)
{
	kmint_t k;
	pthread_mutex_lock(lock);
	k = internal_kmhash_put(h, bc);
	if (k != KMHASH_MAX_SIZE)
		internal_kmhash_put_umi(h->bucks + k, umi);
	pthread_mutex_unlock(lock);

	while (k == KMHASH_MAX_SIZE) {
		if (__sync_bool_compare_and_swap32(&(h->status), KMHASH_IDLE, KMHASH_BUSY)) {
			kmhash_resize(h);
			__sync_val_compare_and_swap32(&(h->status), KMHASH_BUSY, KMHASH_IDLE);
		}
		pthread_mutex_lock(lock);
		k = internal_kmhash_put(h, bc);
		if (k != KMHASH_MAX_SIZE)
			internal_kmhash_put_umi(h->bucks + k, umi);
		pthread_mutex_unlock(lock);
	}
}

// void kmhash_put_bc_umi(struct kmhash_t *h, kmkey_t bc, kmkey_t umi)
// {
// 	kmint_t k;

// 	sem_wrap_wait(&(h->gsem));
// 	k = kmhash_put_bc(h, bc);
// 	if (k < KMHASH_MAX_SIZE)
// 		kmhash_put_umi(h->bucks + k, umi, h->n_workers);
// 	sem_wrap_post(&(h->gsem));

// 	if (k == KMHASH_MAX_SIZE) {
// 		do {
// 			if (__sync_bool_compare_and_swap32(&(h->status), KMHASH_IDLE, KMHASH_BUSY)) {
// 				kmhash_resize(h);
// 				__sync_val_compare_and_swap32(&(h->status), KMHASH_BUSY, KMHASH_IDLE);
// 			}

// 			sem_wrap_wait(&(h->gsem));
// 			k = kmhash_put_bc(h, bc);
// 			if (k < KMHASH_MAX_SIZE)
// 				kmhash_put_umi(h->bucks + k, umi, h->n_workers);
// 			sem_wrap_post(&(h->gsem));
// 		} while (k == KMHASH_MAX_SIZE);
// 	}
// }

struct kmhash_t *init_kmhash(kmint_t size, int n_threads)
{
	struct kmhash_t *h;
	kmint_t i;

	h = calloc(1, sizeof(struct kmhash_t));
	h->size = size;
	__round_up_kmint(h->size);
	h->bucks = calloc(h->size, sizeof(struct kmbucket_t));
	for (i = 0; i < h->size; ++i) {
		h->bucks[i].idx = TOMB_STONE;
		h->bucks[i].umis = NULL;
		pthread_mutex_init(&(h->bucks[i].lock), NULL);
	}
	h->n_workers = n_threads;
	h->locks = calloc(h->n_workers, sizeof(pthread_mutex_t));
	int k;
	for (k = 0; k < h->n_workers; ++k)
		pthread_mutex_init(h->locks + k, NULL);
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
	for (i = 0; i < h->size; ++i) {
		pthread_mutex_destroy(&(h->bucks[i].lock));
		umihash_destroy(h->bucks[i].umis);
	}
	free(h->bucks);
	int k;
	for (k = 0; k < h->n_workers; ++k)
		pthread_mutex_detroy(h->locks + k);
	free(h->locks);
	free(h->pos);
	free(h);
}

