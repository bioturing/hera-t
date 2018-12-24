#ifndef _DQUEUE_H_
#define _DQUEUE_H_

#include <stdint.h>
#include <pthread.h>

#include "semaphore_wrapper.h"

struct dqueue_t {
	int cap;
	void **in_data;
	void **out_data;
	int in_head;
	int in_tail;
	int out_head;
	int out_tail;
	struct sem_wrap_t sem_in;
	struct sem_wrap_t sem_out;
	pthread_mutex_t lock_in_head;
	pthread_mutex_t lock_in_tail;
	pthread_mutex_t lock_out_head;
	pthread_mutex_t lock_out_tail;
};

struct dqueue_t *init_dqueue(int cap);

void dqueue_destroy(struct dqueue_t *q);

void d_enqueue_in(struct dqueue_t *q, void *ptr);

void d_enqueue_out(struct dqueue_t *q, void *ptr);

void *d_dequeue_in(struct dqueue_t *q);

void *d_dequeue_out(struct dqueue_t *q);

#endif
