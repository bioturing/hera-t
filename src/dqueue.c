#include "dqueue.h"
#include <stdlib.h>

struct dqueue_t *init_dqueue(int cap)
{
	struct dqueue_t *ret = malloc(sizeof(struct dqueue_t));
	ret->in_head = ret->out_head = -1;
	ret->in_tail = ret->out_tail = 0;
	ret->cap = cap;
	ret->in_data = malloc(cap * sizeof(void *));
	ret->out_data = malloc(cap * sizeof(void *));
	sem_wrap_init(&ret->sem_in, 0);
	sem_wrap_init(&ret->sem_out, 0);
	pthread_mutex_init(&ret->lock_in_head, NULL);
	pthread_mutex_init(&ret->lock_in_tail, NULL);
	pthread_mutex_init(&ret->lock_out_head, NULL);
	pthread_mutex_init(&ret->lock_out_tail, NULL);
	return ret;
}

void dqueue_destroy(struct dqueue_t *q)
{
	if (!q) return;
	free(q->in_data);
	free(q->out_data);
	pthread_mutex_destroy(&q->lock_in_head);
	pthread_mutex_destroy(&q->lock_in_tail);
	pthread_mutex_destroy(&q->lock_out_head);
	pthread_mutex_destroy(&q->lock_out_tail);
	sem_wrap_destroy(&q->sem_in);
	sem_wrap_destroy(&q->sem_out);
	free(q);
}

void d_enqueue_in(struct dqueue_t *q, void *ptr)
{
	pthread_mutex_lock(&q->lock_in_head);
	++q->in_head;
	if (q->in_head == q->cap)
		q->in_head = 0;
	q->in_data[q->in_head] = ptr;
	sem_wrap_post(&q->sem_in);
	pthread_mutex_unlock(&q->lock_in_head);
}

void d_enqueue_out(struct dqueue_t *q, void *ptr)
{
	pthread_mutex_lock(&q->lock_out_head);
	++q->out_head;
	if (q->out_head == q->cap)
		q->out_head = 0;
	q->out_data[q->out_head] = ptr;
	sem_wrap_post(&q->sem_out);
	pthread_mutex_unlock(&q->lock_out_head);
}

void *d_dequeue_in(struct dqueue_t *q)
{
	void *ret;
	pthread_mutex_lock(&q->lock_in_tail);
	sem_wrap_wait(&q->sem_in);
	ret = q->in_data[q->in_tail];
	++q->in_tail;
	if (q->in_tail == q->cap)
		q->in_tail = 0;
	pthread_mutex_unlock(&q->lock_in_tail);
	return ret;
}

void *d_dequeue_out(struct dqueue_t *q)
{
	void *ret;
	pthread_mutex_lock(&q->lock_out_tail);
	sem_wrap_wait(&q->sem_out);
	ret = q->out_data[q->out_tail];
	++q->out_tail;
	if (q->out_tail == q->cap)
		q->out_tail = 0;
	pthread_mutex_unlock(&q->lock_out_tail);
	return ret;
}

