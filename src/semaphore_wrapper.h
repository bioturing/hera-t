#ifndef __SEMAPHORE_WRAPPER_H__
#define __SEMAPHORE_WRAPPER_H__

#include <stdint.h>

#if defined (__APPLE__) && defined (__MACH__)
#include <dispatch/dispatch.h>
#else
#include <semaphore.h>
#endif

struct sem_wrap_t {
#if defined (__APPLE__) && defined (__MACH__)
	dispatch_semaphore_t sem;
#else
	sem_t sem;
#endif
};

void sem_wrap_init(struct sem_wrap_t *sem, uint32_t value);

void sem_wrap_wait(struct sem_wrap_t *sem);

void sem_wrap_post(struct sem_wrap_t *sem);

void sem_wrap_destroy(struct sem_wrap_t *sem);

#endif /* __SEMAPHORE_WRAPPER_H__ */
