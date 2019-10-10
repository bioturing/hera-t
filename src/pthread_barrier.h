#if defined (__APPLE__) && defined (__MACH__)

#ifndef __PTHREAD_MAC_H__
#define __PTHREAD_MAC_H__

#define HAVE_STRUCT_TIMESPEC
#include <pthread.h>

typedef int pthread_barrierattr_t;

typedef struct
{
	pthread_mutex_t mutex;
	pthread_cond_t cond;
	int count;
	int tripCount;
} pthread_barrier_t;


int pthread_barrier_init(pthread_barrier_t *barrier, const pthread_barrierattr_t *attr, unsigned int count);

int pthread_barrier_destroy(pthread_barrier_t *barrier);

int pthread_barrier_wait(pthread_barrier_t *barrier);


#endif // __PTHREAD_MAC_H__

#endif // (__APPLE__ && __MACH__)
