#ifndef _ATOMIC_H
#define _ATOMIC_H

#if  defined(_WIN32) || defined(_WIN64)
#include <Windows.h>

#define __sync_fetch_and_add32 InterlockedExchangeAdd
#define __sync_val_compare_and_swap32(ptr, a, v) InterlockedCompareExchange(ptr, v, a)
#define __sync_bool_compare_and_swap32(ptr, a, v) (InterlockedCompareExchange(ptr, v, a) == (a))

#define __sync_fetch_and_add64 InterlockedExchangeAdd64
#define __sync_val_compare_and_swap64(ptr, a, v) InterlockedCompareExchange64(ptr, v, a)
#define __sync_bool_compare_and_swap64(ptr, a, v) (InterlockedCompareExchange64(ptr, v, a) == (a))

#else
#define __sync_fetch_and_add8 __sync_fetch_and_add
#define __sync_val_compare_and_swap8(ptr, a, v) __sync_val_compare_and_swap(ptr, a, v)
#define __sync_bool_compare_and_swap8(ptr, a, v) __sync_bool_compare_and_swap(ptr, a, v)

#define __sync_fetch_and_add32 __sync_fetch_and_add
#define __sync_val_compare_and_swap32(ptr, a, v) __sync_val_compare_and_swap(ptr, a, v)
#define __sync_bool_compare_and_swap32(ptr, a, v) __sync_bool_compare_and_swap(ptr, a, v)

#define __sync_fetch_and_add64 __sync_fetch_and_add
#define __sync_val_compare_and_swap64(ptr, a, v) __sync_val_compare_and_swap(ptr, a, v)
#define __sync_bool_compare_and_swap64(ptr, a, v) __sync_bool_compare_and_swap(ptr, a, v)
#endif
#endif
