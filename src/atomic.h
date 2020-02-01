#ifndef _ATOMIC_H
#define _ATOMIC_H

#if  defined(_WIN32) || defined(_WIN64)
#include <Windows.h>
#include <intrin.h>

#define __sync_fetch_and_add32 _InterlockedAdd
#define __sync_val_compare_and_swap32(ptr, a, v) _InterlockedCompareExchange(ptr, v, a)
#define __sync_bool_compare_and_swap32(ptr, a, v) (_InterlockedCompareExchange(ptr, v, a) == (a))

#define __sync_fetch_and_add64 _InterlockedAdd64
#define __sync_val_compare_and_swap64(ptr, a, v) _InterlockedCompareExchange64(ptr, v, a)
#define __sync_bool_compare_and_swap64(ptr, a, v) (_InterlockedCompareExchange64(ptr, v, a) == (a))

#define __sync_fetch_and_add8 _InterlockedAdd8
#define __sync_val_compare_and_swap8(ptr, a, v) _InterlockedCompareExchange8(ptr, v, a)
#define __sync_bool_compare_and_swap8(ptr, a, v) (_InterlockedCompareExchange8(ptr, v, a) == (a))

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

#define atomic_add_and_fetch32(a, x)		__sync_add_and_fetch(a, x)
#define atomic_add_and_fetch64(a, x)		__sync_add_and_fetch(a, x)

#define atomic_val_CAS8(ptr, old_val, new_val)				       \
			__sync_val_compare_and_swap(ptr, old_val, new_val)
#define atomic_val_CAS32(ptr, old_val, new_val)				       \
			__sync_val_compare_and_swap(ptr, old_val, new_val)
#define atomic_val_CAS64(ptr, old_val, new_val)				       \
			__sync_val_compare_and_swap(ptr, old_val, new_val)
#define atomic_bool_CAS8(ptr, old_val, new_val)			       \
			__sync_bool_compare_and_swap(ptr, old_val, new_val)
#define atomic_bool_CAS32(ptr, old_val, new_val)			       \
			__sync_bool_compare_and_swap(ptr, old_val, new_val)
#define atomic_bool_CAS64(ptr, old_val, new_val)			       \
			__sync_bool_compare_and_swap(ptr, old_val, new_val)

#endif
#endif
