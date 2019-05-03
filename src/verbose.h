#ifndef _VERBOSE_H_
#define _VERBOSE_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#if defined(_MSC_VER)
#define __VERBOSE_INFO(tag, fmt, ...) do {				\
	fprintf(stderr, "[" tag "] " fmt, __VA_ARGS__);			\
	fflush(stderr);							\
} while (0) /* VERBOSE_INFO */

#define __VERBOSE_LOG(tag, fmt, ...) do {				\
	fprintf(stderr, "[" tag "] " fmt, __VA_ARGS__);			\
	log_write("[" tag "] " fmt, __VA_ARGS__);			\
} while (0) /* VERBOSE_AND_LOG */

#define __VERBOSE(fmt, ...) do {					\
	fprintf(stderr, fmt, __VA_ARGS__);				\
	fflush(stderr);							\
} while (0) /* VERBOSE */

#if defined(NDEBUG)
#define __DEBUG(fmt, ...) 0
#else
#define __DEBUG(fmt, ...) do {						\
	fprintf(stderr, "[DEBUG] " fmt, __VA_ARGS__);			\
	fflush(stderr);							\
} while (0) /* __DEBUG */
#endif /* NDEBUG */

/* FIXME: Write error handler: i.e. delete temporary files */
#define __ERROR(fmt, ...)do {						\
	fprintf(stderr, "[ERROR] " fmt "\n", __VA_ARGS__);		\
	log_write("[ERROR] " fmt "\n", __VA_ARGS__);			\
	exit(EXIT_FAILURE);						\
} while(0) /* ERROR */

#else /* __MSC_VER */

#define __VERBOSE_INFO(tag, fmt, args...) do {				\
	fprintf(stderr, "[" tag "] " fmt, ##args);			\
	fflush(stderr);							\
} while (0) /* VERBOSE_INFO */

#define __VERBOSE_LOG(tag, fmt, args...) do {				\
	fprintf(stderr, "[" tag "] " fmt, ##args);			\
	log_write("[" tag "] " fmt, ##args);				\
} while (0) /* VERBOSE_AND_LOG */

#define __VERBOSE(fmt, args...) do {					\
	fprintf(stderr, fmt, ##args);					\
	fflush(stderr);							\
} while (0) /* VERBOSE */

#if defined(NDEBUG)
#define __DEBUG(fmt, ...) 0
#else
#define __DEBUG(fmt, args...)	fprintf(stderr, "[DEBUG] " fmt, ##args)
#endif /* NDEBUG */

#define __ERROR(fmt, args...) do {					\
	fprintf(stderr, "[ERROR] " fmt "\n", ##args);			\
	log_write("[ERROR] " fmt "\n", ##args);				\
	exit(EXIT_FAILURE);						\
} while(0) /* ERROR */
#endif /* __MSC_VER */

void init_log(const char *path);

void log_write(const char *fmt, ...);

void close_log();

#endif /* _VERBOSE_H_ */
