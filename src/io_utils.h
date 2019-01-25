#ifndef _IO_UTILS_H_
#define _IO_UTILS_H_

#include <pthread.h>
#include <stdlib.h>
#include <stdio.h>
#include "attribute.h"
#if defined(_MSC_VER)
#include <BaseTsd.h>
#include <Windows.h>
#include <tchar.h>
#include <wchar.h>
typedef SSIZE_T ssize_t;
#endif

#ifdef _WIN32
#define xfopen(file, mode) _xfopen(file, _T(mode))
#define make_dir(path) _make_dir(path)
#else
#endif
/* shared stream struct
 * Usage: multiple threads write data to single file, provide each thread a
 * buffer
 */
#define SFS_BUF_SZ		SIZE_2MB
struct shared_fstream_t {
	FILE *fp;
	pthread_mutex_t *lock;
	char *buf;
	int buf_len;
};

#ifdef _WIN32
int windows_path_convert(TCHAR *path);
#endif

#ifdef _WIN32
FILE *_xfopen(const TCHAR *file_path, const TCHAR *mode);
#else
FILE *_xfopen(const char *file_path, const char *mode);
#endif

/* fflush before close file */
void xfclose(FILE *f);

/* check fread function read enough nmemb */
size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream);

/* check fwrite function write enough nmemb */
size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);

/* auto remove /n character if found */
ssize_t xgetline(char **str, size_t *size, FILE *stream);

/* get time */
double realtime();

/* make directory if is not exist */
#ifdef _WIN32
void _make_dir(const TCHAR *path);
#else
void _make_dir(const char *path);
#endif
/* ------------ shared_stream_t utils ------------ */

/* Init shared stream on n threads */
struct shared_fstream_t *init_shared_stream(const char *path, int n);

/* Destroy shared stream on n threads */
void destroy_shared_stream(struct shared_fstream_t *p, int n);

/* flush buffer to output */
void sfs_flush(struct shared_fstream_t *p);

#endif
