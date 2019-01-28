#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#if defined(_MSC_VER)
#include <time.h>
#include <windows.h>
#include <getopt.h>
#include <BaseTsd.h>
#else
#include <unistd.h>
#include <sys/resource.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#endif /* _MSC_VER */

#include "io_utils.h"
#include "verbose.h"

#ifdef _WIN32
int windows_path_convert(TCHAR *path) {
	if (*path == NULL) {
		fprintf(stderr, "errororrororor");
		return 1;
	}
	size_t plen = _tcslen(path);
	size_t i;
	for (i = 0; i < plen; ++i) {
		if (path[i] == _T('/')) {
			path[i] = _T('\\');
		}
	}
	return 0;
}
#endif

#ifdef _WIN32
FILE *_xfopen(const TCHAR *file_path, const TCHAR *mode)
#else
FILE *_xfopen(const char *file_path, const char *mode)
#endif
{
	FILE *fi = NULL;
#ifndef _MSC_VER
	fi = fopen(file_path, mode);
	if (!fi) {
		if (strcmp(mode, "r") == 0 || strcmp(mode, "rb") == 0)
			__ERROR("Unable to open file [%s] to read", file_path);
		else if (strcmp(mode, "w") == 0 || strcmp(mode, "wb") == 0)
			__ERROR("Unable to open file [%s] to write", file_path);
		else
			__ERROR("Unable to open file [%s]"
				"with unknown mode [%s]", file_path, mode);
	}
#else
	TCHAR *wpath = _tcsdup(file_path);
	//windows_path_convert(wpath);
	fi = _tfopen(wpath, mode);
	if (fi == NULL) {
		if (_tcscmp(mode, _T("r")) == 0 || _tcscmp(mode, _T("rb")) == 0)
			__ERROR("Unable to open file [%ls] to read", file_path);
		else if (_tcscmp(mode, _T("w")) == 0 || _tcscmp(mode, _T("wb")) == 0)
			__ERROR("Unable to open file [%ls] to write", file_path);
		else
			__ERROR("Unable to open file [%ls]"
				"with unknown mode [%ls]", file_path, mode);
	}
	__DEBUG("before lul\n");
	__DEBUG("wpath debug %ls\n", wpath);
	free(wpath);
	__DEBUG("after lul\n");
#endif
	return fi;
}

void xfclose(FILE *f)
{
	if (f) {
		fclose(f);
		f = NULL;
	}
}

size_t xfread(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fread(ptr, size, nmemb, stream);
	if (ret != nmemb)
		__ERROR("fread, wrong file or file is corrupted");
	return ret;
}

size_t xfwrite(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
	size_t ret = fwrite(ptr, size, nmemb, stream);
	if (ret != nmemb)
		__ERROR("fwrite, could not write data to file");
	return ret;
}

ssize_t xgetline(char **str, size_t *size, FILE *stream)
{
	ssize_t ret = -1;
#if defined(_MSC_VER)
	ret = _getline(str, size, stream);
#else
	ret = getline(str, size, stream);
#endif
	if (ret == 0 || ret == -1)
		return ret;
	if ((*str)[ret - 1] == '\n')
		(*str)[--ret] = '\0';
	return ret;
}

#ifdef _WIN32
void _normalize_dir(TCHAR *path)
#else
void _normalize_dir(char *path)
#endif
{
#ifdef _WIN32
	int len = _tcslen(path), i, j;
	for (i = 0; i < len - 1; ) {
		if (path[i] == _T('/') && path[i + 1] == _T('/')) {
			// what is this? memmove?
			for (j = i; j < len; ++j)
				path[j] = path[j + 1];
			--len;
		} else {
			++i;
		}
	}
#else
	int len = strlen(path), i, j;
	for (i = 0; i < len - 1; ) {
		if (path[i] == '/' && path[i + 1] == '/') {
			// what is this? memmove?
			for (j = i; j < len; ++j)
				path[j] = path[j + 1];
			--len;
		} else {
			++i;
		}
	}
#endif
}

#ifdef _WIN32
void _make_dir(const TCHAR *path)
#else
void _make_dir(const char *path)
#endif
{
#ifdef _WIN32
	struct stat st = { 0 };
	if (_tstat(path, &st) == -1) {
		if (_tmkdir(path)) {
			perror("Could not make output directory");
			exit(EXIT_FAILURE);
		}
	}
#else
	struct stat st = {0};
	if (stat(path, &st) == -1) {
		if (mkdir(path, 0700)) {
			perror("Could not make output directory");
			exit(EXIT_FAILURE);
		}
	}
#endif
}

#ifdef _WIN32
size_t _fetch_size(TCHAR **file_path, int n_file)
#else
size_t _fetch_size(char **file_path, int n_file)
#endif
{
	size_t ret = 0;
	FILE *fid;
	int i;
	for (i = 0; i < n_file; ++i) {
#ifdef _WIN32
		fid = xfopen(file_path[i], "rb");
#else
		fid = xfopen(file_path[i], "rb");
#endif
		if (!fid)
			__ERROR("Unable to open file: %s", file_path[i]);
		fseek(fid, 0L, SEEK_END);
		ret += ftell(fid);
		fclose(fid);
	}
	return ret;
}

/* ------------ shared_stream_t utils ------------ */

/* Init shared stream on n threads */
struct shared_fstream_t *init_shared_stream(const char *path, int n)
{
	struct shared_fstream_t *ret;
	FILE *fp;
	pthread_mutex_t *lock;
	fp = xfopen(path, "wb");
	ret = calloc(n, sizeof(struct shared_fstream_t));
	lock = calloc(1, sizeof(pthread_mutex_t));
	pthread_mutex_init(lock, NULL);
	int i;
	for (i = 0; i < n; ++i) {
		ret[i].fp = fp;
		ret[i].lock = lock;
		ret[i].buf = malloc(SFS_BUF_SZ + 1);
		ret[i].buf_len = 0;
	}
	return ret;
}

/* Destroy shared stream on n threads */
void destroy_shared_stream(struct shared_fstream_t *p, int n)
{
	if (!p) return;
	/* Flush remaining buffer */
	int i;
	for (i = 0; i < n; ++i) {
		if (p[i].buf_len) {
			/* Just to make color :D Actually don't need to lock */
			pthread_mutex_lock(p[i].lock);
			xfwrite(p[i].buf, 1, p[i].buf_len, p[i].fp);
			p[i].buf_len = 0;
			pthread_mutex_unlock(p[i].lock);
		}
		free(p[i].buf);
	}
	pthread_mutex_destroy(p->lock);
	free(p->lock);
	xfclose(p->fp);
	free(p);
}

/* flush buffer to output */
void sfs_flush(struct shared_fstream_t *p)
{
	pthread_mutex_lock(p->lock);
	xfwrite(p->buf, 1, p->buf_len, p->fp);
	pthread_mutex_unlock(p->lock);
	p->buf_len = 0;
}
