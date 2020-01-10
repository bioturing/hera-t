//
// Created by BioTuring on 2019-09-27.
//

/*
 * Copyright (c) 2017 rxi
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include <sys/time.h>
#include <sys/resource.h>
#include <sys/uio.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <stdint.h>

#include "log.h"

#define __min(a, b) 		((a) < (b) ? (a) : (b))

#define LOG_PADDING 30

static struct {
	void *udata;
	log_LockFn lock;
	FILE *fp;
	int level;
	int quiet;
	struct rusage usage;
	time_t start_time;
	time_t last_time;
	unsigned int mem_used;
	char *stage;
	int fd;
} L;


static const char *level_names[] = {
	"TRACE", "DEBUG_TECH", "DEBUG", "INFO", "WARN", "ERROR", "FATAL"
};

#ifdef LOG_USE_COLOR
static const char *level_colors[] = {
	"\x1b[94m", "\x1b[36m", "\x1b[36m", "\x1b[32m", "\x1b[33m", "\x1b[31m", "\x1b[35m"
};
#endif


static void lock(void)
{
	if (L.lock) {
		L.lock(L.udata, 1);
	}
}

static void unlock(void)
{
	if (L.lock) {
		L.lock(L.udata, 0);
	}
}


void log_set_udata(void *udata)
{
	L.udata = udata;
}


void log_set_lock(log_LockFn fn)
{
	L.lock = fn;
}


void log_set_fp(FILE *fp)
{
	L.fp = fp;
}


void log_set_level(int level)
{
	L.level = level;
}

void init_logger(int level, const char * file_path)
{
	FILE *fp = fopen(file_path, "w");
	log_set_fp(fp);
	log_set_level(__min(level, LOG_INFO));
	L.start_time = time(NULL);
	L.stage = "General";
	L.fd = fileno(L.fp);
}

void set_log_stage(char *stage)
{
	L.stage = stage;
}

void log_change_file(const char *file_path)
{
	if (L.fp) {
		fclose(L.fp);
	}
	L.fp = fopen(file_path, "w");
}

void close_logger()
{
	fclose(L.fp);
}

void log_set_quiet(int enable)
{
	L.quiet = enable ? 1 : 0;
}

void log_log(int level, const char *file, int line, const char *fmt, ...) {
	/* Get current time */
	time_t usr_time, t;
	t = time(NULL);
	struct tm *lt = localtime(&t);
	char src_code[128];
	unsigned int ru_maxrss;

	usr_time = t - L.start_time;
	struct tm *lt_real = localtime(&usr_time);
	/* Get used time and memory */
	if (level > LOG_TRACE) {
		getrusage(RUSAGE_SELF, &L.usage); /* Get resource usage, only available for UNIX */
		ru_maxrss = L.usage.ru_maxrss;
		L.mem_used = ru_maxrss;
	} else{
		ru_maxrss = L.mem_used;
	}

	/* Cast the file name and line number into a file to easily align */
	sprintf(src_code, "%s:%d", file, line);
	/* Log to file . Always log to file*/
	if (L.fp) {
		char file_buf[8192];
		va_list args;
		char buf[32];
		buf[strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M:%S", lt)] = '\0';
		sprintf(file_buf, "%s %-5s %-40s\t%-20s\t%.2f GB\t", buf, level_names[level], src_code,
			L.stage, (double)ru_maxrss/(1<<20));
		va_start(args, fmt);
		vsprintf(file_buf + strlen(file_buf), fmt, args);
		va_end(args);
		sprintf(file_buf + strlen(file_buf), "\n");
		struct iovec iov;
		iov.iov_base = file_buf;
		iov.iov_len = strlen(file_buf);
		writev(L.fd, &iov, 1);
	}

	if (level < L.level) {
		return;
	}

	/* Acquire lock */
	lock();

	/* Log to stderr */
	if (!L.quiet) {
		va_list args;
		char buf[16];
		buf[strftime(buf, sizeof(buf), "%H:%M:%S", lt_real)] = '\0';
#ifdef LOG_USE_COLOR
		fprintf(
      stderr, "%s %s%-5s\x1b[0m \x1b[90m%-40s\x1b[0m\t%-20s\t%.2f GB\t",
      buf, level_colors[level], level_names[level], src_code,  L.stage, (double)ru_maxrss/(1<<20));
#else
		fprintf(stderr, "%s %-5s %-40s\t%-20s\t%.2fGB\t", buf, level_names[level], src_code,
		        L.stage, (double)ru_maxrss/(1<<20));
#endif
		va_start(args, fmt);
		vfprintf(stderr, fmt, args);
		va_end(args);
		fprintf(stderr, "\n");
		fflush(stderr);
	}

	if (level >= LOG_ERROR) {
		perror("Something went wrong, please check the log file or send it to tan@bioturing.com");
		exit(1);
	}

	/* Release lock */
	unlock();
}
