#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "verbose.h"

static FILE *log_file = NULL;

void init_log(const char *path)
{
	if (!log_file) {
		log_file = fopen(path, "wb");
		if (!log_file)
			__VERBOSE_INFO("WARNING", "Unable to open log file."
					" No log will be written\n");
		else
			__VERBOSE_INFO("INFO", "Log file: [%s]\n", path);
	}
}

void log_write(const char *fmt, ...)
{
	if (!log_file) return;

	va_list args;
	char *buffer = alloca(1024);
	va_start(args, fmt);
	vsprintf(buffer, fmt, args);
	fprintf(log_file, "%s", buffer);
	fflush(log_file);
	va_end(args);
}

void close_log()
{
	if (log_file)
		fclose(log_file);
}
