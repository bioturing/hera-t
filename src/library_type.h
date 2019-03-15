#ifndef _LIBRARY_TYPE_H_
#define _LIBRARY_TYPE_H_

#include <stdio.h>
#include <stdint.h>

// protocol index
#define CHROMIUM3_V2		0
#define CHROMIUM3_V3		1
#define CELSEE			1

struct library_t {
	int16_t bc_len;
	int16_t umi_len;
};

int8_t check_valid_library(const int16_t type);

struct library_t get_library(const int16_t type);

#endif