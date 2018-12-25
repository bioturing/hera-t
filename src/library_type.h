#ifndef _LIBRARY_TYPE_H_
#define _LIBRARY_TYPE_H_

#include <stdio.h>
#include <stdint.h>

// protocol index
#define CHROMIUM_5		0

struct library_t {
        int16_t bc_len;
        int16_t umi_len;
};

int8_t check_valid_library(const int16_t type);

struct library_t get_library(const int16_t type);

#endif