#ifndef _LIBRARY_TYPE_H_
#define _LIBRARY_TYPE_H_

#include <stdio.h>
#include <stdint.h>

// protocol index
#define CHROMIUM_V2		0
#define CHROMIUM_V3		1

struct library_t {
        int16_t bc_pos;
        int16_t bc_len;
        int16_t umi_pos;
        int16_t umi_len;
};

struct library_t get_library(char *platform);

struct library_t init_lib();

#endif