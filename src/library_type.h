#ifndef _LIBRARY_TYPE_H_
#define _LIBRARY_TYPE_H_

#include <stdio.h>
#include <stdint.h>

// protocol index
#define CHROMIUM_5		0

const int16_t n_type = 1;
const int16_t protocol_bc[] = {16};
const int16_t protocol_umi[] = {10};

struct library_t {
        int16_t bc_len;
        int16_t umi_len;
};

struct library_t get_library(const int16_t type)
{
        struct library_t lib;

        lib.bc_len = -1;

        if (type >= 0 && type < n_type){
                lib.bc_len = protocol_bc[type];
                lib.umi_len = protocol_umi[type];
        }

        return lib;
}

#endif