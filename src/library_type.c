#include "library_type.h"

const int16_t n_type = 1;
const int16_t protocol_bc[] = {16};
const int16_t protocol_umi[] = {10};

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