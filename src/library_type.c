#include "library_type.h"

const int16_t n_type = 1;
const struct library_t protocol[] = {
        {16, 10} // 10X Chromium 5'
};

int8_t check_valid_library(const int16_t type)
{
        if (type < 0 || type >= n_type)
                return 0;
        return 1;
}

struct library_t get_library(const int16_t type)
{
        return protocol[type];
}