#include <stdlib.h>
#include "library_type.h"

#if defined(_MSC_VER)
#define __LIB_ERROR(fmt, ...) do {					       \
	fprintf(stderr, "ERROR PARSING OPTIONS: " fmt "\n", __VA_ARGS__);      \
	exit(EXIT_FAILURE);						       \
} while(0)
#define __LIB_WARNING(fmt, ...) do {					       \
	fprintf(stderr, "[WARNING]" fmt "\n", __VA_ARGS__);                    \
	exit(EXIT_FAILURE);						       \
} while(0)
#else 
#define __LIB_ERROR(fmt, args...) do {					       \
	fprintf(stderr, "ERROR PARSING OPTIONS: " fmt "\n", ##args);	       \
	exit(EXIT_FAILURE);						       \
} while(0)
#define __LIB_WARNING(fmt, args...) do {				       \
	fprintf(stderr, "[WARNING]" fmt "\n", ##args);                         \
	exit(EXIT_FAILURE);						       \
} while(0)
#endif

const int16_t n_type = 2;
const struct library_t protocol[] = {
        {16, 10},       // 10X Chromium 3' - v2
        {16, 12}        // 10X Chromium 3' - v3
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