#include "library_type.h"
#include "verbose.h"

const int16_t n_type = 2;
const struct library_t protocol[] = {
	{0, 16, 16, 10},	// 10X Chromium 3' - v2
	{0, 16, 16, 12}	// 10X Chromium 3' - v3
};

struct library_t init_lib()
{
	struct library_t lib;
	lib.bc_pos = lib.umi_pos = -1;
	lib.bc_len = lib.umi_len = 0;

	return lib;
}

struct library_t get_library(char *platform)
{
	int type = -1;

	if (!strcmp(platform, "10x_v2"))
		type = CHROMIUM_V2;
	else if (!strcmp(platform, "10x_v3"))
		type = CHROMIUM_V3;
	if (type == -1)
		__ERROR("Invalid type of platform. Please just input one of these: 10x_v2, 10x_v3\n");

	return protocol[type];
}