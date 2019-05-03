#include "opt.h"
#include "index.h"
#include "single_cell.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

int main(int argc, char *argv[])
{
	print_info();
	if (argc == 1) {
		print_usage();
	} else if (!strcmp(argv[1], "index")) {
		build_index(argc - 2, argv + 2);
	} else if (!strcmp(argv[1], "count")) {
		single_cell(argc - 2, argv + 2);
	} else {
		fprintf(stderr, "Invalid command!\n");
		print_usage();
	}

	fflush(stderr);
	return 0;
}
