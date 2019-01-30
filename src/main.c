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
		build_index(2, argc, argv);
	} else if (!strcmp(argv[1], "count")) {
		single_cell(2, argc, argv);
		// quant(2, argc, argv);
	} else {
		fprintf(stderr, "Invalid command!\n");
		print_usage();
	}

	fflush(stderr);
	return 0;
}
