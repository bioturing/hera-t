#include "opt.h"
#include "index.h"
#include "single_cell.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

// static int check_license()
// {
// 	static const int exp_day = 19;
// 	static const int exp_mon = 7;
// 	static const int exp_year = 2019;
// 	time_t t;
// 	t = time(NULL);
// 	struct tm tm;
// 	localtime_r(&t, &tm);
// 	if (tm.tm_year + 1900 > exp_year)
// 		return 0;
// 	if (tm.tm_year + 1900 == 2019 && tm.tm_mon + 1 > exp_mon)
// 		return 0;
// 	if (tm.tm_year + 1900 == 2019 && tm.tm_mon + 1 == exp_mon && tm.tm_mday >= exp_day)
// 		return 0;
// 	return 1;
// 	// fprintf(stderr, "day = %d; month = %d; year = %d\n", tm.tm_mday, tm.tm_mon, tm.tm_year + 1900);
// }

int main(int argc, char *argv[])
{
	print_info();
	// if (!check_license()) {
	// 	fprintf(stderr, "Your 14-day trial license has expired.\n");
	// 	return 0;
	// }
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
