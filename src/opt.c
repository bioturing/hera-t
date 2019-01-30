#include "attribute.h"
#include "io_utils.h"
#include "library_type.h"
#include "opt.h"
#include "utils.h"
#include "verbose.h"

#if defined(_MSC_VER)
#define __OPT_ERROR(fmt, ...) do {					       \
	fprintf(stderr, "ERROR PARSING OPTIONS: " fmt "\n", __VA_ARGS__);	       \
	exit(EXIT_FAILURE);						       \
} while(0)
#else 
#define __OPT_ERROR(fmt, args...) do {					       \
	fprintf(stderr, "ERROR PARSING OPTIONS: " fmt "\n", ##args);	       \
	exit(EXIT_FAILURE);						       \
} while(0)
#endif

void print_info()
{
	__VERBOSE("Hera-T is a program developed by BioTuring for ScRNA-Seq analysis\n");
	__VERBOSE("Please contact info@bioturing.com if you need further support\n");
	__VERBOSE("This is an academic version, which is free for academic labs\n");
	__VERBOSE("No IP or commercial work can be derived from using this free academic version\n");
	__VERBOSE("If you are using Hera-T for IP or commercial related work, please contact info@bioturing.com to obtain a license\n");
	__VERBOSE("Cite Hera-T paper at: https://doi.org/10.1101/530501\n");
	__VERBOSE("Version: %d.%d.%d\n", PROG_VERSION_MAJOR, PROG_VERSION_MINOR, PROG_VERSION_FIX);
	__VERBOSE("\n");
}

void print_usage()
{
	// print_info();
	__VERBOSE("Usage: ./hera-T <CMD> [options] ...\n");
	__VERBOSE("\n");
	__VERBOSE("Where <CMD> can be one of:\n");
	__VERBOSE("    index            to build Hera-T index\n");
	__VERBOSE("    count            to calculate gene count matrix\n");
	__VERBOSE("\n");
}

void print_index_usage()
{
	// print_info();
	__VERBOSE("To build Hera-T index\n");
	__VERBOSE("\n");
	__VERBOSE("Usage: ./hera-T index -g <path/to/genome_fasta> -t <path/to/gene_gtf> -p <index_prefix> -o <output_folder>\n");
	__VERBOSE("Example: ./hera-T index -g Homo_sapiens.GRCh37.75.dna_sm.primary_assembly.fa -t Homo_sapiens.GRCh37.75.gtf -o index -p grch37\n");
	__VERBOSE("\n");
}

void print_count_usage()
{
	// print_info();
	__VERBOSE("To calculate gene count table\n");
	__VERBOSE("\n");
	__VERBOSE("Usage: ./hera-T count [options] -x <idx_name> -1 <R1> -2 <R2>\n");
	__VERBOSE("Option:\n");
	__VERBOSE("-t\t: Number of threads\n");
	__VERBOSE("-o\t: Output directory name\n");
	__VERBOSE("-p\t: Output file prefix\n");
	__VERBOSE("-l\t: Library types\n");
	__VERBOSE("\t\t%u: 10X-Chromium 3' (v2) protocol\n", CHROMIUM3_V2);
	__VERBOSE("\t\t%u: 10X-Chromium 3' (v3) protocol\n", CHROMIUM3_V3);
	__VERBOSE("Example: ./hera-T count -t 32 -o ./result -x index/grch37 -l 0 -1 lane_0.read_1.fq lane_1.read_1.fq -2 lane_0.read_2.fq lane_1.read_2.fq\n");
	//__VERBOSE("--count-intron\t: Count both exonic and intronic reads\n");
	__VERBOSE("\n");
}

static struct opt_count_t *init_opt_count()
{
	struct opt_count_t *opt;
	opt = calloc(1, sizeof(struct opt_count_t));
	opt->out_dir = ".";
	opt->temp_dir = ".";
	opt->prefix = "HeraT";
	opt->n_threads = 1;
	opt->lib = get_library(0);
	opt->n_threads = 1;
	opt->is_dump_align = 0;
	opt->count_intron = 0;
	opt->left_file = opt->right_file = NULL;
	return opt;
}

static struct opt_index_t *init_opt_index()
{
	struct opt_index_t *opt;
	opt = calloc(1, sizeof(struct opt_index_t));
	opt->genome = NULL;
	opt->gtf = NULL;
	opt->prefix = NULL;
	opt->idx_dir = "./";
	opt->k = 29;
	opt->bwt = 1;
	return opt;
}

static void check_valid_opt_index(struct opt_index_t *opt)
{
	if (!opt->genome)
		__OPT_ERROR("Missing -g argument");

	if (!opt->gtf)
		__OPT_ERROR("Missing -t argument");

	if (!opt->prefix)
		__OPT_ERROR("Missing -p argument");
}

static void check_valid_opt_count(struct opt_count_t *opt)
{
	if (opt->n_threads < 1)
		__OPT_ERROR("Invalid  number of threads: %d", opt->n_threads);

	if (opt->index == NULL)
		__OPT_ERROR("Missing index");

	if (opt->n_files == 0)
		__OPT_ERROR("Missing input files");

	if (opt->left_file == NULL)
		__OPT_ERROR("Missing first segment of read files");

	if (opt->right_file == NULL)
		__OPT_ERROR("Missing second segment of read files");
}

static void opt_check_num(int argc, char **argv)
{
	if (argc < 2 || argv[1][0] == '-')
		__OPT_ERROR("Missing data for option %s", argv[0]);

	int i;
	for (i = 0; argv[1][i]; ++i)
		if (argv[1][i] < '0' || argv[1][i] > '9')
			__OPT_ERROR("Invalid data for option %s: %s", argv[0], argv[1]);
}

static void opt_check_str(int argc, char **argv)
{
	if (argc < 2 || argv[1][0] == '-')
		__OPT_ERROR("Missing data for option %s", argv[0]);
}

static int opt_count_list(int argc, char **argv)
{
	int n;
	for (n = 0; n < argc - 1; ++n) {
		if (argv[n + 1][0] == '-')
			break;
	}
	if (n == 0)
		__OPT_ERROR("Missing data for option %s", argv[0]);
	return n;
}

struct opt_index_t *get_opt_index(int argc, char *argv[])
{
	if (!argc) {
		print_index_usage();
		exit(EXIT_FAILURE);
	}

	struct opt_index_t *opt = init_opt_index();

	int pos = 0;
	while (pos < argc) {
		if (!strcmp(argv[pos], "-t")) {
			opt_check_str(argc - pos, argv + pos);
			opt->gtf = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-p")) {
			opt_check_str(argc - pos, argv + pos);
			opt->prefix = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-g")) {
			opt_check_str(argc - pos, argv + pos);
			opt->genome = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-o")) {
			opt_check_str(argc - pos, argv + pos);
			opt->idx_dir = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-k")) {
			opt_check_num(argc - pos, argv + pos);
			opt->k = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "--no-bwt")) {
			opt->bwt = 0;
			++pos;
		} else if (!strcmp(argv[pos], "-h")) {
			print_index_usage();
		} else {
			__OPT_ERROR("Invalid option %s", argv[pos]);
		}
	}
	check_valid_opt_index(opt);
	int i;
	for (i = 0; opt->prefix[i]; ++i) {
		if ((opt->prefix[i] < 'A' || opt->prefix[i] > 'Z') &&
		    (opt->prefix[i] < 'a' || opt->prefix[i] > 'z') &&
		    (opt->prefix[i] < '0' || opt->prefix[i] > '9') &&
		    !strchr("_-.", opt->prefix[i]))
			__OPT_ERROR("Invalid prefix string %s (included character '%c')", opt->prefix, opt->prefix[i]);
	}

	// normalize_dir(opt->idx_dir);
	make_dir(opt->idx_dir);
	return opt;
}

struct opt_count_t *get_opt_count(int argc, char *argv[])
{
	if (argc == 0) {
		print_count_usage();
		exit(EXIT_FAILURE);
	}

	struct opt_count_t *opt;
	opt = init_opt_count();
	int pos = 0, n;
	while (pos < argc) {
		if (!strcmp(argv[pos], "-x")) {
			opt_check_str(argc - pos, argv + pos);
			opt->index = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-o")) {
			opt_check_str(argc - pos, argv + pos);
			opt->out_dir = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "--temp-dir")) {
			opt_check_str(argc - pos, argv + pos);
			opt->temp_dir = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-t")) {
			opt_check_num(argc - pos, argv + pos);
			opt->n_threads = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "-p")) {
			opt_check_str(argc - pos, argv + pos);
			opt->prefix = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-1")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__OPT_ERROR("Number of files in pair are not equal");
			opt->n_files = n;
			opt->left_file = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-2")) {
			n = opt_count_list(argc - pos, argv + pos);
			if (opt->n_files > 0 && opt->n_files != n)
				__OPT_ERROR("Number of files in pair are not equal");
			opt->n_files = n;
			opt->right_file = argv + pos + 1;
			pos += (n + 1);
		} else if (!strcmp(argv[pos], "-l")) {
			opt_check_num(argc - pos, argv + pos);
			int type = atoi(argv[pos + 1]);
			if (!check_valid_library(type))
				__OPT_ERROR("Invalid library type");
			opt->lib = get_library(type);
			pos += 2;
		} else if (!strcmp(argv[pos], "--dump-align")) {
			opt->is_dump_align = 1;
			++pos;
		}  else if (!strcmp(argv[pos], "--count-intron")) {
			opt->count_intron = 1;
			++pos;
		} else {
			__OPT_ERROR("Invalid option %s", argv[pos]);
		}
	}
	check_valid_opt_count(opt);
	int i;
	for (i = 0; opt->prefix[i]; ++i) {
		if ((opt->prefix[i] < 'A' || opt->prefix[i] > 'Z') &&
		    (opt->prefix[i] < 'a' || opt->prefix[i] > 'z') &&
		    (opt->prefix[i] < '0' || opt->prefix[i] > '9') &&
		    !strchr("_-.", opt->prefix[i]))
			__OPT_ERROR("Invalid prefix string %s (included character '%c')", opt->prefix, opt->prefix[i]);
	}
	// normalize_dir(opt->out_dir);
	make_dir(opt->out_dir);
	return opt;
}
