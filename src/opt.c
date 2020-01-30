#include "opt.h"
#include "verbose.h"
#include "attribute.h"
#include "readTSV.h"

void print_info()
{
	__VERBOSE("Hera-T is a program developed by BioTuring for scRNA-Seq analysis\n");
	__VERBOSE("Please contact info@bioturing.com if you need further support\n");
	__VERBOSE("This is an academic version, which is free for academic labs\n");
	__VERBOSE("No IP or commercial work can be derived from using this free academic version\n");
	__VERBOSE("If you are using Hera-T for IP or commercial related work, please contact\ninfo@bioturing.com to obtain a license\n");
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
	__VERBOSE("Usage: ./hera-T count [options] @meta_file.tsv\n");
	__VERBOSE("Option:\n");
	__VERBOSE("-t\t: Number of threads\n");
	__VERBOSE("-o\t: Output directory name\n");

	__VERBOSE("@meta_file content:\n");
	__VERBOSE("rna_dir\tPath to fastqs directory\n");
	__VERBOSE("rna_name\tSample name (which is the prefix of fastqs)\n");
	__VERBOSE("rna_index\tPath_to_index_directory/index_prefix\n");
	__VERBOSE("cell_adt_dir\tPath to fastqs directory\n");
	__VERBOSE("cell_adt_name\tSample name (which is the prefix of fastqs)\n");
	__VERBOSE("cell_adt_ref\tPath_to_reference.csv\n");
	__VERBOSE("protein_adt_dir\tPath to fastqs directory\n");;
	__VERBOSE("protein_adt_name\tSample name (which is the prefix of fastqs)\n");
	__VERBOSE("protein_adt_ref\tPath_to_reference.csv\n");
	__VERBOSE("crispr_dir\tPath to fastqs directory\n");;
	__VERBOSE("crispr_name\tSample name (which is the prefix of fastqs)\n");
	__VERBOSE("crispr_ref\tPath_to_reference.csv\n");
	__VERBOSE("platform\tSequencing platform (10x_v2, 10x_v3)\n");

	__VERBOSE("\nIf your sequencing platform is not listed, you can define barcode and UMI as follow:\n");
	__VERBOSE("barcode_start_pos\tBarcode start position on read 1\n");
	__VERBOSE("barcode_len\tBarcode length\n");
	__VERBOSE("UMI_start_pos\tUMI start position on read 1\n");
	__VERBOSE("UMI_len\tUMI length\n");
	__VERBOSE("\n");
}

static void check_num(int argc, char **argv)
{
	if (argc < 2 || argv[1][0] == '-')
		__ERROR("Missing data for option %s", argv[0]);

	int i;
	for (i = 0; argv[1][i]; ++i)
		if (argv[1][i] < '0' || argv[1][i] > '9')
			__ERROR("Invalid data for option %s: %s", argv[0], argv[1]);
}

static void check_str(int argc, char **argv)
{
	if (argc < 2 || argv[1][0] == '-')
		__ERROR("Missing data for option %s", argv[0]);
}

/*****************************************
*           OPTION FOR INDEXING          *
*****************************************/

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
		__ERROR("Missing -g argument");

	if (!opt->gtf)
		__ERROR("Missing -t argument");

	if (!opt->prefix)
		__ERROR("Missing -p argument");
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
			check_str(argc - pos, argv + pos);
			opt->gtf = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-p")) {
			check_str(argc - pos, argv + pos);
			opt->prefix = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-g")) {
			check_str(argc - pos, argv + pos);
			opt->genome = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-o")) {
			check_str(argc - pos, argv + pos);
			opt->idx_dir = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-k")) {
			check_num(argc - pos, argv + pos);
			opt->k = atoi(argv[pos + 1]);
			pos += 2;
		} else if (!strcmp(argv[pos], "--no-bwt")) {
			opt->bwt = 0;
			++pos;
		} else if (!strcmp(argv[pos], "-h")) {
			print_index_usage();
		} else {
			__ERROR("Invalid option %s", argv[pos]);
		}
	}

	check_valid_opt_index(opt);

	int i;
	for (i = 0; opt->prefix[i]; ++i) {
		if ((opt->prefix[i] < 'A' || opt->prefix[i] > 'Z') &&
		    (opt->prefix[i] < 'a' || opt->prefix[i] > 'z') &&
		    (opt->prefix[i] < '0' || opt->prefix[i] > '9') &&
		    !strchr("_-.", opt->prefix[i]))
			__ERROR("Invalid prefix string %s (included character '%c')", opt->prefix, opt->prefix[i]);
	}

	make_dir(opt->idx_dir);
	return opt;
}

/*****************************************
*           OPTION FOR COUNTING          *
*****************************************/

static struct opt_count_t *init_opt_count()
{
	struct opt_count_t *opt = calloc(1, sizeof(struct opt_count_t));
	opt->out_dir = ".";
	opt->n_threads = 1;

	opt->lib = init_lib();
	opt->count_intron = 0;
	opt->rna = opt->cell = opt->protein = opt->crispr = NULL;

	return opt;
}

static struct input_t *init_input_data()
{
	struct input_t *input = calloc(1, sizeof(struct input_t));
	
	input->in_dir = input->name = NULL;
	input->n_files = 0;
	input->left_file = input->right_file = NULL;
	input->ref = NULL;

	return input;
}

void destroy_input_data(struct input_t *input)
{
	int i;
	if (input->in_dir)
		free(input->in_dir);
	if (input->name)
		free(input->name);

	if (input->left_file) {
		for (i = 0; i < input->n_files; ++i)
			free(input->left_file[i]);
		free(input->left_file);
	}

	if (input->right_file) {
		for (i = 0; i < input->n_files; ++i)
			free(input->right_file[i]);
		free(input->right_file);
	}

	if (input->ref)
		free(input->ref);

	free(input);
}

int check_valid_input(struct input_t *input, char *type)
{
	if (!input->in_dir && !input->name && !input->ref) {
		destroy_input_data(input);
		return 0;
	}

	if (!input->in_dir)
		__ERROR("Missing input directory for %s\n", type);
	if (!input->name)
		__ERROR("Missing sample name for %s\n", type);
	if (input->ref == NULL)
		__ERROR("Missing reference for %s\n", type);

	parse_dir(input);
	return 1;
}

void check_conflict(struct opt_count_t *opt)
{
	if (opt->lib.bc_pos != -1) {
		__VERBOSE("Already has barcode pos: %d\n", opt->lib.bc_pos);
		__ERROR("Barcode start position is defined twice, please just input only platform or barcode_start_pos\n");
	}
	if (opt->lib.bc_len)
		__ERROR("Barcode length is defined twice, please just input only platform or barcode_len\n");
	if (opt->lib.umi_pos != -1)
		__ERROR("UMI start position is defined twice, please just input only platform or UMI_start_pos\n");
	if (opt->lib.umi_len)
		__ERROR("UMI length is defined twice, please just input only platform or UMI_len\n");
}

void check_valid_lib(struct library_t lib)
{
	if (lib.bc_pos == -1)
		__ERROR("Barcode start position is not defined\n");
	if (lib.bc_len == 0)
		__ERROR("Barcode length is not defined\n");
	if (lib.umi_pos == -1)
		__ERROR("UMI start position is not defined\n");
	if (lib.umi_len == 0)
		__ERROR("UMI length is not defined\n");
}

char *get_string(struct content_t value)
{
	char *ref = malloc(value.l + 1);
	memcpy(ref, value.s, value.l);
	ref[value.l] = '\0';

	return ref;
}

void parse_input_meta(struct opt_count_t *opt, char *meta_file)
{
	struct tsv_t *meta = init_readTSV(meta_file, '\t');
	struct content_t param, value;
	struct input_t *rna = init_input_data();
	struct input_t *cell = init_input_data();
	struct input_t *protein = init_input_data();
	struct input_t *crispr = init_input_data();

	while(get_row_content(meta) > 0) {
		param = get_col_content(meta, 0);
		if (!param.l)
			continue;

		value = get_col_content(meta, 1);
		if (!value.l)
			__ERROR("Missing input value for param %s\n", param.s);
		
		if (!strncmp(param.s, "platform", param.l)) {
			check_conflict(opt);
			opt->lib = get_library(value.s);
		} else if (!strncmp(param.s, "barcode_start_pos", param.l)) {
			opt->lib.bc_pos = atoi(value.s);
		} else if (!strncmp(param.s, "barcode_len", param.l)) {
			opt->lib.bc_len = atoi(value.s);
		} else if (!strncmp(param.s, "UMI_start_pos", param.l)) {
			opt->lib.umi_pos = atoi(value.s);
		} else if (!strncmp(param.s, "UMI_len", param.l)) {
			opt->lib.umi_len = atoi(value.s);
		} else if (!strncmp(param.s, "rna_dir", param.l)) {
			rna->in_dir = get_string(value);
		} else if (!strncmp(param.s, "rna_name", param.l)) {
			rna->name = get_string(value);
		} else if (!strncmp(param.s, "rna_index", param.l)) {
			rna->ref = get_string(value);
		} else if (!strncmp(param.s, "cell_adt_dir", param.l)) {
			cell->in_dir = get_string(value);
		} else if (!strncmp(param.s, "cell_adt_name", param.l)) {
			cell->name = get_string(value);
		} else if (!strncmp(param.s, "cell_adt_ref", param.l)) {
			cell->ref = get_string(value);
		} else if (!strncmp(param.s, "protein_adt_dir", param.l)) {
			protein->in_dir = get_string(value);
		} else if (!strncmp(param.s, "protein_adt_name", param.l)) {
			protein->name = get_string(value);
		} else if (!strncmp(param.s, "protein_adt_ref", param.l)) {
			protein->ref = get_string(value);
		} else if (!strncmp(param.s, "crispr_dir", param.l)) {
			crispr->in_dir = get_string(value);
		} else if (!strncmp(param.s, "crispr_name", param.l)) {
			crispr->name = get_string(value);
		} else if (!strncmp(param.s, "crispr_ref", param.l)) {
			crispr->ref = get_string(value);
		}
	}

	check_valid_lib(opt->lib);

	if (check_valid_input(rna, "gene expression"))
		opt->rna = rna;
	if (check_valid_input(cell, "cell hashing"))
		opt->cell = cell;
	if (check_valid_input(protein, "protein measurement"))
		opt->protein = protein;
	if (check_valid_input(crispr, "CRISPR measurement"))
		opt->crispr = crispr;

	destroy_readTSV(meta);
}

struct opt_count_t *get_opt_count(int argc, char *argv[])
{
	if (argc == 0) {
		print_count_usage();
		exit(EXIT_FAILURE);
	}

	struct opt_count_t *opt = init_opt_count();
	int pos = 0;
	while (pos < argc) {
		if (!strcmp(argv[pos], "-o")) {
			check_str(argc - pos, argv + pos);
			opt->out_dir = argv[pos + 1];
			pos += 2;
		} else if (!strcmp(argv[pos], "-t")) {
			check_num(argc - pos, argv + pos);
			opt->n_threads = atoi(argv[pos + 1]);
			pos += 2;
		} else if (argv[pos][0] == '@') {
			parse_input_meta(opt, argv[pos] + 1);
			++pos;
		} else {
			__ERROR("Invalid option %s", argv[pos]);
		}
	}

	make_dir(opt->out_dir);
	return opt;
}
