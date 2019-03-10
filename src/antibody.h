#ifndef _CELL_HASHING_H_
#define _CELL_HASHING_H_

#include "kmhash.h"
#include "khash.h"
#include "utils.h"

// antibody used for
#define CELL_HASHING            0
#define PROTEIN_QUANT           1   

KHASH_MAP_INIT_INT(bc, int);

struct reference_t {
	khash_t(bc) *h;

	int n_ref;
	int ref_len;
	int *name_idx;
	char *ref_name;
	char *ref_id;
};

struct antibody_lib_t {
        char *whitelist_path;   // Path to whitelist tsv file
	char **left_file;       // Path to read 1 fastqs
	char **right_file;      // Path to read 2 fastqs
	int n_files;

        int8_t type;            // antibody used for
        int8_t trim;           // Number of base must be trimmed from start

        struct reference_t *ref;
};

struct antibody_lib_t *init_antibody_lib(int8_t type);

struct antibody_lib_t *check_valid_cell(struct antibody_lib_t *lib);

struct antibody_lib_t *check_valid_protein(struct antibody_lib_t *lib);

void build_reference(struct antibody_lib_t *lib);

void map_antibody_read(struct read_t *read1, struct read_t *read2,
			 struct worker_bundle_t *bundle);

void print_ref(const char *out_dir, struct reference_t *ref);

#endif