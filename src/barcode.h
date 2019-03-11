#ifndef _BARCODE_H_
#define _BARCODE_H_

#include "attribute.h"
#include "kmhash.h"
#include "opt.h"

void init_barcode(struct gene_info_t *g, struct library_t lib);

void quantification(struct opt_count_t *opt, struct kmhash_t *h);

struct kmhash_t *build_hash_from_cutoff(int n_threads);

void antibody_quant(struct opt_count_t *opt, struct kmhash_t *h,
					struct antibody_lib_t *lib, const char *dir);

#endif
