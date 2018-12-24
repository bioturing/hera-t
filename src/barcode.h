#ifndef _BARCODE_H_
#define _BARCODE_H_

#include "attribute.h"
#include "kmhash.h"
#include "opt.h"
#include "utils.h"

void init_barcode(struct gene_info_t *g);

// void add_bc_umi(int64_t idx, int32_t gene);

// void quantification(const char *out_dir, int n_thread);

void quantification(struct opt_count_t *opt, struct kmhash_t *h);

#endif
