#ifndef _BARCODE_H_
#define _BARCODE_H_

#include "attribute.h"
#include "kmhash.h"
#include "opt.h"

void quantification(struct opt_count_t *opt, struct kmhash_t *h,
			const char *type, int n_refs);

void print_genes(const char *out_dir, const char *type, struct gene_info_t genes);

#endif
