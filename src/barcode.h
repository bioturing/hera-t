#ifndef _BARCODE_H_
#define _BARCODE_H_

#include "attribute.h"
#include "kmhash.h"
#include "opt.h"

struct ref_info_t {
	int n_refs;
	char *ref_text;
	int *ref_iter;
	char *gene_id;
	int *gene_iter;
	int type[4];
};

struct ref_info_t *init_ref_info();

void quantification(struct opt_count_t *opt, struct kmhash_t *h,
			struct ref_info_t *ref);

#endif
