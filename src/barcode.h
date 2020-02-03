#ifndef _BARCODE_H_
#define _BARCODE_H_

#include "attribute.h"
#include "bc_hash.h"
#include "opt.h"

struct ref_info_t {
	int n_refs;
	int n_rna;
	char *ref_text;
	int *text_iter;
	char *ref_id;
	int *id_iter;
	int text_len;
	int id_len;
	int type[4];
};

struct ref_info_t *init_ref_info();

void quantification(struct opt_count_t *opt, struct bc_hash_t *h,
			struct ref_info_t *ref);

#endif
