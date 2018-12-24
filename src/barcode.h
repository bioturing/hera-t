#ifndef _BARCODE_H_
#define _BARCODE_H_

#include "attribute.h"
#include "kmhash.h"
#include "opt.h"
#include "utils.h"

void init_barcode(struct gene_info_t *g, struct library_t lib);

void quantification(struct opt_count_t *opt, struct kmhash_t *h);

#endif
