#ifndef _GENOME_H_
#define _GENOME_H_

#include "attribute.h"
#include "bwt.h"

void genome_init_bwt(struct bwt_t *b, int32_t count_intron);

int genome_map_err(struct read_t *read, int max_err,
		   struct worker_bundle_t *bundle);

#endif
