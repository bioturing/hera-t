#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include "attribute.h"
#include "utils.h"

void alignment_init_hash(const char *path);

void alignment_init_ref_info(struct gene_info_t *g, struct transcript_info_t *t);

void align_chromium_read(struct read_t *read1, struct read_t *read2,
			 struct worker_bundle_t *bundle);

#endif
