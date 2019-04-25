#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include "bundles.h"
#include "library_type.h"

void alignment_init_hash(const char *path);

void alignment_init_ref_info(struct gene_info_t *g, struct transcript_info_t *t);

int align_chromium_read(struct read_t *read, struct bundle_data_t *bundle,
			struct align_stat_t *count);

#endif
