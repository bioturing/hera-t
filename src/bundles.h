#ifndef _BUNDLES_H_
#define _BUNDLES_H_

#include "attribute.h"
#include "utils.h"
#include "interval_tree.h"

struct array_2D_t {
	void *data;
	int len;
	int nrow;
	void **rows;
};

struct bundle_data_t {
	struct raw_alg_t *alg_array;
	struct interval_t *intron_array;
	struct recycle_bin_t *recycle_bin;
	struct seed_t *seed_cons;
	struct array_2D_t *tmp_array;
};

void init_bundle(struct bundle_data_t *bundle);

void reinit_bundle(struct bundle_data_t *bundle);

void destroy_bundle(struct bundle_data_t *bundle);

#endif