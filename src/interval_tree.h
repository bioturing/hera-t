#ifndef _INTERVAL_TREE_
#define _INTERVAL_TREE_

#include "utils.h"

struct node_t {
        uint32_t start;
        uint32_t end;
        uint32_t range[2];
        uint32_t child[2];
        int center;
};

struct leaf_t {
        uint32_t start;
        uint32_t end;
};

struct tree_t {
        struct node_t *node;
        struct leaf_t *leaf;
        uint32_t n_leaf;
        uint32_t n_node;
};

struct interval_t {
        int n;
        int m;
        int *id;
};

void build_interval_tree(const char *out_path, struct leaf_t *leaf, int len);

void load_interval_tree(const char *tree_path, char *strand);

void query_interval(uint32_t qstart, uint32_t qend, int pos,
		    struct interval_t *res, char r_str);

#endif
