#include <stdint.h>
#include <stdio.h>
#include "interval_tree.h"

#define intersect(s, e, qs, qe) ((qs >= s && qe <= qe)? 1: 0)

static struct tree_t *tree;
static char *gene_strand;

void init_tree(struct leaf_t *leaf, int n_leaf)
{
	tree = malloc(sizeof(struct tree_t));

	tree->node = calloc(1, sizeof(struct node_t));
	tree->n_node = 0;

	tree->leaf = leaf;
	tree->n_leaf = n_leaf;
}

void build_tree(int start, int end, int pos)
{
	int i;
	uint32_t mid, center, max;
	
	i = 0;
	mid = (start + end) >> 1;
	center = (tree->leaf[mid].start + tree->leaf[mid].end) >> 1;
	max = tree->leaf[mid].end;

	tree->node[pos].center = center;
	if (start == end){
		tree->node[pos].start = tree->leaf[start].start;
		tree->node[pos].end = tree->leaf[start].end;
		tree->node[pos].range[0] = tree->node[pos].range[1] = start;
		tree->node[pos].child[0] = tree->node[pos].child[1] = 0;
		return;
	}

	for (i = mid + 1; i > start; --i){
		if (tree->leaf[i - 1].end < center)
			break;

		max = __max(max, tree->leaf[i - 1].end);
	}
	tree->node[pos].start = tree->leaf[i].start;
	tree->node[pos].range[0] = i;

	for (i = mid + 1; i <= end; ++i){
		if (tree->leaf[i].start > center)
			break;

		max = __max(max, tree->leaf[i - 1].end);
	}
	tree->node[pos].end = max;
	tree->node[pos].range[1] = i - 1;

	if (tree->node[pos].range[0] > start){
		++tree->n_node;
		tree->node = realloc(tree->node,
				     tree->n_node * sizeof(struct node_t));
		tree->node[pos].child[0] = tree->n_node - 1;
		build_tree(start, tree->node[pos].range[0] - 1, tree->n_node - 1);
	} else {
		tree->node[pos].child[0] = 0;
	}

	if (tree->node[pos].range[1] < end){
		++tree->n_node;
		tree->node = realloc(tree->node,
					tree->n_node * sizeof(struct node_t));
		tree->node[pos].child[1] = tree->n_node - 1;
		build_tree(tree->node[pos].range[1] + 1, end, tree->n_node - 1);
	} else {
		tree->node[pos].child[1] = 0;
	}

	return;
}

void write_tree(const char *path)
{
	FILE *f = fopen(path, "wb");

	fwrite(&tree->n_node, sizeof(int), 1, f);
	fwrite(&tree->n_leaf, sizeof(int), 1, f);
	fwrite(tree->node, sizeof(struct node_t), tree->n_node, f);
	fwrite(tree->leaf, sizeof(struct leaf_t), tree->n_leaf, f);

	fclose(f);
}

void add_result(struct node_t node, struct leaf_t *leaf, uint32_t qstart,
		uint32_t qend, struct interval_t *res, char r_str)
{
	int i;
	for (i = node.range[0]; i <= node.range[1]; ++i){
		if (r_str == gene_strand[i] && 
		    intersect(leaf[i].start, leaf[i].end, qstart, qend)){
			if (res->n == res->m){
				res->m += 64;
				res->id = realloc(res->id, res->m * sizeof(int));
			}

			res->id[res->n] = i;
			++res->n;
		} 
	}	
}

void query_interval(uint32_t qstart, uint32_t qend, int pos,
		    struct interval_t *res, char r_str)
{
	int s, e, next;
	struct node_t node;

	node = tree->node[pos];
	s = node.start;
	e = node.end;

	if (intersect(s, e, qstart, qend))
		add_result(node, tree->leaf, qstart, qend, res, r_str);

	next = node.child[0];
	if (next > 0 && qend < node.center)
		query_interval(qstart, qend, next, res, r_str);

	next = node.child[1];
	if (next > 0 && qstart > node.center)
		query_interval(qstart, qend, next, res, r_str);
}

void build_interval_tree(const char *out_path, struct leaf_t *leaf, int len)
{
	init_tree(leaf, len);
	build_tree(0, tree->n_leaf - 1, 0);
	write_tree(out_path);
	free(tree->node);
	free(tree);
}

void load_interval_tree(const char *tree_path, char *strand)
{
	FILE *f = fopen(tree_path, "rb");
	int ret;

	tree = malloc(sizeof(struct tree_t));
	ret = fread(&tree->n_node, sizeof(int), 1, f);
	ret = fread(&tree->n_leaf, sizeof(int), 1, f);

	tree->node = malloc(tree->n_node * sizeof(struct node_t));
	tree->leaf = malloc(tree->n_leaf * sizeof(struct leaf_t));

	ret = fread(tree->node, sizeof(struct node_t), tree->n_node, f);
	ret = fread(tree->leaf, sizeof(struct leaf_t), tree->n_leaf, f);
	
	fclose(f);

	gene_strand = strand;
}