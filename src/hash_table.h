#ifndef _HASH_TABLE_H_
#define _HASH_TABLE_H_

#include <stdint.h>

struct cons_bucket_t {
	uint32_t id;
	int head;
};

struct cons_build_t {
	struct cons_bucket_t **buckets;
	int *bsize;
	int *pos;
	int npos;
	int l2_size;
	uint32_t mask;
};

struct cons_table_t {
	uint32_t *id;
	int *head;
	int *bpos;
	int *pos;
	int l2_size;
	uint32_t mask;
};

/* CONS HASH */
void init_cons_hash(int size);
int insert_cons_hash(uint64_t id);
void addpos_cons_hash(uint64_t id, int pos);
void recount_cons_hash();
void store_cons_hash(const char *file_path, int kcons);
void load_cons_hash(const char *file_path, int *kcons);
int query_cons_hash(uint64_t id, int **pos);
void free_cons_hash_index();
void free_cons_hash();

#endif /* _HASH_TABLE_H_ */
