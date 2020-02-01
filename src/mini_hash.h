//
// Created by BioTuring on 2020-02-01.
//

#ifndef HERA_T_MINI_HASH_H
#define HERA_T_MINI_HASH_H

#include <stdint.h>

#define EMPTY_SLOT 0xffffffffffffffff
#define TOME_STONE 0x00000000ffffffff

#define __mini_empty(table, i) (table->key[i] == EMPTY_SLOT)
#define INIT_PRIME_INDEX 16

struct mini_hash_t {
    uint64_t *h;
    uint64_t *key;
    uint64_t size;
    uint64_t count;
    uint64_t max_cnt;
    int prime_index;
};

void init_mini_hash(struct mini_hash_t **h_table, uint32_t p_index);
void destroy_mini_hash(struct mini_hash_t *h_table);
uint64_t *mini_put(struct mini_hash_t **h_table, uint64_t data);
uint64_t *mini_put_by_key(struct mini_hash_t *h_table, uint64_t data, uint64_t key);
uint64_t *mini_get(struct mini_hash_t *h_table, uint64_t data);
uint64_t twang_mix64(uint64_t key);
#endif //HERA_T_MINI_HASH_H
