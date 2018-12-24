#ifndef _BWT_H_
#define _BWT_H_

#include <stdint.h>
#include "attribute.h"

#ifdef HERA_64_BIT
#define OCC_INTV_SHIFT		7
#define OCC_INTV		0x80	// = (1 << OCC_INTV_SHIFT)
#define OCC_INTV_MASK		0x7f	// = (OCC_INTV - 1)
#define WORD_BWT_INTV		8	// = (OCC_INTV / 16)
#define WORD_OCC_INTV		8	// = sizeof(bioint_t)
#define WORD_OCC_BWT_INTV	16	// = WORD_BWT_INV + WORD_OCC_INTV
#define SHIFT_OCC_BWT_INTV	4	// = log2(WORD_OCC_BWT_INTV)
#define SA_INTV_SHIFT		5
#define SA_INTV			0x20	// = (1 << SA_INTV_SHIFT)
#define SA_INTV_MASK		0x1f	// = (SA_INTV - 1)
#define MAX_PAC			1000111000111000111ull
#else
#define OCC_INTV_SHIFT		6
#define OCC_INTV		0x40	// = (1 << OCC_INTV_SHIFT)
#define OCC_INTV_MASK		0x3f	// = (OCC_INTV - 1)
#define WORD_BWT_INTV		4	// = (OCC_INTV / 16)
#define WORD_OCC_INTV		4	// = sizeof(bioint_t)
#define WORD_OCC_BWT_INTV	8	// = WORD_BWT_INTV + WORD_OCC_INTV
#define SHIFT_OCC_BWT_INTV	3	// = log2(WORD_OCC_BWT_INTV)
#define SA_INTV_SHIFT		4
#define SA_INTV			0x10	// = (1 << SA_INTV_SHIFT)
#define SA_INTV_MASK		0xf	// = (SA_INTV - 1)
#define MAX_PAC			4000111000u
#endif

/* To calculate actually offset (in 32-bit block) of (k)-th character of compact bwt (bwt + occ checkpoint)
 * 1, Mem for Occ before k (in 8-bit block):
 *                    (k / OCC_INTV + 1) * (4 * sizeof(bioint_t))
 *    Since sizeof(uint32_t) = 32-bit, mem for Occ before k in 32-bit block:
 *                    (k >> OCC_INTV_SHIFT) * sizeof(bioint_t) + sizeof(bioint_t)
 * 2, Mem for BWT before k (in 32-bit block):
 *                    (k / OCC_INTV) * (OCC_INTV / 16) + (k % OCC_INTV) / 16
 *              =     (k >> OCC_INTV_SHIFT) + (k & OCC_INTV_MASK) >> 4
 * S = (k >> OCC_INTV_SHIFT) * (sizeof(bioint_t) + 1) 
 */

#define bwt_pos(b, k) ((b)[((k) >> OCC_INTV_SHIFT << SHIFT_OCC_BWT_INTV) +     \
				WORD_OCC_INTV + (((k) & OCC_INTV_MASK) >> 4)])
#define bwt_occ_intv(b, k) ((b) + ((k) >> OCC_INTV_SHIFT << SHIFT_OCC_BWT_INTV))
#define bwt_char(b, k) (bwt_pos(b, k) >> ((~(k) & 0xf) << 1) & 3)

#define __get_pac(pac, l) ((pac)[(l) >> 2] >> ((~(l) & 3) << 1) & 3)
#define __set_pac(pac, l, c) ((pac)[(l) >> 2] |= (c) << ((~(l) & 3) << 1))

struct bwt_t {
	// Genome sequence
	uint8_t *pac;
	bioint_t seq_len;

	// BWT
	bioint_t primary;
	bioint_t CC[5];
	bioint_t bwt_size;
	uint32_t *bwt;

	// SA
	bioint_t n_sa;
	bioint_t *sa;
};

struct bwt_t *bwt_build_from_fasta(const char *path);

bioint_t bwt_match_exact(struct bwt_t *bwt, const char *str, int len, bioint_t *sa_beg, bioint_t *sa_end);

void bwt_2occ(struct bwt_t *bwt, bioint_t l, bioint_t r, uint8_t c, bioint_t *o_l, bioint_t *o_r);

bioint_t bwt_sa(struct bwt_t *bwt, bioint_t k);

void bwt_dump(const char *path, struct bwt_t *bwt);

void bwt_load(const char *path, struct bwt_t *bwt);

void bwt_destroy(struct bwt_t *p);

#endif
