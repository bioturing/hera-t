#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "bwt.h"
#include "divsufsort64.h"
#include "io_utils.h"
#include "index.h"
#include "kseq.h"
#include "utils.h"
#include "verbose.h"
#include "zlib.h"

KSEQ_INIT(gzFile, gzread)

static inline int occ_32_way(uint64_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = (y ^ (uint64_t)(((c & 2) >> 1) - 1)) >> 1 & (y ^ (uint64_t)((c & 1) - 1)) & 0x5555555555555555ull;
	// count the number of 1s in y
	y = (y & 0x3333333333333333ull) + (y >> 2 & 0x3333333333333333ull);
	return ((y + (y >> 4)) & 0xf0f0f0f0f0f0f0full) * 0x101010101010101ull >> 56;
}

static inline int occ_16_way(uint32_t y, int c)
{
	// reduce nucleotide counting to bits counting
	y = (y ^ (uint32_t)(((c & 2) >> 1) - 1)) >> 1 & (y ^ (uint32_t)((c & 1) - 1)) & 0x55555555u;
	// count the number of 1s in y
	y = (y & 0x33333333u) + (y >> 2 & 0x33333333u);
	return ((y + (y >> 4)) & 0xf0f0f0fu) * 0x1010101u >> 24;
}

bioint_t bwt_occ(struct bwt_t *bwt, bioint_t k, uint8_t c)
{
	if (k == bwt->seq_len) return bwt->CC[c + 1] - bwt->CC[c];
	if (k == (bioint_t)(-1)) return 0;

	bioint_t ret;
	uint32_t *p, *end;
	k -= (k >= bwt->primary);

	// get occ value at checkpoint
	ret = ((bioint_t *)(p = bwt_occ_intv(bwt->bwt, k)))[c];
	p += sizeof(bioint_t);

	// get occ value up to 32-bit block
	end = p + ((k >> 4) - ((k & (~OCC_INTV_MASK)) >> 4));
	for (; p < end; ++p)
		ret += occ_16_way(*p, c);

	// get occ value of last partial 32-bit block
	ret += occ_16_way((*p) & ~((1u << (((~k) & 0xf) << 1)) - 1), c);
	ret -= (bioint_t)(c == 0) * ((~k) & 0xf); // Remove count padded bit

	return ret;
}

void bwt_2occ(struct bwt_t *bwt, bioint_t l, bioint_t r, uint8_t c, bioint_t *o_l, bioint_t *o_r)
{
	bioint_t _l, _r;
	_l = (l >= bwt->primary) ? l - 1 : l;
	_r = (r >= bwt->primary) ? r - 1 : r;
	if ((_l >> OCC_INTV_SHIFT) != (_r >> OCC_INTV_SHIFT) || l == (bioint_t)(-1) || r == (bioint_t)(-1)) {
		*o_l = bwt_occ(bwt, l, c);
		*o_r = bwt_occ(bwt, r, c);
	} else {
		bioint_t m, n, i, j;
		uint32_t *p;
		n = ((bioint_t *)(p = bwt_occ_intv(bwt->bwt, _l)))[c];
		p += sizeof(bioint_t);
		j = _l >> 4 << 4;
		for (i = (_l >> OCC_INTV_SHIFT << OCC_INTV_SHIFT); i < j; i += 0x10, ++p)
			n += occ_16_way(*p, c);
		m = n;
		n += occ_16_way((*p) & ~((1u << (((~_l) & 0xf) << 1)) - 1), c);
		n -= (bioint_t)(c == 0) * ((~_l) & 0xf);
		*o_l = n;

		j = _r >> 4 << 4;
		for (; i < j; i += 0x10, ++p)
			m += occ_16_way(*p, c);
		m += occ_16_way((*p) & ~((1u << (((~_r) & 0xf) << 1)) - 1), c);
		m -= (bioint_t)(c == 0) * ((~_r) & 0xf);
		*o_r = m;
	}
}

static inline bioint_t bwt_invPsi(struct bwt_t *bwt, bioint_t k) // compute inverse CSA
{
	bioint_t x = k - (k > bwt->primary);
	x = bwt_char(bwt->bwt, x);
	x = bwt->CC[x] + bwt_occ(bwt, k, x);
	return k == bwt->primary ? 0 : x;
}

void bwt_cal_sa(struct bwt_t *bwt)
{
	bioint_t i, isa, sa;
	bwt->n_sa = (bwt->seq_len + SA_INTV) / SA_INTV;
	bwt->sa = (bioint_t *)calloc(bwt->n_sa, sizeof(bioint_t));

	isa = 0; sa = bwt->seq_len;
	for (i = 0; i < bwt->seq_len; ++i) {
		if (isa % SA_INTV == 0) bwt->sa[isa / SA_INTV] = sa;
		--sa;
		isa = bwt_invPsi(bwt, isa);
	}
	if (isa % SA_INTV == 0) bwt->sa[isa / SA_INTV] = sa;
	bwt->sa[0] = (bioint_t)(-1);
}

void bwt_construct(struct bwt_t *bwt)
{
	/* Construct bwt when bwt->pac is preloaded */
	if (!bwt || !bwt->pac)
		__ERROR("BWT was not allocated or pac was not preloaded");
	bioint_t i, k, l, n_occ, seq_len;
	uint8_t *pac, *buf;

	seq_len = bwt->seq_len;
	pac = bwt->pac;
	buf = calloc(bwt->seq_len + 1, 1);
	memset(bwt->CC, 0, 5 * sizeof(bioint_t));
	for (i = 0; i < seq_len; ++i) {
		buf[i] = pac[i >> 2] >> (((~i) & 3) << 1) & 3;
		++bwt->CC[1 + buf[i]];
	}
	for (i = 2; i <= 4; ++i)
		bwt->CC[i] += bwt->CC[i - 1];

	bwt->primary = (bioint_t)divbwt64(buf, buf, NULL, (saidx64_t)seq_len);
	bwt->bwt_size = (seq_len + 15) >> 4;
	n_occ = (seq_len + OCC_INTV - 1) / OCC_INTV + 1;
	bwt->bwt_size += n_occ * sizeof(bioint_t);
	bwt->bwt = calloc(bwt->bwt_size, 4);
	bioint_t c[4];
	c[0] = c[1] = c[2] = c[3] = 0;
	for (i = k = l = 0; i < seq_len; ++i) {
		if ((i & OCC_INTV_MASK) == 0) {
			memcpy(bwt->bwt + k, c, sizeof(bioint_t) * 4);
			/* literally, this line should be 
			 * k += sizeof(bioint_t) * 4 / sizeof(uint32_t)
			 * since sizeof(uint32_t) = 4, we can reduce the fraction
			 */
			k += sizeof(bioint_t); 
		}
		if ((i & 15) == 0)
			l = k++;
		bwt->bwt[l] |= (uint32_t)buf[i] << (((~i) & 15) << 1);
		++c[buf[i]];
	}
	memcpy(bwt->bwt + k, c, sizeof(bioint_t) * 4);
	assert(k + sizeof(bioint_t) == bwt->bwt_size);
	free(buf);
}

bioint_t bwt_sa(struct bwt_t *bwt, bioint_t k)
{
	bioint_t sa = 0;
	while (k & SA_INTV_MASK) {
		++sa;
		k = bwt_invPsi(bwt, k);
	}
	return sa + bwt->sa[k / SA_INTV];
}

bioint_t bwt_match_exact(struct bwt_t *bwt, const char *str, int len, bioint_t *sa_beg, bioint_t *sa_end)
{
	bioint_t l, r, o_l, o_r;
	int i;
	uint8_t c;
	l = 0; r = bwt->seq_len;
	for (i = len - 1; i >= 0; --i) {
		c = nt4_table[(int)str[i]];
		if (c > 3)
			return 0;
		bwt_2occ(bwt, l - 1, r, c, &o_l, &o_r);
		l = bwt->CC[c] + o_l + 1;
		r = bwt->CC[c] + o_r;
		if (l > r)
			return 0;
	}
	if (sa_beg) *sa_beg = l;
	if (sa_end) *sa_end = r;
	return (r - l + 1);
}

uint8_t *add_seq(kseq_t *seq, uint8_t *pac, bioint_t *l_pac, uint64_t *m_pac)
{
	int i, c;
	bioint_t l = *l_pac;
	uint64_t m = *m_pac;
	for (i = 0; i < (int)seq->seq.l; ++i) {
		c = nt4_table[(int)seq->seq.s[i]];
		if (c >= 4)
#ifndef _MSC_VER
			c = lrand48() & 3;
#else
			c = rand() & 3;
#endif
		if (l == m) {
			m <<= 1;
			pac = realloc(pac, m / 4);
			memset(pac + l / 4, 0, (m - l) / 4);
		}
		__set_pac(pac, l, c);
		++l;
	}
	*m_pac = m;
	*l_pac = l;
	return pac;
}

struct bwt_t *bwt_build_from_fasta(const char *path)
{
	struct bwt_t *bwt;
	bwt = calloc(1, sizeof(struct bwt_t));
	gzFile fp = gzopen(path, "r");
	kseq_t *seq;
	uint64_t m_pac;
	seq = kseq_init(fp);
#ifndef _MSC_VER
	srand48(11);
#else
	srand(11);
#endif
	m_pac = 0x10000;
	bwt->pac = calloc(m_pac >> 2, 1);
	bwt->seq_len = 0;
	while (kseq_read(seq) >= 0) {
		if (MAX_PAC - (bioint_t)seq->seq.l < bwt->seq_len) {
			if (sizeof(bioint_t) == 4)
				__ERROR("Genome length exceeds 4Bbp. Use Hera 64bit instead");
			else
				__ERROR("Genome length is too big!!! (or something naughty happened)");
		}
		bwt->pac = add_seq(seq, bwt->pac, &(bwt->seq_len), &m_pac);
	}
	m_pac = (bwt->seq_len >> 2) + ((bwt->seq_len & 3) == 0 ? 0 : 1);
	bwt->pac = realloc(bwt->pac, m_pac);
	bwt_construct(bwt);
	bwt_cal_sa(bwt);
	kseq_destroy(seq);
	gzclose(fp);
	return bwt;
}

void bwt_dump(const char *path, struct bwt_t *bwt)
{
	FILE *fp;
	fp = xfopen(path, "wb");
	// packed fasta sequences
	bioint_t pac_len;
	xfwrite(&bwt->seq_len, sizeof(bioint_t), 1, fp);
	pac_len = (bwt->seq_len >> 2) + ((bwt->seq_len & 3) == 0 ? 0 : 1);
	xfwrite(bwt->pac, 1, pac_len, fp);
	// BWT
	xfwrite(&bwt->primary, sizeof(bioint_t), 1, fp);
	xfwrite(bwt->CC + 1, sizeof(bioint_t), 4, fp);
	xfwrite(&bwt->bwt_size, sizeof(bioint_t), 1, fp);
	xfwrite(bwt->bwt, 4, bwt->bwt_size, fp);
	// SA
	xfwrite(&bwt->n_sa, sizeof(bioint_t), 1, fp);
	xfwrite(bwt->sa, sizeof(bioint_t), bwt->n_sa, fp);
	xwfclose(fp);
}

void bwt_load(const char *path, struct bwt_t *bwt)
{
	FILE *fp;
	fp = xfopen(path, "rb");
	// packed fasta sequences
	//__VERBOSE("[DEBUG] Reading fasta pack\n");
	bioint_t pac_len;
	xfread(&bwt->seq_len, sizeof(bioint_t), 1, fp);
	pac_len = (bwt->seq_len >> 2) + ((bwt->seq_len & 3) == 0 ? 0 : 1);
	bwt->pac = malloc(pac_len);
	xfread(bwt->pac, 1, pac_len, fp);
	// BWT
	//__VERBOSE("[DEBUG] Reading bwt binary\n");
	//__VERBOSE("[DEBUG] sizeof(bioint_t) = %d\n", (int)sizeof(bioint_t));
	xfread(&bwt->primary, sizeof(bioint_t), 1, fp);
	bwt->CC[0] = 0;
	xfread(bwt->CC + 1, sizeof(bioint_t), 4, fp);
	xfread(&bwt->bwt_size, sizeof(bioint_t), 1, fp);
	//__VERBOSE("[DEBUG] bwt_size = %u\n", bwt->bwt_size);
	bwt->bwt = malloc(bwt->bwt_size * 4);
	//if (bwt->bwt == NULL)
		//__VERBOSE("FAILLLLLLLLLLLL");
	xfread(bwt->bwt, 4, bwt->bwt_size, fp);
	// SA
	//__VERBOSE("[DEBUG] Reading suffix array\n");
	xfread(&bwt->n_sa, sizeof(bioint_t), 1, fp);
	bwt->sa = malloc(bwt->n_sa * sizeof(bioint_t));
	xfread(bwt->sa, sizeof(bioint_t), bwt->n_sa, fp);
	// TODO: check for integrity
	//__VERBOSE("[DEBUG] Done reading bwt\n");
	fclose(fp);
}

void bwt_destroy(struct bwt_t *p)
{
	if (!p) return;
	free(p->pac);
	free(p->bwt);
	free(p->sa);
	// free(p);
}
