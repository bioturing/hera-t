#ifndef _ATTRIBUTE_H_
#define _ATTRIBUTE_H_

#include <stdint.h>

#define READ_BLOCK		16
#define INDEX_VERSION		1
#define PROG_VERSION_MAJOR	0
#define PROG_VERSION_MINOR	2
#define PROG_VERSION_FIX	0

#define SIZE_1MB		1048576
#define SIZE_2MB		2097152
#define SIZE_4MB		4194304
#define SIZE_16MB		16777216
#define SIZE_128MB		134217728

#define GENE_BIT_LEN		22
#define GENE_MASK		UINT64_C(0x3fffff)

#define NNU			4

#define MAX_PATH		4096

#ifdef HERA_64_BIT
typedef uint64_t bioint_t;
#else
typedef uint32_t bioint_t;
#endif

typedef int64_t r1idx_t;

struct gene_info_t {
	int n;				// Number of genes
	int l_name;			// Gene name block length
	int l_id;			// Gene id block length
	int *chr_idx;			// Chromosome index of each gene
	char *gene_name;		// Genes name
	char *gene_id;			// Genes id
	char *strand;			// Genes strand (0 - forward, 1 - reverse)
};

struct genome_info_t {
	int n;				// Number of chromosomes
	int l_name;			// Chromosome name block length
	bioint_t *chr_len;		// Chromosomes length
	char *chr_name;			// Chromosomes name
	char **seq;			// Chromosomes sequence
};

struct exon_t {
	bioint_t beg;
	bioint_t end;
};

struct transcript_info_t {
	int n;				// Number of transcripts
	int l_id;			// Transcript id block length
	char *tran_id;			// Transcripts id
	int *tran_len;			// Transcripts length
	int *tran_beg;			// Begin of transcripts on [seq] array
	int *gene_idx;			// Gene index of each transcript

	int *idx;			// Which transcript index does one base belong to?
	char *seq;			// Concated sequence of all transcript

	int *n_exon;			// Number of exon of transcripts
	struct exon_t **exons;		// Exon info of transcripts
};

struct read_t {
	char *seq;			// Sequence
	char *name;			// Read name, omit '/1', '/2' and additional barcode/info
	char *qual;			// Base quality
	char *note;			// Comment
	char *rseq;			// Reverse complement of sequence
	char *rqual;			// Reverse string of base quality
	char *info;			// Additional barcode/info in read name
	int len;			// Read length
};

#endif
