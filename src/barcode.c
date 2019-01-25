#include <pthread.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#include "barcode.h"
#include "io_utils.h"
#include "kmhash.h"
#include "verbose.h"
#include "utils.h"

#define CUT_OFF_THRES			0.01

#define __get_gene(x) ((int)((x) & GENE_MASK))
#define __get_umi(x) ((x) >> GENE_BIT_LEN)

struct sc_cell_t {
	uint64_t idx;
	uint32_t cnt_umi;
	struct umi_hash_t *h;
};

struct umi_bundle_t {
	int n_threads;
	int thread_no;
	pthread_mutex_t *lock;
	FILE *fbin;
	uint64_t *n_lines;
};

const uint64_t pow5_r[] = {1ull, 5ull, 25ull, 125ull, 625ull, 3125ull, 15625ull, 78125ull, 390625ull, 1953125ull, 9765625ull, 48828125ull, 244140625ull, 1220703125ull, 6103515625ull, 30517578125ull, 152587890625ull, 762939453125ull, 3814697265625ull, 19073486328125ull, 95367431640625ull, 476837158203125ull, 2384185791015625ull, 11920928955078125ull, 59604644775390625ull, 298023223876953125ull, 1490116119384765625ull};

struct sc_cell_t *CBs;
int32_t n_bc;
int16_t bc_len, umi_len;
struct gene_info_t genes;

void init_barcode(struct gene_info_t *g, struct library_t lib)
{
	extern struct gene_info_t genes;
	memcpy(&genes, g, sizeof(struct gene_info_t));

	bc_len = lib.bc_len;
	umi_len = lib.umi_len;
}

static void insertion_sort(struct sc_cell_t *b, struct sc_cell_t *e)
{
	struct sc_cell_t *i, *j, tmp;
	for (i = b + 1; i < e; ++i) {
		if (i->cnt_umi > (i - 1)->cnt_umi) {
			tmp = *i;
			for (j = i; j > b && tmp.cnt_umi > (j - 1)->cnt_umi; j--)
				*j = *(j - 1);
			*j = tmp;
		}
	}
}

void merge_sort(struct sc_cell_t *a, int l, int r, int m, struct sc_cell_t *tmp)
{
	struct sc_cell_t *a1, *a2;
	int len1, len2, i1, i2, k;
	len1 = m - l;
	len2 = r - m;
	memcpy(tmp, a + l, len1 * sizeof(struct sc_cell_t));

	a1 = tmp;
	a2 = a + m;
	a = a + l;

	i1 = i2 = k = 0;
	while (i1 < len1 && i2 < len2) {
		if (a1[i1].cnt_umi > a2[i2].cnt_umi)
			a[k++] = a1[i1++];
		else
			a[k++] = a2[i2++];
	}

	if (i1 < len1)
		memcpy(a + k, a1 + i1, (len1 - i1) * sizeof(struct sc_cell_t));

	if (i2 < len2)
		memcpy(a + k, a2 + i2, (len2 - i2) * sizeof(struct sc_cell_t));
}

void sort_CBs()
{
	extern int n_bc;
	extern struct sc_cell_t *CBs;

	int i, buck_size, l, r, m;
	struct sc_cell_t *tmp;

	tmp = malloc(n_bc * sizeof(struct sc_cell_t));

	buck_size = 64;
	for (i = 0; i < n_bc; i += buck_size)
		insertion_sort(CBs + i, CBs + __min(i + buck_size, n_bc));

	for (i = buck_size; i < n_bc; i <<= 1) {
		for (l = 0; l < n_bc; l += (i << 1)) {
			r = __min(l + (i << 1), n_bc);
			m = __min(l + i, r);
			merge_sort(CBs, l, r, m, tmp);
		}
	}
	free(tmp);
}

void cut_off_barcode(struct kmhash_t *h)
{
	extern int n_bc;
	extern struct sc_cell_t *CBs;
	kmint_t k;
	int i, j, flag, ch, c;
	uint64_t bc_idx, new_bc, tmp_idx;
	uint32_t cut_off;
	CBs = malloc(h->n_items * sizeof(struct sc_cell_t));
	n_bc = 0;
	for (k = 0; k < h->size; ++k) {
		if (h->bucks[k].idx == TOMB_STONE)
			continue;
		CBs[n_bc].idx = h->bucks[k].idx;
		CBs[n_bc].cnt_umi = h->bucks[k].umis->n_items;
		CBs[n_bc].h = h->bucks[k].umis;
		++n_bc;
	}
	
 	//__VERBOSE("Number of raw barcodes: %d\n", n_bc);

	sort_CBs();

	h->pos = malloc(h->size * sizeof(int));
	
	for (i = 0; i < n_bc; ++i) {
		k = kmhash_get(h, CBs[i].idx);
		h->pos[k] = i;
	}

	for (i = 0; i < n_bc; ++i) {
		bc_idx = CBs[i].idx;
		flag = 0;
		for (j = 0; j < bc_len; ++j) {
			ch = (bc_idx / pow5_r[j]) % 5;
			tmp_idx = bc_idx - pow5_r[j] * ch;
			for (c = 0; c < NNU; ++c) {
				if (ch == c)
					continue;
				new_bc = tmp_idx + pow5_r[j] * c;
				k = kmhash_get(h, new_bc);
				if (k < KMHASH_MAX_SIZE && h->bucks[k].umis->n_items >= CBs[i].cnt_umi) {
					flag = 1;
					break;
				}
			}
			if (flag)
				break;
		}
		if (flag) {
			n_bc = i;
			break;
		}
	}

	cut_off = CUT_OFF_THRES * CBs[0].cnt_umi;
	for (i = 0; i < n_bc; ++i) {
		if (CBs[i].cnt_umi < cut_off) {
			n_bc = i;
			break;
		}
	}
}

static void merge_umi(struct kmhash_t *h, uint64_t dst_bc, uint64_t src_bc)
{
	kmint_t dst_k, src_k, k;
	struct umi_hash_t *dst_umi, *src_umi;
	dst_k = kmhash_get(h, dst_bc);
	src_k = kmhash_get(h, src_bc);
	assert(dst_k != KMHASH_MAX_SIZE);
	assert(src_k != KMHASH_MAX_SIZE);
	dst_umi = h->bucks[dst_k].umis;
	src_umi = h->bucks[src_k].umis;
	for (k = 0; k < src_umi->size; ++k) {
		if (src_umi->bucks[k] == TOMB_STONE)
			continue;
		umihash_put_umi_single(dst_umi, src_umi->bucks[k]);
	}
	umihash_destroy(src_umi);
	h->bucks[src_k].umis = NULL;
}

void correct_barcode(struct kmhash_t *h)
{
	extern int n_bc;
	extern struct sc_cell_t *CBs;
	uint64_t bc_idx, cur_bc, new_bc, tmp_idx;
	int l, i, k, ch, c, cnt_new;
	kmint_t iter;
	l = (int)h->n_items;
	for (i = n_bc; i < l; ++i) {
		bc_idx = CBs[i].idx;
		cnt_new = 0;
		new_bc = 0;
		for (k = 0; k < bc_len; ++k) {
			ch = (bc_idx / pow5_r[k]) % 5;
			tmp_idx = bc_idx - ch * pow5_r[k];
			for (c = 0; c < NNU; ++c) {
				if (ch == c)
					continue;
				cur_bc = tmp_idx + pow5_r[k] * c;
				iter = kmhash_get(h, cur_bc);
				if (iter < KMHASH_MAX_SIZE && h->pos[iter] < n_bc) {
					++cnt_new;
					new_bc = cur_bc;
				}
			}
		}
		if (cnt_new == 1)
			merge_umi(h, new_bc, bc_idx);
	}
}

#ifdef _WIN32
void print_barcodes(const TCHAR *out_dir)
{
	TCHAR out_path[MAX_PATH];
	char *seq;
	FILE *fp;
	int i;

	_tcscpy(out_path, out_dir);
	_tcscat(out_path, "/barcodes.tsv");

	fp = xfopen(out_path, "w");

	for (i = 0; i < n_bc; ++i) {
		seq = num2seq(CBs[i].idx, bc_len);
		fprintf(fp, "%s\n", seq);
		free(seq);
	}

	fclose(fp);
}

void print_genes(const TCHAR *out_dir)
{
	TCHAR out_path[MAX_PATH];
	FILE *fp;
	int i;

	_tcscpy(out_path, out_dir);
	_tcscat(out_path, "/genes.tsv");

	fp = xfopen(out_path, "w");

	for (i = 0; i < genes.n; ++i)
		fprintf(fp, "%s\t%s\n", genes.gene_id + genes.l_id * i,
					genes.gene_name + genes.l_name * i);

	fclose(fp);
}

void write_mtx(const TCHAR *bin_path, const TCHAR *out_dir)
{
	TCHAR out_path[MAX_PATH];
	FILE *fbin, *fmtx;
	_tcscpy(out_path, out_dir);
	_tcscat(out_path, "/matrix.mtx");
	fbin = xfopen(bin_path, "rb");
	fmtx = xfopen(out_path, "wb");
	int m, n, r, c, cnt;
	uint64_t k, i;
	xfread(&n, sizeof(int), 1, fbin);
	xfread(&m, sizeof(int), 1, fbin);
	xfread(&k, sizeof(uint64_t), 1, fbin);
	fprintf(fmtx, "%%%%MatrixMarket matrix coordinate integer general\n");
	fprintf(fmtx, "%%Gene count matrix generated by hera-T version %d.%d.%d\n",
		PROG_VERSION_MAJOR, PROG_VERSION_MINOR, PROG_VERSION_FIX);
	fprintf(fmtx, "%d\t%d\t%llu\n", n, m, (unsigned long long)k);
	for (i = 0; i < k; ++i) {
		xfread(&r, sizeof(int), 1, fbin);
		xfread(&c, sizeof(int), 1, fbin);
		xfread(&cnt, sizeof(int), 1, fbin);
		fprintf(fmtx, "%d\t%d\t%d\n", r, c, cnt);
	}
	fclose(fbin);
	xwfclose(fmtx);
}
#else
void print_barcodes(const char *out_dir)
{
	char out_path[MAX_PATH];
	char *seq;
	FILE *fp;
	int i;

	strcpy(out_path, out_dir);
	strcat(out_path, "/barcodes.tsv");

	fp = xfopen(out_path, "w");

	for (i = 0; i < n_bc; ++i) {
		seq = num2seq(CBs[i].idx, bc_len);
		fprintf(fp, "%s\n", seq);
		free(seq);
	}

	fclose(fp);
}

void print_genes(const char *out_dir)
{
	char out_path[MAX_PATH];
	FILE *fp;
	int i;

	strcpy(out_path, out_dir);
	strcat(out_path, "/genes.tsv");

	fp = xfopen(out_path, "w");

	for (i = 0; i < genes.n; ++i)
		fprintf(fp, "%s\t%s\n", genes.gene_id + genes.l_id * i,
					genes.gene_name + genes.l_name * i);

	fclose(fp);
}

void write_mtx(const char *bin_path, const char *out_dir)
{
	char out_path[MAX_PATH];
	FILE *fbin, *fmtx;
	strcpy(out_path, out_dir);
	strcat(out_path, "/matrix.mtx");
	fbin = xfopen(bin_path, "rb");
	fmtx = xfopen(out_path, "wb");
	int m, n, r, c, cnt;
	uint64_t k, i;
	xfread(&n, sizeof(int), 1, fbin);
	xfread(&m, sizeof(int), 1, fbin);
	xfread(&k, sizeof(uint64_t), 1, fbin);
	fprintf(fmtx, "%%%%MatrixMarket matrix coordinate integer general\n");
	fprintf(fmtx, "%%Gene count matrix generated by hera-T version %d.%d.%d\n",
		PROG_VERSION_MAJOR, PROG_VERSION_MINOR, PROG_VERSION_FIX);
	fprintf(fmtx, "%d\t%d\t%llu\n", n, m, (unsigned long long)k);
	for (i = 0; i < k; ++i) {
		xfread(&r, sizeof(int), 1, fbin);
		xfread(&c, sizeof(int), 1, fbin);
		xfread(&cnt, sizeof(int), 1, fbin);
		fprintf(fmtx, "%d\t%d\t%d\n", r, c, cnt);
	}
	fclose(fbin);
	xwfclose(fmtx);
}
#endif

void correct_umi(struct sc_cell_t *bc)
{
	struct umi_hash_t *h;
	kmint_t i, iter;
	uint64_t umi_idx, new_umi, new_bin;
	int gene, has_Ns, k, ch, c, flag;
	h = bc->h;

	for (i = 0; i < h->size; ++i) {
		if (h->bucks[i] == TOMB_STONE || __get_gene(h->bucks[i]) == genes.n)
			continue;
		umi_idx = __get_umi(h->bucks[i]);
		gene = __get_gene(h->bucks[i]);
		has_Ns = 0;

		for (k = 0; k < umi_len; ++k) {
			ch = (umi_idx / pow5_r[k]) % 5;
			has_Ns += (ch == NNU);
		}

		if (has_Ns) {
			h->bucks[i] = h->bucks[i] >> GENE_BIT_LEN << GENE_BIT_LEN | genes.n;
			continue;
		}
		flag = 0;

		for (k = 0; k < umi_len; ++k) {
			ch = (umi_idx / pow5_r[k]) % 5;
			for (c = 0; c < NNU; ++c) {
				if (c == ch)
					continue;
				new_umi = umi_idx + (c - ch) * pow5_r[k];
				new_bin = (new_umi << GENE_BIT_LEN | gene);
				iter = umihash_get(h, new_bin);
				if (iter != KMHASH_MAX_SIZE) {
					h->bucks[i] = h->bucks[i] >> GENE_BIT_LEN << GENE_BIT_LEN | genes.n;
					flag = 1;
					break;
				}
			}
			if (flag)
				break;
		}
	}
}

void count_genes(struct sc_cell_t *bc, int bc_pos, int *cnt, FILE *fbin,
		 pthread_mutex_t *lock, uint64_t *n_lines)
{
	kmint_t i;
	int k;
	struct umi_hash_t *h;
	h = bc->h;

	memset(cnt, 0, genes.n * sizeof(int));
	for (i = 0; i < h->size; ++i) {
		if (h->bucks[i] == TOMB_STONE || __get_gene(h->bucks[i]) == genes.n)
			continue;
		++cnt[__get_gene(h->bucks[i])];
	}

	pthread_mutex_lock(lock);
	for (k = 1; k <= genes.n; ++k) {
		if (cnt[k - 1]) {
			++*n_lines;
			xfwrite(&k, sizeof(int), 1, fbin);
			xfwrite(&bc_pos, sizeof(int), 1, fbin);
			xfwrite(cnt + (k - 1), sizeof(int), 1, fbin);
		}
	}
	pthread_mutex_unlock(lock);
}

void *umi_worker(void *data)
{
	struct umi_bundle_t *bundle = (struct umi_bundle_t *)data;
	int i, n_threads, thread_no;
	FILE *fbin;
	pthread_mutex_t *lock;
	uint64_t *n_lines;
	int *cnt;

	cnt = malloc(genes.n * sizeof(int));
	n_threads = bundle->n_threads;
	thread_no = bundle->thread_no;
	fbin = bundle->fbin;
	lock = bundle->lock;
	n_lines = bundle->n_lines;

	for (i = 0; i < n_bc; ++i) {
		if (i % n_threads == thread_no) {
			correct_umi(CBs + i);
			count_genes(CBs + i, i + 1, cnt, fbin, lock, n_lines);
		}
	}

	free(cnt);
	pthread_exit(NULL);
}

#ifdef _WIN32
void quantification(struct opt_count_t *opt, struct kmhash_t *h)
{
	cut_off_barcode(h);
	__VERBOSE("Done cutting off barcode\n");
	correct_barcode(h);
	__VERBOSE("Done processing barcode\n");

	print_barcodes(opt->out_dir);
	print_genes(opt->out_dir);

	TCHAR out_path[MAX_PATH];
	_tcscpy(out_path, opt->out_dir);
	_tcscat(out_path, "/matrix.bin");

	FILE *fbin;
	fbin = xfopen(out_path, "wb");

	xfwrite(&genes.n, sizeof(int), 1, fbin);
	xfwrite(&n_bc, sizeof(int), 1, fbin);

	uint64_t n_lines = 0;
	xfwrite(&n_lines, sizeof(uint64_t), 1, fbin);

	pthread_attr_t attr;
	pthread_t *t;
	pthread_mutex_t *lock;
	struct umi_bundle_t *bundles;

	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	lock = calloc(1, sizeof(pthread_mutex_t));
	pthread_mutex_init(lock, NULL);

	t = calloc(opt->n_threads, sizeof(pthread_t));
	bundles = calloc(opt->n_threads, sizeof(struct umi_bundle_t));

	int i;
	for (i = 0; i < opt->n_threads; ++i) {
		bundles[i].lock = lock;
		bundles[i].fbin = fbin;
		bundles[i].thread_no = i;
		bundles[i].n_threads = opt->n_threads;
		bundles[i].n_lines = &n_lines;
		pthread_create(t + i, &attr, umi_worker, bundles + i);
	}

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(t[i], NULL);

	fseek(fbin, 8L, SEEK_SET);
	xfwrite(&n_lines, sizeof(uint64_t), 1, fbin);
	xwfclose(fbin);

	write_mtx(out_path, opt->out_dir);

	free(t);
	free(bundles);
}
#else
void quantification(struct opt_count_t *opt, struct kmhash_t *h)
{
	cut_off_barcode(h);
	__VERBOSE("Done cutting off barcode\n");
	correct_barcode(h);
	__VERBOSE("Done processing barcode\n");

	print_barcodes(opt->out_dir);
	print_genes(opt->out_dir);

	char out_path[MAX_PATH];
	strcpy(out_path, opt->out_dir);
	strcat(out_path, "/matrix.bin");

	FILE *fbin;
	fbin = xfopen(out_path, "wb");

	xfwrite(&genes.n, sizeof(int), 1, fbin);
	xfwrite(&n_bc, sizeof(int), 1, fbin);

	uint64_t n_lines = 0;
	xfwrite(&n_lines, sizeof(uint64_t), 1, fbin);

	pthread_attr_t attr;
	pthread_t *t;
	pthread_mutex_t *lock;
	struct umi_bundle_t *bundles;

	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	lock = calloc(1, sizeof(pthread_mutex_t));
	pthread_mutex_init(lock, NULL);

	t = calloc(opt->n_threads, sizeof(pthread_t));
	bundles = calloc(opt->n_threads, sizeof(struct umi_bundle_t));

	int i;
	for (i = 0; i < opt->n_threads; ++i) {
		bundles[i].lock = lock;
		bundles[i].fbin = fbin;
		bundles[i].thread_no = i;
		bundles[i].n_threads = opt->n_threads;
		bundles[i].n_lines = &n_lines;
		pthread_create(t + i, &attr, umi_worker, bundles + i);
	}

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(t[i], NULL);

	fseek(fbin, 8L, SEEK_SET);
	xfwrite(&n_lines, sizeof(uint64_t), 1, fbin);
	xwfclose(fbin);

	write_mtx(out_path, opt->out_dir);

	free(t);
	free(bundles);
}
#endif

