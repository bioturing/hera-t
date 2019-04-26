#include <pthread.h>

#include "barcode.h"
#include "io_utils.h"
#include "verbose.h"
#include "utils.h"

#define __get_gene(x) ((int)((x) & GENE_MASK))
#define __get_umi(x) ((x) >> GENE_BIT_LEN)

const uint64_t pow5_r[] = {1ull, 5ull, 25ull, 125ull, 625ull, 3125ull, 15625ull, 78125ull, 390625ull, 1953125ull, 9765625ull, 48828125ull, 244140625ull, 1220703125ull, 6103515625ull, 30517578125ull, 152587890625ull, 762939453125ull, 3814697265625ull, 19073486328125ull, 95367431640625ull, 476837158203125ull, 2384185791015625ull, 11920928955078125ull, 59604644775390625ull, 298023223876953125ull, 1490116119384765625ull};

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
	int n_refs;
};

struct sc_cell_t *CBs;
int32_t n_bc;
int16_t bc_len, umi_len;

struct ref_info_t *init_ref_info()
{
	struct ref_info_t *ref = malloc(sizeof(struct ref_info_t));
	ref->n_refs = 0;
	ref->text_iter = calloc(1, sizeof(int));
	ref->ref_text = malloc(1);
	ref->ref_id = malloc(1);
	ref->id_iter = calloc(1, sizeof(int));
	ref->text_len = ref->id_len = 0;

	return ref;
}

void destroy_ref_info(struct ref_info_t *ref)
{
	free(ref->text_iter);
	free(ref->ref_text);
	free(ref->ref_id);
	free(ref->id_iter);
	free(ref);
}

static int compare_cell(const void *a, const void *b)
{
	struct sc_cell_t *c_a = (struct sc_cell_t *) a;
	struct sc_cell_t *c_b = (struct sc_cell_t *) b;
	return (c_b->cnt_umi - c_a->cnt_umi);
}

int *cut_off_barcode(struct kmhash_t *h)
{
	kmint_t k;
	int i, j, flag, ch, c;
	uint64_t bc_idx, new_bc, tmp_idx;

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

	qsort(CBs, n_bc, sizeof(struct sc_cell_t), compare_cell);

	int *index = malloc(h->size * sizeof(int));
	
	for (i = 0; i < n_bc; ++i) {
		k = kmhash_get(h, CBs[i].idx);
		index[k] = i;
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

	return index;
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
}

void correct_barcode(struct kmhash_t *h, int *index)
{
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
				if (iter < KMHASH_MAX_SIZE && index[iter] < n_bc) {
					++cnt_new;
					new_bc = cur_bc;
				}
			}
		}
		if (cnt_new == 1)
			merge_umi(h, new_bc, bc_idx);
	}
}

void print_barcodes(const char *out_path)
{
	char *seq;
	FILE *fp;
	int i;

	fp = xfopen(out_path, "w");

	for (i = 0; i < n_bc; ++i) {
		seq = num2seq(CBs[i].idx, bc_len);
		fprintf(fp, "%s\n", seq);
		free(seq);
	}

	fclose(fp);
}

void print_refs(struct ref_info_t *ref, char *out_path)
{
	FILE *fp;
	int i, start;

	fp = xfopen(out_path, "w");
	start = 0;

	// Print genes
	if (ref->type[0]) {
		for (i = start; i < ref->type[0]; ++i)
			fprintf(fp, "%s\t%s\tGene expression\n",
				ref->ref_id + ref->id_iter[i],
				ref->ref_text + ref->text_iter[i]);
		start = ref->type[0];
	}

	// Print cell hashing
	if (ref->type[1]) {
		for (i = start; i < ref->type[1]; ++i)
			fprintf(fp, "%s\t%s\tCell hashing\n",
				ref->ref_id + ref->id_iter[i],
				ref->ref_text + ref->text_iter[i]);

		start += ref->type[1];
	}

	// Print protein measurement
	if (ref->type[2]) {
		for (i = start; i < ref->type[2]; ++i)
			fprintf(fp, "%s\t%s\tProtein measurement\n",
				ref->ref_id + ref->id_iter[i],
				ref->ref_text + ref->text_iter[i]);

		start += ref->type[2];
	}

	// Print CRISPR capture
	if (ref->type[3]) {
		for (i = start; i < ref->type[3]; ++i)
			fprintf(fp, "%s\t%s\tCRISPR capture\n",
				ref->ref_id + ref->id_iter[i],
				ref->ref_text + ref->text_iter[i]);

		start += ref->type[3];
	}
	fclose(fp);
}

void write_mtx(char *path, int len)
{
	FILE *fbin, *fmtx;
	fbin = xfopen(path, "rb");

	concat_str(path, len, "/matrix.mtx", 11);
	fmtx = xfopen(path, "wb");

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

void correct_umi(struct sc_cell_t *bc, int n_refs)
{
	struct umi_hash_t *h;
	kmint_t i, iter;
	uint64_t umi_idx, new_umi, new_bin;
	int gene, has_Ns, k, ch, c, flag;
	h = bc->h;

	for (i = 0; i < h->size; ++i) {
		if (h->bucks[i] == TOMB_STONE || __get_gene(h->bucks[i]) == n_refs)
			continue;
		umi_idx = __get_umi(h->bucks[i]);
		gene = __get_gene(h->bucks[i]);
		has_Ns = 0;

		for (k = 0; k < umi_len; ++k) {
			ch = (umi_idx / pow5_r[k]) % 5;
			has_Ns += (ch == NNU);
		}

		if (has_Ns) {
			h->bucks[i] = h->bucks[i] >> GENE_BIT_LEN << GENE_BIT_LEN | n_refs;
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
					h->bucks[i] = h->bucks[i] >> GENE_BIT_LEN << GENE_BIT_LEN | n_refs;
					flag = 1;
					break;
				}
			}
			if (flag)
				break;
		}
	}
}

void count_genes(struct sc_cell_t *bc, int bc_pos, int *cnt, int n_refs,
		FILE *fbin, pthread_mutex_t *lock, uint64_t *n_lines)
{
	kmint_t i;
	int k;
	struct umi_hash_t *h;
	h = bc->h;

	memset(cnt, 0, n_refs * sizeof(int));
	for (i = 0; i < h->size; ++i) {
		if (h->bucks[i] == TOMB_STONE || __get_gene(h->bucks[i]) == n_refs)
			continue;
		++cnt[__get_gene(h->bucks[i])];
	}

	pthread_mutex_lock(lock);
	for (k = 1; k <= n_refs; ++k) {
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

	cnt = malloc(bundle->n_refs * sizeof(int));
	n_threads = bundle->n_threads;
	thread_no = bundle->thread_no;
	fbin = bundle->fbin;
	lock = bundle->lock;
	n_lines = bundle->n_lines;

	for (i = 0; i < n_bc; ++i) {
		if (i % n_threads == thread_no) {
			correct_umi(CBs + i, bundle->n_refs);
			count_genes(CBs + i, i + 1, cnt, bundle->n_refs, fbin, lock, n_lines);
		}
	}

	free(cnt);
	pthread_exit(NULL);
}

void do_correct_umi(struct opt_count_t *opt, int n_refs, const char *out_path)
{
	FILE *fbin;
	fbin = xfopen(out_path, "wb");

	xfwrite(&n_refs, sizeof(int), 1, fbin);
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
		bundles[i].n_refs = n_refs;
		pthread_create(t + i, &attr, umi_worker, bundles + i);
	}

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(t[i], NULL);
	
	fseek(fbin, 8L, SEEK_SET);
	xfwrite(&n_lines, sizeof(uint64_t), 1, fbin);
	xwfclose(fbin);

	free(t);
	free(bundles);
}

void quantification(struct opt_count_t *opt, struct kmhash_t *h,
			struct ref_info_t *ref)
{
	bc_len = opt->lib.bc_len;
	umi_len = opt->lib.umi_len;

	int *index = cut_off_barcode(h);
	__VERBOSE("Done cutting off barcode\n");
	correct_barcode(h, index);
	__VERBOSE("Done correcting barcode\n");
	free(index);

	int len = strlen(opt->out_dir);
	char path[len + 15];
	memcpy(path, opt->out_dir, len);

	concat_str(path, len, "/barcodes.tsv", 13);
	print_barcodes(path);

	concat_str(path, len, "/features.tsv", 13);
	print_refs(ref, path);

	concat_str(path, len, "/matrix.bin", 11);
	do_correct_umi(opt, ref->n_refs, path);

	write_mtx(path, len);
	free(CBs);
}