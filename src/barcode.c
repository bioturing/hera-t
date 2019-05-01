#include <pthread.h>

#include "barcode.h"
#include "io_utils.h"
#include "verbose.h"
#include "utils.h"

#define __get_gene(x) ((int)((x) & GENE_MASK))
#define __get_umi(x) ((x) >> GENE_BIT_LEN)

#define CUT_OFF 0.01

const uint64_t pow5_r[] = {1ull, 5ull, 25ull, 125ull, 625ull, 3125ull, 15625ull, 78125ull, 390625ull, 1953125ull, 9765625ull, 48828125ull, 244140625ull, 1220703125ull, 6103515625ull, 30517578125ull, 152587890625ull, 762939453125ull, 3814697265625ull, 19073486328125ull, 95367431640625ull, 476837158203125ull, 2384185791015625ull, 11920928955078125ull, 59604644775390625ull, 298023223876953125ull, 1490116119384765625ull};

struct umi_bundle_t {
	int n_threads;
	int thread_no;
	pthread_mutex_t *lock;
	FILE *fbin;
	uint64_t *n_lines;
	int n_refs;
};

int32_t n_bc;
int16_t bc_len, umi_len;
struct bc_hash_t *bc_hash;

struct ref_info_t *init_ref_info()
{
	struct ref_info_t *ref = malloc(sizeof(struct ref_info_t));
	ref->n_refs = 0;
	ref->text_iter = malloc(1);
	ref->ref_text = malloc(1);
	ref->ref_id = malloc(1);
	ref->id_iter = malloc(1);
	ref->text_len = ref->id_len = 0;
	memset(ref->type, 0, 4 * sizeof(int));

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
	struct umi_hash_t *c_a = (struct umi_hash_t *) a;
	struct umi_hash_t *c_b = (struct umi_hash_t *) b;
	return c_b->type == c_a->type? c_b->count - c_a->count: c_b->type - c_a->type;
}

static int compare_umi(const void *a, const void *b)
{
	struct umi_hash_t *c_a = (struct umi_hash_t *) a;
	struct umi_hash_t *c_b = (struct umi_hash_t *) b;
	if (c_a->type < RNA_PRIOR)
		return c_b->count;
	else if (c_b->type < RNA_PRIOR)
		return -c_a->count;
	else
		return c_b->count - c_a->count;
}

void cut_off_barcode()
{
	int i;
	struct umi_hash_t *umi = bc_hash->umi;

	qsort(umi, n_bc, sizeof(struct umi_hash_t), compare_umi);
	float cut_off = umi[0].count * CUT_OFF;
	for (i = 0; i < n_bc; ++i)
		if ((float) umi[i].count < cut_off || umi[i].type < RNA_PRIOR)
			n_bc = i - 1;
}

void merge_bc(int bc1, int bc2)
{
	khiter_t k;
	struct umi_hash_t *umi = bc_hash->umi;

	if (umi[bc1].type == umi[bc2].type && umi[bc2].count / 2 < umi[bc1].count)
		return;

	khash_t(bc_umi) *h = umi[bc1].h;

	for (k = kh_begin(h); k != kh_end(h); ++k)
		if (kh_exist(h, k))
			add_umi(umi + bc2, kh_key(h, k), kh_value(h, k));
	umi[bc1].type = -bc2;
}

void correct_barcode()
{
	uint64_t bc_idx, new_bc, tmp_idx;
	int i, k, ch, c, max_count, merge_iter, flag, count;
	khiter_t iter;
	struct umi_hash_t *umi = bc_hash->umi;
	khash_t(bc_umi) *h = bc_hash->h;
	int *index = malloc(bc_hash->n_bc * sizeof(int));

	n_bc = bc_hash->n_bc;
	qsort(umi, n_bc, sizeof(struct umi_hash_t), compare_cell);
	for (i = 0; i < n_bc; ++i) {
		k = kh_get(bc_umi, h, umi[i].idx);
		index[kh_value(h, k)] = i;
	}

	for (i = bc_hash->n_bc; i > 0 ; --i) {
		if (umi[i].type < 0)
			continue;

		bc_idx = umi[i].idx;
		max_count = 0;
		for (k = 0; k < bc_len; ++k) {
			flag = 0;
			ch = (bc_idx / pow5_r[k]) % 5;
			tmp_idx = bc_idx - ch * pow5_r[k];
			for (c = 0; c < NNU; ++c) {
				if (ch == c)
					continue;
				new_bc = tmp_idx + pow5_r[k] * c;
				iter = kh_get(bc_umi, h, new_bc);
				if (iter == kh_end(h))
					continue;

				iter = index[kh_value(h, iter)];
				while (umi[iter].type < 0) {
					if (i == -umi[iter].type){
						flag = 1;
						break;
					}
					iter = -umi[iter].type;
				}
				if (flag)
					continue;

				// ensure barcode of rna is prior
				count = umi[iter].count +
					(umi[iter].type & RNA_PRIOR) * umi[0].count;
				if (count > max_count) {
					max_count = count;
					merge_iter = kh_get(bc_umi, h, umi[iter].idx);
				}
			}
		}
		if (max_count)
			merge_bc(i, index[kh_value(h, merge_iter)]);
	}
	free(index);
}

void print_barcodes(const char *out_path)
{
	char *seq;
	FILE *fp;
	int i;

	fp = xfopen(out_path, "w");
	struct umi_hash_t *umi = bc_hash->umi;
	for (i = 0; i < n_bc; ++i) {
		seq = num2seq(umi[i].idx, bc_len);
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

		start = ref->type[1];
	}

	// Print protein measurement
	if (ref->type[2]) {
		for (i = start; i < ref->type[2]; ++i)
			fprintf(fp, "%s\t%s\tProtein measurement\n",
				ref->ref_id + ref->id_iter[i],
				ref->ref_text + ref->text_iter[i]);

		start = ref->type[2];
	}

	// Print CRISPR capture
	if (ref->type[3]) {
		for (i = start; i < ref->type[3]; ++i)
			fprintf(fp, "%s\t%s\tCRISPR capture\n",
				ref->ref_id + ref->id_iter[i],
				ref->ref_text + ref->text_iter[i]);
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

void merge_count(khash_t(bc_umi) *h, khiter_t src, khiter_t dst)
{
	if (kh_value(h, src) > kh_value(h, dst)) {
		kh_value(h, src) += kh_value(h, dst);
		kh_value(h, dst) = -1;
	} else {
		kh_value(h, dst) += kh_value(h, src);
		kh_value(h, src) = -1;
	}
}

void correct_umi(struct umi_hash_t *umi)
{
	uint64_t umi_idx, new_umi, key;
	int i, gene, ch, c, flag;
	khash_t(bc_umi) *h = umi->h;
	khiter_t k, iter;

	for (k = kh_begin(h); k != kh_end(h); ++k) {
		if (!kh_exist(h, k) || kh_value(h, k) < 0)
			continue;

		key = kh_key(h, k);
		umi_idx = __get_umi(key);
		gene = __get_gene(key);
		flag = 0;

		for (i = 0; i < umi_len; ++i) {
			ch = (umi_idx / pow5_r[i]) % 5;
			for (c = 0; c < NNU; ++c) {
				if (c == ch)
					continue;
				new_umi = umi_idx + (c - ch) * pow5_r[i];
				new_umi = (new_umi << GENE_BIT_LEN | gene);
				iter = kh_get(bc_umi, h, new_umi);
				if (iter != kh_end(h) && kh_value(h, iter) > 0) {
					merge_count(h, k, iter);
					flag = 1;
					break;
				}
			}
			if (flag)
				break;
		}
	}
}

void count_genes(struct umi_hash_t *umi, int bc_pos, int *cnt, int n_refs,
		FILE *fbin, pthread_mutex_t *lock, uint64_t *n_lines)
{
	khiter_t k;
	int i;
	khash_t(bc_umi) *h = umi->h;

	memset(cnt, 0, n_refs * sizeof(int));
	for (k = kh_begin(h); k != kh_end(h); ++k)
		if (kh_exist(h, k) && kh_value(h, k) > 0)
			++cnt[__get_gene(kh_key(h, k))];

	pthread_mutex_lock(lock);
	for (i = 1; i <= n_refs; ++i) {
		if (cnt[i - 1]) {
			++*n_lines;
			xfwrite(&i, sizeof(int), 1, fbin);
			xfwrite(&bc_pos, sizeof(int), 1, fbin);
			xfwrite(cnt + (i - 1), sizeof(int), 1, fbin);
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
			correct_umi(bc_hash->umi + i);
			count_genes(bc_hash->umi + i, i + 1, cnt, bundle->n_refs, fbin, lock, n_lines);
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

void print_molecule_info(char *out_path, struct ref_info_t *ref)
{
	char *bc_str, *umi_str;
	FILE *fp;
	int i, r;
	int64_t value, umi_idx;
	khiter_t k;
	khash_t(bc_umi) *h;
	struct umi_hash_t *umi = bc_hash->umi;

	fp = xfopen(out_path, "w");
	fprintf(fp, "barcode\tUMI\treference\tread_count\n");
	for (i = 0; i < n_bc; ++i) {
		bc_str = num2seq(umi[i].idx, bc_len);
		h = umi[i].h;

		for (k = kh_begin(h); k != kh_end(h); ++k) {
			if (!kh_exist(h, k) || kh_value(h, k) < 0)
				continue;
			value = kh_key(h, k);
			umi_idx = __get_umi(value);
			umi_str = num2seq(umi_idx, umi_len);
			r = __get_gene(value);
			fprintf(fp, "%s\t%s\t%s\t%u\n", bc_str, umi_str,
				ref->ref_text + ref->text_iter[r],
				kh_value(h, k));
			free(umi_str);
		}
		free(bc_str);
	}

	fclose(fp);
}

void quantification(struct opt_count_t *opt, struct bc_hash_t *h,
			struct ref_info_t *ref)
{
	bc_len = opt->lib.bc_len;
	umi_len = opt->lib.umi_len;
	bc_hash = h;

	correct_barcode();
	__VERBOSE("Done correcting barcode\n");
	cut_off_barcode();
	__VERBOSE("Done cutting off barcode\n");

	int len = strlen(opt->out_dir);
	char path[len + 20];
	memcpy(path, opt->out_dir, len);

	concat_str(path, len, "/barcodes.tsv", 13);
	print_barcodes(path);

	concat_str(path, len, "/features.tsv", 13);
	print_refs(ref, path);

	concat_str(path, len, "/matrix.bin", 11);
	do_correct_umi(opt, ref->n_refs, path);

	write_mtx(path, len);

	concat_str(path, len, "/molecule_info.tsv", 18);
	print_molecule_info(path, ref);
}