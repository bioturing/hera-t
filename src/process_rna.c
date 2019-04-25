#include "process_rna.h"
#include "genome.h"
#include "alignment.h"
#include "io_utils.h"
#include "hash_table.h"

static struct genome_info_t genome;
static struct gene_info_t genes;
static struct transcript_info_t trans;
static struct bwt_t bwt;

/*****************************************
*                LOAD INDEX              *
*****************************************/

void init_bwt(const char *path, int32_t count_intron)
{
	extern struct bwt_t bwt;
	bwt_load(path, &bwt);
	genome_init_bwt(&bwt, count_intron);
}

void init_ref_index(const char *path)
{
	extern struct genome_info_t genome;
	extern struct gene_info_t genes;
	extern struct transcript_info_t trans;

	FILE *fp = xfopen(path, "rb");
	xfread(&genome.n, sizeof(int), 1, fp);
	xfread(&genome.l_name, sizeof(int), 1, fp);
	genome.chr_len = malloc(genome.n * sizeof(bioint_t));
	genome.chr_name = malloc(genome.n * genome.l_name);
	xfread(genome.chr_len, sizeof(bioint_t), genome.n, fp);
	xfread(genome.chr_name, genome.l_name, genome.n, fp);

	xfread(&genes.n, sizeof(int), 1, fp);
	xfread(&genes.l_name, sizeof(int), 1, fp);
	xfread(&genes.l_id, sizeof(int), 1, fp);
	genes.chr_idx = malloc(genes.n * sizeof(int));
	genes.gene_name = malloc(genes.n * genes.l_name);
	genes.gene_id = malloc(genes.n * genes.l_id);
	genes.strand = malloc(genes.n);
	xfread(genes.chr_idx, sizeof(int), genes.n, fp);
	xfread(genes.gene_name, genes.l_name, genes.n, fp);
	xfread(genes.gene_id, genes.l_id, genes.n, fp);
	xfread(genes.strand, 1, genes.n, fp);

	xfread(&trans.n, sizeof(int), 1, fp);
	xfread(&trans.l_id, sizeof(int), 1, fp);
	trans.tran_id = malloc(trans.n * trans.l_id);
	trans.tran_len = malloc(trans.n * sizeof(int));
	trans.tran_beg = malloc((trans.n + 1) * sizeof(int));
	trans.gene_idx = malloc(trans.n * sizeof(int));
	xfread(trans.tran_id, trans.l_id, trans.n, fp);
	xfread(trans.tran_len, sizeof(int), trans.n, fp);
	xfread(trans.tran_beg, sizeof(int), trans.n + 1, fp);
	xfread(trans.gene_idx, sizeof(int), trans.n, fp);
	trans.idx = malloc(trans.tran_beg[trans.n] * sizeof(int));
	trans.seq = malloc(trans.tran_beg[trans.n]);
	xfread(trans.idx, sizeof(int), trans.tran_beg[trans.n], fp);
	xfread(trans.seq, 1, trans.tran_beg[trans.n], fp);

	trans.n_exon = malloc(trans.n * sizeof(int));
	xfread(trans.n_exon, sizeof(int), trans.n, fp);
	trans.exons = malloc(trans.n * sizeof(struct exon_t *));
	int i;
	for (i = 0; i < trans.n; ++i) {
		trans.exons[i] = malloc(trans.n_exon[i] * sizeof(struct exon_t));
		xfread(trans.exons[i], sizeof(struct exon_t), trans.n_exon[i], fp);
	}
	fclose(fp);

	alignment_init_ref_info(&genes, &trans);
}

void load_rna_index(const char *prefix, int32_t count_intron)
{
	int len = strlen(prefix);
	char path[len + 10];
	memcpy(path, prefix, len);

	concat_str(path, len, ".bwt", 4);
	__VERBOSE("Loading BWT...\n");
	init_bwt(path, count_intron);

	concat_str(path, len, ".info", 5);
	__VERBOSE("Loading transcripts and genes info...\n");
	init_ref_index(path);

	concat_str(path, len, ".hash", 5);
	__VERBOSE("Loading kmer hash table...\n");
	alignment_init_hash(path);
}

void destroy_rna_index()
{
	extern struct bwt_t bwt;
	extern struct genome_info_t genome;
	extern struct gene_info_t genes;
	extern struct transcript_info_t trans;
	free_cons_hash();
	bwt_destroy(&bwt);
	free(trans.seq);
}

void add_rna_ref(struct ref_info_t *ref)
{
	ref->n_refs = genes.n;
	ref->type[0] = ref->n_refs;

	ref->ref_text = realloc(ref->ref_text, genes.n * genes.l_name);
	memcpy(ref->ref_text, genes.gene_name, genes.n * genes.l_name);

	ref->gene_id = realloc(ref->gene_id, genes.n * genes.l_id);
	memcpy(ref->gene_id, genes.gene_id, genes.n * genes.l_id);

	ref->ref_iter = realloc(ref->ref_iter, ref->n_refs * sizeof(int));
	ref->gene_iter = realloc(ref->gene_iter, ref->n_refs * sizeof(int));
	int i;
	for (i = 0; i < ref->n_refs; ++i) {
		ref->ref_iter[i] = i * genes.l_name;
		ref->gene_iter[i] = i * genes.l_id;
	}
}

/*****************************************
*               RNA ALIGNMENT            *
*****************************************/

struct bundle_data_t *bundles;
struct align_stat_t *count;

void init_rna_threads(int n_threads)
{
	int i;
	bundles = malloc(n_threads * sizeof(struct bundle_data_t));
	for (i = 0; i < n_threads; ++i)
		init_bundle(bundles + i);
	count = calloc(n_threads + 1, sizeof(struct align_stat_t));
}

void destroy_rna_threads(int n_threads)
{
	int i;
	for (i = 0; i < n_threads; ++i)
		destroy_bundle(bundles + i);
	free(bundles);
	free(count);
}

int align_rna(struct read_t *read, int thread_num)
{
	struct align_stat_t *c = count + (thread_num + 1);
	struct bundle_data_t *bun = bundles + thread_num;
	reinit_bundle(bun);

	return align_chromium_read(read, bun, c);
}

void update_result(struct align_stat_t *res, struct align_stat_t *add)
{
	res->nread += add->nread;
	res->exon += add->exon;
	res->unmap += add->unmap;
	res->intron += add->intron;
	res->intergenic += add->intergenic;
}

void print_rna_count(int thread_num)
{
	struct align_stat_t *c = count + (thread_num + 1);
	update_result(count, c);
	memset(c, 0, sizeof(struct align_stat_t));
	__VERBOSE("\rNumber of processed reads: %ld", count[0].nread);
}

void print_rna_stat(int n_threads, int count_intron)
{
	int i;
	for (i = 1; i <= n_threads; ++i)
		update_result(count, count + i);
	__VERBOSE("\n");
	__VERBOSE_LOG("INFO", "Total number of reads               : %10ld\n",
			count[0].nread);
	__VERBOSE_LOG("INFO", "Number of exonic mapped reads       : %10ld\n",
			count[0].exon);
	if (count_intron){
		__VERBOSE_LOG("INFO", "Number of intronic reads            : %10ld\n", count[0].intron);
		__VERBOSE_LOG("INFO", "Number of intergenic reads          : %10ld\n", count[0].intergenic);
	} else {
		__VERBOSE_LOG("INFO", "Number of nonexonic reads           : %10ld\n", count[0].intergenic);
	}
	__VERBOSE_LOG("INFO", "Number of unmapped reads            : %10ld\n", count[0].unmap);
}
