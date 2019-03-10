#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>

#include "interval_tree.h"
#include "alignment.h"
#include "attribute.h"
#include "barcode.h"
#include "antibody.h"
#include "bwt.h"
#include "genome.h"
#include "get_buffer.h"
#include "hash_table.h"
#include "io_utils.h"
#include "opt.h"
#include "pthread_barrier.h"
#include "verbose.h"

static struct genome_info_t genome;
static struct gene_info_t genes;
static struct transcript_info_t trans;

static struct bwt_t bwt;

void free_align_data();

void load_index(const char *prefix, int32_t count_intron);

void single_cell_process(struct opt_count_t *opt);

void *align_worker(void *data);

void *producer_worker(void *data);

void update_result(struct align_stat_t *res, struct align_stat_t *add);

void debug_index()
{
	extern struct transcript_info_t trans;
	int i;
	FILE *fp = fopen("debug.fa", "wb");
	char *name, *seq;
	seq = NULL;
	for (i = 0; i < trans.n; ++i) {
		name = trans.tran_id + i * trans.l_id;
		fprintf(fp, ">%s\n", name);
		seq = realloc(seq, trans.tran_len[i] + 1);
		memcpy(seq, trans.seq + trans.tran_beg[i], trans.tran_len[i]);
		seq[trans.tran_len[i]] = '\0';
		fprintf(fp, "%s\n", seq);
	}
	fclose(fp);
	exit(0);
}

void single_cell(int pos, int argc, char *argv[])
{
	struct opt_count_t *opt = get_opt_count(argc - pos, argv + pos);
	char prefix[1024];
	char tmp_dir[1024];
	strcpy(prefix, opt->out_dir); strcat(prefix, "/");
	strcat(prefix, opt->prefix);

	strcpy(tmp_dir, prefix); strcat(tmp_dir, ".log");
	init_log(tmp_dir);

	log_write("VERSION: %d.%d\n", PROG_VERSION_MAJOR, PROG_VERSION_MINOR);
	log_write("COMMAND: ");
	int i;
	for (i = 0; i < argc; ++i)
		log_write("%s ", argv[i]);
	log_write("\n");

	load_index(opt->index, opt->count_intron);

	extern struct gene_info_t genes;
	init_barcode(&genes, opt->lib);

	single_cell_process(opt);
}

void check_some_statistics(struct kmhash_t *h)
{
	__VERBOSE("\n");
	__VERBOSE("Number of barcodes                    : %d\n", h->n_items);
	uint64_t s = 0;
	kmint_t i;
	for (i = 0; i < h->size; ++i) {
		if (h->bucks[i].idx != TOMB_STONE)
			s += h->bucks[i].umis->n_items;
	}
	__VERBOSE("Total number of UMI                   : %lu\n", s);
	__VERBOSE("Mean UMI per barcode                  : %.6f\n", s * 1.0 / h->n_items);
}

void rna_work(struct worker_data_t worker_data)
{
	int i;
	char path[1024];
	struct shared_fstream_t *align_fstream;

	struct opt_count_t *opt = worker_data.opt;
	struct align_stat_t result = worker_data.result;
	struct worker_bundle_t *worker_bundles = worker_data.worker_bundles;
	pthread_t *producer_threads = worker_data.producer_threads;
	pthread_t *worker_threads = worker_data.worker_threads;
	pthread_attr_t *attr = worker_data.attr;

	strcpy(path, opt->out_dir); strcat(path, "/");
	strcat(path, opt->prefix); strcat(path, ".align.tsv");
	if (opt->is_dump_align)
		align_fstream = init_shared_stream(path, opt->n_threads);
	else
		align_fstream = NULL;

	for (i = 0; i < opt->n_threads; ++i) {
		worker_bundles[i].q = worker_data.q;
		worker_bundles[i].bc_table = worker_data.bc_table;
		worker_bundles[i].lock_count = worker_data.lock_count;
		worker_bundles[i].result = &result;
		worker_bundles[i].lib = opt->lib;
		if (opt->is_dump_align)
			worker_bundles[i].align_fstream = align_fstream + i;
		else
			worker_bundles[i].align_fstream = NULL;

		worker_bundles[i].init_bundle = &init_bundle;
		worker_bundles[i].destroy_bundle = &destroy_bundle;
		worker_bundles[i].map_read = &align_chromium_read;
		
		pthread_create(worker_threads + i, attr, align_worker, worker_bundles + i);
	}

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	__VERBOSE("\rNumber of processed reads: %d\n", result.nread);
	__VERBOSE_LOG("INFO", "Total number of reads             : %10u\n", result.nread);
	__VERBOSE_LOG("INFO", "Number of exonic mapped reads       : %10u\n", result.exon);
	if (opt->count_intron){
		__VERBOSE_LOG("INFO", "Number of intronic reads   : %10u\n", result.intron);
		__VERBOSE_LOG("INFO", "Number of intergenic reads   : %10u\n", result.intergenic);
	} else {
		__VERBOSE_LOG("INFO", "Number of nonexonic reads   : %10u\n", result.intergenic);
	}
	__VERBOSE_LOG("INFO", "Number of unmapped reads   : %10u\n", result.unmap);

	destroy_shared_stream(align_fstream, opt->n_threads);
	free_align_data();
}


void antibody_work(struct worker_data_t worker_data)
{
	int i;
	struct opt_count_t *opt = worker_data.opt;
	struct align_stat_t result = worker_data.result;
	struct worker_bundle_t *worker_bundles = worker_data.worker_bundles;
	pthread_t *producer_threads = worker_data.producer_threads;
	pthread_t *worker_threads = worker_data.worker_threads;
	pthread_attr_t *attr = worker_data.attr;

	for (i = 0; i < worker_data.opt->n_threads; ++i) {
		worker_bundles[i].q = worker_data.q;
		worker_bundles[i].bc_table = worker_data.bc_table;
		worker_bundles[i].lock_count = worker_data.lock_count;
		worker_bundles[i].result = &result;
		worker_bundles[i].lib = opt->lib;
		worker_bundles[i].antibody_lib = worker_data.lib;

		worker_bundles[i].init_bundle = NULL;
		worker_bundles[i].destroy_bundle = NULL;
		worker_bundles[i].map_read = &map_antibody_read;
		pthread_create(worker_threads + i, attr, align_worker, worker_bundles + i);
	}

	for (i = 0; i < opt->n_files; ++i)
		pthread_join(producer_threads[i], NULL);

	for (i = 0; i < opt->n_threads; ++i)
		pthread_join(worker_threads[i], NULL);

	__VERBOSE("\rNumber of processed reads: %d\n", result.nread);
	__VERBOSE_LOG("INFO", "Total number of reads             : %10u\n", result.nread);
}

void process_read(struct kmhash_t *bc_table, struct opt_count_t *opt, 
					struct quant_data_t input)
{
	pthread_attr_t attr;
	struct dqueue_t *q;
	struct align_stat_t result;
	struct worker_data_t worker_data;
	struct worker_bundle_t *worker_bundles;
	struct producer_bundle_t *producer_bundles;
	pthread_t *producer_threads, *worker_threads;
	pthread_mutex_t producer_lock, lock_count;
	pthread_barrier_t producer_barrier;
	int n_consumer, i;

	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	// parallel read files
	q = init_dqueue_PE(opt->n_threads * 2); // must always >= n_thread * 2 in order to avoid deadlock
	n_consumer = opt->n_threads * 2;

	producer_bundles = malloc(input.n_files * sizeof(struct producer_bundle_t));
	producer_threads = calloc(input.n_files, sizeof(pthread_t));

	pthread_mutex_init(&producer_lock, NULL);
	pthread_barrier_init(&producer_barrier, NULL, input.n_files);

	for (i = 0; i < input.n_files; ++i) {
		struct gb_pair_data *data = calloc(1, sizeof(struct gb_pair_data));
		gb_pair_init(data, input.left_file[i], input.right_file[i]);

		producer_bundles[i].n_consumer = &n_consumer;
		producer_bundles[i].stream = (void *)data;
		producer_bundles[i].q = q;
		producer_bundles[i].barrier = &producer_barrier;
		producer_bundles[i].lock = &producer_lock;
		pthread_create(producer_threads + i, &attr, producer_worker, producer_bundles + i);
	}

	worker_bundles = malloc(opt->n_threads * sizeof(struct worker_bundle_t));
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));
	memset(&result, 0, sizeof(struct align_stat_t));
	pthread_mutex_init(&lock_count, NULL);

	worker_data.bc_table = bc_table;
	worker_data.opt = opt;
	worker_data.q = q;
	worker_data.worker_bundles = worker_bundles;
	worker_data.worker_threads = worker_threads;
	worker_data.producer_threads = producer_threads;
	worker_data.lock_count = &lock_count;
	worker_data.attr = &attr;
	worker_data.result = result;
	worker_data.lib = input.lib;

	(*input.action)(worker_data);

	free(producer_bundles);
	free(producer_threads);
	free(worker_bundles);
	free(worker_threads);
	dqueue_destroy(q);
}

void single_cell_process(struct opt_count_t *opt)
{
	if (opt->cell_hashing != NULL)
		build_reference(opt->cell_hashing);

	if (opt->protein_quant != NULL)
		build_reference(opt->protein_quant);

	struct kmhash_t *bc_table = init_kmhash(KMHASH_KMHASH_SIZE - 1, opt->n_threads);
	struct quant_data_t data;

	// mRna quantification
	data.n_files = opt->n_files;
	data.left_file = opt->left_file;
	data.right_file = opt->right_file;
	data.lib = NULL;
	data.action = &rna_work;

	process_read(bc_table, opt, data);

	check_some_statistics(bc_table);
	quantification(opt, bc_table);

	kmhash_destroy(bc_table);
	
	// Cell hashing
	if (opt->cell_hashing == NULL && opt->protein_quant == NULL)
		return;

	bc_table = build_hash_from_cutoff(opt->n_threads);

	if (opt->cell_hashing != NULL) {
		__VERBOSE_LOG("INFO", "Map reads for cell hashing\n");
		data.lib = opt->cell_hashing;
		data.n_files = data.lib->n_files;
		data.left_file = data.lib->left_file;
		data.right_file = data.lib->right_file;
		data.action = &antibody_work;

		process_read(bc_table, opt, data);
		antibody_quant(opt, bc_table, opt->cell_hashing);
	}

	if (opt->protein_quant != NULL) {
		__VERBOSE_LOG("INFO", "Map reads for protein measurement\n");
		data.lib = opt->protein_quant;
		data.n_files = data.lib->n_files;
		data.left_file = data.lib->left_file;
		data.right_file = data.lib->right_file;
		data.action = &antibody_work;

		process_read(bc_table, opt, data);
		antibody_quant(opt, bc_table, opt->protein_quant);
	}
}

void *align_worker(void *data)
{
	struct worker_bundle_t *bundle = (struct worker_bundle_t *)data;

	if (bundle->init_bundle != NULL)
		bundle->init_bundle(bundle);

	struct dqueue_t *q = bundle->q;
	struct align_stat_t *global_result = bundle->result;
	pthread_mutex_t *lock_count = bundle->lock_count;
	struct align_stat_t own_result;
	memset(&own_result, 0, sizeof(struct align_stat_t));
	bundle->result = &own_result;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2;
	int input_format;

	while (1) {
		ext_buf = d_dequeue_in(q);
		if (!ext_buf)
			break;
		d_enqueue_out(q, own_buf);
		own_buf = ext_buf;
		pos1 = pos2 = 0;
		buf1 = ext_buf->buf1;
		buf2 = ext_buf->buf2;
		input_format = ext_buf->input_format;
		while (1) {
			rc1 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read1, buf1, &pos1) :
				get_read_from_fa(&read1, buf1, &pos1);

			rc2 = input_format == TYPE_FASTQ ?
				get_read_from_fq(&read2, buf2, &pos2) :
				get_read_from_fa(&read2, buf2, &pos2);


			if (rc1 == READ_FAIL || rc2 == READ_FAIL)
				__ERROR("\nWrong format file\n");

			bundle->map_read(&read1, &read2, bundle);

			if (rc1 == READ_END)
				break;
		}

		pthread_mutex_lock(lock_count);
		update_result(global_result, &own_result);
		pthread_mutex_unlock(lock_count);
		memset(&own_result, 0, sizeof(struct align_stat_t));
	}

	if (bundle->destroy_bundle != NULL)
		bundle->destroy_bundle(bundle);
	free_pair_buffer(own_buf);

	pthread_exit(NULL);
}

void *producer_worker(void *data)
{
	struct producer_bundle_t *bundle = (struct producer_bundle_t *)data;
	struct dqueue_t *q = bundle->q;
	struct gb_pair_data *input_stream = bundle->stream;
	struct pair_buffer_t *own_buf = init_pair_buffer();
	struct pair_buffer_t *external_buf;
	int64_t offset;

	while ((offset = gb_get_pair(input_stream, &own_buf->buf1, &own_buf->buf2)) != -1) {
		own_buf->input_format = input_stream->type;
		external_buf = d_dequeue_out(q);
		d_enqueue_in(q, own_buf);
		own_buf = external_buf;
	}
	free_pair_buffer(own_buf);

	int cur;
	pthread_barrier_wait(bundle->barrier);
	while (1) {
		pthread_mutex_lock(bundle->lock);
		cur = *(bundle->n_consumer);
		if (*(bundle->n_consumer) > 0)
			--*(bundle->n_consumer);
		pthread_mutex_unlock(bundle->lock);
		if (cur == 0)
			break;
		external_buf = d_dequeue_out(q);
		free_pair_buffer(external_buf);
		d_enqueue_in(q, NULL);
	}

	pthread_exit(NULL);
}

void update_result(struct align_stat_t *res, struct align_stat_t *add)
{
	res->nread	+= add->nread;
	res->exon	+= add->exon;
	res->unmap	+= add->unmap;
	res->intron     += add->intron;
	res->intergenic += add->intergenic;
	__VERBOSE("\rNumber of processed reads: %d", res->nread);
}

void init_bwt(const char *path, int32_t count_intron)
{
	extern struct bwt_t bwt;
	bwt_load(path, &bwt);
	genome_init_bwt(&bwt, count_intron);
}

void init_ref_info(const char *path)
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

void load_index(const char *prefix, int32_t count_intron)
{
	char tmp_dir[1024];
	strcpy(tmp_dir, prefix); strcat(tmp_dir, ".bwt"); 
	__VERBOSE("Loading BWT...\n");
	init_bwt(tmp_dir, count_intron);

	strcpy(tmp_dir, prefix); strcat(tmp_dir, ".info");
	__VERBOSE("Loading transcripts and genes info...\n");
	init_ref_info(tmp_dir);

	strcpy(tmp_dir, prefix); strcat(tmp_dir, ".hash");
	__VERBOSE("Loading kmer hash table...\n");
	alignment_init_hash(tmp_dir);

	/*

	strcpy(tmp_dir, prefix); strcat(tmp_dir, ".tree");
	__VERBOSE("Loading gene interval...\n");
	load_interval_tree(tmp_dir, genes.strand);

	*/
}

void free_align_data()
{
	extern struct bwt_t bwt;
	extern struct genome_info_t genome;
	extern struct gene_info_t genes;
	extern struct transcript_info_t trans;
	free_cons_hash();
	bwt_destroy(&bwt);
	free(trans.seq);
}
