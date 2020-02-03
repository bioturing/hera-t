#include <pthread.h>

#include "attribute.h"
#include "verbose.h"
#include "opt.h"
#include "read_fastq.h"
#include "process_rna.h"
#include "process_tag.h"
#include "utils.h"

struct worker_data_t {
	pthread_t *worker_threads;
	struct thread_data_t *thread_data;
	pthread_mutex_t *lock_count;
	pthread_attr_t *attr;
	int n_threads;
};

struct thread_data_t {
	int thread_num;
	struct library_t *lib;
	struct dqueue_t *q;
	struct bc_hash_t *bc_table;
	pthread_mutex_t *lock_count;
	pthread_mutex_t *lock_hash;
	int type;
	int (*map_read)(struct read_t*, int);
	void (*print_count)(int);
};

/*****************************************
*             GENERAL PROCESS            *
*****************************************/

void print_input(struct input_t *input)
{
	__VERBOSE("[Process read]\n Index: %s\n Directory: %s\n Sample name: %s\n",
		input->ref, input->in_dir, input->name);
	__VERBOSE(" Number of pair files: %u\n", input->n_files);

	int i;
	for (i = 0; i < input->n_files; ++i)
		__VERBOSE("\t%s -- %s\n", input->left_file[i], input->right_file[i]);
}

void process_read(struct opt_count_t *opt, struct input_t *input,
		  struct bc_hash_t *bc_table,
		  void (*action)(struct worker_data_t),
		  int type)
{
	print_input(input);

	struct worker_data_t worker_data;
	pthread_t *worker_threads;
	struct thread_data_t *thread_data;
	pthread_attr_t attr;
	pthread_mutex_t lock_count;
	pthread_mutex_t lock_hash;
	struct dqueue_t *q;
	int i;
	
	q = init_read_fastq(opt->n_threads, input);
	worker_threads = calloc(opt->n_threads, sizeof(pthread_t));
	thread_data = calloc(opt->n_threads, sizeof(struct thread_data_t));
	pthread_mutex_init(&lock_count, NULL);
	pthread_mutex_init(&lock_hash, NULL);

	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	worker_data.worker_threads = worker_threads;
	worker_data.thread_data = thread_data;
	worker_data.attr = &attr;
	worker_data.n_threads = opt->n_threads;

	for (i = 0; i < opt->n_threads; ++i) {
		thread_data[i].thread_num = i;
		thread_data[i].lib = &opt->lib;
		thread_data[i].q = q;
		thread_data[i].bc_table = bc_table;
		thread_data[i].type = type;
		thread_data[i].lock_count = &lock_count;
		thread_data[i].lock_hash = &lock_hash;
	}

	(*action)(worker_data);

	free(worker_threads);
	free(thread_data);
	finish_read_fastq(q);
}

void align_read(struct read_t *read1, struct read_t *read2,
		struct thread_data_t *data, int r1_len)
{
	if (read1->len < r1_len)
		__ERROR("Read lenght of %s is shorter than |BC| + |UMI|.\n \
			Expect > %u.\n Receive %u.\n", 
			read1->name, r1_len, read1->len);

	if (check_polyA(read1->seq, r1_len))
		return;

	int ref = data->map_read(read2, data->thread_num);

	if (ref == -1)
		return;

	struct library_t *lib = data->lib;
	int64_t bc_idx = seq2num(read1->seq + lib->bc_pos, lib->bc_len);
	int64_t umi_ref = seq2num(read1->seq + lib->umi_pos, lib->umi_len);

	umi_ref = umi_ref << GENE_BIT_LEN | ref;
	pthread_mutex_lock(data->lock_hash);
	add_bc_umi(data->bc_table, bc_idx, umi_ref, data->type);
	pthread_mutex_unlock(data->lock_hash);
}

void *align_worker(void *thread_data)
{
	struct thread_data_t *data = (struct thread_data_t *)thread_data;

	struct dqueue_t *q = data->q;
	pthread_mutex_t *lock_count = data->lock_count;

	struct read_t read1, read2;
	struct pair_buffer_t *own_buf, *ext_buf;
	own_buf = init_pair_buffer();

	char *buf1, *buf2;
	int pos1, pos2, rc1, rc2;
	int input_format;
	int r1_len = data->lib->bc_len + data->lib->umi_len;

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

			align_read(&read1, &read2, data, r1_len);

			if (rc1 == READ_END)
				break;
		}
		data->print_count(data->thread_num);
	}

	free_pair_buffer(own_buf);
	pthread_exit(NULL);
}

/*****************************************
*               PROCESS RNA              *
*****************************************/

void rna_work(struct worker_data_t worker_data)
{
	int i;
	int n_threads = worker_data.n_threads;
	pthread_t *worker_threads = worker_data.worker_threads;
	pthread_attr_t *attr = worker_data.attr;
	struct thread_data_t *thread_data = worker_data.thread_data;

	init_rna_threads(n_threads);

	for (i = 0; i < n_threads; ++i) {
		thread_data[i].map_read = &align_rna;
		thread_data[i].print_count = &print_rna_count;
		
		pthread_create(worker_threads + i, attr, align_worker, thread_data + i);
	}

	for (i = 0; i < n_threads; ++i)
		pthread_join(worker_threads[i], NULL);
}

void process_rna(struct opt_count_t *opt, struct bc_hash_t *bc_table,
		struct ref_info_t *ref)
{
	load_rna_index(opt->rna->ref, opt->count_intron);
	add_rna_ref(ref);
	process_read(opt, opt->rna, bc_table, &rna_work, RNA_PRIOR);
	print_rna_stat(opt->n_threads, opt->count_intron);
	destroy_rna_index(opt->n_threads);
}

/*****************************************
*               PROCESS TAG              *
*****************************************/

void tag_work(struct worker_data_t worker_data)
{
	int i;
	int n_threads = worker_data.n_threads;
	pthread_t *worker_threads = worker_data.worker_threads;
	pthread_attr_t *attr = worker_data.attr;
	struct thread_data_t *thread_data = worker_data.thread_data;

	init_tag_threads(n_threads);

	for (i = 0; i < n_threads; ++i) {
		thread_data[i].map_read = &align_tag;
		thread_data[i].print_count = &print_tag_count;
		
		pthread_create(worker_threads + i, attr, align_worker, thread_data + i);
	}

	for (i = 0; i < n_threads; ++i)
		pthread_join(worker_threads[i], NULL);
}

void process_tag(struct opt_count_t *opt, struct input_t *input,
		struct bc_hash_t *bc_table, struct ref_info_t *ref, int type)
{
	build_tag_ref(input, ref, type);

	process_read(opt, input, bc_table, &tag_work, TAG_PRIOR);
	print_tag_stat(opt->n_threads);

	destroy_tag_ref();
}


/*****************************************
*                  GENERAL               *
*****************************************/

void single_cell(int argc, char *argv[])
{
	struct opt_count_t *opt = get_opt_count(argc, argv);
	char log_path[strlen(opt->out_dir) + 10];
	strcpy(log_path, opt->out_dir);
	strcat(log_path, "/run.log");
	init_log(log_path);

	log_write("VERSION: %d.%d\n", PROG_VERSION_MAJOR, PROG_VERSION_MINOR);
	log_write("COMMAND: ");
	int i;
	for (i = 0; i < argc; ++i)
		log_write("%s ", argv[i]);
	log_write("\n");

	struct bc_hash_t *bc_table = init_bc_hash();
	struct ref_info_t *ref = init_ref_info();

	if (opt->rna)
		process_rna(opt, bc_table, ref);
	
	if (opt->cell)
		process_tag(opt, opt->cell, bc_table, ref, 1);

	if (opt->protein)
		process_tag(opt, opt->protein, bc_table, ref, 2);

	if (opt->crispr)
		process_tag(opt, opt->crispr, bc_table, ref, 3);

	quantification(opt, bc_table, ref);

	destroy_bc_hash(bc_table);
}