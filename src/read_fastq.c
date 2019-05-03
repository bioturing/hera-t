#include "read_fastq.h"
#include "get_buffer.h"
#include "attribute.h"
#include "utils.h"
#include "pthread_barrier.h"

#define BUF_SIZE 		1048576

struct producer_bundle_t {
	int n_files;
	int thread_no;
	char *in_dir;
	char **left_file;
	char **right_file;
	struct dqueue_t *q;
};

static struct gb_pair_data *input_streams;
static pthread_mutex_t lock;
static pthread_barrier_t barrier;
static int n_consumer, n_producer;
static pthread_t *producer_threads;
struct producer_bundle_t *producer_bundles;

void free_pair_buffer(struct pair_buffer_t *p)
{
	if (!p) return;
	free(p->buf1);
	free(p->buf2);
	free(p);
}

struct pair_buffer_t *init_pair_buffer()
{
	struct pair_buffer_t *ret = malloc(sizeof(struct pair_buffer_t));
	ret->buf1 = malloc(BUF_SIZE + 1);
	ret->buf2 = malloc(BUF_SIZE + 1);
	return ret;
}

struct dqueue_t *init_dqueue_PE(int cap)
{
	struct dqueue_t *ret = init_dqueue(cap);
	struct pair_buffer_t *p;
	int i;
	for (i = 0; i < cap; ++i) {
		p = init_pair_buffer();
		d_enqueue_out(ret, p);
	}
	return ret;
}


void *pair_producer_worker(void *data)
{
	struct producer_bundle_t *bun = (struct producer_bundle_t *)data;
	struct dqueue_t *q = bun->q;
	struct pair_buffer_t *own_buf = init_pair_buffer();
	struct pair_buffer_t *external_buf;
	struct gb_pair_data *stream;
	int i, thread_no, n_files;
	int64_t offset;
	char **left_file, **right_file;
	int len = strlen(bun->in_dir);
	char left_path[MAX_PATH], right_path[MAX_PATH];

	memcpy(left_path, bun->in_dir, len);
	memcpy(right_path, bun->in_dir, len);

	left_file = bun->left_file;
	right_file = bun->right_file;
	thread_no = bun->thread_no;
	n_files = bun->n_files;
	for (i = thread_no; i < n_files; i += n_producer) {
		stream = input_streams + i;
		concat_str(left_path, len, left_file[i], strlen(left_file[i]));
		concat_str(right_path, len, right_file[i], strlen(right_file[i]));
		gb_pair_init(stream, left_path, right_path);
		while ((offset = gb_get_pair(stream, &own_buf->buf1, &own_buf->buf2)) != -1) {
			own_buf->input_format = stream->type;
			external_buf = d_dequeue_out(q);
			d_enqueue_in(q, own_buf);
			own_buf = external_buf;
		}
	}
	free_pair_buffer(own_buf);

	pthread_barrier_wait(&barrier);
	while (1) {
		pthread_mutex_lock(&lock);
		int cur = n_consumer;
		if (n_consumer > 0)
			--n_consumer;
		pthread_mutex_unlock(&lock);
		if (cur == 0)
			break;
		external_buf = d_dequeue_out(q);
		free_pair_buffer(external_buf);
		d_enqueue_in(q, NULL);
	}

	pthread_exit(NULL);
}

struct dqueue_t *init_read_fastq(int n_producer_threads, struct input_t *input)
{
	pthread_attr_t attr;
	struct dqueue_t *q;
	int i;

	pthread_attr_init(&attr);
	pthread_attr_setstacksize(&attr, THREAD_STACK_SIZE);
	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	// parallel read files
	// must always >= n_thread * 2 in order to avoid deadlock
	q = init_dqueue_PE(n_producer_threads * 2); 
	n_consumer = n_producer_threads * 2;
	n_producer = __min(input->n_files, n_producer_threads);
	producer_bundles = malloc(n_producer * sizeof(struct producer_bundle_t));
	producer_threads = calloc(n_producer, sizeof(pthread_t));

	pthread_mutex_init(&lock, NULL);
	pthread_barrier_init(&barrier, NULL, n_producer);

	input_streams = calloc(input->n_files, sizeof(struct gb_pair_data));
	for (i = 0; i < n_producer; ++i) {
		producer_bundles[i].thread_no = i;
		producer_bundles[i].in_dir = input->in_dir;
		producer_bundles[i].n_files = input->n_files;
		producer_bundles[i].left_file = input->left_file;
		producer_bundles[i].right_file = input->right_file;
		producer_bundles[i].q = q;
		pthread_create(producer_threads + i, &attr,
				pair_producer_worker, producer_bundles + i);
	}

	return q;
}

void finish_read_fastq(struct dqueue_t *q)
{
	free(producer_threads);
	free(producer_bundles);
	free(input_streams);
	dqueue_destroy(q);
}
