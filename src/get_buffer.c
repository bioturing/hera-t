#include <string.h>
#include "get_buffer.h"
#include "io_utils.h"
#include "utils.h"
#include "verbose.h"

static int get_format(char *file_path)
{
	char buf[1024];
	gzFile file;
	int n, ret;

	file = gzopen(file_path, "r");
	if (!file)
		__ERROR("Could not open file: %s", file_path);

	n = gzread(file, buf, 1);
	if (!n)
		__ERROR("Could not read file: %s", file_path);

	if (buf[0] == '@')
		ret = TYPE_FASTQ;
	else if (buf[0] == '>')
		ret = TYPE_FASTA;
	else
		__ERROR("Unsupport format of file: %s", file_path);

	gzclose(file);
	return ret;
}

void gb_pair_init(struct gb_pair_data *data, char *file_path1, char *file_path2)
{
	if (!strcmp(file_path1, file_path2))
		__ERROR("Two identical read files");

	data->type = get_format(file_path1);
	if (get_format(file_path2) != data->type)
		__ERROR("Format in two read files are not equal");
	data->offset = 0;

	data->file1.name = file_path1;
	data->file2.name = file_path2;
	data->file1.fi = gzopen(file_path1, "r");
	data->file2.fi = gzopen(file_path2, "r");
	data->file1.buf = malloc(SIZE_1MB + 1);
	data->file2.buf = malloc(SIZE_1MB + 1);
	data->file1.size = data->file2.size = 0;
	data->file1.is_eof = data->file2.is_eof = 0;
	data->finish_flag = 0;
	data->warning_flag = 0;
	data->file1.processed = 0;
	data->file2.processed = 0;
}

void gb_pair_destroy(struct gb_pair_data *data)
{
	free(data->file1.buf);
	free(data->file2.buf);
	gzclose(data->file1.fi);
	gzclose(data->file2.fi);
}

void gb_single_init(struct gb_single_data *data, char *file_path)
{
	data->type = get_format(file_path);
	data->file.name = file_path;
	data->offset = 0;
	data->file.fi = gzopen(file_path, "r");

	data->file.buf = malloc(SIZE_1MB + 1);
	data->file.size = 0;
	data->file.is_eof = 0;
	data->finish_flag = 0;
	data->file.processed = 0;
}

void gb_single_destroy(struct gb_single_data *data)
{
	free(data->file.buf);
	gzclose(data->file.fi);
}

/*
 * return 0 if read is not found
 * otherwise return the next read's position
 */
static int get_nxt_pos(struct gb_file_inf *data, int prev, int type)
{
	int i = prev;
	char *buf = data->buf;
	int size = data->size;
	int id = 0;
	int mod = (type == TYPE_FASTQ ? 3 : 1);

	while (1) {
		for (; buf[i] != '\n' && i < size; ++i);
		if (i == size) {
			if (!data->is_eof) 
				break; 
		} else { 
			++i; 
		}
		id = (id + 1) & mod;
		if (id == 0)
			return i;
		if (i == size)
			break;
	}

	if (size - prev >= SIZE_1MB) {
		__VERBOSE("\n");
		__ERROR("Read is too long from file: %s", data->name);
	}

	return 0;
}

static void load_buffer(struct gb_file_inf *data)
{
	// int n_byte = gzread(data->fi, data->buf + data->size, BUF_SIZE);
	int size = SIZE_1MB - data->size;
	int n_byte = gzread(data->fi, data->buf + data->size, size);
	data->size += n_byte;
	// data->processed = gzoffset(data->fi);
	// if (n_byte < BUF_SIZE)
	// 	data->is_eof = 1;
	if (n_byte < size)
		data->is_eof = 1;
}

static void split_buffer(struct gb_file_inf *data, char *old_buf, int prev)
{
	int padding = data->size - prev;
	/* BUF_SIZE + 1 for case end of file is not '/n' */
	// data->buf = malloc(padding + BUF_SIZE + 1);
	// data->buf = kk
	memcpy(data->buf, old_buf + prev, padding);
	data->size = padding;
	data->buf[padding] = '\0';
	old_buf[prev] = '\0';
}

int gb_get_pair(struct gb_pair_data *data, char **buf1, char **buf2)
{
	// buf1, buf2 must first be initilized
	// *buf1 = *buf2 = NULL;
	if (data->finish_flag)
		return -1;

	int prev1, prev2, new_pos1, new_pos2;
	int cnt, ret;

	load_buffer(&data->file1);
	load_buffer(&data->file2);
	prev1 = prev2 = 0;

	cnt = 0;
	while (1) {
		new_pos1 = get_nxt_pos(&data->file1, prev1, data->type);
		new_pos2 = get_nxt_pos(&data->file2, prev2, data->type);
		if (!new_pos1 || !new_pos2) 
			break;
		++cnt;
		prev1 = new_pos1;
		prev2 = new_pos2;
	}

	/* no more read to get */
	if (!prev1 && !new_pos1 && !new_pos2) {
		data->finish_flag = 1;
		return -1;
	}

	/* read of two files are not equal */
	if (!prev1) {
		data->warning_flag = 1;
		data->finish_flag = 1;
		return -1;
	}

	/* split complete buffer */
	__SWAP(*buf1, data->file1.buf);
	__SWAP(*buf2, data->file2.buf);
	// *buf1 = data->file1.buf;
	// *buf2 = data->file2.buf;
	split_buffer(&data->file1, *buf1, prev1);
	split_buffer(&data->file2, *buf2, prev2);

	ret = data->offset;
	data->offset += cnt;
	return ret;
}

int gb_get_single(struct gb_single_data *data, char **buf)
{
	// *buf = NULL;
	if (data->finish_flag)
		return -1;

	int prev, new_pos;
	int cnt, ret;

	load_buffer(&data->file);
	prev = 0;

	cnt = 0;
	while (1) {
		new_pos = get_nxt_pos(&data->file, prev, data->type);
		if (!new_pos) 
			break;
		++cnt;
		prev = new_pos;
	}

	/* no more read to get */
	if (!prev) {
		data->finish_flag = 1;
		return -1;
	}

	/* split complete buffer */
	__SWAP(*buf, data->file.buf);
	// *buf = data->file.buf;
	split_buffer(&data->file, *buf, prev);
	ret = data->offset;
	data->offset += cnt;
	return ret;
}

void read_destroy(struct read_t *read, int is_buf)
{
	if (!is_buf) {
		free(read->seq);
		free(read->qual);
		free(read->name);
		free(read->note);
	}
	free(read->rseq);
	free(read->rqual);
}

int get_read_from_fq(struct read_t *read, char *buf, int *pos)
{
	int i = *pos, prev, k;

	/* name part */
	prev = i;
	if (buf[i] != '@')
		return READ_FAIL;
	read->info = NULL;
	read->name = buf + prev + 1; /* skip @ character */
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i) {
		if (__is_sep(buf[i]) && read->info == NULL) {
			buf[i] = '\0';
			read->info = buf + i + 1;
			k = i - 2;
		}
	}
	if (buf[i] == '\0')
		return READ_FAIL;
	if (read->info == NULL)
		k = i - 2;

	if (k > 0 &&
	    (strncmp(buf + k, "/1", 2) == 0 || strncmp(buf + k, "/2", 2) == 0)) {
		buf[k] = '\0';
	}

	buf[i++] = '\0';

	/* seq part */
	prev = i;
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i);
	if (buf[i] == '\0')
		return READ_FAIL;
	read->seq = buf + prev;
	read->len = i - prev;
	buf[i++] = '\0';
	if (read->len == 0)
		return READ_FAIL;

	/* optionally part */
	prev = i;
	if (buf[i] != '+')
		return READ_FAIL;
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i);
	if (buf[i] == '\0')
		return READ_FAIL;
	read->note = buf + prev;
	buf[i++] = '\0';

	/* quality part */
	prev = i;
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i);
	if (i - prev != read->len)
		return READ_FAIL;
	read->qual = buf + prev;
	if (buf[i] == '\0')
		return READ_END;
	buf[i++] = '\0';
	if (buf[i] == '\0')
		return READ_END;

	*pos = i;
	return READ_SUCCESS;
}

int get_read_from_fa(struct read_t *read, char *buf, int *pos)
{
	int i = *pos, prev, k;
	read->qual = read->note = NULL;

	/* name part */
	prev = i;
	if (buf[i] != '>')
		return READ_FAIL;
	read->info = NULL;
	read->name = buf + prev + 1; /* skip > character */
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i) {
		if (__is_sep(buf[i]) && read->info == NULL) {
			buf[i] = '\0';
			read->info = buf + i + 1;
			k = i - 2;
		}
	}
	if (buf[i] == '\0')
		return READ_FAIL;
	if (read->info == NULL)
		k = i - 2;

	if (k > 0 &&
	    (strncmp(buf + k, "/1", 2) == 0 || strncmp(buf + k, "/2", 2) == 0)) {
		buf[k] = '\0';
	}

	buf[i++] = '\0';

	/* seq part */
	prev = i;
	for (; buf[i] != '\0' && buf[i] != '\n'; ++i);
	read->seq = buf + prev;
	read->len = i - prev;
	if (read->len == 0)
		return READ_FAIL;
	if (buf[i] == '\0')
		return READ_END;
	buf[i++] = '\0';
	if (buf[i] == '\0')
		return READ_END;

	*pos = i;
	return READ_SUCCESS;
}
