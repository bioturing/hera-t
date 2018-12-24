#ifndef _GET_BUFFER_H_
#define _GET_BUFFER_H_

#include <zlib.h>

#include "attribute.h"

#define TYPE_FASTQ		0
#define	TYPE_FASTA		1

// #define BUF_OK			0
// #define BUF_FAIL		1

struct gb_file_inf {
	gzFile fi;
	int size;
	int is_eof;
	char *buf;
	char *name;
	long processed;
};

/* pair */

struct gb_pair_data {
	struct gb_file_inf file1;
	struct gb_file_inf file2;
	int finish_flag;
	int warning_flag;
	int type;
	int file_id;
	int offset;
};

void gb_pair_init(struct gb_pair_data *data, char *file_path1, char *file_path2);
void gb_pair_destroy(struct gb_pair_data *data);
int gb_get_pair(struct gb_pair_data *data, char **buf1, char **buf2);

/* single */

struct gb_single_data {
	struct gb_file_inf file;
	int finish_flag;
	int type;
	int file_id;
	int offset;
};

void gb_single_init(struct gb_single_data *data, char *file_path);
void gb_single_destroy(struct gb_single_data *data);
int gb_get_single(struct gb_single_data *data, char **buf);

/* get read */

#define READ_SUCCESS		0
#define READ_END		1
#define	READ_FAIL		2

void read_destroy(struct read_t *read, int is_buf);
int get_read_from_fq(struct read_t *read, char *buf, int *pos);
int get_read_from_fa(struct read_t *read, char *buf, int *pos);

#endif /* _GET_BUFFER_H_ */
