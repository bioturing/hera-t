#include <zlib.h>
#include "attribute.h"
#include "bwt.h"
#include "hash_table.h"
#include "index.h"
#include "io_utils.h"
#include "kseq.h"
#include "khash.h"
#include "opt.h"
#include "radix_sort.h"
#include "utils.h"
#include "verbose.h"

#define ex_get_block(p, s, mask) ((p).beg >> (s) & (mask))
#define ex_less_than(x, y) ((x).beg < (y).beg)

RS_IMPL(exon, struct exon_t, 32, 8, ex_less_than, ex_get_block)

KHASH_MAP_INIT_STR(str_int, int)

KSEQ_INIT(gzFile, gzread)

static struct genome_info_t *genome;
static struct gene_info_t *genes;
static struct transcript_info_t *trans;
static khash_t(str_int) *chr_map, *gene_map, *tran_map;
static uint32_t *chr_total_len;

void load_genome_fasta(const char *path)
{
	extern struct genome_info_t *genome;
	extern khash_t(str_int) *chr_map;
	gzFile fp = gzopen(path, "r");
	kseq_t *seq;
	bioint_t *chr_len;
	khint_t k;
	int name_len, n_chr, hash_ret, i;
	char **chr_name;

	seq = kseq_init(fp);
	name_len = 0;
	n_chr = 0;
	chr_name = NULL;
	chr_len = NULL;

	while (kseq_read(seq) >= 0) {
		char *seq_name = strdup(seq->name.s);
		k = kh_put(str_int, chr_map, seq_name, &hash_ret);
		if (hash_ret == 0) { // Duplicated name
			__VERBOSE_LOG("WARNING", "Duplicated sequence name [%s] in genome fasta\n",
					chr_name[n_chr]);
			free(seq_name);
			continue;
		}
		chr_name = realloc(chr_name, (n_chr + 1) * sizeof(char *));
		chr_len = realloc(chr_len, (n_chr + 1) * sizeof(bioint_t));
		genome->seq = realloc(genome->seq, (n_chr + 1) * sizeof(char *));
		genome->seq[n_chr] = malloc(seq->seq.l + 1);
		memcpy(genome->seq[n_chr], seq->seq.s, seq->seq.l);
		genome->seq[n_chr][seq->seq.l] = '\0';
		chr_name[n_chr] = seq_name;
		chr_len[n_chr] = seq->seq.l;
		kh_value(chr_map, k) = n_chr;
		name_len = __max(name_len, (int)seq->name.l + 1);
		++n_chr;
	}

	// index total length
	chr_total_len = calloc(n_chr, sizeof(uint32_t));
	for (i = 1; i < n_chr; ++i)
		chr_total_len[i] = (uint32_t)chr_len[i - 1] + chr_total_len[i - 1];

	genome->n = n_chr;
	genome->chr_len = chr_len;
	genome->l_name = name_len;
	genome->chr_name = calloc(n_chr, name_len);
	for (i = 0; i < n_chr; ++i)
		strcpy(genome->chr_name + i * name_len, chr_name[i]);

	free(chr_name);
	free(chr_len);
	kseq_destroy(seq);
	gzclose(fp);
}

struct gz_reader_t {
	gzFile fp;
	char *buf;
	size_t buf_cap;
	size_t buf_len;
	size_t cur;
	int is_eof;
};

struct gz_reader_t *init_gz_reader(const char *path, size_t cap)
{
	struct gz_reader_t *ret;
	ret = calloc(1, sizeof(struct gz_reader_t));
	ret->fp = gzopen(path, "r");
	ret->buf = malloc(cap);
	ret->buf_cap = cap;
	ret->buf_len = gzread(ret->fp, ret->buf, ret->buf_cap);
	ret->cur = 0;
	ret->is_eof = 0;
	return ret;
}

char *gz_reader_get_line(struct gz_reader_t *fp)
{
	if (fp->is_eof) return NULL;
	size_t cur, len;
	cur = fp->cur;
	len = fp->buf_len;
	char *buf = fp->buf;
	while (buf[cur] != '\n') {
		while (cur < len && buf[cur] != '\n')
			++cur;
		if (cur == len) {
			memmove(buf, buf + fp->cur, fp->buf_len - fp->cur);
			fp->buf_len -= fp->cur;
			cur = fp->buf_len;
			fp->cur = 0;
			if (fp->buf_cap - fp->buf_len)
				fp->buf_len += gzread(fp->fp, buf + fp->buf_len,
						     fp->buf_cap - fp->buf_len);
			len = fp->buf_len;
			if (cur == len) { // No newline until EOF or out of buffer
				if (len == fp->buf_cap) { // out of buffer
					__ERROR("Out of buffer, line is too long (over [%lu])\n",
						fp->buf_cap);
				} else {
					__VERBOSE_LOG("WARNING", "Encountered EOF before newline\n");
					fp->is_eof = 1;
					buf[cur] = '\n';
				}
			}
		}
	}
	buf[cur] = '\0';
	char *ret = buf + fp->cur;
	fp->cur = cur + 1;
	return ret;
}

void close_gz_reader(struct gz_reader_t *fp)
{
	gzclose(fp->fp);
	free(fp->buf);
	free(fp);
}

char **split_string(char *str, const char *tok, int *cnt)
{
	char **ret;
	int n, i, flag;
	n = 0;
	ret = NULL;
	flag = 1;
	for (i = 0; str[i]; ++i) {
		if (strchr(tok, str[i])) {
			str[i] = '\0';
			flag = 1;
		} else {
			if (flag) {
				ret = realloc(ret, (n + 1) * sizeof(char *));
				ret[n++] = str + i;
				flag = 0;
			}
		}
	}
	*cnt = n;
	return ret;
}

static inline char *get_token(char **a, int n, const char *w)
{
	int i;
	for (i = 0; i < n; ++i) {
		if (!strcmp(a[i], w)) {
			if (i + 1 == n)
				return NULL;
			return a[i + 1];
		}
	}
	return NULL;
}

char *normalize_str(const char *str)
{
	int len, i;
	len = strlen(str);
	char *ret = malloc(len + 1);
	len = 0;
	for (i = 0; str[i]; ++i)
		if (str[i] >= 0x20 && str[i] <= 0x7e && str[i] != 0x22 && str[i] != 0x27)
			ret[len++] = str[i];
	ret[len] = '\0';
	ret = realloc(ret, len + 1);
	return ret;
}

void parse_gtf(const char *path, const char *tree_path)
{
	extern struct genome_info_t *genome;
	extern struct gene_info_t *genes;
	extern struct transcript_info_t *trans;
	extern khash_t(str_int) *chr_map;
	extern khash_t(str_int) *gene_map;
	extern khash_t(str_int) *tran_map;
	genes->l_name = genes->l_id = trans->l_id = 0;
	genes->n = trans->n = 0;
	struct gz_reader_t *fp;
	fp = init_gz_reader(path, SIZE_16MB);
	char **w, **ph, *str, *p_gene_id, *p_gene_name, *p_tran_id;
	char **trans_id, **genes_id, **genes_name, *tran_id, *gene_id, *gene_name;
	struct leaf_t *leaf;
	int line, n_w, n_ph, hash_ret;
	int chr_idx, gene_idx, tran_idx, cur_gene, cur_tran;
	int m_trans, m_genes, l_id, l_name;
	char strand;
	bioint_t beg, end;
	khint_t k;
	m_trans = m_genes = 0x10000;
	trans_id = malloc(m_trans * sizeof(char *));
	trans->gene_idx = malloc(m_trans * sizeof(int));
	trans->n_exon = malloc(m_trans * sizeof(int));
	trans->exons = malloc(m_trans * sizeof(struct exon_t *));

	genes_id = malloc(m_genes * sizeof(char *));
	genes_name = malloc(m_genes * sizeof(char *));
	genes->strand = malloc(m_genes * sizeof(char));
	genes->chr_idx = malloc(m_genes * sizeof(int));
	leaf = malloc(m_genes * sizeof(struct leaf_t));

	cur_tran = cur_gene = -1;

	line = 0;
	while ((str = gz_reader_get_line(fp)) != NULL) {
		// fprintf(stderr, "%s\n", str);
		++line;
		if (str[0] == '#')
			continue;
		w = split_string(str, "\t", &n_w);
		if (!w) {
			__VERBOSE_LOG("WARNING",
					"Line [%d]: Ignore blank line\n", line);
			free(w);
			continue;
		}
		if (n_w != 9)
			__ERROR("Line [%d]: Number of fields [%d] is not equal to 9",
				line, n_w);
		k = kh_get(str_int, chr_map, w[0]);
		if (k == kh_end(chr_map)) {
			__VERBOSE_LOG("WARNING", "Line [%d]: Sequence [%s] does not exist in fasta file\n",
					line, w[0]);
			free(w);
			continue;
		}
		chr_idx = kh_value(chr_map, k);

		beg = atoi(w[3]);
		end = atoi(w[4]);
		if (beg > end)
			__ERROR("Line [%d]: Start position [%s] > End position [%s]",
				line, w[3], w[4]);
		strand = strcmp(w[6], "+") == 0 ? 0 : 1;

		ph = split_string(w[8], " ;=", &n_ph);

		if (strcmp(w[2], "gene") && strcmp(w[2], "transcript") &&
		    strcmp(w[2], "exon")) {
			free(ph);
			free(w);
			continue;
		}

		p_gene_id = get_token(ph, n_ph, "gene_id");
		p_gene_name = get_token(ph, n_ph, "gene_name");
		p_tran_id = get_token(ph, n_ph, "transcript_id");

		if (!strcmp(w[2], "exon")) {
			if (p_tran_id == NULL) {
				if (cur_tran == -1)
					__ERROR("Line [%d]: No transcript candidate for exon",
						line);
				__VERBOSE_LOG("WARNING", "Line [%d]: Using transcript id in current transcript level\n",
						line);
				tran_idx = cur_tran;
			} else {
				tran_id = normalize_str(p_tran_id);
				k = kh_get(str_int, tran_map, tran_id);
				if (k == kh_end(tran_map)) {
					__ERROR("Line [%d]: No information about transcript [%s] was declared",
						line, tran_id);
				} else {
					tran_idx = kh_value(tran_map, k);
					free(tran_id); tran_id = NULL;
				}
			}
			if (tran_idx != cur_tran) {
				__VERBOSE_LOG("WARNING", "Line [%d]: Exon is not nested in transcript level\n",
						line);
				cur_tran = tran_idx;
			}

			gene_idx = trans->gene_idx[tran_idx];
			if (strand != genes->strand[gene_idx])
				__ERROR("Line [%d]: Exon strand is not consistent with transcript strand",
					line);
			if (chr_idx != genes->chr_idx[gene_idx])
				__ERROR("Line [%d]: Exon is not on the same chromosome with the transcript", line);
			trans->exons[tran_idx] = realloc(trans->exons[tran_idx],
				(trans->n_exon[tran_idx] + 1) * sizeof(struct exon_t));
			trans->exons[tran_idx][trans->n_exon[tran_idx]].beg = beg;
			trans->exons[tran_idx][trans->n_exon[tran_idx]].end = end;
			++trans->n_exon[tran_idx];
		} else if (!strcmp(w[2], "transcript")) {
			if (p_gene_id == NULL) {
				if (cur_gene == -1)
					__ERROR("Line [%d]: No gene candidate for transcript",
						line);
				__VERBOSE_LOG("WARNING", "Line [%d]: Gene id is not found in attribute field, using gene id in current gene level\n",
						line);
				gene_idx = cur_gene;
			} else {
				gene_id = normalize_str(p_gene_id);
				k = kh_get(str_int, gene_map, gene_id);
				if (k == kh_end(gene_map)) {
					__ERROR("Line [%d]: No information about gene [%s] was declared",
						line, gene_id);
				} else {
					gene_idx = kh_value(gene_map, k);
					free(gene_id); gene_id = NULL;
				}
			}
			if (gene_idx != cur_gene) {
				__VERBOSE_LOG("WARNING", "Line [%d]: Transcript is not nested in gene level\n",
						line);
				cur_gene = gene_idx;
			}

			if (strand != genes->strand[gene_idx])
				__ERROR("Line [%d]: Transcript strand is not consistent with gene strand",
					line);
			if (chr_idx != genes->chr_idx[gene_idx])
				__ERROR("Line [%d]: Transcript is not on the same chromosome with gene",
					line);
			if (p_tran_id == NULL)
				__ERROR("Line [%d]: Transcript id is not found in attribute field",
					line);

			tran_id = normalize_str(p_tran_id);
			k = kh_put(str_int, tran_map, tran_id, &hash_ret);
			if (hash_ret == 0) {
				__VERBOSE_LOG("WARNING",
					"Line [%d]: Duplicated tran id\n", line);
				free(tran_id); tran_id = NULL;
				free(w); w = NULL;
				free(ph); ph = NULL;
				continue;
			}
			kh_value(tran_map, k) = trans->n;
			if (trans->n == m_trans) {
				m_trans <<= 1;
				trans_id = realloc(trans_id,
						m_trans * sizeof(char *));
				trans->gene_idx = realloc(trans->gene_idx,
							m_trans * sizeof(int));
				trans->n_exon = realloc(trans->n_exon,
							m_trans * sizeof(int));
				trans->exons = realloc(trans->exons,
					       m_trans * sizeof(struct exon_t));
			}
			trans_id[trans->n] = tran_id;
			trans->gene_idx[trans->n] = gene_idx;
			trans->n_exon[trans->n] = 0;
			trans->exons[trans->n] = NULL;
			l_id = strlen(tran_id) + 1;
			trans->l_id = __max(trans->l_id, l_id);
			cur_tran = trans->n;
			++trans->n;

		} else if (!strcmp(w[2], "gene")) {
			if (p_gene_id == NULL)
				__ERROR("Line [%d]: Gene id is not found in attribute field",
					line);
			gene_id = normalize_str(p_gene_id);
			if (p_gene_name == NULL) {
				__VERBOSE_LOG("WARNING",
						"Line [%d]: Gene name is not found, using gene id instead\n",
							line);
				p_gene_name = p_gene_id;
			}
			gene_name = normalize_str(p_gene_name);
			k = kh_put(str_int, gene_map, gene_id, &hash_ret);
			if (hash_ret == 0) {
				__VERBOSE_LOG("WARNING",
					"Line [%d]: Duplicated gene id\n",
						line);
				free(gene_id); gene_id = NULL;
				free(gene_name); gene_name = NULL;
				free(w); w = NULL;
				free(ph); ph = NULL;
				continue;
			}
			kh_value(gene_map, k) = genes->n;
			if (genes->n == m_genes) {
				m_genes <<= 1;
				genes_id = realloc(genes_id,
						m_genes * sizeof(char *));
				genes_name = realloc(genes_name,
						m_genes * sizeof(char *));
				genes->strand = realloc(genes->strand,
						m_genes * sizeof(char));
				genes->chr_idx = realloc(genes->chr_idx,
						m_genes * sizeof(int));
				leaf = realloc(leaf,
						m_genes * sizeof(struct leaf_t));
			}
			genes_id[genes->n] = gene_id;
			genes_name[genes->n] = gene_name;
			genes->strand[genes->n] = strand;
			genes->chr_idx[genes->n] = chr_idx;
			leaf[genes->n].start = beg + chr_total_len[chr_idx];
			leaf[genes->n].end = end + chr_total_len[chr_idx];
			l_name = strlen(gene_name) + 1;
			l_id = strlen(gene_id) + 1;
			genes->l_name = __max(genes->l_name, l_name);
			genes->l_id = __max(genes->l_id, l_id);
			cur_gene = genes->n;
			++genes->n;
		} else {
			// Shall never happen
			assert(0);
		}
	}
	genes->gene_name = calloc(genes->n, genes->l_name);
	genes->gene_id = calloc(genes->n, genes->l_id);
	int i;
	for (i = 0; i < genes->n; ++i) {
		strcpy(genes->gene_name + i * genes->l_name, genes_name[i]);
		strcpy(genes->gene_id + i * genes->l_id, genes_id[i]);
		free(genes_name[i]);
		// Have not freed genes_id[i] yet
	}
	free(genes_name);
	free(genes_id);

	trans->tran_id = calloc(trans->n, trans->l_id);
	for (i = 0; i < trans->n; ++i) {
		strcpy(trans->tran_id + i * trans->l_id, trans_id[i]);
		// Have not freed trans_id[i] yet
	}
	free(trans_id);
	close_gz_reader(fp);

	// Build gene interval tree
	// build_interval_tree(tree_path, leaf, genes->n);
	free(leaf);
}

void build_transcript(const char *path)
{
	FILE *fp = xfopen(path, "wb");
	bioint_t kb;
	int m_seq, l, len, i, j, k, n;
	int gene_idx, chr_idx;
	char c, strand;
	char *seq, *tran_id, *gene_id;
	trans->tran_len = malloc(trans->n * sizeof(int));
	trans->tran_beg = malloc((trans->n + 1) * sizeof(int));
	m_seq = 0x10000;
	trans->idx = malloc(m_seq * sizeof(int));
	trans->seq = malloc(m_seq);
	l = 0;
	for (i = 0; i < trans->n; ++i) {
		trans->tran_beg[i] = l;
		n = trans->n_exon[i];
		gene_idx = trans->gene_idx[i];
		tran_id = trans->tran_id + i * trans->l_id;
		gene_id = genes->gene_id + gene_idx * genes->l_id;
		if (n == 0) {
			__VERBOSE_LOG("WARNING", "Transcript [%s] have 0 exon\n", tran_id);
			continue;
		}
		rs_sort(exon, trans->exons[i], trans->exons[i] + n);
		strand = genes->strand[gene_idx];
		chr_idx = genes->chr_idx[gene_idx];
		len = 0;
		for (k = 0; k < n; ++k) {
			if (k + 1 < n &&
			    trans->exons[i][k].end >= trans->exons[i][k + 1].beg)
				__ERROR("Transcript [%s] has exons overlap", tran_id);
			len += (trans->exons[i][k].end - trans->exons[i][k].beg + 1);
		}
		trans->tran_len[i] = len;
		if (l + len + 1 > m_seq) {
			m_seq = l + len + 1;
			__round_up_32(m_seq);
			trans->idx = realloc(trans->idx, m_seq * sizeof(int));
			trans->seq = realloc(trans->seq, m_seq);
		}
		seq = trans->seq + l;
		len = 0;
		if (strand == 0) {
			for (j = 0; j < n; ++j) {
				memcpy(seq + len, genome->seq[chr_idx] + (trans->exons[i][j].beg - 1),
					(trans->exons[i][j].end - trans->exons[i][j].beg + 1));
				len += (trans->exons[i][j].end - trans->exons[i][j].beg + 1);
			}
		} else {
			for (j = n - 1; j >= 0; --j) {
				for (kb = trans->exons[i][j].end; kb >= trans->exons[i][j].beg; --kb) {
					c = genome->seq[chr_idx][kb - 1];
					seq[len++] = rev_nt4_char[nt4_table[(int)c]];
				}
			}
		}
		for (k = 0; k < len; ++k)
			trans->idx[l + k] = i;
		l += len;
		seq[len] = '\0';
		fprintf(fp, ">%s_%s\n", gene_id, tran_id);
		while (len > 0) {
			fprintf(fp, "%.80s\n", seq);
			seq += 80;
			len -= 80;
		}
	}
	trans->idx = realloc(trans->idx, l * sizeof(int));
	trans->seq = realloc(trans->seq, l);
	trans->tran_beg[trans->n] = l;
	xwfclose(fp);
}

void construct_hash(int k_s)
{
	init_cons_hash(27);
	uint64_t kmer;
	int i, k, len, c, last, cnt_kmer = 0, cnt_pos = 0;
	char *seq;
	uint64_t mask = (uint64_t)((1ull << (k_s << 1)) - 1);
	
	for (i = 0; i < trans->n; ++i) {
		seq = trans->seq + trans->tran_beg[i];
		len = trans->tran_len[i];
		kmer = 0;
		last = 0;
		for (k = 0; k < len; ++k) {
			c = nt4_table[(int)seq[k]];
			kmer = (kmer << 2) & mask;
			if (c < 4) {
				kmer |= c;
				++last;
			} else {
				last = 0;
			}
			if (last >= k_s) {
				++cnt_pos;
				cnt_kmer += insert_cons_hash(kmer);
			}
		}
	}
	__VERBOSE_LOG("INFO", "Number of kmer: %d\n", cnt_kmer);
	__VERBOSE_LOG("INFO", "Number of hashing position: %d\n", cnt_pos);

	recount_cons_hash();
	for (i = 0; i < trans->n; ++i) {
		seq = trans->seq + trans->tran_beg[i];
		len = trans->tran_len[i];
		kmer = 0;
		last = 0;
		for (k = 0; k < len; ++k) {
			c = nt4_table[(int)seq[k]];
			kmer = (kmer << 2) & mask;
			if (c < 4) {
				kmer |= c;
				++last;
			} else {
				last = 0;
			}
			if (last >= k_s)
				addpos_cons_hash(kmer, k - k_s + 1 + trans->tran_beg[i]);
		}
	}
}

void dump_info(const char *path)
{
	FILE *fp = xfopen(path, "wb");
	xfwrite(&genome->n, sizeof(int), 1, fp);
	xfwrite(&genome->l_name, sizeof(int), 1, fp);
	xfwrite(genome->chr_len, sizeof(bioint_t), genome->n, fp);
	xfwrite(genome->chr_name, genome->l_name, genome->n, fp);

	xfwrite(&genes->n, sizeof(int), 1, fp);
	xfwrite(&genes->l_name, sizeof(int), 1, fp);
	xfwrite(&genes->l_id, sizeof(int), 1, fp);
	xfwrite(genes->chr_idx, sizeof(int), genes->n, fp);
	xfwrite(genes->gene_name, genes->l_name, genes->n, fp);
	xfwrite(genes->gene_id, genes->l_id, genes->n, fp);
	xfwrite(genes->strand, sizeof(char), genes->n, fp);

	xfwrite(&trans->n, sizeof(int), 1, fp);
	xfwrite(&trans->l_id, sizeof(int), 1, fp);
	xfwrite(trans->tran_id, trans->l_id, trans->n, fp);
	xfwrite(trans->tran_len, sizeof(int), trans->n, fp);
	xfwrite(trans->tran_beg, sizeof(int), trans->n + 1, fp);
	xfwrite(trans->gene_idx, sizeof(int), trans->n, fp);
	xfwrite(trans->idx, sizeof(int), trans->tran_beg[trans->n], fp);
	xfwrite(trans->seq, sizeof(char), trans->tran_beg[trans->n], fp);

	xfwrite(trans->n_exon, sizeof(int), trans->n, fp);
	int i;
	for (i = 0; i < trans->n; ++i)
		xfwrite(trans->exons[i], sizeof(struct exon_t), trans->n_exon[i], fp);
	xwfclose(fp);
}

void free_info()
{
	//  TODO: free info before build hash
}

void build_index(int pos, int argc, char **argv)
{
	char str_dir[1024];
	char idx_name[1024];
	struct opt_index_t *opts = get_opt_index(argc - pos, argv + pos);
	strcpy(idx_name, opts->idx_dir); strcat(idx_name, "/");
	strcat(idx_name, opts->prefix);
	// idx_name = str_concate(opts->idx_dir, opts->prefix);

	strcpy(str_dir, idx_name); strcat(str_dir, ".log");
	init_log(str_dir);
	log_write("INDEX VERSION: %d\n", INDEX_VERSION);
	log_write("COMMAND: ");
	int i;
	for (i = 0; i < argc; ++i)
		log_write("%s ", argv[i]);
	log_write("\n");

	// Build bwt
	if (opts->bwt) {
		__VERBOSE_INFO("INFO", "Building Burrow-Wheeler Transform on genome...\n");
		struct bwt_t *bwt;
		bwt = bwt_build_from_fasta(opts->genome);
		strcpy(str_dir, idx_name); strcat(str_dir, ".bwt");
		bwt_dump(str_dir, bwt);
		bwt_destroy(bwt);
	} else {
		__VERBOSE_INFO("INFO", "Not rebuild BWT\n");
	}

	// Load fasta and gtf
	__VERBOSE_INFO("INFO", "Loading fasta...\n");
	genome = calloc(1, sizeof(struct genome_info_t));
	chr_map = kh_init(str_int);
	load_genome_fasta(opts->genome);

	__VERBOSE_INFO("INFO", "Loading gtf...\n");
	genes = calloc(1, sizeof(struct gene_info_t));
	trans = calloc(1, sizeof(struct transcript_info_t));
	tran_map = kh_init(str_int);
	gene_map = kh_init(str_int);

	strcpy(str_dir, idx_name);
	strcat(str_dir, ".tree");
	parse_gtf(opts->gtf, str_dir);

	strcpy(str_dir, idx_name);
	strcat(str_dir, ".fasta");
	build_transcript(str_dir);

	strcpy(str_dir, idx_name);
	strcat(str_dir, ".info");
	dump_info(str_dir);
	free_info();

	__VERBOSE_INFO("INFO", "Constructing kmer hash for transcript...\n");
	construct_hash(opts->k);
	strcpy(str_dir, idx_name); strcat(str_dir, ".hash");
	store_cons_hash(str_dir, opts->k);
}
