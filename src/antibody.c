#include <string.h>
#include <pthread.h>

#include "antibody.h"
#include "readTSV.h"
#include "utils.h"
#include "barcode.h"
#include "attribute.h"
#include "verbose.h"

const uint64_t _pow5_r[] = {1ull, 5ull, 25ull, 125ull, 625ull, 3125ull, 15625ull, 78125ull, 390625ull, 1953125ull, 9765625ull, 48828125ull, 244140625ull, 1220703125ull, 6103515625ull, 30517578125ull, 152587890625ull, 762939453125ull, 3814697265625ull, 19073486328125ull, 95367431640625ull, 476837158203125ull, 2384185791015625ull, 11920928955078125ull, 59604644775390625ull, 298023223876953125ull, 1490116119384765625ull};

struct antibody_lib_t *init_antibody_lib()
{
	struct antibody_lib_t *lib;

	lib = calloc(1, sizeof(struct antibody_lib_t));
	lib->whitelist_path = NULL;
	lib->n_files = 0;
	lib->left_file = lib->right_file = NULL;
	lib->trim = -1;

	return lib;
}

void destroy_antibody_lib(struct antibody_lib_t *lib)
{
        free(lib);
}

int is_defined(struct antibody_lib_t *lib)
{
        int is_defined = lib->whitelist_path != NULL;
        is_defined &= lib->left_file != NULL;
        is_defined &= lib->right_file != NULL;

        if (!is_defined)
                destroy_antibody_lib(lib);
        return is_defined;
}

struct antibody_lib_t *check_valid_cell(struct antibody_lib_t *lib)
{
        if (!is_defined(lib))
                return NULL;

	if (lib->trim < 0){
		__WARNING("--cell_trim param is not defined.\
                        No base will be trimmed off from read 2.");
                lib->trim = 0;
        }

	if (lib->whitelist_path == NULL)
		__ERROR("Missing --cell_tags argument");

	if (lib->left_file == NULL || lib->right_file == NULL)
		__ERROR("Missing input files for cell hashing");

        return lib;
}

struct antibody_lib_t *check_valid_protein(struct antibody_lib_t *lib)
{
        if (!is_defined(lib))
                return NULL;

	if (lib->trim < 0){
		__WARNING("--protein_trim param is not defined.\
                        No base will be trimmed off from read 2.");
                lib->trim = 0;
        }

	if (lib->whitelist_path == NULL)
		__ERROR("Missing --protein_tags argument");

	if (lib->left_file == NULL || lib->right_file == NULL)
		__ERROR("Missing input files for protein level quantification");

        return lib;
}

struct reference_t *init_reference()
{
	struct reference_t *ref;

	ref = calloc(1, sizeof(struct reference_t));
	ref->h = kh_init(bc);
	ref->n_ref = ref->ref_len = 0;
	ref->name_idx = calloc(1, sizeof(int));
	ref->ref_name = malloc(1);
	ref->ref_id = malloc(1);

	return ref;
};

void destroy_reference(struct reference_t *ref)
{
	kh_destroy(bc, ref->h);
	free(ref->name_idx);
	free(ref->ref_name);
	free(ref->ref_id);
	free(ref);
}

void check_valid_ref(struct content_t ref)
{
	if (!check_valid_nu(ref.s, ref.l))
		__ERROR("Reference tags sequence is not valid %s", ref.s);
}

struct content_t add_ref_name(int id)
{
	struct content_t name;
	char n[100];

	sprintf(n, "ref_%d", id);
	name.s = n;
	name.l = strlen(n);

	return name;
}

void add_ref_hash(struct reference_t *ref, int32_t idx)
{
	khash_t(bc) *h;
	khiter_t k;
	int32_t ret, i, tmp_idx, new_idx, ch, c;

	h = ref->h;
	k = kh_get(bc, h, idx);

	if (k != kh_end(h))
		__ERROR("There is pair of reference tags that is 1-hamming-distance away from each other");

	k = kh_put(bc, h, idx, &ret);
	kh_value(h, k) = ref->n_ref - 1;

	for (i = 0; i < ref->ref_len; ++i) {
		ch = (idx / _pow5_r[i]) % 5;
		tmp_idx = idx - ch * _pow5_r[i];
		for (c = 0; c < NNU; ++c) {
			if (ch == c)
				continue;
			new_idx = tmp_idx + _pow5_r[i] * c;

			k = kh_put(bc, h, new_idx, &ret);
			kh_value(h, k) = ref->n_ref - 1;
		}
	}
}

void add_reference(struct reference_t *ref, struct content_t id,
					 struct content_t name)
{
	if (ref->ref_len > 0 && id.l != ref->ref_len)
		__ERROR("Reference tags are not in the same length.");

	check_valid_ref(id);

	if (name.l == 0)
		name = add_ref_name(ref->n_ref);

	int start;

	// Save id	
	ref->ref_len = id.l;
	start = ref->n_ref * ref->ref_len;
	ref->ref_id = realloc(ref->ref_id, start + id.l);
	memcpy(ref->ref_id + start, id.s, id.l);

	// Save name
	start = ref->n_ref == 0? 0: ref->name_idx[ref->n_ref];
	ref->n_ref += 1;
	ref->ref_name = realloc(ref->ref_name, start + name.l);
	memcpy(ref->ref_name + start, name.s, name.l);

	ref->name_idx = realloc(ref->name_idx, (ref->n_ref + 1) * sizeof(int));
	ref->name_idx[ref->n_ref] = start + name.l;

	// Add to hash
	add_ref_hash(ref, seq2num(id.s, id.l));
}

void build_reference(struct antibody_lib_t *lib)
{
	struct tsv_t *f = init_readTSV(lib->whitelist_path);
	struct reference_t *ref = init_reference();

	while(get_row_content(f) > 0)
		add_reference(ref, get_col_content(f, 0), get_col_content(f, 1));

	destroy_readTSV(f);

	lib->ref = ref;
}

int64_t correct_bc(struct kmhash_t *h, int64_t idx, int len)
{

	int32_t i, k;
	int64_t tmp_idx, new_idx, ch, c;

	k = kmhash_get(h, idx);
	if (k != KMHASH_MAX_SIZE)
		return idx;

	for (i = 0; i < len; ++i) {
		ch = (idx / _pow5_r[i]) % 5;
		tmp_idx = idx - ch * _pow5_r[i];
		for (c = 0; c < NNU; ++c) {
			if (ch == c)
				continue;
			new_idx = tmp_idx + _pow5_r[i] * c;

			k = kmhash_get(h, new_idx);
			if (k != KMHASH_MAX_SIZE)
				return new_idx;
		}
	}

	return -1;	
}

int32_t get_tag_idx(struct antibody_lib_t *lib, struct read_t *read)
{
	if (lib->trim > read->len)
		__ERROR("Number of trimmed based is larger than read length");

	int32_t idx;
	khiter_t k;
	khash_t(bc) *h;

	h = lib->ref->h;
	idx = seq2num(read->seq + lib->trim, lib->ref->ref_len);
	k = kh_get(bc, h, idx);

	if (k == kh_end(h))
		return -1;
	return kh_value(h, k);
}

void map_antibody_read(struct read_t *read1, struct read_t *read2,
			 struct worker_bundle_t *bundle)
{
	if (!read1->name || !read1->seq || !read2->name || !read2->seq)
		return;

	int bc_len, umi_len, r1_len;
	int64_t bc_idx, umi_idx;
	int32_t tag_idx;
	struct kmhash_t *h;

	bc_len = bundle->lib.bc_len;
	umi_len = bundle->lib.umi_len;
	r1_len = bc_len + umi_len;
	if (read1->len != r1_len)
		__ERROR("Read lenght of %s is not consistent with library type.\n Expect %u.\n Receive %u.\n", read1->name, r1_len, read1->len);

	++bundle->result->nread;

	h = bundle->bc_table;
	bc_idx = seq2num(read1->seq, bc_len);
	bc_idx = correct_bc(h, bc_idx, bc_len);

	if (bc_idx == -1)
		return;

	tag_idx = get_tag_idx(bundle->antibody_lib, read2);
	if (tag_idx == -1)
		return;

	umi_idx = seq2num(read1->seq + bc_len, umi_len);
	umi_idx = umi_idx << GENE_BIT_LEN | tag_idx;
	kmhash_put_bc_umi(h, bc_idx, umi_idx);
}

void print_ref(const char *out_dir, struct reference_t *ref)
{
	int i, start, len;
	char id[ref->ref_len + 1];
	char out_path[MAX_PATH];
	FILE *fp;

	strcpy(out_path, out_dir);
	strcat(out_path, "/tags.tsv");

	fp = xfopen(out_path, "w");
	for (i = 0; i < ref->n_ref; ++i) {
		len = ref->ref_len;
		start = i * ref->ref_len;
		memcpy(id, ref->ref_id + start, len);
		id[len] = '\0';

		start = ref->name_idx[i];
		len = ref->name_idx[i + 1] - start;
		char name[len + 1];
		memcpy(name, ref->ref_name + start, len);
		name[len] = '\0';
		fprintf(fp, "%s\t%s\n", id, name);
	}

	fclose(fp);
}
