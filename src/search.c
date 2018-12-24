#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct {
    unsigned int range[4];
    unsigned int child[2];
    unsigned int median;
} Node;

typedef struct {
    unsigned int range[2];
    unsigned int info[2];
} Leaf;

typedef struct {
    Node *node;
    Leaf *leaf;
    char *info;
    unsigned int n_leaf;
    unsigned int n_node;
    unsigned int l_info;
} Tree;

unsigned int n_chr = 25;

Tree *init()
{
    unsigned short i = 0;
    Tree *TREE = (Tree*) calloc(n_chr, sizeof(Tree));
    for (i = 0; i < n_chr; ++i){
        TREE[i].node = (Node*)calloc(1, sizeof(Node));
        TREE[i].leaf = (Leaf*)calloc(1, sizeof(Leaf));
        TREE[i].info = (char*) malloc(10000);
        TREE[i].n_node = TREE[i].n_leaf = TREE[i].l_info = 0;
    }
    return TREE;
}

void free_tree(Tree *TREE)
{
    unsigned int i;
    for (i = 0; i < n_chr; ++i){
        free(TREE[i].node);
        free(TREE[i].leaf);
        free(TREE[i].info);
    }
    free(TREE);
}

void load_file(char *path, Tree *TREE)
{
    FILE *f = fopen(path, "r");
    unsigned int chr, start, end, l;
    while (!feof(f)){
        fscanf(f, "%u\t%u\t%u\t%u\t", &chr, &start, &end, &l);
        if (chr >= n_chr){
            break;
        }

        TREE[chr].info = (char*)realloc(TREE[chr].info, TREE[chr].l_info + l + 10);
        fgets(TREE[chr].info + TREE[chr].l_info, l + 10, f);

        ++TREE[chr].n_leaf;
        TREE[chr].leaf =(Leaf*) realloc(TREE[chr].leaf, TREE[chr].n_leaf*sizeof(Leaf));

        TREE[chr].leaf[TREE[chr].n_leaf - 1].range[0] = start;
        TREE[chr].leaf[TREE[chr].n_leaf - 1].range[1] = end;
        TREE[chr].leaf[TREE[chr].n_leaf - 1].info[0] = TREE[chr].l_info;
        TREE[chr].leaf[TREE[chr].n_leaf - 1].info[1] = l;

        TREE[chr].l_info += l;
    }
}

void build_tree(unsigned int id, unsigned int start,
                        unsigned int end, unsigned int pos, Tree *TREE)
{
    unsigned int i = 0;
    unsigned int mid = (start + end) >> 1;
    unsigned int median = (TREE[id].leaf[mid].range[0] +
                            TREE[id].leaf[mid].range[1]) >> 1;
    unsigned int max = TREE[id].leaf[mid].range[1];

    TREE[id].node[pos].median = median;
    if (start == end){
        TREE[id].node[pos].range[0] = TREE[id].leaf[start].range[0];
        TREE[id].node[pos].range[1] = TREE[id].leaf[start].range[1];
        TREE[id].node[pos].range[2] = TREE[id].node[pos].range[3] = start;
        TREE[id].node[pos].child[0] = TREE[id].node[pos].child[1] = 0;
        return;
    }

    for (i = mid + 1; i > start; --i){
        if (TREE[id].leaf[i - 1].range[1] < median)
            break;
        if (TREE[id].leaf[i - 1].range[1] > max)
            max = TREE[id].leaf[i - 1].range[1];
    }
    TREE[id].node[pos].range[0] = TREE[id].leaf[i].range[0];
    TREE[id].node[pos].range[2] = i;

    for (i = mid + 1; i <= end; ++i){
        if (TREE[id].leaf[i].range[0] > median)
            break;
        if (TREE[id].leaf[i].range[1] > max)
            max = TREE[id].leaf[i].range[1];
    }
    TREE[id].node[pos].range[1] = max;
    TREE[id].node[pos].range[3] = i - 1;

    if (TREE[id].node[pos].range[2] > start){
        ++TREE[id].n_node;
        TREE[id].node =(Node*) realloc(TREE[id].node, TREE[id].n_node*sizeof(Node));
        TREE[id].node[pos].child[0] = TREE[id].n_node - 1;
        build_tree(id, start, 
                    TREE[id].node[pos].range[2] - 1, TREE[id].n_node - 1, TREE);
    } else {
        TREE[id].node[pos].child[0] = 0;
    }

    if (TREE[id].node[pos].range[3] < end){
        ++TREE[id].n_node;
        TREE[id].node = (Node*) realloc(TREE[id].node, TREE[id].n_node*sizeof(Node));
        TREE[id].node[pos].child[1] = TREE[id].n_node - 1;
        build_tree(id, TREE[id].node[pos].range[3] + 1, 
                                end, TREE[id].n_node - 1, TREE);
    } else {
        TREE[id].node[pos].child[1] = 0;
    }

    return;
}

void write_tree(char *path, Tree *TREE)
{
    FILE *f = fopen(path, "wb");
    unsigned int i = 0;
    unsigned int pos[n_chr];
    unsigned int len = n_chr*sizeof(int);

    for (i = 0; i < n_chr; ++i){
        pos[i] = len;
        len += TREE[i].n_node*sizeof(Node) + TREE[i].n_leaf*sizeof(Leaf) +
                TREE[i].l_info + 3*sizeof(int);
    }

    fwrite(pos, sizeof(int), n_chr, f);
    for (i = 0; i < n_chr; ++i){
        fwrite(&TREE[i].n_node, sizeof(int), 1, f);
        fwrite(&TREE[i].n_leaf, sizeof(int), 1, f);
        fwrite(&TREE[i].l_info, sizeof(int), 1, f);
        if (TREE[i].n_leaf > 0){
            fwrite(TREE[i].node, sizeof(Node), TREE[i].n_node, f);
            fwrite(TREE[i].leaf, sizeof(Leaf), TREE[i].n_leaf, f);
            fwrite(TREE[i].info, TREE[i].l_info, 1, f);
        }
    }
    fclose(f);
}

void load_index(char *path, unsigned int chr, Tree *TREE)
{
    FILE *f = fopen(path, "rb");
    unsigned int ret;
    unsigned int pos[n_chr];

    ret = fread(pos, sizeof(int), n_chr, f);
    fseek(f, pos[chr], SEEK_SET);

    ret = fread(&TREE[chr].n_node, sizeof(int), 1, f);
    ret = fread(&TREE[chr].n_leaf, sizeof(int), 1, f);
    ret = fread(&TREE[chr].l_info, sizeof(int), 1, f);

    if (TREE[chr].n_leaf > 0){
        TREE[chr].node = (Node*) calloc(TREE[chr].n_node, sizeof(Node));
        TREE[chr].leaf = (Leaf*) calloc(TREE[chr].n_leaf, sizeof(Leaf));
        TREE[chr].info = (char*) malloc(TREE[chr].l_info);

        ret = fread(TREE[chr].node, sizeof(Node), TREE[chr].n_node, f);
        ret = fread(TREE[chr].leaf, sizeof(Leaf), TREE[chr].n_leaf, f);
        ret = fread(TREE[chr].info, TREE[chr].l_info, 1, f);
    }
    fclose(f);
}

unsigned int in_range(unsigned int start, unsigned int end,
                        unsigned int qstart, unsigned int qend)
{
    if ((start >= qstart && start <= qend) || (end >= qstart && end <= qend))
        return 0;
    else
        return 1;
}

void query(unsigned int chr, unsigned int qstart,
                unsigned int qend, unsigned int pos, Tree *TREE)
{
    unsigned int start = TREE[chr].node[pos].range[0];
    unsigned int end = TREE[chr].node[pos].range[1];
    unsigned int i, next;

    if ((qstart < end && qend > start) || in_range(start, end, qstart, qend) == 0){
        for (i = TREE[chr].node[pos].range[2];
                i <= TREE[chr].node[pos].range[3]; ++i){
            if (in_range(TREE[chr].leaf[i].range[0], TREE[chr].leaf[i].range[1],
                            qstart, qend) == 0){
                char text[TREE[chr].leaf[i].info[1] + 1];
                memcpy(text, TREE[chr].info + TREE[chr].leaf[i].info[0],
                                             TREE[chr].leaf[i].info[1]);
                text[TREE[chr].leaf[i].info[1]] = '\0';
                printf("%u %u %u %s\n", chr, TREE[chr].leaf[i].range[0],
                                            TREE[chr].leaf[i].range[1], text);
            }
        }
    }

    next = TREE[chr].node[pos].child[0];
    if (next > 0 && qstart < TREE[chr].node[pos].median)
        query(chr, qstart, qend, next, TREE);

    next = TREE[chr].node[pos].child[1];
    if (next > 0 && qend > TREE[chr].node[pos].median)
        query(chr, qstart, qend, next, TREE);
}

int main(int argc, char *argv[])
{
    unsigned int i;

    Tree *TREE = init();
    if (strncmp(argv[1], "index", 5) == 0){
        load_file(argv[2], TREE);
        printf("loaded\n");
        for (i = 0; i < n_chr; i++){
            if (TREE[i].n_leaf == 0){
                continue;
            }
            build_tree(i, 0, TREE[i].n_leaf - 1, 0, TREE);
            printf("Builded %u tree (%u Leaf, %u Node)\n", i, TREE[i].n_leaf, TREE[i].n_node);
        }
        write_tree(argv[3], TREE);
    } else {
        i = atoi(argv[3]);
        load_index(argv[2], i, TREE);
        if (TREE[i].n_leaf > 0){
            query(i, atoi(argv[4]), atoi(argv[5]), 0, TREE);
        }
        free(TREE);
    }
    return 0;
}
