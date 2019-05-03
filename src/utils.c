#include "utils.h"

int8_t nt4_table[256] = {
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 0, 4, 1,   4, 4, 4, 2,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   3, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4, 
	4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4,   4, 4, 4, 4
};

char *nt4_char = "ACGTN";
char *rev_nt4_char = "TGCAN";

char *str_concate(const char *str1, const char *str2)
{
	/* TODO: put this to depricated */
	size_t len1 = strlen(str1);
	size_t len2 = strlen(str2);
	char *str3 = malloc(len1 + len2 + 1);
	strcpy(str3, str1);
	strcpy(str3 + len1, str2);
	str3[len1 + len2] = '\0';
	return str3;
}

char *get_rev(const char *seq, int len)
{
	if (seq == NULL)
		return NULL;

	int i, k;
	char *ret = malloc(len + 1);
	for (i = 0, k = len - 1; i < len; ++i, --k)
		ret[i] = seq[k];
	ret[len] = '\0';
	return ret;
}

char *get_rev_complement(const char *seq, int len)
{
	if (seq == NULL)
		return NULL;

	char *ret = malloc(len + 1);
	int i, k;
	for (i = 0, k = len - 1; i < len; ++i, --k)
		ret[i] = rev_nt4_char[nt4_table[(int)seq[k]]];
	ret[len] = '\0';
	return ret;
}

double realtime()
{
	struct timeval tp;
#if defined(_MSC_VER)
	_gettimeofday(&tp,  NULL);
#else
	gettimeofday(&tp,  NULL);
#endif
	return tp.tv_sec + tp.tv_usec * 1e-6;
}

int64_t seq2num(const char *seq, int len)
{
	int64_t ret = 0;
	int i;
	for (i = 0; i < len; ++i)
		ret = ret * 5 + nt4_table[(int)seq[i]];
	return ret;
}

char *num2seq(int64_t num, int len)
{
	char *ret = malloc(len + 1);
	int i;
	for (i = 0; i < len; ++i) {
		ret[len - i - 1] = nt4_char[num % 5];
		num /= 5;
	}
	ret[len] = '\0';
	return ret;
}

void concat_str(char *s1, int l1, char *s2, int l2)
{
	memcpy(s1 + l1, s2, l2);
	s1[l1 + l2] = '\0';
}

int check_valid_nu(const char *seq, int len)
{
	int i;
	for (i = 0; i < len; ++i)
		if (nt4_table[(int)seq[i]] >= NNU)
			return 0;
	return 1;
}

int check_polyA(char *str, int len)
{
	int count[2] = {0, 0};
	int i;
	for (i = 0; i < len; ++i)
		if (str[i] == 'A')
			++count[0];
		if (str[i] == 'T')
			++count[1];
	if (count[0] > len/2 || count[1] > len /2)
		return 1;
	return 0;
}