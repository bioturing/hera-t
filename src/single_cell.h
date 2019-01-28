#ifndef _SINGLE_CELL_H_
#define _SINGLE_CELL_H_

#ifdef _WIN32
void single_cell(int pos, int argc, TCHAR *argv[]);
#else
void single_cell(int pos, int argc, char *argv[]);
#endif

#endif
