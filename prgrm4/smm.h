#ifndef SMM_H_INCLUDED
#define SMM_H_INCLUDED

void matmul (int x, int y, int l, int m, int n);

void rec_matmul (int crow, int ccol,
                 int arow, int acol,
                 int brow, int bcol,
                 int l, int m, int n);
                 //float** a, float** b, float** c, int N);

#endif
