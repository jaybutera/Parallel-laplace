#define SMM_H_INCLUDED

void rec_matmul (int crow, int ccol,
                 int arow, int acol,
                 int brow, int bcol,
                 int l, int m, int n,
                 float** a, float** b, float** c, int N);
