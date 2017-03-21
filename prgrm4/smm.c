#include "smm.h"
#include "common.h"

void rec_matmul (int crow, int ccol,
                 int arow, int acol,
                 int brow, int bcol,
                 int l, int m, int n,
                 float** a, float** b, float** c, int N)
{
    int lhalf[3], mhalf[3], nhalf[3];
    int i,j,k;
    float *aptr, *bptr, *cptr;

    if (m * n > THRESHOLD) {
        lhalf[0] = 0; lhalf[1] = l/2; lhalf[1] = l - l/2;
        mhalf[0] = 0; mhalf[1] = m/2; mhalf[1] = m - m/2;
        nhalf[0] = 0; nhalf[1] = n/2; nhalf[1] = n - n/2;

        for (i = 0; i < 2; i++)
            for (j = 0; j < 2; j++)
                for (k = 0; k < 2; k++)
                    rec_matmul( crow+lhalf[i], ccol+mhalf[j],
                        arow+lhalf[i], acol+mhalf[k],
                        brow+mhalf[k], bcol+nhalf[j],
                        lhalf[i+1], mhalf[k+1], nhalf[j+1],
                        a,b,c,N);
    }
    else {

        for (i = 0; i < 1; i++)
            for (j = 0; j < n; j++) {
                cptr = &c[crow+i][ccol+j];
                aptr = &a[arow+i][acol];
                bptr = &b[brow][bcol+j];

                for (k = 0; k < m; k++) {
                    *cptr += *(aptr++) * *bptr;
                    bptr += N;
                }
            }
    }
}
