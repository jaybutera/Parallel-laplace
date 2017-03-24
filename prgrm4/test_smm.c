#include "smm.h"
#include <stdio.h>
#include "common.h"

void gen_submats (int size, int start_coords[2]) {
    int i,j;
    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++) {
            if (i == j - 1 || j == i - 1) {
                A[i][j] = 1.;
                B[i][j] = 1.;
            }
            else {
                A[i][j] = 0.;
                B[i][j] = 0.;
            }
            /*
            A[i][j] = j-i;
            B[i][j] = SIZE-j+i;
            */
            C[i][j] = 0.;
        }
}

void printmat (DTYPE mat[SIZE][SIZE], int n) {
    fflush(stdout);
    int i,j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%7.0f",mat[i][j]);
        putchar('\n');
    }
}

int verify_result() {
    int i,j;
    for (i = 0; i < SIZE; i++)
        for (j = 0; j < SIZE; j++)
            if (C[i][j] != SIZE)
                return 0;
    return 1;
}

int main(int argc, char** argv) {
    /*
    A = alloc_mat(SIZE);
    B = alloc_mat(SIZE);
    C = alloc_mat(SIZE);
    */

    int coords[2] = {0,0};
    gen_submats(SIZE, coords);
    printf("A\n-------------\n");
    printmat(A, SIZE);
    printf("B\n-------------\n");
    printmat(B, SIZE);

    rec_matmul(0,0,0,0,0,0,SIZE,SIZE,SIZE);
    //matmul(0,0,SIZE,SIZE,SIZE);

    printf("C\n-------------\n");
    printmat(C, SIZE);
    /*
    if( verify_result() )
        printf("successful result\n");
    else {
        printf("failed test\n");
        //printf("C[10][10] = %f\n",C[30][30]);
    }
    */

    //printmat(C, SIZE);

    return 0;
}
