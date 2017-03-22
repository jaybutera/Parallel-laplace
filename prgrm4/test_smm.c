#include "smm.h"
#include <stdio.h>
#include "common.h"

float** alloc_mat (int size) {
    float** A;
    float* Astorage;

    // Allocate space
    Astorage = (float*) malloc((long)size*(long)size * (long)sizeof(float));
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (float**) malloc(size * sizeof(float*));
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i = 0; i < size; i++) {
        A[i] = &Astorage[i * size];
    }

    // TODO: Need to free Astorage

    return A;
}

void gen_submats (float** A, float** B, int size, int start_coords[2]) {
    int i,j;
    for (i = start_coords[0]; i < size; i++)
        for (j = start_coords[1]; j < size; j++) {
            A[i][j] = j-i;
            B[i][j] = SIZE - j + i;
        }
}

void printmat (float** mat, int n) {
    fflush(stdout);
    int i,j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%6.3f",mat[i][j]);
        putchar('\n');
    }
}

int main(int argc, char** argv) {
    float** A = alloc_mat(SIZE);
    float** B = alloc_mat(SIZE);
    float** C = alloc_mat(SIZE);

    int coords[2] = {0,0};
    gen_submats(A,B,SIZE, coords);
    printf("A\n-------------\n");
    printmat(A, SIZE);
    printf("B\n-------------\n");
    printmat(B, SIZE);

    rec_matmul(0,0,0,0,0,0,SIZE,SIZE,SIZE,A,B,C,SIZE);

    return 0;
}
