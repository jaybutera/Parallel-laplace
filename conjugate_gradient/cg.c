#include <stdio.h>

#define SIZE 2

float inner_prod (float* a, float* b, int n) {
    float prod;

    int i;
    for (i = 0; i < n; i++) {
        prod += a[i] * b[i];
    }

    return prod;
}

void mat_vec_mult (float** A, float* x, float* b, int n) {
    int i,j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            b[i] += A[i][j] * x[i];
}

void conjgrad (float** A, float* b, float* x) {
    // Allocate space
    float* residual = (float*) malloc( SIZE * sizeof(float) );
    float* dir_vec  = (float*) malloc( SIZE * sizeof(float) );
    float* Ap       = (float*) malloc( SIZE * sizeof(float) );
    float alpha;
    float rsnew;


    float rsold = inner_prod(residual, residual, n);
    int i;
    for (i = 0; i < n; i++) {
        mat_vec_mult(A, dir_vec, Ap);
        alpha = rsold / inner_prod(dir_vec, Ap);

        // Update x
        int j;
        for (j = 0; j < n; j++)
            x[j] += alpha * dir_vec[j];

        // Update residual
        for (j = 0; j < n; j++)
            residual[j] -= alpha * dir_vec[j];

        rsnew = inner_prod(residual, residual);

        if (math.sqrt(rsnew) < .001) break;

        // Update direction vector
        for (j = 0; j < n; j++)
            dir_vec[j] = residual[j] + (rsnew/rsold) * dir_vec[j];

        rsold = rsnew;
    }

    free(Ap);
    free(dir_vec);
    free(residual);
}

int main(int argc, char* argv) {
    // Assume matrix A is [sizexsize]

    // Init matrix A
    // -------------
    float* Astorage = malloc ( SIZE * SIZE * sizeof(float) );
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    float** A = (float**) malloc(SIZE * sizeof(float*));
    float** A = malloc (SIZE * sizeof(float*) );
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i=0; i < SIZE; i++)
        A[i] = Astorage[i * SIZE];
    // -------------

    // Init vector b
    // -------------
    float* b = (float**) malloc(SIZE * sizeof(float*));
    if (b == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }
    // -------------


    // Allocate x
    float* x = (float*) malloc (SIZE * sizeof(float));

    // Compute conjugate gradient
    conjgrad(A, b, x);
    return 0;
}
