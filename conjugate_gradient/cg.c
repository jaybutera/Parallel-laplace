#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define dtype double

#define SIZE 2

dtype inner_prod (dtype* a, dtype* b, int n) {
    dtype prod = 0;

    int i;
    for (i = 0; i < n; i++) {
        prod += a[i] * b[i];
    }

    return prod;
}

void printb (dtype* b, int n) {
    int i;
    for (i = 0; i < n; i++)
        printf("%8.3f", b[i]);
    printf("\n");
}

void mat_vec_mult (dtype** A, dtype* x, dtype* b, int n) {
    int i,j;
    //printf("b\n---------\n");
    for (i = 0; i < n; i++) {
        b[i] = 0;
        //printb(b,n);
        //printf("\n---------\n");
        for (j = 0; j < n; j++) {
            //printf("b[%d] (%f) += A[%d][%d] (%f) * x[%d] (%f)\n",i,b[i],i,j,A[i][j],j,x[j]);
            b[i] += A[i][j] * x[j];
        }
    }
}

void initA (dtype** A, int n) {
    /*
    dtype counter = 0;
    int i,j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            A[i][j] = counter++;
    */
    A[0][0] = 4;
    A[0][1] = 1;
    A[1][0] = 1;
    A[1][1] = 3;
}

void printA (dtype** A, int n) {
    int i,j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%8.3f", A[i][j]);
        printf("\n");
    }
    printf("\n");
}

void initb (dtype* b, int n) {
    int i,j;
    for (i = 0; i < n; i++)
        b[i] = i+1;
}

void conjgrad (dtype** A, dtype* b, dtype* x, int n) {
    int i;
    // Allocate space
    dtype* residual = (dtype*) malloc( SIZE * sizeof(dtype) );
    dtype* dir_vec  = (dtype*) malloc( SIZE * sizeof(dtype) );
    dtype* Ap       = (dtype*) malloc( SIZE * sizeof(dtype) );
    /*
    for (i = 0; i < n; i++) {
        residual[i] = 0;
        //dir_vec[i] = 0;
        Ap[i] = 0;
    }
    */
    dtype alpha;
    dtype rsnew;


    // Compute initial residual
    mat_vec_mult(A, x, residual, n);
    for (i = 0; i < n; i++) {
        residual[i] = b[i] - residual[i];

        // Residual is the initial search direction
        dir_vec[i] = residual[i];
    }
    //printf("residual\n----------\n");
    //printb(dir_vec,n);

    dtype rsold = inner_prod(residual, residual, n);

    int j;
    //printf("residuals\n------------\n");
    for (i = 0; i < n; i++) {
        mat_vec_mult(A, dir_vec, Ap, n);
        //printf("Ap\n-----------\n");
        //printb(Ap, n);

        //printf("---\np * Ap: %f\n", inner_prod(dir_vec, Ap, n));

        alpha = rsold / inner_prod(dir_vec, Ap, n);
        //printf("\nalpha: %f\n", alpha);

        // Update x
        for (j = 0; j < n; j++)
            x[j] += alpha * dir_vec[j];
        //printf("\nnew x:\n");
        //printb(x,n);

        // Update residual
        for (j = 0; j < n; j++)
            residual[j] -= alpha * Ap[j];
        //printf("residual\n----------\n");
        //printb(residual,n);

        rsnew = inner_prod(residual, residual, n);
        printf("%f\n", sqrt(rsnew));

        if (sqrt(rsnew) < .000001) {
            printf("MINIMAL ERROR, EXIT\n");
            break;
        }

        // Update direction vector
        for (j = 0; j < n; j++)
            dir_vec[j] = residual[j] + (rsnew/rsold) * dir_vec[j];
        //printf("\ndir vec\n----------\n");
        //printb(dir_vec,n);

        rsold = rsnew;
    }

    free(Ap);
    free(dir_vec);
    free(residual);
}

int main(int argc, char** argv) {
    // Assume matrix A is [sizexsize]

    // Init matrix A
    // -------------
    dtype* Astorage = malloc ( SIZE * SIZE * sizeof(dtype) );
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    dtype** A = (dtype**) malloc(SIZE * sizeof(dtype*));
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i=0; i < SIZE; i++)
        A[i] = &Astorage[i * SIZE];

    printf("Init A\n");
    initA(A,SIZE);
    // -------------

    // Init vector b
    // -------------
    dtype* b = (dtype*) malloc(SIZE * sizeof(dtype*));
    if (b == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }
    printf("Init b\n");
    initb(b,SIZE);
    // -------------


    // Init vector x
    // -------------
    dtype* x = (dtype*) malloc (SIZE * sizeof(dtype));
    if (x == NULL) {
        printf("x mem could not allocate\n");
        exit(0);
    }
    for (i = 0; i < SIZE; i++)
        x[i] = 0;
    x[0] = 2;
    x[1] = 1;
    // -------------

    printf("A\n---------\n");
    printA(A, SIZE);
    printf("B\n---------\n");
    printb(b, SIZE);

    // Compute conjugate gradient
    conjgrad(A, b, x, SIZE);
    printf("\nx\n---------\n");
    printb(x,SIZE);

    free(Astorage);
    free(A);
    free(b);
    free(x);

    return 0;
}
