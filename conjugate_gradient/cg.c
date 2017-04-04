#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#define dtype double

#define SIZE 8

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
    //#pragma omp parallel for private(j)
    for (i = 0; i < n; i++) {
        b[i] = 0;
        for (j = 0; j < n; j++) {
            b[i] += A[i][j] * x[j];
        }
    }
}

dtype randDtype () {
    return (dtype)rand() / (dtype)RAND_MAX;
}

void initA (dtype** A, int n) {
    //dtype counter = 0;
    int i,j,k;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i == j) { // Make diagonal dominant
                A[i][j] = 1;
                for (k = 0; k < n; k++)
                    A[i][j] += abs(A[i][k]);
            }
            else // Non diag. is random
                A[i][j] = randDtype();
        }
    }
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

void initb (dtype** A, dtype* b, int n) {
    dtype* tmp_x = (dtype*) malloc( SIZE * sizeof(dtype) );

    printf("Original x\n-----------------\n");
    int i;
    for (i = 0; i < n; i++) {
        tmp_x[i] = randDtype();
        printf("%6.3f", tmp_x[i]);
    }
    printf("\n");

    // Compute b from random x
    mat_vec_mult(A, tmp_x, b, n);
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

    dtype rsold = inner_prod(residual, residual, n);

    int j;
    printf("\nresiduals\n------------\n");
    for (i = 0; i<n; i++) {
        mat_vec_mult(A, dir_vec, Ap, n);

        alpha = rsold / inner_prod(dir_vec, Ap, n);

        // Update x
        for (j = 0; j < n; j++)
            x[j] += alpha * dir_vec[j];

        // Update residual
        for (j = 0; j < n; j++)
            residual[j] -= alpha * Ap[j];

        rsnew = inner_prod(residual, residual, n);

        // Print residual error
        printf("%f\n", sqrt(rsnew));

        if (sqrt(rsnew) < .000001) {
            printf("\n[+] MINIMAL ERROR, EXIT ITERATION %d\n",i);
            break;
        }

        // Update direction vector
        for (j = 0; j < n; j++)
            dir_vec[j] = residual[j] + (rsnew/rsold) * dir_vec[j];

        rsold = rsnew;
    }

    free(Ap);
    free(dir_vec);
    free(residual);
}

int main(int argc, char** argv) {
    // Assume matrix A is [sizexsize]
    srand( time(NULL) );

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

    initA(A,SIZE);
    // -------------

    // Init vector b
    // -------------
    dtype* b = (dtype*) malloc(SIZE * sizeof(dtype*));
    if (b == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }
    initb(A,b,SIZE);
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
    // -------------

    // Time record
    struct timeval start_time, stop_time, elapsed_time;

    /*
    printf("A\n---------\n");
    printA(A, SIZE);
    printf("B\n---------\n");
    printb(b, SIZE);
    */

    //-----------
    // Start time
    //-----------
    gettimeofday(&start_time,NULL);

    // Compute conjugate gradient
    conjgrad(A, b, x, SIZE);
    //printf("\nx\n---------\n");
    //printb(x,SIZE);

    //---------
    // End time
    //---------
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    printf("\nFinal x\n--------------\n");
    printb(x,SIZE);

    // 2N^3/T
    //if (!id) {
        //float GFLOPS = (float)(2.f*ROWS*ROWS*ROWS) / (1000000000.f*(elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0));

        printf("\n--------------------------------\n");
        printf("elapsed time (s): %f\n", ((elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0)));
        //printf("GFLOPS: %f\n", GFLOPS);
    //}

    free(Astorage);
    free(A);
    free(b);
    free(x);

    return 0;
}
