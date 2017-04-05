#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

#define dtype float

#define SIZE 10

dtype inner_prod (dtype* a, dtype* b, int n) {
    // 2*N FLOP
    dtype prod = 0;

    int i;
    #pragma omp parallel for
    //#pragma acc kernels
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

void mat_vec_mult (dtype** A, dtype* x, dtype* b, int n, int m) {
    // 2*N^2 FLOP
    int i,j;
    #pragma omp parallel for private(j)
    //#pragma acc kernels
    for (i = 0; i < n; i++) {
        b[i] = 0;
        //#pragma acc kernels
        for (j = 0; j < m; j++) {
            b[i] += A[i][j] * x[j];
        }
    }
}

dtype randDtype () {
    return (dtype)rand() / (dtype)RAND_MAX;
}

void initA (dtype** A, int band, int n) {
    //dtype counter = 0;
    int i,j,k;
    int m = 2*band+1; // Num columns

    #pragma omp parallel for private(j)
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++) {
            A[i][j] = randDtype();
        }
    }

    // Insert 0s at edge case
    for (i = 0; i < band+1; i++)
        for (j = band+(i+1); j < m; j++) {
            printf("A[%d][%d] = 0\n",i,j);
            A[i][j] = 0;
        }
    for (i = n-1; i > n-(band+1); i--)
        for (j = m-(band+(n-i)); j >= 0; j--) {
            printf("A[%d][%d] = 0\n",i,j);
            A[i][j] = 0;
         }


    // Make diagonal dominant
    #pragma omp parallel for private(k)
    for (i = 0; i < n; i++) {
        A[i][band] = 1;
        for (k = 0; k < m; k++)
            if (k != band)
                A[i][band] += A[i][k];
    }
}

void printA (dtype** A, int n, int m) {
    int i,j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < m; j++)
            printf("%8.3f", A[i][j]);
        printf("\n");
    }
    printf("\n");
}

void initb (dtype** A, dtype* b, int n, int m) {
    dtype* tmp_x = (dtype*) malloc( m * sizeof(dtype) );

    //printf("Original x: \n");
    int i;
    for (i = 0; i < m; i++) {
        tmp_x[i] = randDtype();
        //printf("%6.3f", tmp_x[i]);
    }
    //printf("\n");
    //printf("Original x[0]: %6.3f\n", tmp_x[0]);

    // Compute b from random x
    mat_vec_mult(A, tmp_x, b, n,m);
}

int conjgrad (dtype** A, dtype* b, dtype* x, int n, int m) {
    int i;
    // Allocate space
    dtype* residual = (dtype*) malloc( m * sizeof(dtype) );
    dtype* dir_vec  = (dtype*) malloc( m * sizeof(dtype) );
    dtype* Ap       = (dtype*) malloc( m * sizeof(dtype) );
    dtype alpha;
    dtype rsnew;


    // Compute initial residual
    mat_vec_mult(A, x, residual, n, m);
    for (i = 0; i < m; i++) {
        residual[i] = b[i] - residual[i];

        // Residual is the initial search direction
        dir_vec[i] = residual[i];
    }

    dtype rsold = inner_prod(residual, residual, m);

    int j;
    //printf("\nresiduals\n------------\n");
    for (i = 0; i<n; i++) {
        mat_vec_mult(A, dir_vec, Ap, n,m);

        alpha = rsold / inner_prod(dir_vec, Ap, m);

        // Update x
        for (j = 0; j < m; j++)
            x[j] += alpha * dir_vec[j];

        // Update residual
        for (j = 0; j < m; j++)
            residual[j] -= alpha * Ap[j];

        rsnew = inner_prod(residual, residual, m);

        // Print residual error
        //printf("%f\n", sqrt(rsnew));

        if (sqrt(rsnew) < .000001) {
            printf("\n[+] MINIMAL ERROR, EXIT ITERATION %d\n",i);
            break;
        }

        // Update direction vector
        for (j = 0; j < m; j++)
            dir_vec[j] = residual[j] + (rsnew/rsold) * dir_vec[j];

        rsold = rsnew;
    }

    printf("Final residual: %f\n", sqrt(rsnew));

    free(Ap);
    free(dir_vec);
    free(residual);

    return i;
}

int main(int argc, char** argv) {
    // Assume matrix A is [sizexsize]
    srand( time(NULL) );
    int bandsize = 2; // Init matrix w/ semibandwith of bandsize
    int cols = 2*bandsize+1;

    // Init matrix A
    // -------------
    dtype* Astorage = malloc ( SIZE * cols * sizeof(dtype) );
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
        A[i] = &Astorage[i * cols];

    initA(A,bandsize,SIZE);
    printA(A,SIZE,cols);
    // -------------

    // Init vector b
    // -------------
    dtype* b = (dtype*) malloc(SIZE * sizeof(dtype*));
    if (b == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }
    initb(A,b,SIZE,cols);
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
    int iters = conjgrad(A, b, x, SIZE,cols);
    //printf("\nx\n---------\n");
    //printb(x,SIZE);

    //---------
    // End time
    //---------
    gettimeofday(&stop_time,NULL);
    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    //printf("\nFinal x[0]: %6.3f\n", x[0]);
    //printf("\nFinal x\n--------------\n");
    //printb(x,SIZE);

    // 2N^3/T
    //if (!id) {
        float GFLOPS = (float)((2.f*SIZE*SIZE + 3.f*SIZE) +
                       (float)(iters * 2.f*(SIZE*SIZE + 5.f*SIZE))) /
                       (float)(1000000000.f*(elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0));

        printf("\n--------------------------------\n");
        printf("elapsed time (s): %f\n", ((elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0)));
        printf("GFLOPS: %f\n", GFLOPS);
    //}

    free(Astorage);
    free(A);
    free(b);
    free(x);

    return 0;
}
