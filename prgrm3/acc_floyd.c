#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <sys/time.h>

#define ROWS 2000
#define COLS ROWS

#define BLOCK_LOW(id,p) (id * ROWS) / p
#define BLOCK_OWNER(k,p,n) (((p)*(k+1)-1)/n)
#define MIN(x,y) (x > y) ? y : x

float** A;
float* Astorage;

void file_to_mat(char* filename, int id, int p) {
    FILE* f = fopen(filename, "rb");

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;
    //int start_row = BLOCK_LOW(id, p);
    //int end_row   = BLOCK_LOW(id+1, p) - 1;

    printf("Proc %d reading lines [%d-%d]\n", id, start_row, end_row);

    fseek(f, start_row * COLS * sizeof(float), SEEK_SET);  // Jump to the end of the file
    long offset = (end_row - start_row + 1) * COLS * sizeof(float);

    fread((void*)(Astorage), offset, 1, f);

    fclose(f);
}

void compute_shortest_paths (int id, int p, float** a) {
    int i, j, k;
    int offset;
    int root;
    float* tmp;

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;
    int local_rows = (end_row - start_row)+1;

    tmp = (float*) malloc (COLS * sizeof(float));
    #pragma acc data copy(a)
    for (k = 0; k < ROWS; k++) {
        //if (!id) {
            //clock_t start = clock(), diff;
        //}

        //root = BLOCK_OWNER(k,p,ROWS);

        /*
        if (k == ROWS-1) {
            root = p - 1;
        }
        */

        //if (root == id) {
            //offset = k - BLOCK_LOW(id,p);
            //printf("root %d with k: %d, low block: %d, offset: %d\n", root, k, BLOCK_LOW(id,p), offset);
            ///for (j = 0; j < ROWS; j++)
                //tmp[j] = a[k][j];
        //}

        ///printf("%d| k: %d root: %d local_rows: %d\n", id,k, root, local_rows);
        //fflush(stdout);

        #pragma acc kernels
        for (i = 0; i < local_rows; i++)
            for (j = 0; j < ROWS; j++)
                a[i][j] = MIN(a[i][j], a[i][k]+a[k][j]);

        /*
        if (!id) {
            diff = clock() - start;
            printf("k: %d time: %d\n", k, diff / 1000000);
        }
        */
    }

    free( tmp );
}

int main (int argc, char** argv) {
    int id=0; // Process id
    int num_procs=1;

    // Time record
    struct timeval start_time, stop_time, elapsed_time;

    int start_row = (id * ROWS) / num_procs;
    int end_row   = ((id+1) * ROWS) / num_procs - 1;
    int local_rows = (end_row - start_row)+1;

    Astorage = (float*) malloc(local_rows * COLS * sizeof(float)); if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (float**) malloc(local_rows * sizeof(float*));
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i = 0; i < local_rows; i++) {
        A[i] = &Astorage[i * COLS];
    }

    file_to_mat("mp_mat", id, num_procs);


    //-----------
    // Start time
    //-----------
    gettimeofday(&start_time,NULL);

    compute_shortest_paths(0, 1, A);


    //---------
    // End time
    //---------
    gettimeofday(&stop_time,NULL);

    timersub(&stop_time, &start_time, &elapsed_time); // Unix time subtract routine

    // 2N^3/T
        float GFLOPS = (float)(2.f*ROWS*ROWS*ROWS) / (1000000000.f*(elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0));

        printf("elapsed time (s): %f\n", ((elapsed_time.tv_sec+elapsed_time.tv_usec/1000000.0)));
        printf("GFLOPS: %f\n", GFLOPS);

    // Print matrix
    /*
    for (int i = 0; i < ROWS; i++) {
        for (int j = 0; j < COLS; j++)
            printf("%10.1f ", A[i][j]);
        printf("\n");
    }
    */

    // Dealloc
    free(A);
    free(Astorage);

    return 0;
}
