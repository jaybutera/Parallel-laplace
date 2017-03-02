#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <omp.h>

#define ROWS 10
#define COLS 10

#define BLOCK_LOW(id,p) (id * ROWS) / p
#define BLOCK_OWNER(k,p,n) (k * p) / n
#define MIN(x,y) (x > y) ? y : x

double** A;
double* Astorage;

void file_to_mat(char* filename, int id, int p) {
    FILE* f = fopen(filename, "rb");

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;
    //int start_row = BLOCK_LOW(id, p);
    //int end_row   = BLOCK_LOW(id+1, p) - 1;

    printf("Proc %d reading lines [%d-%d]\n", id, start_row, end_row);

    fseek(f, start_row * COLS * sizeof(double), SEEK_SET);  // Jump to the end of the file
    //long filelen = ftell(f);     // Get the current byte offset in the file
    long offset = (end_row - start_row + 1) * COLS * sizeof(double);
    //printf("process %d offset: %d\n", id, offset);

    fread((void*)(Astorage), offset, 1, f);

    fclose(f);
}

void compute_shortest_paths (int id, int p, double** a) {
    int i, j, k;
    int offset;
    int root;
    double* tmp;

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;
    int local_rows = (end_row - start_row)+1;

    tmp = (double*) malloc (COLS * sizeof(double));
    for (k = 0; k < ROWS; k++) {
        root = BLOCK_OWNER(k,p, ROWS);

        /*
        if (k == ROWS-1)
            root = p-1;
        */

        if (root == id) {
            offset = k - BLOCK_LOW(id,p);
            for (j = 0; j < ROWS; j++)
                tmp[j] = a[offset][j];
        }

        MPI_Bcast(tmp, ROWS, MPI_DOUBLE, root, MPI_COMM_WORLD);
        for (i = 0; i < local_rows; i++)
            for (j = 0; j < ROWS; j++)
                a[i][j] = MIN(a[i][j], a[i][k]+tmp[j]);
    }

    free( tmp );
}

int main (int argc, char** argv) {
    int id; // Process id
    int num_procs;

    // Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    int start_row = (id * ROWS) / num_procs;
    int end_row   = ((id+1) * ROWS) / num_procs - 1;
    int local_rows = (end_row - start_row)+1;

    Astorage = (double*) malloc(local_rows * COLS * sizeof(double));
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (double**) malloc(local_rows * sizeof(double*));
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i = 0; i < local_rows; i++) {
        A[i] = &Astorage[i * COLS];
    }

    file_to_mat("mat_test", id, num_procs);


    compute_shortest_paths(id, num_procs, A);
    /*
    int j;
    for (i = 0; i < local_rows; i++) {
        printf("%d| ", id);

        for (j = 0; j < COLS; j++)
            printf("%f ", A[i][j]);
        printf("\n");
    }
    printf("\n");
    */
    MPI_Barrier(MPI_COMM_WORLD);


    if (id == num_procs-1) {
        // Tmp matrix
        //----------------------------
        double* Astorage_tmp = (double*) malloc(local_rows * COLS * sizeof(double));
        if (Astorage_tmp == NULL) {
            printf("Astorage mem could not allocate\n");
            exit(0);
        }

        double** A_tmp = (double**) malloc(local_rows * sizeof(double*));
        if (A_tmp == NULL) {
            printf("A mem could not allocate\n");
            exit(0);
        }

        int i;
        for (i = 0; i < local_rows; i++) {
            A_tmp[i] = &Astorage_tmp[i * COLS];
        }
        //----------------------------

        MPI_Status status;
        for (id = 0; id < num_procs-1; id++) {
            MPI_Recv(
                     Astorage_tmp,
                     (ROWS / num_procs) * COLS,
                     MPI_DOUBLE,
                     id,
                     0,
                     MPI_COMM_WORLD,
                     &status);

            // Print sub matrix
            int j;
            for (i = 0; i < (ROWS / num_procs ); i++) {
                for (j = 0; j < COLS; j++)
                    printf("%f ", A_tmp[i][j]);
                printf("\n");
            }
        }

        // Print sub matrix
        int j;
        for (i = 0; i < local_rows; i++) {
            for (j = 0; j < COLS; j++)
                printf("%f ", A[i][j]);
            printf("\n");
        }

        free(Astorage_tmp);
        free(A_tmp);
    }
    else {
        printf("Send %d elements from proc %d\n", local_rows * COLS, id);
        MPI_Send(Astorage, local_rows * COLS, MPI_DOUBLE, num_procs-1, 0, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    // Dealloc
    free(A);
    free(Astorage);

    MPI_Finalize();
    return 0;
}
