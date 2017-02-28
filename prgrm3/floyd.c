#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define ROWS 10
#define COLS 10

#define BLOCK_LOW(id,p) (id * ROWS) / p
//#define BLOCK_OWNER(k,p) 

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
    printf("process %d offset: %d\n", id, offset);
    //rewind(f);              // Jump back to beginning of file

    //fread((void*)(Astorage + start_row * COLS), offset, 1, f);
    fread((void*)A[0], offset, 1, f);

    fclose(f);
}

/*
void compute_shortest_paths (int id, int p, double** a, int n) {
    int i, j, k;
    int offset;
    int root;
    int* tmp;

    tmp = (double*) malloc (n * sizeof(double));
    for (k = 0; k < n; k++) {
        root = BLOCK_OWNER(k,p,n);

        if (root == id) {
            offset = k - BLOCK_LOW(id,p);
            for (j = 0; j < n; j++)
                tmp[j] = a[offset][j];
        }

        MPI_Bcast(tmp, n, MPI_DOUBLE, root, MPI_COMM_WORLD);
        for (i = 0; i < BLOCK_SIZE(id,p,n); i++)
            for (j = 0; j < n; j++)
                a[i][j] = min(a[i][j], a[i][k]+tmp[j]);
    }

    free( tmp );
}
*/

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

    Astorage = (double*) malloc((end_row - start_row) * COLS * sizeof(double));
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (double**) malloc(ROWS * sizeof(double*));
    if (A == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i = 0; i < local_rows; i++) {
        A[i] = &Astorage[i * COLS];
    }

    file_to_mat("mat", id, num_procs);

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


    /*
    if (!id) {
        // Print sub matrix
        int j;
        for (i = 0; i < ROWS; i++) {
            for (j = 0; j < COLS; j++)
                printf("%f ", A[i][j]);
            printf("\n");
        }

        MPI_Status status;
        for (i = 1; i < num_procs; i++) {
            MPI_Recv(//A[ BLOCK_LOW(i, num_procs) ],
                     &Astorage,
                     (ROWS / num_procs + 1) * COLS,
                     MPI_DOUBLE,
                     i,
                     0,
                     MPI_COMM_WORLD,
                     &status);
        }
    }
    else {
        printf("Send from proc %d\n", id);
        MPI_Send(&Astorage, local_rows * COLS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }
    */

    // Dealloc
    free(A);
    free(Astorage);

    MPI_Finalize();
    return 0;
}
