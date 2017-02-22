#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define ROWS 10
#define COLS 10

double** A;
double* Astorage;

void file_to_mat(char* filename, int id, int p) {
    FILE* f = fopen(filename, "rb");

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;

    fseek(f, start_row * COLS * sizeof(double), SEEK_SET);  // Jump to the end of the file
    //long filelen = ftell(f);     // Get the current byte offset in the file
    long filelen = (end_row - start_row) * COLS * sizeof(double);
    //rewind(f);              // Jump back to beginning of file

    fread((void*)Astorage, filelen, 1, f);

    fclose(f);
}

int main (int argc, char** argv) {
    int id; // Process id
    int num_procs;

    // Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    Astorage = (double*) malloc(ROWS * COLS * sizeof(double));
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (double**) malloc(ROWS * sizeof(double*));
    if (A == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;

    int i;
    for (i = 0; i < ROWS; i++) {
        A[i] = &Astorage[i * COLS];
    }

    file_to_mat("mat");

    int j;
    for (i = 0; i < start_row; i++) {
        for (j = 0; j < end_row; j++)
            printf("%f ", A[i][j]);
        printf("\n");
    }

    // Dealloc
    free(A);
    free(Astorage);


    fflush(stdout);
    MPI_Finalize();
    return 0;
}
