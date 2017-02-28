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

    printf("Proc %d reading lines [%d-%d]\n", id, start_row, end_row);

    fseek(f, start_row * COLS * sizeof(double), SEEK_SET);  // Jump to the end of the file
    //long filelen = ftell(f);     // Get the current byte offset in the file
    long offset = (end_row - start_row + 1) * COLS * sizeof(double);
    printf("proces %d offset: %d\n", id, offset);
    //rewind(f);              // Jump back to beginning of file

    fread((void*)(Astorage + start_row * COLS), offset, 1, f);

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

    int start_row = (id * ROWS) / num_procs;
    int end_row   = ((id+1) * ROWS) / num_procs - 1;

    int i;
    for (i = 0; i < ROWS; i++) {
        A[i] = &Astorage[i * COLS];
    }

    file_to_mat("mat", id, num_procs);

    if (id) {
        int elem_count = (end_row-start_row)*COLS;
        MPI_Gather( Astorage + (start_row*COLS), elem_count, MPI_DOUBLE, NULL, elem_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    else { // Process 0
        double* Agather = (double*) malloc(ROWS * COLS * sizeof(double));
        if (Agather == NULL) {
            printf("Astorage mem could not allocate\n");
            exit(0);
        }

        double** Ag_ptr = (double**) malloc(ROWS * sizeof(double*));
        if (Ag_ptr == NULL) {
            printf("Astorage mem could not allocate\n");
            exit(0);
        }

        // (start data ptr, num of elements to send,
        int elem_count = (end_row-start_row)*COLS;
        //MPI_Gather( Astorage + (start_row*COLS), elem_count, MPI_DOUBLE, Agather+(start_row*COLS), elem_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Gather( Astorage + (start_row*COLS), elem_count, MPI_DOUBLE, Agather+(start_row*COLS), elem_count, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        int j;
        for (i = 0; i < ROWS; i++) {
            printf("%d| ", id);

            for (j = 0; j < COLS; j++)
                printf("%f ", A[i][j]);
            printf("\n");
        }
        printf("\n");

        free(Ag_ptr);
        free(Agather);
    }

    /*
    int j;
    for (i = start_row; i < end_row+1 ; i++) {
        printf("%d| ", id);

        for (j = 0; j < COLS; j++)
            printf("%f ", A[i][j]);
        printf("\n");
    }
    printf("\n");
    */

    // If processor 0
    /*
    if (!id) {
        int j;
        for (i = 0; i < ROWS; i++) {
            printf("%d| ", id);

            for (j = 0; j < COLS; j++)
                printf("%f ", A[i][j]);
            printf("\n");
        }
        printf("\n");
    }
    */

    // Dealloc
    free(A);
    free(Astorage);

    MPI_Finalize();
    return 0;
}
