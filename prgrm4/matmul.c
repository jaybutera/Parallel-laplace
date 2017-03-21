#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <unistd.h>
#include "smm.h"

#define SIZE 10

#define BLOCK_LEN(p) SIZE * (int)sqrt((double)p)

float** alloc_mat (int size, int p) {
    float** A;
    float* Astorage;
    int b_len = BLOCK_LEN(p);

    // Allocate space
    Astorage = (float*) malloc(b_len*b_len * sizeof(float));
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (float**) malloc(b_len * sizeof(float*));
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i = 0; i < b_len; i++) {
        A[i] = &Astorage[i * b_len];
    }

    // TODO: Need to free Astorage

    return A;
}

void file_to_mat(char* filename, float* A, int id, int p, int* coords, MPI_Comm grid_comm) {
    FILE* f = fopen(filename, "rb");

    int start_row = (id * SIZE) / p;
    int end_row   = ((id+1) * SIZE) / p - 1;

    int i_coord = coords[1];
    int j_coord = coords[0];

    int block_len = BLOCK_LEN(p);//SIZE * (int)sqrt((double)p);
    int start_pos = (i_coord + j_coord) * block_len;

    printf("Proc %d reading lines [%d-%d]\n", id, start_row, end_row);

    fseek(f, start_row * SIZE * sizeof(float), SEEK_SET);  // Jump to the end of the file
    long offset = block_len * sizeof(float);

    // Read process block into memory
    int i;
    for (i = 0; i < block_len; i++) {
        fread((void*)(&A[block_len*i]), offset, 1, f);
    }

    fclose(f);
}

void cannon_mult (float** A, float** B, float** C, int n, int m, int p, int id, MPI_Comm grid_comm) {
    int* coords;
    int i_coord,j_coord;

    // Get process coordinates
    MPI_Cart_coords(grid_comm, id, 2, coords);
    i_coord = coords[1]; // Row
    j_coord = coords[0]; // Col

    int hor_start_id;
    int vert_start_id;
    int my_id;

    // Get rank for process i steps to the left
    MPI_Cart_shift(grid_comm, 1, i_coord, &my_id, &hor_start_id);
    // Get rank for process j steps to up
    MPI_Cart_shift(grid_comm, 0, j_coord, &my_id, &vert_start_id);

    //printf("%d @ [%d,%d]\n", id, coords[1], coords[0]

    //MPI_sendrecv_replace(
}

int main(int argc, char** argv) {
    int id;
    int num_procs;

    // Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    //-------------------------
    // Create virtual topology
    //-------------------------

    int dim_sizes[2];
    int wrap[2];
    int reorder = 1; // Allow MPI to reorder process ranks for optimization
    int grid_id;

    MPI_Comm grid_comm;

    dim_sizes[0] = dim_sizes[1] = (int)sqrt(num_procs);
    wrap[0] = wrap[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes,
            wrap, reorder, &grid_comm);

    // Get process id in grid communicator
    MPI_Comm_rank(grid_comm, &grid_id);

    //-------------------------

    float** A = alloc_mat(SIZE, num_procs);
    float** B = alloc_mat(SIZE, num_procs);
    float** C = alloc_mat(SIZE, num_procs);

    int tmp_c[2] = {0,0};

    // Read in block from file
    file_to_mat("mp_mat", A[0], 0, 1, tmp_c, grid_comm);
    file_to_mat("mp_mat", B[0], 0, 1, tmp_c, grid_comm);

    rec_matmul(0,0,0,0,0,0,SIZE,SIZE,SIZE, A,B,C, SIZE);

    /*
    int i,j;
    for (i = 0; i < b_len; i++) {
        for (j = 0; j < b_len; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
        fflush(stdout);
    }
    */

    free(A);
    free(B);
    free(C);

    return 0;
}
