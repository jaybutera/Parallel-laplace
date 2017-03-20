#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define SIZE 10

#define BLOCK_LOW(id,p) (id * SIZE) / p
#define BLOCK_OWNER(k,p,n) (((p)*(k+1)-1)/n)
#define MIN(x,y) (x > y) ? y : x

/*
typedef struct Matrix {
    float** mat;
    int rows;
    int cols;
}
*/

void file_to_mat(char* filename, float* A, int id, int p, int* coords, MPI_Comm grid_comm) {
    FILE* f = fopen(filename, "rb");

    int start_row = (id * SIZE) / p;
    int end_row   = ((id+1) * SIZE) / p - 1;

    int i_coord = coords[1];
    int j_coord = coords[0];

    int block_size = SIZE * (int)sqrt((double)p);
    int start_pos = (i_coord + j_coord) * block_size;

    printf("Proc %d reading lines [%d-%d]\n", id, start_row, end_row);

    fseek(f, start_row * SIZE * sizeof(float), SEEK_SET);  // Jump to the end of the file
    long offset = block_size * sizeof(float);

    // Read process block into memory
    int i;
    for (i = 0; i < block_size; i++) {
        fread((void*)(&A[block_size*i]), offset, 1, f);
    }

    fclose(f);
}

/*
void matMul(float* A, int a_row, int a_col,
            float* B, int b_row, int b_col,
            int size, float** C)
{
    if (size == 1)
        C[0][0] = A[0][0] * B[0][0];
    else {
        //int new_size = size / 2;
        int n_a_row = a_row / 2;
        int n_a_col = a_col / 2;
        int n_b_row = b_row / 2;
        int n_b_col = b_col / 2;

        matMul(A, n_a_row, n_a_col/2, B, n_b_row/2, n_b_col/2, C);
        matMul(A + (n_a_col * sizeof(float)), a_row/2, a_col/2,
               B + (n_b_rows * sizeof(float)), b_row/2, b_col/2, C);
    }
}

void matMul(Matrix* A, Matrix* B, int size, Matrix* C)
{
    if (size == 1)
        C->mat[0][0] = A->mat[A->rows][A->cols] * B->mat[B->rows][B->cols];

    else {
        int new_size = size / 2;

        A->mat = A->mat[0][0];
        B->mat = B->mat[0][0];
        matMul(A, B, new_size, C);
    }
}

void matSum (float** A, float** B, float**C, int c_row, int c_col) {
}
*/

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

// Auxiliary recursive matrix multiply
//void matMulAux (float** A,

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

    float** A;
    float* Astorage;

    // Allocate space
    Astorage = (float*) malloc(SIZE * SIZE * sizeof(float));
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (float**) malloc(SIZE * sizeof(float*));
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i = 0; i < SIZE; i++) {
        A[i] = &Astorage[i * SIZE];
    }

    int tmp_c[2] = {0,0};

    // Read in block from file
    file_to_mat("mp_mat", Astorage, 0, 1, tmp_c, grid_comm);

    int j;
    for (i = 0; i < SIZE; i++) {
        for (j = 0; j < SIZE; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }

    free(A);
    free(Astorage);

    return 0;
}
