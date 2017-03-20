#include <stdio.h>
#include <stdlib.h>

#define ROWS 10
#define COLS ROWS

#define BLOCK_LOW(id,p) (id * ROWS) / p
#define BLOCK_OWNER(k,p,n) (((p)*(k+1)-1)/n)
#define MIN(x,y) (x > y) ? y : x

typedef struct Matrix {
    float** mat;
    int rows;
    int cols;
}

void file_to_mat(char* filename, float** A, int id, int p) {
    FILE* f = fopen(filename, "rb");

    int start_row = (id * ROWS) / p;
    int end_row   = ((id+1) * ROWS) / p - 1;

    printf("Proc %d reading lines [%d-%d]\n", id, start_row, end_row);

    fseek(f, start_row * COLS * sizeof(float), SEEK_SET);  // Jump to the end of the file
    long offset = (end_row - start_row + 1) * COLS * sizeof(float);

    fread((void*)(A[0]), offset, 1, f);

    fclose(f);
}

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

void cannon_mult (float** A, float** B, float** C, int n, int m, int p) {
}

// Auxiliary recursive matrix multiply
//void matMulAux (float** A,

int main(int argc, char** argv) {
    float** A;
    float* Astorage;

    // Allocate space
    Astorage = (float*) malloc(ROWS * COLS * sizeof(float));
    if (Astorage == NULL) {
        printf("Astorage mem could not allocate\n");
        exit(0);
    }

    A = (float**) malloc(ROWS * sizeof(float*));
    if (A == NULL) {
        printf("A mem could not allocate\n");
        exit(0);
    }

    int i;
    for (i = 0; i < ROWS; i++) {
        A[i] = &Astorage[i * COLS];
    }

    // Read in block from file
    file_to_mat("mp_mat", A, 0, 1);

    int j;
    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < ROWS; j++) {
            printf("%f ", A[i][j]);
        }
        printf("\n");
    }

    //-------------------------
    // Create virtual topology
    //-------------------------

    MPI_Conn grid_comm;
    int dim_sizes[2];
    int wrap[2];
    int reorder = 1; // Allow MPI to reorder process ranks for optimization

    dim_sizes[0] = dim_sizes[1] = SIZE;
    wrap[0] = wrap[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes,
            wrap, reorder, &grid_comm);

    free(A);
    free(Astorage);

    return 0;
}
