#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#include <unistd.h>

#include "smm.h"
#include "common.h"

float* my_malloc (int id, int bytes) {
    float* buffer;

    if ((buffer = (float*)malloc((size_t) bytes)) == NULL) {
        printf("Error: Malloc failed for process %d\n", id);
        fflush(stdout);
        MPI_Abort(MPI_COMM_WORLD, MALLOC_ERROR);
    }

    return buffer;
}

/*
void read_checkerboard_matrix (
        char* s,
        void*** subs,
        void** storage,
        MPI_Datatype dtype,
        MPI_Comm grid_comm)
{
    float* buffer;
    int coords[2];
    int datum_size;
    int dest_id;
    int grid_coord[2];
    int grid_id;
    int grid_period[2];
    int grid_size[2];
    int i,j,k;
    FILE *infileptr;
    float* laddr;
    int local_cols;
    int local_rows;
    float** lptr;
    int p;
    void* raddr;

    float* rptr;
    MPI_Status status;

    MPI_Comm_rank(grid_comm, &grid_id);
    MPI_Comm_size(grid_comm, &p);
    datum_size = get_size (dtype);

    if (grid_id == 0) {
        infileptr = fopen(s, "r");
        //if (infileptr == NULL) *m = 0;
        //else {
        //    fread
        //}
    }

    MPI_Cart_get(grid_comm, 2, grid_size, grid_period, grid_coord);

    local_rows = BLOCK_SIZE(grid_coord[0], grid_size[0], SIZE);
    local_cols = BLOCK_SIZE(grid_coord[1], grid_size[1], SIZE);

    if (!grid_id) {
        printf("read submatrix size (%dx%d)\n",local_rows, local_cols);
    }

    *storage = my_malloc(grid_id, local_rows * local_cols * datum_size);
    *subs = (float**) my_malloc(grid_id, local_rows * PTR_SIZE);
    lptr = (float**) *subs;
    rptr = (float*) *storage;

    for (i = 0; i < local_rows; i++) {
        *(lptr++) = (float*) rptr;
        rptr += local_rows * datum_size;
    }

    if (grid_id == 0)
        buffer = my_malloc(grid_id, SIZE * datum_size);

    for (i = 0; i < grid_size[0]; i++) {
        coords[0] = i;

        for (j = 0; j < BLOCK_SIZE(i, grid_size[0], SIZE); j++) {
            if (grid_id == 0)
                fread(buffer, datum_size, SIZE, infileptr);

            for (k = 0; k < grid_size[1]; k++) {
                coords[1] = k;

                raddr = buffer + BLOCK_LOW(k, grid_size[1], SIZE) * datum_size;

                MPI_Cart_rank(grid_comm, coords, &dest_id);

                if (grid_id == 0) {
                    if (dest_id == 0) {
                        laddr = (*subs)[j];
                        memcpy (laddr, raddr, local_cols * datum_size);
                    }
                    else {
                        MPI_Send(raddr, BLOCK_SIZE(k, grid_size[1], SIZE), dtype, dest_id, 0, grid_comm);
                    }
                }
                else if (grid_id == dest_id) {
                    MPI_Recv((*subs)[j], local_cols, dtype, 0, 0, grid_comm, &status);
                }
            }
        }
    }

    if (grid_id == 0) free (buffer);

}
*/

int get_size (MPI_Datatype t) {
    if (t == MPI_FLOAT) return sizeof(float);
    printf("Error: Unrecognized argument to 'get_size'\n'");
    fflush(stdout);
    MPI_Abort(MPI_COMM_WORLD, TYPE_ERROR);
}

void print_subvector(float* a, int n) {
    int i;

    for (i = 0; i < n; i++) {
        printf("%6.3f", ((float*)a)[i]);
    }
}

void print_checkerboard_matrix (
        float a[SIZE][SIZE],
        MPI_Datatype dtype,
        //int m,
        //int n,
        MPI_Comm grid_comm)
{
    void* buffer;
    int coords[2];
    int datum_size;
    int els;
    int grid_coords[2];
    int grid_id;
    int grid_period[2];
    int grid_size[2];
    int i,j,k;
    void* laddr;
    int local_cols;
    int p;
    int src;
    int n=SIZE,m=SIZE;
    MPI_Status status;

    MPI_Comm_rank(grid_comm, &grid_id);
    MPI_Comm_size(grid_comm, &p);
    datum_size = get_size(dtype);

    MPI_Cart_get (grid_comm, 2, grid_size, grid_period, grid_coords);
    local_cols = BLOCK_SIZE(grid_coords[1], grid_size[1], n);

    if (!grid_id)
        buffer = my_malloc (grid_id, n*datum_size);

    // For each row of the process grid
    for (i = 0; i < grid_size[0]; i++) {
        coords[0] = 1;

        // For each matrix row in a process row
        for (j = 0; j < BLOCK_SIZE(i, grid_size[0], m); j++) {

            // Collect the matrix row on process 0 and print it
            if (!grid_id) {
                for (k = 0; k < grid_size[1]; k++) {
                    coords[1] = k;
                    MPI_Cart_rank(grid_comm, coords, &src);
                    els = BLOCK_SIZE(k, grid_size[1], n);
                    laddr = buffer + BLOCK_LOW(k, grid_size[1], n) * datum_size;

                    if (src == 0) {
                        memcpy( laddr, a[j], els * datum_size);
                    }
                    else {
                        MPI_Recv(laddr, els, dtype, src, 0, grid_comm, &status);
                    }
                }
                print_subvector(buffer, n);
                putchar('\n');
            }
            else if (grid_coords[0] == i) {
                MPI_Send(a[j], local_cols, dtype, 0,0, grid_comm);
            }
        }
    }

    if (!grid_id) {
        free(buffer);
        putchar('\n');
    }
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

void gen_submats (int start_coords[2]) {
    int i,j;
    for (i = start_coords[0]; i < SIZE; i++)
        for (j = start_coords[1]; j < SIZE; j++) {
            A[i][j] = j-i;
            B[i][j] = SIZE - j + i;
        }
}

int main(int argc, char** argv) {
    int global_id;
    int grid_id;
    int num_procs;

    // Init MPI
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &global_id);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    //-------------------------
    // Create virtual topology
    //-------------------------

    int dim_sizes[2];
    int wrap[2];
    int reorder = 1; // Allow MPI to reorder process ranks for optimization
    int coords[2];

    MPI_Comm grid_comm;

    dim_sizes[0] = dim_sizes[1] = (int)sqrt(num_procs);
    wrap[0] = wrap[1] = 1;

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes,
            wrap, reorder, &grid_comm);

    // Get process id in grid communicator
    MPI_Comm_rank(grid_comm, &grid_id);
    // Get coordinates
    MPI_Cart_coords(grid_comm, grid_id, 2, coords);

    //-------------------------

    // Generate matrix block for proccess
    gen_submats(coords);

    print_checkerboard_matrix(A, MPI_FLOAT, grid_comm);

    return 0;
}
