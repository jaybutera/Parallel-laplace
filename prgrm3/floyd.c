#include <stdio.h>
#include <stdlib.h>

#define ROWS 10
#define COLS 10

double** A;
double* Astorage;

void file_to_mat(char* filename) {
    FILE* f = fopen(filename, "rb");

    fseek(f, 0, SEEK_END);  // Jump to the end of the file
    long filelen = ftell(f);     // Get the current byte offset in the file
    rewind(f);              // Jump back to beginning of file

    fread((void*)Astorage, filelen, 1, f);
}

int main (int argc, char** argv) {
    //double** A;
    //double* Astorage;

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

    int i;
    for (i = 0; i < ROWS; i++) {
        A[i] = &Astorage[i * COLS];
    }

    file_to_mat("mat");

    int j;
    for (i = 0; i < ROWS; i++) {
        for (j = 0; j < COLS; j++)
            printf("%f ", A[i][j]);
        printf("\n");
    }

    // Dealloc
    free(A);
    free(Astorage);

    return 0;
}
