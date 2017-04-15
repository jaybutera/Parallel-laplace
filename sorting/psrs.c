#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

#define dtype int
#define SIZE 10

dtype randDtype ();
int partition( int a[], int l, int r);
void quickSort( int a[], int l, int r);
void initArray (int arr[], int n);

int main (int argc, char** argv) {
    srand( time(NULL) );

    // Initialize array
    dtype local_arr[SIZE];
    initArray(local_arr, SIZE);

    // Sort locally
    quickSort(local_arr, 0, SIZE-1);

    int i;
    for (i = 0; i < SIZE; i++)
        printf("%d, ", local_arr[i]);
    printf("\n");

    return 0;
}

dtype randDtype () {
    return rand() % 8196;
}

void initArray (dtype arr[], int n) {
    int i;
    for (i = 0; i < n; i++)
        arr[i] = randDtype();
}

void quickSort( dtype a[], int l, int r) {
   int j;

   if( l < r )
   {
       // divide and conquer
        j = partition( a, l, r);
       quickSort( a, l, j-1);
       quickSort( a, j+1, r);
   }
}



int partition( dtype a[], int l, int r) {
   int pivot, i, j, t;
   pivot = a[l];
   i = l; j = r+1;

   while(1)
   {
       do ++i; while( a[i] <= pivot && i <= r );
       do --j; while( a[j] > pivot );

       if( i >= j ) break;
       t = a[i]; a[i] = a[j]; a[j] = t;
   }

   t = a[l]; a[l] = a[j]; a[j] = t;
   return j;
}

