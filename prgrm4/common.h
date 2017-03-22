#ifndef COMMON_INCLUDE_FILE
#define COMMON_INCLUDE_FILE
// ------------------------

#define SIZE 10
#define THRESHOLD 2
#define PTR_SIZE sizeof(float*)

#define MALLOC_ERROR -2
#define TYPE_ERROR -3

#define BLOCK_LOW(id,p,n) ((id)*n/p)
#define BLOCK_HIGH(id,p,n) (BLOCK_LOW(id+1,p,n)-1)
#define BLOCK_SIZE(id,p,n) BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1
#define BLOCK_OWNER(id,p,n) (p*j+1)-1/n

#define BLOCK_LEN(p) SIZE * (int)sqrt((float)p)

// ------------------------
#endif
