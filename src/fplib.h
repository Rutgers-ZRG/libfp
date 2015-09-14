
#include <cell.h>

#define NC  3 
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

void get_fingerprint_periodic(int nat, int ntyp, int types[], double lat[3][3],
        double rxyz[][3], char *symb[], int natx, char *orb, double cutoff,
        double *sfp[], double *lfp[]);

extern void dsyev( char* jobz, char* uplo, int* n, double* a, int* lda,
                double* w, double* work, int* lwork, int* info );

int get_fp(Cell * cell, int ixyz, int nx, int lseg, int l, double cutoff);

