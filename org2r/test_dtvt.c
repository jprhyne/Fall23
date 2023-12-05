// Testing the functionality of the tvt subroutine
// We need to ensure 
// 1) That the output of TVT is such that T = T*V**T
#ifdef USE_AOCL
    #define SOURCE "AOCL"
#else
    #define SOURCE "REF_LAPACK"
#endif

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>

int main(int argc, char *argv[]) {
    // integer variables
    int info, ldv, ldt, m, n, k, lwork, nb, i, j;
    int workQuery = -1;
    // double variables
    double *A, *As, *T, *work=NULL;
    double normA;
    double elapsed_refL, norm_orth_1, norm_repres_1;
    double zero = 0;
    double one = 1;
    double negOne = -1;
    // struct for help with timing
    struct timeval tp;
    // character variables
    char aChar = 'A';
    char fChar = 'F';
    char lChar = 'L';
    char nChar = 'N';
    char rChar = 'R';
    char tChar = 'T';
    char uChar = 'U';
    // Dummy value that is used when calling fortran
    // functions that take characters from C
    size_t dummy = 0;

    m = 30;
    n = 20;
    ldv = -1;
    ldt = -1;

    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-ldt") == 0) {
            ldt  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-ldv") == 0) {
            ldv  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-m") == 0) {
            m  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-n") == 0) {
            n  = atoi( *(argv + i + 1) );
            i++;
        }
    }

    if( ldv < 0 ) ldv = m;
    if( ldt < 0 ) ldt = m;

    char *source = SOURCE;
    printf("dtvt %s: m = %4d, n = %4d, ldv = %4d, ldt = %4d\n",source, m, n, ldv, ldt);

    // Generate two matrices. one is upper triangular and one is a tall and skinny rectangle
    // T = T*V**T
    // T is n by n
    // V is m by n
    V = (double *) malloc(ldv * n * sizeof(double));
    T = (double *) malloc(ldt * n * sizeof(double));
    for (i = 0; i < ldv * n; i++)


}
