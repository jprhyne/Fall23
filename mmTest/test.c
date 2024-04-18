#include <stdio.h>
#include<stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

int main(int argc, char **argv) {

    // Local params
    int info, lda, ldq, m, n, k, lwork, nb;
    double elapsed_refL, perform_refL;
    struct timeval tp;
    int i, j, version;

    char aChar = 'A';
    char fChar = 'F';
    char uChar = 'U';
    char tChar = 'T';
    char rChar = 'R';
    char nChar = 'N';
    int dummy = 0;

    int negOne = -1;
    double one = 1.0;
    double zero = 0.0;
    double dNegOne = -1.0;

    // Init random seed
    //srand(0);

    // Default parameters to ensure that we hit 
    // the blocked code
    m = 30;
    n = 20;
    k = n/2 + 1;
    nb = 3; // Choose a default nb value that is NOT a factor of k
    lda = -1;
    ldq = -1;
    version = -1; // Default to the system version (vendor usually)
    // Flag that helps facilitate testing with the driver.c file
    bool errorsOnly = false;
    // Flag that helps facilitate timing with the time.sh file
    bool timesOnly = false;

    for(i = 1; i < argc; ++i){
        if( strcmp( *(argv + i), "-ldq") == 0) {
            ldq  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-lda") == 0) {
            lda  = atoi( *(argv + i + 1) );
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
        if( strcmp( *(argv + i), "-nb") == 0) {
            nb  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-k") == 0) {
            k  = atoi( *(argv + i + 1) );
            i++;
        }
        if( strcmp( *(argv + i), "-e") == 0) {
            errorsOnly  = true;
        }
        if( strcmp( *(argv + i), "-t") == 0) {
            timesOnly  = true;
        }
        if( strcmp( *(argv + i), "-v") == 0) {
            version = atoi( *(argv + i + 1) );
            i++;
        }
    }

    if( lda < 0 ) lda = m;
    if( ldq < 0 ) ldq = m;
    if( k > n) k = n;

    if (!errorsOnly && !timesOnly) {
        printf("dgeqrf dorgqr | ");
        printf("version = %4d, ", version);
        printf("m = %4d, ",    m);
        printf("n = %4d, ",    n);
        printf("k = %4d, ",    k);
        printf("nb = %4d, ",    nb);
        printf("lda = %4d, ",lda);
        printf("ldq = %4d, ",ldq);
        printf("\n");
    }

    double *A  = (double *) calloc(m * n,sizeof(double));
    double *B  = (double *) calloc(n * k, sizeof(double));
    double *C  = (double *) calloc(m * k, sizeof(double));

    for(i = 0; i < m * n; ++i)
        *(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    for(i = 0; i < n * k; ++i)
        *(B + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    dgemm_(&nChar, &nChar, &m, &n, &k, &one, A, &m, B, &n, &zero, C, &m, dummy, dummy);

    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    if (!errorsOnly && !timesOnly) {
        printf("| time = %f   GFlop/sec = %f", elapsed_refL, perform_refL);
    } else {
        printf("%10.10e:%10.10e", elapsed_refL, perform_refL);
    }

    printf("\n");

    free( A );
    free( B );
    free( C );

    return 0;
}
