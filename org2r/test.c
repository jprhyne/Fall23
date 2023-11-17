#include <stdio.h>
#include<stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <cblas.h>
#include <lapacke.h>

int main(int argc, char *argv[]) {
    int info, lda, ldq, m, n, k, lwork, nb;
    double *A, *Q, *As, *tau, *work=NULL;
    double normA;
    double elapsed_refL, perform_refL;
    struct timeval tp;
    int i, j;

    m = 30;
    n = 20;
    k = n/2 + 1;
    nb = 3; // Choose a default nb value that is NOT a factor of k
    lda = -1;
    ldq = -1;

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
    }

    if( lda < 0 ) lda = m;
    if( ldq < 0 ) ldq = m;
}
