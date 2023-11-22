/*
 *  This file is used to demonstrate that there is a bug in AOCL DORG2R
 */

#include <stdio.h>
#include<stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

#ifdef USE_MKL_BLAS
    #include <mkl_cblas.h>
#else
    #include <cblas.h>
#endif

#ifdef USE_MKL_LAPACK
    #include <mkl_lapacke.h>
#else
    #include <lapacke.h>
#endif

int main(int argc, char **argv) {

    // Local params
    int info, lda, ldq, m, n, k, lwork, nb;
    double *A, *Q, *As, *tau, *work=NULL;
    double normA;
    double elapsed_refL, perform_refL;
    struct timeval tp;
    int i, j;

    // Init random seed
    //srand(0);

    // Default parameters that demonstrate the bugged behavior.
    // I have found this works for any square matrix

    m = 4;
    n = 4;
    k = n;
    lda = -1;
    ldq = -1;
    // Flag that helps facilitate testing with the driver.c file
    bool errorsOnly = false;
    // Flag that helps facilitate timeing with the time.sh file
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
    }

    if( lda < 0 ) lda = m;
    if( ldq < 0 ) ldq = m;

    if (!errorsOnly && !timesOnly) {
        printf("dgeqrf dorg2r LAPACK | ");
        printf("m = %4d, ",    m);
        printf("n = %4d, ",    n);
        printf("k = %4d, ",    k);
        printf("lda = %4d, ",lda);
        printf("ldq = %4d, ",ldq);
        printf("             ");
    }

    A  = (double *) calloc(lda * k,sizeof(double));
    As = (double *) calloc(lda * k, sizeof(double));
    Q  = (double *) calloc(ldq * n, sizeof(double));

    for(i = 0; i < lda * k; ++i)
        *(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    for(i = 0; i < ldq * n; ++i)
        *(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, As, lda );
    normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, A, lda, work );

    tau = (double *) malloc( k * sizeof(double));

    work = (double *) malloc( 1 * sizeof(double));
    LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, A, lda, tau, work, -1 ); 
    lwork = ((int) work[0]);
    if (lwork < n) lwork = n; 
    free( work );
    work = (double *) malloc( lwork * sizeof(double));

    LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, A, lda, tau, work, lwork ); 
    LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, Q, ldq );
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Calls dorg2r with the matrix Q.
    dorg2r_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);

    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    free( work );
    // fix later. This was from timing everything
    perform_refL = (4.0e+00*(double) m*(double)n*(double)k -2.0e+00*(double)m*(double)k*(double)k - 2.0e+00*(double)n*(double)k*(double)k + 4.0e+00/3.0e+00 *(double)k*(double)k*(double)k)/elapsed_refL /1.0e+9;
    
    double norm_orth_1, norm_repres_1;

    work  = (double *) malloc(n * n * sizeof(double));
    info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
    cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
    norm_orth_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
    free( work );

    work  = (double *) malloc(m * k * sizeof(double));
    info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, Q, ldq, work, m );
    cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, k, (1.0e+00), A, lda, work, m );
    for(i = 0; i < m; ++i)
        for(j = 0; j < k; ++j)
            work[ i+j*m ] -= As[ i+j*lda ];
    norm_repres_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, work, m, NULL );
    norm_repres_1 = norm_repres_1 / normA;
    free( work );

    if (!errorsOnly && !timesOnly) {
        printf("| time = %f   GFlop/sec = %f", elapsed_refL, perform_refL);
        printf("| repres  = %5.1e    ortho = %5.1e ", norm_repres_1, norm_orth_1);
    } else if (errorsOnly && !timesOnly) {
        printf("%10.10e %10.10e", norm_repres_1, norm_orth_1);
    } else {
        printf("%f", elapsed_refL);
    }

    printf("\n");

    free( Q );
    free( A );
    free( As );
    free( tau );

    return 0;
}
