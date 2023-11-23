/*
 *  This file is used to demonstrate that there is a bug in AOCL DORG2R
 *  But this is only the distrubuted version given at https://www.amd.com/en/developer/aocl/dense.html#lapack
 *  on the date of 23 November 2023
 */

// Facilitates the printing of the source LAPACK implimentation we are using for printing out to the user.
#ifdef USE_AOCL
    #define SOURCE "AOCL"
#else
    #define SOURCE "REF_LAPACK"
#endif

#include <stdio.h>
#include<stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <cblas.h>
#include <lapacke.h>

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
    // I have found this works for any square matrix for dorg2r
    // and any time where a call to dorg2r from would result in a square 
    // matrix when calling my implimentation of dorgqr.

    m = 4;
    n = 4;
    k = n;
    for(i = 1; i < argc; ++i){
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
    }

    lda = m;
    ldq = m;
    char *source = SOURCE;
    printf("dgeqrf dorg2r %s | m = %4d, n = %4d, k = %4d, \t", source, m, n, k);

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
    // Computing error and timing information
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

    printf("| time = %f GFlop/sec = %f|\nrepres = %5.1e ortho = %5.1e\n", 
            elapsed_refL, perform_refL, norm_repres_1, norm_orth_1);


    // Now we repeat above, BUT use dorgqr instead of dorg2r. This helps demonstrate our issue
    printf("dgeqrf dorgqr %s | m = %4d, n = %4d, k = %4d, \t", source, m, n, k);
    // Copy the saved version of A back into A
    info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, As, lda, A, lda );
    normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, k, A, lda, work );
    // Free tau
    free( tau );

    tau = (double *) malloc( k * sizeof(double));

    work = (double *) malloc( 1 * sizeof(double));
    LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, A, lda, tau, work, -1 ); 
    lwork = ((int) work[0]);
    LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, k, A, lda, tau, work, -1 );
    if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
    free( work );
    work = (double *) malloc( lwork * sizeof(double));

    LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, k, A, lda, tau, work, lwork ); 
    LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, k, A, lda, Q, ldq );
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Calls dorg2r with the matrix Q.
    dorgqr_(&m, &n, &k, Q, &ldq, tau, work, &lwork, &info);

    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    free( work );
    // Computing error and timing information
    perform_refL = (4.0e+00*(double) m*(double)n*(double)k -2.0e+00*(double)m*(double)k*(double)k - 2.0e+00*(double)n*(double)k*(double)k + 4.0e+00/3.0e+00 *(double)k*(double)k*(double)k)/elapsed_refL /1.0e+9;
    
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

    printf("| time = %f GFlop/sec = %f|\nrepres = %5.1e ortho = %5.1e\n", 
            elapsed_refL, perform_refL, norm_repres_1, norm_orth_1);
    free( Q );
    free( A );
    free( As );
    free( tau );

    return 0;
}
