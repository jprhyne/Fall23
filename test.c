#include <stdio.h>
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
    int info, lda, ldq, m, n, lwork, nb;
    double *A, *Q, *As, *tau, *work=NULL;
    double normA;
    double elapsed_refL, perform_refL;
    struct timeval tp;
    int i, j;

    // Init random seed
    srand(0);

    // Default parameters to ensure that we hit 
    // the blocked code
    m = 30;
    n = 20;
    nb = 7; // Choose a default nb value that is NOT a factor of m nor n
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
    }

    if( lda < 0 ) lda = m;
    if( ldq < 0 ) ldq = m;

    printf("dgeqrf dgorgqr LAPACK | ");
    printf("m = %4d, ",    m);
    printf("n = %4d, ",    n);
    printf("lda = %4d, ",lda);
    printf("ldq = %4d, ",ldq);
    printf("           ");
    printf("  ");

    A  = (double *) calloc(lda * n,sizeof(double));
    As = (double *) calloc(lda * n, sizeof(double));
    Q  = (double *) calloc(ldq * n, sizeof(double));

    for(i = 0; i < lda * n; ++i)
        *(A + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    for(i = 0; i < ldq * n; ++i)
        *(Q + i) = (double)rand() / (double)(RAND_MAX) - 0.5e+00;

    info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, As, lda );
    normA = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, A, lda, work );

    tau = (double *) malloc( n * sizeof(double));

    work = (double *) malloc( 1 * sizeof(double));
    LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, -1 ); 
    lwork = ((int) work[0]);
    LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, A, lda, tau, work, -1 );
    if (lwork < ((int) work[0])) lwork = ((int) work[0]); 
    free( work );
    work = (double *) malloc( lwork * sizeof(double));

    int k = n / 2; // integer division
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    LAPACKE_dgeqrf_work( LAPACK_COL_MAJOR, m, n, A, lda, tau, work, lwork ); 
    LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, A, lda, Q, ldq );
    //LAPACKE_dorgqr_work( LAPACK_COL_MAJOR, m, n, n, Q, ldq, tau, work, lwork );
    // Directly calling the fortran function dorgqr
    //dorgqr_(&m, &n, &n, Q,&ldq, tau, work, &lwork, &info);
    // Directly calling my fortran function my_dorgqr
    my_dorgqr_(&m, &n, &n, &nb, Q, &ldq, tau, work, &lwork, &info);

    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    free( work );

    perform_refL = ( 4.0e+00 * ((double) m) * ((double) n) * ((double) n) - 4.0e+00 / 3.0e+00 * ((double) n) * ((double) n) * ((double) n) )  / elapsed_refL / 1.0e+9 ;
    
    double norm_orth_1, norm_repres_1;

    lwork = n*n;
    work  = (double *) malloc(n * n * sizeof(double));
    info  = LAPACKE_dlaset( LAPACK_COL_MAJOR, 'A', n, n, (0e+00), (1e+00), work, n );
    cblas_dsyrk( CblasColMajor, CblasUpper, CblasTrans, n, m, 1.0e+00, Q, ldq, -1.0e+00, work, n );
    norm_orth_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', n, n, work, n, NULL );
    free( work );

    lwork = m*n;
    work  = (double *) malloc(m * n * sizeof(double));
    info  = LAPACKE_dlacpy_work( LAPACK_COL_MAJOR, 'A', m, n, Q, ldq, work, m );
    cblas_dtrmm( CblasColMajor, CblasRight, CblasUpper, CblasNoTrans, CblasNonUnit, m, n, (1.0e+00), A, lda, work, m );
    for(i = 0; i < m; ++i)
        for(j = 0; j < n; ++j)
            work[ i+j*m ] -= As[ i+j*lda ];
    norm_repres_1 = LAPACKE_dlange_work( LAPACK_COL_MAJOR, 'F', m, n, work, m, NULL );
    norm_repres_1 = norm_repres_1 / normA;
    free( work );

    printf("| time = %f   GFlop/sec = %f", elapsed_refL, perform_refL);

    printf("\n");
    printf("| repres  = %5.1e    ortho = %5.1e ", norm_repres_1, norm_orth_1);

    printf("\n");


    free( Q );
    free( A );
    free( As );
    free( tau );

    return 0;
}
