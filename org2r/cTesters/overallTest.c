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
#include <cblas.h>
#include <lapacke.h>

int main(int argc, char *argv[]) {
    // integer variables
    int info, lda, ldq, m, n, k, lwork, nb, i, j;
    int workQuery = -1;
    // double variables
    double *A, *Q, *As, *tau, *work=NULL;
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
    k = n/2 + 1;
    nb = 3; // Choose a default nb value that is NOT a factor of k
    while (k % nb == 0) {
        nb++;
    }
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

    // allocate memory for the matrices and vectors that we need to use
    A =   (double *) malloc(lda * k * sizeof(double));
    As =  (double *) malloc(lda * k * sizeof(double));
    Q =   (double *) malloc(ldq * n * sizeof(double));
    tau = (double *) malloc(k * sizeof(double));

    // Print to the user what we are doing along with any arguments that are used
    char *source = SOURCE;
    printf("dgeqrf dorg2r %s | m = %4d, n = %d, k = %4d, lda = %4d, ldq = %4d\n", source, m, n, k, lda, ldq);

    // Generate A as a random matrix
    for (i = 0; i < lda * k; i++)
        A[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;
    // Store random data inside Q to ensure that we do not assume anything about Q
    for (i = 0; i < lda * n; i++)
        Q[i] = (double) rand() / (double) (RAND_MAX) - 0.5e+00;

    // Store a copy of A inside As
    dlacpy_(&aChar, &m, &k, A, &lda, As, &lda, dummy);
    // Find the norm of A for use in later accuracy computations
    normA = dlange_(&fChar, &m, &k, A, &lda, NULL, dummy);
    // Create the work array to do workspace queries
    work = (double *) malloc(sizeof(double));
    // Determine how much workspace is needed for our operations
    // Check dgeqrf first
    dgeqrf_(&m, &k, A, &lda, tau, work, &workQuery, &info );
    lwork = work[0];
    // dorg2r uses a workspace of size  "n" (k) so check if dgeqrf requires more
    if (lwork < n)
        lwork = n;

    // reallocate work to be of the right size
    work = (double *) realloc(work, lwork * sizeof(double));

    // Call dgeqrf first
    dgeqrf_(&m, &k, A, &lda, tau, work, &lwork, &info);

    // Copy A into Q for use with dorg2r
    dlacpy_(&aChar, &m, &k, A, &lda, Q, &ldq, dummy);
    
    // Take the current time for use with timing dorg2r
    gettimeofday(&tp, NULL);
    elapsed_refL=-((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Call dorg2r.
    dorg2r_(&m, &n, &k, Q, &ldq, tau, work, &info);

    // Determine how much time has taken
    gettimeofday(&tp, NULL);
    elapsed_refL+=((double)tp.tv_sec+(1.e-6)*tp.tv_usec);

    // Compute the error information
    // reallocate work to be of the right size for testing orthogonality of Q.
    work = (double *)realloc(work, n * n * sizeof(double));

    // Set work to be I
    dlaset_(&aChar, &n, &n, &zero, &one, work, &n, dummy);
    // Compute work = Q**T * Q - I
    dsyrk_(&uChar, &tChar, &n, &m, &one, Q, &ldq, &negOne, work, &n);
    // Compute the norm of work
    norm_orth_1 = dlange_(&fChar, &n, &n, work, &n, NULL, dummy);

    // reallocate work to be able to hold Q
    work = (double *)realloc(work, m * k * sizeof(double));
    // Copy Q into work
    dlacpy_(&aChar, &m, &k, Q, &ldq, work, &m, dummy);

    // Compute work = work * R
    dtrmm_(&rChar, &uChar, &nChar, &nChar, &m, &k, &one, A, &lda, work, &m, dummy, dummy, dummy, dummy);
    
    // Compute work = work - A
    for (i = 0; i < m; i++)
        for (j = 0; j < k; j++)
            work[i + j * m] -= As[i + j * lda];
    // Compute ||A - QR||_F
    norm_repres_1 = dlange_(&fChar, &m, &k, work, &m, NULL, dummy);

    // Compute ||A - Q*R||_F / ||A||_F
    norm_repres_1 /= normA;

    printf("time = %f repres = %5.1e ortho = %5.1e\n", elapsed_refL, norm_repres_1, norm_orth_1);
    // Free memory before terminating 
    free(Q);
    free(As);
    free(A);
    free(tau);
    free(work);
}
