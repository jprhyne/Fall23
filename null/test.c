#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
int main()
{
    int m,n,lda,info;
    // Default values for m,n,lda
    m = 5;
    n = 5;
    lda = 5;
    // Create A as an array of size m*n
    double *A = (double *) malloc(m*n * sizeof(double));
    // Create sigma as an array representing a diagonal matrix
    double *sigma = (double *) malloc(n * sizeof(double));
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1e-5;
    sigma[4] = 1e-7;
    create_dmat_(&m,&n,sigma,A,&lda,&info);
    // Find the svd of A to confirm we did it correctly
    double *S = (double *) malloc(n*sizeof(double));
    char JobU = 'T';
    char JobVT= 'N';
    double workQ[1];
    int lwork = -1;
    info = 0;
    /*
    dgesvd_(&JobU, &JobVT, &m, &n, A, &lda, S, NULL, &m, NULL, &n, workQ, &lwork, &info);
    // allocate the workspace now
    double *work = (double *) malloc(workQ[0] * sizeof(double));
    lwork = (int) workQ[0];
    dgesvd_(&JobU, &JobVT, &m, &n, A, &lda, S, NULL, &m, NULL, &n, work, &lwork, &info);
    */
    // Now we call our nullspace computation routine
    double eps = 1e-6;
    // Copy A into a new matrix to be used for storage
    double *ABak = (double *) malloc(m*n*sizeof(double));
    for (int i = 0; i < m * n; i++)
        ABak[i] = A[i];
    nulls_(&m,&n,&A,&m,&eps,&info);
    double *C = (double *) malloc(m*n*sizeof(double));
    double one = 1.0;
    double zero= 0.0;
    int k = n - info;
    for (
    //dgemm_(&JobU, &JobVT, &m, &n, &k, &one, ABak, &m, A, &m, &zero, C, &m);
    free(A);
    free(ABak);
    free(C);
    free(sigma);
    free(S);
    //free(work);
}
