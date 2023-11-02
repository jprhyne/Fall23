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
    char JobU = 'N';
    char JobVT= 'N';
    double workQ[1];
    int lwork = -1;
    info = 0;
    dgesvd_(&JobU, &JobVT, &m, &n, A, &lda, S, NULL, &m, NULL, &n, workQ, &lwork, &info);
    // allocate the workspace now
    double *work = (double *) malloc(workQ[0] * sizeof(double));
    lwork = (int) workQ[0];
    dgesvd_(&JobU, &JobVT, &m, &n, A, &lda, S, NULL, &m, NULL, &n, work, &lwork, &info);

    for (int i = 0; i < n; i++)
        printf("%lf\n",S[i]);
    printf("\n\n\n");
    for (int i = 0; i < n; i++)
        printf("%lf\n",sigma[i] - S[i]);


    free(A);
    free(sigma);
    free(S);
    free(work);
}
