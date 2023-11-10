*
*  This subroutine will compute the null space of A\in\R^{m by n} and m <=n. 
*  The general overview of how this will work is by 
*  (1) Compute the singular values of A
*  (2) Determine how many 
*
*
      SUBROUTINE NULLS(M,N,A,LDA,EPS,INFO)
*
*        Arguments 
*
         INTEGER M,N,LDA, INFO
         DOUBLE PRECISION EPS
         DOUBLE PRECISION A(LDA,*)
*
*        Local scalars
*
         INTEGER LWORK, I, J, K, NB
         DOUBLE PRECISION S(N), TAU(N)
*        allocatable array to be used for workspaces
*        This will be refactored out once we know all the functions we
*        are going to be calling
         real*8, allocatable :: WORK(:)

*
*        External functions
*
         EXTERNAL DGESVD 
*     .. External Functions ..
      INTEGER           ILAENV
      EXTERNAL          ILAENV

*
*        Argument checking
*
         IF( M.LT.N ) THEN
            INFO = -1
         ELSE IF (M.GT.LDA) THEN
            INFO = -2
         ELSE IF (EPS.LT.0) THEN
            INFO = -3
         ELSE IF (N.LT.0) THEN
            INFO = -4
         END IF

         IF(INFO.NE.0) THEN
            RETURN
         END IF
*
*        Beginning of executable statements
*

*
*        Compute our singular values of A
*        Work query
*
         LWORK=-1
         ALLOCATE(WORK(1))
         CALL DGESVD('No U', 'No V', M, N, A, LDA, S, A, LDA, A, LDA,
     $            WORK, LWORK, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
*
*        Now, we actually call DGESVD to compute the singular values of
*        A.
*
         CALL DGESVD('No U', 'No V', M, N, A, LDA, S, A, LDA, A, LDA,
     $            WORK, LWORK, INFO)
*
*        Determine how many singular values are less than EPS. This will
*        be how many extra columns we want to grab
*
*        TODO: Actually use this value. Right now we just compute it and
*        ignore it by assuming that A is full rank
*
         IF(S(N).LT.EPS) THEN
            DO 10 I=N, 1, -1
               IF(S(I).LT.EPS) THEN
                  K = I
               END IF
   10       CONTINUE
         ELSE
            K = N
         END IF
*        Find the QR decomposition of A represented as householder
*        reflectors. 
*        Workspace query
         LWORK=-1
         CALL DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)

         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))
         CALL DGEQRF(M, N, A, LDA, TAU, WORK, LWORK, INFO)
*        After this loop, we will know the first index k such that the
*        k^th singular value is below epsilon. If K = 0 then A is full
*        rank.
*
         NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
         DEALLOCATE(WORK)
         ALLOCATE(WORK(NB*NB))
         CALL COMPUTEQ2(M,N,K,A,LDA,WORK, NB, TAU)
         INFO = K
      END
