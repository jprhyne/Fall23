*
*  This subroutine will compute the null space of A\in\R^{m by n} and m <=n. 
*  The general overview of how this will work is by 
*  (1) Compute the singular values of A
*  (2) Determine how many 
*
*
      SUBROUTINE(M,N,A,LDA,EPS,INFO)
*
*        Arguments 
*
         INTEGER M,N,LDA, INFO
         DOUBLE PRECISION EPS
         DOUBLE PRECISION A(LDA,*)
*
*        Local scalars
*
         INTEGER LWORK, I, J, K
         DOUBLE PRECISION S(N)
*        allocatable array to be used for workspaces
*        This will be refactored out once we know all the functions we
*        are going to be calling
         real*8, allocatable :: WORK(:)

*
*        External functions
*
         EXTERNAL DGESVD 

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
         DGESVD('No U', 'No V', M, N, A, LDA, S, A, LDA, A, LDA, WORK,
     $            LWORK, INFO)
         DEALLOCATE(WORK)
         LWORK = WORK(1)
         ALLOCATE(WORK(LWORK))
*
*        Now, we actually call DGESVD to compute the singular values of
*        A.
*
         DGESVD('No U', 'No V', M, N, A, LDA, S, A, LDA, A, LDA, WORK,
     $            LWORK, INFO)
*
*        Determine how many singular values are less than EPS. This will
*        be how many extra columns we want to grab
*
*        TODO: Actually use this value. Right now we just compute it and
*        ignore it by assuming that A is full rank
*
         DO 10 I=N, 1, -1
            IF(S(I).LT.EPS) THEN
               K = I
            END IF
   10    CONTINUE
*
*        If A was full rank, we just need to create Q2
*
         DO 10 I=N, 1, -1
         DO 10 I=N, 1, -1
         DO 10 I=N, 1, -1
         asdf

      END
