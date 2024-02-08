* We compute C = T*A**T + ALPHA*C or C = A**T*T + ALPHA * C 
* If SIDE = 'L'/'l' or 'R'/'r' respectively
*     where T is upper triangular and A is rectangular
* Currently do not support any other functionality, but can if desired
*     C is m by n
*     T is m by m
*     A is n by m -> A**T is m by n
      RECURSIVE SUBROUTINE DTRMMOOP(SIDE, DIAG, M, N, A, LDA, T, LDT,
     $                              ALPHA, C, LDC, INFO)
*
*        .. Scalar Arguments ..
*
         INTEGER           M, N, LDA, LDT, LDC, INFO
         CHARACTER         SIDE, DIAG
         DOUBLE PRECISION  ALPHA

*
*        .. Array arguments ..
*
         DOUBLE PRECISION  A(LDA,*), T(LDT,*), C(LDC,*)
*
*        .. Local variables ..
*
         INTEGER           I,J,K
         DOUBLE PRECISION  SUM
*
*        .. Local parameters ..
*
         DOUBLE PRECISION  ONE, ZERO
         PARAMETER(ONE=1.0d+0,ZERO=0.0d+0)
*
*        .. External Functions ..
*
         DOUBLE PRECISION  DDOT
         EXTERNAL          DDOT
*
*        .. External Subroutines ..
*
         EXTERNAL          DAXPY
*        We break down each matrix into the following form
*
*        |-------|   |-------|   |-------|   |-------|**T
*        |C11 C12| = |C11 C12| + |T11 T12| * |A11 A12|
*        |C21 C22|   |C21 C22|   |T21 T22|   |A21 A22|
*        |-------|   |-------|   |-------|   |-------|
*
*
*
*        C_{1,1}\in\R^{k by k}
*        C_{1,2}\in\R^{k by n - k}
*        C_{2,1}\in\R^{m - k by k}
*        C_{2,2}\in\R^{m - k by n - k}
*        
*        We choose K = MIN(M,N)/2
*
*        C_{2,2} = C_{2,2} + T_{2,2}*A_{2,2}**T
*        C_{2,1} = C_{2,1} + T_{2,2}*A_{1,2}**T
*
*        C_{1,1} = C_{1,1} + T_{1,1}*A_{1,1}**T + T_{1,2}*A_{1,2}**T
*        C_{1,2} = C_{1,2} + T_{1,1}*A_{2,1}**T + T_{1,2}*A_{2,2}**T
*
*        C_{1,1} can be broken down into two different operations
*        C_{1,1} = C_{1,1} + T_{1,2}*A_{1,2}**T [GEMM]
*        C_{1,1} = C_{1,1} + T_{1,1}*A_{1,1}**T [DTRMMOOP]
*
*        C_{1,2} can be broken down into two different operations
*        C_{1,2} = C_{1,2} + T_{1,2}*A_{2,2}**T [GEMM]
*        C_{1,2} = C_{1,2} + T_{1,1}*A_{2,1}**T [DTRMMOOP]

*--------------------------------------------------------------------------
*        Begin of executable statements
*--------------------------------------------------------------------------
*        Base case for the left half of the matrix
         IF (M.EQ.0.OR.N.EQ.0) THEN
            RETURN
         END IF
*        Determine if we have T on the left or right
         IF(SIDE.EQ.'L'.OR.SIDE.EQ.'l') GOTO 10
         IF(SIDE.EQ.'R'.OR.SIDE.EQ.'r') GOTO 20
         INFO = -1
         RETURN
*        Base cases
   10    IF (M.EQ.1) THEN
*           In this case, T is 1x1 upper triangular matrix.
*           Therefore, we need to compute C = C + A*T(1,1)
*
*           This special case is done because for some reason, when we go to the
*           10 do loop, we have j=1,1 and this somehow goes to j=2. Not sure
*           why, but this is a workaround. (maybe some gdb issue. will toy
*           around with removing later as a last cleanup step).
*
            ! CALL DAXPY(N, T(1,1), A(1,1), 1, C(1,1), LDC)
            IF (ALPHA.NE.ONE) THEN
               DO J=1, N
                  C(1,J) = ALPHA * C(1,J)
               END DO
            END IF
            CALL DAXPY(N, T(1,1),A(1,1), 1, C(1,1), LDC)
            !IF (N.EQ.1) THEN
            !   C(1,1) = C(1,1) + A(1,1)*T(1,1)
            !ELSE
            !END IF
            RETURN
         ELSE IF (N.EQ.1) THEN
            ! Write quick implementation of dtrmvoop (should be similar to this
            ! case)
*           In this case, we have a C as a column vector and we need to compute
*           a modified matrix vector product (ie a modified form of dtrmv) But
*           also out of place.
*           We accomplish this by computing each element of C through ddot.
            
            DO I = 1,M
               C(I,1) = ALPHA * C(I,1) + DDOT(M-I+1, T(I,I), LDT, 
     $                     A(1,I),LDA)
*               C(I,1) = C(I,1) + A(1, I)
*               IF (I.LT.M) THEN
*                  C(I,1) = C(I,1) + DDOT(M-I, T(I,I+1), LDT, 
*     $                                   A(1,I+1),LDA)
*               END IF
            END DO
            RETURN
         END IF
*        Recursive case
         K = MIN(M,N) / 2
*        Compute C21
         CALL DTRMMOOP(SIDE, DIAG, M - K, K, A(1, K+1), LDA,
     $            T(K + 1, K + 1), LDT, ALPHA, C(K + 1, 1), LDC, INFO)
*        Compute C22
         CALL DTRMMOOP(SIDE, DIAG, M - K, N - K, A(K+1,K+1), LDA,
     $           T(K+1,K+1), LDT, ALPHA, C(K + 1, K + 1), LDC, INFO)
*        Compute C11 part 1
         CALL DTRMMOOP(SIDE, DIAG, K, K, A, LDA, T, LDT, ALPHA, C, LDC, 
     $                  INFO)
*        Compute C11 part 2
         CALL DGEMM('No transpose', 'Transpose', K, K, M - K, ONE, 
     $           T(1, K + 1), LDT, A(1, K + 1), LDA, ONE, C, LDC)
*        Compute C12 part 1
         CALL DTRMMOOP(SIDE, DIAG, K, N-K, A(K + 1, 1), LDA, T, LDT, 
     $           ALPHA, C(1, K + 1), LDC, INFO)
*        Compute C12 part 2
         CALL DGEMM('No transpose', 'Transpose', K, N - K, M - K,
     $           ONE, T(1, K + 1), LDT, A(K + 1, K + 1), LDA,
     $           ONE, C(1, K + 1), LDC)
         INFO = 0
         RETURN
   20    IF (M.EQ.1) THEN
            IF (ALPHA.NE.ONE) THEN
               CALL DSCAL(N, ALPHA, C, LDC)
            END IF
            IF (DIAG.EQ.'U'.OR.DIAG.EQ.'u') THEN
               C(1,1) = C(1,1) + A(1,1)
               DO J = 2, N
                  C(1,J) = C(1,J) + DDOT(J-1, A, 1, T(1,J), 1)
                  C(1,J) = C(1,J) + A(J,1)
               END DO
            ELSE
               DO J = 1, N
                  C(1,J) = C(1,J) + DDOT(J, A, 1, T(1,J), 1)
               END DO
            END IF
            RETURN
         END IF
         IF (N.EQ.1) THEN 
         ! M.EQ.1 for above
            IF (ALPHA.NE.ONE) THEN
               CALL DSCAL(M, ALPHA, C, 1)
            END IF
            IF (DIAG.EQ.'U'.OR.DIAG.EQ.'u') THEN
               CALL DCOPY(M, A, LDA, C, 1)
            ELSE
               CALL DAXPY(M, T(1,1), A(1,1), LDA, C(1,1), 1)
            END IF
            RETURN
         END IF
         K = MIN(M,N) / 2
*        Compute C11
         CALL DTRMMOOP(SIDE, DIAG, K, K, A, LDA, T, LDT, ALPHA, C, 
     $                  LDC, INFO)
*        Compute C21
         CALL DTRMMOOP(SIDE, DIAG, M-K, K, A(1,K+1), LDA, T, LDT,
     $                  ALPHA, C(K+1, 1), LDC, INFO)
*        Compute C12
*        Part 1
         CALL DTRMMOOP(SIDE, DIAG, K, N-K, A(K+1,1), LDA, T(K+1,K+1),
     $                  LDT, ALPHA, C(1, K+1), LDC, INFO)

*        Part 2
         CALL DGEMM('Transpose', 'No transpose', K, N-K, K, ONE, A, LDA,
     $                  T(1,K+1), LDT, ONE, C(1, K+1), LDC)
*        Compute C22
*        Part 1
         CALL DTRMMOOP(SIDE, DIAG, M-K, N-K, A(K+1,K+1), LDA,
     $                  T(K+1,K+1), LDT, ALPHA, C(K+1,K+1), LDC, INFO)

*        Part 2
         CALL DGEMM('Transpose', 'No transpose', M-K, N-K, K, ONE, 
     $                  A(1, K+1), LDA, T(1, K+1), LDT, ONE, C(K+1,K+1),
     $                  LDC)
         INFO = 0
         RETURN
      END SUBROUTINE
