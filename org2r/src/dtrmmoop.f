* This subroutine will be used to essentially do a dtrmm('L','N','N','N',...) but not in
* place on the transpose of matrix A
* We compute C = T*A**T + C
*     where T is upper unit triangular and A is rectangular
* Currently do not support any other functionality, but can if desired
*     C is m by n
*     T is m by m
*     A is n by m -> A**T is m by n
      RECURSIVE SUBROUTINE DTRMMOOP(M, N, A, LDA, T, LDT, C, LDC)
*
*        .. Scalar Arguments ..
*
         INTEGER           M, N, LDA, LDT, LDC

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
*        Base cases
         IF (M.EQ.1) THEN
*           In this case, T is 1x1 upper triangular matrix.
*           Therefore, we need to compute C = C + A*T(1,1)
*
*           This special case is done because for some reason, when we go to the
*           10 do loop, we have j=1,1 and this somehow goes to j=2. Not sure
*           why, but this is a workaround. (maybe some gdb issue. will toy
*           around with removing later as a last cleanup step).
*
            IF (N.EQ.1) THEN
               C(1,1) = C(1,1) + A(1,1)*T(1,1)
            ELSE
               DO 10 J = 1, N
                  C(1,J) = C(1,J) + A(J,1)*T(1,1)
   10          CONTINUE
            END IF
            RETURN
         ELSE IF (N.EQ.1) THEN
*           In this case, we have a C as a column vector and we need to compute
*           a modified matrix vector product (ie a modified form of dtrmv) But
*           also out of place.
*           We accomplish this by computing each element of C through ddot.
            DO 20 I = 1,M
               C(I,1) = C(I,1) + DDOT(M-I+1, T(I,I), LDT, A(1,I),LDA)
*               C(I,1) = C(I,1) + A(1, I)
*               IF (I.LT.M) THEN
*                  C(I,1) = C(I,1) + DDOT(M-I, T(I,I+1), LDT, 
*     $                                   A(1,I+1),LDA)
*               END IF
   20       CONTINUE
            RETURN
         END IF
*        Recursive case
         K = MIN(M,N) / 2
*        Compute C21
         CALL DTRMMOOP(M - K, K, A(1, K+1), LDA, T(K + 1, K + 1), 
     $            LDT, C(K + 1, 1), LDC)
*        Compute C22
         CALL DTRMMOOP(M - K, N - K, A(K+1,K+1), LDA, T(K + 1, K + 1),
     $           LDT, C(K + 1, K + 1), LDC)
*        Compute C11 part 1
         CALL DTRMMOOP(K, K, A, LDA, T, LDT, C, LDC)
*        Compute C11 part 2
         CALL DGEMM('No transpose', 'Transpose', K, K, M - K, ONE, 
     $           T(1, K + 1), LDT, A(1, K + 1), LDA, ONE, C, LDC)
*        Compute C12 part 1
         CALL DTRMMOOP(K, N-K, A(K + 1, 1), LDA, T, LDT, C(1, K + 1),
     $           LDC)
*        Compute C12 part 2
         CALL DGEMM('No transpose', 'Transpose', K, N - K, M - K,
     $           ONE, T(1, K + 1), LDT, A(K + 1, K + 1), LDA,
     $           ONE, C(1, K + 1), LDC)
      END SUBROUTINE 
