* This subroutine will be used to essentially do a dtrmm('L','U','N','N',...) but not in
* place
* We compute C = T*A + ALPHA*C
*     where T is upper unit triangular and A is rectangular
* Currently do not support any other functionality, but can if desired
      SUBROUTINE DTRMMOOP(M, N, A, LDA, T, LDT, ALPHA, C, LDC)
*
*        .. Scalar Arguments ..
*
         INTEGER           M, N, LDA, LDT, LDC

         DOUBLE PRECISION  ALPHA

*
*        .. Array arguments ..
*
         DOUBLE PRECISION  A(LDA,*), T(LDA,*), C(LDA,*)

*
*        .. Local variables ..
*
         INTEGER           I, J, K

         DOUBLE PRECISION  TEMP

*--------------------------------------------------------------------------
*        Begin of executable statements
*--------------------------------------------------------------------------
*        quick return if possible
         IF (M.EQ.0 .OR. N.EQ.0) RETURN
*
*        The formula we are using is for all valid i,j 
*        C(i,j) = sum_1^jA(i,k)*T(k,j) + alpha*C(i,j)
*
         DO 50 J = 1,N
            DO 40 I = 1,M
               TEMP = A(I,J)
               DO 30 K = I+1, M
                  TEMP = TEMP + T(I,K)*A(K,J)
   30          CONTINUE
*              Note: This is done at the end to make comparison with
*              existing functionality easier
*              if alpha is 1, we are adding to the existing element
               IF (ALPHA.EQ.1) THEN
                  TEMP = TEMP + C(I,J)
*              if alpha neither 1 nor 0, we are adding to a non-trivial
*              multiple of the current element
               ELSE IF (ALPHA.NE.0) THEN
                  TEMP = TEMP + ALPHA * C(I,J)
               END IF 
               C(I,J) = TEMP
   40       CONTINUE
   50    CONTINUE

      END SUBROUTINE 
