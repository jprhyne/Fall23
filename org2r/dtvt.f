*     This is a recursive subroutine that will compute T*V**T as part of
*     the dorg2r algorithm.
*
*     Parameters
*     N(in):     Number of columns in the matrix Q
*     Q(in/out): Matrix that will hold V and T as described below.
*        V is input only
*        T is input/output

      SUBROUTINE DTVT(N, Q, LDQ)
*
*        Input paramaters
*
         INTEGER           N, LDQ

         DOUBLE PRECISION  Q(LDQ,*)
*
*        Local variables
*
         INTEGER           K

*
*        Local parameters
*
         DOUBLE PRECISION ONE, ZERO
         ONE = 1.0
         ZERO = 0.0
*
*        Beginning of executable statements
*
*
*        Base case
*
         IF (N.EQ.2) THEN
            Q(1,2) = Q(1,1)*Q(2,1) + Q(1,2)
         ELSE IF (N.EQ.3) THEN
*           Manual computation of 
*              T_{1,2} = T_{1,1}V_{1,2} + T_{1,2}V_{2,2}
            Q(1,3) = Q(1,1)*Q(3,1) + Q(1,2)*Q(3,2) + Q(1,3)
            Q(1,2) = Q(1,1)*Q(2,1) + Q(1,2)
            DTVT(N-1, Q(2,2), LDQ)
         ELSE
*
*        Recursive case
*
            K = N / 2
*           Computes T_{1,2} = T_{1,2}V_{2,2}^\top
            DTRMM('Right', 'Lower Triangular', 'Transpose', 'Unit', K,
      $      N - K, ONE, Q(K+1, K+1), LDQ, Q(1,K+1), LDQ)
*           Compute T_{1,2} = T_{1,2} + T_{1,1}V_{1,2}^\top
         END IF

            
      END SUBROUTINE 
