      RECURSIVE SUBROUTINE LUMM(N, Q, LDQ)
         ! Scalar Arguments
         INTEGER           N, LDQ
         ! Matrix Arguments
         DOUBLE PRECISION  Q(LDQ,*)

         ! Local Variables
         ! Scalars
         INTEGER           K

         ! External Subroutines
         EXTERNAL          DGEMM, DTRMM

         ! Parameters
         DOUBLE PRECISION  ONE
         PARAMETER(ONE=1.0D+0)


         ! Q is of the form
         !
         !     |---------------|
         ! Q = |Q_{1,1} Q_{1,2}|
         !     |Q_{2,1} Q_{2,2}|
         !     |---------------|
         ! k = n/2
         ! Q_{1,1}\in\R^{k\times k}
         ! Q_{1,2}\in\R^{k\times n-k}
         ! Q_{2,1}\in\R^{n-k\times k}
         ! Q_{2,2}\in\R^{n-k\times n-k}
         !
         ! Q is also of the form 
         !     |-|
         ! Q = |U|
         !     |L|
         !     |-|
         ! Due to the breakdown, Q_{1,1} and Q_{2,2} will also have the above
         ! structure
         ! Q_{1,2} is solely a part of U and Q_{2,1} is solely a part of L
         ! U is upper triangular and L is unit lower triangular

         ! U has the form
         ! U = |---------------|
         !     |U_{1,1} U_{1,2}|
         !     |0       U_{2,2}|
         !     |---------------|
         ! L has the form
         ! L = |---------------|
         !     |L_{1,1} 0      |
         !     |L_{2,1} L_{2,2}|
         !     |---------------|
         !
         ! So, we are computing
         ! Q = L*U = 
         !  |--------------------------------------------------|
         !  |L_{1,1}*U_{1,1}  L_{1,1}*U_{1,2}                  |
         !  |L_{2,1}*U_{1,1}  L_{2,1}*U_{1,2} + L_{2,2}*U_{2,2}|
         !  |--------------------------------------------------|
         ! IE
         !  Q_{1,1} = L_{1,1}*U_{1,1} (LUMM)
         !  Q_{1,2} = L_{1,1}*U_{1,2} (TRMM)
         !  Q_{2,1} = L_{2,1}*U_{1,1} (TRMM)
         !  Q_{2,2} = L_{2,1}*U_{1,2} + L_{2,2}*U_{2,2} (LUMM then GEMM)
         !  We compute these from bottom to top

         ! Base case of when N = 1
         IF (N.EQ.1) THEN
            ! This is because we have that L is a 1x1 unit lower triangular
            ! matrix, so Q = L*U = U, which is already stored in Q correctly
            RETURN
         END IF
         K = N / 2

         ! Recursive Case
         ! Compute Q_{2,2} first
         ! Q_{2,2} = L_{2,2}*U_{2,2} (LUMM)
         CALL LUMM(N-K, Q(K + 1, K + 1), LDQ)

         ! Q_{2,2} = L_{2,1}*U_{1,2} + Q_{2,2} (GEMM)
         CALL DGEMM('No transpose', 'No transpose', N-K, N-K, K, 
     $            ONE, Q(K+1,1), LDQ, Q(1,K+1), LDQ, ONE, 
     $            Q(K+1, K+1), LDQ)
         
         ! Compute Q_{2,1}
         ! Q_{2,1} = L_{2,1}*U_{1,1} (TRMM)
         CALL DTRMM('Right', 'Upper', 'No-transpose', 'Non-unit',
     $            N-K, K, ONE, Q, LDQ, Q(K+1,1), LDQ)

         ! Compute Q_{2,1}
         CALL DTRMM('Left', 'Lower', 'No-transpose', 'Unit',
     $            K, N-K, ONE, Q, LDQ, Q(1,K+1), LDQ)
         
         ! Compute Q_{1,1}
         CALL LUMM(K, Q, LDQ)

         ! Done!
      END SUBROUTINE
