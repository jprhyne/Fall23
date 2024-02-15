      RECURSIVE SUBROUTINE MY_DLARFT_UT(M, N, V, LDV, TAU, T, LDT)
         ! Arguments
         ! Scalars
         INTEGER           M, N, LDV, LDT
         ! Matrix 
         DOUBLE PRECISION  V(LDV,*), T(LDT,*), TAU(N)
         ! Local variables
         INTEGER           I,J,K,INFO
         ! Parameters
         DOUBLE PRECISION ONE, NEG_ONE, ZERO, HALF
         PARAMETER(ONE=1.0D+0, HALF=0.5D+0, ZERO = 0.0)
         ! Implementation of the algorithm listed in the following paper
         ! https://www.cs.utexas.edu/users/flame/pubs/p169-joffrain.pdf
         ! Compute T = V^\top V
         V(1,1) = ONE
         DO J = 2, N
            DO I = 1, J - 1
               V(I,J) = ZERO
            END DO
            V(J,J) = ONE
         END DO
               
         CALL DSYRK('Upper', 'Transpose', N, M, ONE, V, LDV, ZERO, 
     $         T, LDT)
c         CALL MY_DST3RK(M, N, V, LDV, ZERO, T, LDT)
         ! Scales the diagonal by 1/2
         CALL DSCAL(N, HALF, T(1,1), LDT + 1)
         ! Replaces T with T^{-1}
         CALL DTRTRI('Upper', 'Non unit', N, T, LDT, INFO)
      END SUBROUTINE
