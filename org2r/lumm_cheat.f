      SUBROUTINE LUMM_CHEAT(N, L, LDL, LDIAG, U, LDU, UDIAG)

         !Arguments
         INTEGER           M, N, LDL, LDU

         DOUBLE PRECISION  L(LDL,*), U(LDU,*)

         CHARACTER         LDIAG, UDIAG

         ! Scalar variables
         INTEGER           I, J
         ! Matrix variables
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:)

         ! Parameters
         DOUBLE PRECISION ONE
         PARAMETER(ONE=1.0D+0)

         ALLOCATE(A(N,N))

         DO I = 1, N
            DO J = 1, N
               A(I,J) = 0
            END DO
         END DO

         DO 10 I = N, 1, -1
            CALL DGER(N - I + 1, N - I + 1, ONE, L(I,I), 1, U(I,I), LDU, 
     $         A(I,I), N)
   10    CONTINUE

         CALL DLACPY('All', N, N, A, N, L, LDL)

         DEALLOCATE(A)

      END SUBROUTINE
