      SUBROUTINE LUMM_CHEAT(N, L, LDL, LDIAG, U, LDU, UDIAG)

         !Arguments
         INTEGER           M, N, LDL, LDU

         DOUBLE PRECISION  L(LDL,*), U(LDU,*)

         CHARACTER         LDIAG, UDIAG

         ! Scalar variables
         INTEGER           I, J
         ! Matrix variables
         DOUBLE PRECISION, ALLOCATABLE :: Ls(:,:), Us(:,:)

         ! Parameters
         DOUBLE PRECISION ONE
         PARAMETER(ONE=1.0D+0)
         ! We overwrite L as we need to choose one.
         ALLOCATE(Ls(N, N))
         ALLOCATE(Us(N, N))

         CALL DLACPY('Lower', N, N, L, LDL, Ls, N)

         IF(LDIAG.EQ.'U'.OR.LDIAG.EQ.'u') THEN
            DO I = 1, N
               Ls(I,I) = 1.0
            END DO
         END IF
         DO I = 1, N-1
            DO J = I+1, N
               Ls(I,J) = 0.0
            END DO
         END DO

         CALL DLACPY('Upper', N, N, U, LDU, Us, N)
         IF(UDIAG.EQ.'U'.OR.UDIAG.EQ.'u') THEN
            DO I = 1, N
               Us(I,I) = 1.0
            END DO
         END IF
         DO I = 2, N
            DO J = 1, I - 1
               Us(I,J) = 0.0
            END DO
         END DO
         ! A = L, B = U
         CALL DTRMM('Left', 'Lower', 'No-transpose', LDIAG, N, N, ONE,
     $      Ls, N, Us, N)

         ! Copy Ls back into L

         CALL DLACPY('All',N, N, Ls, N, L, LDL)

         ! Free memory
         DEALLOCATE(Ls)
         DEALLOCATE(Us)

      END SUBROUTINE
