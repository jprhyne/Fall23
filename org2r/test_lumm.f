      SUBROUTINE TEST_LUMM(N, LDL, LDU)
         INTEGER N, LDL, LDU

         ! Local variables
         DOUBLE PRECISION  NORMF, TMP

         INTEGER           I, J

         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: L(:,:), U(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: Ls(:,:), Us(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:)

         ! External subroutines
         EXTERNAL LUMM_CHEAT, DLACPY, DTRMM

         ! Parameters
         DOUBLE PRECISION  ONE, ZERO
         PARAMETER(ONE=1.0D+0, ZERO=0.0D+0)

         ! Allocate memory
         ALLOCATE(L(LDL,N))
         ALLOCATE(U(LDU,N))
         ALLOCATE(Ls(N,N))
         ALLOCATE(Us(N,N))
         ALLOCATE(A(N,N))

         ! Generate L and U as random matrices
         CALL RANDOM_NUMBER(L)
         CALL RANDOM_NUMBER(U)

         DO I = 1, N
            DO J = 1, N
               IF(I.LT.J) THEN
                  L(I,J) = 0
               ELSE IF (I.GT.J) THEN
                  U(I,J) = 0
               END IF
            END DO
         END DO

         CALL DLACPY('All', N, N, L, LDL, Ls, N)
         CALL DLACPY('All', N, N, U, LDU, Us, N)


         CALL LUMM_CHEAT(N, L, LDL, 'Non-unit', U, LDU, 'Non-unit')

         ! Check that U was not touched
         DO I = 1, N
            DO J = 1, N
               IF (U(I,J).NE.Us(I,J)) THEN
                  WRITE(*,*) "Inconsistency at index i = ", i, "j = ", j
               END IF
            END DO
         END DO

         ! Using Ls and Us compute Ls * Us using TRMM
         ! A = Ls, B = Us
         ! First, since we are treating Us as the 'non-triangular' matrix in 
         ! dtrmm, we need to first enforce it being upper triangular.
         DO I = 2, N
            DO J = 1, I - 1
               Us(I,J) = 0
            END DO
         END DO
         CALL DTRMM('Left', 'Lower', 'No', 'Non-unit', N, N, ONE,
     $      Ls, N, Us, N)

         ! Compute Us - Ls and find the norm thereof
         NORMF=0
         DO I = 1, N
            DO J = 1, N
               TMP = Us(I,J) - L(I,J)
               NORMF = NORMF + TMP * TMP
            END DO
         END DO

         WRITE(*,*) "||Us - L||_F = ", NORMF

         DEALLOCATE(L)
         DEALLOCATE(Ls)
         DEALLOCATE(U)
         DEALLOCATE(Us)
         DEALLOCATE(A)
      END SUBROUTINE
