      SUBROUTINE TEST_DORG2R(M, N, K, LDA, LDQ)
         ! Arguments
         INTEGER  M, N, K, LDA, LDQ

         ! Local variables
         DOUBLE PRECISION  NORMA, NORM_ORTH, NORM_REPRES
         INTEGER           LWORK, I, J, INFO
         ! Local arrays
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), Q(:,:), As(:,:)

         DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: TAU

         DOUBLE PRECISION, ALLOCATABLE :: WORKMAT(:,:), WORK(:), T(:,:)

         ! External Subroutines
         EXTERNAL DLACPY, DGEQRF, DLARFT

         ! External Functions
         DOUBLE PRECISION, EXTERNAL :: DLANGE

         ! Parameters
         DOUBLE PRECISION ONE, ZERO
         INTEGER           NEG_ONE
         PARAMETER(ONE=1.0D+0, ZERO=0, NEG_ONE=-1.0D+0)



         ALLOCATE(A(LDA,K))
         ALLOCATE(As(LDA,K))
         ALLOCATE(Q(LDQ,N))
         ALLOCATE(WORK(1))
         ALLOCATE(TAU(K))
         ALLOCATE(T(N,N))

         CALL DLACPY('All', M, K, A, LDA, As, LDA)
         NORMA = DLANGE('Frobenius', M, K, A, LDA, WORK)

         CALL DGEQRF(M, K, A, LDA, TAU, WORK, NEG_ONE, INFO)
         LWORK = WORK(1)
         DEALLOCATE(WORK)
         ALLOCATE(WORK(LWORK))

         CALL DGEQRF(M, K, A, LDA, TAU, WORK, LWORK, INFO)
         ! Copy into Q
         CALL DLACPY('All', M, K, A, LDA, Q, LDQ)

         DEALLOCATE(WORK)
         ! Compute the triangular factor T
         CALL DLARFT('Forward', 'Column', K, N, Q, LDQ, TAU, T, N)

         ! Copy T into where R was inside Q
         CALL DLACPY('Upper', N, N, T, N, A, LDA)

         ! Now call MY_DORG2R
         CALL MY_DORG2R_CHEAT(M, N, Q, LDQ)
         
         ALLOCATE(WORKMAT(N,N))
         CALL DLASET('A', N, N, ZERO, ONE, WORKMAT, N)
         CALL DSYRK('Upper', 'Transpose', N, M, ONE, Q, LDQ, NEG_ONE,
     $         WORKMAT, N)

         NORM_ORTH = DLANGE('Frobenius',N,N, WORKMAT, N, WORK)

         DEALLOCATE(WORKMAT)
         ALLOCATE(WORKMAT(M,K))
         CALL DLACPY('All', M, K, Q, LDQ, WORKMAT, M)
         CALL DTRMM('Right', 'Upper', 'No-transpose', 'non-unit',
     $      M, K, ONE, A, LDA, WORKMAT, M)

         DO I = 1, M
            DO J = 1, K
               WORKMAT(I,J) = WORKMAT(I,J) - As(I,J)
            END DO
         END DO
         NORM_REPRES = DLANGE('Frobenius', M, K, WORKMAT,
     $      M, WORK)
         NORM_REPRES = NORM_REPRES / NORMA

         DEALLOCATE(WORKMAT)

         WRITE(*,*) "representation norm: ", NORM_REPRES
         WRITE(*,*) "orthogonal norm:     ", NORM_ORTH

         DEALLOCATE(A)
         DEALLOCATE(As)
         DEALLOCATE(Q)
         DEALLOCATE(TAU)
         DEALLOCATE(T)

      END SUBROUTINE
