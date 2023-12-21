      SUBROUTINE TESTDTRMMOOP(M, N, LDA, LDT)
         ! Parameters
         INTEGER M,N,LDA,LDT

         ! Scalar variables
         INTEGER           I, J
         DOUBLE PRECISION  NORM_F, TMP

         ! Matrix variables
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), T(:,:), C(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: As(:,:), ATrmm(:,:), Cs(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: Ts(:,:)

         ! External SUBROUTINES
         EXTERNAL DLACPY, DTRMMOOP, DTRMM

         ! Parameters
         DOUBLE PRECISION ONE
         PARAMETER(ONE=1.0D+0)


         ! Allocate our matrices to be the proper sizes
         ! A,As are n by m, T is m by m, C,Cs are m by n
         ! ATrmm is m by n and ATrmm = A**T
         ALLOCATE(A(LDA,M))
         ALLOCATE(As(LDA,M))
         ALLOCATE(ATrmm(M,N))
         ALLOCATE(T(LDT,M))
         ALLOCATE(Ts(LDT,M))
         ALLOCATE(C(M,N))
         ALLOCATE(Cs(M,N))
         ! Fill A, T, C with random elements
         CALL RANDOM_NUMBER(A)
         CALL RANDOM_NUMBER(T)
         CALL RANDOM_NUMBER(C)

         ! Copy A into As
         CALL DLACPY('All', N, M, A, LDA, As, LDA)

         ! Copy C into Cs
         CALL DLACPY('All', M, N, C, M, Cs, M)

         ! Copy T into Ts
         CALL DLACPY('All', M, M, T, LDT, Ts, LDT)

         ! Copy A**T into ATrmm for use with DTRMM
         DO I = 1, M
            DO J = 1, N
               ATrmm(I,J) = A(J,I)
            END DO
         END DO

         CALL DTRMMOOP(M, N, A, LDA, T, LDT, C, M)

         ! Check that A was not touched
         DO I = 1, N
            DO J = 1, M
               IF (A(I,J).NE.As(I,J)) THEN
                  WRITE(*,*) "Inconsistency in A at index i=", I, 
     $                       " j=", J
               END IF
            END DO
         END DO

         ! Check that T was not touched
         DO I = 1, M
            DO J = 1, M
               IF (T(I,J).NE.Ts(I,J)) THEN
                  WRITE(*,*) "Inconsistency in T at index i=", I, 
     $                       " j=", J
               END IF
            END DO
         END DO
         CALL DTRMM('Left', 'Upper', 'No transpose', 'Unit', M, N, ONE, 
     $               T, LDT, ATrmm, M)

         ! Compute ATrmm += Cs
         DO I = 1,M
            DO J = 1,N
               ATrmm(I,J) = ATrmm(I,J) + Cs(I,J)
            END DO
         END DO

         ! Compute ||ATrmm - C||_F
         NORM_F = 0
         DO I = 1, M
            DO J = 1, N
               TMP = ATrmm(I,J) - C(I,J)
               NORM_F = NORM_F + TMP * TMP
            END DO
         END DO

         WRITE(*,*) "||ATrmm - C||_F = ", NORM_F

         ! Deallocate memory
         DEALLOCATE(A)
         DEALLOCATE(As)
         DEALLOCATE(ATrmm)
         DEALLOCATE(T)
         DEALLOCATE(Ts)
         DEALLOCATE(C)
         DEALLOCATE(Cs)
      END SUBROUTINE
