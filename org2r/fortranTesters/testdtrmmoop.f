      SUBROUTINE TESTDTRMMOOP(M, N, LDA, LDT)
         ! Parameters
         INTEGER M,N,LDA,LDT

         ! Scalar variables
         INTEGER           I, J, K
         DOUBLE PRECISION  NORM_F, TMP

         ! Matrix variables
         DOUBLE PRECISION, ALLOCATABLE :: A(:,:), T(:,:), C(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: As(:,:), ATrmm(:,:), Cs(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: Ts(:,:), B(:,:), Bs(:,:)

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
         CALL DTRMM('Left', 'Upper', 'No transpose', 'Non-unit', M, N,
     $               ONE, T, LDT, ATrmm, M)

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

         ! Now, we will try testing when we pass in matrices that start in the
         ! same place. The way we do this is that we are going to construct a
         ! matrix B as a k by k matrix where k=max(m,n). Then the leading m by m
         ! principle submatrix will be T, and the leading n by m principle
         ! submatrix will be A. This will ensure that fortran doesn't throw a
         ! fit when we call things this way
         K = MAX(M,N)
         ALLOCATE(B(K,K))
         ALLOCATE(BS(K,K))

         CALL RANDOM_NUMBER(B)
         ! regenerate C
         CALL RANDOM_NUMBER(C)
         ! Copy A into As
         CALL DLACPY('All', N, M, B, K, As, LDA)

         ! Copy B into Bs
         CALL DLACPY('All', K, K , B, K, Bs, K)

         ! Copy C into Cs
         CALL DLACPY('All', M, N, C, M, Cs, M)

         ! Copy T into Ts
         CALL DLACPY('All', M, M, B, LDT, Ts, LDT)

         ! Copy A**T into ATrmm for use with DTRMM
         DO I = 1, M
            DO J = 1, N
               ATrmm(I,J) = As(J,I)
            END DO
         END DO

         ! Call our function
         CALL DTRMMOOP(M, N, B, K, B, K, C, M)

         ! Check that all of B was not touched
         DO I = 1, K
            DO J = 1, K
               IF (B(I,J).NE.Bs(I,J)) THEN
                  WRITE(*,*) "Inconsistency in B at index i=", I, 
     $                       " j=", J
               END IF
            END DO
         END DO
         CALL DTRMM('Left', 'Upper', 'No transpose', 'Non-unit', M, N,
     $               ONE, Ts, LDT, ATrmm, M)

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
         DEALLOCATE(B)
         DEALLOCATE(Bs)
      END SUBROUTINE
