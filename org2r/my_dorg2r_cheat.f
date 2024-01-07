      SUBROUTINE MY_DORG2R_CHEAT(M, N, Q, LDQ)
         ! Arguments
         INTEGER           M, N, LDQ

         ! Array arguments
         DOUBLE PRECISION  Q(LDQ,*)

         ! Scalar variables
         INTEGER           I

         ! Matrix variables
         DOUBLE PRECISION, ALLOCATABLE :: T(:,:), V(:,:)

         ! External subroutines
         EXTERNAL DTVT_CHEAT, DTRMM

         ! Parameters
         DOUBLE PRECISION ONE
         PARAMETER(ONE=1.0D+0)
         

         CALL DTVT_CHEAT(N, Q, LDQ)

         ! Now, we have computed T=TV_1^\top

         ! Compute Q = VT
         ! Q_1 = V_1*T ( LUMM_CHEAT lower-upper matrix mult )
         ! Q_2 = V_2*T (TRMM)
         
         ! Copy the 'V' and 'T' into a matrix for use with DTRMM
         ! First allocate the memory
         ALLOCATE(T(N,N))
         ALLOCATE(V(M-N,N))
         ! Copy into T
         CALL DLACPY('Upper', N, N, Q, LDQ, T, N) 
         ! Copy into V
         CALL DLACPY('All', M-N, N, Q(N+1,1), LDQ, V, M-N)

         ! Compute Q_2
         !CALL DTRMM('Right', 'Upper', 'No-transpose', 'Non-unit', M-N, N, ONE,
*     $      T, N, V, M-N)
         CALL DTRMM('R', 'U', 'N', 'N', M-N, N, ONE,
     $      T, N, V, M-N)
         ! Copy V back into Q_2
         CALL DLACPY('All', M-N, N, V, M-N, Q(N+1,1), LDQ)

         ! Reallocate memory for V to be used in computing Q_1
         DEALLOCATE(V)
         ALLOCATE(V(N, N))
         ! Copy lower part of Q_1 into V
         CALL DLACPY('Lower', N, N, Q, LDQ, V, N)
         ! Set the diagonal of V to be 1
         V(1,1) = 1
         DO I = 2, N
            V(I,I) = 1
         END DO
         ! Compute Q_2 first
         !CALL DTRMM('Right', 'Upper', 'No-transpose', 'Non-unit', M-N, N, ONE,
*     $      Q, LDQ, Q(N+1,1), LDQ)

         ! Compute Q_1

         CALL LUMM_CHEAT(N, V, N, 'Unit', T, N, 'Non-unit')

         ! rows, cols, L, LDL, LUnit, U, LDU, UUnit
         !CALL LUMM_CHEAT(N, N, Q, LDQ, 'Unit', Q, LDQ, 'Non-unit')

         ! Compute "I" - Q
         ! M must be greater than or equal to N. (restricting to the "tall and
         ! skinny" case
         DO I = 1, N
            Q(I,I) = 1 - Q(I,I)
         END DO

         ! Now, we should have Q that satisfies the conditions of org2r. IE Q
         ! has orthonormal columns and is a Q associated with A in the QR
         ! decomposition

      END SUBROUTINE
