* This is a "cheat" computation of TV_1**T
* We will do this with the construction of 
*
*         |-|
*     Q = |T|
*         |-|
*         |V|
*         |-|
*
*     V = |---|
*         |V_1|
*         |---|
*         |V_2|
*         |---|
*
      SUBROUTINE DTVT_CHEAT(N, Q, LDQ)
*
*        Input paramaters
*
         INTEGER           N, LDQ

         DOUBLE PRECISION  Q(LDQ,*)
*
*        Local variables
*
         INTEGER           K, L
*
*        Local arrays
*
         DOUBLE PRECISION, ALLOCATABLE :: T(:,:), V(:,:)

*
*        External subroutines
*
         EXTERNAL DLACPY, DTRMM, DTRMMOOP
*
*        Local parameters
*
         DOUBLE PRECISION ONE, ZERO
         PARAMETER(ONE=1.0d+0, ZERO=0.0d+0)

         ! Executable statements
         K = N/2
         ! Break down T and V 
         ! We have that 
         ! T_12 = T_11*V_12^\top + T_12*V_22^\top
         ! T_11 = T_11*V_11
         ! T_22 = T_22 * V_22

         ! Compute T_12
         ! T_12 = T_12*V_22
         ! This is a trmm with a lower triangular
         ! matrix
         CALL DTRMM('Right', 'Lower', 'Transpose', 'Unit', K, N - K,
     $      ONE, Q(K+1,K+1), LDQ, Q(1,K+1), LDQ)
         ! Compute T_12 = T_12 + T_11*V_12^\top
         CALL DTRMMOOP(K, N - K, Q(K,1), LDQ, Q(1,1), LDQ, Q(1,K+1), 
     $      LDQ)

         ! Compute T_11 = T_11*V_11^\top
         ! Allocate T and V
         IF (K.NE.0) THEN

            ALLOCATE(T(K,K))
            ALLOCATE(V(K,K))
            ! Copy T_11 into a matrix T
            CALL DLACPY('Upper', K, K, Q, LDQ, T, K)
            ! Copy V_11 into a matrix V
            CALL DLACPY('Lower', K, K, Q, LDQ, V, K)
            ! Compute T = T*V^\top (TRMM with A = V_11)
            CALL DTRMM('Right', 'Lower', 'Transpose', 'Unit',
     $            K, K, ONE, V, K, T, K)
            ! Copy T back into T_11
            CALL DLACPY('Upper', K, K, T, K, Q, LDQ)
            ! Free the memory
            DEALLOCATE(V)
            DEALLOCATE(T)
         END IF

         ! Compute T_22 = T_22*V_22
         ! Allocate T and V
         ALLOCATE(T(N-K,N-K))
         ALLOCATE(V(N-K,N-K))
         ! Copy T_22 into a matrix T
         CALL DLACPY('Upper', N-K, N-K, Q(K+1,K+1), LDQ, T, N-K)
         ! Copy V_22 into a matrix V
         CALL DLACPY('Lower', N-K, N-K, Q(K+1,K+1), LDQ, V, N-K)
         ! Compute T = T*V^\top (TRMM with B = V_22)
         CALL DTRMM('Right', 'Lower', 'Transpose', 'Unit',
     $         N-K, n-K, ONE, V, N-K, T, N-K)
         ! Copy T back into T_22
         CALL DLACPY('Upper', N-K, N-K, T, N-K, Q(K+1,K+1), LDQ)

         ! Free the memory
         DEALLOCATE(V)
         DEALLOCATE(T)
      END SUBROUTINE 
