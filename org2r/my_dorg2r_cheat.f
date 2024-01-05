      SUBROUTINE MY_DORG2R_CHEAT(M, N, Q, LDQ)
         ! Arguments
         INTEGER           M, N, LDQ

         ! Array arguments
         DOUBLE PRECISION  Q(LDQ,*)

         ! External subroutines
         EXTERNAL DTVT_CHEAT, DTRMM

         CALL DTVT_CHEAT(N, Q, LDQ)

         ! Now, we have computed T=TV_1^\top

         ! Compute Q = VT
         ! Q_1 = V_1*T ( LUMM_CHEAT lower-upper matrix mult )
         ! Q_2 = V_2*T (TRMM)

         ! Compute Q_2 first

         CALL DTRMM('Right', 'Upper', 'No-transpose', 'Non-unit', M-N, N, ONE,
     $      Q, LDQ, Q(N+1,1), LDQ)

         ! Compute Q_1

         ! rows, cols, L, LDL, LUnit, U, LDU, UUnit
         CALL LUMM_CHEAT(N, N, Q, LDQ, 'Unit', Q, LDQ, 'Non-unit')

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
