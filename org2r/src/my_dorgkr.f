      SUBROUTINE MY_DORGKR(M, N, Q, LDQ)
         ! Arguments
         INTEGER           M, N, LDQ

         ! Array arguments
         DOUBLE PRECISION  Q(LDQ,*)

         ! Scalar variables
         INTEGER           I, J

         DOUBLE PRECISION, ALLOCATABLE :: V(:,:)
         ! External subroutines
         EXTERNAL DTRMM, DTVT, LUMM

         ! Parameters
         DOUBLE PRECISION ONE
         PARAMETER(ONE=1.0D+0)
         

         CALL DTVT(N, Q, LDQ)

         ! Now, we have computed T=TV_1^\top

         ! Compute Q = VT
         ! first n rows
         ! Q_1 = V_1*T ( LUMM lower-upper matrix mult )
         ! Q_2 = V_2*T (TRMM)
         
         IF (M.GT.N) THEN
            CALL DTRMM('Right', 'Upper', 'No-transpose', 'Non-unit',
     $         M-N, N, ONE, Q, LDQ, Q(N+1,1), LDQ)
         END IF

         CALL LUMM(N, Q, LDQ)

         ! Compute "I" - Q
         ! M must be greater than or equal to N. (restricting to the "tall and
         ! skinny" case)

         ! Convert to this
         ! Q = -Q
         ! Q += "I"
         DO J = 1, N
            DO I = 1, J-1
               Q(I,J) = -Q(I,J)
            END DO
            Q(J,J) = 1 - Q(J,J)
            DO I = J+1, M
               Q(I,J) = -Q(I,J)
            END DO
         END DO
         ! Now, we should have Q that satisfies the conditions of org2r. IE Q
         ! has orthonormal columns and is a Q associated with A in the QR
         ! decomposition

      END SUBROUTINE
