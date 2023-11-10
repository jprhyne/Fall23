*TODO: Rewrite this to be actually standalone. IE change the references
*to A(I,I) to be something like A(1,1)      
      SUBROUTINE COMPUTEQ2(M, N, K, A, LDA, T, LDT, TAU)
*     .. Scalar Arguments ..
      INTEGER           M, N, K, LDA, LDT

*     .. Array Arguments ..
      DOUBLE PRECISION  A(LDA, *), T(LDT,*), TAU(*)

*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )

*     .. Local Scalars ..
      INTEGER           I, II, IB, NB, KI, KK, JJ, IINFO

*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT

*     .. External Functions ..
      INTEGER           ILAENV
      EXTERNAL          ILAENV

*     Same blocking value as is done from DORGQR
      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )

*     Compute KK, KI like is done in DORGQR
      KI = K - 2 * NB
      KK = K - NB
      I  = KK + 1
      IB = NB
*
*     Form the triangular factor of the block reflector
*     H = H(i) H(i+1) . . . H(i+ib-1)
*
      CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $             A( I, I ), LDA, TAU( I ), T, LDT )
*
*     Apply H to A(i:m,i+ib:n) from the left
*
*
*     C1 := V2**T
*
      DO 36 JJ = I, K
         DO 26 II = K + 1, N
            A( JJ, II ) = A( II, JJ )
   26    CONTINUE
   36 CONTINUE
*
*     C1 := T * C1
*
      CALL DTRMM( 'Left', 'Upper', 'No transpose', 'Non-unit', IB,
     $            N-K,ONE, T, LDT, A(I,I+IB),LDA )
*
*     C2 := C2 - V2 * C1
*
      CALL DGEMM( 'No transpose', 'No transpose', M-IB-KK, N-K, IB,
     $            -ONE, A( I+IB, I ), LDA, A(I,I+IB),LDA, ZERO,
     $            A( I+IB, I+IB ), LDA )
      DO 14 JJ = K + 1, N
         A(JJ,JJ) = 1 + A(JJ,JJ)
   14 CONTINUE
*
*     C1 := -V1 * C1 
*
      CALL DTRMM( 'Left', 'Lower', 'No transpose', 'Unit', IB, N-K,
     $            -ONE, A(I,I), LDA, A(I,I+IB),LDA )
*
*     Apply H to rows i:m of current block
*
      CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), T,
     $             IINFO )

*     First, we compute the part of Q (Q2) that is to the right of K.
*     If A was of rank K, then this is an orthogonal basis of the
*     nullspace of A. Otherwise, this is a subspace of the nullspace of
*     A
      DO 55 I = KI + 1, 1, -NB
         IB = NB
*
*         Form the triangular factor of the block reflector
*        H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                A( I, I ), LDA, TAU( I ), T, LDT )
*
*        Apply H to A(i:m,k+1:n) from the left
*
*        C12 = V2**T * C22
*
         CALL DGEMM( 'Transpose', 'No transpose', IB, N - K,
     $               M-I+1-IB, ONE, A(I+IB,I), LDA, A(I+IB,K+1),
     $               LDA, ZERO, A(I,K+1), LDA)

*
*        C12 := T * C12
*
         CALL DTRMM('Left', 'Upper', 'No transpose', 'Non-unit',
     $              IB, N- K, ONE, T, LDT, A(I,K + 1), LDA)
*
*        C22 := C22 - V2 * C12
*
         CALL DGEMM( 'No transpose', 'No transpose', M-I-IB+1,
     $               N-K, IB, -ONE, A(I+IB,I), LDA, A(I,K+1),
     $               LDA,ONE, A(I+IB, K+1), LDA)
*
*        C12 := -V1 * C12
*
         CALL DTRMM( 'Left', 'Lower', 'No-transpose', 'Unit', IB, 
     $               N-K, -ONE, A(I,I), LDA, A(I,K+1), LDA)
   55 CONTINUE
*     This checks for if K was a perfect multiple of NB
*     so that we only have a special case for the last block when
*     necessary
      IF(I.LT.1) THEN
         IB = I + NB - 1
         I = 1
*
*        Form the triangular factor of the block reflector
*        H = H(i) H(i+1) . . . H(i+ib-1)
*
         CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                A( I, I ), LDA, TAU( I ), T, LDT )
*
*        C12 = V2**T * C22
*
         CALL DGEMM( 'Transpose', 'No transpose', IB, N - K,
     $               M-I+1-IB, ONE, A(I+IB,I), LDA, A(I+IB,K+1),
     $               LDA, ZERO, A(I,K+1), LDA)
*
*        C12 := T * C12
*
         CALL DTRMM('Left', 'Upper', 'No transpose', 'Non-unit',
     $              IB, N- K, ONE, T, LDT, A(I,K + 1), LDA)
*
*        C22 := C22 - V2 * C12
*
         CALL DGEMM( 'No transpose', 'No transpose', M-I-IB+1,
     $               N-K, IB, -ONE, A(I+IB,I), LDA, A(I,K+1),
     $               LDA,ONE, A(I+IB, K+1), LDA)
*
*        C12 := -V1 * C12
*
         CALL DTRMM( 'Left', 'Lower', 'No-transpose', 'Unit', IB, 
     $                  N-K, -ONE, A(I,I), LDA, A(I,K+1), LDA)
      END IF
      END
