*> \brief \b DORGQR
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DORGQR + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dorgqr.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dorgqr.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dorgqr.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DORGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*       .. Scalar Arguments ..
*       INTEGER            INFO, K, LDA, LWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DORGQR generates an M-by-N real matrix Q with orthonormal columns,
*> which is defined as the first N columns of a product of K elementary
*> reflectors of order M
*>
*>       Q  =  H(1) H(2) . . . H(k)
*>
*> as returned by DGEQRF.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix Q. M >= 0.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix Q. M >= N >= 0.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The number of elementary reflectors whose product defines the
*>          matrix Q. N >= K >= 0.
*> \endverbatim
*>
*> \param[in,out] A
*> \verbatim
*>          A is DOUBLE PRECISION array, dimension (LDA,N)
*>          On entry, the i-th column must contain the vector which
*>          defines the elementary reflector H(i), for i = 1,2,...,k, as
*>          returned by DGEQRF in the first k columns of its array
*>          argument A.
*>          On exit, the M-by-N matrix Q.
*> \endverbatim
*>
*> \param[in] LDA
*> \verbatim
*>          LDA is INTEGER
*>          The first dimension of the array A. LDA >= max(1,M).
*> \endverbatim
*>
*> \param[in] TAU
*> \verbatim
*>          TAU is DOUBLE PRECISION array, dimension (K)
*>          TAU(i) must contain the scalar factor of the elementary
*>          reflector H(i), as returned by DGEQRF.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of the array WORK. LWORK >= max(1,N).
*>          For optimum performance LWORK >= N*NB, where NB is the
*>          optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument has an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup doubleOTHERcomputational
*
*  =====================================================================
      SUBROUTINE MY_DORGQR( M, N, K, NB, A, LDA, TAU, WORK, LWORK, INFO)
      IMPLICIT NONE
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARFB, DLARFT, DORG2R, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV

*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*      NB = ILAENV( 1, 'DORGQR', ' ', M, N, K, -1 )
*
*     Debugging purposes
*
      
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DORGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'DORGQR', ' ', M, N, K, -1 ) )
         NX = 0
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'DORGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.EQ.0 ) THEN
        CALL DORG2R( M, N, K, A( 1, 1 ), LDA,
     $                TAU( 1 ), WORK, IINFO )
      END IF
      IF( KK.GT.0 ) THEN
*       First set our matrix to be of the form
*       —————
*       |A|0|
*       —————
*       |A|I|
*       —————
        DO 20 J = KK + 1, N
          DO 10 I = 1, M
            A( I, J ) = ZERO
   10     CONTINUE
          A( J, J ) = ONE
   20   CONTINUE
*       Form the triangular factor of the block reflector
*       H = H(kk+1) H(kk+2) . . . H(k)
*       Note: May need to move to start at A(KK,KK)
        IB = K - KK
        CALL DLARFT( 'Forward', 'Columnwise', M - KK, IB,
     $               A(KK + 1, KK + 1), LDA, TAU(KK + 1), WORK,
     $               LDWORK )

*       Apply H to A(kk+1:m, k:n)
*        CALL DLARFB( 'Left', 'No transpose', 'Forward', 'Columnwise', 
*     $               M - KK, N - K, IB, A( KK+1, KK+1), LDA, WORK,
*     $               LDWORK, A( KK, K), LDA, WORK(IB + 1), LDWORK)
*        SIDE = 'L'
*        TRANS = 'N', so TRANST = 'T'
*        DIRECT = 'F'
*        STOREV = 'C'
*        MT = M - KK
*        NT = N - K
*        K = IB
*        V = A(KK+1,KK+1)
*        note: V(i,j) = A(KK + 1 + i, KK + 1 + j)
*        LDV = LDA
*        T = WORK
*        note: T(i,j) = WORK(i + j * LDWORK)
*        LDT = LDWORK
*        C = A(KK,K)
*        note: C(i,j) = A(KK + i, K + j)
*        LDC = LDA
*        WORK = WORK(IB + 1)
*        note: WORK(i,j) = WORK(IB + 1 + i + j * ldwork)
*        LDWORKT = LDWORK
*
*              Form  H * C  or  H**T * C  where  C = ( C1 )
*                                                    ( C2 )
*
*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
*
*              W := C1**T
*
        DO 15 J = 1, IB
            CALL DCOPY( N - K, A( KK + J, K + 1 ), LDA, 
     $                  WORK(IB + 1 + 1 + J * LDWORK ), 1 )
   15   CONTINUE
*
*              W := W * V1
*
        CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', 
     $              N - KK, IB, ONE, A(KK + 1, KK + 1), LDA, 
     $              WORK(IB + 1), LDWORK )
        IF( M - KK.GT.IB ) THEN
*
*                 W := W + C2**T * V2
*
           CALL DGEMM( 'Transpose', 'No transpose', N - K, IB, 
     $                 M - KK - IB,
     $                 ONE, A( KK + K + 1, K + 1 ), LDA, 
     $                 A( KK + 1 + K + 1, KK + 1 + 1 ), LDA,
     $                 ONE, WORK(IB + 1), LDWORK )
        END IF
*
*              W := W * T**T  or  W * T
*
        CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-unit', 
     $              N - K, IB, ONE, WORK, LDWORK, WORK(IB + 1), LDWORK )
*
*              C := C - V * W**T
*
        IF( M - KK.GT.IB ) THEN
*
*                 C2 := C2 - V2 * W**T
*
           CALL DGEMM( 'No transpose', 'Transpose', M - KK - IB, N - K,
     $                 IB, -ONE, A( KK + 1 + IB + 1, KK + 1 + 1 ), 
     $                 LDA, WORK(IB + 1), LDWORK, ONE,
     $                 A( KK + IB + 1, K + 1 ), LDA )
        END IF
*
*              W := W * V1**T
*
        CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N - K, IB,
     $              ONE, A(KK + 1, KK + 1), LDA, WORK(IB + 1), LDWORK )
*
*              C1 := C1 - W**T
*
        DO 35 J = 1, IB
           DO 25 I = 1, N - K
             A(KK + J, K + I) = A( KK + J, K + I ) - WORK(IB + 1 + I + J
     $                                                    * LDWORK)
   25      CONTINUE
   35   CONTINUE

*       Apply H to rows kk+1:m of the current block        
        CALL DORG2R( M - KK, IB, IB, A(KK+1,KK+1), LDA, TAU(KK + 1),
     $               WORK, IINFO )

         DO 50 I = KI+1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL DLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL DLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL DORG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of DORGQR
*
      END
