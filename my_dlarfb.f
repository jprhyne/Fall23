*> \brief \b DLARFB applies a block reflector or its transpose to a general rectangular matrix.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DLARFB + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dlarfb.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dlarfb.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dlarfb.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
*                          T, LDT, C, LDC, WORK, LDWORK )
*
*       .. Scalar Arguments ..
*       CHARACTER          DIRECT, SIDE, STOREV, TRANS
*       INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*       ..
*       .. Array Arguments ..
*       DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
*      $                   WORK( LDWORK, * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DLARFB applies a real block reflector H or its transpose H**T to a
*> real m by n matrix C, from either the left or the right.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*>
*> \param[in] M
*> \verbatim
*>          M is INTEGER
*>          The number of rows of the matrix C.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of columns of the matrix C.
*> \endverbatim
*>
*> \param[in] K
*> \verbatim
*>          K is INTEGER
*>          The order of the matrix T (= the number of elementary
*>          reflectors whose product defines the block reflector).
*>          M >= K >= 0;
*> \endverbatim
*>
*> \param[in] V
*> \verbatim
*>          V is DOUBLE PRECISION array, dimension
*>                                (LDV,K) if STOREV = 'C'
*>          The matrix V. See Further Details.
*> \endverbatim
*>
*> \param[in] LDV
*> \verbatim
*>          LDV is INTEGER
*>          The leading dimension of the array V.
*>          LDV >= max(1,M);
*> \endverbatim
*>
*> \param[in] T
*> \verbatim
*>          T is DOUBLE PRECISION array, dimension (LDT,K)
*>          The triangular k by k matrix T in the representation of the
*>          block reflector.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= K.
*> \endverbatim
*>
*> \param[in,out] C
*> \verbatim
*>          C is DOUBLE PRECISION array, dimension (LDC,N)
*>          On entry, the m by n matrix C.
*>          On exit, C is overwritten by H*C or H**T*C or C*H or C*H**T.
*> \endverbatim
*>
*> \param[in] LDC
*> \verbatim
*>          LDC is INTEGER
*>          The leading dimension of the array C. LDC >= max(1,M).
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is DOUBLE PRECISION array, dimension (LDWORK,K)
*> \endverbatim
*>
*> \param[in] LDWORK
*> \verbatim
*>          LDWORK is INTEGER
*>          The leading dimension of the array WORK.
*>          LDWORK >= max(1,N);
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
*> \ingroup doubleOTHERauxiliary
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The shape of the matrix V and the storage of the vectors which define
*>  the H(i) is best illustrated by the following example with n = 5 and
*>  k = 3. The elements equal to 1 are not stored; the corresponding
*>  array elements are modified but restored on exit. The rest of the
*>  array is not used.
*>
*>               V = (  1       )
*>                   ( v1  1    )
*>                   ( v1 v2  1 )
*>                   ( v1 v2 v3 )
*>                   ( v1 v2 v3 )
*>
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE MY_DLARFB( M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DGEMM, DTRMM
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
*
*              Form  H * C  or  H**T * C  where  C = ( C1 = 0)
*                                                    ( C2 = I)
*
*              W := C**T * V  =  (C1**T * V1 + C2**T * V2)  (stored in WORK)
*                             =  (V2)
*
*              W := C1**T
*
*              Copies C1^T into Work      
*              Note: C1 is all 0s so we probably don't need to do this
*------------------------------------------------------------------
* Removed below
*------------------------------------------------------------------
*               DO 10 J = 1, K
*                  CALL DCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
*   10          CONTINUE
*
*              W := W * V1. Note: W is 0
*
*               CALL DTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
*     $                     K, ONE, V, LDV, WORK, LDWORK )
*               IF( M.GT.K ) THEN
*
**                 W := C2**T * V2
*                 W := V2
*
*                CALL DGEMM( 'Transpose', 'No transpose', N, K, M-K,
*     $                      ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV,
*     $                      ZERO, WORK, LDWORK )
*------------------------------------------------------------------
* Removed above
*------------------------------------------------------------------
*------------------------------------------------------------------
* Added below
*------------------------------------------------------------------
*                 CALL DLACPY('All', N, K, V(K+1,1), LDV,WORK,LDWORK)
*------------------------------------------------------------------
* Added above
*------------------------------------------------------------------
*------------------------------------------------------------------
* Removed below
*------------------------------------------------------------------
*               END IF
*------------------------------------------------------------------
* Removed above
*------------------------------------------------------------------
*
*              W := W * T**T  or  W * T
*
*               CALL DTRMM( 'Right', 'Upper', 'Transpose', 'Non-unit', N,
*     $                     K,ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W**T
*
*------------------------------------------------------------------
* Removed below
*------------------------------------------------------------------
*               IF( M.GT.K ) THEN
*------------------------------------------------------------------
* Removed above
*------------------------------------------------------------------
*
*                 C2 := C2 - V2 * W**T
*
*                CALL DGEMM( 'No transpose', 'Transpose', M-K, N, K,
*     $                      -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE,
*     $                      C( K+1, 1 ), LDC )
*               END IF
*
*              W := W * V1**T
*
*               CALL DTRMM( 'Right', 'Lower', 'Transpose', 'Unit', N, K,
*     $                     ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W**T
*
*
*              C1 := -W**T
*
*------------------------------------------------------------------
* Modified below
*------------------------------------------------------------------
*               DO 30 J = 1, K
*                  DO 20 I = 1, N
*                     C( J, I ) = -WORK( I, J )
*   20             CONTINUE
*   30          CONTINUE
*------------------------------------------------------------------
* Modified above
*------------------------------------------------------------------
      RETURN
*
*     End of DLARFB
*
      END
