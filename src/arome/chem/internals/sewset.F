!     ######spl
      SUBROUTINE SEWSET (N, ITOL, RTOL, ATOL, YCUR, EWT)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################################################
      REAL RTOL, ATOL, YCUR, EWT
      INTEGER N, ITOL
      DIMENSION RTOL(*), ATOL(*), YCUR(N), EWT(N)
!-----------------------------------------------------------------------
! Call sequence input -- N, ITOL, RTOL, ATOL, YCUR
! Call sequence output -- EWT
! COMMON block variables accessed -- None
!
! Subroutines/functions called by SEWSET.. None
!-----------------------------------------------------------------------
! This subroutine sets the error weight vector EWT according to
!     EWT(i) = RTOL(i)*abs(YCUR(i)) + ATOL(i),  i = 1,...,N,
! with the subscript on RTOL and/or ATOL possibly replaced by 1 above,
! depending on the value of ITOL.
!-----------------------------------------------------------------------
      INTEGER I
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
      IF (LHOOK) CALL DR_HOOK('SEWSET',0,ZHOOK_HANDLE)
      GO TO (10, 20, 30, 40), ITOL
 10   CONTINUE
      DO 15 I = 1, N
 15     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(1)
      IF (LHOOK) CALL DR_HOOK('SEWSET',1,ZHOOK_HANDLE)
      RETURN
 20   CONTINUE
      DO 25 I = 1, N
 25     EWT(I) = RTOL(1)*ABS(YCUR(I)) + ATOL(I)
      IF (LHOOK) CALL DR_HOOK('SEWSET',1,ZHOOK_HANDLE)
      RETURN
 30   CONTINUE
      DO 35 I = 1, N
 35     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(1)
      IF (LHOOK) CALL DR_HOOK('SEWSET',1,ZHOOK_HANDLE)
      RETURN
 40   CONTINUE
      DO 45 I = 1, N
 45     EWT(I) = RTOL(I)*ABS(YCUR(I)) + ATOL(I)
      IF (LHOOK) CALL DR_HOOK('SEWSET',1,ZHOOK_HANDLE)
      RETURN
      END SUBROUTINE SEWSET
