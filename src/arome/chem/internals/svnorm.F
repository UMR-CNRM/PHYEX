!     ######spl
      FUNCTION SVNORM (N, V, W)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##############################
      REAL SVNORM
      REAL V, W
      INTEGER N
      DIMENSION V(N), W(N)
!-----------------------------------------------------------------------
! Call sequence input -- N, V, W
! Call sequence output -- None
! COMMON block variables accessed -- None
!
! Subroutines/functions called by SVNORM.. None
!-----------------------------------------------------------------------
! This function routine computes the weighted root-mean-square norm
! of the vector of length N contained in the array V, with weights
! contained in the array W of length N..
!   SVNORM = sqrt( (1/N) * sum( V(i)*W(i) )**2 )
!-----------------------------------------------------------------------
      REAL SUM
      INTEGER I
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('SVNORM',0,ZHOOK_HANDLE)
      SUM = 0.0E0
      DO 10 I = 1, N
 10     SUM = SUM + (V(I)*W(I))**2
      SVNORM = SQRT(SUM/REAL(N))
      IF (LHOOK) CALL DR_HOOK('SVNORM',1,ZHOOK_HANDLE)
      RETURN
      ENDFUNCTION SVNORM
