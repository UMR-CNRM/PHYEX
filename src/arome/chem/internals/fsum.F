!     ######spl
      FUNCTION fsum(n,x)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      INTEGER n
      REAL x(n)

! function value:
      REAL fsum

! local:
      INTEGER i
!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FSUM',0,ZHOOK_HANDLE)
      fsum = 0.
      DO 10, i = 1, n
         fsum=fsum+x(i)
   10 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('FSUM',1,ZHOOK_HANDLE)
      RETURN
      ENDFUNCTION FSUM
