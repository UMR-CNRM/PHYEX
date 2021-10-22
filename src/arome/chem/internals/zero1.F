!     ######spl
      subroutine zero1(x,m)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      implicit none
      integer i, m
      real x(m)
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ZERO1',0,ZHOOK_HANDLE)
      do 1 i = 1, m
         x(i) = 0.
 1    continue
      IF (LHOOK) CALL DR_HOOK('ZERO1',1,ZHOOK_HANDLE)
      return
      endsubroutine zero1
