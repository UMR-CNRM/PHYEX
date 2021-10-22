!     ######spl
      subroutine zero2(x,m,n)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      implicit none
! m,n : dimensions of x, exactly as specified in the calling program
      integer i, j, m, n
      real x(m,n)
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ZERO2',0,ZHOOK_HANDLE)
      do 1 j = 1, n
         do 2 i = 1, m
            x(i,j) = 0.
 2       continue
 1    continue
      IF (LHOOK) CALL DR_HOOK('ZERO2',1,ZHOOK_HANDLE)
      return
      endsubroutine zero2
