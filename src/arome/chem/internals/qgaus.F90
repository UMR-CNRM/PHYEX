!     ######spl
      subroutine qgaus(func,a,b,ss)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      IMPLICIT NONE
      real a, b, ss, func
      INTEGER :: J
      external func
      real dx,xm,xr,w(5),x(5)
      data w/.2955242247,.2692667193,.2190863625,.1494513491,.066671344/
      data x/.1488743389,.4333953941,.6794095682,.8650633666,.973906529/
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('QGAUS',0,ZHOOK_HANDLE)
      xm = 0.5*(b+a)
      xr = 0.5*(b-a)
      ss = 0.
      do 11 j = 1, 5
            dx = xr*x(j)
            ss = ss + w(j)*(func(xm+dx)+func(xm-dx))
 11   continue
      ss = xr*ss
      IF (LHOOK) CALL DR_HOOK('QGAUS',1,ZHOOK_HANDLE)
      return
      endsubroutine qgaus
