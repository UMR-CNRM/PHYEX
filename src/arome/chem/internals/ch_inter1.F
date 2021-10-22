!     ######spl
      SUBROUTINE ch_inter1(ng,xg,yg, n,x,y)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)

! output:
      REAL yg(ng)

! local:
      REAL slope
      INTEGER jsave, i, j
!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CH_INTER1',0,ZHOOK_HANDLE)
      jsave = 1
      DO 20, i = 1, ng
         yg(i) = 0.
         j = jsave
   10    CONTINUE
            IF ((x(j) .GT. xg(i)) .OR. (xg(i) .GE. x(j+1))) THEN
               j = j+1
               IF (j .LE. n-1) GOTO 10
!        ---- end of loop 10 ----
            ELSE
               slope = (y(j+1)-y(j)) / (x(j+1)-x(j))
               yg(i) = y(j) + slope * (xg(i) - x(j))
               jsave = j
             ENDIF
   20 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('CH_INTER1',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE CH_INTER1
