!     ######spl
      SUBROUTINE inter3(ng,xg,yg, n,x,y, FoldIn)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      INTEGER n, ng
      REAL xg(ng)
      REAL x(n), y(n)

      INTEGER FoldIn

! output:
      REAL yg(ng)

! local:
      REAL a1, a2, sum
      REAL tail
      INTEGER jstart, i, j, k
!_______________________________________________________________________

! check whether flag given is legal
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('INTER3',0,ZHOOK_HANDLE)
      IF ((FoldIn .NE. 0) .AND. (FoldIn .NE. 1)) THEN
         WRITE(0,*) '>>> ERROR (inter3) <<<  Value for FOLDIN invalid. '
         WRITE(0,*) '                        Must be 0 or 1'
         STOP
      ENDIF

! do interpolation

      jstart = 1

      DO 30, i = 1, ng - 1

         yg(i) = 0.
         sum = 0.
         j = jstart

         IF (j .LE. n-1) THEN

   20      CONTINUE

             IF (x(j+1) .LT. xg(i)) THEN
                jstart = j
                j = j+1
                IF (j .LE. n-1) GO TO 20
             ENDIF

   25      CONTINUE

             IF ((x(j) .LE. xg(i+1)) .AND. (j .LE. n-1)) THEN

                a1 = MAX(x(j),xg(i))
                a2 = MIN(x(j+1),xg(i+1))

                sum = sum + y(j) * (a2-a1)/(x(j+1)-x(j))
                j = j+1
                GO TO 25

             ENDIF

           yg(i) = sum

         ENDIF

   30 CONTINUE


! if wanted, integrate data "overhang" and fold back into last bin

      IF (FoldIn .EQ. 1) THEN

         j = j-1
         a1 = xg(ng)     ! upper limit of last interpolated bin
         a2 = x(j+1)     ! upper limit of last input bin considered

!        do folding only if grids don't match up and there is more input
         IF ((a2 .GT. a1) .OR. (j+1 .LT. n)) THEN
           tail = y(j) * (a2-a1)/(x(j+1)-x(j))
           DO k = j+1, n-1
              tail = tail + y(k) * (x(k+1)-x(k))
           ENDDO
           yg(ng-1) = yg(ng-1) + tail
         ENDIF

      ENDIF
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('INTER3',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE INTER3
