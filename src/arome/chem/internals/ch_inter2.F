!     ######spl
      SUBROUTINE ch_inter2(ng,xg,yg,n,x,y,ierr)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      INTEGER ng, n
      REAL x(n), y(n), xg(ng)

! output:
      REAL yg(ng)

! local:
      REAL area, xgl, xgu
      REAL darea, slope
      REAL a1, a2, b1, b2
      INTEGER ngintv
      INTEGER i, k, jstart
      INTEGER ierr
!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CH_INTER2',0,ZHOOK_HANDLE)
      ierr = 0

!  test for correct ordering of data, by increasing value of x

      DO 10, i = 2, n
         IF (x(i) .LE. x(i-1)) THEN
            ierr = 1
            WRITE(*,*)'data not sorted'
            IF (LHOOK) CALL DR_HOOK('CH_INTER2',1,ZHOOK_HANDLE)
            RETURN
         ENDIF
   10 CONTINUE

      DO i = 2, ng
        IF (xg(i) .LE. xg(i-1)) THEN
           ierr = 2
          WRITE(0,*) '>>> ERROR (inter2) <<<  xg-grid not sorted!'
          IF (LHOOK) CALL DR_HOOK('CH_INTER2',1,ZHOOK_HANDLE)
          RETURN
        ENDIF
      ENDDO

! check for xg-values outside the x-range

      IF ( (x(1) .GT. xg(1)) .OR. (x(n) .LT. xg(ng)) ) THEN
          WRITE(0,*) '>>> ERROR (inter2) <<<  Data do not span '//      &
     &               'grid.  '
          WRITE(0,*) '                        Use ADDPNT to '//         &
     &               'expand data and re-run.'
          STOP
      ENDIF

!  find the integral of each grid interval and use this to
!  calculate the average y value for the interval
!  xgl and xgu are the lower and upper limits of the grid interval

      jstart = 1
      ngintv = ng - 1
      DO 50, i = 1,ngintv

! initalize:

            area = 0.0
            xgl = xg(i)
            xgu = xg(i+1)

!  discard data before the first grid interval and after the
!  last grid interval
!  for internal grid intervals, start calculating area by interpolating
!  between the last point which lies in the previous interval and the
!  first point inside the current interval

            k = jstart
            IF (k .LE. n-1) THEN

!  if both points are before the first grid, go to the next point
   30         CONTINUE
                IF (x(k+1) .LE. xgl) THEN
                   jstart = k - 1
                   k = k+1
                   IF (k .LE. n-1) GO TO 30
                ENDIF


!  if the last point is beyond the end of the grid, complete and go to the next
!  grid
   40         CONTINUE
                 IF ((k .LE. n-1) .AND. (x(k) .LT. xgu)) THEN

                    jstart = k-1

! compute x-coordinates of increment

                    a1 = MAX(x(k),xgl)
                    a2 = MIN(x(k+1),xgu)

!  if points coincide, contribution is zero

                    IF (x(k+1).EQ.x(k)) THEN
                       darea = 0.e0
                    ELSE
                       slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                       b1 = y(k) + slope*(a1 - x(k))
                       b2 = y(k) + slope*(a2 - x(k))
                       darea = (a2 - a1)*(b2 + b1)/2.
                    ENDIF


!  find the area under the trapezoid from a1 to a2

                    area = area + darea

! go to next point

                    k = k+1
                    GO TO 40

                ENDIF

            ENDIF

!  calculate the average y after summing the areas in the interval
            yg(i) = area/(xgu - xgl)

   50 CONTINUE
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('CH_INTER2',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE CH_INTER2
