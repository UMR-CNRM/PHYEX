!     ######spl
      SUBROUTINE addpnt ( x, y, ld, n, xnew, ynew )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! calling parameters

      INTEGER ld, n
      REAL x(ld), y(ld)
      REAL xnew, ynew
      INTEGER ierr

! local variables

      INTEGER insert
      INTEGER i


! initialize error flag

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('ADDPNT',0,ZHOOK_HANDLE)
      ierr = 0

! check n<ld to make sure x will hold another point

      IF (n .GE. ld) THEN
         WRITE(0,*) '>>> ERROR (ADDPNT) <<<  Cannot expand array '
         WRITE(0,*) '                        All elements used.'
         STOP
      ENDIF

      insert = 1
      i = 2

! check, whether x is already sorted.
! also, use this loop to find the point at which xnew needs to be inserted
! into vector x, if x is sorted.

 10   CONTINUE
      IF (i .LT. n) THEN
        IF (x(i) .LT. x(i-1)) THEN
           WRITE(0,*) '>>> ERROR (ADDPNT) <<<  x-data must be '//       &
     &                'in ascending order!'
           STOP
        ELSE
           IF (xnew .GT. x(i)) insert = i + 1
        ENDIF
        i = i+1
        GOTO 10
      ENDIF

! if <xnew,ynew> needs to be appended at the end, just do so,
! otherwise, insert <xnew,ynew> at position INSERT

      IF ( xnew .GT. x(n) ) THEN

         x(n+1) = xnew
         y(n+1) = ynew

      ELSE

! shift all existing points one index up

         DO i = n, insert, -1
           x(i+1) = x(i)
           y(i+1) = y(i)
         ENDDO

! insert new point

         x(insert) = xnew
         y(insert) = ynew

      ENDIF

! increase total number of elements in x, y

      n = n+1

      IF (LHOOK) CALL DR_HOOK('ADDPNT',1,ZHOOK_HANDLE)
      ENDSUBROUTINE ADDPNT
