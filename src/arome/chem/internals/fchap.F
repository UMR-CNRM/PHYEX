!     ######spl
      FUNCTION fchap(zeta)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      REAL zeta

! function value:
      REAL fchap

! internal:
      REAL rm
      INTEGER j
      REAL y(21)
!_______________________________________________________________________

      DATA y /                                                          &
     &     3.800,4.055,4.348,4.687,5.083,5.551,6.113,                   &
     &     6.799,7.650,8.732,10.144,12.051,14.730,18.686,               &
     &     24.905,35.466,55.211,96.753,197.,485.,1476./
!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FCHAP',0,ZHOOK_HANDLE)
      j = MIN( MAX(INT(zeta)+1,75),95 )
      rm = FLOAT(j)

      fchap = y(j-75) + (y(j-74) - y(j-75))*(zeta - (rm-1.))
!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('FCHAP',1,ZHOOK_HANDLE)
      RETURN
      ENDFUNCTION FCHAP
