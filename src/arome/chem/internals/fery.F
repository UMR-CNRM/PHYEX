!     ######spl
      FUNCTION fery(w)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK


      IMPLICIT NONE

! input:
      REAL w

! function value:
      REAL fery
!_______________________________________________________________________

      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('FERY',0,ZHOOK_HANDLE)
      IF (w .LT. 250.) THEN
          fery = 1.
! outside the ery spectrum range
      ELSEIF ((w .GE. 250.) .AND. (w .LT. 298)) THEN
          fery = 1.
      ELSEIF ((w .GE. 298.) .AND. (w .LT. 328.)) THEN
          fery = 10.**( 0.094*(298.-w) )
      ELSEIF ((w .GE. 328.) .AND. (w .LT. 400.)) THEN
          fery = 10.**( 0.015*(139.-w) )
      ELSE
         fery = 1.E-36
! outside the ery spectrum range
      ENDIF

!_______________________________________________________________________

      IF (LHOOK) CALL DR_HOOK('FERY',1,ZHOOK_HANDLE)
      RETURN
      ENDFUNCTION FERY
