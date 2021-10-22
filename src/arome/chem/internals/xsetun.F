!     ######spl
      SUBROUTINE XSETUN (LUN)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #######################
!-----------------------------------------------------------------------
! This routine resets the logical unit number for messages.
!
! Subroutines called by XSETUN.. None
! Function routines called by XSETUN.. LUNSAV
!-----------------------------------------------------------------------
      INTEGER LUN, JUNK, LUNSAV
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('XSETUN',0,ZHOOK_HANDLE)
      IF (LUN .GT. 0) JUNK = LUNSAV (LUN,.TRUE.)
      IF (LHOOK) CALL DR_HOOK('XSETUN',1,ZHOOK_HANDLE)
      RETURN
!----------------------- End of Subroutine XSETUN ----------------------
      ENDSUBROUTINE XSETUN
