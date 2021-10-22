!     ######spl
      SUBROUTINE XSETF (MFLAG)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ########################
!-----------------------------------------------------------------------
! This routine resets the print control flag MFLAG.
!
! Subroutines called by XSETF.. None
! Function routines called by XSETF.. MFLGSV
!-----------------------------------------------------------------------
      INTEGER MFLAG, JUNK, MFLGSV
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('XSETF',0,ZHOOK_HANDLE)
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = MFLGSV (MFLAG,.TRUE.)
      IF (LHOOK) CALL DR_HOOK('XSETF',1,ZHOOK_HANDLE)
      RETURN
      ENDSUBROUTINE XSETF
