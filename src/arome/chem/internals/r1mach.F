!     ######spl
      FUNCTION R1MACH (IDUM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###########################
      REAL R1MACH
      INTEGER IDUM
!-----------------------------------------------------------------------
! This routine computes the unit roundoff of the machine.
! This is defined as the smallest positive machine number
! u such that  1.0 + u .ne. 1.0
!
! Subroutines/functions called by R1MACH.. None
!-----------------------------------------------------------------------
      REAL U, COMP
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('R1MACH',0,ZHOOK_HANDLE)
      U = 1.0E0
 10   U = U*0.5E0
      COMP = 1.0E0 + U
      IF (COMP .NE. 1.0E0) GO TO 10
      R1MACH = U*2.0E0
      IF (LHOOK) CALL DR_HOOK('R1MACH',1,ZHOOK_HANDLE)
      RETURN
      ENDFUNCTION R1MACH
