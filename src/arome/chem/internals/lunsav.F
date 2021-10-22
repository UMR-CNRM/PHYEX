!     ######spl
      FUNCTION LUNSAV (IVALUE, ISET)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      INTEGER LUNSAV
      LOGICAL ISET
      INTEGER IVALUE
!-----------------------------------------------------------------------
! LUNSAV saves and recalls the parameter LUNIT which is the logical
! unit number to which error messages are printed.
!
! Saved local variable..
!
!  LUNIT   = Logical unit number for messages.
!            The default is 6 (machine-dependent).
!
! On input..
!
!   IVALUE = The value to be set for the LUNIT parameter,
!            if ISET is .TRUE. .
!
!   ISET   = Logical flag to indicate whether to read or write.
!            If ISET=.TRUE., the LUNIT parameter will be given
!            the value IVALUE.  If ISET=.FALSE., the LUNIT
!            parameter will be unchanged, and IVALUE is a dummy
!            parameter.
!
! On return..
!
!   The (old) value of the LUNIT parameter will be returned
!   in the function value, LUNSAV.
!
! This is a modification of the SLATEC library routine J4SAVE.
!
! Subroutines/functions called by LUNSAV.. None
!-----------------------------------------------------------------------
      INTEGER LUNIT
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE LUNIT
      DATA LUNIT/6/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('LUNSAV',0,ZHOOK_HANDLE)
      LUNSAV = LUNIT
      IF (ISET) LUNIT = IVALUE
      IF (LHOOK) CALL DR_HOOK('LUNSAV',1,ZHOOK_HANDLE)
      RETURN
!----------------------- End of Function LUNSAV ------------------------
      ENDFUNCTION LUNSAV
