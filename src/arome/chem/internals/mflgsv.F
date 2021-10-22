!     ######spl
      FUNCTION MFLGSV (IVALUE, ISET)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
      INTEGER MFLGSV
      LOGICAL ISET
      INTEGER IVALUE
!-----------------------------------------------------------------------
! MFLGSV saves and recalls the parameter MESFLG which controls the
! printing of the error messages.
!
! Saved local variable..
!
!   MESFLG = Print control flag..
!            1 means print all messages (the default).
!            0 means no printing.
!
! On input..
!
!   IVALUE = The value to be set for the MESFLG parameter,
!            if ISET is .TRUE. .
!
!   ISET   = Logical flag to indicate whether to read or write.
!            If ISET=.TRUE., the MESFLG parameter will be given
!            the value IVALUE.  If ISET=.FALSE., the MESFLG
!            parameter will be unchanged, and IVALUE is a dummy
!            parameter.
!
! On return..
!
!   The (old) value of the MESFLG parameter will be returned
!   in the function value, MFLGSV.
!
! This is a modification of the SLATEC library routine J4SAVE.
!
! Subroutines/functions called by MFLGSV.. None
!-----------------------------------------------------------------------
      INTEGER MESFLG
!-----------------------------------------------------------------------
! The following Fortran-77 declaration is to cause the values of the
! listed (local) variables to be saved between calls to this integrator.
!-----------------------------------------------------------------------
      SAVE MESFLG
      DATA MESFLG/1/
!
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('MFLGSV',0,ZHOOK_HANDLE)
      MFLGSV = MESFLG
      IF (ISET) MESFLG = IVALUE
      IF (LHOOK) CALL DR_HOOK('MFLGSV',1,ZHOOK_HANDLE)
      RETURN
!----------------------- End of Function MFLGSV ------------------------
      ENDFUNCTION MFLGSV
