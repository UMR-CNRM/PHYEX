!!    ##############################
      SUBROUTINE AROINI_NSV(KSV,KSV_CHEMBEG, KSV_CHEMEND, KSV_AERBEG, KSV_AEREND, &
                            KSV_DSTBEG, KSV_DSTEND,KSV_DSTDEPBEG, KSV_DSTDEPEND,&
                            KSV_CO2)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!!    ##############################
!!
!!*** *AROINI_MNHC*
!!
!!    PURPOSE
!!    -------
!        initialize nsv
!!
!!
!!    AUTHOR
!!    ------
!!    P. Tulet  *CNRM / GMEI* 
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/04/05
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_NSV
USE MODD_CH_AEROSOL, ONLY : LORILAM, CAERONAMES, LVARSIGI, LVARSIGJ, NM6_AER
USE MODD_CH_MNHC_n,  ONLY : LUSECHEM
USE MODD_CH_M9,      ONLY : CNAMES
USE MODD_DUST
IMPLICIT NONE 
!
!*       0.   Declarations of arguments
!
INTEGER, INTENT(IN)  :: KSV,KSV_CHEMBEG, KSV_CHEMEND, KSV_AERBEG, KSV_AEREND, &
                        KSV_DSTBEG, KSV_DSTEND, KSV_CO2, &
                        KSV_DSTDEPBEG, KSV_DSTDEPEND


INTEGER :: JN, ICO2
!-------------------------------------------------------------------------------
!
! Initialize NSV and CSV (names of chemical aerosols and dusts species)
! 
!
! Initialize NSV
! 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AROINI_NSV',0,ZHOOK_HANDLE)
NSV=KSV  ! for LIMA case
IF (LUSECHEM) THEN
  NSV_CHEMBEG     = KSV_CHEMBEG
  NSV_CHEMEND     = KSV_CHEMEND
  NSV_CHEM = KSV_CHEMEND - KSV_CHEMBEG + 1
  NSV=NSV_CHEM
ELSE
! force First index to be superior to last index
! in order to create a null section
  NSV_CHEMBEG     = 1
  NSV_CHEMEND     = 0
  NSV_CHEM        = 0
END IF

IF (LORILAM) THEN
  NSV_AERBEG     = KSV_AERBEG
  NSV_AEREND     = KSV_AEREND
  NSV_AER = KSV_AEREND - KSV_AERBEG + 1
  NM6_AER = 0
  NSV=NSV_CHEM+NSV_AER 
  IF (LVARSIGI) NM6_AER = 1
  IF (LVARSIGJ) NM6_AER = NM6_AER + 1
ELSE
! force First index to be superior to last index
! in order to create a null section
  NSV_AERBEG     = 1
  NSV_AEREND     = 0
  NSV_AER        = 0 
  NM6_AER = 0 
END IF

IF (LDUST) THEN
  NSV_DSTBEG     = KSV_DSTBEG
  NSV_DSTEND     = KSV_DSTEND
  NSV_DST = KSV_DSTEND - KSV_DSTBEG + 1
  NSV=NSV_CHEM+NSV_AER+NSV_DST
ELSE
! force First index to be superior to last index
! in order to create a null section
  NSV_DSTBEG     = 1
  NSV_DSTEND     = 0
  NSV_DST        = 0
END IF 

IF ((LDUST).AND.(LDEPOS_DST(1))) THEN
  NSV_DSTDEPBEG     = KSV_DSTDEPBEG
  NSV_DSTDEPEND     = KSV_DSTDEPEND
  NSV_DSTDEP = KSV_DSTDEPEND - KSV_DSTDEPBEG + 1
  NSV  = NSV_DST + NSV_DSTDEP + NSV_CHEM + NSV_AER 
ELSE
! force First index to be superior to last index
! in order to create a null section
  NSV_DSTDEPBEG     = 1
  NSV_DSTDEPEND     = 0
  NSV_DSTDEP        = 0
ENDIF

NSV_CO2     = KSV_CO2
ICO2 = 0
IF (NSV_CO2 > 0) THEN
  ICO2 = 1
  NSV  = NSV_DST + NSV_DSTDEP + NSV_CHEM + NSV_AER + ICO2
ENDIF
!
! Initialize CSV
! 
IF ((LUSECHEM).OR.(LORILAM).OR.(LDUST).OR.(LDEPOS_DST(1))) THEN
  IF (.NOT.ASSOCIATED(CSV)) THEN
    ALLOCATE(CSV(NSV))
  ELSE 
    DEALLOCATE(CSV)
    ALLOCATE(CSV(NSV))
  ENDIF

  DO JN=NSV_CHEMBEG,NSV_CHEMEND
   CSV(JN)(1:6) = CNAMES(JN-NSV_CHEMBEG+1)(1:6)
  END DO

  DO JN=NSV_AERBEG,NSV_AEREND
   CSV(JN)(1:6) = CAERONAMES(JN-NSV_AERBEG+1)(1:6)
  END DO

  DO JN=NSV_DSTBEG,NSV_DSTEND
   CSV(JN)(1:6) = CDUSTNAMES(JN-NSV_DSTBEG+1)(1:6)
  END DO

  DO JN=NSV_DSTDEPBEG,NSV_DSTDEPEND
   CSV(JN)(1:6) = CDEDSTNAMES(JN-NSV_DSTDEPBEG+1)(1:6)
  END DO
ELSE
  IF (.NOT.ASSOCIATED(CSV)) THEN
    ALLOCATE(CSV(0))
  ELSE 
    DEALLOCATE(CSV)
    ALLOCATE(CSV(0))
  ENDIF
ENDIF
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AROINI_NSV',1,ZHOOK_HANDLE)
END SUBROUTINE AROINI_NSV
