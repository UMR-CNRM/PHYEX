!     ######spl
     SUBROUTINE CH_INI_ORILAM(PM, PSIG0, PRG0, PN0,PCTOTG, PCTOTA, PCCTOT, &
                              PSEDA, POM, PRHOP0, PAERO, PCHEM, PRV, PDENAIR,&
                              PPRESSURE,  PTEMP, PRC, PFRAC, PMI, GSCHEME)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ##########################################
!!
!!    PURPOSE
!!    -------
!!     initialize the aerosol variables (vectorwise) by calling NNARES
!!
!!    METHOD
!!    -------
!!      call the solver with zero coag/growth/cond terms
!!    then only ares should be active and we won't need to recode everyting
!!    here ;-)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre Tulet (GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!! 
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODI_CH_AER_SOLV
USE MODI_CH_AER_TRANS
USE MODD_CH_AEROSOL
USE MODD_CH_M9, ONLY : CNAMES
!
!*       0.     DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PM, POM 
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PFRAC
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PMI
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCHEM, PAERO
REAL, DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
CHARACTER(LEN=10),      INTENT(IN)    :: GSCHEME
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(SIZE(PM,1),(JPMODE)*3)    :: ZDMINTRA, ZDMINTER, ZDMCOND
REAL, DIMENSION(SIZE(PM,1),JPMODE)        :: ZMASK, ZSOLORG

INTEGER  :: JJ, JI
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_INI_ORILAM',0,ZHOOK_HANDLE)
POM(:,:)      = 0.
PFRAC(:,:)    = 0.
PSEDA(:,:)    = 0.
ZDMINTRA(:,:) = 0.
ZDMINTER(:,:) = 0.
ZDMCOND(:,:)  = 0.
ZSOLORG(:,:)  = 0.
ZMASK(:,:)    = 1.
!
! Initialization of constants
XPI =  2.*ASIN(1.)
XBOLTZ = 1.380658E-23
XAVOGADRO = 6.0221367E+23
XG = 9.80665
XP00 = 1.E5
XMD = 28.9644E-3
XRD = XAVOGADRO * XBOLTZ / XMD
XCPD = 7.* XRD /2.
!
! Moments index
NM0(1) = 1
NM3(1) = 2
NM6(1) = 3
NM0(2) = 4
NM3(2) = 5
NM6(2) = 6
!
! Aerosol Density
! Cf Ackermann (all to black carbon except water)
XRHOI(:) = 1.8e3
XRHOI(JP_AER_H2O) = 1.0e3   ! water
!
DO JJ=1,NSP+NCARB+NSOA
  XFAC(JJ)=(4./3.)*XPI*XRHOI(JJ)*1.e-9
ENDDO
!
! verify that all array elements are defined
DO JI = 1, SIZE(XRHOI)
  IF (XRHOI(JI) .LE. 0.0) THEN
    PRINT *, 'CH_AER_MOD_INIT ERROR: density for species ', JI, ' not defined'
    STOP 'CH_AER_MOD_INIT ERROR: density not defined'
  END IF
ENDDO
!
! Index gas scheme <=> Index Orilam
!
JP_CH_SO42M  = 0 ! unuse in many schemes
!
DO JJ=1,SIZE(CNAMES)
! for heterogeneous chemistry
IF (CNAMES(JJ) == "O3") JP_CH_O3  = JJ
IF (CNAMES(JJ) == "SO2") JP_CH_SO2  = JJ
IF (CNAMES(JJ) == "SO42M") JP_CH_SO42M  = JJ
IF (CNAMES(JJ) == "H2O2") JP_CH_H2O2  = JJ

! Inorganics
IF (CNAMES(JJ) == "HNO3") JP_CH_HNO3  = JJ
IF (CNAMES(JJ) == "NH3")  JP_CH_NH3   = JJ
IF ((CNAMES(JJ) == "H2SO4").OR.(CNAMES(JJ) == "SULF")) JP_CH_H2SO4 = JJ

! SOA group 1
IF (CNAMES(JJ) == "URG1") JP_CH_URG1 = JJ
IF (CNAMES(JJ) == "UR21") JP_CH_UR21 = JJ
IF (CNAMES(JJ) == "UR28") JP_CH_UR28 = JJ

! SOA group 2
IF (CNAMES(JJ) == "URG2") JP_CH_URG2 = JJ
IF (CNAMES(JJ) == "RPG2") JP_CH_RPG2 = JJ
IF (CNAMES(JJ) == "RP18") JP_CH_RP18 = JJ
IF (CNAMES(JJ) == "UR29") JP_CH_UR29 = JJ
IF (CNAMES(JJ) == "UR30") JP_CH_UR30 = JJ
IF (CNAMES(JJ) == "RP13") JP_CH_RP13 = JJ
IF (CNAMES(JJ) == "RP17") JP_CH_RP17 = JJ

! SOA group 3
IF (CNAMES(JJ) == "RPG3") JP_CH_RPG3 = JJ
IF (CNAMES(JJ) == "RPR9") JP_CH_RPR9 = JJ
IF (CNAMES(JJ) == "RP12") JP_CH_RP12 = JJ

! SOA group 4
IF (CNAMES(JJ) == "URG4") JP_CH_URG4 = JJ
IF (CNAMES(JJ) == "UR8")  JP_CH_UR8  = JJ ! only for MPMPO (for PUN it is group 10)
IF (CNAMES(JJ) == "UR3")  JP_CH_UR3  = JJ
IF (CNAMES(JJ) == "UR23") JP_CH_UR23 = JJ

! SOA group 5
IF (CNAMES(JJ) == "UR17") JP_CH_UR17 = JJ
IF (CNAMES(JJ) == "AP7")  JP_CH_AP7 = JJ
IF (CNAMES(JJ) == "UR7")  JP_CH_UR7  = JJ ! only for MPMPO (for PUN it is group 10)
IF (CNAMES(JJ) == "RPR3") JP_CH_RPR3 = JJ ! only for PUN   (for MPMPO it is not a SOA precursor)

! SOA group 6
IF (CNAMES(JJ) == "URG6") JP_CH_URG6 = JJ
IF (CNAMES(JJ) == "ARAC") JP_CH_ARAC = JJ
IF (CNAMES(JJ) == "UR22") JP_CH_UR22 = JJ ! only for PUN (for MPMPO it is not a SOA precursor)
IF (CNAMES(JJ) == "UR31") JP_CH_UR31 = JJ
IF (CNAMES(JJ) == "AP1")  JP_CH_AP1  = JJ
IF (CNAMES(JJ) == "AP6")  JP_CH_AP6  = JJ

! SOA group 7
IF (CNAMES(JJ) == "URG7") JP_CH_URG7 = JJ
IF (CNAMES(JJ) == "RPG7") JP_CH_RPG7 = JJ
IF (CNAMES(JJ) == "RPR7") JP_CH_RPR7 = JJ
IF (CNAMES(JJ) == "RPR4") JP_CH_RPR4 = JJ ! only for PUN (for MPMPO it is not a SOA precursor)
IF (CNAMES(JJ) == "RP14") JP_CH_RP14 = JJ ! only for PUN (for MPMPO it is not a SOA precursor)
IF (CNAMES(JJ) == "RP19") JP_CH_RP19 = JJ ! only for PUN (for MPMPO it is not a SOA precursor)
IF (CNAMES(JJ) == "ADAC") JP_CH_ADAC = JJ
IF (CNAMES(JJ) == "UR2")  JP_CH_UR2  = JJ
IF (CNAMES(JJ) == "UR14") JP_CH_UR14 = JJ
IF (CNAMES(JJ) == "UR27") JP_CH_UR27 = JJ

! SOA group 8
IF (CNAMES(JJ) == "URG8") JP_CH_URG8 = JJ
IF (CNAMES(JJ) == "UR19") JP_CH_UR19 = JJ ! only for MPMPO (for PUN it is not a SOA precursor)
IF (CNAMES(JJ) == "UR11") JP_CH_UR11 = JJ
IF (CNAMES(JJ) == "UR15") JP_CH_UR15 = JJ
IF (CNAMES(JJ) == "AP10") JP_CH_AP10 = JJ

!  SOA group 9
IF (CNAMES(JJ) == "URG9") JP_CH_URG9 = JJ
IF (CNAMES(JJ) == "UR20") JP_CH_UR20 = JJ
IF (CNAMES(JJ) == "UR34") JP_CH_UR34 = JJ
IF (CNAMES(JJ) == "AP11") JP_CH_AP11 = JJ
IF (CNAMES(JJ) == "AP12") JP_CH_AP12 = JJ
IF (CNAMES(JJ) == "UR26") JP_CH_UR26 = JJ 

!  SOA group 10
IF (CNAMES(JJ) == "URG10") JP_CH_URG10 = JJ
IF (CNAMES(JJ) == "PAN8")  JP_CH_PAN8  = JJ  ! only for PUN (for MPMPO it is not a SOA precursor)
IF (CNAMES(JJ) == "UR5")   JP_CH_UR5   = JJ
IF (CNAMES(JJ) == "UR6")   JP_CH_UR6   = JJ
IF (CNAMES(JJ) == "UR7")   JP_CH_UR7   = JJ
IF (CNAMES(JJ) == "UR8")   JP_CH_UR8   = JJ
IF (CNAMES(JJ) == "AP8")   JP_CH_AP8   = JJ

END DO
!
!
!*      0.4    initialization aerosol solveur
CALL CH_AER_TRANS(0, PM, PSIG0, PRG0, PN0, PRHOP0, PAERO,&
                  PCHEM, PCTOTG, PCTOTA, PCCTOT, PFRAC, PMI, ZMASK, GSCHEME )

CALL CH_AER_SOLV(PM,PSIG0, PRG0, PN0, PCTOTG, PCTOTA, PCCTOT, &
                 ZDMINTRA,ZDMINTER,ZDMCOND,PSEDA,0.,POM,&
                 PRV, PDENAIR, PPRESSURE, PTEMP, PRC, 0., ZSOLORG)

CALL CH_AER_TRANS(1, PM, PSIG0, PRG0, PN0, PRHOP0, PAERO,&
                  PCHEM, PCTOTG, PCTOTA, PCCTOT, PFRAC, PMI, ZMASK, GSCHEME)

!
IF (LHOOK) CALL DR_HOOK('CH_INI_ORILAM',1,ZHOOK_HANDLE)
RETURN
!
END SUBROUTINE CH_INI_ORILAM
