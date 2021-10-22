!     ######spl
SUBROUTINE CH_AER_DRIVER(PM, PSIG0, PRG0, PN0, PCTOTG, PCTOTA,&
                           PCCTOT, PDTACT, PSEDA,&
                           PMU, PLAMBDA, PRHOP0, POM, PSO4RAT,  &
                           PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PMASK,&
                           PTIME, PSOLORG)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!#####################################################################################
!!
!!    PURPOSE
!!    -------
!!
!!    compute the right hand side of the moment equations
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Vincent Crassier (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_AER_COAG
USE MODI_CH_AER_GROWTH
USE MODI_CH_AER_SOLV
!!
!!  IMPLICIT ARGUMENTS
!!  ------------------
!!
USE MODD_CH_AEROSOL
!
!
IMPLICIT NONE
!  Declaration arguments
REAL,                                   INTENT(IN)    :: PDTACT, PTIME
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PRHOP0, POM
REAL,                 DIMENSION(:),     INTENT(INOUT) :: PLAMBDA, PMU, PSO4RAT
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PM
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSIG0, PRG0, PN0
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PSEDA
REAL,                 DIMENSION(:,:),   INTENT(INOUT) :: PMASK
REAL,                 DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL,                 DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
!
!  Declarations variables internes
!
INTEGER                       :: II

! Variables utilisees pour le tranfert de moment de chaque espece
! pour la condensation 
!----------------------------------------------------------------

REAL, DIMENSION(SIZE(PM,1),(JPMODE)*3)  :: ZDMINTRA,ZDMINTER,ZDMCOND

REAL                          :: ZGASMW       ! Molecular weight of background
                                           ! gas (g/mol) 
REAL, DIMENSION(SIZE(PM,1))   :: ZPGAS        ! background gas pressure (Pa)
REAL, DIMENSION(SIZE(PM,1))   :: ZRH,PSAT            ! Relative humidity
REAL                          :: ZDT          ! Pas de temps
REAL, DIMENSION(SIZE(PM,1))   :: ZPKM, ZPKH2O

!-----------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_DRIVER',0,ZHOOK_HANDLE)
ZDT=PDTACT

!******************************************************
!      Thermodynamic variables initialization
!         from Meso-NHC
!******************************************************
  
ZPKM(:) = 1E-3*PDENAIR(:) * 6.0221367E+23 / 28.9644
ZPKH2O(:) = ZPKM(:)*1.6077*PRV(:)
PSAT(:)=0.611*EXP(17.2694*(PTEMP(:)-273.16)/(PTEMP(:)-35.86))
PSAT(:)=PSAT(:)*1000.
ZRH(:)=(ZPKH2O(:)/(ZPKM(:)*1.6077))*PPRESSURE(:)/&
      &(0.622+(ZPKH2O(:)/(ZPKM(:)*1.6077)))/PSAT(:)

ZPGAS(:)=PPRESSURE(:)
ZGASMW=29.

!******************************************************
!      calculate gas viscosity and mean free path
!******************************************************
PMU(:)=0.003661*PTEMP(:)
PMU(:)=.0066164*PMU(:)*sqrt(PMU(:))/(PTEMP(:)+114.d0)

PLAMBDA(:)=PMU(:)/PDENAIR(:)*sqrt(1.89d-4*ZGASMW/PTEMP(:))*1.e6

CALL CH_AER_COAG(PM, PSIG0, PRG0, PN0,ZDMINTRA,ZDMINTER,&
                 PTEMP,PMU,PLAMBDA,PRHOP0)


CALL CH_AER_GROWTH(PM, PSIG0, PRG0, ZDMCOND,PDENAIR,ZGASMW,&
                     ZPGAS,PTEMP,ZRH,POM,PSO4RAT,PDTACT)

DO II=1,JPMODE
ZDMINTRA(:,NM0(II)) = ZDMINTRA(:,NM0(II)) * PMASK(:,II)
ZDMINTRA(:,NM3(II)) = ZDMINTRA(:,NM3(II)) * PMASK(:,II)
ZDMINTRA(:,NM6(II)) = ZDMINTRA(:,NM6(II)) * PMASK(:,II)
ZDMINTER(:,NM0(II)) = ZDMINTER(:,NM0(II)) * PMASK(:,II)
ZDMINTER(:,NM3(II)) = ZDMINTER(:,NM3(II)) * PMASK(:,II)
ZDMINTER(:,NM6(II)) = ZDMINTER(:,NM6(II)) * PMASK(:,II)
ZDMCOND(:,NM0(II)) = ZDMCOND(:,NM0(II)) * PMASK(:,II)
ZDMCOND(:,NM3(II)) = ZDMCOND(:,NM3(II)) * PMASK(:,II)
ZDMCOND(:,NM6(II)) = ZDMCOND(:,NM6(II)) * PMASK(:,II)
POM(:,II)          = POM(:,II)  * PMASK(:,II)
PSEDA(:,NM0(II))   = PSEDA(:,NM0(II))   * PMASK(:,II)
PSEDA(:,NM3(II))   = PSEDA(:,NM3(II))   * PMASK(:,II)
PSEDA(:,NM6(II))   = PSEDA(:,NM6(II))   * PMASK(:,II)
END DO


CALL CH_AER_SOLV(PM, PSIG0, PRG0, PN0,PCTOTG, PCTOTA, PCCTOT, &
                   ZDMINTRA,ZDMINTER,ZDMCOND,PSEDA,ZDT,POM,&
                   PRV, PDENAIR, PPRESSURE, PTEMP, PRC, PTIME, PSOLORG)


IF (LHOOK) CALL DR_HOOK('CH_AER_DRIVER',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_DRIVER
