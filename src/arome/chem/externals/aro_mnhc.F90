!    ################################################### 
      SUBROUTINE ARO_MNHC(PRSVSIN, PRHODREF, PTSTEP,      &
                          PTHT, PABST,                  &
                          PRT, PLAT, PLON,              &
                          PALB_UV, PZS, PZENITH, PZZ,   &
                          KYEAR, KMONTH, KDAY, PTIME,   &
                          KLON,KLEV,NSV, KRR,           &
                          KTCOUNT, KLUOUT,NDIAG,PPEZDIAG,PRSVS )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ###################################################
!!
!!*** *ARO_MNHC*  monitor of the chemical module of MesoNH-C
!!
!!    PURPOSE
!!    -------
!!       The purpose of this subroutine is to control the chemical module
!!    i.e. to pass the meteorological parameters from MesoNH to its chemical
!!    part and to call the different subroutines (calculation of rate constants,
!!    photolysis rates, stiff solver,..)
!!
!!    METHOD
!!    ------
!!       The calculation  of the chemical terms is performed using a loop
!!    over all spatial dimensions. 
!!
!!       For each single grid point, all necessary meteorological parameters are
!!    passed into the chemical core system (variable TZM). This variable is
!!    then passed on to the subroutines that calculate the reaction and
!!    photolysis rates. Then the chemical solver is called. As the chemistry
!!    part works with different units than MesoNH (MesoNH uses mixing ratio,
!!    the chemisty part uses molec/cm3) some unit conversion is also performed.
!!
!!       Temporal integration is performed over a double timestep 2*PTSTEP
!!    (except in the case of a cold start). If the timestep of MesoNH
!!    is too large for the chemical solver, several smaller steps can
!!    be taken using the NCH_SUBSTEPS parameter.
!!    "SPLIT"  : from PRSVSIN the scalar variable at t+dt is calculated and
!!               given as input to the solver; the result is written 
!!               into PRSVS; this corresponds to applying first only dynamics
!!               and then only chemistry; this option assures positivity, but
!!               degrades the order of the temporal integration.
!!               In fact, an overhead of a factor two is produced here.
!!               A future solution will be to calculate the dynamics
!!               of the scalar variables not using leapfrog, but forward
!!               temporal integration.
!!     Vectorization Mask : Need ISVECMASK for photolysis (input(IN))
!!     ISVECNMASK = Nb of vector 
!!     IDT1       = Nb of points for the first dimension  (for arome = size of vector)
!!     IDT2       = Nb of points for the second dimension (for arome = 1)
!!     IDT3       = Nb of points in vertical
!!
!!
!!
!!    REFERENCE
!!    ---------
!!    Book 1, 2, 3 of MesoNH-chemistry
!!
!!    AUTHOR
!!    ------
!!    P. Tulet  *CNRM / GMEI* and contributors of MesoNH-C (K. Shure, C. Mari, V. Crassier)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 10/11/04
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_METEO_TRANS
USE MODI_CH_SET_RATES
USE MODI_CH_SET_PHOTO_RATES
USE MODI_CH_SOLVER_n
USE MODI_CH_AER_SEDIM_n 
USE MODI_CH_ORILAM
USE MODI_CH_INI_ORILAM
USE MODI_CH_AER_EQM_CORMASS
!USE MODI_CH_AER_SURF
USE MODI_CH_AER_INIT
USE MODI_CH_UPDATE_JVALUES_n
!
!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
USE MODD_CH_M9_SCHEME, ONLY: JP_CO, JP_O3, JP_NO, JP_NO2, JP_SO2, JP_CH4
USE MODD_CH_M9,     ONLY: NEQ,            &! number of prognostic chem. species
                          NMETEOVARS,     &! number of meteorological variables
                          CNAMES,         &! names of the chem. species
                          METEOTRANSTYPE   ! type for meteo . transfer
!!
USE MODD_CH_INIT_JVALUES, ONLY: JPJVMAX    ! number of photolysis reactions in TUVMAIN
USE MODD_CH_MNHC_n, ONLY: NCH_SUBSTEPS
                  ! ZDTSOLVER = PTSTEP/NCH_SUBSTEPS
USE MODD_CH_MNHC_n, ONLY: LCH_INIT_FIELD
USE MODD_CH_MNHC_n, ONLY: CCH_VEC_METHOD
USE MODD_CH_SOLVER_n
USE MODD_CH_MNHC_n, ONLY: LCH_TUV_ONLINE, CCH_TUV_LOOKUP, CCH_TUV_CLOUDS,  &
                          XCH_TUV_ALBNEW, XCH_TUV_DOBNEW, XCH_TUV_TUPDATE
USE MODD_CH_MNHC_n, ONLY: CCH_SCHEME

!!
USE MODD_CONF,      ONLY: CCONF, & ! Configuration of models (START/RESTART)
                          NVERB    ! Level of informations on output-listing

!!
USE MODD_CST,       ONLY: XAVOGADRO, &! Avogadro number
                          XMD,&         ! Molar mass of dry air
                          XP00, XRD, XCPD, XPI
!
USE MODD_NSV,       ONLY: NSV_CHEMBEG, NSV_CHEMEND, NSV_CHEM, &
                          NSV_AERBEG, NSV_AEREND, NSV_AER
USE MODD_GRID_AROME,ONLY: XLAT0, XLON0
USE MODD_CH_AEROSOL
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
INTEGER,   INTENT(IN)   :: KLON     !NPROMA under CPG 
INTEGER,   INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,   INTENT(IN)   :: NSV     ! Number of passive scalar
INTEGER,   INTENT(IN)   :: KRR      ! Number of moist variables
REAL, DIMENSION(KLON,1,KLEV,NSV),INTENT(IN) :: PRSVSIN       ! source of scalar variable
REAL, DIMENSION(KLON,1,KLEV),    INTENT(IN)    :: PRHODREF    ! iteration count
REAL,    INTENT(IN)    :: PTSTEP      ! time step of MesoNH
REAL, DIMENSION(KLON,1,KLEV),    INTENT(IN)    :: PTHT, PABST ! theta and pressure at t
REAL, DIMENSION(KLON,1,KLEV,KRR),INTENT(IN)    :: PRT         ! moist variables at t
REAL, DIMENSION(KLON,1),           INTENT(IN)    :: PLAT, PLON  ! lat, lon of each points
REAL, DIMENSION(KLON,1),           INTENT(IN)    :: PALB_UV, PZS, PZENITH 
REAL, DIMENSION(KLON,1,KLEV),    INTENT(IN)    :: PZZ
INTEGER, INTENT(IN)    :: KYEAR       ! Current Year
INTEGER, INTENT(IN)    :: KMONTH      ! Current Month
INTEGER, INTENT(IN)    :: KDAY        ! Current Day
REAL,    INTENT(IN)    :: PTIME       ! Current time in second
INTEGER, INTENT(IN)    :: KTCOUNT     ! iteration count
INTEGER, INTENT(IN)    :: KLUOUT      ! unit for output listing count
INTEGER, INTENT(IN)   :: NDIAG   ! nb of diagnostics
REAL, DIMENSION(KLON,KLEV,NDIAG),INTENT(INOUT)   :: PPEZDIAG !O [-] Diagnostics table
REAL, DIMENSION(KLON,1,KLEV,NSV),INTENT(OUT) :: PRSVS       ! source of scalar variable

!
!
!*      0.2    declarations of local variables
!
INTEGER :: JM,JI,JK,JL,JK1  ! loop counters
REAL    :: ZDTSOLVER        ! timestep for the solver
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCHEM, ZOLDCHEM, ZNEWCHEM 
REAL, DIMENSION(:,:), ALLOCATABLE :: ZAERO
        ! arrays for parameter passage to solver
!
REAL    :: ZDEN2MOL
        !  ZDEN2MOL = 6.0221367E+23 * 1E-6 / 28.9644E-3
        !  conversion factor density to mol/cm3
        !  n_molec (moelc./cm3):  M = 1E-6*RHO(kg/m3) * XAVOGADRO / XMD
!
TYPE(METEOTRANSTYPE), DIMENSION(:), ALLOCATABLE  :: TZM 
        ! meteo variables to be transferred into CCS
!
INTEGER                :: JN, JSV ! loop index for SV
INTEGER                :: ILAT
!
!-------------------------------------------------------------------------------
!   variables for the vectorization
!
INTEGER, SAVE :: ISVECNMASK
INTEGER, SAVE :: ISVECNPT = 0
INTEGER, SAVE, DIMENSION(:,:), ALLOCATABLE :: ISVECMASK
!
INTEGER :: IDTI,IDTJ,IDTK
INTEGER :: IDT1,IDT2,IDT3
INTEGER             :: IIU  ! Upper dimension in x direction
INTEGER             :: IJU  ! Upper dimension in y direction
INTEGER             :: IKU  ! Upper dimension in z direction
INTEGER             :: IIB  ! indice I Beginning in x direction
INTEGER             :: IJB  ! indice J Beginning in y direction
INTEGER             :: IKB  ! indice K Beginning in z direction
INTEGER             :: IIE  ! indice I End       in x direction
INTEGER             :: IJE  ! indice J End       in y direction
INTEGER             :: IKE  ! indice K End       in z direction

!-------------------------------------------------------------------------------
!
! variables for TUV
!
REAL, DIMENSION(KLON,1,KLEV,JPJVMAX)  :: ZJVALUES
!--------------------------Arguments vs Local--------------------------------------------------
INTEGER            :: IDAY, IMONTH, IYEAR
REAL               :: ZTIME, ZTIMERAD, ZUT, ZDATE, ZAD, ZTSIDER, ZA1, ZA2
REAL, DIMENSION(KLON,1,KLEV)  :: ZRHODREF, ZZZ, ZTHT, ZABST
REAL, DIMENSION(KLON,1)         :: ZLAT, ZLON, ZTUT
REAL, DIMENSION(KLON,1,KLEV,KRR)  :: ZRT
REAL, DIMENSION(KLON,1)             :: ZZENITH, ZZS, ZALB_UV
REAL, DIMENSION(KLON,1)             :: ZCOSZEN, ZSOLANG
REAL                                :: ZSINDEL, ZCOSDEL, ZDECSOL
INTEGER, DIMENSION(0:11) :: IBIS, INOBIS ! Cumulative number of days per month
                                         ! for bissextile and regular years
!-------------------------------------------------------------------------------
! variables for the aerosol module
REAL, ALLOCATABLE, DIMENSION(:,:)   :: ZM, ZSIG0, ZN0, ZSOLORG, ZMASS, &
                                       ZRG0, ZCTOTG, ZSEDA, ZFRAC, ZMI  ! work array for aerosols
REAL, ALLOCATABLE, DIMENSION(:,:,:) :: ZCTOTA, ZCCTOT
                                         ! first dimension is vectorization,
                                         ! second dim. are the modes*moments
REAL, ALLOCATABLE, DIMENSION(:) :: ZRV, ZDENAIR, ZPRESSURE, ZTEMP, ZRC
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRHOP0, ZOM
REAL, DIMENSION(:),   ALLOCATABLE :: ZLAMBDA, ZMU, ZSO4RAT

REAL,DIMENSION(KLON,1,KLEV,NSV)  :: ZSVT
REAL,DIMENSION(KLON,1,KLEV,JPIN) :: ZXSEDA
REAL,DIMENSION(KLON,1,JPIN)        :: ZVDEPAERO
REAL,DIMENSION(KLON,KLEV,NDIAG)  :: ZPEZDIAG 
!
!-------------------------------------------------------------------------------
!
!
!*       1.    PREPARE MONITOR
!              ---------------
!     Arome is vectorized upon horizontal (I,1,K)
!     I = horizontal vector (not necessary equal to west/east nb of points)
!     K = nb of points in vertical and physical domain upon vertical is 2=>KLEV+1

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_MNHC',0,ZHOOK_HANDLE)
ZJVALUES(:,:,:,:) = 0.
ZVDEPAERO(:,:,:)  = 0.

ZZENITH(:,:) = PZENITH(:,:)
ZTIME        = PTIME
ZRT(:,:,:,:) = PRT(:,:,:,:)
ZTHT(:,:,:)  = PTHT(:,:,:)
ZABST(:,:,:) = PABST(:,:,:)
ZALB_UV(:,:) = PALB_UV(:,:)
ZZS(:,:)     = PZS(:,:)
ZZZ(:,:,:)   = PZZ(:,:,:)
ZRHODREF(:,:,:) = PRHODREF(:,:,:)
IDAY   = KDAY
IMONTH = KMONTH
IYEAR  = KYEAR

IIU = KLON
IJU = 1
IKU = KLEV
IIB=1
IJB=1
IKB=1
IIE = IIU
IJE = 1
IKE = IKU

IDTI=(IIE-IIB+1)
IDTJ=(IJE-IJB+1)
IDTK=(IKE-IKB+1)

ISVECNPT=IDTI*IDTJ*IDTK
ISVECNMASK=1
IDT1=IDTI
IDT2=IDTJ
IDT3=IDTK

! Possible to change the vector size only on NPROMA
! to do it uses:
!     ISVECNPT=IDTI*IDTJ
!     ISVECNMASK=IDTK
!     IDT1=IDTI
!     IDT2=IDTJ
!     IDT3=1 

IF (.NOT. ALLOCATED(ISVECMASK)) ALLOCATE(ISVECMASK(6,ISVECNMASK))

!
!**********************************
! Compute mask boundaries
!**********************************
!
  ISVECMASK(1,1)=IIB
  ISVECMASK(2,1)=IIB+IDT1-1
  ISVECMASK(3,1)=IJB
  ISVECMASK(4,1)=IJB+IDT2-1
  ISVECMASK(5,1)=IKB
  ISVECMASK(6,1)=IKB+IDT3-1

 IF (ISVECNMASK .GE. 2) THEN
    DO JI=2,ISVECNMASK
    ISVECMASK(1,JI)=ISVECMASK(1,JI-1)+IDT1-IIB
    ISVECMASK(3,JI)=ISVECMASK(3,JI-1)
    ISVECMASK(5,JI)=ISVECMASK(5,JI-1)
!
    ISVECMASK(3,JI)=ISVECMASK(3,JI)+(ISVECMASK(1,JI)/IDTI)*IDT2-IJB
    ISVECMASK(5,JI)=ISVECMASK(5,JI)+(ISVECMASK(3,JI)/IDTJ)*IDT3-IKB
!
    ISVECMASK(1,JI)=ISVECMASK(1,JI)-IDTI*(ISVECMASK(1,JI)/IDTI)+IIB
    ISVECMASK(3,JI)=ISVECMASK(3,JI)-IDTJ*(ISVECMASK(3,JI)/IDTJ)+IJB
    ISVECMASK(5,JI)=ISVECMASK(5,JI)-IDTK*(ISVECMASK(5,JI)/IDTK)+IKB
!
    ISVECMASK(2,JI)=ISVECMASK(1,JI)+IDT1-1
    ISVECMASK(4,JI)=ISVECMASK(3,JI)+IDT2-1
    ISVECMASK(6,JI)=ISVECMASK(5,JI)+IDT3-1
  END DO

 END IF

!*       1.1   calculate timestep variables
!
ZDTSOLVER = PTSTEP / NCH_SUBSTEPS
!
PRSVS(:,:,:,NSV_CHEMBEG:NSV_CHEMEND) = MAX(PRSVSIN(:,:,:,NSV_CHEMBEG:NSV_CHEMEND), 1.E-80) 

IF ((KTCOUNT == 1).AND.(LCH_INIT_FIELD)) THEN
! Bulk initialization; to be improve 
   PRSVS(:,:,:,:) = MAX(PRSVS(:,:,:,:),1E-20 / PTSTEP)
   PRSVS(:,:,:,JP_CO)  = 50E-9 / PTSTEP  ! CO initialize to 50 ppb (low value)
   PRSVS(:,:,:,JP_O3)  = 30E-9 / PTSTEP  ! O3 initialize to 30 ppb (low value)
   PRSVS(:,:,:,JP_NO)  = 1E-11 / PTSTEP  ! NO initialize to 1E-3 ppb (low value)
   PRSVS(:,:,:,JP_NO2) = 1E-11 / PTSTEP  ! NO2 initialize to 1E-3 ppb (low value)
   PRSVS(:,:,:,JP_SO2) = 1E-11 / PTSTEP  ! SO2 initialize to 1E-3 ppb (low value)
   PRSVS(:,:,:,JP_CH4) = 1700E-9 / PTSTEP ! CH4 mean value
! increase ozone over tropopause
PRSVS(:,:,11,JP_O3) = 5.*PRSVS(:,:,11,JP_O3)
PRSVS(:,:,9:10,JP_O3) = 15.*PRSVS(:,:,9:10,JP_O3)
PRSVS(:,:,1:8,JP_O3) = 40.*PRSVS(:,:,1:8,JP_O3)
! increase primary species in the BLM
PRSVS(:,:,IKE-8:IKE,JP_NO) = 10.*PRSVS(:,:,IKE-8:IKE,JP_NO)
PRSVS(:,:,IKE-8:IKE,JP_NO2) = 10.*PRSVS(:,:,IKE-8:IKE,JP_NO2)
PRSVS(:,:,IKE-8:IKE,JP_SO2) = 10.*PRSVS(:,:,IKE-8:IKE,JP_SO2)
END IF

IF (LORILAM) THEN
 IF (.NOT.ALLOCATED(XRHOI))  ALLOCATE(XRHOI(NSP+NCARB+NSOA))
 IF (.NOT.ALLOCATED(XFAC))  ALLOCATE(XFAC(NSP+NCARB+NSOA))
 ! Density: all to black carbon except water
 XRHOI(:) = 1.8e3
 XRHOI(JP_AER_H2O) = 1.0e3   ! water
 ! Moments index
 NM0(1) = 1
 NM3(1) = 2
 NM6(1) = 3 
 NM0(2) = 4 
 NM3(2) = 5 
 NM6(2) = 6

 DO JN=1,NSP+NCARB+NSOA
  XFAC(JN)=(4./3.)*3.14292654*XRHOI(JN)*1.e-9
 ENDDO

ZSVT(:,:,:,:) = PRSVS(:,:,:,:) * PTSTEP
ZSVT(:,:,:,NSV_AERBEG:NSV_AEREND) = MAX(ZSVT(:,:,:,NSV_AERBEG:NSV_AEREND), 1.E-80)

 IF ((KTCOUNT == 1).AND.(LAERINIT)) THEN
 CALL CH_AER_INIT(ZSVT(:,:,:,NSV_CHEMBEG:NSV_CHEMEND),&
                  ZSVT(:,:,:,NSV_AERBEG:NSV_AEREND),  &
                  PRHODREF )
 ENDIF

!
!*       2.5   sedimentation of aerosols tendency (ZXSEDA) 
!
 ZXSEDA(:,:,:,:)=0.
 IF ((LSEDIMAERO).AND.(KTCOUNT .GT. 1)) THEN
 CALL CH_AER_SEDIM_n(PTSTEP, ZSVT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_AERBEG:NSV_AEREND),&
                    PTHT(IIB:IIE,IJB:IJE,IKB:IKE), PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),&
                    PABST(IIB:IIE,IJB:IJE,IKB:IKE), ZVDEPAERO(IIB:IIE,IJB:IJE,:),    &
                    PZZ(IIB:IIE,IJB:IJE,IKB:IKE), ZXSEDA(IIB:IIE,IJB:IJE,IKB:IKE,:)) 
 END IF

PRSVS(:,:,:,:) = ZSVT(:,:,:,:) / PTSTEP

ENDIF
!
ALLOCATE(TZM(ISVECNPT))
ALLOCATE(ZCHEM(ISVECNPT,NEQ))
ALLOCATE(ZOLDCHEM(ISVECNPT,NEQ))
ALLOCATE(ZNEWCHEM(ISVECNPT,NEQ))
IF (LORILAM) THEN
ALLOCATE(ZAERO(ISVECNPT,NSV_AER))
ALLOCATE(ZM(ISVECNPT,JPIN))
ALLOCATE(ZSEDA(ISVECNPT,JPIN))
ALLOCATE(ZRHOP0(ISVECNPT,JPMODE))
ALLOCATE(ZSIG0(ISVECNPT,JPMODE))
ALLOCATE(ZRG0(ISVECNPT,JPMODE))
ALLOCATE(ZN0(ISVECNPT,JPMODE))
ALLOCATE(ZCTOTA(ISVECNPT,NSP+NCARB+NSOA,JPMODE))
ALLOCATE(ZMASS(ISVECNPT,JPMODE))
ALLOCATE(ZCCTOT(ISVECNPT,NSP+NCARB+NSOA,JPMODE))
ALLOCATE(ZCTOTG(ISVECNPT,NSP+NCARB+NSOA))
ALLOCATE(ZMU(ISVECNPT))
ALLOCATE(ZLAMBDA(ISVECNPT))
ALLOCATE(ZOM(ISVECNPT,JPMODE))
ALLOCATE(ZSO4RAT(ISVECNPT))
ALLOCATE(ZRV(ISVECNPT))
ALLOCATE(ZRC(ISVECNPT))
ALLOCATE(ZPRESSURE(ISVECNPT))
ALLOCATE(ZTEMP(ISVECNPT))
ALLOCATE(ZDENAIR(ISVECNPT))
ALLOCATE(ZFRAC(ISVECNPT,NEQ))
ALLOCATE(ZMI(ISVECNPT,NSP+NCARB+NSOA))
ALLOCATE(ZSOLORG(ISVECNPT,NSP+NCARB+NSOA))
ZSEDA(:,:)=0.
ZSOLORG(:,:) = 0.
! initialisation of aerosol molecular weight
ZMI(:,:) = 250.
ZMI(:,JP_AER_SO4)  = 98.
ZMI(:,JP_AER_NO3)  = 63.
ZMI(:,JP_AER_NH3)  = 17.
ZMI(:,JP_AER_H2O)  = 18.
IF (NSOA .EQ. 10) THEN
  ZMI(:,JP_AER_SOA1) = 88. 
  ZMI(:,JP_AER_SOA2) = 180.
  ZMI(:,JP_AER_SOA3) = 1.5374857E3
  ZMI(:,JP_AER_SOA4) = 1.9586780E3
  ZMI(:,JP_AER_SOA5) = 195.
  ZMI(:,JP_AER_SOA6) = 195.
  ZMI(:,JP_AER_SOA7) = 165.
  ZMI(:,JP_AER_SOA8) = 195.
  ZMI(:,JP_AER_SOA9) = 270.
  ZMI(:,JP_AER_SOA10) = 210.
END IF
END IF


!
!-------------------------------------------------------------------------------

ZLAT(:,:) = PLAT(:,:)*(XPI/180.)
ZLON(:,:) = PLON(:,:)*(XPI/180.)
!
ZTIMERAD     = ZTIME + 0.5* PTSTEP

ZUT       = MOD( 24.0+MOD(ZTIMERAD/3600.,24.0),24.0 )

INOBIS(:) = (/0,31,59,90,120,151,181,212,243,273,304,334/)
IBIS(0) = INOBIS(0)
DO JI=1,11
  IBIS(JI) = INOBIS(JI)+1
END DO
IF( MOD(IYEAR,4).EQ.0 ) THEN
  ZDATE = FLOAT(IDAY +   IBIS(IMONTH-1)) - 1
  ZAD = 2.0*XPI*ZDATE/366.0
ELSE
  ZDATE = FLOAT(IDAY + INOBIS(IMONTH-1)) - 1
  ZAD = 2.0*XPI*ZDATE/365.0
END IF
ZDECSOL = 0.006918-0.399912*COS(ZAD)   +0.070257*SIN(ZAD)    &
         -0.006758*COS(2.*ZAD)+0.000907*SIN(2.*ZAD) &
         -0.002697*COS(3.*ZAD)+0.00148 *SIN(3.*ZAD)
ZSINDEL = SIN(ZDECSOL)
ZCOSDEL = COS(ZDECSOL)

ZA1 = (1.00554*ZDATE- 6.28306)*(XPI/180.0)
ZA2 = (1.93946*ZDATE+23.35089)*(XPI/180.0)
ZTSIDER = (7.67825*SIN(ZA1)+10.09176*SIN(ZA2)) / 60.0

ZTUT(:,:) = ZUT - ZTSIDER + ZLON(:,:)*((180./XPI)/15.0)
ZSOLANG(:,:) = (ZTUT(:,:)-12.0)*15.0*(XPI/180.)          ! hour angle in radians
ZCOSZEN(:,:) = SIN(ZLAT(:,:))*ZSINDEL +                 &! Cosine of the zenithal
               COS(ZLAT(:,:))*ZCOSDEL*COS(ZSOLANG(:,:))  !       solar angle
ZZENITH(:,:) = ACOS(ZCOSZEN(:,:))
!
ILAT = 1

  CALL CH_UPDATE_JVALUES_n(KLUOUT, ZZENITH, ZRT,  &
    &   ZALB_UV, ZZS, ZZZ, XLAT0, XLON0,          &
    &   KLON, ILAT, KLEV, KRR,                  &
    &   IDAY, IMONTH, IYEAR, ZTIME,               &
    &   LCH_TUV_ONLINE, CCH_TUV_CLOUDS,           &
    &   XCH_TUV_ALBNEW, XCH_TUV_DOBNEW, ZRHODREF, ZJVALUES,&
    &   IIB,IIE,IJB,IJE,IIU,IJU,10)

!
!*       3.    INTEGRATE OVER ALL GRID POINTS
!              -------------------------------
!
!
DO JL=1,ISVECNMASK
                     
!
!*       3.1   transfer chemical species for solver
!              and convert from part/part to molec./cm3
  IDTI=ISVECMASK(2,JL)-ISVECMASK(1,JL)+1
  IDTJ=ISVECMASK(4,JL)-ISVECMASK(3,JL)+1
  IDTK=ISVECMASK(6,JL)-ISVECMASK(5,JL)+1
!
ZDEN2MOL = 1E-6 * XAVOGADRO / XMD


DO JN = 1, NSV_CHEM
!Vectorization:
!ocl novrec
!cdir nodep
  DO JM=0,ISVECNPT-1
    JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
    JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
    ZCHEM(JM+1,JN) = PRSVS(JI,1,JK,NSV_CHEMBEG+JN-1) * PTSTEP * &
                                   ZDEN2MOL * ZRHODREF(JI,1,JK) 
  END DO
END DO


ZCHEM(:,:) = MAX(ZCHEM(:,:), 0.)

IF (LORILAM) THEN
ZRV(:) = 0.
ZRC(:) = 0.
!ocl novrec
!cdir nodep
DO JM=0,ISVECNPT-1
    JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
    JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
! Sedimentation tendancy
    ZSEDA(JM+1,:) = ZXSEDA(JI,1,JK,:)
! Pressure (Pa)
    ZPRESSURE(JM+1) = PABST(JI,1,JK)
! Air density (kg/m3)
    ZDENAIR(JM+1) = PRHODREF(JI,1,JK)
! Temperature (K)
    ZTEMP(JM+1)   = PTHT(JI,1,JK)*((PABST(JI,1,JK)/XP00)**(XRD/XCPD))
! Water vapor (kg/kg)
    ZRV(JM+1)    = PRT(JI,1,JK,1)
! Cloud vapor (kg/kg)
    ZRC(JM+1)  = PRT(JI,1,JK,2)
!  END DO
!
END DO
DO JN = 1, NSV_AER
!Vectorization:
!ocl novrec
  DO JM=0,ISVECNPT-1
    JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
    JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
    ZAERO(JM+1,JN) = PRSVS(JI,1,JK,NSV_AERBEG+JN-1) * PTSTEP * &
                                  ZDEN2MOL * ZRHODREF(JI,1,JK)
  END DO
END DO

 IF (KTCOUNT == 1) THEN
!*       3.4 initialize aerosol parameters and moments of 0th,
!            6th, aerosol surface and aerosol diameter order

          CALL CH_INI_ORILAM(ZM, ZSIG0, ZRG0, ZN0, ZCTOTG, ZCTOTA, ZCCTOT,&
                               ZSEDA, ZOM, ZRHOP0, ZAERO, ZCHEM, ZRV, ZDENAIR, &
                               ZPRESSURE, ZTEMP, ZRC, ZFRAC, ZMI, CCH_SCHEME)
!          
 END IF
ZCHEM(:,JP_CH_H2SO4) =  ZAERO(:,JP_CH_SO4i) + ZAERO(:,JP_CH_SO4j)
END IF


!
!*       3.2   transfer meteo data into chemical core system
!
  CALL CH_METEO_TRANS(JL, ZRHODREF, ZRT, ZTHT, ZABST, &
                     &   ISVECNPT, ISVECMASK, TZM, IDAY,&
                     &   IMONTH, IYEAR, PLAT, PLON, XLAT0, XLON0,&
                     &   .TRUE., .TRUE.,  KLUOUT)

!*       3.3   calculate reaction and photolysis rates
!
  ZCHEM(:,:) = MAX(ZCHEM(:,:), 1E-40)
  CALL CH_SET_RATES (ZTIME, ZCHEM, TZM, 1, KLUOUT, NVERB, ISVECNPT, NEQ)

  CALL CH_SET_PHOTO_RATES &
    (ZTIME, ZCHEM, JL, TZM, 1, KLUOUT, NVERB, ISVECNPT, ISVECMASK, NEQ, ZJVALUES)
!
!
!*       3.4   solve chemical system for the timestep of the monitor
!
  ZNEWCHEM(:,:) = 0.
  ZOLDCHEM(:,:) =  ZCHEM(:,:)

  DO JM = 1, NCH_SUBSTEPS
    CALL CH_SOLVER_n &
          (ZTIME, ZDTSOLVER, ZCHEM, ZNEWCHEM, NEQ, ISVECNPT, 1)
    ZCHEM(:,:) = MAX(ZNEWCHEM(:,:), 1E-40)
  END DO

IF (LORILAM) THEN
!
  ZSO4RAT(:) = (ZNEWCHEM(:,JP_CH_H2SO4) - ZOLDCHEM(:,JP_CH_H2SO4) ) / PTSTEP
!
!*       3.6   solve aerosol system 
  CALL CH_ORILAM(ZAERO,ZCHEM, ZM, ZSIG0, ZRG0, ZN0,  ZCTOTG, ZCTOTA, ZCCTOT, &
                 PTSTEP, ZSEDA, ZMU, ZLAMBDA, ZRHOP0, ZOM, ZSO4RAT,      &
                 ZRV, ZDENAIR,ZPRESSURE, ZTEMP, ZRC, ZFRAC, ZMI, ZTIME,      &
                 CCH_SCHEME,ZSOLORG)
!
!*       3.7   return result to MesoNH scalar variables
!
DO JN = 1, NSV_AER
!Vectorization:
!ocl novrec
!cdir nodep
  DO JM=0,ISVECNPT-1
      JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
      JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
      PRSVS(JI,1,JK,NSV_AERBEG+JN-1) = ZAERO(JM+1,JN) / &
                          (PTSTEP * ZDEN2MOL * ZRHODREF(JI,1,JK))
   END DO
 END DO
IF (NDIAG .GE. JPMODE*2) THEN
ZPEZDIAG(:,:,:) = 0.
ZMASS(:,:) = 0.
DO JN = 1, JPMODE
  DO JSV=1,NSP+NCARB+NSOA
     ZMASS(:,JN) = ZMASS(:,JN)+ZCTOTA(:,JSV,JN)
  ENDDO
!Vectorization:
!ocl novrec
!cdir nodep
  DO JM=0,ISVECNPT-1
    JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
    JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
    ZPEZDIAG(JI,JK,JN) = ZN0(JM+1,JN) 
    IF (NDIAG .GT. (JPMODE+1)) ZPEZDIAG(JI,JK,JPMODE+JN) = ZRG0(JM+1,JN) 
    IF (NDIAG .GT. (JPMODE*2+1)) ZPEZDIAG(JI,JK,JPMODE*2+JN) = ZMASS(JM+1,JN) 
  END DO
END DO
PPEZDIAG(:,:,:) = 0.
DO JN = 1, NDIAG
  DO JK=1,KLEV
!    JK1=KLEV-JK
    PPEZDIAG(:,JK,JN) = ZPEZDIAG(:,JK,JN)
  ENDDO
ENDDO
END IF
!
END IF
!

!*       3.5   return result to Arome scalar variables
!
DO JN = 1, NSV_CHEM
!Vectorization:
!ocl novrec
!cdir nodep
  DO JM=0,ISVECNPT-1
    JI=JM-IDTI*(JM/IDTI)+ISVECMASK(1,JL)
    JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+ISVECMASK(5,JL)
    PRSVS(JI,1,JK,NSV_CHEMBEG+JN-1) = ZCHEM(JM+1,JN) / &
                          (PTSTEP * ZDEN2MOL * ZRHODREF(JI,1,JK))
  END DO
END DO

END DO
!
!
!*       4.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!   no mass correction (mnh surcouche not included in arome)
PRSVS(:,:,:,NSV_CHEMBEG:NSV_CHEMEND) = MAX(0., PRSVS(:,:,:,NSV_CHEMBEG:NSV_CHEMEND))
PRSVS(:,:,:,NSV_AERBEG:NSV_AEREND) = MAX(0., PRSVS(:,:,:,NSV_AERBEG:NSV_AEREND))
!
!----------------------------------------------------------------------
!
!
DEALLOCATE(TZM)
DEALLOCATE(ZCHEM)
DEALLOCATE(ZNEWCHEM)
DEALLOCATE(ZOLDCHEM)
IF (LORILAM) THEN
  DEALLOCATE(ZAERO)
  DEALLOCATE(ZM)
  DEALLOCATE(ZSEDA)
  DEALLOCATE(ZN0)
  DEALLOCATE(ZRG0)
  DEALLOCATE(ZSIG0)
  DEALLOCATE(ZRHOP0)
  DEALLOCATE(ZCTOTA)
  DEALLOCATE(ZMASS)
  DEALLOCATE(ZCCTOT)
  DEALLOCATE(ZCTOTG)
  DEALLOCATE(ZMU)
  DEALLOCATE(ZLAMBDA)
  DEALLOCATE(ZOM)
  DEALLOCATE(ZSO4RAT)
  DEALLOCATE(ZRV)
  DEALLOCATE(ZRC)
  DEALLOCATE(ZPRESSURE)
  DEALLOCATE(ZTEMP)
  DEALLOCATE(ZDENAIR)
  DEALLOCATE(ZFRAC)
  DEALLOCATE(ZMI)
  DEALLOCATE(ZSOLORG)
END IF

!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_MNHC',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_MNHC
