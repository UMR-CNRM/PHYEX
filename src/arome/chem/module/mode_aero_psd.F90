!     ######spl
MODULE MODE_AERO_PSD
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ########################
!!
!! MODULE DUST PSD (Particle Size Distribution)
!! Purpose: Contains subroutines to convert from transported variables (ppp)
!! to understandable aerosol variables, e.g. #/m3, kg/m3, sigma, R_{n}

USE MODD_CH_AEROSOL
!
IMPLICIT NONE
! 
CONTAINS
! 
!    ############################################################
  SUBROUTINE PPP2AERO(             &
       PSVT                         & !I [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PSIG3D                     & !O [-] standard deviation of aerosol distribution
       , PRG3D                      & !O [um] number median diameter of aerosol distribution
       , PN3D                       & !O [#/m3] number concentration of aerosols
       , PCTOTA                     & !O [ug/m3] mass of each aerosol compounds
       , PM3D                       & !moments 0, 3 and 6
       )
!    ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the three moments M0, M3 and M6 given in ppp into
!!    Values which can be understood more easily (R, sigma, N, M)
!! 
!!    CALLING STRUCTURE NOTE: OPTIONAL VARIABLES
!!    -------
!!    CALL PPP2AEROS(PSVT, PRHODREF, PSIG3D=SIGVAR,  &
!!       PRG3D=RVAR, PN3D=NVAR, PM3D=ZM)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Alf Grini (CNRM)
!!
!!    EXTERNAL
!!    --------
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,    DIMENSION(:,:,:,:),  INTENT(IN)     :: PSVT     !I [#/molec_{air}] first moment
                                                             !I [molec_{aer}/molec_{air} 3rd moment
                                                             !I [um6/molec_{air}*(cm3/m3)] 6th moment
REAL,    DIMENSION(:,:,:),    INTENT(IN)     :: PRHODREF !I [kg/m3] density of air

REAL,    DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)   :: PSIG3D   !O [-] standard deviation
REAL,    DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)   :: PRG3D    !O [um] number median diameter
REAL,    DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)   :: PN3D     !O [#/m3] number concentration
REAL,    DIMENSION(:,:,:,:,:),OPTIONAL, INTENT(OUT)   :: PCTOTA   !O [ug/m3] mass of each component
REAL,    DIMENSION(:,:,:,:),  OPTIONAL, INTENT(OUT)   :: PM3D     !O moments 0,3 and 6 
!
!*      0.2    declarations local variables
!
REAL, DIMENSION(NSP+NCARB+NSOA)      :: ZMI                 ! [kg/mol] molar weight of aerosol
REAL                                 :: ZRGMIN              ! [um] minimum radius accepted
REAL,DIMENSION(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3),SIZE(PSVT,4)) :: ZSV ! [aerosol concentration]
REAL,DIMENSION(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3)) :: ZSIGMA   ! [-] standard deviation
REAL,DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NSP+NCARB+NSOA,JPMODE):: ZCTOTA
REAL,DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),JPMODE*3) :: ZM

REAL,DIMENSION(NSP+NCARB+NSOA)       :: ZFAC                ! M3 / mass conversion factor
REAL                                 :: ZDEN2MOL
REAL,DIMENSION(JPMODE*3)             :: ZPMIN               ! [aerosol units] minimum values for N, sigma, M
INTEGER                              :: JJ, JN              ! [idx] loop counters
REAL    :: ZINIRADIUSI, ZINIRADIUSJ
!
!-------------------------------------------------------------------------------
!
!        1.    initialisation 

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_AERO_PSD:PPP2AERO',0,ZHOOK_HANDLE)
IF (CRGUNIT=="MASS") THEN
  ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
  ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
  ZINIRADIUSI = XINIRADIUSI 
  ZINIRADIUSJ = XINIRADIUSJ
END IF



    !Get minimum values possible
    ZPMIN(1) = XN0IMIN
    ZRGMIN = XCOEFRADIMIN * ZINIRADIUSI
    ZPMIN(2) = ZPMIN(1) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGIMIN)**2) 
    ZPMIN(3) = ZPMIN(1) * (ZRGMIN**6)*EXP(18. * LOG(XSIGIMIN)**2)

    ZPMIN(4) = XN0JMIN
    ZRGMIN = XCOEFRADJMIN * ZINIRADIUSJ
    ZPMIN(5) = ZPMIN(4) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGJMIN)**2) 
    ZPMIN(6) = ZPMIN(4) * (ZRGMIN**6)*EXP(18. * LOG(XSIGJMIN)**2)

    !Set density of aerosol, here dust (kg/m3)
    ! Aerosol Density
! Cf Ackermann (all to black carbon except water)

    !Set molecular weightn g/mol 
    ZMI(:) = 250.
    ZMI(JP_AER_SO4)  = 98.
    ZMI(JP_AER_NO3)  = 63.
    ZMI(JP_AER_NH3)  = 17.
    ZMI(JP_AER_H2O)  = 18.
    IF (NSOA .EQ. 10) THEN
    ZMI(JP_AER_SOA1) = 88. 
    ZMI(JP_AER_SOA2) = 180.
    ZMI(JP_AER_SOA3) = 1.5374857E2
    ZMI(JP_AER_SOA4) = 1.9586780E2
    ZMI(JP_AER_SOA5) = 195.
    ZMI(JP_AER_SOA6) = 195.
    ZMI(JP_AER_SOA7) = 165.
    ZMI(JP_AER_SOA8) = 195.
    ZMI(JP_AER_SOA9) = 270.
    ZMI(JP_AER_SOA10) = 210.
    END IF

! conversion into mol.cm-3
ZDEN2MOL = 1E-6 * XAVOGADRO  / XMD
!
DO JJ=1, SIZE(PSVT,4)
  ZSV(:,:,:,JJ) =  PSVT(:,:,:,JJ) * ZDEN2MOL * PRHODREF(:,:,:)
ENDDO
!
DO JJ=1,NSP+NCARB+NSOA
  ZFAC(JJ)=(4./3.)*XPI*XRHOI(JJ)*1.e-9
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       2    transfer aerosol mass from gas to aerosol variables
!               (and conversion of mol.cm-3 --> microgram/m3)
!
ZCTOTA(:,:,:,:,:) = 0.
! aerosol phase
  ZCTOTA(:,:,:,JP_AER_SO4,1) = ZSV(:,:,:,JP_CH_SO4i)*ZMI(JP_AER_SO4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SO4,2) = ZSV(:,:,:,JP_CH_SO4j)*ZMI(JP_AER_SO4)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_NO3,1) = ZSV(:,:,:,JP_CH_NO3i)*ZMI(JP_AER_NO3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_NO3,2) = ZSV(:,:,:,JP_CH_NO3j)*ZMI(JP_AER_NO3)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_NH3,1) = ZSV(:,:,:,JP_CH_NH3i)*ZMI(JP_AER_NH3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_NH3,2) = ZSV(:,:,:,JP_CH_NH3j)*ZMI(JP_AER_NH3)/6.0221367E+11
!
! water
  ZCTOTA(:,:,:,JP_AER_H2O,1) = ZSV(:,:,:,JP_CH_H2Oi)*ZMI(JP_AER_H2O)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_H2O,2) = ZSV(:,:,:,JP_CH_H2Oj)*ZMI(JP_AER_H2O)/6.0221367E+11
!
! primary organic carbon
  ZCTOTA(:,:,:,JP_AER_OC,1) = ZSV(:,:,:,JP_CH_OCi)*ZMI(JP_AER_OC)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_OC,2) = ZSV(:,:,:,JP_CH_OCj)*ZMI(JP_AER_OC)/6.0221367E+11
!
! primary black carbon
  ZCTOTA(:,:,:,JP_AER_BC,1) = ZSV(:,:,:,JP_CH_BCi)*ZMI(JP_AER_BC)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_BC,2) = ZSV(:,:,:,JP_CH_BCj)*ZMI(JP_AER_BC)/6.0221367E+11
!
IF (NSOA .EQ. 10) THEN
  ZCTOTA(:,:,:,JP_AER_SOA1,1) = ZSV(:,:,:,JP_CH_SOA1i)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA1,2) = ZSV(:,:,:,JP_CH_SOA1j)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA2,1) = ZSV(:,:,:,JP_CH_SOA2i)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA2,2) = ZSV(:,:,:,JP_CH_SOA2j)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA3,1) = ZSV(:,:,:,JP_CH_SOA3i)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA3,2) = ZSV(:,:,:,JP_CH_SOA3j)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA4,1) = ZSV(:,:,:,JP_CH_SOA4i)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA4,2) = ZSV(:,:,:,JP_CH_SOA4j)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA5,1) = ZSV(:,:,:,JP_CH_SOA5i)*ZMI(JP_AER_SOA5)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA5,2) = ZSV(:,:,:,JP_CH_SOA5j)*ZMI(JP_AER_SOA5)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_SOA6,1) = ZSV(:,:,:,JP_CH_SOA6i)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA6,2) = ZSV(:,:,:,JP_CH_SOA6j)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA7,1) = ZSV(:,:,:,JP_CH_SOA7i)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA7,2) = ZSV(:,:,:,JP_CH_SOA7j)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA8,1) = ZSV(:,:,:,JP_CH_SOA8i)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA8,2) = ZSV(:,:,:,JP_CH_SOA8j)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA9,1) = ZSV(:,:,:,JP_CH_SOA9i)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA9,2) = ZSV(:,:,:,JP_CH_SOA9j)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA10,1) = ZSV(:,:,:,JP_CH_SOA10i)*ZMI(JP_AER_SOA10)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA10,2) = ZSV(:,:,:,JP_CH_SOA10j)*ZMI(JP_AER_SOA10)/6.0221367E+11
END IF
!
!-------------------------------------------------------------------------------
!
!*       3    calculate moment 3 from total aerosol mass
!
ZM(:,:,:,2) = 0.
ZM(:,:,:,5) = 0.
DO JJ = 1,NSP+NCARB+NSOA
    ZM(:,:,:,2) = ZM(:,:,:,2)+ZCTOTA(:,:,:,JJ,1)/ZFAC(JJ) !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
    ZM(:,:,:,5) = ZM(:,:,:,5)+ZCTOTA(:,:,:,JJ,2)/ZFAC(JJ) !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       4    set  moment 0 
!
   ZM(:,:,:,1)=   MAX(ZSV(:,:,:,JP_CH_M0i) * 1E+6, 1.E-10) ! molec_{aer}/m3_{air}
   ZM(:,:,:,4)=   MAX(ZSV(:,:,:,JP_CH_M0j) * 1E+6, 1.E-10) ! molec_{aer}/m3_{air}
!
!-------------------------------------------------------------------------------
!
!*       5    set moment 6  ==> um6_{aer}/m3_{air}
!
IF (LVARSIGI) THEN ! set M6 variable standard deviation
  ZM(:,:,:,3) = MAX(ZSV(:,:,:,JP_CH_M6i), 1E-80)

  ZSIGMA(:,:,:)=ZM(:,:,:,2)**2/(ZM(:,:,:,1)*ZM(:,:,:,3))
  ZSIGMA(:,:,:)=MIN(1-1E-10,ZSIGMA(:,:,:))
  ZSIGMA(:,:,:)=MAX(1E-10,ZSIGMA(:,:,:))
  ZSIGMA(:,:,:)= LOG(ZSIGMA(:,:,:))
  ZSIGMA(:,:,:)= EXP(1./3.*SQRT(-ZSIGMA(:,:,:)))
  WHERE (ZSIGMA(:,:,:) > XSIGIMAX)
   ZSIGMA(:,:,:) =  XSIGIMAX
  END WHERE
  WHERE (ZSIGMA(:,:,:) < XSIGIMIN)
   ZSIGMA(:,:,:) =  XSIGIMIN
  END WHERE
  ZM(:,:,:,3) = ZM(:,:,:,1) &
          * ( (ZM(:,:,:,2)/ZM(:,:,:,1))**(1./3.)  &
          * exp(-(3./2.)*log(ZSIGMA(:,:,:))**2))**6 &
          * exp(18.*log(ZSIGMA(:,:,:))**2)

 IF(PRESENT(PSIG3D)) PSIG3D(:,:,:,1) = ZSIGMA(:,:,:)

ELSE ! fixed standard deviation
 ZM(:,:,:,3) = ZM(:,:,:,1) &
          * ( (ZM(:,:,:,2)/ZM(:,:,:,1))**(1./3.)  &
          * exp(-(3./2.)*log(XINISIGI)**2))**6 &
          * exp(18.*log(XINISIGI)**2)

 IF(PRESENT(PSIG3D)) PSIG3D(:,:,:,1) = XINISIGI
END IF

IF (LVARSIGJ) THEN ! set M6 variable standard deviation
  ZM(:,:,:,6) = MAX(ZSV(:,:,:,JP_CH_M6j), 1E-80)

  ZSIGMA(:,:,:)=ZM(:,:,:,5)**2/(ZM(:,:,:,4)*ZM(:,:,:,6))
  ZSIGMA(:,:,:)=MIN(1-1E-10,ZSIGMA(:,:,:))
  ZSIGMA(:,:,:)=MAX(1E-10,ZSIGMA(:,:,:))
  ZSIGMA(:,:,:)= LOG(ZSIGMA(:,:,:))
  ZSIGMA(:,:,:)= EXP(1./3.*SQRT(-ZSIGMA(:,:,:)))
  WHERE (ZSIGMA(:,:,:) > XSIGJMAX)
    ZSIGMA(:,:,:) =  XSIGJMAX
  END WHERE
  WHERE (ZSIGMA(:,:,:) < XSIGJMIN)
    ZSIGMA(:,:,:) =  XSIGJMIN
  END WHERE

  ZM(:,:,:,6) = ZM(:,:,:,4) &
          * ( (ZM(:,:,:,5)/ZM(:,:,:,4))**(1./3.)  &
          * exp(-(3./2.)*log(ZSIGMA(:,:,:))**2))**6 &
          * exp(18.*log(ZSIGMA(:,:,:))**2)

  IF(PRESENT(PSIG3D)) PSIG3D(:,:,:,2) = ZSIGMA(:,:,:)

ELSE ! fixed standard deviation
  ZM(:,:,:,6) = ZM(:,:,:,4) &
          * ( (ZM(:,:,:,5)/ZM(:,:,:,4))**(1./3.)  &
          * exp(-(3./2.)*log(XINISIGJ)**2))**6 &
          * exp(18.*log(XINISIGJ)**2)

  IF(PRESENT(PSIG3D)) PSIG3D(:,:,:,2) = XINISIGJ
END IF


!-------------------------------------------------------------------------------
!
!*       6    calculate modal parameters from moments
!
DO JN=1,JPMODE
  IF(PRESENT(PN3D)) PN3D(:,:,:,JN) = ZM(:,:,:,NM0(JN))

  IF(PRESENT(PRG3D)) PRG3D(:,:,:,JN)=(ZM(:,:,:,NM3(JN))**4. &
              / (ZM(:,:,:,NM6(JN))*ZM(:,:,:,NM0(JN))**3.))**(1./6.)

ENDDO
!
IF(PRESENT(PCTOTA)) PCTOTA(:,:,:,:,:) = ZCTOTA(:,:,:,:,:)
IF(PRESENT(PM3D)) PM3D(:,:,:,:) = ZM(:,:,:,:)
!
!
IF (LHOOK) CALL DR_HOOK('MODE_AERO_PSD:PPP2AERO',1,ZHOOK_HANDLE)
END SUBROUTINE PPP2AERO
!
!   ############################################################
  SUBROUTINE AERO2PPP(             &
       PSVT                         & !IO [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PSIG3D                     & !I [-] standard deviation of aerosol distribution
       , PRG3D                      & !I [um] number median diameter of aerosol distribution
       )
!!   ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the aerosol Mass, RG and SIGMA in the  three moments M0, M3 and M6 given in ppp 
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Alf Grini (CNRM)
!!
!!    EXTERNAL
!!    --------
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,       DIMENSION(:,:,:,:),  INTENT(INOUT)     :: PSVT      !I [#/molec_{air}] first moment
                                                                !I [molec_{aer}/molec_{air} 3rd moment
                                                                !I [um6/molec_{air}*(cm3/m3)] 6th moment
REAL,       DIMENSION(:,:,:),  INTENT(IN)      :: PRHODREF !I [kg/m3] density of air

REAL,       DIMENSION(:,:,:,:),  INTENT(IN)     :: PSIG3D   !O [-] standard deviation
REAL,       DIMENSION(:,:,:,:),  INTENT(IN)     :: PRG3D    !O [um] number median diameter
!
!
!*      0.2    declarations local variables
!
REAL,DIMENSION(NSP+NCARB+NSOA)       :: ZMI                 ! [kg/mol] molar weight of aerosol
REAL                                 :: ZRGMIN              ! [um] minimum radius accepted
REAL,DIMENSION(:,:,:,:), ALLOCATABLE :: ZM                  ! [aerosol units] local array which goes to output later
REAL,DIMENSION(JPMODE*3)             :: ZPMIN               ! [aerosol units] minimum values for N, sigma, M
REAL,DIMENSION(SIZE(PSVT,1),SIZE(PSVT,2),SIZE(PSVT,3),NSP+NCARB+NSOA,JPMODE):: ZCTOTA

REAL,DIMENSION(NSP+NCARB+NSOA)       :: ZFAC                ! M3 / mass conversion factor
INTEGER                              :: JJ                  ! [idx] loop counters
REAL                                 :: ZDEN2MOL
REAL    :: ZINIRADIUSI, ZINIRADIUSJ
!
!-------------------------------------------------------------------------------
!
!        1.    initialisation 

!Calculations here are for one mode only
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_AERO_PSD:AERO2PPP',0,ZHOOK_HANDLE)
IF (CRGUNIT=="MASS") THEN
    ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
    ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
    ZINIRADIUSI = XINIRADIUSI 
    ZINIRADIUSJ = XINIRADIUSJ
END IF


!Get minimum values possible
ZPMIN(1) = XN0IMIN 
ZRGMIN = XCOEFRADIMIN * ZINIRADIUSI
ZPMIN(2) = ZPMIN(1) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGIMIN)**2) 
ZPMIN(3) = ZPMIN(1) * (ZRGMIN**6)*EXP(18. * LOG(XSIGIMIN)**2)
ZPMIN(4) = XN0JMIN 
ZRGMIN = XCOEFRADJMIN * ZINIRADIUSJ
ZPMIN(5) = ZPMIN(4) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGJMIN)**2) 
ZPMIN(6) = ZPMIN(4) * (ZRGMIN**6)*EXP(18. * LOG(XSIGJMIN)**2)

ALLOCATE (ZM(SIZE(PSVT,1), SIZE(PSVT,2), SIZE(PSVT,3), JPMODE*3))
    
    !Set density of aerosol, here dust (kg/m3)
    ! Aerosol Density
! Cf Ackermann (all to black carbon except water)

    !Set molecular weightn g/mol 
    ZMI(:) = 250.
    ZMI(JP_AER_SO4)  = 98.
    ZMI(JP_AER_NO3)  = 63.
    ZMI(JP_AER_NH3)  = 17.
    ZMI(JP_AER_H2O)  = 18.
    IF (NSOA .EQ. 10) THEN
    ZMI(JP_AER_SOA1) = 88. 
    ZMI(JP_AER_SOA2) = 180.
    ZMI(JP_AER_SOA3) = 1.5374857E2
    ZMI(JP_AER_SOA4) = 1.9586780E2
    ZMI(JP_AER_SOA5) = 195.
    ZMI(JP_AER_SOA6) = 195.
    ZMI(JP_AER_SOA7) = 165.
    ZMI(JP_AER_SOA8) = 195.
    ZMI(JP_AER_SOA9) = 270.
    ZMI(JP_AER_SOA10) = 210.
    END IF

! conversion into mol.cm-3
ZDEN2MOL = 1E-6 * XAVOGADRO  / XMD
DO JJ=1, SIZE(PSVT, 4)
  PSVT(:,:,:,JJ) =  PSVT(:,:,:,JJ) * ZDEN2MOL * PRHODREF(:,:,:)
ENDDO
!
DO JJ=1,NSP+NCARB+NSOA
  ZFAC(JJ)=(4./3.)*XPI*XRHOI(JJ)*1.e-9
ENDDO
!
!
!*       2    transfer aerosol mass from gas to aerosol variables
!               (and conversion of mol.cm-3 --> microgram/m3)
!
ZCTOTA(:,:,:,:,:) = 0.
! aerosol phase
  ZCTOTA(:,:,:,JP_AER_SO4,1) = PSVT(:,:,:,JP_CH_SO4i)*ZMI(JP_AER_SO4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SO4,2) = PSVT(:,:,:,JP_CH_SO4j)*ZMI(JP_AER_SO4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_NO3,1) = PSVT(:,:,:,JP_CH_NO3i)*ZMI(JP_AER_NO3)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_NO3,2) = PSVT(:,:,:,JP_CH_NO3j)*ZMI(JP_AER_NO3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_NH3,1) = PSVT(:,:,:,JP_CH_NH3i)*ZMI(JP_AER_NH3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_NH3,2) = PSVT(:,:,:,JP_CH_NH3j)*ZMI(JP_AER_NH3)/6.0221367E+11
!
! water
  ZCTOTA(:,:,:,JP_AER_H2O,1) = PSVT(:,:,:,JP_CH_H2Oi)*ZMI(JP_AER_H2O)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_H2O,2) = PSVT(:,:,:,JP_CH_H2Oj)*ZMI(JP_AER_H2O)/6.0221367E+11
!
! primary organic carbon
  ZCTOTA(:,:,:,JP_AER_OC,1) = PSVT(:,:,:,JP_CH_OCi)*ZMI(JP_AER_OC)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_OC,2) = PSVT(:,:,:,JP_CH_OCj)*ZMI(JP_AER_OC)/6.0221367E+11
!
! primary black carbon
  ZCTOTA(:,:,:,JP_AER_BC,1) = PSVT(:,:,:,JP_CH_BCi)*ZMI(JP_AER_BC)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_BC,2) = PSVT(:,:,:,JP_CH_BCj)*ZMI(JP_AER_BC)/6.0221367E+11
!
 IF (NSOA .EQ. 10) THEN
  ZCTOTA(:,:,:,JP_AER_SOA1,1) = PSVT(:,:,:,JP_CH_SOA1i)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA1,2) = PSVT(:,:,:,JP_CH_SOA1j)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA2,1) = PSVT(:,:,:,JP_CH_SOA2i)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA2,2) = PSVT(:,:,:,JP_CH_SOA2j)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA3,1) = PSVT(:,:,:,JP_CH_SOA3i)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA3,2) = PSVT(:,:,:,JP_CH_SOA3j)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA4,1) = PSVT(:,:,:,JP_CH_SOA4i)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA4,2) = PSVT(:,:,:,JP_CH_SOA4j)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA5,1) = PSVT(:,:,:,JP_CH_SOA5i)*ZMI(JP_AER_SOA5)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA5,2) = PSVT(:,:,:,JP_CH_SOA5j)*ZMI(JP_AER_SOA5)/6.0221367E+11

  ZCTOTA(:,:,:,JP_AER_SOA6,1) = PSVT(:,:,:,JP_CH_SOA6i)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA6,2) = PSVT(:,:,:,JP_CH_SOA6j)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA7,1) = PSVT(:,:,:,JP_CH_SOA7i)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA7,2) = PSVT(:,:,:,JP_CH_SOA7j)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA8,1) = PSVT(:,:,:,JP_CH_SOA8i)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA8,2) = PSVT(:,:,:,JP_CH_SOA8j)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA9,1) = PSVT(:,:,:,JP_CH_SOA9i)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA9,2) = PSVT(:,:,:,JP_CH_SOA9j)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA10,1) = PSVT(:,:,:,JP_CH_SOA10i)*ZMI(JP_AER_SOA10)/6.0221367E+11
  ZCTOTA(:,:,:,JP_AER_SOA10,2) = PSVT(:,:,:,JP_CH_SOA10j)*ZMI(JP_AER_SOA10)/6.0221367E+11
  END IF

!
!*       3    calculate moment 3 from total aerosol mass
!
  ZM(:,:,:,2) = 0.
  ZM(:,:,:,5) = 0.
  DO JJ = 1,NSP+NCARB+NSOA
    ZM(:,:,:,2) = ZM(:,:,:,2)+ZCTOTA(:,:,:,JJ,1)/ZFAC(JJ)
    ZM(:,:,:,5) = ZM(:,:,:,5)+ZCTOTA(:,:,:,JJ,2)/ZFAC(JJ)
  ENDDO
  ZM(:,:,:,2) = MAX(ZM(:,:,:,2), ZPMIN(2))
  ZM(:,:,:,5) = MAX(ZM(:,:,:,5), ZPMIN(5))
!
!
!*       4    calculate moment 0 from dispersion and mean radius
!
   ZM(:,:,:,1)= ZM(:,:,:,2)/ &
              ( (PRG3D(:,:,:,1)**3)*EXP(4.5 * LOG(PSIG3D(:,:,:,1))**2) )
   ZM(:,:,:,4)= ZM(:,:,:,5)/ &
              ( (PRG3D(:,:,:,2)**3)*EXP(4.5 * LOG(PSIG3D(:,:,:,2))**2) )
!

!*       5    calculate moment 6 from dispersion and mean radius
!
 ZM(:,:,:,3) = ZM(:,:,:,1)*(PRG3D(:,:,:,1)**6) * &
               EXP(18 *(LOG(PSIG3D(:,:,:,1)))**2)
 ZM(:,:,:,6) = ZM(:,:,:,4)*(PRG3D(:,:,:,2)**6) * &
               EXP(18 *(LOG(PSIG3D(:,:,:,2)))**2)

!*       6    return to ppp
!
PSVT(:,:,:,JP_CH_M0i) = ZM(:,:,:,1) * 1E-6 
PSVT(:,:,:,JP_CH_M0j) = ZM(:,:,:,4) * 1E-6

IF (LVARSIGI) PSVT(:,:,:,JP_CH_M6i) = ZM(:,:,:,3) 
IF (LVARSIGJ) PSVT(:,:,:,JP_CH_M6j) = ZM(:,:,:,6)

DO JJ=1,SIZE(PSVT,4)
  PSVT(:,:,:,JJ) =  PSVT(:,:,:,JJ) / (ZDEN2MOL * PRHODREF(:,:,:))
ENDDO

DEALLOCATE(ZM)
!
IF (LHOOK) CALL DR_HOOK('MODE_AERO_PSD:AERO2PPP',1,ZHOOK_HANDLE)
END SUBROUTINE AERO2PPP
!
!    ############################################################
  SUBROUTINE PPP2AERO1D(            &
       PSVT                         & !I [ppp] input scalar variables (moment of distribution)
       , PRHODREF                   & !I [kg/m3] density of air       
       , PSIG1D                     & !O [-] standard deviation of aerosol distribution
       , PRG1D                      & !O [um] number median diameter of aerosol distribution
       , PN1D                       & !O [#/m3] number concentration of aerosols
       , PCTOTA                     & !O [ug/m3] mass of each aerosol compounds
       , PM1D                       & !moments 0, 3 and 6
       )
!    ############################################################
!
!!
!!    PURPOSE
!!    -------
!!    Translate the three moments M0, M3 and M6 given in ppp into
!!    Values which can be understood more easily (R, sigma, N, M)
!! 
!!    CALLING STRUCTURE NOTE: OPTIONAL VARIABLES
!!    -------
!!    CALL PPP2AERO1D(PSVT, PRHODREF, PSIG1D=SIGVAR,  &
!!       PRG1D=RVAR, PN1D=NVAR, PM1D=ZM)
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Pierre TULET (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Alf Grini (CNRM)
!!
!!    EXTERNAL
!!    --------
!!
    IMPLICIT NONE
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,    DIMENSION(:,:),  INTENT(IN)     :: PSVT     !I [#/molec_{air}] first moment
                                                             !I [molec_{aer}/molec_{air} 3rd moment
                                                             !I [um6/molec_{air}*(cm3/m3)] 6th moment
REAL,    DIMENSION(:),    INTENT(IN)     :: PRHODREF !I [kg/m3] density of air

REAL,    DIMENSION(:,:),  OPTIONAL, INTENT(OUT)   :: PSIG1D   !O [-] standard deviation
REAL,    DIMENSION(:,:),  OPTIONAL, INTENT(OUT)   :: PRG1D    !O [um] number median diameter
REAL,    DIMENSION(:,:),  OPTIONAL, INTENT(OUT)   :: PN1D     !O [#/m3] number concentration
REAL,    DIMENSION(:,:,:),OPTIONAL, INTENT(OUT)   :: PCTOTA   !O [ug/m3] mass of each component
REAL,    DIMENSION(:,:),  OPTIONAL, INTENT(OUT)   :: PM1D     !O moments 0,3 and 6 
!
!*      0.2    declarations local variables
!
REAL, DIMENSION(NSP+NCARB+NSOA)      :: ZMI                 ! [kg/mol] molar weight of aerosol
REAL                                 :: ZRGMIN              ! [um] minimum radius accepted
REAL,DIMENSION(SIZE(PSVT,1), SIZE(PSVT,2)) :: ZSV ! [aerosol concentration]
REAL,DIMENSION(SIZE(PSVT,1)) :: ZSIGMA   ! [-] standard deviation
REAL,DIMENSION(SIZE(PSVT,1),NSP+NCARB+NSOA,JPMODE):: ZCTOTA
REAL,DIMENSION(SIZE(PSVT,1),JPMODE*3) :: ZM

REAL,DIMENSION(NSP+NCARB+NSOA)       :: ZFAC                ! M3 / mass conversion factor
REAL                                 :: ZDEN2MOL
REAL,DIMENSION(JPMODE*3)             :: ZPMIN               ! [aerosol units] minimum values for N, sigma, M
INTEGER                              :: JJ, JN              ! [idx] loop counters
REAL    :: ZINIRADIUSI, ZINIRADIUSJ
!
!-------------------------------------------------------------------------------
!
!        1.    initialisation 

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODE_AERO_PSD:PPP2AERO1D',0,ZHOOK_HANDLE)
IF (CRGUNIT=="MASS") THEN
  ZINIRADIUSI = XINIRADIUSI * EXP(-3.*(LOG(XINISIGI))**2)
  ZINIRADIUSJ = XINIRADIUSJ * EXP(-3.*(LOG(XINISIGJ))**2)
ELSE
  ZINIRADIUSI = XINIRADIUSI 
  ZINIRADIUSJ = XINIRADIUSJ
END IF



    !Get minimum values possible
    ZPMIN(1) = XN0IMIN
    ZRGMIN = XCOEFRADIMIN * ZINIRADIUSI
    ZPMIN(2) = ZPMIN(1) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGIMIN)**2) 
    ZPMIN(3) = ZPMIN(1) * (ZRGMIN**6)*EXP(18. * LOG(XSIGIMIN)**2)

    ZPMIN(4) = XN0JMIN
    ZRGMIN = XCOEFRADJMIN * ZINIRADIUSJ
    ZPMIN(5) = ZPMIN(4) * (ZRGMIN**3)*EXP(4.5 * LOG(XSIGJMIN)**2) 
    ZPMIN(6) = ZPMIN(4) * (ZRGMIN**6)*EXP(18. * LOG(XSIGJMIN)**2)

    !Set density of aerosol, here dust (kg/m3)
    ! Aerosol Density
! Cf Ackermann (all to black carbon except water)

    !Set molecular weightn g/mol 
    ZMI(:) = 250.
    ZMI(JP_AER_SO4)  = 98.
    ZMI(JP_AER_NO3)  = 63.
    ZMI(JP_AER_NH3)  = 17.
    ZMI(JP_AER_H2O)  = 18.
    IF (NSOA .EQ. 10) THEN
    ZMI(JP_AER_SOA1) = 88. 
    ZMI(JP_AER_SOA2) = 180.
    ZMI(JP_AER_SOA3) = 1.5374857E3
    ZMI(JP_AER_SOA4) = 1.9586780E3
    ZMI(JP_AER_SOA5) = 195.
    ZMI(JP_AER_SOA6) = 195.
    ZMI(JP_AER_SOA7) = 165.
    ZMI(JP_AER_SOA8) = 195.
    ZMI(JP_AER_SOA9) = 270.
    ZMI(JP_AER_SOA10) = 210.
    END IF

! conversion into mol.cm-3
ZDEN2MOL = 1E-6 * XAVOGADRO  / XMD
!
DO JJ=1, SIZE(PSVT,2)
  ZSV(:,JJ) =  PSVT(:,JJ) * ZDEN2MOL * PRHODREF(:)
ENDDO
!
DO JJ=1,NSP+NCARB+NSOA
  ZFAC(JJ)=(4./3.)*XPI*XRHOI(JJ)*1.e-9
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       2    transfer aerosol mass from gas to aerosol variables
!               (and conversion of mol.cm-3 --> microgram/m3)
!
ZCTOTA(:,:,:) = 0.
! aerosol phase
  ZCTOTA(:,JP_AER_SO4,1) = ZSV(:,JP_CH_SO4i)*ZMI(JP_AER_SO4)/6.0221367E+11
  ZCTOTA(:,JP_AER_SO4,2) = ZSV(:,JP_CH_SO4j)*ZMI(JP_AER_SO4)/6.0221367E+11

  ZCTOTA(:,JP_AER_NO3,1) = ZSV(:,JP_CH_NO3i)*ZMI(JP_AER_NO3)/6.0221367E+11
  ZCTOTA(:,JP_AER_NO3,2) = ZSV(:,JP_CH_NO3j)*ZMI(JP_AER_NO3)/6.0221367E+11

  ZCTOTA(:,JP_AER_NH3,1) = ZSV(:,JP_CH_NH3i)*ZMI(JP_AER_NH3)/6.0221367E+11
  ZCTOTA(:,JP_AER_NH3,2) = ZSV(:,JP_CH_NH3j)*ZMI(JP_AER_NH3)/6.0221367E+11
!
! water
  ZCTOTA(:,JP_AER_H2O,1) = ZSV(:,JP_CH_H2Oi)*ZMI(JP_AER_H2O)/6.0221367E+11
  ZCTOTA(:,JP_AER_H2O,2) = ZSV(:,JP_CH_H2Oj)*ZMI(JP_AER_H2O)/6.0221367E+11
!
! primary organic carbon
  ZCTOTA(:,JP_AER_OC,1) = ZSV(:,JP_CH_OCi)*ZMI(JP_AER_OC)/6.0221367E+11
  ZCTOTA(:,JP_AER_OC,2) = ZSV(:,JP_CH_OCj)*ZMI(JP_AER_OC)/6.0221367E+11
!
! primary black carbon
  ZCTOTA(:,JP_AER_BC,1) = ZSV(:,JP_CH_BCi)*ZMI(JP_AER_BC)/6.0221367E+11
  ZCTOTA(:,JP_AER_BC,2) = ZSV(:,JP_CH_BCj)*ZMI(JP_AER_BC)/6.0221367E+11
!
IF (NSOA .EQ. 10) THEN
  ZCTOTA(:,JP_AER_SOA1,1) = ZSV(:,JP_CH_SOA1i)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA1,2) = ZSV(:,JP_CH_SOA1j)*ZMI(JP_AER_SOA1)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA2,1) = ZSV(:,JP_CH_SOA2i)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA2,2) = ZSV(:,JP_CH_SOA2j)*ZMI(JP_AER_SOA2)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA3,1) = ZSV(:,JP_CH_SOA3i)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA3,2) = ZSV(:,JP_CH_SOA3j)*ZMI(JP_AER_SOA3)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA4,1) = ZSV(:,JP_CH_SOA4i)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA4,2) = ZSV(:,JP_CH_SOA4j)*ZMI(JP_AER_SOA4)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA5,1) = ZSV(:,JP_CH_SOA5i)*ZMI(JP_AER_SOA5)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA5,2) = ZSV(:,JP_CH_SOA5j)*ZMI(JP_AER_SOA5)/6.0221367E+11

  ZCTOTA(:,JP_AER_SOA6,1) = ZSV(:,JP_CH_SOA6i)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA6,2) = ZSV(:,JP_CH_SOA6j)*ZMI(JP_AER_SOA6)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA7,1) = ZSV(:,JP_CH_SOA7i)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA7,2) = ZSV(:,JP_CH_SOA7j)*ZMI(JP_AER_SOA7)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA8,1) = ZSV(:,JP_CH_SOA8i)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA8,2) = ZSV(:,JP_CH_SOA8j)*ZMI(JP_AER_SOA8)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA9,1) = ZSV(:,JP_CH_SOA9i)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA9,2) = ZSV(:,JP_CH_SOA9j)*ZMI(JP_AER_SOA9)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA10,1) = ZSV(:,JP_CH_SOA10i)*ZMI(JP_AER_SOA10)/6.0221367E+11
  ZCTOTA(:,JP_AER_SOA10,2) = ZSV(:,JP_CH_SOA10j)*ZMI(JP_AER_SOA10)/6.0221367E+11
END IF
!
!-------------------------------------------------------------------------------
!
!*       3    calculate moment 3 from total aerosol mass
!
ZM(:,2) = 0.
ZM(:,5) = 0.
DO JJ = 1,NSP+NCARB+NSOA
    ZM(:,2) = ZM(:,2)+ZCTOTA(:,JJ,1)/ZFAC(JJ) !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
    ZM(:,5) = ZM(:,5)+ZCTOTA(:,JJ,2)/ZFAC(JJ) !==>um3_{aer}/m3_{air} (volume ==> 3rd moment)
ENDDO
!
!-------------------------------------------------------------------------------
!
!*       4    set  moment 0 
!
   ZM(:,1)=   MAX(ZSV(:,JP_CH_M0i) * 1E+6, 1.E-10) ! molec_{aer}/m3_{air}
   ZM(:,4)=   MAX(ZSV(:,JP_CH_M0j) * 1E+6, 1.E-10) ! molec_{aer}/m3_{air}
!
!-------------------------------------------------------------------------------
!
!*       5    set moment 6  ==> um6_{aer}/m3_{air}
!
IF (LVARSIGI) THEN ! set M6 variable standard deviation
  ZM(:,3) = MAX(ZSV(:,JP_CH_M6i), 1E-80)

  ZSIGMA(:)=ZM(:,2)**2/(ZM(:,1)*ZM(:,3))
  ZSIGMA(:)=MIN(1-1E-10,ZSIGMA(:))
  ZSIGMA(:)=MAX(1E-10,ZSIGMA(:))
  ZSIGMA(:)= LOG(ZSIGMA(:))
  ZSIGMA(:)= EXP(1./3.*SQRT(-ZSIGMA(:)))
  WHERE (ZSIGMA(:) > XSIGIMAX)
   ZSIGMA(:) =  XSIGIMAX
  END WHERE
  WHERE (ZSIGMA(:) < XSIGIMIN)
   ZSIGMA(:) =  XSIGIMIN
  END WHERE
  ZM(:,3) = ZM(:,1) &
          * ( (ZM(:,2)/ZM(:,1))**(1./3.)  &
          * exp(-(3./2.)*log(ZSIGMA(:))**2))**6 &
          * exp(18.*log(ZSIGMA(:))**2)

 IF(PRESENT(PSIG1D)) PSIG1D(:,1) = ZSIGMA(:)

ELSE ! fixed standard deviation
 ZM(:,3) = ZM(:,1) &
          * ( (ZM(:,2)/ZM(:,1))**(1./3.)  &
          * exp(-(3./2.)*log(XINISIGI)**2))**6 &
          * exp(18.*log(XINISIGI)**2)

 IF(PRESENT(PSIG1D)) PSIG1D(:,1) = XINISIGI
END IF

IF (LVARSIGJ) THEN ! set M6 variable standard deviation
  ZM(:,6) = MAX(ZSV(:,JP_CH_M6j), 1E-80)

  ZSIGMA(:)=ZM(:,5)**2/(ZM(:,4)*ZM(:,6))
  ZSIGMA(:)=MIN(1-1E-10,ZSIGMA(:))
  ZSIGMA(:)=MAX(1E-10,ZSIGMA(:))
  ZSIGMA(:)= LOG(ZSIGMA(:))
  ZSIGMA(:)= EXP(1./3.*SQRT(-ZSIGMA(:)))
  WHERE (ZSIGMA(:) > XSIGJMAX)
    ZSIGMA(:) =  XSIGJMAX
  END WHERE
  WHERE (ZSIGMA(:) < XSIGJMIN)
    ZSIGMA(:) =  XSIGJMIN
  END WHERE

  ZM(:,6) = ZM(:,4) &
          * ( (ZM(:,5)/ZM(:,4))**(1./3.)  &
          * exp(-(3./2.)*log(ZSIGMA(:))**2))**6 &
          * exp(18.*log(ZSIGMA(:))**2)

  IF(PRESENT(PSIG1D)) PSIG1D(:,2) = ZSIGMA(:)

ELSE ! fixed standard deviation
  ZM(:,6) = ZM(:,4) &
          * ( (ZM(:,5)/ZM(:,4))**(1./3.)  &
          * exp(-(3./2.)*log(XINISIGJ)**2))**6 &
          * exp(18.*log(XINISIGJ)**2)

  IF(PRESENT(PSIG1D)) PSIG1D(:,2) = XINISIGJ
END IF


!-------------------------------------------------------------------------------
!
!*       6    calculate modal parameters from moments
!
DO JN=1,JPMODE
  IF(PRESENT(PN1D)) PN1D(:,JN) = ZM(:,NM0(JN))

  IF(PRESENT(PRG1D)) PRG1D(:,JN)=(ZM(:,NM3(JN))**4. &
              / (ZM(:,NM6(JN))*ZM(:,NM0(JN))**3.))**(1./6.)

ENDDO
!
IF(PRESENT(PCTOTA)) PCTOTA(:,:,:) = ZCTOTA(:,:,:)
IF(PRESENT(PM1D)) PM1D(:,:) = ZM(:,:)
!
IF (LHOOK) CALL DR_HOOK('MODE_AERO_PSD:PPP2AERO1D',1,ZHOOK_HANDLE)
END SUBROUTINE PPP2AERO1D
!
!
END MODULE MODE_AERO_PSD
