!     ######spl
      SUBROUTINE AERMR_NC ( PHYEX, PTSTEP, PRHODREF, PEXNREF,PPABST, &
                            PCLDFR, PTHT, PRVT, PRCT, KNAERO, &
                            PAER_MR, PSSAT, PAER_NC, PCCN_NC, PIFN_NC)

      USE YOMHOOK,  ONLY : LHOOK, DR_HOOK, JPHOOK
!     ##########################################################################
!
!!**** * -  Calculate the number of condensation nuclei
!!
!!    PURPOSE
!!    -------
!!      Calculate the number of nuclei from the mixing ratio
!!
!!**  METHOD
!!    ------
!!    Considering a log normal distribution for the particles of every specie
!!    the mean cubic radius is used for the caltulation.
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!      Module MODD_CST
!!          XPI
!!          XP00               ! Reference pressure
!!          XRD,XRV            ! Gaz  constant for dry air, vapor
!!          XMD,XMV            ! Molecular weight for dry air, vapor
!!          XCPD               ! Cpd (dry air)
!!          XCL                ! Cl (liquid)
!!          XCI                ! Ci (solid)
!!          XTT                ! Triple point temperature
!!          XLVTT              ! Vaporization heat constant
!!          XALPW,XBETAW,XGAMW ! Constants for saturation vapor pressure
!!                               function over liquid water
!!          XALPI,XBETAI,XGAMI ! Constants for saturation vapor pressure
!!                               function over solid ice
!!    AUTHOR
!!    ------
!!    D. Martin-Perez
!!
!!    MODIFICATIONS
!!    -------------
!!    Original    03.2018
!!    02.2021: Inclusion of IFN
!!    28.06.2022: Modification of SSat to include effect of coarse Sea Salt
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAMETERS
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_AEROSOL_PROP
USE MODD_NRT_AEROSOLS, ONLY: XSSMINLO,XSSFACSS,XIFNMINSIZE,XCLDROPMIN
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PHYEX_t),            INTENT(IN) :: PHYEX
REAL,                     INTENT(IN)    :: PTSTEP  ! Time step
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR  ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! Cloud water m.r. at t
!
INTEGER,                  INTENT(IN)    :: KNAERO   ! Number of aerosol species
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PAER_MR  ! Aerosol mixing ratio (kg/kg)
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSSAT    ! Supersaturation
REAL, DIMENSION(:,:,:,:), INTENT(  OUT) :: PAER_NC  ! Aerosol number concentration (m-3)
REAL, DIMENSION(:,:,:,:), INTENT(  OUT) :: PCCN_NC  ! Cloud condensation nuclei number concentration (m-3) per species
REAL, DIMENSION(:,:,:),   INTENT(  OUT) :: PIFN_NC  ! Ice freezing nuclei number concentration (m-3)

! ****** Used for aerosol . ****
! Weighted concentration
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3),KNAERO) :: &
      ZCONC_AER                                                              ! Weighted concentration
REAL, DIMENSION(SIZE(PEXNREF,1),SIZE(PEXNREF,2),SIZE(PEXNREF,3))        :: &
      ZT,                                                                  & ! Temperature
      ZAER_NCCOARSS, ZAER_NCSS, ZAER_NCDU
!
! Working variables
!
REAL :: ZZT        ! Temperature
REAL :: ZZZSSW     ! Supersaturation
REAL :: ZACTAE     ! Activated nuclei
REAL :: ZPPABST    ! pressure
!
REAL :: ZDKA, ZDKB ! Kohler parameters (a,b)
REAL :: ZDKR       ! Critical radius
REAL :: ZDKRMIN    ! Minimum aerosol radius for saturation
REAL :: ZRTINCLOUD

! ****** end aerosol *****
!
!*       0.2   Declarations of local variables :
!
REAL :: ZCCNAE   ! Condensation nuclei
REAL :: ZAERPR, ZR3AE
REAL :: ZERFRAE, ZERFIFN
REAL :: ZWVDIF,ZRA,ZRD,ZW1,ZSSC,ZDELTASS

INTEGER :: IC,JC
INTEGER :: JI,JJ,JK
INTEGER :: IIB,IIE,IJE,IJB
INTEGER :: IKT,IKTB,IKTE
INTEGER :: JCCN

!
!-------------------------------------------------------------------------------
!
!*       1.1     COMPUTE THE LOOP BOUNDS
!               -----------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AERMR_NC',0,ZHOOK_HANDLE)

IIB=1+JPHEXT
IIE=SIZE(PRHODREF,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PRHODREF,2) - JPHEXT
IKT=SIZE(PRHODREF,3)
IKTB=1+JPVEXT
IKTE=IKT-JPVEXT

!
!*       1.2     COMPUTE SOME CONSTANT PARAMETERS
!

! Initialize Variables 

ZDKRMIN=2.0E-8  ! Minimum aerosol radius to be activated (m)
ZERFIFN=0.0

PAER_NC(:,:,:,:)=0.0
PCCN_NC(:,:,:,:)=0.0
PIFN_NC(:,:,:)=0.0

ZCONC_AER(:,:,:,:)=0.0

ZRTINCLOUD=PHYEX%RAIN_ICE_DESCRN%XRTMIN(2)

! Temperature
ZT(:,:,:) = PTHT(:,:,:) * ( PPABST(:,:,:) / XP00 ) ** (XRD/XCPD)

! Compute the number of aerosol nuclei from the mixing ratio.
DO JC=1,KNAERO
   !     rm**3 exp(4.5 (ln s)**2 )
   ZR3AE=XR3CORR(JC)*1.0E-18
   ZAERPR=XRHOAE(JC)*ZR3AE
   PAER_NC(:,:,:,JC)=(0.75*PAER_MR(:,:,:,JC)*PRHODREF(:,:,:))/(XPI*ZAERPR)
ENDDO

! Total number of aerosol nuclei per specie
ZCONC_AER(:,:,:,:)=MAX(PAER_NC(:,:,:,:),0.0)

! Estimation of the reduction of the supersaturation due to the presence coarse sea salt
ZWVDIF=0.211E-4  ! Water vapor diffusivity (m**2/s)
ZRA=5.0E-6       ! reference aerosol radius (m)
ZAER_NCCOARSS(:,:,:)=0.0
ZDELTASS=0.03E-2
IF ( NCOARSEASALT > 0 ) THEN
   DO JK = IKTB, IKTE
      DO JI = IIB,IIE
         DO JJ = IJE,IJB
            IF (PSSAT(JI,JJ,JK)>0.0) THEN
               ZZT = ZT(JI,JJ,JK)
               ZDKA=0.32538E-6/ZZT
               DO IC=1,NCOARSEASALT
                  JC=XCOARSEASALT(IC)
                  ZAER_NCCOARSS(JI,JJ,JK)=ZAER_NCCOARSS(JI,JJ,JK)+ &
                          & ZCONC_AER(JI,JJ,JK,JC)
                  ZDKB=XKHYGROS(JC)*(ZRA)**3.0
               ENDDO
               ! Critical radius and supersaturation
               ZRD=(3.0*ZDKB/ZDKA)**0.5
               ZSSC=(4.0*ZDKA**3.0/(27.0*ZDKB))**0.5
               ZW1=MAX(0.0,4.0*XPI*ZWVDIF*(XSSFACSS*PTSTEP)*ZRD*((XSSMINLO-ZDELTASS)-ZSSC))
               PSSAT(JI,JJ,JK)=PSSAT(JI,JJ,JK)-MIN(ZDELTASS, &
                        & ZW1*(MIN(XCLDROPMIN,ZAER_NCCOARSS(JI,JJ,JK))))
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF

! *****************************************************
! Compute total number of IFN

DO IC=1,NIFN
   JC=XIFNI(IC)
   IF ( JC > KNAERO ) CYCLE
   ! Number of c.n. bigger than XIFNMINSIZE micrometers
   ZERFIFN=ERFPDF(XIFNMINSIZE,XMRAE(JC),XSDAE(JC))
   PIFN_NC(:,:,:)=PIFN_NC(:,:,:)+PAER_NC(:,:,:,JC)*(XERFUP(JC)-ZERFIFN)/(XERFUP(JC)-XERFDOWN(JC))
ENDDO


! *****************************************************
! Compute Active CCN
!

! Number of aerosol species to become CCN


!  1. Consider grid points with Cloud water content greater than ZRTINCLOUD

DO JK = IKTB,IKTE
   DO JI = IIB,IIE
      DO JJ = IJE,IJB
         IF ( PRCT(JI,JJ,JK)>ZRTINCLOUD ) THEN
            ZZT = ZT(JI,JJ,JK)
            ZZZSSW=PSSAT(JI,JJ,JK)
            ZPPABST=PPABST(JI,JJ,JK)

            ZDKA=0.32538/ZZT

!  2. Calculate activated cloud condensation nuclei
            DO JC=1,NCCN
               ZACTAE=0.0
               JCCN=XACCNI(JC)
               IF ( JCCN > KNAERO ) CYCLE
               ZCCNAE=0.0
               ZERFRAE=0.0
               ZDKR=0.0
               ZDKB=XKHYGROS(JCCN)
               ZDKR=MAX(ZDKA/(6.75*ZDKB)**(0.33)/&
                       & (ZZZSSW)**(0.66),ZDKRMIN)
               ZCCNAE=ZCONC_AER(JI,JJ,JK,JCCN)/(XERFUP(JCCN)-XERFDOWN(JCCN))
               ZERFRAE=MAX(XERFDOWN(JCCN),0.5*&
                       & erf(0.7071*LOG(ZDKR/XMRAE(JCCN))/LOG(XSDAE(JCCN))))
               ZACTAE=MAX(ZCCNAE*(XERFUP(JCCN)-ZERFRAE),0.0)
! Total cloud condensation nuclei
               PCCN_NC(JI,JJ,JK,JC)=MAX(ZACTAE,0.0)
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Total aerosol nuclei
PAER_NC(:,:,:,:)=MAX(PAER_NC(:,:,:,:),0.0)

! End of calculation of Active CCN
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('AERMR_NC',1,ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
!
CONTAINS

!
!

FUNCTION ERFPDF(PRAD,PMODRAD,PSIGMA) RESULT (PERFPDF)
!
! 0.5*erf(log(PRAD/Rg)/(sqrt(2.0))*log(sigma)
!
!
  IMPLICIT NONE
!
  REAL   :: PRAD
  REAL   :: PMODRAD
  REAL   :: PSIGMA
  REAL   :: PERFPDF
!
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('AERMR_NC:ERFPDF',0,ZHOOK_HANDLE)
!
  PERFPDF = 0.5*erf(0.7071*LOG(PRAD/PMODRAD)/LOG(PSIGMA))
!
  IF (LHOOK) CALL DR_HOOK('AERMR_NC:ERFPDF',1,ZHOOK_HANDLE)
!
END FUNCTION ERFPDF
!
!
END SUBROUTINE AERMR_NC
