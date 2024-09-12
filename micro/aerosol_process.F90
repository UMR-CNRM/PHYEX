!     ######spl
  SUBROUTINE AEROSOL_PROCESS(PHYEX,KKA,KKU,KKL,KNAERO,PTSTEP,&
        &  PPABST,PTHT,PDZZ,PRHODREF,PRCT,PRIT,&
        &  PCLDFR,PFPR, &
        &  PLSM,PUGST,PVGST, &
        &  PAERMR,PAERTEND)

      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ##########################################################################
!
!!    Aerosol removal processes
!!
!!    PURPOSE
!!    -------
!!     Calculate the scavenging and deposition of aerosols.
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!     Module MODD_CST
!!          XPI
!!          XG           : Gravity constant
!!          XP00         : Reference pressure
!!          XRD          : Gaz  constant for dry air
!!          XCPD         : Cpd (dry air)
!!          XRHOWL       : density of liquid water
!!          XRHOLI       : density of ice water
!!     Module MODD_RAIN_ICE_DESCR
!!          XRTMIN       : Min values allowed for mixing ratios
!!     Module MODD_AEROSOL_PROP
!!          XVDRY        : Dry Deposition Velocity
!!          XFRICS       : Fracion of aerosols in aquaeous phase
!!                       : for in-cloud scavenging
!!
!!    REFERENCE
!!    ---------
!!       Mocrette et al. 2009 and Redy, 2019.
!!
!!    AUTHOR
!!    ------
!!       Daniel Martin Perez (AEMET)
!!
!!    MODIFICATIONS
!!    -------------
!!       Original:    2018
!!       D. Martin :  27.03.2023 : Adition of sea salt emission.
!!       D. Martin :  28.06.2022 : Modification of scavenging, addition of 
!!                                 snow fluxes.
!!
!-------------------------------------------------------------------------
!!
!*      0. DECLARATIONS
!          ------------

USE MODD_PARAMETERS,     ONLY : JPHEXT, JPVEXT
USE MODD_CST,            ONLY : XG, XPI, XRHOLW, XRHOLI, XP00, XRD, XCPD
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_AEROSOL_PROP,   ONLY : XVDRY, XFRICS, NCOAR, XCOARSE, XSDAE, XRHOAE, XMRAE, &
         & XERFUP, XERFDOWN, XKHYGROS, NDRYDEP, XDRYDEP, NSEASALT, XSEASALT, &
         & NDUST, XDUST, XINTSS
USE MODD_NRT_AEROSOLS,   ONLY : LAERDRDEP, LAERSSEM
!
IMPLICIT NONE
!
TYPE(PHYEX_t),            INTENT(IN) :: PHYEX
INTEGER,                  INTENT(IN)    :: KKA      ! near ground array index
INTEGER,                  INTENT(IN)    :: KKU      ! uppest atmosphere array index
INTEGER,                  INTENT(IN)    :: KKL      ! vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)    :: KNAERO   ! number of aerosol species
REAL,                     INTENT(IN)    :: PTSTEP   ! Double Time step
                                                    ! (single if cold start)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST   ! Pressure
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PDZZ     ! Layer thickness (m)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCLDFR   ! Convective Mass Flux Cloud fraction
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PFPR     ! upper-air precipitation fluxes
REAL, DIMENSION(:),       INTENT(IN)    :: PLSM     ! Land/Sea Mask
REAL, DIMENSION(:),       INTENT(IN)    :: PUGST    ! 10m u-wind gust
REAL, DIMENSION(:),       INTENT(IN)    :: PVGST    ! 10m v-wind gust
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PAERMR
REAL, DIMENSION(:,:,:,:), INTENT(  OUT) :: PAERTEND ! Aerosol tendencies: emission +  scavenging coefficient rate
!
!*       0.2  declaration of local variables
!
INTEGER :: IIB           !  Define the domain where is
INTEGER :: IIE           !  the microphysical sources have to be computed
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB           !
INTEGER :: IKE           !
INTEGER :: JC, IC
INTEGER :: JK, JI, JJ

REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3),KNAERO) :: &
      ZW1, ZW2, ZW3, ZW4,                                                     & ! work arrays
      ZAER_MR,                                                                &
      ZAER_DEP                                                                  ! Aerosol deposition
REAL, DIMENSION(SIZE(PRHODREF,1),SIZE(PRHODREF,2),SIZE(PRHODREF,3))        :: &
      ZT,                                                                     & ! Temperature
      ZVELSED,                                                                & ! Velocity sedimentation
      ZMFREEPATH,                                                             & ! mean free path
      ZCUNN,                                                                  & ! Cunningham factor for sedimentation
      ZPFPRR, ZPFPRS

REAL, DIMENSION(SIZE(PRHODREF,1),3)                                        :: &
      ZWS_SS                                                                    ! work array

! Sink and source at surface of aerosol in kg/m2
REAL, DIMENSION(SIZE(PRHODREF,1))                                          :: &
      ZAERSNKDU,                                                              & ! Desert dust sink
      ZAERSNKSS,                                                              & ! Sea salt sink
      ZAERSRCSS                                                                 ! Sea salt source

! XRTMIN = Minimum value for the mixing ratio
REAL :: ZINVTSTEP
REAL :: ZPFMIN     ! Minimum Flux precipitation
REAL :: ZEF_RAIN, ZEF_SNOW
REAL :: ZRAY_RAIN, ZRAY_SNOW
REAL :: ZCLFK
REAL :: ZVDRY
REAL :: ZVISC
REAL :: ZAIRMAS
REAL :: ZBOLTZMANN
REAL :: ZAER_R
REAL :: Z10WIND

REAL :: ZDKA, ZDKB, ZDKR, ZDKRMIN, ZFRICS, ZERFRAE

!-------------------------------------------------------------------------------
!
!*       1. Parameters for aerosol wet deposition
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AEROSOL_PROCESS',0,ZHOOK_HANDLE)


IIB=1+JPHEXT
IIE=SIZE(PDZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=SIZE(PDZZ,2) - JPHEXT
IKB=KKA+JPVEXT*KKL
IKE=KKU-JPVEXT*KKL
!
!*       1.2     COMPUTE SOME CONSTANT PARAMETERS
!
ZINVTSTEP=1./PTSTEP

!
!*       2.    compute the fluxes
!

!ZPFMIN=1.0E-12
ZPFMIN=0.0
ZEF_RAIN=0.001
ZEF_SNOW=0.01
ZRAY_RAIN=1.0E-3
ZRAY_SNOW=1.0E-3
ZVISC=1.78E-5      ! Air viscosity (Pa*s)
ZAIRMAS=28.97E-3   ! Air mean molar mass(Kg*mol-1)
ZBOLTZMANN=8.31445 ! (J*K-1*mol-1)
ZDKRMIN=2.0E-8

ZAER_MR(:,:,:,:)=MAX(PAERMR(:,:,:,:),0.0)

ZW1(:,:,:,:) = 0.
ZW2(:,:,:,:) = 0.
ZW3(:,:,:,:) = 0.
ZW4(:,:,:,:) = 0.

ZWS_SS(:,:)=0.0

ZT(:,:,:) = PTHT(:,:,:) * ( PPABST(:,:,:) / XP00 ) ** (XRD/XCPD)

ZAERSNKDU(:)=0.0
ZAERSNKSS(:)=0.0
ZAERSRCSS(:)=0.0

! ************************************************
! AEROSOL REMOVAL PROCESSES
! ************************************************

! Rain and snow fluxes initialized with sedimentation fluxes

DO JI = IIB,IIE
   DO JJ = IJE,IJB
      DO JK = IKE+1 , IKB, -1*KKL
         ZPFPRR(JI,JJ,JK)=MAX(PFPR(JI,JJ,JK,3),0.0)
         ZPFPRS(JI,JJ,JK)=MAX(PFPR(JI,JJ,JK,5),0.0)
      ENDDO
   ENDDO
ENDDO


! Aerosol wet deposition
! (Loop from top to bottom KKL=-1 ARO)
DO JK = IKE+1 , IKB, -1*KKL
   DO JI = IIB,IIE
      DO JJ = IJE,IJB
         ! inside cloud
         IF  (ZPFPRR(JI,JJ,JK) > ZPFMIN .OR. &
            & ZPFPRS(JI,JJ,JK) > ZPFMIN) THEN
            ! inside cloud
            ZCLFK=PCLDFR(JI,JJ,JK)
            ZDKA=0.32538/ZT(JI,JJ,JK)
            IF ((PRCT(JI,JJ,JK) > PHYEX%RAIN_ICE_DESCRN%XRTMIN(2) .OR. &
               & PRIT(JI,JJ,JK) > PHYEX%RAIN_ICE_DESCRN%XRTMIN(4)) .AND. ZCLFK > 0.1) THEN
               DO JC=1,KNAERO
                  ZW1(JI,JJ,JK,JC)=((ZPFPRR(JI,JJ,JK+KKL)-ZPFPRR(JI,JJ,JK))+&
                  & (ZPFPRS(JI,JJ,JK+KKL)-ZPFPRS(JI,JJ,JK)))/ &
                  & (PRCT(JI,JJ,JK)+PRIT(JI,JJ,JK))/ &
                  & (PRHODREF(JI,JJ,JK)*PDZZ(JI,JJ,JK)*ZCLFK)
                  ZDKB=XKHYGROS(JC)
                  IF ( ZDKB == 0.0 ) THEN
                     ZFRICS=XFRICS(JC)
                  ELSE
                     ZDKR=MAX(ZDKA/(6.75*ZDKB)**(0.33)/(0.05E-2)**(0.66),ZDKRMIN)
                     ZERFRAE=MAX(XERFDOWN(JC),0.5*erf(0.7071*LOG(ZDKR/XMRAE(JC))/LOG(XSDAE(JC))))
                     ZFRICS=MIN((XERFUP(JC)-ZERFRAE)/(XERFUP(JC)-XERFDOWN(JC)),XFRICS(JC))
                  ENDIF
                  ZW1(JI,JJ,JK,JC)=ZCLFK*ZFRICS*ZAER_MR(JI,JJ,JK,JC)*MIN(MAX(-1.0,ZW1(JI,JJ,JK,JC)),0.0)
               ENDDO
            ELSE
               ! below cloud (rain)
               DO JC=1,KNAERO
                  ZW2(JI,JJ,JK,JC)=MIN(-0.75*ZPFPRR(JI,JJ,JK)*ZEF_RAIN/&
                  & (ZRAY_RAIN*XRHOLW)- &
                  & 0.75*ZPFPRS(JI,JJ,JK)*ZEF_SNOW/(ZRAY_SNOW*XRHOLI),0.0)* &
                  & ZAER_MR(JI,JJ,JK,JC)
               ENDDO
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDDO

! Aerosol dry deposition
! Surface layer
!
IF (LAERDRDEP) THEN
   DO JI = IIB,IIE
      DO JJ = IJE,IJB
         IF (PRCT(JI,JJ,IKB) < PHYEX%RAIN_ICE_DESCRN%XRTMIN(2) .OR. &
            & (ZPFPRR(JI,JJ,IKB) < ZPFMIN .AND. &
            & ZPFPRS(JI,JJ,IKB) < ZPFMIN)) THEN
            DO IC=1,NDRYDEP ! only for those aerosol with parametrized emission
               JC=XDRYDEP(IC)
               IF ( JC > KNAERO ) CYCLE
               ZVDRY=(1.0-PLSM(JI))*XVDRY(JC,1)+PLSM(JI)*XVDRY(JC,2)
               ZW3(JI,JJ,IKB,JC)=-0.01*ZVDRY*ZAER_MR(JI,JJ,IKB,JC)/&
                  & PDZZ(JI,JJ,IKB)
            ENDDO
         ENDIF
      ENDDO
   ENDDO
ENDIF

ZAER_DEP(:,:,:,:)=ZW1(:,:,:,:)+ZW2(:,:,:,:)+ZW3(:,:,:,:)


! Dry sedimentation only for coarse modes (Redy, 2019)
ZMFREEPATH(:,:,:)=ZVISC/PPABST(:,:,:)*SQRT(0.5*XPI*ZBOLTZMANN*ZT(:,:,:)/ZAIRMAS)

DO IC=1,NCOAR
   JC=XCOARSE(IC)
   IF ( JC > KNAERO ) CYCLE
   ZAER_R=XMRAE(JC)*EXP(0.5*(LOG(XSDAE(JC))**2))*1.0E-6
   ZCUNN(:,:,:)=1.0+ZMFREEPATH(:,:,:)/ZAER_R*&
           & (1.257+0.4*EXP(-1.1*ZAER_R/ZMFREEPATH(:,:,:)))
   ZVELSED(:,:,:)=MIN(0.22*(ZAER_R**2)*XRHOAE(JC)*XG*ZCUNN(:,:,:)/ZVISC, &
           & PDZZ(:,:,:)*ZINVTSTEP)

   DO JK = IKE+1 , IKB, -1*KKL
      ZW4(:,:,JK,JC)=(ZAER_MR(:,:,JK+KKL,JC)*ZVELSED(:,:,JK+KKL)-&
              & ZAER_MR(:,:,JK,JC)*ZVELSED(:,:,JK))/PDZZ(:,:,JK)
   ENDDO

ENDDO

! Sum all the aerosol removal contributions

PAERTEND(:,:,:,:)=MAX(ZAER_DEP(:,:,:,:)+ZW4(:,:,:,:),-0.8*ZAER_MR(:,:,:,:)*ZINVTSTEP)

! Aerosol sinks on the surface (kg/m2/s)

DO JI = IIB,IIE
   ZAERSNKDU(JI)=0.0
   DO JJ = IJE,IJB
      DO JK = IKE+1 , IKB, -1*KKL
         DO IC=1,NSEASALT  ! Sea Salt species
            JC=XSEASALT(IC)
            ! diagnostics of all aerosol deposition in kg/m2/s
            ZAERSNKSS(JI)=ZAERSNKSS(JI)+(MIN(ZW1(JI,JJ,JK,JC),0.0)+ &
                    & MIN(ZW2(JI,JJ,JK,JC),0.0))*PRHODREF(JI,JJ,JK)*&
                    & PDZZ(JI,JJ,JK)
         ENDDO
         DO IC=1,NDUST  ! Dust species
            JC=XDUST(IC)
            ! diagnostics of all aerosol deposition in kg/m2/s
            ZAERSNKDU(JI)=ZAERSNKDU(JI)+(MIN(ZW1(JI,JJ,JK,JC),0.0)+ &
                    & MIN(ZW2(JI,JJ,JK,JC),0.0))*PRHODREF(JI,JJ,JK)*&
                    & PDZZ(JI,JJ,JK)
         ENDDO
      ENDDO
      DO IC=1,NSEASALT  ! Sea Salt species
         JC=XSEASALT(IC)
         ZAERSNKSS(JI)=ZAERSNKSS(JI)+(MIN(ZW3(JI,JJ,IKB,JC),0.0)+&
                  & MIN(ZW4(JI,JJ,IKB,JC),0.0))*&
                  & PRHODREF(JI,JJ,IKB)*PDZZ(JI,JJ,IKB)
      ENDDO
      DO IC=1,NDUST  ! Dust species
         JC=XDUST(IC)
         ZAERSNKDU(JI)=ZAERSNKDU(JI)+(MIN(ZW3(JI,JJ,IKB,JC),0.0)+&
                  & MIN(ZW4(JI,JJ,IKB,JC),0.0))*&
                  & PRHODREF(JI,JJ,IKB)*PDZZ(JI,JJ,IKB)
      ENDDO
   ENDDO
   ZAERSNKDU(JI)=MIN(ZAERSNKDU(JI),0.0)
   ZAERSNKSS(JI)=MIN(ZAERSNKSS(JI),0.0)
ENDDO


!*********************************
!  AEROSOL EMISSIONS
!*********************************
!
! Sea salt production (Monahan,1986)
!
IF (LAERSSEM) THEN
   DO JI = IIB,IIE
      ZAERSRCSS(JI)=0.0
      IF ( PLSM(JI)<0.5 ) THEN
         Z10WIND=SQRT(PUGST(JI)**2+PVGST(JI)**2)
         DO IC=1,NSEASALT  ! Sea Salt species
            JC=XSEASALT(IC)
            ZWS_SS(JI,JC)=1.373*Z10WIND**(3.41)*XINTSS(JC)* &
                  & (XRHOAE(JC)/8.0)*(4.0*XPI/3.0)*1.0E-18
            ZAERSRCSS(JI)=ZAERSRCSS(JI)+ZWS_SS(JI,JC)
            DO JJ = IJE,IJB
               PAERTEND(JI,JJ,IKB,JC)=PAERTEND(JI,JJ,IKB,JC)+ZWS_SS(JI,JC)/&
                     & (PRHODREF(JI,JJ,IKB)*PDZZ(JI,JJ,IKB))
            ENDDO
         ENDDO
      ENDIF
   ENDDO
ENDIF


IF (LHOOK) CALL DR_HOOK('AEROSOL_PROCESS',1,ZHOOK_HANDLE)
!
END SUBROUTINE AEROSOL_PROCESS
