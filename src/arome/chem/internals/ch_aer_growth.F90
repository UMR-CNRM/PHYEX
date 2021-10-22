!     ######spl
     SUBROUTINE CH_AER_GROWTH(PM,PSIG0, PRG0, PDMCOND,PDENAIR,&
                                PGASMW,PPGAS,PTGAS,PRH, POM,&
                                PSO4RAT, PDT)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ##############################################
!!
!!   PURPOSE
!!   -------
!!
!!   This routine computes the rate of change due to condensation
!!   and homogene nucleation
!!
!!*************************************************************
!!
!!   Sans test pour savoir si toute la vapeur est utilisee
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
!!    Tulet P.  ajout nucleation Kulmala, 1998
!!
!*************************************************************
! Entry variables:
!
! PM(JPIN)       -Array of moments
! ZT          -Present time in the scheme
!
!*************************************************************
! Exit variables:
!
! ZCOEFM(JPIN)  -Array of moment variation due to condensation
!            and homogeneous nucleation
!
!*************************************************************
! Variables used during the condensation calculation
!
! ZALPHA    - accomodation coefficient
! ZCBAR     - kinetic velocity of vapor molecules (m/s)
! ZDV       - vapor diffusivity (m2/s)
! ZPSIT     - size-independant component of the growth law
! ZSATUR    - saturation ratio of condensed species
!*************************************************************
! Variables used during nucleation calculation
!
! ZCCRIT     -Critical concentration for production of new 
!            particles (kg/m3)
! ZC0        -Initial monomer concentration (kg/m3)
! ZG0        -Critical cluster number
! ZP         -Rate of gas phase production of sulfuric acid
!            concentration C (kg/m3)
! ZSURTEN    -Surface tension (N/m)
! ZTHETA     -Dimensionless surface energy
! ZW0        -Dimensionless energy barrier to nucleation
! ************************************************************
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
USE MODD_CH_AEROSOL
USE MODI_CH_AER_NUCL
!!
IMPLICIT NONE
!
!
!*       0.1   Declarations of arguments
!
REAL, DIMENSION(:,:), INTENT(INOUT) :: PM
REAL, DIMENSION(:,:), INTENT(INOUT) :: PDMCOND, POM
REAL, DIMENSION(:),   INTENT(IN)    :: PDENAIR,PPGAS,PTGAS
REAL, DIMENSION(:),   INTENT(INOUT) :: PRH, PSO4RAT
REAL, DIMENSION(:,:), INTENT(INOUT) :: PSIG0, PRG0
REAL, INTENT(INOUT) :: PGASMW
REAL, INTENT(IN)    :: PDT
!
!*       0.2   Declarations of local variables
!
INTEGER :: JI,JJ

REAL, DIMENSION(SIZE(PM,1),JPMODE,2) :: ZRIK
REAL, DIMENSION(SIZE(PM,1)) :: ZRG,ZLN2S

REAL, DIMENSION(SIZE(PM,1)) :: ZRIKNC,ZRIKFM
REAL, DIMENSION(SIZE(PM,1)) :: ZSIGGAS,ZSIGAIR
REAL, DIMENSION(SIZE(PM,1)) :: ZSIG
REAL, DIMENSION(SIZE(PM,1)) :: ZCBAR
REAL, DIMENSION(SIZE(PM,1),(JPMODE)*3) :: ZMOM
REAL, DIMENSION(SIZE(PM,1)) :: ZCCRIT
REAL, DIMENSION(SIZE(PM,1)) :: ZDTD,ZTINF,ZCSO4SS
REAL, DIMENSION(SIZE(PM,1)) :: ZDMDT,ZDNDT,ZDM3DT,ZDM6DT
REAL, DIMENSION(SIZE(PM,1)) :: ZAL, ZJA, ZSULF

REAL :: ZDV,ZALPHA

REAL :: ZMSO4

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_GROWTH',0,ZHOOK_HANDLE)
ZALPHA=0.05
ZMSO4 = 98.
PDMCOND(:,:)=0.d0
!
!-------------------------------------------------------------------------------
!
! Pour l'instant seul H2SO4 peut nucleer, d'une part
! de facon homogene (creation de nouvelles particules)
! d'autre part de facon heterogene (sur les particules
! d'aerosol deja existante)

!*******************************************************
! Compute the binary diffusivity of the gaseous species
!*******************************************************
ZSIGAIR(:)=(PGASMW/1000.*3./(6.023e23*4.*XPI*PDENAIR(:)))**(1./3.)

ZSIGGAS(:)=(ZMSO4/1000.*3./(6.023e23*4.*XPI*XRHOI(2)))**(1./3.)
ZSIG(:)=(ZSIGGAS(:)+ZSIGAIR(:))/2.
ZCBAR(:)=SQRT(8.*PTGAS(:)*8.31441/(XPI*ZMSO4*1.e-3))
ZDV=0.08e-4

!*************************
! Compute the Omega terms
!*************************

DO JI=1,JPMODE
     
  ZRG(:)=PRG0(:,JI)*1.e-6
  ZLN2S(:)=PSIG0(:,JI)**2

  DO JJ=1,6

    ZMOM(:,JJ)=PM(:,NM0(JI))*ZRG(:)**JJ*exp(real(JJ)**2*ZLN2S(:)/2.)

  ENDDO
  

  ZRIKFM(:)=XPI*ZALPHA*ZCBAR(:)/8.*ZMOM(:,2)
  ZRIKNC(:)=XPI*ZDV/2.*ZMOM(:,1)

  ZRIK(:,JI,1)=ZRIKFM(:)*(ZRIKNC(:)/(ZRIKFM(:)+ZRIKNC(:)))
     
  ZRIKFM(:)=XPI*ZALPHA*ZCBAR(:)/8.*ZMOM(:,5)
  ZRIKNC(:)=XPI*ZDV/2.*ZMOM(:,4)

  ZRIK(:,JI,2)=ZRIKFM(:)*(ZRIKNC(:)/(ZRIKFM(:)+ZRIKNC(:)))

ENDDO
POM(:,1)=(ZRIK(:,1,1)/(ZRIK(:,1,1)+ZRIK(:,2,1)))
POM(:,2)=(ZRIK(:,2,1)/(ZRIK(:,1,1)+ZRIK(:,2,1)))
!
IF (CNUCLEATION=='KERMINEN') THEN
!******************************************************************
! Debut de la partie nucleation homogene en utilisant l'approche de
! Kerminen et Wexler (1994)
!******************************************************************

  ZCCRIT(:)=0.16*exp(0.1*PTGAS(:)-3.5*PRH(:)-27.7)
!KS: suppress nucleation
!  ZCCRIT(:)=1E20
  
!  ZTINF, the time constant for particles to condense onto
!  existing particles is given by Tinf=1/(dM3i/dt+dm3j/dt)
!  where M3i and M3j are in third moment par CM3

   ZDTD(:)=8.*(ZRIK(:,1,1)+ZRIK(:,2,1))
   ZTINF(:)=1./ZDTD(:)

   ZCSO4SS(:)=PSO4RAT(:)*(ZMSO4/6.0221367E+11)*ZTINF(:)

  DO JI = 1,SIZE(PM,1)

    IF (ZCSO4SS(JI) <= ZCCRIT(JI))  THEN  !No nucleation
      ZDNDT(JI)=0.
      ZDMDT(JI)=0.
      ZDM6DT(JI)=0.
    ELSE                        !Calculate nucleation

!  Nucleation of particles from excess mass concentration of sulfuric acid
!  above critical mass concentration of sulfuric acid
!  and condensation of the remaining mass
 
     ZDMDT(JI)=ZDTD(JI)*(ZCSO4SS(JI)-ZCCRIT(JI))
! Les nouvelles particules fraichement crees sont inclues dans le mode
! d'aitken avec les parametres d'initialisation au niveau de la distribution
     ZDNDT(JI)=ZDMDT(JI)*1.e-18/XFAC(JP_AER_SO4)/((0.0025e-6)**3*exp(9./2.*log(1.5)**2))
     ZDM6DT(JI)=ZDNDT(JI)*(0.0025)**6*exp(18.*log(1.5)**2)
!   write(*,*) 'Nucleation: ','DNDT= ',ZDNDT,' DM6DT= ',ZDM6DT

    ENDIF

  ENDDO
ELSE
     ZDNDT(:)=0.
     ZDMDT(:)=0.
     ZDM6DT(:)=0.
ENDIF
!
IF (CNUCLEATION=='KULMALA') THEN
! compute nucleation rate
!
  ZSULF(:) = PSO4RAT(:) * PDT
!
  CALL CH_AER_NUCL(PRH,PTGAS,ZSULF,ZJA,ZAL,SIZE(PSO4RAT,1))
!
! new mass in molec.cm-3.s-1
  ZDMDT(:)= ZAL(:)*ZJA(:) 
! convert into microgram.m-3.s-1
  ZDMDT(:)= ZDMDT(:) * ZMSO4/6.0221367E+11
!
! Les nouvelles particules fraichement crees sont inclues dans le mode
! d'aitken avec les parametres d'initialisation au niveau de la distribution
!
  ZDNDT(:)  = ZDMDT(:)/(XFAC(JP_AER_SO4)*(PRG0(:,1)**3)*EXP(4.5 * PSIG0(:,1)**2))
  ZDM6DT(:) = ZDNDT(:)*(PRG0(:,1)**6*EXP(18.*PSIG0(:,1)**2))

ELSE
  ZDNDT(:)=0.
  ZDMDT(:)=0.
  ZDM6DT(:)=0.
ENDIF
!
! condensation des sulfates
ZDM3DT(:)=PSO4RAT(:)*(ZMSO4/6.0221367E+11)/XFAC(JP_AER_SO4)*1.e-18
!
! Enlever la quantite de 3e moment deja consommee pour la nucleation homogene
ZDM3DT(:)=ZDM3DT(:)-ZDMDT(:)/XFAC(JP_AER_SO4)*1.e-18
!
!
PDMCOND(:,1)=ZDNDT(:)
PDMCOND(:,2)=ZDMDT(:)/XFAC(JP_AER_SO4)
PDMCOND(:,3)=ZDM6DT(:)
!
DO JI=1,JPMODE
  PDMCOND(:,NM3(JI))=PDMCOND(:,NM3(JI))+ZDM3DT(:)*POM(:,JI)*1.e18
  PDMCOND(:,NM6(JI))=PDMCOND(:,NM6(JI))+ZDM3DT(:)*POM(:,JI)*ZRIK(:,JI,2)/ZRIK(:,JI,1)*1.e36
    
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CH_AER_GROWTH',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_GROWTH
