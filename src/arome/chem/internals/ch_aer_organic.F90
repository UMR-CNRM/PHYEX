!     ######spl
     SUBROUTINE CH_AER_ORGANIC(PCTOTG, PCTOTA, PRV, PDENAIR, PPRESSURE, PTEMP, PRC, POM,&
                               PCCTOT,PSIG0, PRG0, PDT, PSOLORG)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!###########################################################################################
!!
!!   PURPOSE
!!   -------
!!   solve the organic thermodynamic balance , if use CACM our ReLACS2 chemical
!!   scheme
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    P. Tulet  (GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CH_AEROSOL
USE MODI_CH_AER_MPMPO
USE MODI_CH_AER_PUN
!USE MODI_CH_AER_DIFF

IMPLICIT NONE
!!
!!  Declaration arguments
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PCTOTG, POM
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PCTOTA, PCCTOT
REAL, DIMENSION(:),     INTENT(IN)    :: PRV, PDENAIR, PPRESSURE, PTEMP, PRC
REAL, DIMENSION(:,:),   INTENT(IN)    :: PSIG0, PRG0
REAL,                   INTENT(IN)    :: PDT
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSOLORG


!!
!!  Declaration variables internes
!! 
INTEGER :: II,IJ, JJ
REAL, DIMENSION(SIZE(PCTOTA,1),NSP+NCARB+NSOA) :: ZTOT, ZTOTG, ZTOTNEW, ZTOTGNEW
REAL, DIMENSION(SIZE(PCTOTA,1),NSP+NCARB+NSOA) :: ZDEL
REAL, DIMENSION(SIZE(PCTOTA,1),NSP+NCARB+NSOA) :: ZDIFF
REAL, DIMENSION(SIZE(PCTOTA,1)) :: ZL3DIFF
REAL, DIMENSION(SIZE(PCTOTA,1)) :: ZLWC, ZRH
REAL, DIMENSION(SIZE(PCTOTA,1)) :: ZPROTON, ZNEUTRAL
REAL, DIMENSION(SIZE(PCTOTA,1)) :: ZPSAT, ZPKM, ZPKH2O
INTEGER, SAVE :: ZZZ=0


! INPUTS
! PCTOTG : array of gas phase in microgram / m3
! PCTOTA : array of aerosol phase in microgram / m3 (last index represents
! lognormal modes)
! PCCTOT : percentage of aerosol mass composition 
! For both PCTOTG, PCTOTA and PCCTOT : First dimension is the data of simulation domain
!                                      Second dimension represent the aerosol species
!                                      1 => NSP : mineral aerosol  species
!                                      NSP+1 => NSP + NCARB : primary aerosol species
!                                      NSP + NCARB +1 => NSP + NCARB + NSOA : secondary organic aerosol
!                                      Third dimension (only for PCTOTA and PCCTOT): lognormal mode
! PRESSURE : pressure in hPa
! PDENAIR  : air density
! PTEMP    : air temperature in K
! PRC      : cloud mixing ration in kg/kg (not use now)
! PRV      : vapor mixing ration in kg/kg
! PSIG0    : log(sigma) of the distribution 
! PRG0     : Mean radius of the lognormal distribution in number of particles 
! PDT      : Time step (in s)

!*****************************************************************
!*****************************************************************
! SOLVEUR DE L'EQUILIBRE CHIMIQUE ORGANIC
!*****************************************************************
!*****************************************************************
! Option  flag to turn on type B oa module if == 1
!          (type a is iterated, but no need to rerun type b)
! Compute relative humidity
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_ORGANIC',0,ZHOOK_HANDLE)
ZPKM(:) = 1E-3*PDENAIR(:) * 6.0221367E+23 / 28.9644
ZPKH2O(:) = ZPKM(:)*1.6077*PRV(:)
ZPSAT(:)=0.611*EXP(17.2694*(PTEMP(:)-273.16)/(PTEMP(:)-35.86))
ZPSAT(:)=ZPSAT(:)*1000.
ZRH(:)=(ZPKH2O(:)/(ZPKM(:)*1.6077))*PPRESSURE(:)/&
     &(0.622+(ZPKH2O(:)/(ZPKM(:)*1.6077)))/ZPSAT(:)
ZRH(:) = MIN(0.95, MAX(ZRH(:), .1)) ! until 0.95 thermodynamic code is not valid


! Mass need to be positive
 PCTOTA(:,:,:)= MAX (PCTOTA(:,:,:),1.E-40)
 PCTOTG(:,:)  = MAX (PCTOTG(:,:),1.E-40)


!******************************************************************
! Calcul de la repartition des differentes especes entre les modes
! pour pouvoir conserver celle ci apres l'equilibre chimique
!******************************************************************

ZTOT(:,:)=0.
DO II=1,NSP+NCARB+NSOA
  ZTOTNEW(:,II)=0.
  ZTOT(:,II)=PCTOTA(:,II,1)+PCTOTA(:,II,2)
  ZTOTG(:,II)=PCTOTG(:,II)
ENDDO
  ZTOTNEW(:,:) = ZTOT(:,:)

ZLWC(:) = ZTOT(:,JP_AER_H2O)

ZPROTON(:) = 1E-10

! Neutralisation H2SO4 + NH3 => H2SO4NH3
ZNEUTRAL(:) =  PCTOTA(:,JP_AER_SO4,1) / 98. - PCTOTA(:,JP_AER_NH3,1) / 17. + &
               PCTOTA(:,JP_AER_SO4,2) / 98. - PCTOTA(:,JP_AER_NH3,2) / 17.

WHERE (ZNEUTRAL(:) .GT. 0.)  ! SO4 > NH3 / excess of H2SO4 => SO4-- + 2*H+
   ZNEUTRAL(:) = 2.*ZNEUTRAL(:)
ENDWHERE


WHERE(ZLWC(:) .GT. 1.E-20) 
ZPROTON(:) = ZNEUTRAL(:) +  PCTOTA(:,JP_AER_NO3,1)/63. + PCTOTA(:,JP_AER_NO3,2)/63.
ZPROTON(:) = 1E-10 + ZPROTON(:) / ZLWC(:)  ! proton concentration in mole / g water
ENDWHERE
ZPROTON(:) = MIN(MAX(ZPROTON(:), 1E-17),1E-4)

PSOLORG(:,:) = 0.
IF (TRIM(CORGANIC) == 'PUN') THEN
CALL CH_AER_PUN(ZTOT, ZTOTG, PTEMP, ZRH, ZLWC, ZPROTON, ZTOTNEW, ZTOTGNEW)

ELSE IF (TRIM(CORGANIC) == 'MPMPO') THEN

CALL CH_AER_MPMPO(ZTOT, ZTOTG, PTEMP, ZRH, ZLWC, ZPROTON,ZTOTNEW,ZTOTGNEW,PSOLORG)

ELSE
IF (ZZZ == 0) THEN
PRINT *,' WARNING WARNING WARNING WARNING WARNING WARNING'
PRINT *,' PAS D EQUILIBRE THERMODYNAMIQUE ENTRE LES ORGANIQUES'
PRINT *,' WARNING WARNING WARNING WARNING WARNING WARNING'
ZZZ=1
END IF
END IF


! Calcul de la nouvelle composition chimique
! de chacun des modes apres equilibre chimique
!---------------------------------------------

ZDEL(:,:)=0.
! Concentration des especes 'totales' presentes dans l'aerosol
! on repartie la masse en fonction de la surface integree des aerosols (omega)
ZDEL(:,1:NSP+NCARB+NSOA)=ZTOTNEW(:,1:NSP+NCARB+NSOA)-ZTOT(:,1:NSP+NCARB+NSOA)

DO II=1,JPMODE
 DO IJ=NSP+NCARB+1,NSP+NCARB+NSOA
  PCTOTA(:,IJ,II)=MAX(1.E-60,PCTOTA(:,IJ,II)+ZDEL(:,IJ)*POM(:,II))
 ENDDO
 PCTOTA(:,JP_AER_H2O,II)=MAX(1.E-60,PCTOTA(:,JP_AER_H2O,II)+ZDEL(:,JP_AER_H2O)*POM(:,II))
ENDDO
 DO IJ=NSP+NCARB+1,NSP+NCARB+NSOA
 PCTOTG(:,IJ)=ZTOTGNEW(:,IJ)
 ENDDO


IF (LHOOK) CALL DR_HOOK('CH_AER_ORGANIC',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_ORGANIC
