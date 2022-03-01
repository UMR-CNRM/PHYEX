MODULE MODE_ICECLOUD
IMPLICIT NONE
CONTAINS
SUBROUTINE ICECLOUD  &
!   Input :
     & ( NP,PP,PZ,PDZ,PT,PR,PTSTEP,PPBLH,PWCLD,XW2D, &
!   Output :
     &  SIFRC,SSIO,SSIU,W2D,RSI)

  USE PARKIND1, ONLY : JPRB
  USE YOMHOOK , ONLY : LHOOK, DR_HOOK
  USE MODD_CST,ONLY : XCPD,XCPV,XLVTT,XLSTT,XG,XRD,XEPSILO
  USE MODE_TIWMX, ONLY: ESATW, ESATI
  USE MODE_QSATMX_TAB
!  USE MODI_TIWMX
  IMPLICIT NONE
!-----------------------------------------------------------------------
!
! Purpose:
! calculate subgridscale fraction of supersaturation with respect to ice.
! Method:
! Assume a linear distubution of relative humidity and let the variability 
! of humidity be a function of model level thickness.
! (Also a function of of humidity itself in the boundary layer)
!     Interface:    subroutine ICECLOUD  is called
!     ------------  from subroutine 'rain_ice'
!
!     variable        type         content
!     ========        ====         =======
!
!     INPUT  arguments  (arguments d'entree)
!----------------------------------------------
!     PP        : pressure at model level (Pa)
!     PZ        : model level height (m)
!     PDZ       : model level thickness (m)
!     PT        : temperature (K)
!     PR        : model level humidity mixing ratio (kg/kg)
!     PTSTEP    : timestep
!     PPBLH     : plantetary layer height (m) (negative value means unknown)
!     PWCLD     : water and / mixed phase cloud cover (negative means unknown)
!     XW2D      : quota between ice crystal concentration between dry and wet
!                 part of a gridbox

!     OUTPUT  arguments  (arguments d'sortie)
!---------------------------------------------
!     SIFRC     : subgridscale fraction with supersaturation with respect to ice.
!     SSIO      : Super-saturation with respect to ice in the  
!                 supersaturated fraction
!     SSIU      : Sub-saturation with respect to ice in the sub-saturated 
!                 fraction
!     W2D       : Factor used to get consistncy between the mean value of 
!                 the gridbox and parts of the gridbox 
!     RSI       : Saturation mixing ratio over ice

INTEGER, INTENT(IN) ::     NP
REAL,  INTENT(IN)  ::      PP(NP)
REAL,  INTENT(IN)  ::      PZ(NP)
REAL,  INTENT(IN)  ::      PDZ(NP)
REAL,  INTENT(IN)  ::      PT(NP)
REAL,  INTENT(IN)  ::      PR(NP)
REAL,  INTENT(IN)  ::      PTSTEP
REAL,  INTENT(IN)  ::      PPBLH
REAL,  INTENT(IN)  ::      PWCLD(NP)
REAL,  INTENT(IN)  ::      XW2D

!     OUTPUT  arguments  (arguments d'sortie)
!---------------------------------------------
REAL,  INTENT(OUT) ::      SIFRC(NP)
REAL,  INTENT(OUT) ::      SSIO(NP)
REAL,  INTENT(OUT) ::      SSIU(NP)
REAL,  INTENT(OUT) ::      W2D(NP)
REAL,  INTENT(OUT) ::      RSI(NP)

!     Working variables:
REAL :: ZSIGMAX,ZSIGMAY,ZSIGMAZ,ZXDIST,ZYDIST,&
     & ZRSW,ZRHW,ZRHIN,ZDRHDZ,ZZ,ZRHDIST,ZRHLIM, &
     & ZRHDIF,ZWCLD,ZI2W,ZRHLIMICE,ZRHLIMINV,ZA,ZRHI,ZR
INTEGER :: JK

!     External function
!REAL :: QSATMX_TAB

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ICECLOUD',0,ZHOOK_HANDLE)

ZSIGMAX=3.E-4         ! assumed rh variation in x axis direction
ZSIGMAY=ZSIGMAX            ! assumed rh variation in y axis direction
ZSIGMAZ=1.6E-2 ! t5

!ZXDIST=DTHETA*110000.
ZXDIST=2500.
   ! gridsize in  x axis (m) Avoid too low
   ! since the model has a tendency to become
   ! drier at high horizontal resolution
   ! due to stronger vertical velocities.
ZYDIST=ZXDIST          ! gridsize in  y axis (m)

DO JK = 1, NP
   ZR = MAX(0.,PR(JK)*PTSTEP)
   SIFRC(JK) = 0.
   ZA = ZR*PP(JK)/(XEPSILO + ZR)
   ZRHW = ZA/ESATW(PT(JK))  
   RSI(JK) = QSATMX_TAB(PP(JK),PT(JK),1.)
   ZRHI = ZA/ESATI(PT(JK))
   ZI2W =  ESATW(PT(JK))/ESATI(PT(JK))

   SSIU(JK) = MIN(ZI2W,ZRHI)
   SSIO(JK) = SSIU(JK)
   W2D(JK) = 1.

   IF (PT(JK)>273.1 .OR. ZR<=0. .OR. ESATI(PT(JK)) >= PP(JK)*0.5) THEN
      SSIU(JK) = SSIU(JK) - 1.
      SSIO(JK) = SSIU(JK)
      IF(PWCLD(JK)>=0.) SIFRC(JK) = PWCLD(JK)
      CYCLE
   ENDIF


   ZRHIN = MAX(0.05, MIN(1.,ZRHW))

   ZDRHDZ=ZRHIN*XG /(PT(JK)*XRD)*  &
        &     ( XEPSILO*XLVTT/(XCPD*PT(JK)) - 1.) ! correct
!              &     ( ZEPSILO*XLSTT/(XCPD*PT) -1.)  ! incorrect 
!          more exact
!          assumed rh variation in the z axis (rh/m) in the pbl .
!          Also possible to just use
!         zdrhdz=4.2e-4_jprb ! rh/m !

   ZZ=0.
   IF(PPBLH < 0. )THEN ! Assume boundary layer height is not available 
      ZZ = MIN(1.,MAX(0.,PZ(JK)*0.001))
   ELSE
      IF(PZ(JK) > 35. .AND. PZ(JK) > PPBLH) ZZ = 1.
   ENDIF

!        1.6e-2 rh/m means variations is of order 0.5 for a 1km dept.
!        sigmaz=4e-2 ! EO 140 lev.


!        Compute rh-variation is x,y,z direction as approxmately
!        independent, exept for the z variation in the pbl, where rh is
!        assumed to be fairly constantly increasing with height

   ZRHDIST = SQRT( ZXDIST*ZSIGMAX**2 + ZYDIST*ZSIGMAY**2 +  &
        &         (1.-ZZ)* (PDZ(JK)*ZDRHDZ)**2 + ZZ*PDZ(JK)*ZSIGMAZ**2)
!         z-variation of rh in the pbl    z-variation of rh outside the pbl
!         Safety for very coarse vertical resolution:
   IF(ZZ > 0.1) ZRHDIST = ZRHDIST/(1.+ZRHDIST)

!!!! Note ZRHDIST is with respect to water ! !!!!!!!!!!!!

   ZRHLIM = MAX(0.5,MIN(0.99,1.-0.5*ZRHDIST))

   IF(PWCLD(JK) < 0.)THEN
   !  Assume water/mixed-phase cloud cover from e.g. 
   ! statistical cloud scheme is not available
      ZRHDIF = (1. - ZRHW)/(1.0-ZRHLIM)
      ZRHDIF =  1. - SQRT(MAX(0.,ZRHDIF))
      ZWCLD = MIN(1.,MAX(ZRHDIF,0.0))
   ELSE
      ZWCLD = PWCLD(JK)
! possible to backwards compute a critical relative humity consitent with 
!  input cloudcover:
!   IF(PWCLD < 0.99 .AND. PWCLD > 0.01) ZRHLIM= 1. - (1.-ZRHW)/(1.-PWCLD)**2
   ENDIF

   SIFRC(JK) = ZWCLD

!              relation rhlim with respect to water to that of ice:
!ZRHLIMICE = MAX(ZRHDMIN*ZI2W,1.+ ZI2W*( ZRHLIM - 1.))
   ZRHLIMICE = 1.+ ZI2W*( ZRHLIM - 1.)

   IF(ZRHLIM <= 0.999)THEN

   !              compute a 1/(1-rhlim) constistant with  lstmp(i,k):
      ZRHLIMINV = 1./(1. - ZRHLIMICE)
      ZRHDIF = (ZRHI - ZRHLIMICE)*ZRHLIMINV

      IF(ZWCLD==0.)THEN
         SIFRC(JK) = MIN(1.,0.5*MAX(0.,ZRHDIF))
      ELSE
         ZA = 1. -  1./ZI2W
         SIFRC(JK) = MIN(1.,ZA*0.5/ (1. - ZRHLIM))
         SIFRC(JK) = MIN(1.,ZWCLD + SIFRC(JK))
      ENDIF
   ENDIF

   IF(SIFRC(JK) > 0.01) THEN
      SSIU(JK) = SIFRC(JK) + ZRHLIMICE*(1.-SIFRC(JK))
      SSIO(JK) = (ZRHI - (1.- SIFRC(JK))*SSIU(JK))/SIFRC(JK)
   ELSE
      SIFRC(JK) = 0.! to aviod mismatch with output variables
      ZA = MIN(0.,ZRHI-ZRHLIMICE)
      SSIU(JK) = MAX(0.,SIFRC(JK) + ZRHLIMICE*(1.-SIFRC(JK)) + 2.*ZA )
   ENDIF
   SSIO(JK) = MIN(ZI2W,SSIO(JK))
   SSIU(JK) = MAX(0.,SSIU(JK))

! Transform from relative humidity to degree of saturation:
   SSIU(JK) = SSIU(JK) - 1.
   SSIO(JK) = SSIO(JK) - 1.

   IF (XW2D > 1.) W2D(JK) = 1./(1. - SIFRC(JK) + XW2D*SIFRC(JK))

ENDDO

IF (LHOOK) CALL DR_HOOK('ICECLOUD',1,ZHOOK_HANDLE)
END SUBROUTINE ICECLOUD
END MODULE MODE_ICECLOUD
