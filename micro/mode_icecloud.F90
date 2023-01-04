MODULE MODE_ICECLOUD
IMPLICIT NONE
CONTAINS
SUBROUTINE ICECLOUD  &
!   Input :
     & ( D,PP,PZ,PDZ,PT,PR,PTSTEP,PPBLH,PWCLD,XW2D, &
!   Output :
     &  SIFRC,SSIO,SSIU,W2D,RSI)

  USE PARKIND1, ONLY : JPRB
  USE YOMHOOK , ONLY : LHOOK, DR_HOOK
  USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
  USE MODD_CST,ONLY : XCPD,XCPV,XLVTT,XLSTT,XG,XRD,XEPSILO
  USE MODE_TIWMX, ONLY: ESATW, ESATI
  USE MODE_QSATMX_TAB, ONLY: QSATMX_TAB
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

TYPE(DIMPHYEX_t), INTENT(IN)    :: D
REAL,  INTENT(IN)  ::      PP(D%NIJT)
REAL,  INTENT(IN)  ::      PZ(D%NIJT)
REAL,  INTENT(IN)  ::      PDZ(D%NIJT)
REAL,  INTENT(IN)  ::      PT(D%NIJT)
REAL,  INTENT(IN)  ::      PR(D%NIJT)
REAL,  INTENT(IN)  ::      PTSTEP
REAL,  INTENT(IN)  ::      PPBLH
REAL,  INTENT(IN)  ::      PWCLD(D%NIJT)
REAL,  INTENT(IN)  ::      XW2D

!     OUTPUT  arguments  (arguments d'sortie)
!---------------------------------------------
REAL,  INTENT(OUT) ::      SIFRC(D%NIJT)
REAL,  INTENT(OUT) ::      SSIO(D%NIJT)
REAL,  INTENT(OUT) ::      SSIU(D%NIJT)
REAL,  INTENT(OUT) ::      W2D(D%NIJT)
REAL,  INTENT(OUT) ::      RSI(D%NIJT)

!     Working variables:
REAL :: ZSIGMAX,ZSIGMAY,ZSIGMAZ,ZXDIST,ZYDIST,&
     & ZRSW,ZRHW,ZRHIN,ZDRHDZ,ZZ,ZRHDIST,ZRHLIM, &
     & ZRHDIF,ZWCLD,ZI2W,ZRHLIMICE,ZRHLIMINV,ZA,ZRHI,ZR
INTEGER :: JIJ, IIJB, IIJE

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ICECLOUD',0,ZHOOK_HANDLE)
!
IIJB=D%NIJB
IIJE=D%NIJE
!
ZSIGMAX=3.E-4         ! assumed rh variation in x axis direction
ZSIGMAY=ZSIGMAX            ! assumed rh variation in y axis direction
ZSIGMAZ=1.E-2

!ZXDIST=DTHETA*110000.
ZXDIST=2500.
   ! gridsize in  x axis (m) Avoid too low
   ! since the model has a tendency to become
   ! drier at high horizontal resolution
   ! due to stronger vertical velocities.
ZYDIST=ZXDIST          ! gridsize in  y axis (m)

DO JIJ = IIJB, IIJE
   ZR = MAX(0.,PR(JIJ)*PTSTEP)
   SIFRC(JIJ) = 0.
   ZA = ZR*PP(JIJ)/(XEPSILO + ZR)
   ZRHW = ZA/ESATW(PT(JIJ))  
   RSI(JIJ) = QSATMX_TAB(PP(JIJ),PT(JIJ),1.)
   ZRHI = ZA/ESATI(PT(JIJ))
   ZI2W =  ESATW(PT(JIJ))/ESATI(PT(JIJ))

   SSIU(JIJ) = MIN(ZI2W,ZRHI)
   SSIO(JIJ) = SSIU(JIJ)
   W2D(JIJ) = 1.

   IF (PT(JIJ)>273.1 .OR. ZR<=0. .OR. ESATI(PT(JIJ)) >= PP(JIJ)*0.5) THEN
      SSIU(JIJ) = SSIU(JIJ) - 1.
      SSIO(JIJ) = SSIU(JIJ)
      IF(PWCLD(JIJ)>=0.) SIFRC(JIJ) = PWCLD(JIJ)
      CYCLE
   ENDIF


   ZRHIN = MAX(0.05, MIN(1.,ZRHW))

   ZDRHDZ=ZRHIN*XG /(PT(JIJ)*XRD)*  &
        &     ( XEPSILO*XLVTT/(XCPD*PT(JIJ)) - 1.) ! correct
!              &     ( ZEPSILO*XLSTT/(XCPD*PT) -1.)  ! incorrect 
!          more exact
!          assumed rh variation in the z axis (rh/m) in the pbl .
!          Also possible to just use
!         zdrhdz=4.2e-4_jprb ! rh/m !

   ZZ=0.
   IF(PPBLH < 0. )THEN ! Assume boundary layer height is not available 
      ZZ = MIN(1.,MAX(0.,PZ(JIJ)*0.001))
   ELSE
      IF(PZ(JIJ) > 35. .AND. PZ(JIJ) > PPBLH) ZZ = 1.
   ENDIF

!        1.6e-2 rh/m means variations is of order 0.5 for a 1km dept.
!        sigmaz=4e-2 ! EO 140 lev.


!        Compute rh-variation is x,y,z direction as approxmately
!        independent, exept for the z variation in the pbl, where rh is
!        assumed to be fairly constantly increasing with height

   ZRHDIST = SQRT( ZXDIST*ZSIGMAX**2 + ZYDIST*ZSIGMAY**2 +  &
        &         (1.-ZZ)* (PDZ(JIJ)*ZDRHDZ)**2 + ZZ*PDZ(JIJ)*ZSIGMAZ**2)
!         z-variation of rh in the pbl    z-variation of rh outside the pbl
!         Safety for very coarse vertical resolution:
   IF(ZZ > 0.1) ZRHDIST = ZRHDIST/(1.+ZRHDIST)

!!!! Note ZRHDIST is with respect to water ! !!!!!!!!!!!!

   ZRHLIM = MAX(0.5,MIN(0.99,1.-0.5*ZRHDIST))

   IF(PWCLD(JIJ) < 0.)THEN
   !  Assume water/mixed-phase cloud cover from e.g. 
   ! statistical cloud scheme is not available
      ZRHDIF = (1. - ZRHW)/(1.0-ZRHLIM)
      ZRHDIF =  1. - SQRT(MAX(0.,ZRHDIF))
      ZWCLD = MIN(1.,MAX(ZRHDIF,0.0))
   ELSE
      ZWCLD = PWCLD(JIJ)
! possible to backwards compute a critical relative humity consitent with 
!  input cloudcover:
!   IF(PWCLD < 0.99 .AND. PWCLD > 0.01) ZRHLIM= 1. - (1.-ZRHW)/(1.-PWCLD)**2
   ENDIF

   SIFRC(JIJ) = ZWCLD

!              relation rhlim with respect to water to that of ice:
!ZRHLIMICE = MAX(ZRHDMIN*ZI2W,1.+ ZI2W*( ZRHLIM - 1.))
   ZRHLIMICE = 1.+ ZI2W*( ZRHLIM - 1.)

   IF(ZRHLIM <= 0.999)THEN

   !              compute a 1/(1-rhlim) constistant with  lstmp(i,k):
      ZRHLIMINV = 1./(1. - ZRHLIMICE)
      ZRHDIF = (ZRHI - ZRHLIMICE)*ZRHLIMINV

      IF(ZWCLD==0.)THEN
         SIFRC(JIJ) = MIN(1.,0.5*MAX(0.,ZRHDIF))
      ELSE
         ZA = 1. -  1./ZI2W
         SIFRC(JIJ) = MIN(1.,ZA*0.5/ (1. - ZRHLIM))
         SIFRC(JIJ) = MIN(1.,ZWCLD + SIFRC(JIJ))
      ENDIF
   ENDIF

   IF(SIFRC(JIJ) > 0.01) THEN
      SSIU(JIJ) = SIFRC(JIJ) + ZRHLIMICE*(1.-SIFRC(JIJ))
      SSIO(JIJ) = (ZRHI - (1.- SIFRC(JIJ))*SSIU(JIJ))/SIFRC(JIJ)
   ELSE
      SIFRC(JIJ) = 0.! to aviod mismatch with output variables
      ZA = MIN(0.,ZRHI-ZRHLIMICE)
      SSIU(JIJ) = MAX(0.,SIFRC(JIJ) + ZRHLIMICE*(1.-SIFRC(JIJ)) + 2.*ZA )
   ENDIF
   SSIO(JIJ) = MIN(ZI2W,SSIO(JIJ))
   SSIU(JIJ) = MAX(0.,SSIU(JIJ))

! Transform from relative humidity to degree of saturation:
   SSIU(JIJ) = SSIU(JIJ) - 1.
   SSIO(JIJ) = SSIO(JIJ) - 1.

   IF (XW2D > 1.) W2D(JIJ) = 1./(1. - SIFRC(JIJ) + XW2D*SIFRC(JIJ))

ENDDO

IF (LHOOK) CALL DR_HOOK('ICECLOUD',1,ZHOOK_HANDLE)
END SUBROUTINE ICECLOUD
END MODULE MODE_ICECLOUD
