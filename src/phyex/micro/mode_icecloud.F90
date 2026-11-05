MODULE MODE_ICECLOUD
IMPLICIT NONE
CONTAINS
SUBROUTINE ICECLOUD  &
!   Input :
     & ( D,CST,ICEP,PP,PZ,PDZ,PT,PR,PTSTEP,PPBLH,PWCLD,XW2D, &
!   Output :
     &  SIFRC,SSIO,SSIU,W2D,RSI)

  USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
  USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
  USE MODD_CST,ONLY : CST_t
  USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
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
TYPE(CST_t),      INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t), INTENT(IN) :: ICEP
REAL,  INTENT(IN)  ::      PP(D%NIJT,D%NKT)
REAL,  INTENT(IN)  ::      PZ(D%NIJT,D%NKT)
REAL,  INTENT(IN)  ::      PDZ(D%NIJT,D%NKT)
REAL,  INTENT(IN)  ::      PT(D%NIJT,D%NKT)
REAL,  INTENT(IN)  ::      PR(D%NIJT,D%NKT)
REAL,  INTENT(IN)  ::      PTSTEP
REAL,  INTENT(IN)  ::      PPBLH
REAL,  INTENT(IN)  ::      PWCLD(D%NIJT,D%NKT)
REAL,  INTENT(IN)  ::      XW2D

!     OUTPUT  arguments  (arguments d'sortie)
!---------------------------------------------
REAL,  INTENT(OUT) ::      SIFRC(D%NIJT,D%NKT)
REAL,  INTENT(OUT) ::      SSIO(D%NIJT,D%NKT)
REAL,  INTENT(OUT) ::      SSIU(D%NIJT,D%NKT)
REAL,  INTENT(OUT) ::      W2D(D%NIJT,D%NKT)
REAL,  INTENT(OUT) ::      RSI(D%NIJT,D%NKT)

!     Working variables:
REAL :: ZSIGMAX,ZSIGMAY,ZSIGMAZ,ZXDIST,ZYDIST,&
     & ZRHW,ZRHIN,ZDRHDZ,ZZ,ZRHDIST,ZRHLIM, &
     & ZRHDIF,ZWCLD,ZI2W,ZRHLIMICE,ZRHLIMINV,ZA,ZRHI,ZR
INTEGER :: JIJ, JK, IIJB, IIJE, IKTB, IKTE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ICECLOUD',0,ZHOOK_HANDLE)
!
IIJB=D%NIJB
IIJE=D%NIJE
IKTB=D%NKTB
IKTE=D%NKTE
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
DO JK = IKTB, IKTE
 DO JIJ = IIJB, IIJE
   ZR = MAX(0.,PR(JIJ,JK)*PTSTEP)
   SIFRC(JIJ,JK) = 0.
   ZA = ZR*PP(JIJ,JK)/(CST%XEPSILO + ZR)
   ZRHW = ZA/ESATW(ICEP%TIWMX, PT(JIJ,JK))  
   RSI(JIJ,JK) = QSATMX_TAB(CST,ICEP%TIWMX,PP(JIJ,JK),PT(JIJ,JK),1.)
   ZRHI = ZA/ESATI(ICEP%TIWMX, PT(JIJ,JK))
   ZI2W =  ESATW(ICEP%TIWMX, PT(JIJ,JK))/ESATI(ICEP%TIWMX, PT(JIJ,JK))

   SSIU(JIJ,JK) = MIN(ZI2W,ZRHI)
   SSIO(JIJ,JK) = SSIU(JIJ,JK)
   W2D(JIJ,JK) = 1.
   IF (PT(JIJ,JK)>273.1 .OR. ZR<=0. .OR. ESATI(ICEP%TIWMX, PT(JIJ,JK)) >= PP(JIJ,JK)*0.5) THEN
      SSIU(JIJ,JK) = SSIU(JIJ,JK) - 1.
      SSIO(JIJ,JK) = SSIU(JIJ,JK)
      IF(PWCLD(JIJ,JK)>=0.) SIFRC(JIJ,JK) = PWCLD(JIJ,JK)
   ELSE

   ZRHIN = MAX(0.05, MIN(1.,ZRHW))

   ZDRHDZ=ZRHIN*CST%XG /(PT(JIJ,JK)*CST%XRD)*  &
        &     ( CST%XEPSILO*CST%XLVTT/(CST%XCPD*PT(JIJ,JK)) - 1.) ! correct
!              &     ( ZEPSILO*XLSTT/(CST%XCPD*PT) -1.)  ! incorrect 
!          more exact
!          assumed rh variation in the z axis (rh/m) in the pbl .
!          Also possible to just use
!         zdrhdz=4.2e-4_jprb ! rh/m !

   ZZ=0.
   IF(PPBLH < 0. )THEN ! Assume boundary layer height is not available 
      ZZ = MIN(1.,MAX(0.,PZ(JIJ,JK)*0.001))
   ELSE
      IF(PZ(JIJ,JK) > 35. .AND. PZ(JIJ,JK) > PPBLH) ZZ = 1.
   ENDIF

!        1.6e-2 rh/m means variations is of order 0.5 for a 1km dept.
!        sigmaz=4e-2 ! EO 140 lev.


!        Compute rh-variation is x,y,z direction as approxmately
!        independent, exept for the z variation in the pbl, where rh is
!        assumed to be fairly constantly increasing with height

   ZRHDIST = SQRT( ZXDIST*ZSIGMAX**2 + ZYDIST*ZSIGMAY**2 +  &
        &         (1.-ZZ)* (PDZ(JIJ,JK)*ZDRHDZ)**2 + ZZ*PDZ(JIJ,JK)*ZSIGMAZ**2)
!         z-variation of rh in the pbl    z-variation of rh outside the pbl
!         Safety for very coarse vertical resolution:
   IF(ZZ > 0.1) ZRHDIST = ZRHDIST/(1.+ZRHDIST)

!!!! Note ZRHDIST is with respect to water ! !!!!!!!!!!!!

   ZRHLIM = MAX(0.5,MIN(0.99,1.-0.5*ZRHDIST))

   IF(PWCLD(JIJ,JK) < 0.)THEN
   !  Assume water/mixed-phase cloud cover from e.g. 
   ! statistical cloud scheme is not available
      ZRHDIF = (1. - ZRHW)/(1.0-ZRHLIM)
      ZRHDIF =  1. - SQRT(MAX(0.,ZRHDIF))
      ZWCLD = MIN(1.,MAX(ZRHDIF,0.0))
   ELSE
      ZWCLD = PWCLD(JIJ,JK)
! possible to backwards compute a critical relative humity consitent with 
!  input cloudcover:
!   IF(PWCLD < 0.99 .AND. PWCLD > 0.01) ZRHLIM= 1. - (1.-ZRHW)/(1.-PWCLD)**2
   ENDIF

   SIFRC(JIJ,JK) = ZWCLD

!              relation rhlim with respect to water to that of ice:
!ZRHLIMICE = MAX(ZRHDMIN*ZI2W,1.+ ZI2W*( ZRHLIM - 1.))
   ZRHLIMICE = 1.+ ZI2W*( ZRHLIM - 1.)

   IF(ZRHLIM <= 0.999)THEN

   !              compute a 1/(1-rhlim) constistant with  lstmp(i,k):
      ZRHLIMINV = 1./(1. - ZRHLIMICE)
      ZRHDIF = (ZRHI - ZRHLIMICE)*ZRHLIMINV

      IF(ZWCLD==0.)THEN
         SIFRC(JIJ,JK) = MIN(1.,0.5*MAX(0.,ZRHDIF))
      ELSE
         ZA = 1. -  1./ZI2W
         SIFRC(JIJ,JK) = MIN(1.,ZA*0.5/ (1. - ZRHLIM))
         SIFRC(JIJ,JK) = MIN(1.,ZWCLD + SIFRC(JIJ,JK))
      ENDIF
   ENDIF

   IF(SIFRC(JIJ,JK) > 0.01) THEN
      SSIU(JIJ,JK) = SIFRC(JIJ,JK) + ZRHLIMICE*(1.-SIFRC(JIJ,JK))
      SSIO(JIJ,JK) = (ZRHI - (1.- SIFRC(JIJ,JK))*SSIU(JIJ,JK))/SIFRC(JIJ,JK)
   ELSE
      SIFRC(JIJ,JK) = 0.! to aviod mismatch with output variables
      ZA = MIN(0.,ZRHI-ZRHLIMICE)
      SSIU(JIJ,JK) = MAX(0.,SIFRC(JIJ,JK) + ZRHLIMICE*(1.-SIFRC(JIJ,JK)) + 2.*ZA )
   ENDIF
   SSIO(JIJ,JK) = MIN(ZI2W,SSIO(JIJ,JK))
   SSIU(JIJ,JK) = MAX(0.,SSIU(JIJ,JK))

! Transform from relative humidity to degree of saturation:
   SSIU(JIJ,JK) = SSIU(JIJ,JK) - 1.
   SSIO(JIJ,JK) = SSIO(JIJ,JK) - 1.

   IF (XW2D > 1.) W2D(JIJ,JK) = 1./(1. - SIFRC(JIJ,JK) + XW2D*SIFRC(JIJ,JK))
 END IF
 ENDDO
ENDDO

IF (LHOOK) CALL DR_HOOK('ICECLOUD',1,ZHOOK_HANDLE)
END SUBROUTINE ICECLOUD
END MODULE MODE_ICECLOUD
