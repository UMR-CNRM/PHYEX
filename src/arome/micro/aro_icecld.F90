SUBROUTINE ARO_ICECLD  &
!   Input :
     & ( PP,PZ,PDZ,PT,PR,PPBLH,PWCLD,XW2D, &
!   Output :
     &  SIFRC,SSIO,SSIU,W2D,RSI)


  USE PARKIND1, ONLY : JPRB
  USE YOMHOOK , ONLY : LHOOK, DR_HOOK
  USE MODD_CST,ONLY : XCPD,XCPV,XLVTT,XLSTT,XG,XRD,XTT,XMD,XMV,XEPSILO
  IMPLICIT NONE
!-----------------------------------------------------------------------
!
! Purpose:
! calculate subgridscale fraction of supersaturation with respect to ice.
! Method:
! Assume a linear distubution of relative humidity and let the variability 
! of humidity be a function of model level thickness.
! (Also a function of of humidity itself in the boundary layer)
!     Interface:    subroutine ARO_ICECLD  is called
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




REAL,  INTENT(IN)  ::      PP
REAL,  INTENT(IN)  ::      PZ
REAL,  INTENT(IN)  ::      PDZ
REAL,  INTENT(IN)  ::      PT
REAL,  INTENT(IN)  ::      PR
REAL,  INTENT(IN)  ::      PPBLH
REAL,  INTENT(IN)  ::      PWCLD
REAL,  INTENT(IN)  ::      XW2D

!     OUTPUT  arguments  (arguments d'sortie)
!---------------------------------------------
REAL,  INTENT(OUT) ::      SIFRC
REAL,  INTENT(OUT) ::      SSIO
REAL,  INTENT(OUT) ::      SSIU
REAL,  INTENT(OUT) ::      W2D
REAL,  INTENT(OUT) ::      RSI
!     Working variables:

 REAL  :: ZSIGMAX,ZSIGMAY,ZSIGMAZ,ZFICE,ZXDIST, ZYDIST,&
      & ZRSW,ZRHW,ZRHIN,ZDRHDZ,ZZ,ZRHDIST ,ZRHLIM, &
      & ZRHDIF,ZWCLD ,ZI2W,ZRHLIMICE,ZRHLIMINV,ZA,ZRHI

!     ==================================================================
!     1. Declarations.
!     ==================================================================
!     1.1 MODULES USED
!-----------------------------------------------------------------------
 REAL AROQSATMX

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_ICECLD',0,ZHOOK_HANDLE)

SIFRC = 0.
ZFICE = 0.
ZRSW= AROQSATMX(PP,PT,ZFICE)
ZRHW= PR/ZRSW
ZFICE=1.
RSI= AROQSATMX(PP,PT,ZFICE)
ZRHI= PR/RSI
ZI2W = ZRSW/RSI
SSIU=ZRHI
SSIO=SSIU
W2D  = 1.

IF(PT>XTT.OR. PR<=0.)THEN
   SSIU=SSIU-1
   SSIO=SSIU
   IF(PWCLD>=0.)SIFRC=PWCLD
   IF (LHOOK) CALL DR_HOOK('ARO_ICECLD',1,ZHOOK_HANDLE)
   RETURN
ENDIF

ZSIGMAX=3.E-4         ! assumed rh variation in x axis direction
ZFICE=0.                    ! fraction of ice
ZSIGMAY=ZSIGMAX            ! assumed rh variation in y axis direction
ZSIGMAZ=1.6E-2 ! t5

!ZXDIST=MAX(0.10,DTHETA)*110000.
!ZXDIST=DTHETA*110000.
ZXDIST=2500. 
   ! gridsize in  x axis (m) Avoid too low
   ! since the model has a tendency to become
   ! drier at high horizontal resolution
   ! due to stronger vertical velocities.
ZYDIST=ZXDIST          ! gridsize in  y axis (m)




ZRHIN = MAX(0.05, MIN(1.,PR/ZRSW))

ZDRHDZ=ZRHIN*XG /(PT*XRD)*  &
              &     ( XEPSILO*XLVTT/(XCPD*PT) -1.) ! correct
!              &     ( XEPSILO*XLSTT/(XCPD*PT) -1.)  ! incorrect but currently used
!          more exact
!          assumed rh variation in the z axis (rh/m) in the pbl .
!          Also possible to just use
!         zdrhdz=4.2e-4_jprb ! rh/m !
ZZ=0.


IF(PPBLH < 0. )THEN ! Assume boundary layer height is not available 
   ZZ= MIN(1.,MAX(0.,PZ*0.001))
ELSE
   IF(PZ > 35. .AND. PZ > PPBLH)ZZ= 1.
ENDIF

!        1.6e-2 rh/m means variations is of order 0.5 for a 1km dept.
!        sigmaz=4e-2 ! EO 140 lev.


!        Compute rh-variation is x,y,z direction as approxmately
!        independent, exept for the z variation in the pbl, where rh is
!        assumed to be fairly constantly increasing with height
if(ZXDIST*ZSIGMAX**2 + ZYDIST*ZSIGMAY**2 +  &
     &           (1.-ZZ)* (PDZ* ZDRHDZ)**2 + ZZ*PDZ* ZSIGMAZ**2 < 0.)then
write(*,*)'in ARO_ICECLD: PDZ ZDRHDZ,expression=',PDZ, ZDRHDZ&
& ,ZXDIST*ZSIGMAX**2 + ZYDIST*ZSIGMAY**2 +  &
     &           (1.-ZZ)* (PDZ* ZDRHDZ)**2 + ZZ*PDZ* ZSIGMAZ**2 
endif
ZRHDIST = SQRT( ZXDIST*ZSIGMAX**2 + ZYDIST*ZSIGMAY**2 +  &
     &           (1.-ZZ)* (PDZ* ZDRHDZ)**2 + ZZ*PDZ* ZSIGMAZ**2)
!         z-variation of rh in the pbl    z-variation of rh outside the pbl
!         Safety for very coarse vertical resolution:
IF(ZZ > 0.1) ZRHDIST = ZRHDIST/(1.+ZRHDIST)

!!!! Note ZRHDIST is with respect to water ! !!!!!!!!!!!!


ZRHLIM = MAX(0.5, MIN(0.99,1. - 0.5*ZRHDIST))


IF(PWCLD < 0.)THEN
   !  Assume water/mixed-phase cloud cover from e.g. 
   ! statistical cloud scheme is not avialabe
   ZRHDIF = (1. - ZRHW)/(1.0-ZRHLIM)
   ZRHDIF =  1. - SQRT(MAX(0.,ZRHDIF))
   ZWCLD = MIN(1.,MAX(ZRHDIF,0.0))
ELSE
   ZWCLD = PWCLD
! possible to backwards compute a critical relative humity consitent with 
!  input cloudcover:
!   IF(PWCLD < 0.99 .AND. PWCLD > 0.01) ZRHLIM= 1. - (1.-ZRHW)/(1.-PWCLD)**2
ENDIF

SIFRC = ZWCLD

!              relation rhlim with respect to water to that of ice:
!ZRHLIMICE = MAX(ZRHDMIN*ZI2W,1.+ ZI2W*( ZRHLIM - 1.))
ZRHLIMICE = 1.+ ZI2W*( ZRHLIM - 1.)

IF(ZRHLIM <= 0.999)THEN

   !              compute a 1/(1-rhlim) constistant with  lstmp(i,k):
   ZRHLIMINV = 1./(1. - ZRHLIMICE)
   ZRHDIF = (ZRHI - ZRHLIMICE)*ZRHLIMINV

   IF(ZWCLD==0.)THEN
      SIFRC = MIN(1.,0.5*MAX(0.,ZRHDIF))
   ELSE
      ZA =1. -  1./ZI2W
      SIFRC =MIN(1.,ZA*0.5/ (1. - ZRHLIM))
      SIFRC = MIN(1.,ZWCLD + SIFRC)
   ENDIF
ENDIF



IF(SIFRC > 0.01) THEN
   SSIU = SIFRC + ZRHLIMICE*(1.-SIFRC)
   SSIO = (ZRHI - (1.- SIFRC)*SSIU)/SIFRC
ELSE
   SIFRC=0.! to aviod mismatch with output variables
   ZA = MIN(0.,ZRHI-ZRHLIMICE)
   SSIU = MAX(0.,SIFRC + ZRHLIMICE*(1.-SIFRC) + 2*ZA )
ENDIF


! Transform from relative humidity to degree of saturation:
SSIU = SSIU - 1.
SSIO = SSIO - 1.

IF (XW2D > 1.) W2D = 1./(1. - SIFRC + XW2D*SIFRC)

IF (LHOOK) CALL DR_HOOK('ARO_ICECLD',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_ICECLD
