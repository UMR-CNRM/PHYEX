!     ######spl
          SUBROUTINE COMPUTE_ENTR_DETR(KK,KKB,KKE,KKL,OTEST,OTESTLCL,&
                            HFRAC_ICE,PFRAC_ICE,PRHODREF,&
                            PPRE_MINUS_HALF,&
                            PPRE_PLUS_HALF,PZZ,PDZZ,&
                            PTHVM,PTHLM,PRTM,PW_UP2,PTH_UP,&
                            PTHL_UP,PRT_UP,PLUP,&
                            PRC_UP,PRI_UP,PTHV_UP,&
                            PRSAT_UP,PRC_MIX,PRI_MIX,      &
                            PENTR,PDETR,PENTR_CLD,PDETR_CLD,&
                            PBUO_INTEG_DRY,PBUO_INTEG_CLD,&
                            PPART_DRY)

          USE PARKIND1, ONLY : JPRB
          USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!         #############################################################

!!
!!***COMPUTE_ENTR_DETR* - calculates caracteristics of the updraft or downdraft
!!                       using model of the EDMF scheme 
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to compute entrainement and
!!      detrainement at one level of the updraft
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     REFERENCE
!!     ---------
!!       Book 1 of Meso-NH documentation (chapter Convection)
!!       
!!
!!     AUTHOR
!!     ------
!!    J.Pergaud : 2009
!!
!!    MODIFICATIONS
!!    -------------
!!      Y.Seity (06/2010) Bug correction
!!      V.Masson (09/2010) Optimization
!!      S. Riette april 2011 : ice added, protection against zero divide by Yves Bouteloup
!!                             protection against too big ZPART_DRY, interface modified
!!      S. Riette Jan 2012: support for both order of vertical levels
!!      P.Marguinaud Jun 2012: fix uninitialized variable
!!      P.Marguinaud Nov 2012: fix gfortran bug
!!      S. Riette Apr 2013: bugs correction, rewriting (for optimisation) and
!!                          improvement of continuity at the condensation level
!!      S. Riette Nov 2013: protection against zero divide for min value of dry PDETR
!!      R. El Khatib 29-Apr-2019 portability fix : compiler may get confused by embricked WHERE statements
!!                          eventually breaking tests with NaN initializations at compile time.
!!                          Replace by IF conditions and traditional DO loops can only improve the performance.
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!                         
USE MODD_CST
!
USE MODD_CMFSHALL
!
USE MODI_TH_R_FROM_THL_RT_1D 

USE MODE_THERMO

IMPLICIT NONE
!
!                         
!*                    1.1  Declaration of Arguments
!
!
INTEGER,                INTENT(IN)   :: KK
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
LOGICAL,DIMENSION(:),   INTENT(IN)   :: OTEST ! test to see if updraft is running
LOGICAL,DIMENSION(:),   INTENT(IN)   :: OTESTLCL !test of condensation 
CHARACTER*1,            INTENT(IN)   :: HFRAC_ICE ! frac_ice can be compute using
                                              ! Temperature (T) or prescribed
                                              ! (Y)
REAL, DIMENSION(:), INTENT(IN)      :: PFRAC_ICE ! fraction of ice
!
!    prognostic variables at t- deltat
!
REAL, DIMENSION(:),     INTENT(IN) ::  PRHODREF  !rhodref
REAL, DIMENSION(:),     INTENT(IN) ::  PPRE_MINUS_HALF ! Pressure at flux level KK
REAL, DIMENSION(:),     INTENT(IN) ::  PPRE_PLUS_HALF ! Pressure at flux level KK+KKL
REAL, DIMENSION(:,:),   INTENT(IN) ::  PZZ       !  Height at the flux point
REAL, DIMENSION(:,:),   INTENT(IN) ::  PDZZ       !  metrics coefficient
REAL, DIMENSION(:,:),   INTENT(IN) ::  PTHVM      ! ThetaV environment 

!
!   thermodynamical variables which are transformed in conservative var.
!
REAL, DIMENSION(:,:), INTENT(IN)     ::  PTHLM     ! Thetal
REAL, DIMENSION(:,:), INTENT(IN)     ::  PRTM      ! total mixing ratio 
REAL, DIMENSION(:,:), INTENT(IN)     ::  PW_UP2    ! Vertical velocity^2
REAL, DIMENSION(:),   INTENT(IN)     ::  PTH_UP,PTHL_UP,PRT_UP  ! updraft properties
REAL, DIMENSION(:),   INTENT(IN)     ::  PLUP      ! LUP compute from the ground
REAL, DIMENSION(:),   INTENT(IN)     ::  PRC_UP,PRI_UP   ! Updraft cloud content
REAL, DIMENSION(:),   INTENT(IN)     ::  PTHV_UP ! Thetav of updraft
REAL, DIMENSION(:),   INTENT(IN)     ::  PRSAT_UP ! Mixing ratio at saturation in updraft
REAL, DIMENSION(:),   INTENT(INOUT)  ::  PRC_MIX, PRI_MIX      ! Mixture cloud content
REAL, DIMENSION(:),   INTENT(OUT)    ::  PENTR     ! Mass flux entrainment of the updraft
REAL, DIMENSION(:),   INTENT(OUT)    ::  PDETR     ! Mass flux detrainment of the updraft
REAL, DIMENSION(:),   INTENT(OUT)    ::  PENTR_CLD ! Mass flux entrainment of the updraft in cloudy part
REAL, DIMENSION(:),   INTENT(OUT)    ::  PDETR_CLD ! Mass flux detrainment of the updraft in cloudy part
REAL, DIMENSION(:),   INTENT(OUT)    ::  PBUO_INTEG_DRY, PBUO_INTEG_CLD! Integral Buoyancy
REAL, DIMENSION(:),   INTENT(OUT)    ::  PPART_DRY ! ratio of dry part at the transition level
!
!
!                       1.2  Declaration of local variables
!
!

! Variables for cloudy part
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZKIC, ZKIC_F2  ! fraction of env. mass in the muxtures
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZEPSI,ZDELTA   ! factor entrainment detrainment
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZEPSI_CLOUD    ! factor entrainment detrainment
REAL                           :: ZCOEFFMF_CLOUD ! factor for compputing entr. detr.
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZMIXTHL,ZMIXRT ! Thetal and rt in the mixtures
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZTHMIX         ! Theta and Thetav  of mixtures
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZRVMIX,ZRCMIX,ZRIMIX ! mixing ratios in mixtures
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZTHVMIX, ZTHVMIX_F2 ! Theta and Thetav  of mixtures
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZTHV_UP_F2     ! thv_up at flux point kk+kkl
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZRSATW, ZRSATI ! working arrays (mixing ratio at saturation)
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZTHV           ! theta V of environment at the bottom of cloudy part  
REAL                           :: ZKIC_INIT      !Initial value of ZKIC
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZCOTHVU              ! Variation of Thvup between bottom and top of cloudy part

! Variables for dry part
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZFOESW, ZFOESI       ! saturating vapor pressure
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZDRSATODP            ! d.Rsat/dP
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZT                   ! Temperature
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZWK                  ! Work array

! Variables for dry and cloudy parts
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZCOEFF_MINUS_HALF,&  ! Variation of Thv between mass points kk-kkl and kk
                                  ZCOEFF_PLUS_HALF     ! Variation of Thv between mass points kk and kk+kkl
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZPRE                 ! pressure at the bottom of the cloudy part
REAL, DIMENSION(SIZE(PTHVM,1)) :: ZG_O_THVREF
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZFRAC_ICE            ! fraction of ice
REAL                           :: ZRVORD               ! RV/RD
REAL, DIMENSION(SIZE(PTHLM,1)) :: ZDZ_STOP,&           ! Exact Height of the LCL above flux level KK
                                  ZTHV_MINUS_HALF,&    ! Thv at flux point(kk)  
                                  ZTHV_PLUS_HALF,&     ! Thv at flux point(kk+kkl)
                                  ZDZ                  ! Delta Z used in computations
INTEGER :: JI

!----------------------------------------------------------------------------------
                        
!                1.3 Initialisation
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('COMPUTE_ENTR_DETR',0,ZHOOK_HANDLE)

  
  ZRVORD   = XRV / XRD   !=1.607
  ZG_O_THVREF(:)=XG/PTHVM(:,KK)
  ZCOEFFMF_CLOUD=XENTR_MF * XG / XCRAD_MF
  
  ZFRAC_ICE(:)=PFRAC_ICE(:) ! to not modify fraction of ice
 
  ZPRE(:)=PPRE_MINUS_HALF(:)
  ZMIXTHL(:)=0.1
  ZMIXRT(:)=0.1

!                1.4 Estimation of PPART_DRY
  WHERE(OTEST)
    WHERE(OTESTLCL)
      !No dry part when condensation level is reached
      PPART_DRY(:)=0.
      ZDZ_STOP(:)=0.
      ZPRE(:)=PPRE_MINUS_HALF(:)
    ELSEWHERE
      !Temperature at flux level KK
      ZT(:)=PTH_UP(:)*(PPRE_MINUS_HALF(:)/XP00) ** (XRD/XCPD)
      !Saturating vapor pressure at flux level KK
      ZFOESW(:) = MIN(EXP( XALPW - XBETAW/ZT(:) - XGAMW*LOG(ZT(:))  ), 0.99*PPRE_MINUS_HALF(:))
      ZFOESI(:) = MIN(EXP( XALPI - XBETAI/ZT(:) - XGAMI*LOG(ZT(:))  ), 0.99*PPRE_MINUS_HALF(:))
      !Computation of d.Rsat / dP (partial derivations with respect to P and T
      !and use of T=Theta*(P/P0)**(R/Cp) to transform dT into dP with theta_up
      !constant at the vertical)
      ZDRSATODP(:)=(XBETAW/ZT(:)-XGAMW)*(1-ZFRAC_ICE(:))+(XBETAI/ZT(:)-XGAMI)*ZFRAC_ICE(:)
      ZDRSATODP(:)=((XRD/XCPD)*ZDRSATODP(:)-1.)*PRSAT_UP(:)/ &
                  &(PPRE_MINUS_HALF(:)-(ZFOESW(:)*(1-ZFRAC_ICE(:)) + ZFOESI(:)*ZFRAC_ICE(:)))
      !Use of d.Rsat / dP and pressure at flux level KK to find pressure (ZPRE)
      !where Rsat is equal to PRT_UP
      ZPRE(:)=PPRE_MINUS_HALF(:)+(PRT_UP(:)-PRSAT_UP(:))/ZDRSATODP(:)
      !Fraction of dry part (computed with pressure and used with heights, no
      !impact found when using log function here and for pressure on flux levels
      !computation)
      PPART_DRY(:)=MAX(0., MIN(1., (PPRE_MINUS_HALF(:)-ZPRE(:))/(PPRE_MINUS_HALF(:)-PPRE_PLUS_HALF(:))))
      !Height above flux level KK of the cloudy part
      ZDZ_STOP(:) = (PZZ(:,KK+KKL)-PZZ(:,KK))*PPART_DRY(:)
    ENDWHERE
  ENDWHERE

!               1.5 Gradient and flux values of thetav
  IF(KK/=KKB)THEN
    ZCOEFF_MINUS_HALF(:)=((PTHVM(:,KK)-PTHVM(:,KK-KKL))/PDZZ(:,KK))
    ZTHV_MINUS_HALF(:) = PTHVM(:,KK) - ZCOEFF_MINUS_HALF(:)*0.5*(PZZ(:,KK+KKL)-PZZ(:,KK))
  ELSE
    ZCOEFF_MINUS_HALF(:)=0.
    ZTHV_MINUS_HALF(:) = PTHVM(:,KK)
  ENDIF
  ZCOEFF_PLUS_HALF(:)  = ((PTHVM(:,KK+KKL)-PTHVM(:,KK))/PDZZ(:,KK+KKL))
  ZTHV_PLUS_HALF(:)  = PTHVM(:,KK) + ZCOEFF_PLUS_HALF(:)*0.5*(PZZ(:,KK+KKL)-PZZ(:,KK))

!               2  Dry part computation:
!                  Integral buoyancy and computation of PENTR and PDETR for dry part
!               --------------------------------------------------------------------

  DO JI=1,SIZE(PTHLM,1)
    IF(OTEST(JI)) THEN
      IF(PPART_DRY(JI)>0.) THEN
        !Buoyancy computation in two parts to use change of gradient of theta v of environment
        !Between flux level KK and min(mass level, bottom of cloudy part)
        ZDZ(JI)=MIN(ZDZ_STOP(JI),(PZZ(JI,KK+KKL)-PZZ(JI,KK))*0.5)
        PBUO_INTEG_DRY(JI) = ZG_O_THVREF(JI)*ZDZ(JI)*&
                    (0.5 * (  - ZCOEFF_MINUS_HALF(JI))*ZDZ(JI)  &
                      - ZTHV_MINUS_HALF(JI) + PTHV_UP(JI) )
  
        !Between mass flux KK and bottom of cloudy part (if above mass flux)
        ZDZ(JI)=MAX(0., ZDZ_STOP(JI)-(PZZ(JI,KK+KKL)-PZZ(JI,KK))*0.5)
        PBUO_INTEG_DRY(JI) = PBUO_INTEG_DRY(JI) + ZG_O_THVREF(JI)*ZDZ(JI)*&
                    (0.5 * (  - ZCOEFF_PLUS_HALF(JI))*ZDZ(JI) &
                      - PTHVM(JI,KK) + PTHV_UP(JI) )
  
        !Entr//Detr. computation
        IF (PBUO_INTEG_DRY(JI)>=0.) THEN
          PENTR(JI) = 0.5/(XABUO-XBENTR*XENTR_DRY)*&
                     LOG(1.+ (2.*(XABUO-XBENTR*XENTR_DRY)/PW_UP2(JI,KK))* &
                     PBUO_INTEG_DRY(JI))
          PDETR(JI) = 0.
        ELSE
          PENTR(JI) = 0.
          PDETR(JI) = 0.5/(XABUO)*&
                     LOG(1.+ (2.*(XABUO)/PW_UP2(JI,KK))* &
                     (-PBUO_INTEG_DRY(JI)))
        ENDIF
        PENTR(JI) = XENTR_DRY*PENTR(JI)/(PZZ(JI,KK+KKL)-PZZ(JI,KK))    
        PDETR(JI) = XDETR_DRY*PDETR(JI)/(PZZ(JI,KK+KKL)-PZZ(JI,KK))
        !Minimum value of detrainment
        ZWK(JI)=PLUP(JI)-0.5*(PZZ(JI,KK)+PZZ(JI,KK+KKL))
        ZWK(JI)=SIGN(MAX(1., ABS(ZWK(JI))), ZWK(JI)) ! ZWK must not be zero
        PDETR(JI) = MAX(PPART_DRY(JI)*XDETR_LUP/ZWK(JI), PDETR(JI))
      ELSE
        !No dry part, consation reached (OTESTLCL)
        PBUO_INTEG_DRY(JI) = 0.
        PENTR(JI)=0.
        PDETR(JI)=0.
      ENDIF
    ELSE
      !No dry part, consation reached (OTESTLCL)
      PBUO_INTEG_DRY(JI) = 0.
      PENTR(JI)=0.
      PDETR(JI)=0.
    ENDIF
  ENDDO

!               3  Wet part computation
!               -----------------------

!               3.1 Integral buoyancy for cloudy part

  ! Compute theta_v of updraft at flux level KK+KKL                   
  !MIX variables are used to avoid declaring new variables
  !but we are dealing with updraft and not mixture
  ZRCMIX(:)=PRC_UP(:)
  ZRIMIX(:)=PRI_UP(:)
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,ZFRAC_ICE,&
               PPRE_PLUS_HALF,PTHL_UP,PRT_UP,&
               ZTHMIX,ZRVMIX,ZRCMIX,ZRIMIX,&
               ZRSATW, ZRSATI)
  ZTHV_UP_F2(:) = ZTHMIX(:)*(1.+ZRVORD*ZRVMIX(:))/(1.+PRT_UP(:))

  ! Integral buoyancy for cloudy part
  WHERE(OTEST)
    WHERE(PPART_DRY(:)<1.)
      !Gradient of Theta V updraft over the cloudy part, assuming that thetaV updraft don't change
      !between flux level KK and bottom of cloudy part
      ZCOTHVU(:)=(ZTHV_UP_F2(:)-PTHV_UP(:))/((PZZ(:,KK+KKL)-PZZ(:,KK))*(1-PPART_DRY(:)))

      !Computation in two parts to use change of gradient of theta v of environment
      !Between bottom of cloudy part (if under mass level) and mass level KK
      ZDZ(:)=MAX(0., 0.5*(PZZ(:,KK+KKL)-PZZ(:,KK))-ZDZ_STOP(:))
      PBUO_INTEG_CLD(:) = ZG_O_THVREF(:)*ZDZ(:)*&
                (0.5*( ZCOTHVU(:) - ZCOEFF_MINUS_HALF(:))*ZDZ(:) &
                  - (PTHVM(:,KK)-ZDZ(:)*ZCOEFF_MINUS_HALF(:)) + PTHV_UP(:) )

      !Between max(mass level, bottom of cloudy part) and flux level KK+KKL
      ZDZ(:)=(PZZ(:,KK+KKL)-PZZ(:,KK))-MAX(ZDZ_STOP(:),0.5*(PZZ(:,KK+KKL)-PZZ(:,KK)))
      PBUO_INTEG_CLD(:) = PBUO_INTEG_CLD(:)+ZG_O_THVREF(:)*ZDZ(:)*&
                          (0.5*( ZCOTHVU(:) - ZCOEFF_PLUS_HALF(:))*ZDZ(:)&
                  - (PTHVM(:,KK)+(0.5*((PZZ(:,KK+KKL)-PZZ(:,KK)))-ZDZ(:))*ZCOEFF_PLUS_HALF(:)) +&
                  PTHV_UP(:) )

    ELSEWHERE
      !No cloudy part
      PBUO_INTEG_CLD(:)=0.
    ENDWHERE
  ELSEWHERE
    !No cloudy part
    PBUO_INTEG_CLD(:)=0.
  ENDWHERE

!               3.2 Critical mixed fraction for KK+KKL flux level (ZKIC_F2) and
!                   for bottom of cloudy part (ZKIC), then a mean for the cloudy part
!                   (put also in ZKIC)
!
!                   computation by estimating unknown  
!                   T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
!                   We determine the zero crossing of the linear curve
!                   evaluating the derivative using ZMIXF=0.1
                
  ZKIC_INIT=0.1  ! starting value for critical mixed fraction for CLoudy Part

  !  Compute thetaV of environment at the bottom of cloudy part
  !    and cons then non cons. var. of mixture at the bottom of cloudy part

  !   JI computed to avoid KKL(KK-KKL) being < KKL*KKB
  JI=KKL*MAX(KKL*(KK-KKL),KKL*KKB)

  WHERE(OTEST .AND. PPART_DRY(:)>0.5)
    ZDZ(:)=ZDZ_STOP(:)-0.5*(PZZ(:,KK+KKL)-PZZ(:,KK))
    ZTHV(:)= PTHVM(:,KK)+ZCOEFF_PLUS_HALF(:)*ZDZ(:)
    ZMIXTHL(:) = ZKIC_INIT * &
                 (PTHLM(:,KK)+ZDZ(:)*(PTHLM(:,KK+KKL)-PTHLM(:,KK))/PDZZ(:,KK+KKL)) + &
                 (1. - ZKIC_INIT)*PTHL_UP(:)
    ZMIXRT(:)  = ZKIC_INIT * &
                 (PRTM(:,KK)+ZDZ(:)*(PRTM(:,KK+KKL)-PRTM(:,KK))/PDZZ(:,KK+KKL)) +   &
                 (1. - ZKIC_INIT)*PRT_UP(:)
  ELSEWHERE(OTEST)
    ZDZ(:)=0.5*(PZZ(:,KK+KKL)-PZZ(:,KK))-ZDZ_STOP(:)
    ZTHV(:)= PTHVM(:,KK)-ZCOEFF_MINUS_HALF(:)*ZDZ(:)
    ZMIXTHL(:) = ZKIC_INIT * &
                 (PTHLM(:,KK)-ZDZ(:)*(PTHLM(:,KK)-PTHLM(:,JI))/PDZZ(:,KK)) + &
                 (1. - ZKIC_INIT)*PTHL_UP(:)
    ZMIXRT(:)  = ZKIC_INIT * &
                 (PRTM(:,KK)-ZDZ(:)*(PRTM(:,KK)-PRTM(:,JI))/PDZZ(:,KK)) + &
                 (1. - ZKIC_INIT)*PRT_UP(:)
  ENDWHERE
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,ZFRAC_ICE,&
               ZPRE,ZMIXTHL,ZMIXRT,&
               ZTHMIX,ZRVMIX,PRC_MIX,PRI_MIX,&
               ZRSATW, ZRSATI)
  ZTHVMIX(:) = ZTHMIX(:)*(1.+ZRVORD*ZRVMIX(:))/(1.+ZMIXRT(:))

  !  Compute cons then non cons. var. of mixture at the flux level KK+KKL  with initial ZKIC
  ZMIXTHL(:) = ZKIC_INIT * 0.5*(PTHLM(:,KK)+PTHLM(:,KK+KKL))+(1. - ZKIC_INIT)*PTHL_UP(:)
  ZMIXRT(:)  = ZKIC_INIT * 0.5*(PRTM(:,KK)+PRTM(:,KK+KKL))+(1. - ZKIC_INIT)*PRT_UP(:)
  CALL TH_R_FROM_THL_RT_1D(HFRAC_ICE,ZFRAC_ICE,&
               PPRE_PLUS_HALF,ZMIXTHL,ZMIXRT,&
               ZTHMIX,ZRVMIX,PRC_MIX,PRI_MIX,&
               ZRSATW, ZRSATI)
  ZTHVMIX_F2(:) = ZTHMIX(:)*(1.+ZRVORD*ZRVMIX(:))/(1.+ZMIXRT(:))

  !Computation of mean ZKIC over the cloudy part
  WHERE (OTEST)
    ! Compute ZKIC at the bottom of cloudy part
    ! Thetav_up at bottom is equal to Thetav_up at flux level KK
    WHERE (ABS(PTHV_UP(:)-ZTHVMIX(:))<1.E-10)
      ZKIC(:)=1.
    ELSEWHERE
      ZKIC(:) = MAX(0.,PTHV_UP(:)-ZTHV(:))*ZKIC_INIT /  &  
                   (PTHV_UP(:)-ZTHVMIX(:))
    ENDWHERE
    ! Compute ZKIC_F2 at flux level KK+KKL
    WHERE (ABS(ZTHV_UP_F2(:)-ZTHVMIX_F2(:))<1.E-10)
      ZKIC_F2(:)=1.
    ELSEWHERE
      ZKIC_F2(:) = MAX(0.,ZTHV_UP_F2(:)-ZTHV_PLUS_HALF(:))*ZKIC_INIT /  &  
                   (ZTHV_UP_F2(:)-ZTHVMIX_F2(:))
    ENDWHERE
    !Mean ZKIC over the cloudy part
    ZKIC(:)=MAX(MIN(0.5*(ZKIC(:)+ZKIC_F2(:)),1.),0.)
  ENDWHERE

!               3.3 Integration of PDF
!                   According to Kain and Fritsch (1990), we replace delta Mt
!                   in eq. (7) and (8) using eq. (5). Here we compute the ratio
!                   of integrals without computing delta Me

  !Constant PDF
  !For this PDF, eq. (5) is delta Me=0.5*delta Mt
  WHERE(OTEST)
    ZEPSI(:) = ZKIC(:)**2. !integration multiplied by 2
    ZDELTA(:) = (1.-ZKIC(:))**2. !idem
  ENDWHERE

  !Triangular PDF
  !Calculus must be verified before activating this part, but in this state,
  !results on ARM case are almost identical
  !For this PDF, eq. (5) is also delta Me=0.5*delta Mt
  !WHERE(OTEST)
  !  !Integration multiplied by 2
  !  WHERE(ZKIC<0.5)
  !    ZEPSI(:)=8.*ZKIC(:)**3/3.
  !    ZDELTA(:)=1.-4.*ZKIC(:)**2+8.*ZKIC(:)**3/3.
  !  ELSEWHERE
  !    ZEPSI(:)=5./3.-4*ZKIC(:)**2+8.*ZKIC(:)**3/3.
  !    ZDELTA(:)=8.*(1.-ZKIC(:))**3/3.
  !  ENDWHERE
  !ENDWHERE

!               3.4 Computation of PENTR and PDETR
  WHERE (OTEST)
    ZEPSI_CLOUD=MIN(ZDELTA,ZEPSI)
    PENTR_CLD(:) = (1.-PPART_DRY(:))*ZCOEFFMF_CLOUD*PRHODREF(:)*ZEPSI_CLOUD(:)
    PDETR_CLD(:) = (1.-PPART_DRY(:))*ZCOEFFMF_CLOUD*PRHODREF(:)*ZDELTA(:)
    PENTR(:) = PENTR(:)+PENTR_CLD(:)
    PDETR(:) = PDETR(:)+PDETR_CLD(:)
  ELSEWHERE
    PENTR_CLD(:) = 0.
    PDETR_CLD(:) = 0.
  ENDWHERE

IF (LHOOK) CALL DR_HOOK('COMPUTE_ENTR_DETR',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_ENTR_DETR
