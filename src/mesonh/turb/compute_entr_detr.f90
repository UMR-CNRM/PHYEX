!MNH_LIC Copyright 2009-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODI_COMPUTE_ENTR_DETR
!    ##############################
!
INTERFACE
!
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

!INTEGER,                INTENT(IN)   :: KK
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
LOGICAL,DIMENSION(:),   INTENT(IN)   :: OTEST ! test to see if updraft is running
LOGICAL,DIMENSION(:),   INTENT(IN)   :: OTESTLCL !test of condensation 
CHARACTER(len=1),       INTENT(IN)   :: HFRAC_ICE ! frac_ice can be compute using
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
END SUBROUTINE COMPUTE_ENTR_DETR

END INTERFACE
!
END MODULE MODI_COMPUTE_ENTR_DETR
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
!!      S. Riette & J. Escobar (11/2013) : remove div by 0 on real*4 case
!!      P.Marguinaud Jun 2012: fix uninitialized variable
!!      P.Marguinaud Nov 2012: fix gfortran bug
!!      S. Riette Apr 2013: bugs correction, rewriting (for optimisation) and
!!                          improvement of continuity at the condensation level
!!      S. Riette Nov 2013: protection against zero divide for min value of dry PDETR
!!      R.Honnert Oct 2016 : Update with AROME
!  P. Wautelet 08/02/2019: bugfix: compute ZEPSI_CLOUD only once and only when it is needed
!  P. Wautelet 10/02/2021: bugfix: initialized PPART_DRY everywhere
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!                         
USE MODD_CST
USE MODD_PARAM_MFSHALL_n
USE MODD_PARAMETERS, ONLY: XUNDEF

USE MODE_THERMO

USE MODI_TH_R_FROM_THL_RT_1D

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
CHARACTER(len=1),       INTENT(IN)   :: HFRAC_ICE ! frac_ice can be compute using
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
INTEGER :: JI,JLOOP

!----------------------------------------------------------------------------------
                        
!                1.3 Initialisation
!                ------------------

  
  ZRVORD   = XRV / XRD   !=1.607
  ZG_O_THVREF(:)=XG/PTHVM(:,KK)
  ZCOEFFMF_CLOUD=XENTR_MF * XG / XCRAD_MF
  
  ZFRAC_ICE(:)=PFRAC_ICE(:) ! to not modify fraction of ice
 
  ZPRE(:)=PPRE_MINUS_HALF(:)
  ZMIXTHL(:)=0.1
  ZMIXRT(:)=0.1

  !Initialize PPART_DRY everywhere to prevent access to non-initialized values
  ! (intent(out) arrays have undefined values at subroutine entry)
  PPART_DRY(:) = XUNDEF

!                1.4 Estimation of PPART_DRY
  DO JLOOP=1,SIZE(OTEST)
    IF(OTEST(JLOOP) .AND. OTESTLCL(JLOOP)) THEN
      !No dry part when condensation level is reached
      PPART_DRY(JLOOP)=0.
      ZDZ_STOP(JLOOP)=0.
      ZPRE(JLOOP)=PPRE_MINUS_HALF(JLOOP)
    ELSE IF (OTEST(JLOOP) .AND. .NOT. OTESTLCL(JLOOP)) THEN
      !Temperature at flux level KK
      ZT(JLOOP)=PTH_UP(JLOOP)*(PPRE_MINUS_HALF(JLOOP)/XP00) ** (XRD/XCPD)
      !Saturating vapor pressure at flux level KK
      ZFOESW(JLOOP) = MIN(EXP( XALPW - XBETAW/ZT(JLOOP) - XGAMW*LOG(ZT(JLOOP))  ), 0.99*PPRE_MINUS_HALF(JLOOP))
      ZFOESI(JLOOP) = MIN(EXP( XALPI - XBETAI/ZT(JLOOP) - XGAMI*LOG(ZT(JLOOP))  ), 0.99*PPRE_MINUS_HALF(JLOOP))
      !Computation of d.Rsat / dP (partial derivations with respect to P and T
      !and use of T=Theta*(P/P0)**(R/Cp) to transform dT into dP with theta_up
      !constant at the vertical)
      ZDRSATODP(JLOOP)=(XBETAW/ZT(JLOOP)-XGAMW)*(1-ZFRAC_ICE(JLOOP))+(XBETAI/ZT(JLOOP)-XGAMI)*ZFRAC_ICE(JLOOP)
      ZDRSATODP(JLOOP)=((XRD/XCPD)*ZDRSATODP(JLOOP)-1.)*PRSAT_UP(JLOOP)/ &
                  &(PPRE_MINUS_HALF(JLOOP)-(ZFOESW(JLOOP)*(1-ZFRAC_ICE(JLOOP)) + ZFOESI(JLOOP)*ZFRAC_ICE(JLOOP)))
      !Use of d.Rsat / dP and pressure at flux level KK to find pressure (ZPRE)
      !where Rsat is equal to PRT_UP
      ZPRE(JLOOP)=PPRE_MINUS_HALF(JLOOP)+(PRT_UP(JLOOP)-PRSAT_UP(JLOOP))/ZDRSATODP(JLOOP)
      !Fraction of dry part (computed with pressure and used with heights, no
      !impact found when using log function here and for pressure on flux levels
      !computation)
      PPART_DRY(JLOOP)=MAX(0., MIN(1., (PPRE_MINUS_HALF(JLOOP)-ZPRE(JLOOP))/(PPRE_MINUS_HALF(JLOOP)-PPRE_PLUS_HALF(JLOOP))))
      !Height above flux level KK of the cloudy part
      ZDZ_STOP(JLOOP) = (PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))*PPART_DRY(JLOOP)
    END IF
  END DO

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

DO JLOOP=1,SIZE(OTEST) 
 IF (OTEST(JLOOP) .AND. PPART_DRY(JLOOP)>0.) THEN
    ZDZ(JLOOP)=MIN(ZDZ_STOP(JLOOP),(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))*0.5)
    PBUO_INTEG_DRY(JLOOP) = ZG_O_THVREF(JLOOP)*ZDZ(JLOOP)*&
                (0.5 * (  - ZCOEFF_MINUS_HALF(JLOOP))*ZDZ(JLOOP)  &
                  - ZTHV_MINUS_HALF(JLOOP) + PTHV_UP(JLOOP) )

    !Between mass flux KK and bottom of cloudy part (if above mass flux)
    ZDZ(JLOOP)=MAX(0., ZDZ_STOP(JLOOP)-(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))*0.5)
    PBUO_INTEG_DRY(JLOOP) = PBUO_INTEG_DRY(JLOOP) + ZG_O_THVREF(JLOOP)*ZDZ(JLOOP)*&
                (0.5 * (  - ZCOEFF_PLUS_HALF(JLOOP))*ZDZ(JLOOP) &
                  - PTHVM(JLOOP,KK) + PTHV_UP(JLOOP) )
    IF (PBUO_INTEG_DRY(JLOOP)>=0.) THEN
      PENTR(JLOOP) = 0.5/(XABUO-XBENTR*XENTR_DRY)*&
                 LOG(1.+ (2.*(XABUO-XBENTR*XENTR_DRY)/PW_UP2(JLOOP,KK))* &
                 PBUO_INTEG_DRY(JLOOP))
      PDETR(JLOOP) = 0.  
    ELSE
      PENTR(JLOOP) = 0.
      PDETR(JLOOP) = 0.5/(XABUO)*&
                 LOG(1.+ (2.*(XABUO)/PW_UP2(JLOOP,KK))* &
                 (-PBUO_INTEG_DRY(JLOOP)))
    ENDIF
    PENTR(JLOOP) = XENTR_DRY*PENTR(JLOOP)/(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))    
    PDETR(JLOOP) = XDETR_DRY*PDETR(JLOOP)/(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))
    !Minimum value of detrainment
    ZWK(JLOOP)=PLUP(JLOOP)-0.5*(PZZ(JLOOP,KK)+PZZ(JLOOP,KK+KKL))
    ZWK(JLOOP)=SIGN(MAX(1., ABS(ZWK(JLOOP))), ZWK(JLOOP)) ! ZWK must not be zero
    PDETR(JLOOP) = MAX(PPART_DRY(JLOOP)*XDETR_LUP/ZWK(JLOOP), PDETR(JLOOP))  
 ELSE
    !No dry part, consation reached (OTESTLCL)
    PBUO_INTEG_DRY(JLOOP) = 0.
    PENTR(JLOOP)=0.
    PDETR(JLOOP)=0.
 END IF
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
  DO JLOOP=1,SIZE(OTEST)
    IF(OTEST(JLOOP) .AND. PPART_DRY(JLOOP)<1.) THEN
      !Gradient of Theta V updraft over the cloudy part, assuming that thetaV updraft don't change
      !between flux level KK and bottom of cloudy part
      ZCOTHVU(JLOOP)=(ZTHV_UP_F2(JLOOP)-PTHV_UP(JLOOP))/((PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))*(1-PPART_DRY(JLOOP)))

      !Computation in two parts to use change of gradient of theta v of environment
      !Between bottom of cloudy part (if under mass level) and mass level KK
      ZDZ(JLOOP)=MAX(0., 0.5*(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))-ZDZ_STOP(JLOOP))
      PBUO_INTEG_CLD(JLOOP) = ZG_O_THVREF(JLOOP)*ZDZ(JLOOP)*&
              (0.5*( ZCOTHVU(JLOOP) - ZCOEFF_MINUS_HALF(JLOOP))*ZDZ(JLOOP) &
                - (PTHVM(JLOOP,KK)-ZDZ(JLOOP)*ZCOEFF_MINUS_HALF(JLOOP)) + PTHV_UP(JLOOP) )

      !Between max(mass level, bottom of cloudy part) and flux level KK+KKL
      ZDZ(JLOOP)=(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))-MAX(ZDZ_STOP(JLOOP),0.5*(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK)))
      PBUO_INTEG_CLD(JLOOP) = PBUO_INTEG_CLD(JLOOP)+ZG_O_THVREF(JLOOP)*ZDZ(JLOOP)*&
                        (0.5*( ZCOTHVU(JLOOP) - ZCOEFF_PLUS_HALF(JLOOP))*ZDZ(JLOOP)&
                - (PTHVM(JLOOP,KK)+(0.5*((PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK)))-ZDZ(JLOOP))*ZCOEFF_PLUS_HALF(JLOOP)) +&
                PTHV_UP(JLOOP) )

    ELSE
      !No cloudy part
      PBUO_INTEG_CLD(JLOOP)=0.
    END IF
  END DO

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
  DO JLOOP=1,SIZE(OTEST)
    IF(OTEST(JLOOP) .AND. PPART_DRY(JLOOP)>0.5) THEN
      ZDZ(JLOOP)=ZDZ_STOP(JLOOP)-0.5*(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))
      ZTHV(JLOOP)= PTHVM(JLOOP,KK)+ZCOEFF_PLUS_HALF(JLOOP)*ZDZ(JLOOP)
      ZMIXTHL(JLOOP) = ZKIC_INIT * &
                 (PTHLM(JLOOP,KK)+ZDZ(JLOOP)*(PTHLM(JLOOP,KK+KKL)-PTHLM(JLOOP,KK))/PDZZ(JLOOP,KK+KKL)) + &
                 (1. - ZKIC_INIT)*PTHL_UP(JLOOP)
      ZMIXRT(JLOOP)  = ZKIC_INIT * &
                 (PRTM(JLOOP,KK)+ZDZ(JLOOP)*(PRTM(JLOOP,KK+KKL)-PRTM(JLOOP,KK))/PDZZ(JLOOP,KK+KKL)) +   &
                 (1. - ZKIC_INIT)*PRT_UP(JLOOP)
     ELSEIF(OTEST(JLOOP)) THEN
      ZDZ(JLOOP)=0.5*(PZZ(JLOOP,KK+KKL)-PZZ(JLOOP,KK))-ZDZ_STOP(JLOOP)
      ZTHV(JLOOP)= PTHVM(JLOOP,KK)-ZCOEFF_MINUS_HALF(JLOOP)*ZDZ(JLOOP)
      ZMIXTHL(JLOOP) = ZKIC_INIT * &
                 (PTHLM(JLOOP,KK)-ZDZ(JLOOP)*(PTHLM(JLOOP,KK)-PTHLM(JLOOP,JI))/PDZZ(JLOOP,KK)) + &
                 (1. - ZKIC_INIT)*PTHL_UP(JLOOP)
      ZMIXRT(JLOOP)  = ZKIC_INIT * &
                 (PRTM(JLOOP,KK)-ZDZ(JLOOP)*(PRTM(JLOOP,KK)-PRTM(JLOOP,JI))/PDZZ(JLOOP,KK)) + &
                 (1. - ZKIC_INIT)*PRT_UP(JLOOP)
    ENDIF
  ENDDO
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
  DO JLOOP=1,SIZE(OTEST) 
    IF (OTEST(JLOOP)) THEN
    ! Compute ZKIC at the bottom of cloudy part
    ! Thetav_up at bottom is equal to Thetav_up at flux level KK
      IF (ABS(PTHV_UP(JLOOP)-ZTHVMIX(JLOOP))<1.E-10) THEN
        ZKIC(JLOOP)=1.
      ELSE
        ZKIC(JLOOP) = MAX(0.,PTHV_UP(JLOOP)-ZTHV(JLOOP))*ZKIC_INIT /  &  
                   (PTHV_UP(JLOOP)-ZTHVMIX(JLOOP))
      END IF
      ! Compute ZKIC_F2 at flux level KK+KKL
      IF (ABS(ZTHV_UP_F2(JLOOP)-ZTHVMIX_F2(JLOOP))<1.E-10) THEN
        ZKIC_F2(JLOOP)=1.
      ELSE
        ZKIC_F2(JLOOP) = MAX(0.,ZTHV_UP_F2(JLOOP)-ZTHV_PLUS_HALF(JLOOP))*ZKIC_INIT /  &  
                   (ZTHV_UP_F2(JLOOP)-ZTHVMIX_F2(JLOOP))
      END IF
      !Mean ZKIC over the cloudy part
      ZKIC(JLOOP)=MAX(MIN(0.5*(ZKIC(JLOOP)+ZKIC_F2(JLOOP)),1.),0.)
    END IF
  END DO


!               3.3 Integration of PDF
!                   According to Kain and Fritsch (1990), we replace delta Mt
!                   in eq. (7) and (8) using eq. (5). Here we compute the ratio
!                   of integrals without computing delta Me

  !Constant PDF
  !For this PDF, eq. (5) is delta Me=0.5*delta Mt
  DO JLOOP=1,SIZE(OTEST)
    IF(OTEST(JLOOP)) THEN
      ZEPSI(JLOOP) = ZKIC(JLOOP)**2. !integration multiplied by 2
      ZDELTA(JLOOP) = (1.-ZKIC(JLOOP))**2. !idem
    ENDIF
  ENDDO

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
  DO JLOOP=1,SIZE(OTEST)
    IF(OTEST(JLOOP)) THEN
      ZEPSI_CLOUD(JLOOP)=MIN(ZDELTA(JLOOP),ZEPSI(JLOOP))
      PENTR_CLD(JLOOP) = (1.-PPART_DRY(JLOOP))*ZCOEFFMF_CLOUD*PRHODREF(JLOOP)*ZEPSI_CLOUD(JLOOP)
      PDETR_CLD(JLOOP) = (1.-PPART_DRY(JLOOP))*ZCOEFFMF_CLOUD*PRHODREF(JLOOP)*ZDELTA(JLOOP)
      PENTR(JLOOP) = PENTR(JLOOP)+PENTR_CLD(JLOOP)
      PDETR(JLOOP) = PDETR(JLOOP)+PDETR_CLD(JLOOP)
    ELSE
      PENTR_CLD(JLOOP) = 0.
      PDETR_CLD(JLOOP) = 0.
    ENDIF
  ENDDO 

END SUBROUTINE COMPUTE_ENTR_DETR  
