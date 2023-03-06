!MNH_LIC Copyright 2002-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_VER_SV_CORR
IMPLICIT NONE
CONTAINS
SUBROUTINE TURB_VER_SV_CORR(D,CST,CSTURB,TLES,KRR,KRRL,KRRI,OOCEAN, &
                      PDZZ,KSV,KSV_LGBEG,KSV_LGEND,ONOMIXLG,        &
                      OBLOWSNOW,OCOMPUTE_SRC,PRSNOW,                &
                      PTHLM,PRM,PTHVREF,                            &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PPHI3,PPSI3,  &
                      PWM,PSVM,                                     &
                      PTKEM,PLM,PLEPS,PPSI_SV                       )
!     ###############################################################
!
!
!!****  *TURB_VER_SV_CORR* -compute the subgrid Sv2 and SvThv terms
!!
!!    PURPOSE
!!    -------
!!
!!            
!!    EXTERNAL
!!    --------
!!
!!      FUNCTIONs ETHETA and EMOIST  :  
!!            allows to compute:
!!            - the coefficients for the turbulent correlation between
!!            any variable and the virtual potential temperature, of its 
!!            correlations with the conservative potential temperature and 
!!            the humidity conservative variable:
!!            -------              -------              -------
!!            A' Thv'  =  ETHETA   A' Thl'  +  EMOIST   A' Rnp'  
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson               * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       October  29, 2002
!!      JP Pinty       Feb      20, 2003 Add PFRAC_ICE
!!--------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPVEXT_TURB
USE MODD_LES, ONLY: TLES_t
!
USE MODE_SHUMAN_PHY, ONLY:  MZF_PHY
USE MODE_GRADIENT_M_PHY, ONLY : GZ_M_W_PHY
USE MODE_EMOIST, ONLY: EMOIST
USE MODE_ETHETA, ONLY: ETHETA
USE MODI_LES_MEAN_SUBGRID_PHY
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   ::  D
TYPE(CST_t),            INTENT(IN)   ::  CST
TYPE(CSTURB_t),         INTENT(IN)   ::  CSTURB
TYPE(TLES_t),           INTENT(INOUT)::  TLES         ! modd_les structure
INTEGER,                INTENT(IN)   ::  KSV, KSV_LGBEG, KSV_LGEND ! number of scalar variables
LOGICAL,                INTENT(IN)   ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                INTENT(IN)   ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                INTENT(IN)   ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
LOGICAL,                INTENT(IN)   ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and
REAL,                   INTENT(IN)   ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
INTEGER,                INTENT(IN)   ::  KRRL         ! number of liquid var.
INTEGER,                INTENT(IN)   ::  KRRI         ! number of ice var.
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PDZZ
                                                      ! Metric coefficients
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHLM        ! potential temperature at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN) ::  PRM          ! Mixing ratios at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTHVREF      ! reference Thv
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(MERGE(D%NIJT,0,OCOMPUTE_SRC),&
                MERGE(D%NKT,0,OCOMPUTE_SRC)), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PPHI3        ! Inv.Turb.Sch.for temperature
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PPSI3        ! Inv.Turb.Sch.for humidity
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PWM          ! w at time t
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(D%NIJT,D%NKT,KSV), INTENT(IN) ::  PPSI_SV      ! Inv.Turb.Sch.for scalars
                            ! cumulated sources for the prognostic variables
!
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(D%NIJT,D%NKT)  :: ZA, ZFLXZ, &
              ZWORK1,ZWORK2,ZWORK3! working var. for shuman operators (array syntax)
!
REAL :: ZCSV          !constant for the scalar flux
!
INTEGER             :: JIJ,JK,JSV          ! loop counters
INTEGER             :: IIJB, IIJE, IKT
!
REAL :: ZTIME1, ZTIME2
!
REAL :: ZCSVD  = 1.2  ! constant for scalar variance dissipation
REAL :: ZCTSVD = 2.4  ! constant for temperature - scalar covariance dissipation
REAL :: ZCQSVD = 2.4  ! constant for humidity - scalar covariance dissipation
!----------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TURB_VER_SV_CORR',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
CALL SECOND_MNH(ZTIME1)
!
IF(OBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables          
   ZCSV= CSTURB%XCHF/PRSNOW
ELSE
   ZCSV= CSTURB%XCHF
ENDIF
!
DO JSV=1,KSV
  !
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  !
  ! variance Sv2
  !
  IF (TLES%LLES_CALL) THEN
    ! approximation: diagnosed explicitely (without implicit term)
    CALL GZ_M_W_PHY(D,PSVM(:,:,JSV),PDZZ,ZWORK1)
    CALL MZF_PHY(D,ZFLXZ,ZWORK2)
    CALL MZF_PHY(D,PWM,ZWORK3)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZFLXZ(IIJB:IIJE,1:IKT) =  PPSI_SV(IIJB:IIJE,1:IKT,JSV)*ZWORK1(IIJB:IIJE,1:IKT)**2
    ZFLXZ(IIJB:IIJE,1:IKT) = ZCSV / ZCSVD * PLM(IIJB:IIJE,1:IKT) * PLEPS(IIJB:IIJE,1:IKT) * ZWORK2(IIJB:IIJE,1:IKT)
    ZWORK1(IIJB:IIJE,1:IKT) = -2.*ZCSVD*SQRT(PTKEM(IIJB:IIJE,1:IKT))*ZFLXZ(IIJB:IIJE,1:IKT)/PLEPS(IIJB:IIJE,1:IKT)
    ZWORK2(IIJB:IIJE,1:IKT) = ZWORK3(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK1, TLES%X_LES_SUBGRID_DISS_Sv2(:,:,:,JSV) )
    CALL LES_MEAN_SUBGRID_PHY(D,TLES,ZWORK2, TLES%X_LES_RES_W_SBG_Sv2(:,:,:,JSV) )
  END IF
  !
  ! covariance ThvSv
  !
  IF (TLES%LLES_CALL) THEN
    ! approximation: diagnosed explicitely (without implicit term)
    CALL ETHETA(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM,OOCEAN,OCOMPUTE_SRC,ZA)
    !
    CALL GZ_M_W_PHY(D,PTHLM,PDZZ,ZWORK1)
    CALL GZ_M_W_PHY(D,PSVM(:,:,JSV),PDZZ,ZWORK2)
    !
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZFLXZ(IIJB:IIJE,1:IKT)= ( CSTURB%XCSHF * PPHI3(IIJB:IIJE,1:IKT) + ZCSV * PPSI_SV(IIJB:IIJE,1:IKT,JSV) ) &
                  *  ZWORK1(IIJB:IIJE,1:IKT) *  ZWORK2(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !
    CALL MZF_PHY(D,ZFLXZ,ZWORK3)
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZFLXZ(IIJB:IIJE,1:IKT)= PLM(IIJB:IIJE,1:IKT) * PLEPS(IIJB:IIJE,1:IKT) / (2.*ZCTSVD) * ZWORK3(IIJB:IIJE,1:IKT)
    ZWORK1(IIJB:IIJE,1:IKT) = ZA(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
    ZWORK2(IIJB:IIJE,1:IKT) = -CST%XG/PTHVREF(IIJB:IIJE,1:IKT)/3.*ZA(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)    
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    !
    CALL LES_MEAN_SUBGRID_PHY(D,TLES, ZWORK1, TLES%X_LES_SUBGRID_SvThv(:,:,:,JSV) )
    CALL LES_MEAN_SUBGRID_PHY(D,TLES, ZWORK2, TLES%X_LES_SUBGRID_SvPz(:,:,:,JSV), .TRUE.)
    !
    IF (KRR>=1) THEN
      CALL EMOIST(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM,OOCEAN,ZA)
      !
      CALL GZ_M_W_PHY(D,PRM(:,:,1),PDZZ,ZWORK1)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZFLXZ(IIJB:IIJE,1:IKT)= ( ZCSV * PPSI3(IIJB:IIJE,1:IKT) + ZCSV * PPSI_SV(IIJB:IIJE,1:IKT,JSV) )  &
                    *  ZWORK1(IIJB:IIJE,1:IKT) *  ZWORK2(IIJB:IIJE,1:IKT)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      CALL MZF_PHY(D,ZFLXZ,ZWORK3)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZFLXZ(IIJB:IIJE,1:IKT)= PLM(IIJB:IIJE,1:IKT) * PLEPS(IIJB:IIJE,1:IKT) / (2.*ZCQSVD) * ZWORK3(IIJB:IIJE,1:IKT)
      ZWORK1(IIJB:IIJE,1:IKT) = ZA(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
      ZWORK2(IIJB:IIJE,1:IKT) = -CST%XG/PTHVREF(IIJB:IIJE,1:IKT)/3.*ZA(IIJB:IIJE,1:IKT)*ZFLXZ(IIJB:IIJE,1:IKT)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      CALL LES_MEAN_SUBGRID_PHY(D,TLES, ZWORK1, TLES%X_LES_SUBGRID_SvThv(:,:,:,JSV) , .TRUE.)
      CALL LES_MEAN_SUBGRID_PHY(D,TLES, ZWORK2, TLES%X_LES_SUBGRID_SvPz(:,:,:,JSV), .TRUE.)
    END IF
  END IF
  !
END DO   ! end of scalar loop 
!
CALL SECOND_MNH(ZTIME2)
IF(TLES%LLES_CALL) TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('TURB_VER_SV_CORR',1,ZHOOK_HANDLE)
END SUBROUTINE TURB_VER_SV_CORR
END MODULE MODE_TURB_VER_SV_CORR
