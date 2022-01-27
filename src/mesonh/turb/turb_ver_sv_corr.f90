!MNH_LIC Copyright 2002-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    #################### 
     MODULE MODI_TURB_VER_SV_CORR
!    ####################
!
INTERFACE 
!
      SUBROUTINE TURB_VER_SV_CORR(KKA,KKU,KKL,KRR,KRRL,KRRI,                    &
                      PDZZ,                                         &
                      PTHLM,PRM,PTHVREF,                            &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PPHI3,PPSI3,  &
                      PWM,PSVM,                                     &
                      PTKEM,PLM,PLEPS,PPSI_SV                       )
!
INTEGER,                 INTENT(IN)  ::  KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                 INTENT(IN)  ::  KKL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
INTEGER,                INTENT(IN)   ::  KRRL         ! number of liquid var.
INTEGER,                INTENT(IN)   ::  KRRI         ! number of ice var.
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHLM        ! potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM          ! Mixing ratios at t-Delta t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF      ! reference Thv
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PPHI3        ! Inv.Turb.Sch.for temperature
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PPSI3        ! Inv.Turb.Sch.for humidity
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PWM          ! w at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PPSI_SV      ! Inv.Turb.Sch.for scalars
                            ! cumulated sources for the prognostic variables
!
!
END SUBROUTINE TURB_VER_SV_CORR
!
END INTERFACE
!
END MODULE MODI_TURB_VER_SV_CORR
!
!
!     ###############################################################
      SUBROUTINE TURB_VER_SV_CORR(KKA,KKU,KKL,KRR,KRRL,KRRI,        &
                      PDZZ,                                         &
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
USE MODD_CST
USE MODD_CTURB
USE MODD_PARAMETERS
USE MODD_LES
USE MODD_CONF
USE MODD_NSV, ONLY : NSV,NSV_LGBEG,NSV_LGEND
USE MODD_BLOWSNOW
!
!
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_GRADIENT_M
USE MODI_SHUMAN 
USE MODI_EMOIST
USE MODI_ETHETA
USE MODI_LES_MEAN_SUBGRID
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!*      0.1  declarations of arguments
!
!
!
INTEGER,                 INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                 INTENT(IN)  :: KKL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
INTEGER,                INTENT(IN)   ::  KRRL         ! number of liquid var.
INTEGER,                INTENT(IN)   ::  KRRI         ! number of ice var.
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ
                                                      ! Metric coefficients
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHLM        ! potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM          ! Mixing ratios at t-Delta t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF      ! reference Thv
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PPHI3        ! Inv.Turb.Sch.for temperature
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PPSI3        ! Inv.Turb.Sch.for humidity
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PWM          ! w at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PPSI_SV      ! Inv.Turb.Sch.for scalars
                            ! cumulated sources for the prognostic variables
!
!
!
!
!*       0.2  declaration of local variables
!
!
REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2),SIZE(PSVM,3))  ::  &
       ZA, ZFLXZ
!
REAL :: ZCSV          !constant for the scalar flux
!
INTEGER             :: JSV          ! loop counters
!
REAL :: ZTIME1, ZTIME2
!
REAL :: ZCSVD  = 1.2  ! constant for scalar variance dissipation
REAL :: ZCTSVD = 2.4  ! constant for temperature - scalar covariance dissipation
REAL :: ZCQSVD = 2.4  ! constant for humidity - scalar covariance dissipation
!----------------------------------------------------------------------------
!
CALL SECOND_MNH(ZTIME1)
!
IF(LBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables          
   ZCSV= XCHF/XRSNOW
ELSE
   ZCSV= XCHF
ENDIF
!
DO JSV=1,NSV
  !
  IF (LNOMIXLG .AND. JSV >= NSV_LGBEG .AND. JSV<= NSV_LGEND) CYCLE
  !
  ! variance Sv2
  !
  IF (LLES_CALL) THEN
    ! approximation: diagnosed explicitely (without implicit term)
    ZFLXZ(:,:,:) =  PPSI_SV(:,:,:,JSV)*GZ_M_W(KKA,KKU,KKL,PSVM(:,:,:,JSV),PDZZ)**2
    ZFLXZ(:,:,:) = ZCSV / ZCSVD * PLM * PLEPS * MZF(ZFLXZ(:,:,:) )
    CALL LES_MEAN_SUBGRID( -2.*ZCSVD*SQRT(PTKEM)*ZFLXZ/PLEPS, X_LES_SUBGRID_DISS_Sv2(:,:,:,JSV) )
    CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLXZ, X_LES_RES_W_SBG_Sv2(:,:,:,JSV) )
  END IF
  !
  ! covariance ThvSv
  !
  IF (LLES_CALL) THEN
    ! approximation: diagnosed explicitely (without implicit term)
    ZA(:,:,:)   =  ETHETA(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM)
    ZFLXZ(:,:,:)= ( XCSHF * PPHI3 + ZCSV * PPSI_SV(:,:,:,JSV) )              &
                  *  GZ_M_W(KKA,KKU,KKL,PTHLM,PDZZ)                          &
                  *  GZ_M_W(KKA,KKU,KKL,PSVM(:,:,:,JSV),PDZZ)
    ZFLXZ(:,:,:)= PLM * PLEPS / (2.*ZCTSVD) * MZF(ZFLXZ)
    CALL LES_MEAN_SUBGRID( ZA*ZFLXZ, X_LES_SUBGRID_SvThv(:,:,:,JSV) ) 
    CALL LES_MEAN_SUBGRID( -XG/PTHVREF/3.*ZA*ZFLXZ, X_LES_SUBGRID_SvPz(:,:,:,JSV), .TRUE.)
    !
    IF (KRR>=1) THEN
      ZA(:,:,:)   =  EMOIST(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM)
      ZFLXZ(:,:,:)= ( XCHF * PPSI3 + ZCSV * PPSI_SV(:,:,:,JSV) )             &
                    *  GZ_M_W(KKA,KKU,KKL,PRM(:,:,:,1),PDZZ)                 &
                    *  GZ_M_W(KKA,KKU,KKL,PSVM(:,:,:,JSV),PDZZ)
      ZFLXZ(:,:,:)= PLM * PLEPS / (2.*ZCQSVD) * MZF(ZFLXZ)
      CALL LES_MEAN_SUBGRID( ZA*ZFLXZ, X_LES_SUBGRID_SvThv(:,:,:,JSV) , .TRUE.)
      CALL LES_MEAN_SUBGRID( -XG/PTHVREF/3.*ZA*ZFLXZ, X_LES_SUBGRID_SvPz(:,:,:,JSV), .TRUE.)
    END IF
  END IF
  !
END DO   ! end of scalar loop 
!
CALL SECOND_MNH(ZTIME2)
XTIME_LES = XTIME_LES + ZTIME2 - ZTIME1
!----------------------------------------------------------------------------
!
END SUBROUTINE TURB_VER_SV_CORR
