!MNH_LIC Copyright 2002-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_TURB_HOR_SV_CORR
IMPLICIT NONE
CONTAINS
      SUBROUTINE TURB_HOR_SV_CORR(D,CST,CSTURB,TURBN,TLES,KSV,KSV_LGBEG,KSV_LGEND,&
                      KRR,KRRL,KRRI,OOCEAN,OCOMPUTE_SRC,OBLOWSNOW,   &
                      ONOMIXLG,O2D,                                  &
                      PDXX,PDYY,PDZZ,PDZX,PDZY,PRSNOW,               &
                      PLM,PLEPS,PTKEM,PTHVREF,                       &
                      PTHLM,PRM,                                     &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,               &
                      PWM,PSVM                                       )
!     ################################################################
!
!
!!****  *TURB_HOT_SV_CORR*  computes subgrid Sv2 and SvThv terms
!!
!!    PURPOSE
!!    -------
!!
!!     see TURB_HOR
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!
!!      V. Masson               * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!   Original  06/11/02
!!      JP Pinty       Feb 20, 2003 Add PFRAC_ICE
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST, ONLY: CST_t
USE MODD_CTURB, ONLY : CSTURB_t
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_PARAMETERS
USE MODD_LES, ONLY: TLES_t
!
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_SHUMAN 
USE MODI_LES_MEAN_SUBGRID
USE MODE_EMOIST, ONLY: EMOIST
USE MODE_ETHETA, ONLY: ETHETA
!
USE MODI_SECOND_MNH
!
IMPLICIT NONE
!
!
!*       0.1  declaration of arguments
!
!
!
TYPE(DIMPHYEX_t),         INTENT(IN)    ::  D
TYPE(CST_t),              INTENT(IN)    ::  CST
TYPE(CSTURB_t),           INTENT(IN)    ::  CSTURB
TYPE(TURB_t),             INTENT(IN)    :: TURBN
TYPE(TLES_t),             INTENT(INOUT) :: TLES          ! modd_les structure
INTEGER,                  INTENT(IN)    ::  KRR          ! number of moist var.
INTEGER,                  INTENT(IN)    ::  KRRL         ! number of liquid var.
INTEGER,                  INTENT(IN)    ::  KRRI         ! number of ice var.
INTEGER,                  INTENT(IN)    ::  KSV,KSV_LGBEG,KSV_LGEND ! number of sv var.
LOGICAL,                  INTENT(IN)    ::  OOCEAN       ! switch for Ocean model version
LOGICAL,                  INTENT(IN)    ::  OCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables
LOGICAL,                  INTENT(IN)    ::  OBLOWSNOW    ! switch to activate pronostic blowing snow
LOGICAL,                  INTENT(IN)    ::  ONOMIXLG     ! to use turbulence for lagrangian variables (modd_conf)
LOGICAL,                  INTENT(IN)    ::  O2D          ! Logical for 2D model version (modd_conf)
REAL,                     INTENT(IN)    ::  PRSNOW       ! Ratio for diffusion coeff. scalar (blowing snow)
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PDXX, PDYY, PDZZ, PDZX, PDZY 
                                                         ! Metric coefficients
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLM          ! mixing length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLEPS        ! dissipative length
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTKEM        ! tke
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHVREF      ! reference Thv
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PTHLM        ! potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PRM          ! Mixing ratios at t-Delta t
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:),   INTENT(IN)    ::  PWM          ! w at t-1
REAL, DIMENSION(:,:,:,:), INTENT(IN)    ::  PSVM         ! scalar var. at t-1
!
!
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(SIZE(PSVM,1),SIZE(PSVM,2),SIZE(PSVM,3))       &
                                     :: ZFLX, ZA
!
INTEGER             :: JSV          ! loop counter
INTEGER             :: IKU
!
REAL :: ZTIME1, ZTIME2
!
REAL :: ZCSVD  = 1.2  ! constant for scalar variance dissipation
REAL :: ZCTSVD = 2.4  ! constant for temperature - scalar covariance dissipation
REAL :: ZCQSVD = 2.4  ! constant for humidity - scalar covariance dissipation
!
REAL :: ZCSV          !constant for the scalar flux 
! ---------------------------------------------------------------------------
!
IKU=SIZE(PTKEM,3)
CALL SECOND_MNH(ZTIME1)
!
IF(OBLOWSNOW) THEN
! See Vionnet (PhD, 2012) for a complete discussion around the value of the Schmidt number for blowing snow variables        
   ZCSV= TURBN%XCHF/PRSNOW 
ELSE
   ZCSV= TURBN%XCHF
ENDIF
!
DO JSV=1,KSV
!
  IF (ONOMIXLG .AND. JSV >= KSV_LGBEG .AND. JSV<= KSV_LGEND) CYCLE
  !
  ! variance Sv2
  !
  IF (TLES%LLES_CALL) THEN
    IF (.NOT. O2D) THEN
      ZFLX(:,:,:) =  ZCSV / ZCSVD * PLM(:,:,:) * PLEPS(:,:,:) *   &
         (  GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)**2             &
          + GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)**2 )
    ELSE
      ZFLX(:,:,:) =  ZCSV / ZCSVD * PLM(:,:,:) * PLEPS(:,:,:) *   &
            GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)**2
    END IF
    CALL LES_MEAN_SUBGRID( -2.*ZCSVD*SQRT(PTKEM)*ZFLX/PLEPS, &
                           TLES%X_LES_SUBGRID_DISS_Sv2(:,:,:,JSV), .TRUE. )
    CALL LES_MEAN_SUBGRID( MZF(PWM)*ZFLX, TLES%X_LES_RES_W_SBG_Sv2(:,:,:,JSV), .TRUE. )
  END IF
  !
  ! covariance SvThv
  !
  IF (TLES%LLES_CALL) THEN
    CALL ETHETA(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM,OOCEAN,OCOMPUTE_SRC,ZA)
    IF (.NOT. O2D) THEN
      ZFLX(:,:,:)=  PLM(:,:,:) * PLEPS(:,:,:)                                          &
          *  (  GX_M_M(PTHLM,PDXX,PDZZ,PDZX) * GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)  &
              + GY_M_M(PTHLM,PDYY,PDZZ,PDZY) * GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)  &
             ) * (TURBN%XCSHF+ZCSV) / (2.*ZCTSVD)
    ELSE
      ZFLX(:,:,:)=  PLM(:,:,:) * PLEPS(:,:,:)                                          &
              * GX_M_M(PTHLM,PDXX,PDZZ,PDZX) * GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)  &
              * (TURBN%XCSHF+ZCSV) / (2.*ZCTSVD)
    END IF
    CALL LES_MEAN_SUBGRID( ZA*ZFLX, TLES%X_LES_SUBGRID_SvThv(:,:,:,JSV) , .TRUE.)
    CALL LES_MEAN_SUBGRID( -CST%XG/PTHVREF/3.*ZA*ZFLX, TLES%X_LES_SUBGRID_SvPz(:,:,:,JSV), .TRUE. )
    !
    IF (KRR>=1) THEN
      CALL  EMOIST(D,CST,KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PAMOIST,PSRCM,OOCEAN,ZA)
      IF (.NOT. O2D) THEN
        ZFLX(:,:,:)=  PLM(:,:,:) * PLEPS(:,:,:)                                                 &
            *  (  GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX) * GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)  &
                + GY_M_M(PRM(:,:,:,1),PDYY,PDZZ,PDZY) * GY_M_M(PSVM(:,:,:,JSV),PDYY,PDZZ,PDZY)  &
               ) * (TURBN%XCHF+ZCSV) / (2.*ZCQSVD)
      ELSE
        ZFLX(:,:,:)=  PLM(:,:,:) * PLEPS(:,:,:)                                                 &
                * GX_M_M(PRM(:,:,:,1),PDXX,PDZZ,PDZX) * GX_M_M(PSVM(:,:,:,JSV),PDXX,PDZZ,PDZX)  &
                * (TURBN%XCHF+ZCSV) / (2.*ZCQSVD)
      END IF
      CALL LES_MEAN_SUBGRID( ZA*ZFLX, TLES%X_LES_SUBGRID_SvThv(:,:,:,JSV) , .TRUE.)
      CALL LES_MEAN_SUBGRID( -CST%XG/PTHVREF/3.*ZA*ZFLX, TLES%X_LES_SUBGRID_SvPz(:,:,:,JSV), .TRUE. )
    END IF
  END IF
!
END DO    ! end loop JSV
!
CALL SECOND_MNH(ZTIME2)
TLES%XTIME_LES = TLES%XTIME_LES + ZTIME2 - ZTIME1
!
END SUBROUTINE TURB_HOR_SV_CORR
END MODULE MODE_TURB_HOR_SV_CORR

