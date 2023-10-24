!MNH_LIC Copyright 2022-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!     ##########################
      MODULE MODI_COMPUTE_LAMBDA
!     ##########################
!
INTERFACE
      SUBROUTINE COMPUTE_LAMBDA (KID, KMOMENT, KSIZE, &
                                 PRHO, PRTMIN, PRX, PCX, PLBDX)
!
INTEGER,                INTENT(IN)    :: KID      ! nb of the hydrometeor category
INTEGER,                INTENT(IN)    :: KMOMENT  ! nb of moments of the microphysics scheme
INTEGER,                INTENT(IN)    :: KSIZE
REAL,                   INTENT(IN)    :: PRTMIN
REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRHO     ! reference density
REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRX      ! Mixing ratio
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PCX      ! Nb concentration
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PLBDX    ! Slope parameter of the distribution
!
END SUBROUTINE COMPUTE_LAMBDA
END INTERFACE
END MODULE MODI_COMPUTE_LAMBDA
!
!
! #########################################################
  SUBROUTINE COMPUTE_LAMBDA (KID, KMOMENT, KSIZE, &
                             PRHO, PRTMIN, PRX, PCX, PLBDX)
! #########################################################
!
! Purpose : compute lambda, the slope parameter of the distribution
!           - for 1-moment species: lbda_x = [(rho r_x) / (a_x C_x G(b/alpha))]^(1/(x-b))
!           - for 2-moment species: lbda_x = [(rho r_x) / (a_x N_x G(b/alpha))]^(-1/b)
!
!   AUTHOR
!   ------
!     C. Barthe    * LAERO *
!
!   MODIFICATIONS
!   -------------
!     Original    June 2022
!     C. Barthe    12/07/23  adapt the code for LIMA2
!
!------------------------------------------------------------------
!
!*      0.      DECLARATIONS
!               ------------
!
USE MODD_PARAM_n, ONLY: CCLOUD
USE MODD_RAIN_ICE_DESCR_n, ONLY: XLBC_I=>XLBC, XLBR_I=>XLBR, XLBI_I=>XLBI, XLBS_I=>XLBS, XLBG_I=>XLBG, XLBH_I=>XLBH, &
                                 XLBEXC_I=>XLBEXC, XLBEXR_I=>XLBEXR, XLBEXI_I=>XLBEXI, XLBEXS_I=>XLBEXS,             &
                                 XLBEXG_I=>XLBEXG, XLBEXH_I=>XLBEXH,                                                 &
                                 XLBDAS_MAX_I=>XLBDAS_MAX,                                                           &
                                 XCCR_I=>XCCR, XCCS_I=>XCCS, XCCG_I=>XCCG, XCCH_I=>XCCH,                             &
                                 XCXS_I=>XCXS, XCXG_I=>XCXG, XCXH_I=>XCXH
USE MODD_ELEC_DESCR,       ONLY: XCXR_I=>XCXR
USE MODD_PARAM_LIMA_WARM,  ONLY: XLBC_L=>XLBC, XLBR_L=>XLBR, XLBEXC_L=>XLBEXC, XLBEXR_L=>XLBEXR, &
                                 XCXR_L=>XCXR, XCCR_L=>XCCR
USE MODD_PARAM_LIMA_COLD,  ONLY: XLBI_L=>XLBI, XLBS_L=>XLBS, XLBEXI_L=>XLBEXI, XLBEXS_L=>XLBEXS, &
                                 XLBDAS_MAX_L=>XLBDAS_MAX,                                       &
                                 XCXS_L=>XCXS, XCCS_L=>XCCS
USE MODD_PARAM_LIMA_MIXED, ONLY: XLBG_L=>XLBG, XLBH_L=>XLBH, XLBEXG_L=>XLBEXG, XLBEXH_L=>XLBEXH, &
                                 XCXG_L=>XCXG, XCCG_L=>XCCG, XCXH_L=>XCXH, XCCH_L=>XCCH
!
USE MODI_MOMG
!
IMPLICIT NONE
!
!*      0.1   Declaration of dummy arguments
!
INTEGER,                INTENT(IN)    :: KID      ! nb of the hydrometeor category
INTEGER,                INTENT(IN)    :: KMOMENT  ! nb of moments of the microphysics scheme
INTEGER,                INTENT(IN)    :: KSIZE
REAL,                   INTENT(IN)    :: PRTMIN
REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRHO     ! reference density
REAL, DIMENSION(KSIZE), INTENT(IN)    :: PRX      ! Mixing ratio
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PCX      ! Nb concentration
REAL, DIMENSION(KSIZE), INTENT(INOUT) :: PLBDX    ! Slope parameter of the distribution
!
!*     0.2    Declaration of local variables
!
REAL :: ZRTMIN, ZLBX, ZLBEX, ZLBDAX_MAX, ZCCX, ZCXX
!
!---------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
ZRTMIN = PRTMIN
!
IF (KID == 2) THEN
  IF (CCLOUD == 'LIMA' .AND. KMOMENT == 2) THEN
    ZLBX  = XLBC_L
    ZLBEX = XLBEXC_L
!  ELSE
!    print*, 'ERROR: the computation of lambda_c is not available if c is 1-moment species'
  END IF
ELSE IF (KID == 3) THEN
  IF (CCLOUD == 'LIMA') THEN
    ZLBX  = XLBR_L
    ZLBEX = XLBEXR_L
    IF (KMOMENT == 1) THEN
      ZCCX = XCCR_L
      ZCXX = XCXR_L
    END IF  
  ELSE IF (CCLOUD(1:3) == 'ICE') THEN
    ZLBX  = XLBR_I
    ZLBEX = XLBEXR_I
    ZCCX  = XCCR_I
    ZCXX  = XCXR_I
  ELSE
    PRINT*, 'ERROR: something wrong with the computation of lambda_r'
  END IF
ELSE IF (KID == 4) THEN
  IF (CCLOUD == 'LIMA') THEN
    ZLBX  = XLBI_L
    ZLBEX = XLBEXI_L
  ELSE IF (CCLOUD(1:3) == 'ICE') THEN
    ZLBX  = XLBI_I
    ZLBEX = XLBEXI_I
  ELSE
    PRINT*, 'ERROR: something wrong with the computation of lambda_i'
  END IF
ELSE IF (KID == 5) THEN
  IF (CCLOUD == 'LIMA') THEN
    ZLBX  = XLBS_L
    ZLBEX = XLBEXS_L
    ZLBDAX_MAX = XLBDAS_MAX_L
    IF (KMOMENT == 1) THEN
      ZCCX = XCCS_L
      ZCXX = XCXS_L
    END IF
  ELSE IF (CCLOUD(1:3) == 'ICE') THEN
    ZLBX  = XLBS_I
    ZLBEX = XLBEXS_I
    ZLBDAX_MAX = XLBDAS_MAX_I
    ZCCX  = XCCS_I
    ZCXX  = XCXS_I
  ELSE
    PRINT*, 'ERROR: something wrong with the computation of lambda_s'
  END IF
ELSE IF (KID == 6) THEN
  IF (CCLOUD == 'LIMA') THEN ! .AND. KMOMENT == 1) THEN
    ZLBX  = XLBG_L
    ZLBEX = XLBEXG_L
    IF (KMOMENT == 1) THEN
      ZCCX = XCCG_L
      ZCXX = XCXG_L
    END IF
  ELSE IF (CCLOUD(1:3) == 'ICE') THEN
    ZLBX  = XLBG_I
    ZLBEX = XLBEXG_I
    ZCCX  = XCCG_I
    ZCXX  = XCXG_I
  ELSE
    PRINT*, 'ERROR: something wrong with the computation of lambda_g'
  END IF        
ELSE IF (KID == 7) THEN
  IF (CCLOUD == 'LIMA') THEN ! .AND. KMOMENT == 1) THEN
    ZLBX  = XLBH_L
    ZLBEX = XLBEXH_L
    IF (KMOMENT == 1) THEN
      ZCCX = XCCH_L
      ZCXX = XCXH_L
    END IF
  ELSE IF (CCLOUD(1:3) == 'ICE') THEN
    ZLBX  = XLBH_I
    ZLBEX = XLBEXH_I
    ZCCX  = XCCH_I
    ZCXX  = XCXH_I
  ELSE
    PRINT*, 'ERROR: something wrong with the computation of lambda_h'
  END IF
END IF
!
PLBDX(:) = 0. !1.E10
!
!
!*       2.     COMPUTE LBDA_x FOR 2-MOMENT SPECIES
!               -----------------------------------
!
IF (KMOMENT == 2) THEN
  WHERE (PRX(:) > ZRTMIN .AND. PCX(:) > 0.)
    PLBDX(:) = (ZLBX * PCX(:) / PRX(:))**ZLBEX
  END WHERE
  IF (KID == 5) PLBDX(:) = MIN(ZLBDAX_MAX, PLBDX(:))
!
!
!*       3.     COMPUTE LBDA_x and N_x FOR 1-MOMENT SPECIES
!               -------------------------------------------
!
ELSE IF (KMOMENT == 1) THEN
!
!*       3.1    Special case of cloud droplets
!
  IF (KID == 2) THEN
!    print*, 'computation of lambda_c in 1-moment configuration not treated'
!
!*       3.2    Special case of ice crystals
!
  ELSE IF (KID == 4) THEN
! formulation utilisee dans rain_ice_fast_ri
    WHERE (PRX(:) > ZRTMIN .AND. PCX(:) > 0.0)
      PLBDX(:) = ZLBX * (PRHO(:) * PRX(:) / PCX(:))**ZLBEX
    ENDWHERE
!
!*       3.3    Special case of snow
!
  ELSE IF (KID == 5) THEN
! limitation of lbdas
    WHERE (PRX(:) > ZRTMIN)
      PLBDX(:) = MIN(200000., ZLBX * (PRHO(:) * PRX(:))**ZLBEX)
      PCX(:)   = ZCCX * PLBDX(:)**ZCXX / PRHO(:)
    ENDWHERE
!
!*       3.4    Computation for all other hydrometeors
!
  ELSE
    WHERE (PRX(:) > ZRTMIN)
      PLBDX(:) = ZLBX * (PRHO(:) * PRX(:))**ZLBEX
      PCX(:)   = ZCCX * PLBDX(:)**ZCXX / PRHO(:)
    ENDWHERE
  END IF
END IF
!
END  SUBROUTINE COMPUTE_LAMBDA
