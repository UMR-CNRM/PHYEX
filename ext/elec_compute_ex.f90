!MNH_LIC Copyright 2022-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
!     ###########################
      MODULE MODI_ELEC_COMPUTE_EX
!     ###########################
!
INTERFACE
      SUBROUTINE ELEC_COMPUTE_EX (KID, KMOMENT, KSIZE, &
                                  PDUM, PRHO, PRTMIN,  &
                                  PRX, PQX, PEX, PLBDX, PCX)
!
INTEGER,                          INTENT(IN)    :: KID      ! nb of the hydrometeor category
INTEGER,                          INTENT(IN)    :: KMOMENT  ! nb of moments of the microphysics scheme
INTEGER,                          INTENT(IN)    :: KSIZE
REAL,                             INTENT(IN)    :: PDUM     ! =1. if mixing ratio
                                                            ! =timestep if source
REAL,                             INTENT(IN)    :: PRTMIN
REAL, DIMENSION(KSIZE),           INTENT(IN)    :: PRHO     ! reference density
REAL, DIMENSION(KSIZE),           INTENT(IN)    :: PQX      ! Electric charge
REAL, DIMENSION(KSIZE),           INTENT(IN)    :: PRX      ! Mixing ratio
REAL, DIMENSION(KSIZE),           INTENT(INOUT) :: PEX      ! e coef of the q-D relation
REAL, DIMENSION(KSIZE), OPTIONAL, INTENT(IN)    :: PLBDX    ! Slope parameter of the distribution
REAL, DIMENSION(KSIZE), OPTIONAL, INTENT(IN)    :: PCX      ! Nb concentration
!
END SUBROUTINE ELEC_COMPUTE_EX
END INTERFACE
END MODULE MODI_ELEC_COMPUTE_EX
!
!
! #######################################################
  SUBROUTINE ELEC_COMPUTE_EX (KID, KMOMENT, KSIZE,     &
                              PDUM, PRHO, PRTMIN,      &
                              PRX, PQX, PEX, PLBDX, PCX )
! #######################################################
!
! Purpose : update the parameter e_x in the relation q_x = e_x d**f_x
!           e_x = q_x/(N_x *  M(f_x))
!
!   AUTHOR
!   ------
!     C. Barthe    * LAERO *
!
!   MODIFICATIONS
!   -------------
!     Original    June 2022
!
!------------------------------------------------------------------
!
!*      0.      DECLARATIONS
!               ------------
!
USE MODD_PARAM_n,          ONLY : CCLOUD
USE MODD_ELEC_PARAM,       ONLY : XECMAX, XERMAX, XEIMAX, XESMAX, XEGMAX, XEHMAX, &
                                  XFQUPDC, XFQUPDR, XFQUPDI, XEXFQUPDI, XFQUPDS, XFQUPDG, XFQUPDH
USE MODD_ELEC_DESCR,       ONLY : XCXR, XFC, XFR, XFI, XFS, XFG, XFH
USE MODD_RAIN_ICE_DESCR_n, ONLY : XCXS_I=>XCXS, XCXG_I=>XCXG, XCXH_I=>XCXH
USE MODD_PARAM_LIMA_COLD,  ONLY : XCXS_L=>XCXS
USE MODD_PARAM_LIMA_MIXED, ONLY : XCXG_L=>XCXG, XCXH_L=>XCXH, XALPHAH, XNUH
USE MODD_PARAM_LIMA,       ONLY : XALPHAC, XALPHAR, XALPHAI, XALPHAS, XALPHAG, &
                                  XNUC, XNUR, XNUI, XNUS, XNUG
!
USE MODI_MOMG
!
IMPLICIT NONE
!
!*      0.1   Declaration of dummy arguments
!
INTEGER,                          INTENT(IN)    :: KID      ! nb of the hydrometeor category
INTEGER,                          INTENT(IN)    :: KMOMENT  ! nb of moments of the microphysics scheme
INTEGER,                          INTENT(IN)    :: KSIZE
REAL,                             INTENT(IN)    :: PDUM     ! =1. if mixing ratio
                                                            ! =timestep if source
REAL,                             INTENT(IN)    :: PRTMIN
REAL, DIMENSION(KSIZE),           INTENT(IN)    :: PRHO     ! reference density
REAL, DIMENSION(KSIZE),           INTENT(IN)    :: PQX      ! Electric charge
REAL, DIMENSION(KSIZE),           INTENT(IN)    :: PRX      ! Mixing ratio
REAL, DIMENSION(KSIZE),           INTENT(INOUT) :: PEX      ! e coef of the q-D relation
REAL, DIMENSION(KSIZE), OPTIONAL, INTENT(IN)    :: PLBDX    ! Slope parameter of the distribution
REAL, DIMENSION(KSIZE), OPTIONAL, INTENT(IN)    :: PCX      ! Nb concentration
!
!*     0.2    Declaration of local variables
!
REAL :: ZRTMIN, ZFX, ZCX, ZEXMAX, ZFQUPDX, ZALPHAX, ZNUX
!
!---------------------------------------------------------------------------------
!
!*       1.     PRELIMINARIES
!               -------------
!
ZRTMIN = PRTMIN / PDUM
PEX(:) = 0.
!
IF (KID == 2) THEN       ! parameters for cloud droplets
  ZFX     = XFC
  ZEXMAX  = XECMAX
  IF (CCLOUD(1:3) == 'ICE') ZFQUPDX = XFQUPDC
  IF (CCLOUD == 'LIMA') THEN
    ZALPHAX = XALPHAC
    ZNUX    = XNUC
  END IF
ELSE IF (KID == 3) THEN  ! parameters for raindrops
  ZFX     = XFR
  ZCX     = XCXR
  ZEXMAX  = XERMAX
  IF (CCLOUD(1:3) == 'ICE') ZFQUPDX = XFQUPDR
  IF (CCLOUD == 'LIMA' .AND. KMOMENT == 2) THEN
    ZALPHAX = XALPHAR
    ZNUX    = XNUR
  END IF
ELSE IF (KID == 4) THEN  ! parameters for ice crystals
  ZFX     = XFI
  ZEXMAX  = XEIMAX
  IF (CCLOUD(1:3) == 'ICE') ZFQUPDX = XFQUPDI
  IF (CCLOUD == 'LIMA' .AND. KMOMENT == 2) THEN
    ZALPHAX = XALPHAI
    ZNUX    = XNUI
  END IF
ELSE IF (KID == 5) THEN  ! parameters for snow/aggregates
  ZFX     = XFS
  ZEXMAX  = XESMAX
  ZFQUPDX = XFQUPDS
  IF (CCLOUD == 'LIMA' .AND. KMOMENT == 2) THEN
    ZALPHAX = XALPHAS
    ZNUX    = XNUS
  ELSE IF (CCLOUD == 'LIMA' .AND. KMOMENT == 1) THEN
    ZCX   = XCXS_L
  ELSE IF (CCLOUD(1:3) == 'ICE') THEN
    ZCX   = XCXS_I
  END IF
ELSE IF (KID == 6) THEN  ! parameters for graupel
  ZFX     = XFG
  ZEXMAX  = XEGMAX
  ZFQUPDX = XFQUPDG
  IF (CCLOUD == 'LIMA' .AND. KMOMENT == 2) THEN
    ZALPHAX = XALPHAG
    ZNUX    = XNUG
  ELSE IF (CCLOUD == 'LIMA' .AND. KMOMENT == 1) THEN
    ZCX   = XCXG_L
  ELSE IF (CCLOUD(1:3) == 'ICE') THEN
    ZCX   = XCXG_I
  END IF
ELSE IF (KID == 7) THEN  ! parameters for hail
  ZFX     = XFH
  ZEXMAX  = XEHMAX
  ZFQUPDX = XFQUPDH
  IF (CCLOUD == 'LIMA' .AND. KMOMENT == 2) THEN
    ZALPHAX = XALPHAH
    ZNUX    = XNUH
  ELSE IF (CCLOUD == 'LIMA' .AND. KMOMENT == 1) THEN
    ZCX   = XCXH_L
  ELSE IF (CCLOUD(1:3) == 'ICE') THEN
    ZCX   = XCXH_I
  END IF  
END IF
!
IF (CCLOUD == 'LIMA') THEN
  IF (KID == 2) THEN
    ZALPHAX = XALPHAC
    ZNUX    = XNUC
  ELSE IF (KID == 3) THEN
    ZALPHAX = XALPHAR
    ZNUX    = XNUR
  ELSE IF (KID == 4) THEN
    ZALPHAX = XALPHAI
    ZNUX    = XNUI       
  ELSE IF (KID == 5) THEN
    ZALPHAX = XALPHAS
    ZNUX    = XNUS       
  ELSE IF (KID == 6) THEN
    ZALPHAX = XALPHAG
    ZNUX    = XNUG       
  ELSE IF (KID == 7) THEN
    ZALPHAX = XALPHAH
    ZNUX    = XNUH          
  END IF
END IF
!
!
!*       2.     UPDATE E_x FOR 2-MOMENT SPECIES
!               -------------------------------
!
IF (KMOMENT == 2) THEN
  WHERE (PRX(:) > ZRTMIN .AND. PCX(:) > 0.0)
    PEX(:) = PDUM * PRHO(:) * PQX(:) * PLBDX(:)**ZFX / (PCX(:) * MOMG(ZALPHAX,ZNUX,ZFX))       
  ENDWHERE
!
!
!*       3.     UPDATE E_x FOR 1-MOMENT SPECIES
!               -------------------------------
!
ELSE IF (KMOMENT == 1) THEN
!
!*       3.1    Special case of cloud droplets
!
  IF (KID == 2) THEN
    WHERE (PRX(:) > ZRTMIN)
      PEX(:) = PDUM * PRHO(:) * PQX(:) / ZFQUPDX
      PEX(:) = SIGN( MIN(ABS(PEX(:)), ZEXMAX), PEX(:))
    ENDWHERE
!
!*       3.2    Special case of ice crystals
!
  ELSE IF (KID == 4) THEN
    WHERE (PRX(:) > ZRTMIN .AND. PCX(:) > 0.0)
      PEX(:) = PDUM * PRHO(:) * PQX(:) /                      &
               ((PCX**(1 - XEXFQUPDI)) * ZFQUPDX * (PRHO(:) * &
               PDUM * PRX(:))**XEXFQUPDI)
      PEX(:) = SIGN( MIN(ABS(PEX(:)), ZEXMAX), PEX(:))
    ENDWHERE
!
!*       3.3    Computation for all other hydrometeors
!
  ELSE
    WHERE (PRX(:) > ZRTMIN .AND. PLBDX(:) > 0.)
      PEX(:) = PDUM * PRHO(:) * PQX(:) / (ZFQUPDX * PLBDX(:)**(ZCX - ZFX))
      PEX(:) = SIGN( MIN(ABS(PEX(:)), ZEXMAX), PEX(:))
    ENDWHERE
  END IF
END IF
!
END  SUBROUTINE ELEC_COMPUTE_EX
                                                              
