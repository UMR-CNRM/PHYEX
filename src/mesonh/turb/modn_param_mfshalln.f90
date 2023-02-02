!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    MODIFICATIONS
!!    -------------
!!       10/16 R.Honnert Update with AROME
!!       01/2019 R.Honnert add parameters for the reduction of mass-flux surface closure with resolution
!-----------------------------------------------------------------
!     ###########################
      MODULE MODN_PARAM_MFSHALL_n
!     ###########################
!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAM_MFSHALL_n, ONLY: &
         XIMPL_MF_n => XIMPL_MF, &
         CMF_UPDRAFT_n => CMF_UPDRAFT, &
         CMF_CLOUD_n => CMF_CLOUD, &
         LMIXUV_n => LMIXUV, &
         LMF_FLX_n => LMF_FLX, &  
         XALP_PERT_n => XALP_PERT, &
         XABUO_n => XABUO, &
         XBENTR_n => XBENTR, &
         XBDETR_n => XBDETR, &
         XCMF_n => XCMF, &
         XENTR_MF_n => XENTR_MF, &
         XCRAD_MF_n => XCRAD_MF, &
         XENTR_DRY_n => XENTR_DRY, &
         XDETR_DRY_n => XDETR_DRY, &
         XDETR_LUP_n => XDETR_LUP, &
         XKCF_MF_n => XKCF_MF, &
         XKRC_MF_n => XKRC_MF, &
         XTAUSIGMF_n => XTAUSIGMF, &
         XPRES_UV_n => XPRES_UV, &
         XALPHA_MF_n => XALPHA_MF, &
         XSIGMA_MF_n => XSIGMA_MF, &
         XFRAC_UP_MAX_n => XFRAC_UP_MAX, &
         XA1_n => XA1, &
         XB_n => XB, &
         XC_n => XC, &
         XBETA1_n => XBETA1, &
         XR_n => XR, &
         XLAMBDA_MF_n => XLAMBDA_MF, &
         LGZ_n => LGZ, &
         XGZ_n => XGZ
!
IMPLICIT NONE
!
REAL             ,SAVE  :: XIMPL_MF  
CHARACTER (LEN=4),SAVE  :: CMF_UPDRAFT
CHARACTER (LEN=4),SAVE  :: CMF_CLOUD
LOGICAL   ,SAVE  :: LMIXUV   
LOGICAL   ,SAVE  :: LMF_FLX           
!
REAL,SAVE          :: XALP_PERT   ! coefficient for the perturbation of
                                ! theta_l and r_t at the first level of 
                                ! the updraft
REAL,SAVE          ::    XABUO    ! coefficient of the buoyancy term in the w_up equation
REAL,SAVE          ::    XBENTR   ! coefficient of the entrainment term in the w_up equation
REAL,SAVE          ::    XBDETR   ! coefficient of the detrainment term in the w_up equation
REAL,SAVE          ::    XCMF     ! coefficient for the mass flux at the first level 
                                ! of the updraft (closure)
REAL,SAVE          :: XENTR_MF    ! entrainment constant (m/Pa) = 0.2 (m) 
REAL,SAVE          :: XCRAD_MF    ! cloud radius in cloudy part
REAL,SAVE          :: XENTR_DRY   ! coefficient for entrainment in dry part 
REAL,SAVE          :: XDETR_DRY   ! coefficient for detrainment in dry part
REAL,SAVE          :: XDETR_LUP   ! coefficient for detrainment in dry part
REAL,SAVE          :: XKCF_MF     ! coefficient for cloud fraction
REAL,SAVE          :: XKRC_MF     ! coefficient for convective rc
REAL,SAVE          :: XTAUSIGMF
REAL,SAVE          :: XPRES_UV    ! coefficient for pressure term in wind
                                  ! mixing

REAL,SAVE          :: XALPHA_MF   ! coefficient for cloudy fraction
REAL,SAVE          :: XSIGMA_MF   ! coefficient for sigma computation

REAL,SAVE          :: XFRAC_UP_MAX! maximum Updraft fraction
!
!
! Tuning variables for RHCJ10 updraft :
!
REAL,SAVE          :: XA1 
REAL,SAVE          :: XB
REAL,SAVE          :: XC  
REAL,SAVE          :: XBETA1
!
! Tuning variables for RAHA updraft :
!
REAL,SAVE          :: XR
REAL,SAVE          :: XLAMBDA_MF
!
! Tuning variables for Grey Zone updraft :
!
LOGICAL,SAVE       :: LGZ
REAL,SAVE          :: XGZ
!
NAMELIST/NAM_PARAM_MFSHALLn/XIMPL_MF,CMF_UPDRAFT,CMF_CLOUD,LMIXUV,LMF_FLX,&
                            XALP_PERT,XABUO,XBENTR,XBDETR,XCMF,XENTR_MF,&
                            XCRAD_MF,XENTR_DRY,XDETR_DRY,XDETR_LUP,XKCF_MF,&
                            XKRC_MF,XTAUSIGMF,XPRES_UV,XALPHA_MF,XSIGMA_MF,&
                            XFRAC_UP_MAX,XA1,XB,XC,XBETA1,XR,XLAMBDA_MF,LGZ,XGZ


!
CONTAINS
!
SUBROUTINE INIT_NAM_PARAM_MFSHALLn
   XIMPL_MF = XIMPL_MF_n
   CMF_UPDRAFT = CMF_UPDRAFT_n
   CMF_CLOUD = CMF_CLOUD_n
   LMIXUV = LMIXUV_n
   LMF_FLX = LMF_FLX_n
   XALP_PERT = XALP_PERT_n
   XABUO = XABUO_n
   XBENTR = XBENTR_n
   XBDETR = XBDETR_n
   XCMF = XCMF_n
   XENTR_MF = XENTR_MF_n
   XCRAD_MF = XCRAD_MF_n
   XENTR_DRY = XENTR_DRY_n
   XDETR_DRY = XDETR_DRY_n
   XDETR_LUP = XDETR_LUP_n
   XKCF_MF = XKCF_MF_n
   XKRC_MF = XKRC_MF_n
   XTAUSIGMF = XTAUSIGMF_n
   XPRES_UV = XPRES_UV_n
   XALPHA_MF = XALPHA_MF_n
   XSIGMA_MF = XSIGMA_MF_n
   XFRAC_UP_MAX = XFRAC_UP_MAX_n
   XA1 = XA1_n
   XB = XB_n
   XC = XC_n
   XBETA1 = XBETA1_n
   XR = XR_n
   XLAMBDA_MF = XLAMBDA_MF_n
   LGZ = LGZ_n
   XGZ = XGZ_n

END SUBROUTINE INIT_NAM_PARAM_MFSHALLn

SUBROUTINE UPDATE_NAM_PARAM_MFSHALLn
   XIMPL_MF_n = XIMPL_MF
   CMF_UPDRAFT_n = CMF_UPDRAFT
   CMF_CLOUD_n = CMF_CLOUD
   LMIXUV_n = LMIXUV
   LMF_FLX_n = LMF_FLX
   XALP_PERT_n = XALP_PERT
   XABUO_n = XABUO
   XBENTR_n = XBENTR
   XBDETR_n = XBDETR
   XCMF_n = XCMF
   XENTR_MF_n = XENTR_MF
   XCRAD_MF_n = XCRAD_MF
   XENTR_DRY_n = XENTR_DRY
   XDETR_DRY_n = XDETR_DRY
   XDETR_LUP_n = XDETR_LUP
   XKCF_MF_n = XKCF_MF
   XKRC_MF_n = XKRC_MF
   XTAUSIGMF_n = XTAUSIGMF
   XPRES_UV_n = XPRES_UV
   XALPHA_MF_n = XALPHA_MF
   XSIGMA_MF_n = XSIGMA_MF
   XFRAC_UP_MAX_n = XFRAC_UP_MAX
   XA1_n = XA1
   XB_n = XB
   XC_n = XC
   XBETA1_n = XBETA1
   XR_n = XR
   XLAMBDA_MF_n = XLAMBDA_MF
   LGZ_n = LGZ
   XGZ_n = XGZ

END SUBROUTINE UPDATE_NAM_PARAM_MFSHALLn

END MODULE MODN_PARAM_MFSHALL_n
