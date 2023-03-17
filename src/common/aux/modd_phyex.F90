MODULE MODD_PHYEX
!
!> @file 
!!    MODD_PHYEX - decalration of the PHYEX structure gathering all the parametrisation strucutres of PHYEX
!! 
!!    PURPOSE 
!!    ------- 
!!       The purpose of this declarative module is to declare the 
!!       the PHYEX type that allows to gather all the different structures used
!!       by the paramteriation available in PHYEX
!! 
!!    AUTHOR
!!    ------
!!      S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    Mar 2023
!!
!-------------------------------------------------------------------------------
!
USE MODD_CST, ONLY: CST_t
USE MODD_PARAM_ICE_n, ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_CLOUDPAR_N, ONLY: CLOUDPAR_t
USE MODD_PARAM_MFSHALL_N, ONLY: PARAM_MFSHALL_t
USE MODD_TURB_n, ONLY: TURB_t
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_NEB_n, ONLY: NEB_t
USE MODD_MISC, ONLY: MISC_t
!
IMPLICIT NONE
!
TYPE PHYEX_t
  ! Structures for the different parametrisations
  TYPE(CST_t)            :: CST              !< Physical constants
  TYPE(PARAM_ICE_t)      :: PARAM_ICEN       !< Control parameters for microphysics
  TYPE(RAIN_ICE_DESCR_t) :: RAIN_ICE_DESCRN  !< Microphysical descriptive constants
  TYPE(RAIN_ICE_PARAM_t) :: RAIN_ICE_PARAMN  !< Microphysical factors
  TYPE(CLOUDPAR_t)       :: CLOUDPARN        !< Some other microphysical values
  TYPE(PARAM_MFSHALL_t)  :: PARAM_MFSHALLN   !< Mass flux scheme free parameters
  TYPE(CSTURB_t)         :: CSTURB           !< Turbulence scheme constants
  TYPE(TURB_t)           :: TURBN            !< Turbulence scheme constants set by namelist
  TYPE(NEB_t)            :: NEBN             !< Cloud scheme constants
  !
  ! Supplementary strucuture to hold model specific values
  TYPE(MISC_t)           :: MISC             !< Model specific values
END TYPE PHYEX_t
!
END MODULE MODD_PHYEX
