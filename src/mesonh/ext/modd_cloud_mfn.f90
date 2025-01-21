!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     #################### 
      MODULE MODD_CLOUD_MF_n
!     ####################
!
!!****  *MODD_CLOUD_MF* - declaration of diagnostic variables of
!!           cloud creating using mass flux convection scheme            
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the 
!     characteritics of cloud
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------ 
!!      None  
!!
!!    REFERENCE
!!    ---------
!!      
!!      
!!
!!    AUTHOR
!!    ------
!!      J. Pergaud * Meteo-France *	
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CLOUD_MF_t
!
REAL, DIMENSION(:,:,:), POINTER :: XCF_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRC_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRI_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XHLC_HRC_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XHLC_HCF_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XHLI_HRI_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XHLI_HCF_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XWEIGHT_MF_CLOUD=>NULL()
!
END TYPE CLOUD_MF_t

TYPE(CLOUD_MF_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CLOUD_MF_MODEL

REAL, DIMENSION(:,:,:), POINTER :: XCF_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRC_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XRI_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XHLC_HRC_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XHLC_HCF_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XHLI_HRI_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XHLI_HCF_MF=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XWEIGHT_MF_CLOUD=>NULL()

CONTAINS

SUBROUTINE CLOUD_MF_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
CLOUD_MF_MODEL(KFROM)%XCF_MF=>XCF_MF
CLOUD_MF_MODEL(KFROM)%XRC_MF=>XRC_MF
CLOUD_MF_MODEL(KFROM)%XRI_MF=>XRI_MF
CLOUD_MF_MODEL(KFROM)%XHLC_HRC_MF=>XHLC_HRC_MF
CLOUD_MF_MODEL(KFROM)%XHLC_HCF_MF=>XHLC_HCF_MF
CLOUD_MF_MODEL(KFROM)%XHLI_HRI_MF=>XHLI_HRI_MF
CLOUD_MF_MODEL(KFROM)%XHLI_HCF_MF=>XHLI_HCF_MF
CLOUD_MF_MODEL(KFROM)%XWEIGHT_MF_CLOUD=>XWEIGHT_MF_CLOUD
!
! Current model is set to model KTO
XCF_MF=>CLOUD_MF_MODEL(KTO)%XCF_MF
XRC_MF=>CLOUD_MF_MODEL(KTO)%XRC_MF
XRI_MF=>CLOUD_MF_MODEL(KTO)%XRI_MF
XHLC_HRC_MF=>CLOUD_MF_MODEL(KTO)%XHLC_HRC_MF
XHLC_HCF_MF=>CLOUD_MF_MODEL(KTO)%XHLC_HCF_MF
XHLI_HRI_MF=>CLOUD_MF_MODEL(KTO)%XHLI_HRI_MF
XHLI_HCF_MF=>CLOUD_MF_MODEL(KTO)%XHLI_HCF_MF
XWEIGHT_MF_CLOUD=>CLOUD_MF_MODEL(KTO)%XWEIGHT_MF_CLOUD

END SUBROUTINE CLOUD_MF_GOTO_MODEL

END MODULE MODD_CLOUD_MF_n
