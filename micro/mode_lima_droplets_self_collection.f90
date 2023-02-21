!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_DROPLETS_SELF_COLLECTION
  IMPLICIT NONE
CONTAINS
!     ######################################################################
  SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION (LDCOMPUTE,                      &
                                            PRHODREF,                       &
                                            PCCT, PLBDC3,                   &
                                            P_CC_SELF                       )
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the self-collection of cloud droplets rate
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,      ONLY : XCTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : XSELFC
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT     ! Cloud water C. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC3   ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CC_SELF
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PCCT)) :: ZW ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Cloud droplets self collection
!	        ------------------------------
!
!
P_CC_SELF(:)=0.
!
WHERE( PCCT(:)>XCTMIN(2) .AND. LDCOMPUTE(:) )
   ZW(:) = XSELFC*(PCCT(:)/PLBDC3(:))**2 * PRHODREF(:) ! analytical integration
   P_CC_SELF(:) = - ZW(:)
END WHERE
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION
END MODULE MODE_LIMA_DROPLETS_SELF_COLLECTION
