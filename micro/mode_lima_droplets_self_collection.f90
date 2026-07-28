!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_DROPLETS_SELF_COLLECTION
  IMPLICIT NONE
CONTAINS
!     ######################################################################
  SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION (LIMAP, LIMAW, KSIZE, ODCOMPUTE,               &
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
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,              INTENT(IN)    :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT     ! Cloud water C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC3   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_SELF
!
!*       0.2   Declarations of local variables :
!
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
REAL, DIMENSION(SIZE(PCCT)) :: ZW ! work arrays
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Cloud droplets self collection
!               ------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_SELF_COLLECTION', 0, ZHOOK_HANDLE)
P_CC_SELF(:)=0.
!
WHERE( PCCT(:)>LIMAP%XCTMIN(2) .AND. ODCOMPUTE(:) )
   ZW(:) = LIMAW%XSELFC*(PCCT(:)/PLBDC3(:))**2 * PRHODREF(:) ! analytical integration
   P_CC_SELF(:) = - ZW(:)
END WHERE
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_SELF_COLLECTION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPLETS_SELF_COLLECTION
END MODULE MODE_LIMA_DROPLETS_SELF_COLLECTION
