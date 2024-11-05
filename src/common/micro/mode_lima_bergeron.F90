!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_LIMA_BERGERON
  IMPLICIT NONE
  CONTAINS
!     #############################################################
    SUBROUTINE LIMA_BERGERON( LIMAP, LIMAC, KSIZE, ODCOMPUTE,                  &
                              PRCT, PRIT, PCIT, PLBDI,           &
                              PSSIW, PAI, PCJ, PLVFACT, PLSFACT, &
                              P_TH_BERFI, P_RC_BERFI             )
!     #############################################################
!
!!    PURPOSE
!!    -------
!!      Compute the Bergeron-Findeisen process rate
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
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN
USE MODD_PARAM_LIMA_COLD,  ONLY : XDI, X0DEPI, X2DEPI
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_t
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDI   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PSSIW   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PAI     ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCJ     ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLVFACT ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_BERFI
TYPE(PARAM_LIMA_COLD_t),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_t),INTENT(IN)::LIMAP
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_BERFI
!
!-------------------------------------------------------------------------------
!
! Bergeron-Findeisen process
!
P_TH_BERFI(:) = 0.0
P_RC_BERFI(:) = 0.0
!
WHERE( (PRCT(:)>XRTMIN(2)) .AND. (PRIT(:)>XRTMIN(4)) .AND. (PCIT(:)>XCTMIN(4)) .AND. ODCOMPUTE(:))
   P_RC_BERFI(:) = - ( PSSIW(:) / PAI(:) ) * PCIT(:) *        &
        ( X0DEPI/PLBDI(:)+X2DEPI*PCJ(:)*PCJ(:)/PLBDI(:)**(XDI+2.0) )
   P_TH_BERFI(:) = - P_RC_BERFI(:)*(PLSFACT(:)-PLVFACT(:))
END WHERE
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_BERGERON
END MODULE MODE_LIMA_BERGERON
