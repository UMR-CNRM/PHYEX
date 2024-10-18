!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_HAIL_DEPOSITION
  IMPLICIT NONE
CONTAINS
!     ###########################################################################
  SUBROUTINE LIMA_HAIL_DEPOSITION (KSIZE, ODCOMPUTE, PRHODREF,                 &
                                   PRHT, PCHT, PSSI, PLBDH, PAI, PCJ, PLSFACT, &
                                   P_TH_DEPH, P_RH_DEPH                        )
!     ###########################################################################
!
!!    PURPOSE
!!    -------
!!      Deposition of water vapour on graupel
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
!       M. Taufour              07/2022 add concentration for snow, graupel, hail        
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN
USE MODD_PARAM_LIMA_MIXED, ONLY : X0DEPH, XEX0DEPH, X1DEPH, XEX1DEPH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHT     ! hail mr
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCHT     ! hail conc
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PSSI     ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDH    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PAI      ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCJ      ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT  ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_DEPH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RH_DEPH
!
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Deposition of vapour on graupel
!	        -------------------------------
!
P_TH_DEPH(:) = 0.0
P_RH_DEPH(:) = 0.0
WHERE ( PRHT(:)>XRTMIN(7) .AND. PCHT(:)>XCTMIN(7) .AND. ODCOMPUTE(:) )
   P_RH_DEPH(:) = PSSI(:) / PAI(:) * PCHT(:) *                      &
                ( X0DEPH*PLBDH(:)**XEX0DEPH + X1DEPH*PCJ(:)*PLBDH(:)**XEX1DEPH )
   P_TH_DEPH(:) = P_RH_DEPH(:)*PLSFACT(:)
END WHERE
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_HAIL_DEPOSITION
END MODULE MODE_LIMA_HAIL_DEPOSITION
