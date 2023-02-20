!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_GRAUPEL_DEPOSITION
  IMPLICIT NONE
CONTAINS
!     ###########################################################################
  SUBROUTINE LIMA_GRAUPEL_DEPOSITION (LDCOMPUTE, PRHODREF,                        &
                                      PRGT, PCGT, PSSI, PLBDG, PAI, PCJ, PLSFACT, &
                                      P_TH_DEPG, P_RG_DEPG                        )
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
USE MODD_PARAM_LIMA_MIXED, ONLY : X0DEPG, XEX0DEPG, X1DEPG, XEX1DEPG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT     ! graupel mr
REAL, DIMENSION(:),   INTENT(IN)    :: PCGT     ! graupel conc
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI     ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDG    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PAI      ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ      ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT  ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_DEPG
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RG_DEPG
!
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Deposition of vapour on graupel
!	        -------------------------------
!
P_TH_DEPG(:) = 0.0
P_RG_DEPG(:) = 0.0
WHERE ( PRGT(:)>XRTMIN(6) .AND. PCGT(:)>XCTMIN(6) .AND. LDCOMPUTE(:) )
   P_RG_DEPG(:) = PSSI(:) / PAI(:) * PCGT(:) *                      &
                ( X0DEPG*PLBDG(:)**XEX0DEPG + X1DEPG*PCJ(:)*PLBDG(:)**XEX1DEPG )
   P_TH_DEPG(:) = P_RG_DEPG(:)*PLSFACT(:)
END WHERE
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_GRAUPEL_DEPOSITION
END MODULE MODE_LIMA_GRAUPEL_DEPOSITION
