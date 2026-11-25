!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!--------------------------------------------------------------------------
! ******* TO BE INCLUDED IN THE *CONTAINS* OF A SUBROUTINE OR FUNCTION ******
! => Don't use drHook !!!
! "gamma.h" must be included at the same time (MOMG uses GAMMA)
!
!     #############################################################
      ELEMENTAL FUNCTION MOMG(PALPHA, PNU, PP)  RESULT(PMOMG)
!     #############################################################
!
!!****  *MOMG* -  Moment of the generalized gamma law
!!
!!    PURPOSE
!!    -------
!!      Compute: G(p) = Gamma(nu + p/alpha) / Gamma(nu)
!!                    = M(p) * lambda^p
!!
!!    AUTHOR
!!    ------
!!      C. Barthe    * Laboratoire de l'Atmosphere et des Cyclones *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   26 Nov. 2009
!!      S. Riette  Nov. 2026: elemental function in momg.h, replaces momg.F90
!!                 (MOMG_X0D, MOMG_X1D) and the copies in several routines
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE



!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN) :: PALPHA ! first shape parameter of the dimensionnal distribution
REAL, INTENT(IN) :: PNU    ! second shape parameter of the dimensionnal distribution
REAL, INTENT(IN) :: PP     ! order of the moment
REAL             :: PMOMG  ! result: moment of order PP
!
!-------------------------------------------------------------------------------
!
PMOMG = GAMMA(PNU+PP/PALPHA) / GAMMA(PNU)
!
END FUNCTION MOMG
