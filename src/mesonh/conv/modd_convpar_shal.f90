!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######spl
      MODULE MODD_CONVPAR_SHAL
!     ########################
!
!!****  *MODD_CONVPAR_SHAL* - Declaration of convection constants 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to declare  the 
!!      constants in the deep convection parameterization.    
!!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_CONVPAR_SHAL)
!!          
!!    AUTHOR
!!    ------
!!      P. Bechtold   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96                      
!!   Last modified  04/10/98
!!      E. Bazile   05/05/09
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
!
REAL, SAVE :: XA25        ! 25 km x 25 km reference grid area
!
REAL, SAVE :: XCRAD       ! cloud radius 
REAL, SAVE :: XCTIME_SHAL ! convective adjustment time
REAL, SAVE :: XCDEPTH     ! minimum necessary cloud depth
REAL, SAVE :: XCDEPTH_D   ! maximum allowed cloud thickness
REAL, SAVE :: XDTPERT     ! add small Temp perturb. at LCL
REAL, SAVE :: XATPERT     ! Parameter for temp Perturb
REAL, SAVE :: XBTPERT     ! Parameter for temp Perturb
                          ! (XATPERT* TKE/Cp + XBTPERT) * XDTPERT
REAL, SAVE :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)  
!
REAL, SAVE :: XZLCL       ! maximum allowed allowed height 
                          ! difference between departure level and surface
REAL, SAVE :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL, SAVE :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL, SAVE :: XNHGAM      ! accounts for non-hydrost. pressure 
			  ! in buoyancy term of w equation
                          ! = 2 / (1+gamma)
REAL, SAVE :: XTFRZ1      ! begin of freezing interval
REAL, SAVE :: XTFRZ2      ! end of freezing interval
!
!
REAL, SAVE :: XSTABT      ! factor to assure stability in  fractional time
                          ! integration, routine CONVECT_CLOSURE
REAL, SAVE :: XSTABC      ! factor to assure stability in CAPE adjustment,
                          !  routine CONVECT_CLOSURE
REAL, SAVE :: XAW,XBW     ! Parameters for WLCL = XAW * W + XBW 
LOGICAL, SAVE :: LLSMOOTH ! Default=TRUE but not necessary 
!
END MODULE MODD_CONVPAR_SHAL 
