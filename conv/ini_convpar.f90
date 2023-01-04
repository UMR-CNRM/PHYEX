!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 init 2006/05/18 13:07:25
!-----------------------------------------------------------------
!       #####################
        MODULE MODI_INI_CONVPAR
!       #####################
!
INTERFACE
!
SUBROUTINE INI_CONVPAR
END SUBROUTINE INI_CONVPAR
!
END INTERFACE
!
END MODULE MODI_INI_CONVPAR
!
!
!
!     ######################
      SUBROUTINE INI_CONVPAR
!     ######################
!
!!****  *INI_CONVPAR * - routine to initialize the constants modules 
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize  the constants
!     stored in  modules MODD_CONVPAR, MODD_CST, MODD_CONVPAREXT.
!      
!
!!**  METHOD
!!    ------
!!      The deep convection constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CONVPAR, routine INI_CONVPAR)
!!      
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Last modified  15/04/98 adapted for ARPEGE
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAR
!
IMPLICIT NONE
!  
!-------------------------------------------------------------------------------
!
!*       1.    Set the thermodynamical and numerical constants for
!              the deep convection parameterization
!              ---------------------------------------------------
!
!
XA25     = 625.E6    ! 25 km x 25 km reference grid area
!
XCRAD    = 1500.     ! cloud radius 
XCDEPTH  = 2.5E3      ! minimum necessary cloud depth
XENTR    = 0.03      ! entrainment constant (m/Pa) = 0.2 (m)  
!
XZLCL    = 3.5E3     ! maximum allowed allowed height 
                     ! difference between the surface and the LCL
XZPBL    = 60.E2     ! minimum mixed layer depth to sustain convection
XWTRIG   = 6.00      ! constant in vertical velocity trigger
!
!
XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure 
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
XTFRZ1   = 268.16    ! begin of freezing interval
XTFRZ2   = 248.16    ! end of freezing interval
!
XRHDBC   = 0.9       ! relative humidity below cloud in downdraft

XRCONV   = 0.015     ! constant in precipitation conversion 
XSTABT   = 0.75      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
XUSRDPTH = 165.E2    ! pressure thickness used to compute updraft
                     ! moisture supply rate for downdraft
XMELDPTH = 100.E2    ! layer (Pa) through which precipitation melt is
                     ! allowed below downdraft
XUVDP    = 0.7       ! constant for pressure perturb in momentum transport
!
!
END SUBROUTINE INI_CONVPAR 
