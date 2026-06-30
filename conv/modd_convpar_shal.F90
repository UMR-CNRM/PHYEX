!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
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

!$ACDC methods 

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
TYPE CONVPAR_SHAL
REAL :: XA25        ! 25 km x 25 km reference grid area
!
REAL :: XCRAD       ! cloud radius
REAL :: XCTIME_SHAL ! convective adjustment time
REAL :: XCDEPTH     ! minimum necessary cloud depth
REAL :: XCDEPTH_D   ! maximum allowed cloud thickness
REAL :: XDTPERT     ! add small Temp perturb. at LCL
REAL :: XATPERT     ! Parameter for temp Perturb
REAL :: XBTPERT     ! Parameter for temp Perturb
                    ! (XATPERT* TKE/Cp + XBTPERT) * XDTPERT
REAL :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)
!
REAL :: XZLCL       ! maximum allowed allowed height
                    ! difference between departure level and surface
REAL :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL :: XNHGAM      ! accounts for non-hydrost. pressure
                    ! in buoyancy term of w equation
                    ! = 2 / (1+gamma)
REAL :: XTFRZ1      ! begin of freezing interval
REAL :: XTFRZ2      ! end of freezing interval
!
!
REAL :: XSTABT      ! factor to assure stability in  fractional time
                    ! integration, routine CONVECT_CLOSURE
REAL :: XSTABC      ! factor to assure stability in CAPE adjustment,
                    !  routine CONVECT_CLOSURE
REAL :: XAW,XBW     ! Parameters for WLCL = XAW * W + XBW
LOGICAL :: LLSMOOTH ! Default=TRUE but not necessary
END TYPE CONVPAR_SHAL
!Keep global variables for parts of the code not ported to the type yet
REAL :: XA25        ! 25 km x 25 km reference grid area
!
REAL :: XCRAD       ! cloud radius
REAL :: XCTIME_SHAL ! convective adjustment time
REAL :: XCDEPTH     ! minimum necessary cloud depth
REAL :: XCDEPTH_D   ! maximum allowed cloud thickness
REAL :: XDTPERT     ! add small Temp perturb. at LCL
REAL :: XATPERT     ! Parameter for temp Perturb
REAL :: XBTPERT     ! Parameter for temp Perturb
                    ! (XATPERT* TKE/Cp + XBTPERT) * XDTPERT
REAL :: XENTR       ! entrainment constant (m/Pa) = 0.2 (m)
!
REAL :: XZLCL       ! maximum allowed allowed height
                    ! difference between departure level and surface
REAL :: XZPBL       ! minimum mixed layer depth to sustain convection
REAL :: XWTRIG      ! constant in vertical velocity trigger
!
!
REAL :: XNHGAM      ! accounts for non-hydrost. pressure
                    ! in buoyancy term of w equation
                    ! = 2 / (1+gamma)
REAL :: XTFRZ1      ! begin of freezing interval
REAL :: XTFRZ2      ! end of freezing interval
!
!
REAL :: XSTABT      ! factor to assure stability in  fractional time
                    ! integration, routine CONVECT_CLOSURE
REAL :: XSTABC      ! factor to assure stability in CAPE adjustment,
                    !  routine CONVECT_CLOSURE
REAL :: XAW,XBW     ! Parameters for WLCL = XAW * W + XBW
LOGICAL :: LLSMOOTH ! Default=TRUE but not necessary
!
INTERFACE INI_CONVPAR_SHAL
  MODULE PROCEDURE INI_CONVPAR_SHAL0
  MODULE PROCEDURE INI_CONVPAR_SHAL1
END INTERFACE



CONTAINS


SUBROUTINE INI_CONVPAR_SHAL1(CVP_SHAL)

USE YOMHOOK , ONLY : LHOOK, JPHOOK, DR_HOOK

TYPE (CONVPAR_SHAL), INTENT(OUT) :: CVP_SHAL

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('INI_CONVPAR_SHAL1',0,ZHOOK_HANDLE)

CVP_SHAL%XA25        = XA25
CVP_SHAL%XCRAD       = XCRAD
CVP_SHAL%XCTIME_SHAL = XCTIME_SHAL
CVP_SHAL%XCDEPTH     = XCDEPTH
CVP_SHAL%XCDEPTH_D   = XCDEPTH_D
CVP_SHAL%XDTPERT     = XDTPERT
CVP_SHAL%XATPERT     = XATPERT
CVP_SHAL%XBTPERT     = XBTPERT
CVP_SHAL%XENTR       = XENTR
CVP_SHAL%XZLCL       = XZLCL
CVP_SHAL%XZPBL       = XZPBL
CVP_SHAL%XWTRIG      = XWTRIG
CVP_SHAL%XNHGAM      = XNHGAM
CVP_SHAL%XTFRZ1      = XTFRZ1
CVP_SHAL%XTFRZ2      = XTFRZ2
CVP_SHAL%XSTABT      = XSTABT
CVP_SHAL%XSTABC      = XSTABC
CVP_SHAL%XAW         = XAW
CVP_SHAL%XBW         = XBW
CVP_SHAL%LLSMOOTH    = LLSMOOTH

IF (LHOOK) CALL DR_HOOK('INI_CONVPAR_SHAL1',1,ZHOOK_HANDLE)

END SUBROUTINE INI_CONVPAR_SHAL1

SUBROUTINE INI_CONVPAR_SHAL0
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ###########################
!
!!****  *INI_CONVPAR * - routine to initialize the constants modules 
!!
!!    PURPOSE
!!    -------
!!       The purpose of this routine is to initialize  the constants
!!     stored in  modules MODD_CONVPAR_SHAL
!!      
!!
!!**  METHOD
!!    ------
!!      The shallow convection constants are set to their numerical values 
!!     
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CONVPAR_SHAL   : contains deep convection constants
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (module MODD_CONVPAR_SHAL, routine INI_CONVPAR)
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
!!                  05/05/09 E. Bazile
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('INI_CONVPAR_SHAL0',0,ZHOOK_HANDLE)
XA25     = 625.E6    ! 25 km x 25 km reference grid area
!
XCRAD       = 50.    ! cloud radius 
XCTIME_SHAL = 10800. ! convective adjustment time
XCDEPTH     = 0.5E3  ! minimum necessary shallow cloud depth
XCDEPTH_D   = 2.5E3  ! maximum allowed shallow cloud depth
XDTPERT     = .2     ! add small Temp perturbation at LCL
XATPERT     = 0.     ! 0.=original scheme , recommended = 1000. 
XBTPERT     = 1.     ! 1.=original scheme , recommended = 0.
!
XENTR    = 0.02      ! entrainment constant (m/Pa) = 0.2 (m)  
!
XZLCL    = 0.5E3     ! maximum allowed allowed height 
                     ! difference between the DPL and the surface
XZPBL    = 40.E2     ! minimum mixed layer depth to sustain convection
!
!
XNHGAM   = 1.3333    ! accounts for non-hydrost. pressure 
                     ! in buoyancy term of w equation
                     ! = 2 / (1+gamma)
XTFRZ1   = 268.16    ! begin of freezing interval
XTFRZ2   = 248.16    ! end of freezing interval
!

XSTABT   = 0.75      ! factor to assure stability in  fractional time
                     ! integration, routine CONVECT_CLOSURE
XSTABC   = 0.95      ! factor to assure stability in CAPE adjustment,
                     !  routine CONVECT_CLOSURE
XAW      = 0.        ! 0.= Original scheme , 1 = recommended 
XBW      = 1.        ! 1.= Original scheme , 0 = recommended
LLSMOOTH = .TRUE.
!
!
IF (LHOOK) CALL DR_HOOK('INI_CONVPAR_SHAL0',1,ZHOOK_HANDLE)
END SUBROUTINE INI_CONVPAR_SHAL0
!
END MODULE MODD_CONVPAR_SHAL
