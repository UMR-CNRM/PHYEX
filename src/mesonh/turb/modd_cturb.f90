!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/23 10:10:13
!-----------------------------------------------------------------
!      #######################
          MODULE MODD_CTURB
!      #######################
!
!!****   *MODD_CTURB*  - declaration of the turbulent scheme constants
!!
!!     PURPOSE
!!     -------
!        The purpose of this declarative module is to declare the 
!      turbulence scheme constants.
!
!!
!!**   IMPLICIT ARGUMENTS
!!     ------------------
!!       NONE
!!
!!     REFERENCE
!!     ---------
!!       Book 2 of Meso-NH documentation (MODD_CTURB)
!!       Book 1 of Meso-NH documentation (Chapter Turbulence)
!!
!!     AUTHOR
!!     ------
!1       Joan Cuxart         * INM  and Meteo-France *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original            08/08/94
!!     Nov 06, 2002 (V. Masson)  add XALPSBL and XASBL
!!     May 06                    Remove EPS
!!     Jan 2019     (Q. Rodier)  Remove XASBL
!----------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE
!
REAL,SAVE :: XCMFS        ! constant for the momentum flux due to shear   
REAL,SAVE :: XCMFB        ! constant for the momentum flux due to buoyancy
REAL,SAVE :: XCSHF        ! constant for the sensible heat flux 
REAL,SAVE :: XCHF         ! constant for the humidity flux 
REAL,SAVE :: XCTV         ! constant for the temperature variance
REAL,SAVE :: XCHV         ! constant for the humidity variance
REAL,SAVE :: XCHT1        ! first ct. for the humidity-temperature correlation
REAL,SAVE :: XCHT2        ! second ct. for the humidity-temperature correlation
!
REAL,SAVE :: XCPR1        ! first ct. for the turbulent Prandtl numbers
REAL,SAVE :: XCPR2        ! second ct. for the turbulent Prandtl numbers
REAL,SAVE :: XCPR3        ! third ct. for the turbulent Prandtl numbers
REAL,SAVE :: XCPR4        ! fourth ct. for the turbulent Prandtl numbers
REAL,SAVE :: XCPR5        ! fifth ct. for the turbulent Prandtl numbers
!
REAL,SAVE :: XCET         ! constant into the transport term of the TKE eq.
REAL,SAVE :: XCED         ! constant into the dissipation term of the TKE eq.
!
REAL,SAVE :: XCDP         ! ct. for the production term in the dissipation eq.
REAL,SAVE :: XCDD         ! ct. for the destruction term in the dissipation eq.
REAL,SAVE :: XCDT         ! ct. for the transport term in the dissipation eq.
!
REAL,SAVE :: XTKEMIN      ! mimimum value for the TKE
REAL,SAVE :: XRM17        ! Rodier et al 2017 constant in shear term for mixing length
!
REAL,SAVE :: XLINI        ! initial value for BL mixing length
REAL,SAVE :: XLINF        ! to prevent division by zero in the BL algorithm
!
REAL,SAVE :: XALPSBL      ! constant linking TKE and friction velocity in the SBL
!
REAL,SAVE :: XCEP         ! Constant for wind pressure-correlations
REAL,SAVE :: XA0          ! Constant a0 for wind pressure-correlations
REAL,SAVE :: XA2          ! Constant a2 for wind pressure-correlations
REAL,SAVE :: XA3          ! Constant a3 for wind pressure-correlations
REAL,SAVE :: XA5          ! Constant a5 for temperature pressure-correlations
REAL,SAVE :: XCTD         ! Constant for temperature and vapor dissipation
REAL,SAVE :: XCTP         ! Constant for temperature and vapor pressure-correlations
!
REAL,SAVE :: XPHI_LIM     ! Threshold value for Phi3 and Psi3
REAL,SAVE :: XSBL_O_BL    ! SBL height / BL height ratio
REAL,SAVE :: XFTOP_O_FSURF! Fraction of surface (heat or momentum) flux used to define top of BL
!
END MODULE MODD_CTURB
