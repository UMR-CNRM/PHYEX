!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############
      MODULE MODD_CST      
!     ###############
!
!!****  *MODD_CST* - declaration of Physic constants 
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the 
!     Physics constants.    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_CST)
!!          
!!    AUTHOR
!!    ------
!!      V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    16/05/94  
!!      J. Stein    02/01/95  add xrholw                    
!!      J.-P. Pinty 13/12/95  add XALPI,XBETAI,XGAMI
!!      J. Stein    25/07/97  add XTH00                    
!!      V. Masson   05/10/98  add XRHOLI
!!      C. Mari     31/10/00  add NDAYSEC
!!      V. Masson   01/03/03  add conductivity of ice
!!      J.Escobar : 10/2017 : for real*4 , add XMNH_HUGE_12_LOG
!!      J.L. Redelsperger 03/2021  add constants for ocean penetrating solar
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
REAL,SAVE :: XPI                ! Pi
!
REAL,SAVE :: XDAY,XSIYEA,XSIDAY ! day duration, sideral year duration,
                                ! sideral day duration
!
REAL,SAVE :: XKARMAN            ! von karman constant
REAL,SAVE :: XLIGHTSPEED        ! light speed
REAL,SAVE :: XPLANCK            ! Planck constant
REAL,SAVE :: XBOLTZ             ! Boltzman constant 
REAL,SAVE :: XAVOGADRO          ! Avogadro number
!
REAL,SAVE :: XRADIUS,XOMEGA     ! Earth radius, earth rotation
REAL,SAVE :: XG                 ! Gravity constant
!
REAL,SAVE :: XP00               ! Reference pressure
REAL,SAVE :: XP00OCEAN          ! Reference pressure for ocean model
REAL,SAVE :: XRH00OCEAN         ! Reference density for ocean model
!
REAL,SAVE :: XSTEFAN,XI0        ! Stefan-Boltzman constant, solar constant
!
REAL,SAVE :: XMD,XMV            ! Molar mass of dry air and molar mass of vapor
REAL,SAVE :: XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
REAL,SAVE :: XEPSILO            ! XMV/XMD
REAL,SAVE :: XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
REAL,SAVE :: XRHOLW             ! Volumic mass of liquid water
REAL,SAVE :: XCL,XCI            ! Cl (liquid), Ci (ice)
REAL,SAVE :: XTT                ! Triple point temperature
REAL,SAVE :: XLVTT              ! Vaporization heat constant
REAL,SAVE :: XLSTT              ! Sublimation heat constant
REAL,SAVE :: XLMTT              ! Melting heat constant
REAL,SAVE :: XESTT              ! Saturation vapor pressure  at triple point
                                ! temperature  
REAL,SAVE :: XALPW,XBETAW,XGAMW ! Constants for saturation vapor 
                                !  pressure  function 
REAL,SAVE :: XALPI,XBETAI,XGAMI ! Constants for saturation vapor
                                !  pressure  function over solid ice
REAL,SAVE :: XCONDI             ! thermal conductivity of ice (W m-1 K-1)
REAL,SAVE :: XALPHAOC           ! thermal expansion coefficient for ocean (K-1)
REAL,SAVE :: XBETAOC             ! Haline contraction coeff for ocean (S-1)
REAL,SAVE :: XTH00              ! reference value  for the potential temperature
REAL,SAVE :: XTH00OCEAN         ! Ref value for pot temp in ocean model
REAL,SAVE :: XSA00OCEAN         ! Ref value for SAlinity in ocean model 
REAL,SAVE :: XROC=0.69! 3 coeffs for SW penetration in  Ocean (Hoecker et al)
REAL,SAVE :: XD1=1.1
REAL,SAVE :: XD2=23.
! Values used in SURFEX CMO
!REAL,SAVE :: XROC=0.58
!REAL,SAVE :: XD1=0.35
!REAL,SAVE :: XD2=23.

REAL,SAVE :: XRHOLI             ! Volumic mass of liquid water
!
INTEGER, SAVE :: NDAYSEC        ! Number of seconds in a day
!
!
!   Some machine precision value depending of real4/8 use  
!
REAL,SAVE     :: XMNH_TINY          ! minimum real on this machine
REAL,SAVE     :: XMNH_TINY_12       ! sqrt(minimum real on this machine)
REAL,SAVE     :: XMNH_EPSILON       ! minimum space with 1.0
REAL,SAVE     :: XMNH_HUGE          ! maximum real on this machine
REAL,SAVE     :: XMNH_HUGE_12_LOG   ! maximum log(sqrt(real)) on this machine

REAL,SAVE     :: XEPS_DT            ! default value for DT test 
REAL,SAVE     :: XRES_FLAT_CART     ! default     flat&cart residual tolerance
REAL,SAVE     :: XRES_OTHER         ! default not flat&cart residual tolerance
REAL,SAVE     :: XRES_PREP          ! default     prep      residual tolerance

!
END MODULE MODD_CST
