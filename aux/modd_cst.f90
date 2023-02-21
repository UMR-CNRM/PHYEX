!MNH_LIC Copyright 1994-2023 CNRS, Meteo-France and Universite Paul Sabatier
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
!!      R. El Khatib 04/08/14 add pre-computed quantities
!!      J.Escobar : 10/2017 : for real*4 , add XMNH_HUGE_12_LOG
!  J.L. Redelsperger 03/2021: add constants for ocean penetrating solar
!  S. Riette      01/2022: introduction of a structure
!  P. Wautelet 20/05/2022: add RASTA cloud radar wavelength
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 

REAL, PARAMETER :: XLAM_CRAD = 3.154E-3 ! RASTA cloud radar wavelength (m) <=> 95.04 GHz

TYPE CST_t
  REAL :: XPI                ! Pi
  !
  REAL :: XDAY,XSIYEA,XSIDAY ! day duration, sideral year duration, sideral day duration
  !
  REAL :: XKARMAN            ! von karman constant
  REAL :: XLIGHTSPEED        ! light speed
  REAL :: XPLANCK            ! Planck constant
  REAL :: XBOLTZ             ! Boltzman constant
  REAL :: XAVOGADRO          ! Avogadro number
  !
  REAL :: XRADIUS,XOMEGA     ! Earth radius, earth rotation
  REAL :: XG                 ! Gravity constant
  !
  REAL :: XP00               ! Reference pressure
  REAL :: XP00OCEAN          ! Reference pressure for ocean model
  REAL :: XRH00OCEAN         ! Reference density for ocean model
  !
  REAL :: XSTEFAN,XI0        ! Stefan-Boltzman constant, solar constant
  !
  REAL :: XMD,XMV            ! Molar mass of dry air and molar mass of vapor
  REAL :: XRD,XRV            ! Gaz constant for dry air, gaz constant for vapor
  REAL :: XEPSILO            ! XMV/XMD
  REAL :: XCPD,XCPV          ! Cpd (dry air), Cpv (vapor)
  REAL :: XRHOLW             ! Volumic mass of liquid water
  REAL :: XCL,XCI            ! Cl (liquid), Ci (ice)
  REAL :: XTT                ! Triple point temperature
  REAL :: XLVTT              ! Vaporization heat constant
  REAL :: XLSTT              ! Sublimation heat constant
  REAL :: XLMTT              ! Melting heat constant
  REAL :: XESTT              ! Saturation vapor pressure  at triple point temperature
  REAL :: XALPW,XBETAW,XGAMW ! Constants for saturation vapor pressure function
  REAL :: XALPI,XBETAI,XGAMI ! Constants for saturation vapor pressure  function over solid ice
  REAL :: XCONDI             ! thermal conductivity of ice (W m-1 K-1)
  REAL :: XALPHAOC           ! thermal expansion coefficient for ocean (K-1)
  REAL :: XBETAOC            ! Haline contraction coeff for ocean (S-1)
  REAL :: XTH00              ! reference value  for the potential temperature
  REAL :: XTH00OCEAN         ! Ref value for pot temp in ocean model
  REAL :: XSA00OCEAN         ! Ref value for SAlinity in ocean model
  REAL :: XROC=0.69! 3 coeffs for SW penetration in  Ocean (Hoecker et al)
  REAL :: XD1=1.1
  REAL :: XD2=23.
  ! Values used in SURFEX CMO
  !REAL :: XROC=0.58
  !REAL :: XD1=0.35
  !REAL :: XD2=23.

  REAL :: XRHOLI             ! Volumic mass of ice
  !
  INTEGER :: NDAYSEC         ! Number of seconds in a day
  !
  REAL :: RDSRV              !  XRD/XRV
  REAL :: RDSCPD             !  XRD/XCPD
  REAL :: RINVXP00           !  1./XP00
  !
  !   Some machine precision value depending of real4/8 use
  !
  REAL :: XMNH_TINY          ! minimum real on this machine
  REAL :: XMNH_TINY_12       ! sqrt(minimum real on this machine)
  REAL :: XMNH_EPSILON       ! minimum space with 1.0
  REAL :: XMNH_HUGE          ! maximum real on this machine
  REAL :: XMNH_HUGE_12_LOG   ! maximum log(sqrt(real)) on this machine

  REAL :: XEPS_DT            ! default value for DT test
  REAL :: XRES_FLAT_CART     ! default     flat&cart residual tolerance
  REAL :: XRES_OTHER         ! default not flat&cart residual tolerance
  REAL :: XRES_PREP          ! default     prep      residual tolerance
END TYPE CST_t

TYPE(CST_t), TARGET, SAVE :: CST

REAL, POINTER :: XPI=>NULL()
REAL, POINTER :: XDAY=>NULL(), XSIYEA=>NULL(), XSIDAY=>NULL()
REAL, POINTER :: XKARMAN=>NULL()
REAL, POINTER :: XLIGHTSPEED=>NULL()
REAL, POINTER :: XPLANCK=>NULL()
REAL, POINTER :: XBOLTZ=>NULL()
REAL, POINTER :: XAVOGADRO=>NULL()
REAL, POINTER :: XRADIUS=>NULL(), XOMEGA=>NULL()
REAL, POINTER :: XG=>NULL()
REAL, POINTER :: XP00=>NULL()
REAL, POINTER :: XP00OCEAN=>NULL()
REAL, POINTER :: XRH00OCEAN=>NULL()
REAL, POINTER :: XSTEFAN=>NULL(), XI0=>NULL()
REAL, POINTER :: XMD=>NULL(), XMV=>NULL()
REAL, POINTER :: XRD=>NULL(), XRV=>NULL()
REAL, POINTER :: XEPSILO=>NULL()
REAL, POINTER :: XCPD=>NULL(), XCPV=>NULL()
REAL, POINTER :: XRHOLW=>NULL()
REAL, POINTER :: XCL=>NULL(), XCI=>NULL()
REAL, POINTER :: XTT=>NULL()
REAL, POINTER :: XLVTT=>NULL()
REAL, POINTER :: XLSTT=>NULL()
REAL, POINTER :: XLMTT=>NULL()
REAL, POINTER :: XESTT=>NULL()
REAL, POINTER :: XALPW=>NULL(), XBETAW=>NULL(), XGAMW=>NULL()
REAL, POINTER :: XALPI=>NULL(), XBETAI=>NULL(), XGAMI=>NULL()
REAL, POINTER :: XCONDI=>NULL()
REAL, POINTER :: XALPHAOC=>NULL()
REAL, POINTER :: XBETAOC=>NULL()
REAL, POINTER :: XTH00=>NULL()
REAL, POINTER :: XTH00OCEAN=>NULL()
REAL, POINTER :: XSA00OCEAN=>NULL()
REAL, POINTER :: XROC=>NULL()
REAL, POINTER :: XD1=>NULL()
REAL, POINTER :: XD2=>NULL()
REAL, POINTER :: XRHOLI=>NULL()
INTEGER, POINTER :: NDAYSEC=>NULL()
REAL, POINTER :: RDSRV=>NULL()
REAL, POINTER :: RDSCPD=>NULL()
REAL, POINTER :: RINVXP00=>NULL()
REAL, POINTER :: XMNH_TINY=>NULL()
REAL, POINTER :: XMNH_TINY_12=>NULL()
REAL, POINTER :: XMNH_EPSILON=>NULL()
REAL, POINTER :: XMNH_HUGE=>NULL()
REAL, POINTER :: XMNH_HUGE_12_LOG=>NULL()
REAL, POINTER :: XEPS_DT=>NULL()
REAL, POINTER :: XRES_FLAT_CART=>NULL()
REAL, POINTER :: XRES_OTHER=>NULL()
REAL, POINTER :: XRES_PREP=>NULL()
!
CONTAINS

SUBROUTINE CST_ASSOCIATE()
  IMPLICIT NONE
  XPI=>CST%XPI
  XDAY=>CST%XDAY
  XSIYEA=>CST%XSIYEA
  XSIDAY=>CST%XSIDAY
  XKARMAN=>CST%XKARMAN
  XLIGHTSPEED=>CST%XLIGHTSPEED
  XPLANCK=>CST%XPLANCK
  XBOLTZ=>CST%XBOLTZ
  XAVOGADRO=>CST%XAVOGADRO
  XRADIUS=>CST%XRADIUS
  XOMEGA=>CST%XOMEGA
  XG=>CST%XG
  XP00=>CST%XP00
  XP00OCEAN=>CST%XP00OCEAN
  XRH00OCEAN=>CST%XRH00OCEAN
  XSTEFAN=>CST%XSTEFAN
  XI0=>CST%XI0
  XMD=>CST%XMD
  XMV=>CST%XMV
  XRD=>CST%XRD
  XRV=>CST%XRV
  XEPSILO=>CST%XEPSILO
  XCPD=>CST%XCPD
  XCPV=>CST%XCPV
  XRHOLW=>CST%XRHOLW
  XCL=>CST%XCL
  XCI=>CST%XCI
  XTT=>CST%XTT
  XLVTT=>CST%XLVTT
  XLSTT=>CST%XLSTT
  XLMTT=>CST%XLMTT
  XESTT=>CST%XESTT
  XALPW=>CST%XALPW
  XBETAW=>CST%XBETAW
  XGAMW=>CST%XGAMW
  XALPI=>CST%XALPI
  XBETAI=>CST%XBETAI
  XGAMI=>CST%XGAMI
  XCONDI=>CST%XCONDI
  XALPHAOC=>CST%XALPHAOC
  XBETAOC=>CST%XBETAOC
  XTH00=>CST%XTH00
  XTH00OCEAN=>CST%XTH00OCEAN
  XSA00OCEAN=>CST%XSA00OCEAN
  XROC=>CST%XROC
  XD1=>CST%XD1
  XD2=>CST%XD2
  XRHOLI=>CST%XRHOLI
  NDAYSEC=>CST%NDAYSEC
  RDSRV=>CST%RDSRV
  RDSCPD=>CST%RDSCPD
  RINVXP00=>CST%RINVXP00
  XMNH_TINY=>CST%XMNH_TINY
  XMNH_TINY_12=>CST%XMNH_TINY_12
  XMNH_EPSILON=>CST%XMNH_EPSILON
  XMNH_HUGE=>CST%XMNH_HUGE
  XMNH_HUGE_12_LOG=>CST%XMNH_HUGE_12_LOG
  XEPS_DT=>CST%XEPS_DT
  XRES_FLAT_CART=>CST%XRES_FLAT_CART
  XRES_OTHER=>CST%XRES_OTHER
  XRES_PREP=>CST%XRES_PREP
END SUBROUTINE CST_ASSOCIATE
!
END MODULE MODD_CST

