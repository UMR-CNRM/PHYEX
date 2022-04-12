!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! masdev4_7 BUG1 2007/06/15 17:47:18
!-----------------------------------------------------------------
!      ########################
       MODULE MODI_INI_RADAR
!      ########################
!
INTERFACE
      SUBROUTINE INI_RADAR (HPRISTINE_ICE  )
!
CHARACTER (LEN=4), INTENT(IN) :: HPRISTINE_ICE ! Indicator of ice crystal characteristics
!
!
END SUBROUTINE INI_RADAR
!
END INTERFACE
!
END MODULE MODI_INI_RADAR
!     ###########################################################
      SUBROUTINE INI_RADAR ( HPRISTINE_ICE )
!     ###########################################################
!
!!****  *INI_RADAR * 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used to 
!!    compute radar reflectivity (radar_rain_ice.f90 or radar_simulator.f90)
!!    for DIAG after PREP_REAL_CASE with AROME file (CCLOUD=NONE)
!!
!!**  METHOD
!!    ------
!!      The constants useful to radar are initialized to their 
!!      numerical values as in ini_rain_ice.f90 for ICE3 
!!
!!    EXTERNAL
!!    --------
!!      GAMMA    :  gamma function
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XP00                 ! Reference pressure
!!        XRD                  ! Gaz constant for dry air
!!        XRHOLW               ! Liquid water density
!!      Module MODD_RAIN_ICE_DESCR
!!
!!
!!    AUTHOR
!!    ------
!!      G. TANGUY * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    27/10/2009
!!      P.Scheffknecht 22/04/2015: test missing on already allocated XRTMIN 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_RAIN_ICE_DESCR
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
CHARACTER (LEN=4), INTENT(IN)       :: HPRISTINE_ICE    ! Indicator of ice crystal caracteristics
!
!-------------------------------------------------------------------------------
!
!
!
!*      1.1    Raindrop characteristics
!
!
!
XAR = (XPI/6.0)*XRHOLW
XBR = 3.0
XCR = 842.
XDR = 0.8
XCCR = 8.E6   
!
!*       1.2    Ice crystal characteristics
!
!
SELECT CASE (HPRISTINE_ICE)
  CASE('PLAT')
    XAI = 0.82      ! Plates
    XBI = 2.5       ! Plates 
    XC_I = 800.     ! Plates
    XDI = 1.0       ! Plates
  CASE('COLU')
    XAI = 2.14E-3   ! Columns
    XBI = 1.7       ! Columns
    XC_I = 2.1E5    ! Columns
    XDI = 1.585     ! Columns
  CASE('BURO')
    XAI = 44.0      ! Bullet rosettes
    XBI = 3.0       ! Bullet rosettes
    XC_I = 4.3E5    ! Bullet rosettes
    XDI = 1.663     ! Bullet rosettes
END SELECT
!
!
!*       1.3    Snowflakes/aggregates characteristics
!
!
XAS = 0.02
XBS = 1.9
XCS = 5.1
XDS = 0.27
XCCS = 5.0
XCXS = 1.0
!
!*      1.4    Graupel/Frozen drop characteristics
!
!
XAG = 19.6 
XBG = 2.8 
XCG = 124. 
XDG = 0.66 
XCCG = 5.E5
XCXG = -0.5
!
!*       1.5    Hailstone characteristics
!
!
XAH = 470.
XBH = 3.0
XCH = 207.
XDH = 0.64
XCCH = 4.E4 
XCXH = -1.0 
!
!-------------------------------------------------------------------------------
!
!*       2.     DIMENSIONAL DISTRIBUTIONS OF THE SPECIES
!	            ----------------------------------------
!
!*       2.1    Raindrops distribution
!
XALPHAR = 1.0  ! Exponential law
XNUR    = 1.0  ! Exponential law
!
!*       2.2    Ice crystal distribution
!
XALPHAI = 3.0  ! Gamma law for the ice crystal volume
XNUI    = 3.0  ! Gamma law with little dispersion
!
XALPHAS = 1.0  ! Exponential law
XNUS    = 1.0  ! Exponential law
!
XALPHAG = 1.0  ! Exponential law
XNUG    = 1.0  ! Exponential law
!
XALPHAH = 1.0  ! Gamma law
XNUH    = 8.0  ! Gamma law with little dispersion
!
!*       2.3    Constants for shape parameter
!
XLBEXR = 1.0/(-1.0-XBR)
XLBR   = ( XAR*XCCR*MOMG(XALPHAR,XNUR,XBR) )**(-XLBEXR)
!
XLBEXI = 1.0/(-XBI)
XLBI   = ( XAI*MOMG(XALPHAI,XNUI,XBI) )**(-XLBEXI)
!
XLBEXS = 1.0/(XCXS-XBS)
XLBS   = ( XAS*XCCS*MOMG(XALPHAS,XNUS,XBS) )**(-XLBEXS)
!
XLBEXG = 1.0/(XCXG-XBG)
XLBG   = ( XAG*XCCG*MOMG(XALPHAG,XNUG,XBG) )**(-XLBEXG)
!
XLBEXH = 1.0/(XCXH-XBH)
XLBH   = ( XAH*XCCH*MOMG(XALPHAH,XNUH,XBH) )**(-XLBEXH)
!
!*       2.4    Minimal values allowed for the mixing ratios
! ICE3
IF(.NOT.ALLOCATED(XRTMIN)) ALLOCATE( XRTMIN(6) )
!
XRTMIN(1) = 1.0E-20
XRTMIN(2) = 1.0E-20
XRTMIN(3) = 1.0E-20
XRTMIN(4) = 1.0E-20
XRTMIN(5) = 1.0E-15
XRTMIN(6) = 1.0E-15

!
CONTAINS
!
!------------------------------------------------------------------------------
!
  FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
!
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law
!
  USE MODI_GAMMA

  IMPLICIT NONE

  REAL     :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL     :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL     :: PP     ! order of the moment
  REAL     :: PMOMG  ! result: moment of order ZP

!------------------------------------------------------------------------------


  PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)

  END FUNCTION MOMG

!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INI_RADAR


