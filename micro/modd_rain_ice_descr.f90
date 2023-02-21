!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##########################
      MODULE MODD_RAIN_ICE_DESCR
!     ##########################
!
!!****  *MODD_RAIN_ICE_DESCR* - declaration of the microphysical descriptive
!!                              constants for use in the warm and cold schemes.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the microphysical
!     constants. This includes the descriptive parameters for the raindrop and
!     the ice crystal habits and the parameters relevant of the dimensional
!     distributions.
!
!         m(D)    = XAx * D**XBx      : Mass-MaxDim relationship
!         v(D)    = XCx * D**XDx      : Fallspeed-MaxDim relationship
!         N(Lbda) = XCCx * Lbda**XCXx : NumberConc-Slopeparam relationship
!         XF0x, XF1x, XF2x            : Ventilation factors
!         XC1x                        : Shape parameter for deposition
!
!              and
!
!         XALPHAx, XNUx                        : Generalized GAMMA law
!         Lbda = XLBx * (r_x*rho_dref)**XLBEXx : Slope parameter of the
!                                                distribution law
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_RAIN_ICE_DESCR)
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/95
!!       J.-P. Pinty   29/11/02 add ICE4
!!       C. LAC     26/01/2012 : suppression de XCONC qui n'était pas utilisé
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
TYPE RAIN_ICE_DESCR_t
REAL :: XCEXVT               ! air density fall speed correction
!
REAL :: XAC,XBC,XCC,XDC                          ! Cloud droplet  charact.
REAL :: XAR,XBR,XCR,XDR,XCCR     ,XF0R,XF1R,XC1R ! Raindrop       charact.
REAL :: XAI,XBI,XC_I,XDI          ,XF0I,XF2I,XC1I ! Cloud ice      charact.
REAL :: XAS,XBS,XCS,XDS,XCCS,XCXS,XF0S,XF1S,XC1S ! Snow/agg.      charact.
REAL :: XAG,XBG,XCG,XDG,XCCG,XCXG,XF0G,XF1G,XC1G ! Graupel        charact.
REAL :: XAH,XBH,XCH,XDH,XCCH,XCXH,XF0H,XF1H,XC1H ! Hail           charact.
!
REAL :: XALPHAC,XNUC,XALPHAC2,XNUC2, XLBEXC      ! Cloud droplet  distribution parameters
REAL,DIMENSION(2) :: XLBC ! Cloud droplet distribution parameters
REAL :: XALPHAR,XNUR,XLBEXR,XLBR ! Raindrop       distribution parameters
REAL :: XALPHAI,XNUI,XLBEXI,XLBI ! Cloud ice      distribution parameters
REAL :: XALPHAS,XNUS,XLBEXS,XLBS,XNS ! Snow/agg.      distribution parameters
REAL :: XALPHAG,XNUG,XLBEXG,XLBG ! Graupel        distribution parameters
REAL :: XALPHAH,XNUH,XLBEXH,XLBH ! Hail           distribution parameters
!
REAL :: XFVELOS            ! factor for snow fall speed after Thompson (2008)
REAL :: XTRANS_MP_GAMMAS    ! coefficient to convert lambdas for gamma function
REAL :: XLBDAR_MAX,XLBDAS_MIN,XLBDAS_MAX,XLBDAG_MAX ! Max values allowed for the shape
                                              ! parameters (rain,snow,graupeln)
!
REAL,DIMENSION(:),ALLOCATABLE :: XRTMIN ! Min values allowed for the mixing ratios
REAL :: XCONC_SEA   ! Diagnostic concentration of droplets over sea
REAL :: XCONC_LAND  !  Diagnostic concentration of droplets over land
REAL :: XCONC_URBAN ! Diagnostic concentration of droplets over urban area
END TYPE RAIN_ICE_DESCR_t
!
TYPE(RAIN_ICE_DESCR_t), SAVE, TARGET :: RAIN_ICE_DESCR
!
REAL,DIMENSION(:),POINTER :: XLBC=>NULL(), XRTMIN=>NULL()
REAL, POINTER :: XCEXVT => NULL(), &
                 XAC => NULL(), &
                 XBC => NULL(), &
                 XCC => NULL(), &
                 XDC => NULL(), &
                 XAR => NULL(), &
                 XBR => NULL(), &
                 XCR => NULL(), &
                 XDR => NULL(), &
                 XCCR => NULL(), &
                 XF0R => NULL(), &
                 XF1R => NULL(), &
                 XC1R => NULL(), &
                 XAI => NULL(), &
                 XBI => NULL(), &
                 XC_I => NULL(), &
                 XDI => NULL(), &
                 XF0I => NULL(), &
                 XF2I => NULL(), &
                 XC1I => NULL(), &
                 XAS => NULL(), &
                 XBS => NULL(), &
                 XCS => NULL(), &
                 XDS => NULL(), &
                 XCCS => NULL(), &
                 XCXS => NULL(), &
                 XF0S => NULL(), &
                 XF1S => NULL(), &
                 XC1S => NULL(), &
                 XAG => NULL(), &
                 XBG => NULL(), &
                 XCG => NULL(), &
                 XDG => NULL(), &
                 XCCG => NULL(), &
                 XCXG => NULL(), &
                 XF0G => NULL(), &
                 XF1G => NULL(), &
                 XC1G => NULL(), &
                 XAH => NULL(), &
                 XBH => NULL(), &
                 XCH => NULL(), &
                 XDH => NULL(), &
                 XCCH => NULL(), &
                 XCXH => NULL(), &
                 XF0H => NULL(), &
                 XF1H => NULL(), &
                 XC1H => NULL(), &
                 XALPHAC => NULL(), &
                 XNUC => NULL(), &
                 XALPHAC2 => NULL(), &
                 XNUC2 => NULL(), &
                 XLBEXC => NULL(), &
                 XALPHAR => NULL(), &
                 XNUR => NULL(), &
                 XLBEXR => NULL(), &
                 XLBR => NULL(), &
                 XALPHAI => NULL(), &
                 XNUI => NULL(), &
                 XLBEXI => NULL(), &
                 XLBI => NULL(), &
                 XALPHAS => NULL(), &
                 XNUS => NULL(), &
                 XNS => NULL(), &
                 XLBEXS => NULL(), &
                 XLBS => NULL(), &
                 XALPHAG => NULL(), &
                 XNUG => NULL(), &
                 XLBEXG => NULL(), &
                 XLBG => NULL(), &
                 XALPHAH => NULL(), &
                 XNUH => NULL(), &
                 XLBEXH => NULL(), &
                 XLBH => NULL(), &
                 XLBDAR_MAX => NULL(), &
                 XLBDAS_MAX => NULL(), &
                 XLBDAG_MAX => NULL(), &
                 XCONC_SEA => NULL(), &
                 XCONC_LAND => NULL(), &
                 XCONC_URBAN => NULL(), &
                 XFVELOS => NULL(), &
                 XTRANS_MP_GAMMAS => NULL(), &
                 XLBDAS_MIN => NULL()
!
CONTAINS
SUBROUTINE RAIN_ICE_DESCR_ASSOCIATE()
  IMPLICIT NONE
  XCEXVT => RAIN_ICE_DESCR%XCEXVT
  XAC => RAIN_ICE_DESCR%XAC
  XBC => RAIN_ICE_DESCR%XBC
  XCC => RAIN_ICE_DESCR%XCC
  XDC => RAIN_ICE_DESCR%XDC
  XAR => RAIN_ICE_DESCR%XAR
  XBR => RAIN_ICE_DESCR%XBR
  XCR => RAIN_ICE_DESCR%XCR
  XDR => RAIN_ICE_DESCR%XDR
  XCCR => RAIN_ICE_DESCR%XCCR
  XF0R => RAIN_ICE_DESCR%XF0R
  XF1R => RAIN_ICE_DESCR%XF1R
  XC1R => RAIN_ICE_DESCR%XC1R
  XAI => RAIN_ICE_DESCR%XAI
  XBI => RAIN_ICE_DESCR%XBI
  XC_I => RAIN_ICE_DESCR%XC_I
  XDI => RAIN_ICE_DESCR%XDI
  XF0I => RAIN_ICE_DESCR%XF0I
  XF2I => RAIN_ICE_DESCR%XF2I
  XC1I => RAIN_ICE_DESCR%XC1I
  XAS => RAIN_ICE_DESCR%XAS
  XBS => RAIN_ICE_DESCR%XBS
  XCS => RAIN_ICE_DESCR%XCS
  XDS => RAIN_ICE_DESCR%XDS
  XCCS => RAIN_ICE_DESCR%XCCS
  XCXS => RAIN_ICE_DESCR%XCXS
  XF0S => RAIN_ICE_DESCR%XF0S
  XF1S => RAIN_ICE_DESCR%XF1S
  XC1S => RAIN_ICE_DESCR%XC1S
  XAG => RAIN_ICE_DESCR%XAG
  XBG => RAIN_ICE_DESCR%XBG
  XCG => RAIN_ICE_DESCR%XCG
  XDG => RAIN_ICE_DESCR%XDG
  XCCG => RAIN_ICE_DESCR%XCCG
  XCXG => RAIN_ICE_DESCR%XCXG
  XF0G => RAIN_ICE_DESCR%XF0G
  XF1G => RAIN_ICE_DESCR%XF1G
  XC1G => RAIN_ICE_DESCR%XC1G
  XAH => RAIN_ICE_DESCR%XAH
  XBH => RAIN_ICE_DESCR%XBH
  XCH => RAIN_ICE_DESCR%XCH
  XDH => RAIN_ICE_DESCR%XDH
  XCCH => RAIN_ICE_DESCR%XCCH
  XCXH => RAIN_ICE_DESCR%XCXH
  XF0H => RAIN_ICE_DESCR%XF0H
  XF1H => RAIN_ICE_DESCR%XF1H
  XC1H => RAIN_ICE_DESCR%XC1H
  XALPHAC => RAIN_ICE_DESCR%XALPHAC
  XNUC => RAIN_ICE_DESCR%XNUC
  XALPHAC2 => RAIN_ICE_DESCR%XALPHAC2
  XNUC2 => RAIN_ICE_DESCR%XNUC2
  XLBEXC => RAIN_ICE_DESCR%XLBEXC
  XALPHAR => RAIN_ICE_DESCR%XALPHAR
  XNUR => RAIN_ICE_DESCR%XNUR
  XLBEXR => RAIN_ICE_DESCR%XLBEXR
  XLBR => RAIN_ICE_DESCR%XLBR
  XALPHAI => RAIN_ICE_DESCR%XALPHAI
  XNUI => RAIN_ICE_DESCR%XNUI
  XLBEXI => RAIN_ICE_DESCR%XLBEXI
  XLBI => RAIN_ICE_DESCR%XLBI
  XALPHAS => RAIN_ICE_DESCR%XALPHAS
  XNUS => RAIN_ICE_DESCR%XNUS
  XLBEXS => RAIN_ICE_DESCR%XLBEXS
  XLBS => RAIN_ICE_DESCR%XLBS
  XALPHAG => RAIN_ICE_DESCR%XALPHAG
  XNUG => RAIN_ICE_DESCR%XNUG
  XLBEXG => RAIN_ICE_DESCR%XLBEXG
  XLBG => RAIN_ICE_DESCR%XLBG
  XALPHAH => RAIN_ICE_DESCR%XALPHAH
  XNUH => RAIN_ICE_DESCR%XNUH
  XLBEXH => RAIN_ICE_DESCR%XLBEXH
  XLBH => RAIN_ICE_DESCR%XLBH
  XLBDAR_MAX => RAIN_ICE_DESCR%XLBDAR_MAX
  XLBDAS_MAX => RAIN_ICE_DESCR%XLBDAS_MAX
  XLBDAG_MAX => RAIN_ICE_DESCR%XLBDAG_MAX
  XCONC_SEA => RAIN_ICE_DESCR%XCONC_SEA
  XCONC_LAND => RAIN_ICE_DESCR%XCONC_LAND
  XCONC_URBAN => RAIN_ICE_DESCR%XCONC_URBAN
  XNS => RAIN_ICE_DESCR%XNS
  XFVELOS => RAIN_ICE_DESCR%XFVELOS
  XTRANS_MP_GAMMAS => RAIN_ICE_DESCR%XTRANS_MP_GAMMAS
  XLBDAS_MIN => RAIN_ICE_DESCR%XLBDAS_MIN
END SUBROUTINE
!
SUBROUTINE RAIN_ICE_DESCR_ALLOCATE(KRR)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: KRR
  ALLOCATE(RAIN_ICE_DESCR%XRTMIN(KRR))
  XRTMIN=>RAIN_ICE_DESCR%XRTMIN
  XLBC=>RAIN_ICE_DESCR%XLBC
END SUBROUTINE RAIN_ICE_DESCR_ALLOCATE
!
SUBROUTINE RAIN_ICE_DESCR_DEALLOCATE()
  IMPLICIT NONE
  XRTMIN=>NULL()
  DEALLOCATE(RAIN_ICE_DESCR%XRTMIN)
END SUBROUTINE RAIN_ICE_DESCR_DEALLOCATE
!
END MODULE MODD_RAIN_ICE_DESCR
