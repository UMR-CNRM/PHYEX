!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ##########################
      MODULE MODD_RAIN_ICE_DESCR_n
!     ##########################
!> @file
!!****  *MODD_RAIN_ICE_DESCR_n* - declaration of the microphysical descriptive
!!                                constants for use in the warm and cold schemes.
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
USE MODD_PARAMETERS, ONLY: JPMODELMAX
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
TYPE(RAIN_ICE_DESCR_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: RAIN_ICE_DESCR_MODEL
TYPE(RAIN_ICE_DESCR_t), POINTER, SAVE :: RAIN_ICE_DESCRN => NULL()
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
SUBROUTINE RAIN_ICE_DESCR_GOTO_MODEL(KFROM, KTO)
!! This subroutine associate all the pointers to the right component of
!! the right strucuture. A value can be accessed through the structure RAIN_ICE_DESCRN
!! or through the strucuture RAIN_ICE_DESCR_MODEL(KTO) or directly through these pointers.
IMPLICIT NONE
INTEGER, INTENT(IN) :: KFROM, KTO
!
IF(.NOT. ASSOCIATED(RAIN_ICE_DESCRN, RAIN_ICE_DESCR_MODEL(KTO))) THEN    
  !
  RAIN_ICE_DESCRN => RAIN_ICE_DESCR_MODEL(KTO) 
  !
  XCEXVT => RAIN_ICE_DESCRN%XCEXVT
  XAC => RAIN_ICE_DESCRN%XAC
  XBC => RAIN_ICE_DESCRN%XBC
  XCC => RAIN_ICE_DESCRN%XCC
  XDC => RAIN_ICE_DESCRN%XDC
  XAR => RAIN_ICE_DESCRN%XAR
  XBR => RAIN_ICE_DESCRN%XBR
  XCR => RAIN_ICE_DESCRN%XCR
  XDR => RAIN_ICE_DESCRN%XDR
  XCCR => RAIN_ICE_DESCRN%XCCR
  XF0R => RAIN_ICE_DESCRN%XF0R
  XF1R => RAIN_ICE_DESCRN%XF1R
  XC1R => RAIN_ICE_DESCRN%XC1R
  XAI => RAIN_ICE_DESCRN%XAI
  XBI => RAIN_ICE_DESCRN%XBI
  XC_I => RAIN_ICE_DESCRN%XC_I
  XDI => RAIN_ICE_DESCRN%XDI
  XF0I => RAIN_ICE_DESCRN%XF0I
  XF2I => RAIN_ICE_DESCRN%XF2I
  XC1I => RAIN_ICE_DESCRN%XC1I
  XAS => RAIN_ICE_DESCRN%XAS
  XBS => RAIN_ICE_DESCRN%XBS
  XCS => RAIN_ICE_DESCRN%XCS
  XDS => RAIN_ICE_DESCRN%XDS
  XCCS => RAIN_ICE_DESCRN%XCCS
  XCXS => RAIN_ICE_DESCRN%XCXS
  XF0S => RAIN_ICE_DESCRN%XF0S
  XF1S => RAIN_ICE_DESCRN%XF1S
  XC1S => RAIN_ICE_DESCRN%XC1S
  XAG => RAIN_ICE_DESCRN%XAG
  XBG => RAIN_ICE_DESCRN%XBG
  XCG => RAIN_ICE_DESCRN%XCG
  XDG => RAIN_ICE_DESCRN%XDG
  XCCG => RAIN_ICE_DESCRN%XCCG
  XCXG => RAIN_ICE_DESCRN%XCXG
  XF0G => RAIN_ICE_DESCRN%XF0G
  XF1G => RAIN_ICE_DESCRN%XF1G
  XC1G => RAIN_ICE_DESCRN%XC1G
  XAH => RAIN_ICE_DESCRN%XAH
  XBH => RAIN_ICE_DESCRN%XBH
  XCH => RAIN_ICE_DESCRN%XCH
  XDH => RAIN_ICE_DESCRN%XDH
  XCCH => RAIN_ICE_DESCRN%XCCH
  XCXH => RAIN_ICE_DESCRN%XCXH
  XF0H => RAIN_ICE_DESCRN%XF0H
  XF1H => RAIN_ICE_DESCRN%XF1H
  XC1H => RAIN_ICE_DESCRN%XC1H
  XALPHAC => RAIN_ICE_DESCRN%XALPHAC
  XNUC => RAIN_ICE_DESCRN%XNUC
  XALPHAC2 => RAIN_ICE_DESCRN%XALPHAC2
  XNUC2 => RAIN_ICE_DESCRN%XNUC2
  XLBEXC => RAIN_ICE_DESCRN%XLBEXC
  XALPHAR => RAIN_ICE_DESCRN%XALPHAR
  XNUR => RAIN_ICE_DESCRN%XNUR
  XLBEXR => RAIN_ICE_DESCRN%XLBEXR
  XLBR => RAIN_ICE_DESCRN%XLBR
  XALPHAI => RAIN_ICE_DESCRN%XALPHAI
  XNUI => RAIN_ICE_DESCRN%XNUI
  XLBEXI => RAIN_ICE_DESCRN%XLBEXI
  XLBI => RAIN_ICE_DESCRN%XLBI
  XALPHAS => RAIN_ICE_DESCRN%XALPHAS
  XNUS => RAIN_ICE_DESCRN%XNUS
  XLBEXS => RAIN_ICE_DESCRN%XLBEXS
  XLBS => RAIN_ICE_DESCRN%XLBS
  XALPHAG => RAIN_ICE_DESCRN%XALPHAG
  XNUG => RAIN_ICE_DESCRN%XNUG
  XLBEXG => RAIN_ICE_DESCRN%XLBEXG
  XLBG => RAIN_ICE_DESCRN%XLBG
  XALPHAH => RAIN_ICE_DESCRN%XALPHAH
  XNUH => RAIN_ICE_DESCRN%XNUH
  XLBEXH => RAIN_ICE_DESCRN%XLBEXH
  XLBH => RAIN_ICE_DESCRN%XLBH
  XLBDAR_MAX => RAIN_ICE_DESCRN%XLBDAR_MAX
  XLBDAS_MAX => RAIN_ICE_DESCRN%XLBDAS_MAX
  XLBDAG_MAX => RAIN_ICE_DESCRN%XLBDAG_MAX
  XCONC_SEA => RAIN_ICE_DESCRN%XCONC_SEA
  XCONC_LAND => RAIN_ICE_DESCRN%XCONC_LAND
  XCONC_URBAN => RAIN_ICE_DESCRN%XCONC_URBAN
  XNS => RAIN_ICE_DESCRN%XNS
  XFVELOS => RAIN_ICE_DESCRN%XFVELOS
  XTRANS_MP_GAMMAS => RAIN_ICE_DESCRN%XTRANS_MP_GAMMAS
  XLBDAS_MIN => RAIN_ICE_DESCRN%XLBDAS_MIN
ENDIF
END SUBROUTINE RAIN_ICE_DESCR_GOTO_MODEL
!
SUBROUTINE RAIN_ICE_DESCR_ALLOCATE(KRR)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: KRR
  ALLOCATE(RAIN_ICE_DESCRN%XRTMIN(KRR))
  XRTMIN=>RAIN_ICE_DESCRN%XRTMIN
  XLBC=>RAIN_ICE_DESCRN%XLBC
END SUBROUTINE RAIN_ICE_DESCR_ALLOCATE
!
SUBROUTINE RAIN_ICE_DESCR_DEALLOCATE()
  IMPLICIT NONE
  XRTMIN=>NULL()
  DEALLOCATE(RAIN_ICE_DESCRN%XRTMIN)
END SUBROUTINE RAIN_ICE_DESCR_DEALLOCATE
!
END MODULE MODD_RAIN_ICE_DESCR_n
