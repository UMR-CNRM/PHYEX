!MNH_LIC Copyright 2000-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODD_RAIN_C2R2_DESCR
!     ###########################
!
!!****  *MODD_RAIN_C2R2_DESCR* - declaration of the microphysical descriptive
!!                               constants for use in the warm scheme.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the microphysical
!     constants. This includes the descriptive parameters for the raindrop 
!     and the parameters relevant of the dimensional distributions.
!
!         m(D)    = XAx * D**XBx      : Mass-MaxDim relationship
!         v(D)    = XCx * D**XDx      : Fallspeed-MaxDim relationship
!         XF0x, XF1x, XF2x            : Ventilation factors
!
!              and
!
!         XALPHAx, XNUx                        : Generalized GAMMA law 
!         XLBx                                 : Slope parameter of the 
!                                                distribution law
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_RAIN_C2R2_DESCR)
!!          
!!    AUTHOR
!!    ------
!!	J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/11/2000
!!       J.-P. Pinty   29/11/02 add cloud doplet fall speed parameters
!!
!-------------------------------------------------------------------------------
USE MODD_PARAMETERS, ONLY: JPSVNAMELGTMAX
IMPLICIT NONE
!
!*       0.   DECLARATIONS
!             ------------
!
REAL,SAVE :: XCEXVT                    ! air density fall speed correction
!
REAL,SAVE :: XAR,XBR,XCR,XDR,XF0R,XF1R,     & ! Raindrop       charact.
	     XAC,XBC,XCC,XDC,XF0C,XF2C,XC1C   ! Cloud droplet  charact.
!
REAL,DIMENSION(:),SAVE,ALLOCATABLE :: XRTMIN
                                       ! Min values of the mixing ratios
REAL,DIMENSION(:),SAVE,ALLOCATABLE :: XCTMIN
                                       ! Min values of the drop concentrations
REAL,SAVE ::  XLBC, XLBEXC,          & ! shape parameters of the cloud droplets
	      XLBR, XLBEXR             ! shape parameters of the raindrops
!
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(4),PARAMETER &
                                   :: C2R2NAMES=(/'CCCN  ','CCLOUD','CRAIN ','SUPSAT'/)
                                       ! basenames of the SV articles stored
                                       ! in the binary files
!
END MODULE MODD_RAIN_C2R2_DESCR
!
!
