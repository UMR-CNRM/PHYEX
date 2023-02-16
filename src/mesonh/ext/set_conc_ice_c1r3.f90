!MNH_LIC Copyright 2001-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #############################
       MODULE MODI_SET_CONC_ICE_C1R3
!      #############################
!
INTERFACE
!
      SUBROUTINE SET_CONC_ICE_C1R3 (PRHODREF,PRT,PSVT)
!
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:,:),  INTENT(INOUT) :: PRT    ! microphysical mixing ratios
!
REAL,  DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT   ! microphys. concentrations
!
!
END SUBROUTINE SET_CONC_ICE_C1R3
!
END INTERFACE
!
END MODULE MODI_SET_CONC_ICE_C1R3
!
!     ##########################################################
      SUBROUTINE SET_CONC_ICE_C1R3 (PRHODREF,PRT,PSVT)
!     ##########################################################
!
!!****  *SET_CONC_ICE_C1R3 * - initialize the ice crystal
!!                   concentration for a RESTArt simulation of the C1R3 scheme
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the pristine ice crystal
!!    concentrations when the cloud ice mixing ratios are only available.
!!      This routine is used to initialize the small ice crystal concentrations
!!    using the r_i of a previous ICE3 run but also to compute the LB tendencies
!!    in ONE_WAY$n in case of grid-nesting when the optional argument PTIME is
!!    set (a C3R5 run embedded in a ICE3 run).
!!
!!**  METHOD
!!    ------
!!      The method uses the contact nucleation formulation of Meyers as a rough
!!    estimate (a function of the temperature). A limiting value of XCONCI_MAX
!!    is also assumed in the case of very cold temperatures
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_ICE_C1R3_DESCR, ONLY : XRTMIN, XCTMIN
!!      Module MODD_ICE_C1R3_PARAM, ONLY : XCONCI_INI
!!      Module MODD_CONF,           ONLY : NVERB
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine SET_CONC_ICE_C1R3 )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    15/04/01
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,            ONLY : XRHOLI
USE MODD_CONF,           ONLY : NVERB
USE MODD_ICE_C1R3_DESCR, ONLY : XRTMIN, XCTMIN
USE MODD_ICE_C1R3_PARAM, ONLY : XCONCI_MAX, XNUC_CON, XEXTT_CON, XEX_CON
USE MODD_LUNIT_n,        ONLY : TLUOUT
USE MODD_RAIN_ICE_DESCR, ONLY : XAI, XBI
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF   ! Reference density
!
REAL, DIMENSION(:,:,:,:),  INTENT(INOUT) :: PRT    ! microphysical mixing ratios
!
REAL,  DIMENSION(:,:,:,:), INTENT(INOUT) :: PSVT   ! microphys. concentrations
!
!*       0.2   Declarations of local variables :
!
INTEGER    :: IRESP   ! Return code of FM routines
INTEGER    :: ILUOUT  ! Logical unit number of output-listing
!
!
!-------------------------------------------------------------------------------
!*       1.    RETRIEVE LOGICAL UNIT NUMBER
!              ----------------------------
!
ILUOUT = TLUOUT%NLU
!
!*       2.    INITIALIZATION
!              --------------
!
! Assume the ice crystal concentration according to the
! contact nucleation formulation of Meyers et al. (1992)
!
WHERE ( PRT(:,:,:,4) > XRTMIN(4) )
  PSVT(:,:,:,4) = MIN( PRHODREF(:,:,:) /                                     &
                             ( XRHOLI * XAI*(10.E-06)**XBI * PRT(:,:,:,4) ), &
                      XCONCI_MAX )
  PSVT(:,:,:,5) = 0.0
END WHERE
WHERE ( PRT(:,:,:,4) <= XRTMIN(4) )
  PRT(:,:,:,4)  = 0.0
  PSVT(:,:,:,4) = 0.0
  PSVT(:,:,:,5) = 0.0
END WHERE
IF( NVERB >= 5 ) THEN
  WRITE (UNIT=ILUOUT,FMT=*) "!INI_MODEL$n: The cloud ice concentration has "
  WRITE (UNIT=ILUOUT,FMT=*) "been roughly initialised to a value of 1 per liter"
END IF
!
END SUBROUTINE SET_CONC_ICE_C1R3
