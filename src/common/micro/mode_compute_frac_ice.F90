!MNH_LIC Copyright 2006-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODE_COMPUTE_FRAC_ICE
!    ############################
!
!!****  *COMPUTE_FRAC_ICE* - computes ice fraction
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Julien PERGAUD      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         13/03/06
!!      S. Riette        April 2011 optimisation
!!      S. Riette        08/2016 add option O
!!      P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!!      R. El Khatib 24-Aug-2021 Optimization by cache re-use + assume data is contiguous
!!      S. Riette Jan-2022 Merge the 3 procedures in one module + array shape declaration
!       R. El Khatib / S. Riette Jan-2022 written as an elemental subroutine
!
! 
!


!****************** Don't use drHook !!!



IMPLICIT NONE
CONTAINS

ELEMENTAL SUBROUTINE COMPUTE_FRAC_ICE(HFRAC_ICE, PFRAC_ICE, PT, KERR)
!! --------------------------------------------------------------------------
!       0. DECLARATIONS
!          ------------
!
USE MODD_NEB, ONLY : XTMINMIX, XTMAXMIX
USE MODD_CST, ONLY : XTT
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER(LEN=1),  INTENT(IN)    :: HFRAC_ICE ! scheme to use
REAL,              INTENT(IN)    :: PT        ! Temperature
REAL,              INTENT(INOUT) :: PFRAC_ICE ! Ice fraction (1 for ice only, 0 for liquid only)
INTEGER,           INTENT(OUT)   :: KERR      ! Error code in return
!
!               0.2  declaration of local variables
! 
!------------------------------------------------------------------------
!                1. Compute FRAC_ICE
!
KERR=0
SELECT CASE(HFRAC_ICE)
  CASE ('T') !using Temperature
    PFRAC_ICE = MAX( 0., MIN(1., (( XTMAXMIX - PT ) / ( XTMAXMIX - XTMINMIX )) ) ) ! freezing interval
  CASE ('O') !using Temperature with old formulae
    PFRAC_ICE = MAX( 0., MIN(1., (( XTT - PT ) / 40.) ) ) ! freezing interval
  CASE ('N') !No ice
    PFRAC_ICE = 0.
  CASE ('S') !Same as previous
    ! (almost) nothing to do
    PFRAC_ICE = MAX( 0., MIN(1., PFRAC_ICE ) )
  CASE DEFAULT
    KERR=1
END SELECT
!
END SUBROUTINE COMPUTE_FRAC_ICE
!
END MODULE MODE_COMPUTE_FRAC_ICE
