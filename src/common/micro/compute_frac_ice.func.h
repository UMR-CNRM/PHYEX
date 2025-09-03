!MNH_LIC Copyright 2006-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
      ELEMENTAL SUBROUTINE COMPUTE_FRAC_ICE(CST,HFRAC_ICE,NEBN,PFRAC_ICE,PT)

! ******* TO BE INCLUDED IN THE *CONTAINS* OF A SUBROUTINE, IN ORDER TO EASE AUTOMATIC INLINING ******
! => Don't use drHook !!!
!
!!****  *COMPUTE_FRAC_ICE* - computes ice fraction
!
!!    AUTHOR
!!    ------
!!      Julien PERGAUD      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         13/03/06
!!      S. Riette        April 2011 optimisation
!!      S. Riette        08/2016 add option O
!!      R. El Khatib     12-Aug-2021 written as a include file
!!      S. Riette        Sept 2025, remove KERR argument
!
!! --------------------------------------------------------------------------
USE MODD_NEB_n, ONLY : NEB_t
USE MODD_CST, ONLY : CST_t
!
IMPLICIT NONE
!
TYPE(CST_t),      INTENT(IN)    :: CST
CHARACTER(LEN=1), INTENT(IN)    :: HFRAC_ICE       ! scheme to use
TYPE(NEB_t),      INTENT(IN)    :: NEBN
REAL,             INTENT(IN)    :: PT              ! temperature
REAL,             INTENT(INOUT) :: PFRAC_ICE       ! Ice fraction (1 for ice only, 0 for liquid only)
!
!------------------------------------------------------------------------

!                1. Compute FRAC_ICE
!
SELECT CASE(HFRAC_ICE)
  CASE ('T') !using Temperature
    PFRAC_ICE = MAX( 0., MIN(1., (( NEBN%XTMAXMIX - PT ) / ( NEBN%XTMAXMIX - NEBN%XTMINMIX )) ) ) ! freezing interval
  CASE ('O') !using Temperature with old formulae
    PFRAC_ICE = MAX( 0., MIN(1., (( CST%XTT - PT ) / 40.) ) ) ! freezing interval
  CASE ('N') !No ice
    PFRAC_ICE = 0.
  CASE ('S') !Same as previous
    ! (almost) nothing to do
    PFRAC_ICE = MAX( 0., MIN(1., PFRAC_ICE ) )
END SELECT

END SUBROUTINE COMPUTE_FRAC_ICE
