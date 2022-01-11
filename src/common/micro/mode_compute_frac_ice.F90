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
!
! 
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
INTERFACE COMPUTE_FRAC_ICE 
  MODULE PROCEDURE COMPUTE_FRAC_ICE1D, COMPUTE_FRAC_ICE2D, COMPUTE_FRAC_ICE3D
END INTERFACE COMPUTE_FRAC_ICE

CONTAINS

!!       ==========
!!       3D version
!!       ==========
!
SUBROUTINE COMPUTE_FRAC_ICE3D(KIT, KJT, KKT, HFRAC_ICE, PFRAC_ICE, PT)
IMPLICIT NONE
INTEGER,                              INTENT(IN)    :: KIT, KJT, KKT
CHARACTER(LEN=1),                     INTENT(IN)    :: HFRAC_ICE ! scheme to use
REAL, TARGET, DIMENSION(KIT,KJT,KKT), INTENT(IN)    :: PT        ! Temperature
REAL, TARGET, DIMENSION(KIT,KJT,KKT), INTENT(INOUT) :: PFRAC_ICE ! Ice fraction (1 for ice only, 0 for liquid only)
REAL, POINTER, DIMENSION(:) :: ZT, ZFRAC_ICE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE3D',0,ZHOOK_HANDLE)
!
ZT(1:SIZE(PT))=>PT
ZFRAC_ICE(1:SIZE(PFRAC_ICE))=>PFRAC_ICE
CALL COMPUTE_FRAC_ICE1D(SIZE(PT), HFRAC_ICE,ZFRAC_ICE,ZT)
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE3D',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_FRAC_ICE3D
!
!
!!       ==========
!!       2D version
!!       ==========
!
SUBROUTINE COMPUTE_FRAC_ICE2D(KIT, KJT, HFRAC_ICE, PFRAC_ICE, PT)
IMPLICIT NONE
INTEGER,                          INTENT(IN)    :: KIT, KJT
CHARACTER(LEN=1),                 INTENT(IN)    :: HFRAC_ICE ! scheme to use
REAL, TARGET, DIMENSION(KIT,KJT), INTENT(IN)    :: PT        ! Temperature
REAL, TARGET, DIMENSION(KIT,KJT), INTENT(INOUT) :: PFRAC_ICE ! Ice fraction (1 for ice only, 0 for liquid only)
REAL, POINTER, DIMENSION(:) :: ZT, ZFRAC_ICE
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE2D',0,ZHOOK_HANDLE)
!
ZT(1:SIZE(PT))=>PT
ZFRAC_ICE(1:SIZE(PFRAC_ICE))=>PFRAC_ICE
CALL COMPUTE_FRAC_ICE1D(SIZE(PT), HFRAC_ICE,ZFRAC_ICE,ZT)
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE2D',1,ZHOOK_HANDLE)
END SUBROUTINE COMPUTE_FRAC_ICE2D
!
!
!!       ==========
!!       1D version
!!       ==========
!
SUBROUTINE COMPUTE_FRAC_ICE1D(KIT, HFRAC_ICE, PFRAC_ICE, PT)
!! --------------------------------------------------------------------------
!       0. DECLARATIONS
!          ------------
!
USE MODD_NEB, ONLY : XTMINMIX, XTMAXMIX
USE MODD_CST, ONLY : XTT
USE MODE_MSG, ONLY : PRINT_MSG, NVERB_FATAL
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
INTEGER,                      INTENT(IN)    :: KIT
CHARACTER(LEN=1),             INTENT(IN)    :: HFRAC_ICE ! scheme to use
REAL, TARGET, DIMENSION(KIT), INTENT(IN)    :: PT        ! Temperature
REAL, TARGET, DIMENSION(KIT), INTENT(INOUT) :: PFRAC_ICE ! Ice fraction (1 for ice only, 0 for liquid only)
!
!               0.2  declaration of local variables
! 
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE1D',0,ZHOOK_HANDLE)
!------------------------------------------------------------------------
!                1. Compute FRAC_ICE
!
SELECT CASE(HFRAC_ICE)
  CASE ('T') !using Temperature
    PFRAC_ICE(:) = MAX( 0., MIN(1., (( XTMAXMIX - PT(:) ) / ( XTMAXMIX - XTMINMIX )) ) ) ! freezing interval
  CASE ('O') !using Temperature with old formulae
    PFRAC_ICE(:) = MAX( 0., MIN(1., (( XTT - PT(:) ) / 40.) ) ) ! freezing interval
  CASE ('N') !No ice
    PFRAC_ICE(:) = 0.
  CASE ('S') !Same as previous
    ! (almost) nothing to do
    PFRAC_ICE(:) = MAX( 0., MIN(1., PFRAC_ICE(:) ) )
  CASE DEFAULT
    CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'COMPUTE_FRAC_ICE', 'invalid option for HFRAC_ICE='//HFRAC_ICE)
END SELECT
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_FRAC_ICE1D',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_FRAC_ICE1D
!
END MODULE MODE_COMPUTE_FRAC_ICE
