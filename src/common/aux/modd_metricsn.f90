!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/06/27 14:20:29
!-----------------------------------------------------------------
!     #####################
      MODULE MODD_METRICS_n
!     #####################
!
!!****  *MODD_METRICS$n* - metric coefficients
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the 
!     metric coefficients. 
!    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!	P. Jabouille   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/04/99                      
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE METRICS_t
  REAL, DIMENSION(:,:,:), POINTER :: XDXX=>NULL(),XDZX=>NULL(), &
                                  XDYY=>NULL(),XDZY=>NULL(),XDZZ=>NULL()
                                            !metric coefficients
END TYPE METRICS_t

TYPE(METRICS_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: METRICS_MODEL

REAL, DIMENSION(:,:,:), POINTER :: XDXX=>NULL(),XDZX=>NULL(), &
                                  XDYY=>NULL(),XDZY=>NULL(),XDZZ=>NULL()

CONTAINS

SUBROUTINE METRICS_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
METRICS_MODEL(KFROM)%XDXX=>XDXX
METRICS_MODEL(KFROM)%XDZX=>XDZX
METRICS_MODEL(KFROM)%XDYY=>XDYY
METRICS_MODEL(KFROM)%XDZY=>XDZY
METRICS_MODEL(KFROM)%XDZZ=>XDZZ
!
! Current model is set to model KTO
XDXX=>METRICS_MODEL(KTO)%XDXX
XDZX=>METRICS_MODEL(KTO)%XDZX
XDYY=>METRICS_MODEL(KTO)%XDYY
XDZY=>METRICS_MODEL(KTO)%XDZY
XDZZ=>METRICS_MODEL(KTO)%XDZZ

END SUBROUTINE METRICS_GOTO_MODEL

END MODULE MODD_METRICS_n
