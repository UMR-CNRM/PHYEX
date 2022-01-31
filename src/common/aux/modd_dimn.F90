!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_DIM_n
!     ##################
!
!!****  *MODD_DIM$n* - declaration of dimensions
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the dimensions 
!     of the data arrays.   
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_DIMn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!          
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94     
!!      Modifications 13/08/98 (V. Ducrocq) // NIINF .. NJSUP are no more used in the init part                
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE DIM_t
  INTEGER :: NIMAX,NJMAX,NKMAX  !  Dimensions respectively  in x , 
                              ! y ,  z directions of the physical sub-domain.
  INTEGER :: NIMAX_ll,NJMAX_ll  !  Dimensions respectively  in x and y
                                   ! directions of the physical domain
  INTEGER :: NIINF, NISUP       !  Lower bound and upper bound of the arrays 
                                   ! in x direction 
  INTEGER :: NJINF, NJSUP       !  Lower bound and upper bound of the arrays 
                                   ! in y direction
!
END TYPE DIM_t

TYPE(DIM_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: DIM_MODEL

INTEGER, POINTER :: NIMAX=>NULL(),NJMAX=>NULL(),NKMAX=>NULL()
INTEGER, POINTER :: NIMAX_ll=>NULL(),NJMAX_ll=>NULL()
INTEGER, POINTER :: NIINF=>NULL(), NISUP=>NULL()
INTEGER, POINTER :: NJINF=>NULL(), NJSUP=>NULL()

CONTAINS

SUBROUTINE DIM_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
!
! Current model is set to model KTO
NIMAX=>DIM_MODEL(KTO)%NIMAX
NJMAX=>DIM_MODEL(KTO)%NJMAX
NKMAX=>DIM_MODEL(KTO)%NKMAX
NIMAX_ll=>DIM_MODEL(KTO)%NIMAX_ll
NJMAX_ll=>DIM_MODEL(KTO)%NJMAX_ll
NIINF=>DIM_MODEL(KTO)%NIINF
NISUP=>DIM_MODEL(KTO)%NISUP
NJINF=>DIM_MODEL(KTO)%NJINF
NJSUP=>DIM_MODEL(KTO)%NJSUP

END SUBROUTINE DIM_GOTO_MODEL

END MODULE MODD_DIM_n
