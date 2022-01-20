!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############
      MODULE MODD_REF
!     ###############
!
!!****  *MODD_REF* - declaration of reference state  profile
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the vertical
!     profile of  the reference state, used for the anelastic 
!     approximation. 
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_REF)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   07/06/94   
!!                    07/13 (C.Lac) Add LBOUSS
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!          
REAL,SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: XRHODREFZ ! rhod(z) for reference
                                             ! state without orography
REAL,SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: XTHVREFZ  ! Thetav(z) for reference
                                             ! state without orography    
REAL,SAVE                            :: XEXNTOP   ! Exner function at model top 
!
! For coupled A-O case
REAL,SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: XRHODREFZO! rhod(z) for ocean ref state in coupled mode
REAL,SAVE, DIMENSION(:), ALLOCATABLE, TARGET :: XTHVREFZO !Thetav(z) for ocean ref state in coupled mode
REAL,SAVE                            :: XEXNTOPO   ! Exner function at ocean  model top in coupled mode
!
LOGICAL, SAVE                        :: LBOUSS    ! Boussinesq approximation
LOGICAL, SAVE   ::LCOUPLES ! AUTOCOUPLED ATMS-OCEAN LES VERSION
! 
END MODULE MODD_REF
