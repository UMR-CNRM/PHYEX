!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
      MODULE MODD_NEB
!     #############################
!
!!****  *MODD_NEB* - Declaration of nebulosity constants
!!
!!    PURPOSE
!!    -------
!!      The purpose of this declarative module is to declare some
!!      constants for nebulosity calculation
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!
!!
!!    AUTHOR
!!    ------
!!       S. Riette (Meteo France)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    24 Aug 2011
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
REAL,SAVE          :: XTMINMIX   ! minimum temperature of mixed phase
REAL,SAVE          :: XTMAXMIX   ! maximum temperature of mixed phase
!
!
END MODULE MODD_NEB
