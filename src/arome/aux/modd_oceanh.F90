!MNH_LIC Copyright 2021-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      MODULE MODD_OCEANH
!     #################
!
!!****  *MODD_OCEAN* - declaration of variables used in ocean version
!!
!!    PURPOSE
!!    -------
!       Declarative module for the variables
!!      at interface for OCEAN LES MESONH version including auto-coupling O-A LES
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    AUTHOR
!!    ------
!!    JL Redelsperger LOPS
!!   
!!    MODIFICATIONS
!!    -------------
!!      Original 03/2021
!
!*       0.   DECLARATIONS
!             ------------
!
!USE MODD_TYPE_DATE
!
IMPLICIT NONE
!
!*            fields for Sea Sfc FORCINGs
!             ------------------
!
INTEGER,          SAVE                  :: NFRCLT     ! number of sea surface forcings PLUS 1
INTEGER,          SAVE                  :: NINFRT     ! Interval in second between forcings
!TYPE (DATE_TIME), SAVE, DIMENSION(:), ALLOCATABLE :: TFRCLT ! date/time of sea surface forcings
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE ::  XSSUFL,XSSVFL,XSSTFL,XSSOLA ! Time evol Flux U V T Solar_Rad at sea surface
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XSSUFL_XY,XSSVFL_XY,XSSTFL_XY! XY flux shape
REAL, SAVE, DIMENSION(:), ALLOCATABLE :: XSSUFL_T,XSSVFL_T,XSSTFL_T,XSSOLA_T ! given time forcing fluxes
!
END MODULE MODD_OCEANH
