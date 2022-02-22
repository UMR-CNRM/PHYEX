!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_GRID_n
!     ##################
!
!!****  *MODD_GRID$n* - declaration of grid variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare  the variables
!     describing the grid. 
!    
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_GRIDn)
!!      Technical Specifications Report of the Meso-NH (chapters 2 and 3)
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/05/94                      
!!      J. Stein    15/11/95  add the slope angle
!!      V. Ducrocq   13/08/98  // : add XLATOR_ll and XLONOR_ll       
!!      V. Masson   nov 2004  supress XLATOR,XLONOR,XLATOR_ll,XLONOR_ll
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

REAL, DIMENSION(:,:),  POINTER :: XLON=>NULL(),XLAT=>NULL() ! Longitude and latitude  
REAL, DIMENSION(:),    POINTER :: XXHAT=>NULL()             ! Position x in the conformal or cartesian plane
REAL, DIMENSION(:),    POINTER :: XYHAT=>NULL()             ! Position y in the conformal or cartesian plane
REAL, DIMENSION(:),    POINTER :: XDXHAT=>NULL()            ! horizontal stretching in x
REAL, DIMENSION(:),    POINTER :: XDYHAT=>NULL()            ! horizontal stretching in y
REAL, DIMENSION(:,:),  POINTER :: XMAP=>NULL()              ! Map factor 
REAL, DIMENSION(:,:),  POINTER :: XZS=>NULL()               ! orography
REAL, DIMENSION(:,:,:),POINTER :: XZZ=>NULL()               ! height z 
REAL,                  POINTER :: XZTOP=>NULL()             ! model top (m)
REAL, DIMENSION(:),    POINTER :: XZHAT=>NULL()             ! height level without orography
REAL, DIMENSION(:,:),  POINTER :: XDIRCOSXW=>NULL(),XDIRCOSYW=>NULL(),XDIRCOSZW=>NULL() ! director cosinus of the normal
                                                                                        ! to the ground surface
REAL, DIMENSION(:,:),  POINTER  :: XCOSSLOPE=>NULL()         ! cosinus of the angle between i and the slope vector
REAL, DIMENSION(:,:),  POINTER  :: XSINSLOPE=>NULL()         ! sinus   of the angle between i and the slope vector
! Quantities for SLEVE vertical coordinate
LOGICAL,               POINTER  :: LSLEVE=>NULL()            ! Logical for SLEVE coordinate
REAL,                  POINTER  :: XLEN1=>NULL()             ! Decay scale for smooth topography
REAL,                  POINTER  :: XLEN2=>NULL()             ! Decay scale for small-scale topography deviation
REAL, DIMENSION(:,:),  POINTER  :: XZSMT=>NULL()             ! smooth orography for SLEVE coordinate

END MODULE MODD_GRID_n
