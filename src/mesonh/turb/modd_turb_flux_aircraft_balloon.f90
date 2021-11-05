!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     ######################################
      MODULE MODD_TURB_FLUX_AIRCRAFT_BALLOON
!     ######################################
!
!!****  *MODD_CVERT* - Declares work arrays for vertical cross-sections
!!
!!    PURPOSE
!!    -------
!       For vertical cross-sections only, this declarative module declares 
!     the arrays containing the sea-level altitudes and the model topography 
!     of the oblique cross-section points.     
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!
!!     Book2 of the TRACE volume of the Meso-NH user manual
!!     (MODD_CVERT) 
!!       
!!    AUTHOR
!!    ------
!!      P.Lacarrere
!!
!!    MODIFICATIONS
!!    -------------
!!      Original       18/09/06                      
!!        
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!
IMPLICIT NONE
!
REAL,DIMENSION(:,:,:)  ,ALLOCATABLE,SAVE   :: XTHW_FLUX !sensible flux 
REAL,DIMENSION(:,:,:)  ,ALLOCATABLE,SAVE   :: XRCW_FLUX !Latent flux
REAL,DIMENSION(:,:,:,:),ALLOCATABLE,SAVE   :: XSVW_FLUX !turb scalar flux
!
END MODULE MODD_TURB_FLUX_AIRCRAFT_BALLOON
