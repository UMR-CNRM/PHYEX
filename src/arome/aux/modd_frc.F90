!MNH_LIC Copyright 1996-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############
      MODULE MODD_FRC
!     ###############
!
!!***  *MODD_FRC -  Declarative module for the forcing fields
!!
!!    PURPOSE
!!    -------
!       This module contains NFRC 1D-arrays used by FORCING (geostrophic wind
!     components, large scale vertical wind, theta and humidity profiles when
!     the relaxation option is used,large scale theta and humidity gradients
!     and the translation speed of the domain of simulation.
!     The following control parameters are used by FORCING:
!     - LGEOST_UV_FRC and LGEOST_TH_FRC
!     - LTEND_THRV_FRC and LTEND_UV_FRC
!     - LVERT_MOTION_FRC  
!     - LRELAX_THRV_FRC, LRELAX_UV_FRC and LRELAX_UVMEAN_FRC using:
!         XRELAX_TIME_FRC, XRELAX_HEIGHT_FRC and CRELAX_HEIGHT_TYPE
!     - LTRANS
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_FRC)
!!      
!!
!!    AUTHOR
!!    ------
!!	    Marc Georgelin Labo d'aerologie
!!
!!    MODIFICATIONS
!!    -------------
!!      Original 29/07/96 
!!      29/07/96 (Pinty&Suhre) revised
!!      18/11/96 J.-P. Pinty   addition of the translation
!!      27/01/98 P. Bechtold   use tendency forcing
!!                             add SST and surface pressure forcing
!!      01/2004  V. Masson     surface externalization: removes SST forcing
!!                   09/2017 Q.Rodier add LTEND_UV_FRC
!!      03/2021 JL Redelsperger Parameters defining sfc forcing shape for idealized ocean case
!!      06/2021 F. Couvreux    add LRELAX_UVMEAN_FRC
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
!USE MODD_TYPE_DATE
!
IMPLICIT NONE
!
!*            fields for FORCING
!             ------------------
!
INTEGER,          SAVE                  :: NFRC     ! number of forcing profiles
!TYPE (DATE_TIME), SAVE, DIMENSION(:), ALLOCATABLE :: TDTFRC ! date of
                                                    !  each forcing profile
!
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XUFRC,   &! geostrophic wind 
					                       XVFRC,   &! components U and V
					                       XWFRC     ! large scale vertical wind
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XTHFRC,  &! large scale TH profile
					                       XRVFRC    ! large scale RV profile
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XGXTHFRC,&! large scale TH gradient
                                           XGYTHFRC  ! along the X and Y axis
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XTENDTHFRC,&! large scale TH tendency
                                           XTENDRVFRC  ! large scale RV tendency
REAL, SAVE                              :: XUTRANS, &! horizontal components of
                                           XVTRANS   !        a constant
                                                     !   Galilean TRANSlation
REAL, SAVE, DIMENSION(:), ALLOCATABLE   :: XPGROUNDFRC! surf. pressure 
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XTENDUFRC   ! large scale U tendency
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XTENDVFRC   ! large scale V tendency
!
!*            control parameters for FORCING
!             ------------------------------
!
LOGICAL, SAVE     :: LGEOST_UV_FRC      ! enables geostrophic wind term
LOGICAL, SAVE     :: LGEOST_TH_FRC      ! enables thermal wind advection
LOGICAL, SAVE     :: LTEND_THRV_FRC     ! enables tendency forcing
LOGICAL, SAVE     :: LTEND_UV_FRC       ! enables tendency forcing of the wind
LOGICAL, SAVE     :: LVERT_MOTION_FRC   ! enables prescribed a forced vertical
					                    ! transport for all prognostic variables
LOGICAL, SAVE     :: LRELAX_THRV_FRC    ! enables temp. and humidity relaxation
LOGICAL, SAVE     :: LRELAX_UV_FRC      ! enables  horizontal wind relaxation applied to the full wind field
LOGICAL, SAVE     :: LRELAX_UVMEAN_FRC  ! enables  horizontal wind relaxation applied to the horiz. avg. wind
!
REAL,    SAVE     :: XRELAX_TIME_FRC    ! e-folding time for relaxation 
REAL,    SAVE     :: XRELAX_HEIGHT_FRC  ! height below which relaxation
                                        ! is never applied
CHARACTER(len=4), SAVE :: CRELAX_HEIGHT_TYPE ! "THGR" relax. above maximal dTH/dz
					                    ! (but always above XRELAX_HEIGHT_FRC)
					                    ! "FIXE" relax. above XRELAX_HEIGHT_FRC
!
LOGICAL, SAVE     :: LTRANS             ! enables a Galilean translation of the
                                        !         domain of simulation
LOGICAL, SAVE     :: LPGROUND_FRC       ! enables surf. pressure forcing
!
LOGICAL, SAVE     :: LDEEPOC            ! activates sfc forcing for ideal ocean deep conv 
REAL,    SAVE     :: XCENTX_OC          ! center of sfc forc for ideal ocean
REAL,    SAVE     :: XRADX_OC           ! radius of sfc forc for ideal ocean
REAL,    SAVE     :: XCENTY_OC          ! center of sfc forc for ideal ocean
REAL,    SAVE     :: XRADY_OC           ! radius of sfc forc for ideal ocean
!
END MODULE MODD_FRC
