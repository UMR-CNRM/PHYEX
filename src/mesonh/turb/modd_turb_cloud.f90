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
      MODULE MODD_TURB_CLOUD
!     ##################
!
!!****  *MODD_TURB_CLOUD* - declaration of parameters for cloud mixing length
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the
!     variables that may be set by namelist for the cloud mixing length
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
!!	M. Tomasini     *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    September, 2004
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
! 
INTEGER,SAVE            :: NMODEL_CLOUD    ! model number where the modification                                ! of the mixing length in the clouds is computed
CHARACTER (LEN=4),SAVE  :: CTURBLEN_CLOUD  ! type of length in the clouds
                                     ! 'DEAR' Deardorff mixing length
                                     ! 'BL89' Bougeault and Lacarrere scheme
                                     ! 'DELT' length = ( volum) ** 1/3
REAL,SAVE               :: XCOEF_AMPL_SAT  ! saturation of the amplification coefficient
REAL,SAVE               :: XCEI_MIN  ! minimum threshold for the instability index CEI
                                     !(beginning of the amplification)
REAL,SAVE               :: XCEI_MAX  ! maximum threshold for the instability index CEI
                                     !(beginning of the saturation of the amplification)
REAL,SAVE,DIMENSION(:,:,:), ALLOCATABLE  :: XCEI ! Cloud Entrainment instability
                                                 ! index to emphasize localy 
                                                 ! turbulent fluxes
!
END MODULE MODD_TURB_CLOUD
