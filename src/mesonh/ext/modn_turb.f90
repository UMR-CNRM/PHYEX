!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modn 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ###################
      MODULE MODN_TURB
!     ###################
!
!!****  *MODN_TURB* - declaration of namelist NAM_TURB
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_TURB
!     which concern the parameters of the turbulence scheme for all models
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	    V. Masson                  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    November 2005
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_CTURB
!
IMPLICIT NONE
!
NAMELIST/NAM_TURB/XPHI_LIM, XSBL_O_BL, XFTOP_O_FSURF
!
END MODULE MODN_TURB
