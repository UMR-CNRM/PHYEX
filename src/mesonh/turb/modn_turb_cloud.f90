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
!     ##################
      MODULE MODN_TURB_CLOUD
!     ##################
!
!!****  *MODN_TURB_CLOUD* - declaration of namelist NAM_TURB_CLOUD
!!
!!    PURPOSE
!!    -------
!       The purpose of this module is to specify  the namelist NAM_TURB_CLOUD
!     which concern the parameters of the cloud mixing length for a given model.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_TURB_CLOUD
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	M. Tomasini    *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    September, 2004
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_TURB_CLOUD
!
IMPLICIT NONE
!
NAMELIST/NAM_TURB_CLOUD/NMODEL_CLOUD, CTURBLEN_CLOUD, &
                        XCOEF_AMPL_SAT, XCEI_MIN, XCEI_MAX
!
END MODULE MODN_TURB_CLOUD
