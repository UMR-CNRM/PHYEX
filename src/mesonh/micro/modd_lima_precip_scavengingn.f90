!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################################
      MODULE MODD_LIMA_PRECIP_SCAVENGING_n
!     ####################################
!
!!****  *MODD_PRECIP_SCAVENGING$n* - declaration of scavenged aerosols 
!!                                   precipitating fields
!!
!!    PURPOSE
!!    -------
!       Stores the INstantaneous and ACcumulated PRecipitating fields of 
!!      scavenged aerosol by rain
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_PARAMETERS, ONLY: JPMODELMAX
!
IMPLICIT NONE
!
REAL, DIMENSION(:,:), POINTER :: XINPAP=>NULL(), XACPAP=>NULL()
                                         ! Instant and cumul of ground
                                         ! precipitation fields of Scavenged
                                         ! Aerosol Particles

CONTAINS

SUBROUTINE LIMA_PRECIP_SCAVENGING_GOTO_MODEL(KFROM, KTO)
  INTEGER, INTENT(IN) :: KFROM, KTO
END SUBROUTINE LIMA_PRECIP_SCAVENGING_GOTO_MODEL
!
!
END MODULE MODD_LIMA_PRECIP_SCAVENGING_n
