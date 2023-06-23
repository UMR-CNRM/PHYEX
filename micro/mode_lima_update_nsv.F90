!MNH_LIC Copyright 1995-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_LIMA_UPDATE_NSV
!> @file
!! This module contains code to update the NSV module variable relative to the LIMA scheme
!
IMPLICIT NONE
!
CONTAINS
!
SUBROUTINE LIMA_UPDATE_NSV(LDINIT, KMI, KSV, CDCLOUD, LDUPDATE)
!!*** *LIMA_UPDATE_NSV* - update modd_nsv values realtive to LIMA
!!
!!*   PURPOSE
!!    -------
!!    The modd_nsv values relative to the LIMA scheme are initialised (if LDINIT is .TRUE)
!!    according to the micromisics scheme used.
!!    If LDUPDATE is .TRUE., the scalar values of modd_nsv module receive the values
!!    assigned for the KMI model.
!!
!!*   METHOD
!!    ------
!!
!!*   EXTERNAL
!!    --------
!!
!!*   IMPLICIT ARGUMENTS
!!    ------------------
!!    MODD_NSV, MODD_PARAM_LIMA
!!
!!*   REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    S. Riette
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    April 2023
!-------------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!       ---------------
!
USE MODD_NSV, ONLY: NSV_LIMA_BEG_A, NSV_LIMA_END_A, NSV_LIMA_A, &
                  & NSV_LIMA_NC_A, NSV_LIMA_NR_A, NSV_LIMA_NI_A, NSV_LIMA_NS_A, NSV_LIMA_NG_A, NSV_LIMA_NH_A, &
                  & NSV_LIMA_CCN_FREE_A, NSV_LIMA_CCN_ACTI_A, NSV_LIMA_SCAVMASS_A, &
                  & NSV_LIMA_IFN_FREE_A, NSV_LIMA_IFN_NUCL_A, NSV_LIMA_IMM_NUCL_A, &
                  & NSV_LIMA_HOM_HAZE_A, NSV_LIMA_SPRO_A, &
                  & NSV_LIMA_BEG, NSV_LIMA_END, NSV_LIMA, &
                  & NSV_LIMA_NC, NSV_LIMA_NR, NSV_LIMA_NI, NSV_LIMA_NS, NSV_LIMA_NG, NSV_LIMA_NH, &
                  & NSV_LIMA_CCN_FREE, NSV_LIMA_CCN_ACTI, NSV_LIMA_SCAVMASS, &
                  & NSV_LIMA_IFN_FREE, NSV_LIMA_IFN_NUCL, NSV_LIMA_IMM_NUCL, &
                  & NSV_LIMA_HOM_HAZE, NSV_LIMA_SPRO
USE MODD_PARAM_LIMA,      ONLY: NMOD_CCN, LSCAV, LAERO_MASS, &
                                NMOD_IFN, NMOD_IMM, LHHONI, &
                                LSPRO,  &
                                NMOM_C, NMOM_R, NMOM_I, NMOM_S, NMOM_G, NMOM_H
!
!* 0.1. Declaration of arguments
!       ------------------------
!
IMPLICIT NONE
LOGICAL,          INTENT(IN)    :: LDINIT   !< .TRUE. to fill the different NSV_LIMA_*_A arrays
INTEGER,          INTENT(IN)    :: KMI      !< model number
INTEGER,          INTENT(INOUT) :: KSV      !< IN: Initial value to use when filling the NSV_LIMA_*_A arrays; 
                                            !! OUT: Final value after having filled the arrays
CHARACTER(LEN=4), INTENT(IN)    :: CDCLOUD  !< Cloud scheme
LOGICAL,          INTENT(IN)    :: LDUPDATE !< .TRUE. to goto model
!
!*      1. INITIALISATION
!       -----------------
!
IF(LDINIT) THEN
  IF(CDCLOUD=='LIMA') THEN
    KSV = KSV+1
    NSV_LIMA_BEG_A(KMI) = KSV
    ! Nc
    IF (NMOM_C.GE.2) THEN
      NSV_LIMA_NC_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Nr
    IF (NMOM_R.GE.2) THEN
      NSV_LIMA_NR_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! CCN
    IF (NMOD_CCN .GT. 0) THEN
      NSV_LIMA_CCN_FREE_A(KMI) = KSV
      KSV = KSV + NMOD_CCN
      NSV_LIMA_CCN_ACTI_A(KMI) = KSV
      KSV = KSV + NMOD_CCN
    END IF
    ! Scavenging
    IF (LSCAV .AND. LAERO_MASS) THEN
      NSV_LIMA_SCAVMASS_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Ni
    IF (NMOM_I.GE.2) THEN
      NSV_LIMA_NI_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Ns
    IF (NMOM_S.GE.2) THEN
      NSV_LIMA_NS_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Ng
    IF (NMOM_G.GE.2) THEN
      NSV_LIMA_NG_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Nh
    IF (NMOM_H.GE.2) THEN
      NSV_LIMA_NH_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! IFN
    IF (NMOD_IFN .GT. 0) THEN
      NSV_LIMA_IFN_FREE_A(KMI) = KSV
      KSV = KSV + NMOD_IFN
      NSV_LIMA_IFN_NUCL_A(KMI) = KSV
      KSV = KSV + NMOD_IFN
    END IF
    ! IMM
    IF (NMOD_IMM .GT. 0) THEN
      NSV_LIMA_IMM_NUCL_A(KMI) = KSV
      KSV = KSV + MAX(1,NMOD_IMM)
    END IF
    ! Homogeneous freezing of CCN
    IF (LHHONI) THEN
      NSV_LIMA_HOM_HAZE_A(KMI) = KSV
      KSV = KSV + 1
    END IF
    ! Supersaturation
    IF (LSPRO) THEN
      NSV_LIMA_SPRO_A(KMI) = KSV
      KSV = KSV + 1
    END IF
    !
    ! End and total variables
    !
    KSV = KSV - 1
    NSV_LIMA_END_A(KMI) = KSV
    NSV_LIMA_A(KMI) = NSV_LIMA_END_A(KMI) - NSV_LIMA_BEG_A(KMI) + 1
  ELSE
    NSV_LIMA_A(KMI)    = 0
    !
    ! force First index to be superior to last index
    ! in order to create a null section
    !
    NSV_LIMA_BEG_A(KMI) = 1
    NSV_LIMA_END_A(KMI) = 0
  ENDIF
ENDIF
!
!*      2. UPDATE
!       ---------
!
IF(LDUPDATE) THEN
  NSV_LIMA          = NSV_LIMA_A(KMI)
  NSV_LIMA_BEG      = NSV_LIMA_BEG_A(KMI)
  NSV_LIMA_END      = NSV_LIMA_END_A(KMI)
  NSV_LIMA_NC       = NSV_LIMA_NC_A(KMI)
  NSV_LIMA_NR       = NSV_LIMA_NR_A(KMI)
  NSV_LIMA_CCN_FREE = NSV_LIMA_CCN_FREE_A(KMI)
  NSV_LIMA_CCN_ACTI = NSV_LIMA_CCN_ACTI_A(KMI)
  NSV_LIMA_SCAVMASS = NSV_LIMA_SCAVMASS_A(KMI)
  NSV_LIMA_NI       = NSV_LIMA_NI_A(KMI)
  NSV_LIMA_NS       = NSV_LIMA_NS_A(KMI)
  NSV_LIMA_NG       = NSV_LIMA_NG_A(KMI)
  NSV_LIMA_NH       = NSV_LIMA_NH_A(KMI)
  NSV_LIMA_IFN_FREE = NSV_LIMA_IFN_FREE_A(KMI)
  NSV_LIMA_IFN_NUCL = NSV_LIMA_IFN_NUCL_A(KMI)
  NSV_LIMA_IMM_NUCL = NSV_LIMA_IMM_NUCL_A(KMI)
  NSV_LIMA_HOM_HAZE = NSV_LIMA_HOM_HAZE_A(KMI)
  NSV_LIMA_SPRO     = NSV_LIMA_SPRO_A(KMI)
ENDIF
!
END SUBROUTINE LIMA_UPDATE_NSV
!
END MODULE MODE_LIMA_UPDATE_NSV
