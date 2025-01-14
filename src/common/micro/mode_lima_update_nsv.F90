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
SUBROUTINE LIMA_UPDATE_NSV(LIMAP, TNSV, ODINIT, KMI, KSV, HDCLOUD, ODUPDATE)
!!*** *LIMA_UPDATE_NSV* - update modd_nsv values realtive to LIMA
!!
!!*   PURPOSE
!!    -------
!!    The modd_nsv values relative to the LIMA scheme are initialised (if ODINIT is .TRUE)
!!    according to the micromisics scheme used.
!!    If ODUPDATE is .TRUE., the scalar values of modd_nsv module receive the values
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
USE MODD_NSV, ONLY: NSV_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
!* 0.1. Declaration of arguments
!       ------------------------
!
IMPLICIT NONE
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(NSV_T),      INTENT(INOUT) :: TNSV
LOGICAL,          INTENT(IN)    :: ODINIT   !< .TRUE. to fill the different NSV_LIMA_*_A arrays
INTEGER,          INTENT(IN)    :: KMI      !< model number
INTEGER,          INTENT(INOUT) :: KSV      !< IN: Initial value to use when filling the NSV_LIMA_*_A arrays; 
                                            !! OUT: Final value after having filled the arrays
CHARACTER(LEN=4), INTENT(IN)    :: HDCLOUD  !< Cloud scheme
LOGICAL,          INTENT(IN)    :: ODUPDATE !< .TRUE. to goto model
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!*      1. INITIALISATION
!       -----------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_UPDATE_NSV', 0, ZHOOK_HANDLE)
IF(ODINIT) THEN
  IF(HDCLOUD=='LIMA') THEN
    KSV = KSV+1
    TNSV%NSV_LIMA_BEG_A(KMI) = KSV
    ! Nc
    IF (LIMAP%NMOM_C.GE.2) THEN
      TNSV%NSV_LIMA_NC_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Nr
    IF (LIMAP%NMOM_R.GE.2) THEN
      TNSV%NSV_LIMA_NR_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! CCN
    IF (LIMAP%NMOD_CCN .GT. 0) THEN
      TNSV%NSV_LIMA_CCN_FREE_A(KMI) = KSV
      KSV = KSV + LIMAP%NMOD_CCN
      TNSV%NSV_LIMA_CCN_ACTI_A(KMI) = KSV
      KSV = KSV + LIMAP%NMOD_CCN
    END IF
    ! Scavenging
    IF (LIMAP%LSCAV .AND. LIMAP%LAERO_MASS) THEN
      TNSV%NSV_LIMA_SCAVMASS_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Ni
    IF (LIMAP%NMOM_I.GE.2) THEN
      TNSV%NSV_LIMA_NI_A(KMI) = KSV
      IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
        KSV = KSV+1
      ELSE
        KSV = KSV + LIMAP%NNB_CRYSTAL_SHAPE
      END IF
    END IF
    ! Ns
    IF (LIMAP%NMOM_S.GE.2) THEN
      TNSV%NSV_LIMA_NS_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Ng
    IF (LIMAP%NMOM_G.GE.2) THEN
      TNSV%NSV_LIMA_NG_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! Nh
    IF (LIMAP%NMOM_H.GE.2) THEN
      TNSV%NSV_LIMA_NH_A(KMI) = KSV
      KSV = KSV+1
    END IF
    ! IFN
    IF (LIMAP%NMOD_IFN .GT. 0) THEN
      TNSV%NSV_LIMA_IFN_FREE_A(KMI) = KSV
      KSV = KSV + LIMAP%NMOD_IFN
      TNSV%NSV_LIMA_IFN_NUCL_A(KMI) = KSV
      KSV = KSV + LIMAP%NMOD_IFN
    END IF
    ! IMM
    IF (LIMAP%NMOD_IMM .GT. 0) THEN
      TNSV%NSV_LIMA_IMM_NUCL_A(KMI) = KSV
      KSV = KSV + MAX(1,LIMAP%NMOD_IMM)
    END IF
    ! Homogeneous freezing of CCN
    IF (LIMAP%LHHONI) THEN
      TNSV%NSV_LIMA_HOM_HAZE_A(KMI) = KSV
      KSV = KSV + 1
    END IF
    ! Supersaturation
    IF (LIMAP%LSPRO) THEN
      TNSV%NSV_LIMA_SPRO_A(KMI) = KSV
      KSV = KSV + 1
    END IF
    !
    ! End and total variables
    !
    KSV = KSV - 1
    TNSV%NSV_LIMA_END_A(KMI) = KSV
    TNSV%NSV_LIMA_A(KMI) = TNSV%NSV_LIMA_END_A(KMI) - TNSV%NSV_LIMA_BEG_A(KMI) + 1
  ELSE
    TNSV%NSV_LIMA_A(KMI)    = 0
    !
    ! force First index to be superior to last index
    ! in order to create a null section
    !
    TNSV%NSV_LIMA_BEG_A(KMI) = 1
    TNSV%NSV_LIMA_END_A(KMI) = 0
  ENDIF
ENDIF
!
!*      2. UPDATE
!       ---------
!
IF(ODUPDATE) THEN
  TNSV%NSV_LIMA          = TNSV%NSV_LIMA_A(KMI)
  TNSV%NSV_LIMA_BEG      = TNSV%NSV_LIMA_BEG_A(KMI)
  TNSV%NSV_LIMA_END      = TNSV%NSV_LIMA_END_A(KMI)
  TNSV%NSV_LIMA_NC       = TNSV%NSV_LIMA_NC_A(KMI)
  TNSV%NSV_LIMA_NR       = TNSV%NSV_LIMA_NR_A(KMI)
  TNSV%NSV_LIMA_CCN_FREE = TNSV%NSV_LIMA_CCN_FREE_A(KMI)
  TNSV%NSV_LIMA_CCN_ACTI = TNSV%NSV_LIMA_CCN_ACTI_A(KMI)
  TNSV%NSV_LIMA_SCAVMASS = TNSV%NSV_LIMA_SCAVMASS_A(KMI)
  TNSV%NSV_LIMA_NI       = TNSV%NSV_LIMA_NI_A(KMI)
  TNSV%NSV_LIMA_NS       = TNSV%NSV_LIMA_NS_A(KMI)
  TNSV%NSV_LIMA_NG       = TNSV%NSV_LIMA_NG_A(KMI)
  TNSV%NSV_LIMA_NH       = TNSV%NSV_LIMA_NH_A(KMI)
  TNSV%NSV_LIMA_IFN_FREE = TNSV%NSV_LIMA_IFN_FREE_A(KMI)
  TNSV%NSV_LIMA_IFN_NUCL = TNSV%NSV_LIMA_IFN_NUCL_A(KMI)
  TNSV%NSV_LIMA_IMM_NUCL = TNSV%NSV_LIMA_IMM_NUCL_A(KMI)
  TNSV%NSV_LIMA_HOM_HAZE = TNSV%NSV_LIMA_HOM_HAZE_A(KMI)
  TNSV%NSV_LIMA_SPRO     = TNSV%NSV_LIMA_SPRO_A(KMI)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('LIMA_UPDATE_NSV', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_UPDATE_NSV
!
END MODULE MODE_LIMA_UPDATE_NSV       
