!MNH_LIC Copyright 1995-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_STEPPING
IMPLICIT NONE

INTERFACE

SUBROUTINE ICE4_STEPPING(CST, PARAMI, ICEP, ICED, BUCONF, &
                        &KPROMA, KSIZE , PTSTEP, &
                        &KRR, OSAVE_MICRO, LDMICRO, OELEC, &
                        &PEXN, PRHODREF, PPABST, PCIT, PCLDFR, &
                        &PHLC_HCF, PHLC_HRC, PHLI_HCF, PHLI_HRI, &
                        &PTHS, PRS, PRREVAV, PRAINFR, PSIGS, PTHT, PRT, &
                        &PICLDFR, PZZZ, PCONC3D, PSSIO, PSSIU, PIFR, &
                        &PBUDGETS, PLATHAM_IAGGS)

USE MODD_BUDGET,         ONLY: TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_FIELDS_ADDRESS, ONLY : & ! common fields adress
      & IRV,     & ! Water vapor
      & IRC,     & ! Cloud water
      & IRR,     & ! Rain water
      & IRI,     & ! Pristine ice
      & IBUNUM,       & ! Number of tendency terms
      & IBUNUM_MR,    & ! Number of tendency terms expressed as mixing ratio changes
      & IBUNUM_EXTRA, & ! Number of extra tendency terms
      & IRREVAV,      & ! Index for the evaporation tendency
      & IBUEXTRAIND
! Index indirection
!
USE MODE_ICE4_TENDENCIES, ONLY: ICE4_TENDENCIES
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
INTEGER,                  INTENT(IN)    :: KPROMA ! cache-blocking factor for microphysic loop
INTEGER,                  INTENT(IN)    :: KSIZE  ! Case r_x>0 locations
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL,                  INTENT(IN)    :: OSAVE_MICRO   ! if true, save the microphysical tendencies
LOGICAL, DIMENSION(KPROMA), INTENT(IN)  :: LDMICRO
LOGICAL,                  INTENT(IN)    :: OELEC         ! if true, cloud electricity is activated
!
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PEXN    ! Exner function
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PRHODREF! Reference density
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PPABST

REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PCIT

REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PCLDFR ! Cloud fraction
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PICLDFR ! Ice cloud fraction
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PZZZ    ! Model level height
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PCONC3D ! Cloud croplet number concentration
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the supersaturated fraction
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the  subsaturated fraction
REAL,    DIMENSION(KPROMA),                     INTENT(IN)    :: PIFR    ! Ratio cloud ice moist part to dry part
REAL,    DIMENSION(MERGE(KPROMA,0,PARAMI%CSUBG_AUCV_RC=='ADJU' .OR. PARAMI%CSUBG_AUCV_RI=='ADJU')), INTENT(INOUT) :: PHLC_HRC, &
                                                                                                                  & PHLC_HCF, &
                                                                                                                  & PHLI_HRI, &
                                                                                                                  & PHLI_HCF
REAL,    DIMENSION(MERGE(KPROMA,0,PARAMI%LEXT_TEND)),   INTENT(IN) :: PTHS !To take into acount external tendencies inside the splitting
REAL,    DIMENSION(MERGE(KPROMA,0,PARAMI%LEXT_TEND),7), INTENT(IN) :: PRS !To take into acount external tendencies inside the splitting
REAL,    DIMENSION(KPROMA),                     INTENT(OUT)   :: PRREVAV
REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PRAINFR
REAL,    DIMENSION(MERGE(KPROMA,0,PARAMI%CSUBG_AUCV_RC=='PDF ' .AND. PARAMI%CSUBG_PR_PDF=='SIGM')), INTENT(IN) :: PSIGS
REAL,    DIMENSION(KPROMA),                     INTENT(INOUT) :: PTHT
REAL,    DIMENSION(KPROMA,7),                   INTENT(INOUT) :: PRT !Packed variables
REAL,    DIMENSION(MERGE(KPROMA,0,BUCONF%LBU_ENABLE .OR. OSAVE_MICRO), &
                   MERGE(IBUNUM-IBUNUM_EXTRA,0,BUCONF%LBU_ENABLE .OR. OSAVE_MICRO)),INTENT(OUT)   :: PBUDGETS
REAL,    DIMENSION(MERGE(KPROMA,0,OELEC)),      INTENT(IN)    :: PLATHAM_IAGGS ! E Function to simulate
                                                                               ! enhancement of IAGGS
!
!
END SUBROUTINE ICE4_STEPPING

END INTERFACE

END MODULE MODI_ICE4_STEPPING
