!MNH_LIC Copyright 1995-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_PACK
IMPLICIT NONE

INTERFACE

SUBROUTINE ICE4_PACK(D, CST, PARAMI, ICEP, ICED, BUCONF, &
                    &PTSTEP, &
                    &KRR, OSAVE_MICRO, LDMICRO, OELEC, &
                    &PEXN, PRHODREF, PPABST, PCIT, PCLDFR, &
                    &PHLC_HCF, PHLC_HRC, PHLI_HCF, PHLI_HRI, &
                    &PTHS, PRS, PRREVAV, PRAINFR, PSIGS, PTHT, PRT, &
                    &PICLDFR, PZZZ, PCONC3D, PSSIO, PSSIU, PIFR, &
                    &PBUDGETS, PLATHAM_IAGGS)
!     ######################################################################
!
!!****  * -  compute the explicit microphysical sources
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to pack arrays to compute
!!      the microphysics tendencies
!!
!!
!!    METHOD
!!    ------
!!      Pack arrays by chuncks
!!
!!
!!    MODIFICATIONS
!!    -------------
!!     R. El Khatib 28-Apr-2023 Fix (and re-enable) the cache-blocking mechanism on top of phyex
!!     S. Riette Sept 23: all 3D arrays are suppressed from ice4_stepping
!  -----------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_BUDGET,         ONLY: TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE MODD_FIELDS_ADDRESS, ONLY : & ! common fields adress
      & IBUNUM,       & ! Number of tendency terms
      & IBUNUM_EXTRA

!
IMPLICIT NONE

TYPE(DIMPHYEX_t),         INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
REAL,                     INTENT(IN)    :: PTSTEP  ! Double Time step (single if cold start)
INTEGER,                  INTENT(IN)    :: KRR     ! Number of moist variable
LOGICAL,                  INTENT(IN)    :: OSAVE_MICRO  ! If true, microphysical tendencies are saved
LOGICAL,                  INTENT(IN)    :: OELEC        ! if true, cloud electricity is activated
LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: LDMICRO ! mask to limit computation
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCLDFR  ! Cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HCF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PHLI_HRI
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(D%NIJT,D%NKT,KRR),   INTENT(INOUT) :: PRS    ! m.r. source
!
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRREVAV! Rain evap profile
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PRAINFR !Precipitation fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSIGS   ! Sigma_s at t
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PICLDFR
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PZZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PCONC3D
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSSIO   ! Super-saturation with respect to ice in the supersaturated fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PSSIU   ! Sub-saturation with respect to ice in the  subsaturated fraction 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)    :: PIFR    ! Ratio cloud ice moist part to dry part 
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(INOUT) :: PTHT
REAL, DIMENSION(D%NIJT,D%NKT,7), INTENT(INOUT) :: PRT
REAL, DIMENSION(MERGE(D%NIJT,0,OSAVE_MICRO .OR. BUCONF%LBU_ENABLE), &
                MERGE(D%NKT,0,OSAVE_MICRO .OR. BUCONF%LBU_ENABLE), &
                MERGE(IBUNUM-IBUNUM_EXTRA,0,OSAVE_MICRO .OR. BUCONF%LBU_ENABLE)), &
                                          INTENT(OUT) :: PBUDGETS  ! Microphysical tendencies
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC),MERGE(D%NKT,0,OELEC)), &
                                          INTENT(IN)    :: PLATHAM_IAGGS  ! E Function to simulate
                                                                          ! enhancement of IAGGS
END SUBROUTINE ICE4_PACK

END INTERFACE

END MODULE MODI_ICE4_PACK
