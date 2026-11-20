!MNH_LIC Copyright 2014-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! ******* TO BE INCLUDED IN THE *CONTAINS* OF A SUBROUTINE, IN ORDER TO EASE AUTOMATIC INLINING ******
! => Don't use drHook !!!
!
! Functions of temperature (saturation pressures and their derivatives, ice crystal and
! water droplet growth functions...) read in the lookup tables of MODD_TIWMX (filled by
! INI_TIWMX), the tables being given as argument (ICEP%TIWMX)
!
! Author: K. I. Ivarsson *SMHI* (20/11/14), formerly in mode_tiwmx.F90
! Modifications: S. Riette Nov. 2026: elemental functions in tiwmx.func.h
!
  ELEMENTAL REAL FUNCTION ESATW(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    ESATW = TIWMX%ESTABW(NINT(XNDEGR*TT))
  END FUNCTION ESATW

  ELEMENTAL REAL FUNCTION DESDTW(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    DESDTW = TIWMX%DESTABW(NINT(XNDEGR*TT))
  END FUNCTION

  ELEMENTAL REAL FUNCTION ESATI(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    ESATI = TIWMX%ESTABI(NINT(XNDEGR*TT))
  END FUNCTION

  ELEMENTAL REAL FUNCTION DESDTI(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    DESDTI = TIWMX%DESTABI(NINT(XNDEGR*TT))
  END FUNCTION

! Water droplet function:
  ELEMENTAL REAL FUNCTION AA2W(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    AA2W = TIWMX%A2WTAB(NINT(XNDEGR*TT))
  END FUNCTION

! Ice crystal function
  ELEMENTAL REAL FUNCTION AA2(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    AA2 = TIWMX%A2TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Meyers IN concentration function:
  ELEMENTAL REAL FUNCTION AM3(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    AM3 = TIWMX%AM3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Fletchers IN concentration function:
  ELEMENTAL REAL FUNCTION AF3(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    AF3 = TIWMX%AF3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Ice crystal function
  ELEMENTAL REAL FUNCTION BB3(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    BB3 = TIWMX%BB3TAB(NINT(XNDEGR*TT))
  END FUNCTION

! Water droplet function:
  ELEMENTAL REAL FUNCTION BB3W(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    BB3W = TIWMX%BB3WTAB(NINT(XNDEGR*TT))
  END FUNCTION

! Function for IN concentration reduction between 0 and -25 C:
  ELEMENTAL REAL FUNCTION REDIN(TIWMX, TT)
    USE MODD_TIWMX, ONLY: XNDEGR, TIWMX_t
    IMPLICIT NONE

!$acc routine seq

    TYPE(TIWMX_t),   INTENT(IN) :: TIWMX
    REAL,INTENT(IN) :: TT 
    REDIN = TIWMX%REDINTAB(NINT(XNDEGR*TT))
  END FUNCTION REDIN
  
