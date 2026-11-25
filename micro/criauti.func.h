!MNH_LIC Copyright 1994-2026 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
ELEMENTAL SUBROUTINE CRIAUTI(PRCRIAUTI, PT0CRIAUTI, PACRIAUTI, PBCRIAUTI)
!!
!!**  PURPOSE
!!    -------
!!      Compute the A and B coefficients of the 10**(A*T+B) law used for the
!!      ice->snow autoconversion threshold (cirrus law), from the threshold
!!      value PRCRIAUTI (possibly perturbed, e.g. by SPP) and the temperature
!!      PT0CRIAUTI at which the cirrus law starts to be used.
!!      The second point needed to determine the law is fixed (-40 degC, 1.25E-6).
!!
!!    AUTHOR
!!    ------
!!      U. Andrae (2026)
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!

!
IMPLICIT NONE
!
REAL,    INTENT(IN)    :: PRCRIAUTI   ! Constants for pristine ice autoconversion
REAL,    INTENT(IN)    :: PT0CRIAUTI  ! Temp degC at which cirrus law starts to be used
REAL,    INTENT(OUT)   :: PACRIAUTI   ! A Coef. for cirrus law
REAL,    INTENT(OUT)   :: PBCRIAUTI   ! B Coef. for cirrus law
!
REAL :: ZTCRI0, ZCRI0
!
! Second point to determine 10**(aT+b) law
!
ZTCRI0=-40.0
ZCRI0=1.25E-6
!
PBCRIAUTI=-( LOG10(PRCRIAUTI) - LOG10(ZCRI0)*PT0CRIAUTI/ZTCRI0 )&
            *ZTCRI0/(PT0CRIAUTI-ZTCRI0)
PACRIAUTI=(LOG10(ZCRI0)-PBCRIAUTI)/ZTCRI0
!
END SUBROUTINE CRIAUTI
