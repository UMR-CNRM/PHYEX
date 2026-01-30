!MNH_LIC Copyright 2022-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      MODULE MODD_DIMPHYEX
!     ####################
!
!!****  *MODD_DIMPHYEX* - declaration of dimensions for the physics
!!
!!    PURPOSE
!!    -------
!       Declaration of array dimensions used by the physics
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
!!      S. Riette, Météo-France
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    January 2022
!
!-----------------------------------------------------------------
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE
!
TYPE DIMPHYEX_t
  ! 
  !On x direction
  INTEGER :: NIT ! Array total dimension
  INTEGER :: NIB ! First inner mass point index
  INTEGER :: NIE ! Last inner mass point index
  !
  !On y direction
  INTEGER :: NJT ! Array total dimension
  INTEGER :: NJB ! First inner mass point index
  INTEGER :: NJE ! Last inner mass point index
  !
  !On z direction
  !Ordering can be different depending on the host model
  INTEGER :: NKL  ! Order of the vertical levels
                  !  1: as for Méso-NH, levels are numbered from ground to space
                  ! -1: as for AROME, levels are numbered from space to ground
  INTEGER :: NKT  ! Array total dimension
  INTEGER :: NKLES ! Number of vertical levels for LES diagnostics
  INTEGER :: NKA  ! Near ground array index (is an unphysical level if JPVEXT!=0)
  INTEGER :: NKU  ! Uppest atmosphere array index (is an unphysical level if JPVEXT!=0)
  INTEGER :: NKB  ! Near ground physical array index (e.g. equal to 1+JPVEXT if NKL==1)
  INTEGER :: NKE  ! Uppest physical atmosphere array index (e.g. equal to 1+JPVEXT if NKL==-1)
  INTEGER :: NKTB ! Smaller index of the physical domain (equals to MIN(NKB, NKE)=1+JPVEXT)
  INTEGER :: NKTE ! Greater index of the physical domain (equals to MAX(NKB, NKE)=NKT-JPVEXT)
  !Explanations about the different values. To loop on:
  !* all (including non physical) levels from ground to top of atm: DO JK=NKA, NKU, KKL
  !* all (including non physical) levels from top of atm to ground: DO JK=NKU, NKA, -KKL
  !* physical levels only from ground to top of atm: DO JK=NKB, NKE, KKL
  !* physical levels only from top of atm to ground: DO JK=NKE, NKB, -KKL
  !* all (including non physical) following the array ordering: DO JK=1, NKT
  !* physical levels only following the array ordering: DO JK=NKTB, NKTE
  INTEGER :: NIBC  ! Computational indices used in DO LOOP
  INTEGER :: NJBC  ! = NIB/NJC/NIE/NJE in all schemes
  INTEGER :: NIEC  ! except in turbulence where external HALO points must be
  INTEGER :: NJEC  ! included so NIBC=NJBC=1 and NIEC/NJEC=NIT/NJT
  INTEGER :: NIJT  ! NIT*NJT for horizontal packing
  INTEGER :: NIJB  ! First horizontal inner mass point index
  INTEGER :: NIJE  ! Last horizontal inner mass point index
  !
  INTEGER :: NLESMASK ! Number of LES masks
  INTEGER :: NLES_TIMES ! Number of LES time data storage
!
END TYPE DIMPHYEX_t
!
END MODULE MODD_DIMPHYEX

