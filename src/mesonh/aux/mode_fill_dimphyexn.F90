!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_FILL_DIMPHYEX
IMPLICIT NONE
CONTAINS
SUBROUTINE FILL_DIMPHYEX(YDDIMPHYEX, KIT, KJT, KKT)
!     #########################
!
!!
!!    PURPOSE
!!    -------
!       This subroutine computes the dimensions according to the running model.
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
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPVEXT
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t), INTENT(OUT) :: YDDIMPHYEX ! Structure to fill in
INTEGER, INTENT(IN) :: KIT, KJT, KKT ! Array dimensions

!
!*       0.2  declaration of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FILL_DIMPHYEX', 0, ZHOOK_HANDLE)
!
YDDIMPHYEX%NIT=KIT
YDDIMPHYEX%NJT=KJT
CALL GET_INDICE_ll(YDDIMPHYEX%NIB, YDDIMPHYEX%NJB,&
                  &YDDIMPHYEX%NIE, YDDIMPHYEX%NJE)
!
YDDIMPHYEX%NKL=1
YDDIMPHYEX%NKT=KKT
YDDIMPHYEX%NKA=1
YDDIMPHYEX%NKU=KKT
YDDIMPHYEX%NKB=1+JPVEXT
YDDIMPHYEX%NKE=KKT-JPVEXT
YDDIMPHYEX%NKTB=1+JPVEXT
YDDIMPHYEX%NKTE=KKT-JPVEXT
!
IF (LHOOK) CALL DR_HOOK('FILL_DIMPHYEX', 1, ZHOOK_HANDLE)
!
END SUBROUTINE FILL_DIMPHYEX
END MODULE MODE_FILL_DIMPHYEX
