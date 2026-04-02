!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_FILL_DIMPHYEX
IMPLICIT NONE
CONTAINS
SUBROUTINE FILL_DIMPHYEX(YDDIMPHYEX, KIT, KJT, KKT, KVEXT, KIB, KIE)
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
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t), INTENT(OUT) :: YDDIMPHYEX ! Structure to fill in
INTEGER, INTENT(IN) :: KIT, KJT, KKT ! Array dimensions
INTEGER, INTENT(IN) :: KVEXT ! Number of unphysical points at each end of the vertical axis
INTEGER, INTENT(IN) :: KIB ! Index of the first horizontal point to consider
INTEGER, INTENT(IN) :: KIE ! Index of the last horizontal point to consider

!
!*       0.2  declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('FILL_DIMPHYEX', 0, ZHOOK_HANDLE)
!
YDDIMPHYEX%NIT=KIT
YDDIMPHYEX%NIB=KIB
YDDIMPHYEX%NIE=KIE !used to be KIT before considering KPROMA blocs
!
!In AROME, KJT is always 1
YDDIMPHYEX%NJT=KJT
YDDIMPHYEX%NJB=1
YDDIMPHYEX%NJE=KJT
!
YDDIMPHYEX%NIJT=KIT*KJT
YDDIMPHYEX%NIJB=KIB
YDDIMPHYEX%NIJE=KIE
!
YDDIMPHYEX%NKL=-1
YDDIMPHYEX%NKT=KKT
YDDIMPHYEX%NKA=KKT
YDDIMPHYEX%NKU=1
YDDIMPHYEX%NKB=KKT-KVEXT
YDDIMPHYEX%NKE=1+KVEXT
YDDIMPHYEX%NKLES=KKT-2*KVEXT
YDDIMPHYEX%NKTB=1+KVEXT
YDDIMPHYEX%NKTE=KKT-KVEXT
!
YDDIMPHYEX%NIBC=1
YDDIMPHYEX%NJBC=1
YDDIMPHYEX%NIEC=KIE
YDDIMPHYEX%NJEC=KJT
!
YDDIMPHYEX%NLESMASK = 0 ! never used in AROME
IF (LHOOK) CALL DR_HOOK('FILL_DIMPHYEX', 1, ZHOOK_HANDLE)
!
END SUBROUTINE FILL_DIMPHYEX
END MODULE MODE_FILL_DIMPHYEX
