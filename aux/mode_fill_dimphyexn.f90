!MNH_LIC Copyright 2022-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_FILL_DIMPHYEX
IMPLICIT NONE
CONTAINS
SUBROUTINE FILL_DIMPHYEX( YDDIMPHYEX, KIT, KJT, KKT, LTURB, KLES_TIMES, KLES_K )
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
!  P. Wautelet 04/10/2023: bugfix: set NKLES correctly
!
!-----------------------------------------------------------------
!*       0.   DECLARATIONS
!             ------------
!
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_LES, ONLY : LLES_NEB_MASK, LLES_CORE_MASK, LLES_CS_MASK, LLES_MY_MASK, &
                     NLES_MASKS_USER
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t), INTENT(OUT) :: YDDIMPHYEX ! Structure to fill in
INTEGER, INTENT(IN) :: KIT, KJT, KKT ! Array dimensions
LOGICAL, INTENT(IN), OPTIONAL :: LTURB ! Flag to replace array dimensions I/JB and I/JE to the full array size
                                       ! needed if computation in HALO points (e.g. in turbulence)
INTEGER, INTENT(IN), OPTIONAL :: KLES_TIMES  ! number of LES computations in time
INTEGER, INTENT(IN), OPTIONAL :: KLES_K      ! number of vertical levels for LES diagnostics

LOGICAL :: YTURB
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
YDDIMPHYEX%NIJT=KIT*KJT
CALL GET_INDICE_ll(YDDIMPHYEX%NIB, YDDIMPHYEX%NJB,&
                  &YDDIMPHYEX%NIE, YDDIMPHYEX%NJE)
!
YDDIMPHYEX%NIJB=1
YDDIMPHYEX%NIJE=KIT*KJT
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
IF(PRESENT(LTURB)) THEN
  YTURB=LTURB
ELSE
  YTURB=.FALSE. ! Default value
END IF
!
IF(YTURB) THEN
  YDDIMPHYEX%NIBC=1
  YDDIMPHYEX%NJBC=1
  YDDIMPHYEX%NIEC=YDDIMPHYEX%NIT
  YDDIMPHYEX%NJEC=YDDIMPHYEX%NJT
ELSE
  YDDIMPHYEX%NIBC=YDDIMPHYEX%NIB
  YDDIMPHYEX%NJBC=YDDIMPHYEX%NJB
  YDDIMPHYEX%NIEC=YDDIMPHYEX%NIE
  YDDIMPHYEX%NJEC=YDDIMPHYEX%NJE
END IF
IF (LHOOK) CALL DR_HOOK('FILL_DIMPHYEX', 1, ZHOOK_HANDLE)
!
YDDIMPHYEX%NLESMASK=1
YDDIMPHYEX%NLES_TIMES=0
IF (PRESENT(KLES_TIMES)) THEN
  YDDIMPHYEX%NLES_TIMES = KLES_TIMES
END IF
YDDIMPHYEX%NKLES=0
IF (PRESENT(KLES_K)) THEN
  YDDIMPHYEX%NKLES = KLES_K
END IF
IF (LLES_MY_MASK) YDDIMPHYEX%NLESMASK = YDDIMPHYEX%NLESMASK + NLES_MASKS_USER
IF (LLES_NEB_MASK) YDDIMPHYEX%NLESMASK = YDDIMPHYEX%NLESMASK + 2
IF (LLES_CORE_MASK) YDDIMPHYEX%NLESMASK = YDDIMPHYEX%NLESMASK + 2
IF (LLES_CS_MASK) YDDIMPHYEX%NLESMASK = YDDIMPHYEX%NLESMASK + 3
!
END SUBROUTINE FILL_DIMPHYEX
END MODULE MODE_FILL_DIMPHYEX
