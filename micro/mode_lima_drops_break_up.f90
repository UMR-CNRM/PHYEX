!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_LIMA_DROPS_BREAK_UP
  IMPLICIT NONE
CONTAINS
!     ##########################################
  SUBROUTINE LIMA_DROPS_BREAK_UP (LIMAP, LIMAW, KSIZE, ODCOMPUTE,  &
                                  PCRT, PRRT, &
                                  P_CR_BRKU,  &
                                  PB_CR       )
!     ##########################################
!
!!
!!    PURPOSE
!!    -------
!!      Numerical filter to prevent drops from growing too much
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
!!
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                   INTENT(IN)    :: KSIZE
LOGICAL, DIMENSION(KSIZE), INTENT(IN)    :: ODCOMPUTE  
!
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PCRT             !
REAL, DIMENSION(KSIZE),    INTENT(IN)    :: PRRT             !
!
REAL, DIMENSION(KSIZE),    INTENT(OUT)   :: P_CR_BRKU        ! Concentration change (#/kg)
REAL, DIMENSION(KSIZE),    INTENT(INOUT) :: PB_CR            ! Cumulated concentration change (#/kg)
!
!*       0.2   Declarations of local variables :
!
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
REAL,    DIMENSION(SIZE(PCRT)) :: ZWLBDR,ZWLBDR3
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!              SPONTANEOUS BREAK-UP (NUMERICAL FILTER)
!              ---------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPS_BREAK_UP', 0, ZHOOK_HANDLE)
P_CR_BRKU(:)=0.
!
ZWLBDR3(:) = 1.E30
ZWLBDR(:) = 1.E10
WHERE ( PRRT(:)>LIMAP%XRTMIN(3) .AND. PCRT(:)>LIMAP%XCTMIN(3) .AND. ODCOMPUTE(:) )
   ZWLBDR3(:) = LIMAW%XLBR * PCRT(:) / PRRT(:)
   ZWLBDR(:)  = ZWLBDR3(:)**LIMAW%XLBEXR
END WHERE
WHERE (ZWLBDR(:)<(LIMAW%XACCR1/LIMAW%XSPONBUD1) .AND. ODCOMPUTE(:))
   P_CR_BRKU(:) = PCRT(:)*( MAX((1.+LIMAW%XSPONCOEF2*(LIMAW%XACCR1/ZWLBDR(:)-LIMAW%XSPONBUD1)**2),&
                                                     (LIMAW%XACCR1/ZWLBDR(:)/LIMAW%XSPONBUD3)**3) -1. )
END WHERE
!
PB_CR(:) = PB_CR(:) + P_CR_BRKU(:)
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPS_BREAK_UP', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPS_BREAK_UP
END MODULE MODE_LIMA_DROPS_BREAK_UP
