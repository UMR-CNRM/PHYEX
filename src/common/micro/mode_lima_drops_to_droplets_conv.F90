!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_LIMA_DROPS_TO_DROPLETS_CONV
  IMPLICIT NONE
CONTAINS
!     ######################################################################
  SUBROUTINE LIMA_DROPS_TO_DROPLETS_CONV (LIMAP, D, CST, PRHODREF, PRCT, PRRT, PCCT, PCRT, &
                                          P_RR_CVRC, P_CR_CVRC    )
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      Conversion of rain drops into cloud droplets if mean volume diameter < 82µm
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * CNRM *
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
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_T
USE MODD_CST,             ONLY : CST_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(DIMPHYEX_T),         INTENT(IN)    :: D
TYPE(CST_T),              INTENT(IN)    :: CST
!
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)    :: PRHODREF! Cloud water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
!
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(IN)    :: PCRT    ! Rain water C. at t
!
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(OUT)   :: P_RR_CVRC
REAL, DIMENSION(D%NIJT,D%NKT),    INTENT(OUT)   :: P_CR_CVRC
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2)) :: ZDR
!
LOGICAL, DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2)) :: ZMASKR, ZMASKC
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPS_TO_DROPLETS_CONV', 0, ZHOOK_HANDLE)
P_RR_CVRC(:,:) = 0.
P_CR_CVRC(:,:) = 0.
!
ZDR(:,:) = 9999.
ZMASKR(:,:) = PRRT(:,:).GT.LIMAP%XRTMIN(3) .AND. PCRT(:,:).GT.LIMAP%XCTMIN(3)
ZMASKC(:,:) = PRCT(:,:).GT.LIMAP%XRTMIN(2) .AND. PCCT(:,:).GT.LIMAP%XCTMIN(2)
WHERE(ZMASKR(:,:))
   ZDR(:,:)=(6.*PRRT(:,:)/CST%XPI/CST%XRHOLW/PCRT(:,:))**0.33
END WHERE
!
! Transfer all drops in droplets if out of cloud and Dr<82microns
!
WHERE( ZMASKR(:,:) .AND. .NOT.ZMASKC(:,:) .AND. ZDR(:,:).LT.82.E-6)
   P_RR_CVRC(:,:) = -PRRT(:,:)
   P_CR_CVRC(:,:) = -PCRT(:,:)
END WHERE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPS_TO_DROPLETS_CONV', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPS_TO_DROPLETS_CONV
END MODULE MODE_LIMA_DROPS_TO_DROPLETS_CONV
