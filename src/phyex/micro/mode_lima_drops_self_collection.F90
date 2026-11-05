!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_DROPS_SELF_COLLECTION
  IMPLICIT NONE
CONTAINS
!     #############################################################
  SUBROUTINE LIMA_DROPS_SELF_COLLECTION (LIMAP, LIMAW, KSIZE, ODCOMPUTE,    &
                                         PRHODREF,            &
                                         PCRT, PLBDR, PLBDR3, &
                                         P_CR_SCBU            )
!     #############################################################
!
!!    PURPOSE
!!    -------
!!      Compute the self-collection and physical break-up of rain drops
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
USE MODD_PARAM_LIMA_WARM, ONLY:PARAM_LIMA_WARM_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,              INTENT(IN)    :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF  ! Reference Exner function
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT      ! Rain drops C. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR     ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR3    ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CR_SCBU
!
!*       0.2   Declarations of local variables :
!
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
REAL, DIMENSION(SIZE(PCRT)) :: &
                                           ZW1, & ! work arrays
                                           ZW2, &
                                           ZW3, &
                                           ZW4, &
                                           ZSCBU
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!*       1.     Rain drops self-collection and break-up
!               ---------------------------------------
!
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPS_SELF_COLLECTION', 0, ZHOOK_HANDLE)
P_CR_SCBU(:)=0.
!
ZW4(:) =0.
!
WHERE( PCRT(:)>LIMAP%XCTMIN(3) .AND. ODCOMPUTE(:) )
   ZW4(:)  = LIMAW%XACCR1 / PLBDR(:)                ! Mean diameter
END WHERE
ZSCBU(:)=1.
WHERE (ZW4(:)>=LIMAW%XSCBU_EFF1 .AND. PCRT(:)>LIMAP%XCTMIN(3) .AND. ODCOMPUTE(:)) &
     ZSCBU(:) = EXP(LIMAW%XSCBUEXP1*(ZW4(:)-LIMAW%XSCBU_EFF1))            ! coalescence efficiency
WHERE (ZW4(:)>=LIMAW%XSCBU_EFF2 .AND. ODCOMPUTE(:)) ZSCBU(:) = 0.0  ! Break-up
!
ZW1(:) = 0.0
ZW2(:) = 0.0
ZW3(:) = 0.0
!
WHERE ( PCRT(:)>LIMAP%XCTMIN(3) .AND. ZW4(:)>1.E-4 .AND. ODCOMPUTE(:))  ! analytical integration
   ZW1(:) = LIMAW%XSCBU2 * PCRT(:)**2 / PLBDR3(:)                        ! D>100 10-6 m
   ZW3(:) = ZW1(:)*ZSCBU(:)
END WHERE
!
WHERE ( PCRT(:)>LIMAP%XCTMIN(3) .AND. ZW4(:)<=1.E-4 .AND. ODCOMPUTE(:))
   ZW2(:) = LIMAW%XSCBU3 *(PCRT(:) / PLBDR3(:))**2                       ! D<100 10-6 m
   ZW3(:) = ZW2(:)
END WHERE
!
P_CR_SCBU(:) = - ZW3(:) * PRHODREF(:)
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPS_SELF_COLLECTION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPS_SELF_COLLECTION
END MODULE MODE_LIMA_DROPS_SELF_COLLECTION
