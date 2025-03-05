!MNH_LIC Copyright 2018-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_DROPLETS_ACCRETION
  IMPLICIT NONE
CONTAINS
!     #####################################################################
  SUBROUTINE LIMA_DROPLETS_ACCRETION (LIMAP, LIMAW, KSIZE, ODCOMPUTE,               &
                                      PRHODREF,                       &
                                      PRCT, PRRT, PCCT, PCRT,         &
                                      PLBDC, PLBDC3, PLBDR, PLBDR3,   &
                                      P_RC_ACCR, P_CC_ACCR            )
!     #####################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the accretion of cloud droplets by rain drops
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             15/03/2018 
!!
!       Delbeke/Vie     03/2022 : KHKO option
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
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRRT     ! Rain m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCCT     ! Cloud water conc. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCRT     ! Rain conc. at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDC3   ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDR3   ! 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_ACCR
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_ACCR
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRCT))    :: ZW1, ZW2, ZW3, ZW4 ! work arrays
TYPE(PARAM_LIMA_WARM_T),INTENT(IN)::LIMAW
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
LOGICAL, DIMENSION(SIZE(PRCT)) :: GACCR
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
!
!
!
!*       1. Accretion of cloud droplets on rain drops
!           --------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_ACCRETION', 0, ZHOOK_HANDLE)
P_RC_ACCR(:) = 0.0
P_CC_ACCR(:) = 0.0
!
ZW1(:) = 0.0
ZW2(:) = 0.0
ZW3(:) = 0.0
ZW4(:) = 0.0
!
!
!
IF ( LIMAP%LKHKO ) THEN
!
   GACCR(:) = PRRT(:)>LIMAP%XRTMIN(3) .AND. &
              PCRT(:)>LIMAP%XCTMIN(3) .AND. &
              PRCT(:)>LIMAP%XRTMIN(2) .AND. &
              PCCT(:)>LIMAP%XCTMIN(2)
!
   WHERE ( GACCR(:) )
!
      ZW1(:) = 67.0 * ( PRCT(:) * PRRT(:) )**1.15
      P_RC_ACCR(:) = - ZW1(:)
!
      ZW2(:) = ZW1(:) * PCCT(:) / PRCT(:)
      P_CC_ACCR(:) = - ZW2(:)
!
   END WHERE
!
ELSE IF (LIMAP%NMOM_C.EQ.1 .AND. LIMAP%NMOM_R.EQ.1) THEN
   GACCR(:) = PRRT(:)>LIMAP%XRTMIN(3) .AND. PCRT(:)>LIMAP%XCTMIN(3) .AND. &
              PRCT(:)>LIMAP%XRTMIN(2) .AND. PCCT(:)>LIMAP%XCTMIN(2)
   WHERE ( GACCR(:) )
      P_RC_ACCR(:) = - LIMAW%XFCACCR * PRCT(:)    &
                   * PLBDR(:)**LIMAW%XEXCACCR     &
                   * PRHODREF(:)**(-LIMAP%XCEXVT)
   END WHERE
ELSE
!
   WHERE( PRCT(:)>LIMAP%XRTMIN(2) .AND. PCCT(:)>LIMAP%XCTMIN(2) .AND. PRRT(:)>LIMAP%XRTMIN(3) &
        .AND. PCRT(:)>LIMAP%XCTMIN(3) .AND. ODCOMPUTE(:) )
      ZW2(:) = MAX( 0.0,LIMAW%XLAUTR*PRHODREF(:)*PRCT(:)*(LIMAW%XAUTO1/PLBDC(:)**4-LIMAW%XLAUTR_THRESHOLD) ) ! L 
      ZW4(:) = LIMAW%XACCR1/PLBDR(:)
   END WHERE
!
   GACCR(:) = ODCOMPUTE(:)      .AND. &
              PRRT(:)>LIMAP%XRTMIN(3) .AND. &
              PCRT(:)>LIMAP%XCTMIN(3) .AND. &
              PRCT(:)>LIMAP%XRTMIN(2) .AND. &
              PCCT(:)>LIMAP%XCTMIN(2) .AND. &
              (PRRT(:)>1.2*ZW2(:)/PRHODREF(:) .OR. &
                ZW4(:)>=MAX(LIMAW%XACCR2,LIMAW%XACCR3/(LIMAW%XACCR4/PLBDC(:)-LIMAW%XACCR5)) )
!
! Accretion for D>100 10-6 m
   WHERE( GACCR(:).AND.(ZW4(:)>1.E-4) )
      ZW3(:) = MIN(PLBDC3(:) / PLBDR3(:),1.E15)
      ZW1(:) = ( PCCT(:)*PCRT(:) / PLBDC3(:) )*PRHODREF(:)
      ZW2(:) = ZW1(:)*(LIMAW%XACCR_CLARGE1+LIMAW%XACCR_CLARGE2*ZW3(:))
!
      P_CC_ACCR(:) = - ZW2(:)
!
      ZW1(:) = ( ZW1(:) / PLBDC3(:) )
      ZW2(:) = ZW1(:)*(LIMAW%XACCR_RLARGE1+LIMAW%XACCR_RLARGE2*ZW3(:))
!
      P_RC_ACCR(:) = - ZW2(:)
   END WHERE
!
! Accretion for D<100 10-6 m
   WHERE( GACCR(:).AND.(ZW4(:)<=1.E-4) )
      ZW3(:) = MIN(PLBDC3(:) / PLBDR3(:), 1.E8)
      ZW1(:) = ( PCCT(:)*PCRT(:) / PLBDC3(:) )*PRHODREF(:)
      ZW1(:) =  ZW1(:)/PLBDC3(:)
   
      ZW3(:) = ZW3(:)**2
      ZW2(:) = ZW1(:)*(LIMAW%XACCR_CSMALL1+LIMAW%XACCR_CSMALL2*ZW3(:))
!
      P_CC_ACCR(:) = - ZW2(:)
!
      ZW1(:) = ZW1(:) / PLBDC3(:)
      ZW2(:) = ZW1(:)*(LIMAW%XACCR_RSMALL1+LIMAW%XACCR_RSMALL2*ZW3(:))
!
      P_RC_ACCR(:) = - ZW2(:)
   END WHERE
!
END IF
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_DROPLETS_ACCRETION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_DROPLETS_ACCRETION
END MODULE MODE_LIMA_DROPLETS_ACCRETION
