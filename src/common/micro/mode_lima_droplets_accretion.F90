!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_DROPLETS_ACCRETION
  IMPLICIT NONE
CONTAINS
!     #####################################################################
  SUBROUTINE LIMA_DROPLETS_ACCRETION (LDCOMPUTE,                      &
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
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, LKHKO, NMOM_C, NMOM_R, XCEXVT
USE MODD_PARAM_LIMA_WARM, ONLY : XLAUTR, XAUTO1, XLAUTR_THRESHOLD, &
                                 XACCR4, XACCR5, XACCR3, XACCR2, XACCR1, &
                                 XACCR_CLARGE1, XACCR_CLARGE2, XACCR_RLARGE1, XACCR_RLARGE2, &
                                 XACCR_CSMALL1, XACCR_CSMALL2, XACCR_RSMALL1, XACCR_RSMALL2, &
                                 XFCACCR, XEXCACCR
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF ! Reference Exner function
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT     ! Rain m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCCT     ! Cloud water conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PCRT     ! Rain conc. at t
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDC3   ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDR3   ! 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RC_ACCR
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CC_ACCR
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRCT))    :: ZW1, ZW2, ZW3, ZW4 ! work arrays
LOGICAL, DIMENSION(SIZE(PRCT)) :: GACCR
!
!-------------------------------------------------------------------------------
!
!
!
!*       1. Accretion of cloud droplets on rain drops
!   	 --------------------------------------------
!
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
IF ( LKHKO ) THEN
!
   GACCR(:) = PRRT(:)>XRTMIN(3) .AND. &
              PCRT(:)>XCTMIN(3) .AND. &
              PRCT(:)>XRTMIN(2) .AND. &
              PCCT(:)>XCTMIN(2)
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
ELSE IF (NMOM_C.EQ.1 .AND. NMOM_R.EQ.1) THEN
   GACCR(:) = PRRT(:)>XRTMIN(3) .AND. PCRT(:)>XCTMIN(3) .AND. &
              PRCT(:)>XRTMIN(2) .AND. PCCT(:)>XCTMIN(2)
   WHERE ( GACCR(:) )
      P_RC_ACCR(:) = - XFCACCR * PRCT(:)    &
                   * PLBDR(:)**XEXCACCR     &
                   * PRHODREF(:)**(-XCEXVT)
   END WHERE
ELSE
!
   WHERE( PRCT(:)>XRTMIN(2) .AND. PCCT(:)>XCTMIN(2) .AND. PRRT(:)>XRTMIN(3) .AND. PCRT(:)>XCTMIN(3) .AND. LDCOMPUTE(:) )
      ZW2(:) = MAX( 0.0,XLAUTR*PRHODREF(:)*PRCT(:)*(XAUTO1/PLBDC(:)**4-XLAUTR_THRESHOLD) ) ! L 
      ZW4(:) = XACCR1/PLBDR(:)
   END WHERE
!
   GACCR(:) = LDCOMPUTE(:)      .AND. &
              PRRT(:)>XRTMIN(3) .AND. &
              PCRT(:)>XCTMIN(3) .AND. &
              PRCT(:)>XRTMIN(2) .AND. &
              PCCT(:)>XCTMIN(2) .AND. &
              (PRRT(:)>1.2*ZW2(:)/PRHODREF(:) .OR. &
                ZW4(:)>=MAX(XACCR2,XACCR3/(XACCR4/PLBDC(:)-XACCR5)) )
!
! Accretion for D>100 10-6 m
   WHERE( GACCR(:).AND.(ZW4(:)>1.E-4) )
      ZW3(:) = MIN(PLBDC3(:) / PLBDR3(:),1.E15)
      ZW1(:) = ( PCCT(:)*PCRT(:) / PLBDC3(:) )*PRHODREF(:)
      ZW2(:) = ZW1(:)*(XACCR_CLARGE1+XACCR_CLARGE2*ZW3(:))
!
      P_CC_ACCR(:) = - ZW2(:)
!
      ZW1(:) = ( ZW1(:) / PLBDC3(:) )
      ZW2(:) = ZW1(:)*(XACCR_RLARGE1+XACCR_RLARGE2*ZW3(:))
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
      ZW2(:) = ZW1(:)*(XACCR_CSMALL1+XACCR_CSMALL2*ZW3(:))
!
      P_CC_ACCR(:) = - ZW2(:)
!
      ZW1(:) = ZW1(:) / PLBDC3(:)
      ZW2(:) = ZW1(:)*(XACCR_RSMALL1+XACCR_RSMALL2*ZW3(:))
!
      P_RC_ACCR(:) = - ZW2(:)
   END WHERE
!
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPLETS_ACCRETION
END MODULE MODE_LIMA_DROPLETS_ACCRETION
