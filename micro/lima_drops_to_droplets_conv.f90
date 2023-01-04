!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      #################################
       MODULE MODI_LIMA_DROPS_TO_DROPLETS_CONV
!      #################################
!
INTERFACE
      SUBROUTINE LIMA_DROPS_TO_DROPLETS_CONV (PRHODREF, PRCT, PRRT, PCCT, PCRT, &
                                              P_RR_CVRC, P_CR_CVRC    )
!
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
!
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PCRT    ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),    INTENT(INOUT) :: P_RR_CVRC
REAL, DIMENSION(:,:,:),    INTENT(INOUT) :: P_CR_CVRC
!
END SUBROUTINE LIMA_DROPS_TO_DROPLETS_CONV
END INTERFACE
END MODULE MODI_LIMA_DROPS_TO_DROPLETS_CONV
!
!     ######################################################################
      SUBROUTINE LIMA_DROPS_TO_DROPLETS_CONV (PRHODREF, PRCT, PRRT, PCCT, PCRT, &
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
USE MODD_CST,             ONLY : XPI, XRHOLW 
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN
USE MODD_PARAM_LIMA_WARM, ONLY : XLBR, XLBEXR, XLBC, XLBEXC, &
                                 XACCR1, XACCR3, XACCR4, XACCR5
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
!
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PCCT    ! Cloud water C. at t
REAL, DIMENSION(:,:,:),    INTENT(IN)    :: PCRT    ! Rain water C. at t
!
REAL, DIMENSION(:,:,:),    INTENT(INOUT) :: P_RR_CVRC
REAL, DIMENSION(:,:,:),    INTENT(INOUT) :: P_CR_CVRC
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3)) :: ZDR
!
LOGICAL, DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3)) :: ZMASKR, ZMASKC
!
REAL :: ZFACT
!
!
!
!-------------------------------------------------------------------------------
!
P_RR_CVRC(:,:,:) = 0.
P_CR_CVRC(:,:,:) = 0.
!
ZDR(:,:,:) = 9999.
ZMASKR(:,:,:) = PRRT(:,:,:).GT.XRTMIN(3) .AND. PCRT(:,:,:).GT.XCTMIN(3)
ZMASKC(:,:,:) = PRCT(:,:,:).GT.XRTMIN(2) .AND. PCCT(:,:,:).GT.XCTMIN(2)
WHERE(ZMASKR(:,:,:))
   ZDR(:,:,:)=(6.*PRRT(:,:,:)/XPI/XRHOLW/PCRT(:,:,:))**0.33
END WHERE
!
! Transfer all drops in droplets if out of cloud and Dr<82microns
!
WHERE( ZMASKR(:,:,:) .AND. .NOT.ZMASKC(:,:,:) .AND. ZDR(:,:,:).LT.82.E-6)
   P_RR_CVRC(:,:,:) = -PRRT(:,:,:)
   P_CR_CVRC(:,:,:) = -PCRT(:,:,:)
END WHERE
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPS_TO_DROPLETS_CONV
