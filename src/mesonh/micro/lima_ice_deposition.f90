!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!      #####################
       MODULE MODI_LIMA_ICE_DEPOSITION
!      #####################
!
INTERFACE
      SUBROUTINE LIMA_ICE_DEPOSITION (PTSTEP, LDCOMPUTE,                 &
                                      PRHODREF, PSSI, PAI, PCJ, PLSFACT, &
                                      PRIT, PCIT, PLBDI,                 &
                                      P_TH_DEPI, P_RI_DEPI,              &
                                      P_RI_CNVS, P_CI_CNVS               )
!
REAL,                 INTENT(IN)    :: PTSTEP
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PAI  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT  ! abs. pressure at time t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI    ! Graupel m.r. at t 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_DEPI
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_DEPI
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_CNVS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_CNVS
!
END SUBROUTINE LIMA_ICE_DEPOSITION
END INTERFACE
END MODULE MODI_LIMA_ICE_DEPOSITION
!
!     ##########################################################################
SUBROUTINE LIMA_ICE_DEPOSITION (PTSTEP, LDCOMPUTE,                        &
                                PRHODREF,  PSSI, PAI, PCJ, PLSFACT,       &
                                PRIT, PCIT, PLBDI,                        &
                                P_TH_DEPI, P_RI_DEPI,                     &
                                P_RI_CNVS, P_CI_CNVS                      )
!     ##########################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources
!!    for slow cold processes :
!!      - conversion of snow to ice
!!      - deposition of vapor on snow
!!      - conversion of ice to snow (Harrington 1995)
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
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XALPHAI, XALPHAS, XNUI, XNUS 
USE MODD_PARAM_LIMA_COLD, ONLY : XCXS, XCCS, &
                                 XLBDAS_MAX, XDSCNVI_LIM, XLBDASCNVI_MAX,     &
                                 XC0DEPSI, XC1DEPSI, XR0DEPSI, XR1DEPSI,      &
                                 XSCFAC, X1DEPS, X0DEPS, XEX1DEPS, XEX0DEPS,  &
                                 XDICNVS_LIM, XLBDAICNVS_LIM,                 &
                                 XC0DEPIS, XC1DEPIS, XR0DEPIS, XR1DEPIS,      &
                                 XCOLEXIS, XAGGS_CLARGE1, XAGGS_CLARGE2,      &
                                 XAGGS_RLARGE1, XAGGS_RLARGE2,                &
                                 XDI, X0DEPI, X2DEPI

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PAI  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ  ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT  ! abs. pressure at time t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI    ! Graupel m.r. at t 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_DEPI
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_DEPI
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_CNVS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_CNVS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMICRO ! Computations only where necessary
REAL,    DIMENSION(SIZE(PRHODREF)) :: ZZW, ZZW2, ZZX ! Work array
!
!
!-------------------------------------------------------------------------------
!
P_TH_DEPI(:) = 0.
P_RI_DEPI(:) = 0.
P_RI_CNVS(:) = 0.
P_CI_CNVS(:) = 0.
!
! Physical limitations
!
!
! Looking for regions where computations are necessary
!
GMICRO(:) = LDCOMPUTE(:) .AND. PRIT(:)>XRTMIN(4)
!
!
WHERE( GMICRO )
!
!
!*       2.2    Deposition of water vapor on r_i: RVDEPI
!        -----------------------------------------------
!
!
   ZZW(:) = 0.0
   WHERE ( (PRIT(:)>XRTMIN(4)) .AND. (PCIT(:)>XCTMIN(4)) )
      ZZW(:) = ( PSSI(:) / PAI(:) ) * PCIT(:) *        &
        ( X0DEPI/PLBDI(:)+X2DEPI*PCJ(:)*PCJ(:)/PLBDI(:)**(XDI+2.0) )
   END WHERE
!
   P_RI_DEPI(:) = ZZW(:)
!!$   P_TH_DEPI(:) = P_RI_DEPI(:) * PLSFACT(:)
!
!!$   PA_TH(:) = PA_TH(:) + P_TH_DEPI(:)
!!$   PA_RV(:) = PA_RV(:) - P_RI_DEPI(:) 
!!$   PA_RI(:) = PA_RI(:) + P_RI_DEPI(:) 
!
!
!*       2.3    Conversion of pristine ice to r_s: RICNVS
!        ------------------------------------------------
!
!
   ZZW(:) = 0.0
   ZZW2(:) = 0.0
   WHERE ( (PLBDI(:)<XLBDAICNVS_LIM) .AND. (PCIT(:)>XCTMIN(4)) &
                                     .AND. (PSSI(:)>0.0)       )
      ZZW(:) = (PLBDI(:)*XDICNVS_LIM)**(XALPHAI)
      ZZX(:) = ( PSSI(:)/PAI(:) )*PCIT(:) * (ZZW(:)**XNUI) *EXP(-ZZW(:))
!
      ZZW(:) = ( XR0DEPIS + XR1DEPIS*PCJ(:) )*ZZX(:)                             
!
      ZZW2(:) = ZZW(:) * (XC0DEPIS+XC1DEPIS*PCJ(:)) / (XR0DEPIS+XR1DEPIS*PCJ(:))
   END WHERE
!
P_RI_CNVS(:) = - ZZW(:)
P_CI_CNVS(:) = - ZZW2(:)
!
!
END WHERE
!
!
END SUBROUTINE LIMA_ICE_DEPOSITION
