!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      #####################
       MODULE MODI_LIMA_ICE_SNOW_DEPOSITION
!      #####################
!
INTERFACE
      SUBROUTINE LIMA_ICE_SNOW_DEPOSITION (PTSTEP, LDCOMPUTE,                 &
                                           PRHODREF, PSSI, PAI, PCJ, PLSFACT, &
                                           PRIT, PRST, PCIT, PLBDI, PLBDS,    &
                                           P_RI_CNVI, P_CI_CNVI,              &
                                           P_TH_DEPS, P_RS_DEPS,              &
                                           P_RI_CNVS, P_CI_CNVS,              &
                                           PA_TH, PA_RV, PA_RI, PA_CI, PA_RS  )
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
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI    ! Graupel m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS    ! Graupel m.r. at t 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DEPS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVS
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RV
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RS
!
END SUBROUTINE LIMA_ICE_SNOW_DEPOSITION
END INTERFACE
END MODULE MODI_LIMA_ICE_SNOW_DEPOSITION
!
!     ##########################################################################
SUBROUTINE LIMA_ICE_SNOW_DEPOSITION (PTSTEP, LDCOMPUTE,                        &
                                     PRHODREF,  PSSI, PAI, PCJ, PLSFACT,       &
                                     PRIT, PRST, PCIT, PLBDI, PLBDS,           &
                                     P_RI_CNVI, P_CI_CNVI,                     &
                                     P_TH_DEPS, P_RS_DEPS,                     &
                                     P_RI_CNVS, P_CI_CNVS,                     &
                                     PA_TH, PA_RV, PA_RI, PA_CI, PA_RS         )
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
                                 XAGGS_RLARGE1, XAGGS_RLARGE2  

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
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI    ! Graupel m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS    ! Graupel m.r. at t 
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVI
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_DEPS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RS_DEPS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RI_CNVS
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CI_CNVS
!
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RV
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_CI
REAL, DIMENSION(:),   INTENT(INOUT) :: PA_RS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMICRO ! Computations only where necessary
REAL,    DIMENSION(SIZE(PRHODREF)) :: ZZW, ZZW2, ZZX ! Work array
!
!
!-------------------------------------------------------------------------------
!
P_RI_CNVI(:) = 0.
P_CI_CNVI(:) = 0.
P_TH_DEPS(:) = 0.
P_RS_DEPS(:) = 0.
P_RI_CNVS(:) = 0.
P_CI_CNVS(:) = 0.
!
! Physical limitations
!
!
! Looking for regions where computations are necessary
!
GMICRO(:) = .FALSE.
GMICRO(:) = LDCOMPUTE(:) .AND. &
     (PRIT(:)>XRTMIN(4)      .OR.  &
      PRST(:)>XRTMIN(5))
!
!
WHERE( GMICRO )
!
!*       2.1    Conversion of snow to r_i: RSCNVI
!        ----------------------------------------
!
!
   ZZW2(:) = 0.0
   ZZW(:) = 0.0
   WHERE ( PLBDS(:)<XLBDASCNVI_MAX .AND. (PRST(:)>XRTMIN(5)) &
                                   .AND. (PSSI(:)<0.0)       )
      ZZW(:) = (PLBDS(:)*XDSCNVI_LIM)**(XALPHAS)
      ZZX(:) = ( -PSSI(:)/PAI(:) ) * (XCCS*PLBDS(:)**XCXS)/PRHODREF(:) * (ZZW(:)**XNUS) * EXP(-ZZW(:))
!
      ZZW(:) = ( XR0DEPSI+XR1DEPSI*PCJ(:) )*ZZX(:)
!
      ZZW2(:) = ZZW(:)*( XC0DEPSI+XC1DEPSI*PCJ(:) )/( XR0DEPSI+XR1DEPSI*PCJ(:) )
   END WHERE
!
   P_RI_CNVI(:) = ZZW(:)
   P_CI_CNVI(:) = ZZW2(:)
!
   PA_RI(:) = PA_RI(:) + P_RI_CNVI(:)
   PA_CI(:) = PA_CI(:) + P_CI_CNVI(:)
   PA_RS(:) = PA_RS(:) - P_RI_CNVI(:)
!
!
!*       2.2    Deposition of water vapor on r_s: RVDEPS
!        -----------------------------------------------
!
!
   ZZW(:) = 0.0
   WHERE ( (PRST(:)>XRTMIN(5)) )
      ZZW(:) = ( PSSI(:)/(PAI(:))/PRHODREF(:) ) * &
           ( X0DEPS*PLBDS(:)**XEX0DEPS + X1DEPS*PCJ(:)*PLBDS(:)**XEX1DEPS )
      ZZW(:) =    ZZW(:)*(0.5+SIGN(0.5,ZZW(:))) - ABS(ZZW(:))*(0.5-SIGN(0.5,ZZW(:)))
   END WHERE
!
   P_RS_DEPS(:) = ZZW(:)
   P_TH_DEPS(:) = P_RS_DEPS(:) * PLSFACT(:)
!
   PA_TH(:) = PA_TH(:) + P_TH_DEPS(:)
   PA_RV(:) = PA_RV(:) - P_RS_DEPS(:) 
   PA_RS(:) = PA_RS(:) + P_RS_DEPS(:) 
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
PA_RI(:) = PA_RI(:) + P_RI_CNVS(:)
PA_CI(:) = PA_CI(:) + P_CI_CNVS(:)
PA_RS(:) = PA_RS(:) - P_RI_CNVS(:)
!
!
END WHERE
!
!
END SUBROUTINE LIMA_ICE_SNOW_DEPOSITION
