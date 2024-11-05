!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_ICE_DEPOSITION
  IMPLICIT NONE
CONTAINS
!     ##########################################################################
  SUBROUTINE LIMA_ICE_DEPOSITION (LIMAP, LIMAC, KSIZE, PTSTEP, ODCOMPUTE,                 &
                                  PRHODREF, PT,  PSSI, PAI, PCJ, PLSFACT,   &
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
!!      B. Vié               30/08/2021      Disable CNVS if LSNOW=F  
!!      B. Vie                  03/2022   Add option for 1-moment pristine ice
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY : XTT
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_t
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,              INTENT(IN)    :: KSIZE
REAL,                 INTENT(IN)    :: PTSTEP
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT  ! abs. pressure at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PSSI  ! abs. pressure at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PAI  ! abs. pressure at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCJ  ! abs. pressure at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT  ! abs. pressure at time t
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDI    ! Graupel m.r. at t 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_DEPI
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_DEPI
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_CNVS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_CNVS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMICRO ! Computations only where necessary
TYPE(PARAM_LIMA_COLD_t),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_t),INTENT(IN)::LIMAP
REAL,    DIMENSION(SIZE(PRHODREF)) :: ZZW, ZZW2, ZZX, ZCRIAUTI ! Work array
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
GMICRO(:) = ODCOMPUTE(:) .AND. PRIT(:)>LIMAP%XRTMIN(4)
!
!
IF (LIMAP%NMOM_I.EQ.1) THEN
   WHERE( GMICRO )
!
!*       Conversion of pristine ice to r_s: RICNVS
!        -----------------------------------------
!
      ZCRIAUTI(:)=MIN(0.2E-4,10**(0.06*(PT(:)-XTT)-3.5))
      ZZW(:) = 0.0
      WHERE ( (PRIT(:)>LIMAP%XRTMIN(4)))
         ZZW(:)   = 1.E-3 * EXP( 0.015*(PT(:)-XTT) ) * MAX( PRIT(:)-ZCRIAUTI(:),0.0 )
      END WHERE
!
      P_RI_CNVS(:) = - ZZW(:)
   END WHERE
ELSE
   WHERE( GMICRO )
!
!*       Deposition of water vapor on r_i: RVDEPI
!        ----------------------------------------
!
      ZZW(:) = 0.0
      WHERE ( (PRIT(:)>LIMAP%XRTMIN(4)) .AND. (PCIT(:)>LIMAP%XCTMIN(4)) )
         ZZW(:) = ( PSSI(:) / PAI(:) ) * PCIT(:) *        &
              ( LIMAC%X0DEPI/PLBDI(:)+LIMAC%X2DEPI*PCJ(:)*PCJ(:)/PLBDI(:)**(LIMAC%XDI+2.0) )
      END WHERE
      P_RI_DEPI(:) = ZZW(:)
!
!*       Conversion of pristine ice to r_s: RICNVS
!        -----------------------------------------
!
      ZZW(:) = 0.0
      ZZW2(:) = 0.0
      WHERE ( (PLBDI(:)<LIMAC%XLBDAICNVS_LIM) .AND. (PCIT(:)>LIMAP%XCTMIN(4)) &
                                        .AND. (PSSI(:)>0.0)       )
         ZZW(:) = (PLBDI(:)*LIMAC%XDICNVS_LIM)**(LIMAP%XALPHAI)
         ZZX(:) = (PSSI(:)/PAI(:))*PCIT(:) * (ZZW(:)**LIMAP%XNUI) *EXP(-ZZW(:))
!
         ZZW(:) = (LIMAC%XR0DEPIS + LIMAC%XR1DEPIS*PCJ(:))*ZZX(:)                             
!
         ZZW2(:) = ZZW(:) * (LIMAC%XC0DEPIS+LIMAC%XC1DEPIS*PCJ(:)) / (LIMAC%XR0DEPIS+LIMAC%XR1DEPIS*PCJ(:))
      END WHERE
      P_RI_CNVS(:) = - ZZW(:)
      P_CI_CNVS(:) = - ZZW2(:)
   END WHERE
END IF
!
IF (LIMAP%NMOM_S.EQ.0) THEN
   P_RI_CNVS(:) = 0.
   P_CI_CNVS(:) = 0.
END IF
!
END SUBROUTINE LIMA_ICE_DEPOSITION
END MODULE MODE_LIMA_ICE_DEPOSITION
