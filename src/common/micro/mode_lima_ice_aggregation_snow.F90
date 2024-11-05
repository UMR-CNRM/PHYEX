!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_ICE_AGGREGATION_SNOW
  IMPLICIT NONE
CONTAINS
!     #######################################################################
  SUBROUTINE LIMA_ICE_AGGREGATION_SNOW (LIMAP, LIMAC, KSIZE, ODCOMPUTE,                     &
                                        PT, PRHODREF,                         &
                                        PRIT, PRST, PCIT, PCST, PLBDI, PLBDS, &
                                        PLATHAM_IAGGS,                        &
                                        P_RI_AGGS, P_CI_AGGS                  )
!     #######################################################################
!
!!    PURPOSE
!!    -------
!!      Compute the aggregation of pristine ice on snow/aggregates
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
!  J. Wurtz       03/2022: new snow characteristics
!  B. Vie         03/2022: Add option for 1-moment pristine ice
!  M. Taufour     07/2022: add concentration for snow, graupel, hail        
!  C. Barthe       06/2023: add Latham effect (Efield) for IAGGS
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY : XTT
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XCEXVT, NMOM_I, XNUS, XALPHAS, XCEXVT
USE MODD_PARAM_LIMA_COLD, ONLY : XBI, XCOLEXIS, XAGGS_CLARGE1, XAGGS_CLARGE2, &
                                 XAGGS_RLARGE1, XAGGS_RLARGE2, XFIAGGS, XFVELOS, XEXIAGGS
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_t
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCST
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDI 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDS 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLATHAM_IAGGS
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_AGGS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_AGGS
!
!*       0.2   Declarations of local variables :
!
TYPE(PARAM_LIMA_COLD_t),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_t),INTENT(IN)::LIMAP
REAL, DIMENSION(SIZE(PRIT)) :: ZZW1, ZZW2, ZZW3 ! work arrays
!
!-------------------------------------------------------------------------------
!
!
!*       2.4    Aggregation of r_i on r_s: CIAGGS and RIAGGS
!        ---------------------------------------------------
!
ZZW1(:) = 0.
ZZW2(:) = 0.
ZZW3(:) = 0.
!
P_RI_AGGS(:) = 0.
P_CI_AGGS(:) = 0.
!
!
IF (NMOM_I.EQ.1) THEN
   WHERE ( PRIT(:)>XRTMIN(4) .AND. PRST(:)>XRTMIN(5) .AND. ODCOMPUTE(:) )
      ZZW1(:) = XFIAGGS * EXP( XCOLEXIS*(PT(:)-XTT) ) &
                        * PLATHAM_IAGGS(:)            &
                        * PRIT(:)                     &
                        * PCST(:) * (1+(XFVELOS/PLBDS(:))**XALPHAS)**(-XNUS+XEXIAGGS/XALPHAS) &
                        * PRHODREF(:)**(-XCEXVT+1.) &
                        * PLBDS(:)**XEXIAGGS
!
      P_RI_AGGS(:) = - ZZW1(:)
   END WHERE
ELSE
   WHERE ( PRIT(:)>XRTMIN(4) .AND. PRST(:)>XRTMIN(5) .AND. &
           PCIT(:)>XCTMIN(4) .AND. PCST(:)>XCTMIN(5) .AND. ODCOMPUTE(:) )
      ZZW1(:) = (PLBDI(:) / PLBDS(:))**3
      ZZW2(:) = PCIT(:)*PCST(:)*EXP(XCOLEXIS*(PT(:)-XTT))*PRHODREF(:) / (PLBDI(:)**3)
      ZZW3(:) = ZZW2(:)*(XAGGS_CLARGE1+XAGGS_CLARGE2*ZZW1(:))
!
      P_CI_AGGS(:) = - ZZW3(:)
!
      ZZW2(:) = ZZW2(:) / PLBDI(:)**XBI
      ZZW2(:) = ZZW2(:)*(XAGGS_RLARGE1+XAGGS_RLARGE2*ZZW1(:))
!
      P_RI_AGGS(:) = - ZZW2(:)
   END WHERE
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ICE_AGGREGATION_SNOW
END MODULE MODE_LIMA_ICE_AGGREGATION_SNOW
