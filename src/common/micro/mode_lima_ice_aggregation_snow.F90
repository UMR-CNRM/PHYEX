!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_ICE_AGGREGATION_SNOW
  IMPLICIT NONE
CONTAINS
!     #######################################################################
  SUBROUTINE LIMA_ICE_AGGREGATION_SNOW (CST, LIMAP, LIMAC, KSIZE, ODCOMPUTE,                     &
                                        PT, PRHODREF,                         &
                                        PRIT, PRST, PCIT, PCST, PLBDI, PLBDS, &
                                        PCIT_SHAPE, PRIT_SHAPE, PLBDAI_SHAPE, &
                                        PLATHAM_IAGGS,                        &
                                        P_RI_AGGS, P_CI_AGGS,                 &
                                        P_SHCI_AGGS                           )
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
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
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
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN) :: PCIT_SHAPE
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN) :: PRIT_SHAPE
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN) :: PLBDAI_SHAPE
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLATHAM_IAGGS
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_AGGS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_AGGS
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(OUT):: P_SHCI_AGGS
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRIT)) :: ZZW1, ZZW2, ZZW3 ! work arrays
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: ISH
!
!-------------------------------------------------------------------------------
!
!
!*       2.4    Aggregation of r_i on r_s: CIAGGS and RIAGGS
!        ---------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE_AGGREGATION_SNOW', 0, ZHOOK_HANDLE)
ZZW1(:) = 0.
ZZW2(:) = 0.
ZZW3(:) = 0.
!
P_RI_AGGS(:) = 0.
P_CI_AGGS(:) = 0.
!
!
IF (LIMAP%NMOM_I.EQ.1) THEN
   WHERE ( PRIT(:)>LIMAP%XRTMIN(4) .AND. PRST(:)>LIMAP%XRTMIN(5) .AND. ODCOMPUTE(:) )
      ZZW1(:) = LIMAC%XFIAGGS * EXP( LIMAC%XCOLEXIS*(PT(:)-CST%XTT) ) &
                        * PLATHAM_IAGGS(:)            &
                        * PRIT(:)                     &
                        * PCST(:) * (1+(LIMAC%XFVELOS/PLBDS(:))**LIMAP%XALPHAS)**(-LIMAP%XNUS+LIMAC%XEXIAGGS/LIMAP%XALPHAS) &
                        * PRHODREF(:)**(-LIMAP%XCEXVT+1.) &
                        * PLBDS(:)**LIMAC%XEXIAGGS
!
      P_RI_AGGS(:) = - ZZW1(:)
   END WHERE
ELSE
   IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
      WHERE ( PRIT(:)>LIMAP%XRTMIN(4) .AND. PRST(:)>LIMAP%XRTMIN(5) .AND. &
           PCIT(:)>LIMAP%XCTMIN(4) .AND. PCST(:)>LIMAP%XCTMIN(5) .AND. ODCOMPUTE(:) )
         ZZW1(:) = (PLBDI(:) / PLBDS(:))**3
         ZZW2(:) = PCIT(:)*PCST(:)*EXP(LIMAC%XCOLEXIS*(PT(:)-CST%XTT))*PRHODREF(:) / (PLBDI(:)**3)
         ZZW3(:) = ZZW2(:)*(LIMAC%XAGGS_CLARGE1+LIMAC%XAGGS_CLARGE2*ZZW1(:))
!
         P_CI_AGGS(:) = - ZZW3(:)
!
         ZZW2(:) = ZZW2(:) / PLBDI(:)**LIMAC%XBI
         ZZW2(:) = ZZW2(:)*(LIMAC%XAGGS_RLARGE1+LIMAC%XAGGS_RLARGE2*ZZW1(:))
!
         P_RI_AGGS(:) = - ZZW2(:)
      END WHERE
   ELSE
      DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
         WHERE ( PRIT(:) > LIMAP%XRTMIN(4) .AND. PRST(:) > LIMAP%XRTMIN(5) .AND. PCST(:) > LIMAP%XCTMIN(5) .AND. &
              PCIT_SHAPE(:,ISH) > 0. .AND. PRIT_SHAPE(:,ISH) > 0. .AND. &
              ODCOMPUTE(:) )
            ZZW1(:) = (PLBDAI_SHAPE(:,ISH) / PLBDS(:))**3
            ZZW2(:) = PCIT_SHAPE(:,ISH) * PCST(:) * EXP(LIMAC%XCOLEXIS*(PT(:)-CST%XTT)) * PRHODREF(:) / &
                 (PLBDAI_SHAPE(:,ISH)**3)
            ZZW3(:) = ZZW2(:) * (LIMAC%XAGGS_CLARGE1 + LIMAC%XAGGS_CLARGE2 * ZZW1(:))
!
            P_SHCI_AGGS(:,ISH) = - ZZW3(:)
!
            ZZW2(:) = ZZW2(:) / PLBDAI_SHAPE(:,ISH)**LIMAC%XBI_SHAPE(ISH)
            ZZW2(:) = ZZW2(:) * (LIMAC%XAGGS_RLARGE1 + LIMAC%XAGGS_RLARGE2 * ZZW1(:))
!
            P_RI_AGGS(:) = P_RI_AGGS(:) - ZZW2(:)
         END WHERE
      END DO
   END IF
END IF
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE_AGGREGATION_SNOW', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_ICE_AGGREGATION_SNOW
END MODULE MODE_LIMA_ICE_AGGREGATION_SNOW
