!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_ICE_AGGREGATION_SNOW
  IMPLICIT NONE
CONTAINS
! #############################################################################
  SUBROUTINE LIMA_ICE_AGGREGATION_SNOW (LDCOMPUTE,                            &
                                        PT, PRHODREF,                         &
                                        PRIT, PRST, PCIT, PCST, PLBDI, PLBDS, &
                                        PLATHAM_IAGGS,                        &
                                        P_RI_AGGS, P_CI_AGGS,                 &
                                        PCIT_SHAPE, PRIT_SHAPE, PLBDAI_SHAPE, &
                                        P_SHCI_AGGS                           )
! #############################################################################
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
!  C. Barthe      06/2023: add Latham effect (Efield) for IAGGS
!  C. Barthe      01/2024: add several shapes for ice crystals
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY : XTT
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XCEXVT, NMOM_I, XNUS, XALPHAS, XCEXVT, &
                                 LCRYSTAL_SHAPE, NNB_CRYSTAL_SHAPE
USE MODD_PARAM_LIMA_COLD, ONLY : XBI, XCOLEXIS, XAGGS_CLARGE1, XAGGS_CLARGE2, &
                                 XAGGS_RLARGE1, XAGGS_RLARGE2, XFIAGGS, XFVELOS, XEXIAGGS, &
                                 XBI_SHAPE
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PT
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT
REAL, DIMENSION(:),   INTENT(IN)    :: PRST
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT
REAL, DIMENSION(:),   INTENT(IN)    :: PCST
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI 
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDS 
REAL, DIMENSION(:),   INTENT(IN)    :: PLATHAM_IAGGS
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_AGGS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_AGGS
!
REAL, DIMENSION(:,:), INTENT(IN)    :: PCIT_SHAPE
REAL, DIMENSION(:,:), INTENT(IN)    :: PRIT_SHAPE
REAL, DIMENSION(:,:), INTENT(IN)    :: PLBDAI_SHAPE
REAL, DIMENSION(:,:), INTENT(OUT)   :: P_SHCI_AGGS
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PRIT)) :: ZZW1, ZZW2, ZZW3 ! work arrays
INTEGER :: JSH
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
   WHERE ( PRIT(:)>XRTMIN(4) .AND. PRST(:)>XRTMIN(5) .AND. LDCOMPUTE(:) )
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
!++cb++ 31/01/24
   IF (.NOT. LCRYSTAL_SHAPE) THEN
     WHERE ( PRIT(:)>XRTMIN(4) .AND. PRST(:)>XRTMIN(5) .AND. &
             PCIT(:)>XCTMIN(4) .AND. PCST(:)>XCTMIN(5) .AND. LDCOMPUTE(:) )
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
   ELSE
     DO JSH = 1, NNB_CRYSTAL_SHAPE
       WHERE ( PRIT(:) > XRTMIN(4) .AND. PRST(:) > XRTMIN(5) .AND. PCST(:) > XCTMIN(5) .AND. &
               PCIT_SHAPE(:,JSH) > 0. .AND. PRIT_SHAPE(:,JSH) > 0. .AND. &
               LDCOMPUTE(:) )
          ZZW1(:) = (PLBDAI_SHAPE(:,JSH) / PLBDS(:))**3
          ZZW2(:) = PCIT_SHAPE(:,JSH) * PCST(:) * EXP(XCOLEXIS*(PT(:)-XTT)) * PRHODREF(:) / &
                   (PLBDAI_SHAPE(:,JSH)**3)
          ZZW3(:) = ZZW2(:) * (XAGGS_CLARGE1 + XAGGS_CLARGE2 * ZZW1(:))
!
          P_SHCI_AGGS(:,JSH) = - ZZW3(:)
!
          ZZW2(:) = ZZW2(:) / PLBDAI_SHAPE(:,JSH)**XBI_SHAPE(JSH)
          ZZW2(:) = ZZW2(:) * (XAGGS_RLARGE1 + XAGGS_RLARGE2 * ZZW1(:))
!
!++cb++ 21/02/24 different de ce qui est fait pour une seule forme car on ajoute les contributions pour chaque forme
          P_RI_AGGS(:) = P_RI_AGGS(:) - ZZW2(:)
       END WHERE
     END DO
   END IF
!--cb--
END IF
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ICE_AGGREGATION_SNOW
END MODULE MODE_LIMA_ICE_AGGREGATION_SNOW
