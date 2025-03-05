!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_ICE_DEPOSITION
  IMPLICIT NONE
CONTAINS
!     ##########################################################################
  SUBROUTINE LIMA_ICE_DEPOSITION (CST, LIMAP, LIMAC, KSIZE, PTSTEP, ODCOMPUTE, &
                                  PRHODREF, PT,  PSSI, PAI, PCJ, PLSFACT,   &
                                  PRIT, PCIT, PCIT_SHAPE, PLBDI,            &
                                  P_TH_DEPI, P_RI_DEPI, P_SHCI_HACH,        &
                                  P_RI_CNVS, P_CI_CNVS, P_SHCI_CNVS         )
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
USE MODD_PARAM_LIMA_COLD, ONLY:PARAM_LIMA_COLD_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
USE MODE_LIMA_SHAPE_COMPUTE_LBDA, ONLY : LIMA_SHAPE_COMPUTE_LBDA
USE MODE_LIMA_CHANGE_SHAPE, ONLY : LIMA_CHANGE_SHAPE
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_COLD_T),INTENT(IN)::LIMAC
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(CST_T),INTENT(IN)::CST
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
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN) :: PCIT_SHAPE ! Ice crystal C. at t for each shape
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLBDI    ! Graupel m.r. at t 
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_DEPI
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_DEPI
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE),   INTENT(OUT)   :: P_SHCI_HACH
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RI_CNVS
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CI_CNVS
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE),   INTENT(OUT)   :: P_SHCI_CNVS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMICRO ! Computations only where necessary
REAL,    DIMENSION(SIZE(PRHODREF)) :: ZZW, ZZW2, ZZX, ZCRIAUTI ! Work array
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRIT_SHAPE
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLBDAI_SHAPE
INTEGER                           :: ISH
INTEGER                           :: ISIZE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE_DEPOSITION', 0, ZHOOK_HANDLE)
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
IF (LIMAP%NMOM_I.GE.2) GMICRO(:) = GMICRO(:) .AND. PCIT(:)>LIMAP%XCTMIN(4)
!
ZZW(:) = 0.0
ZZW2(:) = 0.0
IF (LIMAP%NMOM_I.EQ.1) THEN
!
!*       Conversion of pristine ice to r_s: RICNVS
!        -----------------------------------------
!
   WHERE( GMICRO )
      ZCRIAUTI(:)=MIN(0.2E-4,10**(0.06*(PT(:)-CST%XTT)-3.5))
      ZZW(:)   = 1.E-3 * EXP( 0.015*(PT(:)-CST%XTT) ) * MAX( PRIT(:)-ZCRIAUTI(:),0.0 )
      P_RI_CNVS(:) = - ZZW(:)
   END WHERE
ELSE
   IF (.NOT. LIMAP%LCRYSTAL_SHAPE) THEN
!
!*       Deposition of water vapor on r_i: RVDEPI
!        ----------------------------------------
!
      WHERE( GMICRO )
         ZZW(:) = ( PSSI(:) / PAI(:) ) * PCIT(:) *        &
              ( LIMAC%X0DEPI/PLBDI(:)+LIMAC%X2DEPI*PCJ(:)*PCJ(:)/PLBDI(:)**(LIMAC%XDI+2.0) )
         P_RI_DEPI(:) = ZZW(:)
      END WHERE
!
!*       Conversion of pristine ice to r_s: RICNVS
!        -----------------------------------------
!
      ZZW(:) = 0.0
      ZZW2(:) = 0.0
      WHERE ( GMICRO .AND. PLBDI(:)<LIMAC%XLBDAICNVS_LIM .AND. PSSI(:)>0.0 )
         ZZW(:) = (PLBDI(:)*LIMAC%XDICNVS_LIM)**(LIMAP%XALPHAI)
         ZZX(:) = (PSSI(:)/PAI(:))*PCIT(:) * (ZZW(:)**LIMAP%XNUI) *EXP(-ZZW(:))
!
         ZZW(:) = (LIMAC%XR0DEPIS + LIMAC%XR1DEPIS*PCJ(:))*ZZX(:)                             
!
         ZZW2(:) = ZZW(:) * (LIMAC%XC0DEPIS+LIMAC%XC1DEPIS*PCJ(:)) / (LIMAC%XR0DEPIS+LIMAC%XR1DEPIS*PCJ(:))
         P_RI_CNVS(:) = - ZZW(:)
         P_CI_CNVS(:) = - ZZW2(:)
      END WHERE
   ELSE ! LCRYSTAL_SHAPE
!
!*       Deposition of water vapor on r_i: RVDEPI
!        ----------------------------------------
      ALLOCATE(ZRIT_SHAPE(SIZE(PRHODREF),LIMAP%NNB_CRYSTAL_SHAPE))
      ALLOCATE(ZLBDAI_SHAPE(SIZE(PRHODREF),LIMAP%NNB_CRYSTAL_SHAPE))
      ! compute lambdai and ri per shape
      CALL LIMA_SHAPE_COMPUTE_LBDA (LIMAP, LIMAC, KSIZE, PRIT, PCIT_SHAPE, &
           ZRIT_SHAPE, ZLBDAI_SHAPE)
      !
      ! compute the deposition rate per shape
      DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
         ZZW(:) = 0.
         WHERE ( GMICRO(:) .AND. (ZRIT_SHAPE(:,ISH)>LIMAP%XRTMIN(4)) .AND. (PCIT_SHAPE(:,ISH)>LIMAP%XCTMIN(4)) )
            ZZW(:) = (PSSI(:) / PAI(:)) * PCIT_SHAPE(:,ISH) *   &
                 (LIMAC%X0DEPI_SHAPE(ISH) / ZLBDAI_SHAPE(:,ISH) + &
                 LIMAC%X2DEPI_SHAPE(ISH) * PCJ(:) * PCJ(:) /     &
                 ZLBDAI_SHAPE(:,ISH)**(LIMAC%XDI_SHAPE(ISH)+2.0))
         END WHERE
         P_RI_DEPI(:) = P_RI_DEPI(:) + ZZW(:)
      END DO
      !
      ! plates and columns can change shape
      CALL LIMA_CHANGE_SHAPE (CST, LIMAP, KSIZE, PT, LIMAP%XCTMIN, P_RI_DEPI, PCIT_SHAPE, P_SHCI_HACH)
      !
      DEALLOCATE(ZRIT_SHAPE)
!
!*       Conversion of pristine ice to r_s: RICNVS
!        -----------------------------------------
      DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
         ZZW(:) = 0.0
         ZZW2(:) = 0.0
         WHERE ( GMICRO .AND. (ZLBDAI_SHAPE(:,ISH) < LIMAC%XLBDAICNVS_LIM) &
              .AND. (PCIT_SHAPE(:,ISH) > LIMAP%XCTMIN(4))        &
              .AND. (PSSI(:) > 0.0) )
            ZZW(:) = (ZLBDAI_SHAPE(:,ISH) * LIMAC%XDICNVS_LIM)**(LIMAP%XALPHAI)
            ZZX(:) = (PSSI(:) / PAI(:)) * PCIT_SHAPE(:,ISH) * (ZZW(:)**LIMAP%XNUI) * EXP(-ZZW(:))
            ZZW(:) = (LIMAC%XR0DEPIS_SHAPE(ISH) + LIMAC%XR1DEPIS_SHAPE(ISH) * PCJ(:)) * ZZX(:)                             
            !
            ZZW2(:) = ZZW(:) * (LIMAC%XC0DEPIS_SHAPE(ISH) + LIMAC%XC1DEPIS_SHAPE(ISH) * PCJ(:)) / &
                 (LIMAC%XR0DEPIS_SHAPE(ISH) + LIMAC%XR1DEPIS_SHAPE(ISH) * PCJ(:))
         END WHERE
         P_RI_CNVS(:)       = P_RI_CNVS(:) - ZZW(:)
         P_SHCI_CNVS(:,ISH) = - ZZW2(:)
      END DO
      DEALLOCATE(ZLBDAI_SHAPE)
   ENDIF ! LCRYSTAL_SHAPE
END IF
!
IF (LIMAP%NMOM_S.EQ.0) THEN
   P_RI_CNVS(:) = 0.
   P_CI_CNVS(:) = 0.
END IF
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE_DEPOSITION', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_ICE_DEPOSITION
END MODULE MODE_LIMA_ICE_DEPOSITION
