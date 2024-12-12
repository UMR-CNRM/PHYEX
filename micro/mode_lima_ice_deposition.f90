!MNH_LIC Copyright 2018-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_ICE_DEPOSITION
  IMPLICIT NONE
CONTAINS
!     ##########################################################################
  SUBROUTINE LIMA_ICE_DEPOSITION (PTSTEP, LDCOMPUTE,                        &
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
!!      C. Barthe               01/2024   Add different ice crystal shapes
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, XCTMIN, XALPHAI, XNUI, &
                                 NMOM_I, NMOM_S, &
                                 LCRYSTAL_SHAPE, NNB_CRYSTAL_SHAPE 
USE MODD_PARAM_LIMA_COLD, ONLY : XDICNVS_LIM, XLBDAICNVS_LIM,         &
                                 XC0DEPIS, XC1DEPIS, XR0DEPIS, XR1DEPIS,      &
                                 XDI, X0DEPI, X2DEPI, &
                                 XDI_SHAPE,                      &
                                 X0DEPI_SHAPE, X2DEPI_SHAPE,     &
                                 XC0DEPIS_SHAPE, XC1DEPIS_SHAPE, &
                                 XR0DEPIS_SHAPE, XR1DEPIS_SHAPE
USE MODD_CST,             ONLY : XTT
!
USE MODE_LIMA_SHAPE_COMPUTE_LBDA, ONLY : LIMA_SHAPE_COMPUTE_LBDA
USE MODE_LIMA_CHANGE_SHAPE,       ONLY : LIMA_CHANGE_SHAPE 
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:),   INTENT(IN)    :: PT   ! abs. pressure at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PSSI ! sursat. / ice at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PAI  ! thermodyn. fct Ai
REAL, DIMENSION(:),   INTENT(IN)    :: PCJ  ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PLSFACT  ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Ice crystal C. at t
REAL, DIMENSION(:,:), INTENT(IN)    :: PCIT_SHAPE ! Ice crystal C. at t for each shape
!
REAL, DIMENSION(:),   INTENT(IN)    :: PLBDI    ! Ice crystal slope param. 
!
REAL, DIMENSION(:),   INTENT(OUT)   :: P_TH_DEPI
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_DEPI
REAL, DIMENSION(:),   INTENT(OUT)   :: P_RI_CNVS
REAL, DIMENSION(:),   INTENT(OUT)   :: P_CI_CNVS
REAL, DIMENSION(:,:), INTENT(OUT)   :: P_SHCI_HACH
REAL, DIMENSION(:,:), INTENT(OUT)   :: P_SHCI_CNVS
!
!*       0.2   Declarations of local variables :
!
LOGICAL, DIMENSION(SIZE(PRHODREF)) :: GMICRO ! Computations only where necessary
REAL,    DIMENSION(SIZE(PRHODREF)) :: ZZW, ZZW2, ZZX, ZCRIAUTI ! Work array
!
! for ice crystal shapes
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRIT_SHAPE
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLBDAI_SHAPE
INTEGER                           :: JSH
INTEGER                           :: ISIZE
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
IF (NMOM_I.EQ.1) THEN
   WHERE( GMICRO )
!
!*       Conversion of pristine ice to r_s: RICNVS
!        -----------------------------------------
!
      ZCRIAUTI(:)=MIN(0.2E-4,10**(0.06*(PT(:)-XTT)-3.5))
      ZZW(:) = 0.0
      WHERE ( (PRIT(:)>XRTMIN(4)))
         ZZW(:)   = 1.E-3 * EXP( 0.015*(PT(:)-XTT) ) * MAX( PRIT(:)-ZCRIAUTI(:),0.0 )
      END WHERE
!
      P_RI_CNVS(:) = - ZZW(:)
   END WHERE
ELSE
!++cb++
!
!*       Deposition of water vapor on r_i: RVDEPI
!        ----------------------------------------
!
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    ZZW(:) = 0.
    WHERE ( GMICRO .AND. (PRIT(:) > XRTMIN(4)) .AND. (PCIT(:) > XCTMIN(4)) )
      ZZW(:) = (PSSI(:) / PAI(:)) * PCIT(:) *        &
               (X0DEPI / PLBDI(:) + X2DEPI * PCJ(:) * PCJ(:) / PLBDI(:)**(XDI+2.0))
    END WHERE
    P_RI_DEPI(:) = ZZW(:)
  ELSE
    ALLOCATE(ZRIT_SHAPE(SIZE(PRHODREF),NNB_CRYSTAL_SHAPE))
    ALLOCATE(ZLBDAI_SHAPE(SIZE(PRHODREF),NNB_CRYSTAL_SHAPE))
    ! compute lambdai and ri per shape
    ISIZE = SIZE(PRHODREF)
    CALL LIMA_SHAPE_COMPUTE_LBDA (ISIZE, PRIT, PCIT_SHAPE, &
                                  ZRIT_SHAPE, ZLBDAI_SHAPE)
    !
    ! compute the deposition rate per shape
    DO JSH = 1, NNB_CRYSTAL_SHAPE
      ZZW(:) = 0.
      WHERE (GMICRO(:) .AND. (ZRIT_SHAPE(:,JSH) > XRTMIN(4)) .AND. (PCIT_SHAPE(:,JSH) > XCTMIN(4)))
        ZZW(:) = (PSSI(:) / PAI(:)) * PCIT_SHAPE(:,JSH) *   &
                 (X0DEPI_SHAPE(JSH) / ZLBDAI_SHAPE(:,JSH) + &
                  X2DEPI_SHAPE(JSH) * PCJ(:) * PCJ(:) /     &
                  ZLBDAI_SHAPE(:,JSH)**(XDI_SHAPE(JSH)+2.0))
      END WHERE
      P_RI_DEPI(:) = P_RI_DEPI(:) + ZZW(:)
    END DO
    !
    ! plates and columns can change shape
    CALL LIMA_CHANGE_SHAPE (PT, XCTMIN, P_RI_DEPI, PCIT_SHAPE, P_SHCI_HACH)
    !
    DEALLOCATE(ZRIT_SHAPE)
  END IF
!
!
!*       Conversion of pristine ice to r_s: RICNVS
!        -----------------------------------------
!
  IF (.NOT. LCRYSTAL_SHAPE) THEN
    ZZW(:) = 0.0
    ZZW2(:) = 0.0
    WHERE ( GMICRO .AND. (PLBDI(:) < XLBDAICNVS_LIM) .AND. (PCIT(:) > XCTMIN(4)) &
                                                     .AND. (PSSI(:) > 0.0)     )
      ZZW(:) = (PLBDI(:) * XDICNVS_LIM)**(XALPHAI)
      ZZX(:) = (PSSI(:) / PAI(:)) * PCIT(:) * (ZZW(:)**XNUI) * EXP(-ZZW(:))
      ZZW(:) = (XR0DEPIS + XR1DEPIS * PCJ(:)) * ZZX(:)                             
!
      ZZW2(:) = ZZW(:) * (XC0DEPIS + XC1DEPIS * PCJ(:)) / (XR0DEPIS + XR1DEPIS * PCJ(:))
    END WHERE
    P_RI_CNVS(:) = - ZZW(:)
    P_CI_CNVS(:) = - ZZW2(:)
  ELSE
    DO JSH = 1, NNB_CRYSTAL_SHAPE
      ZZW(:) = 0.0
      ZZW2(:) = 0.0
      WHERE ( GMICRO .AND. (ZLBDAI_SHAPE(:,JSH) < XLBDAICNVS_LIM) &
                     .AND. (PCIT_SHAPE(:,JSH) > XCTMIN(4))        &
                     .AND. (PSSI(:) > 0.0) )
        ZZW(:) = (ZLBDAI_SHAPE(:,JSH) * XDICNVS_LIM)**(XALPHAI)
        ZZX(:) = (PSSI(:) / PAI(:)) * PCIT_SHAPE(:,JSH) * (ZZW(:)**XNUI) * EXP(-ZZW(:))
        ZZW(:) = (XR0DEPIS_SHAPE(JSH) + XR1DEPIS_SHAPE(JSH) * PCJ(:)) * ZZX(:)                             
        !
        ZZW2(:) = ZZW(:) * (XC0DEPIS_SHAPE(JSH) + XC1DEPIS_SHAPE(JSH) * PCJ(:)) / &
                           (XR0DEPIS_SHAPE(JSH) + XR1DEPIS_SHAPE(JSH) * PCJ(:))
      END WHERE
      P_RI_CNVS(:)       = P_RI_CNVS(:) - ZZW(:)
      P_SHCI_CNVS(:,JSH) = - ZZW2(:)
    END DO
    DEALLOCATE(ZLBDAI_SHAPE)
  END IF

!   WHERE( GMICRO )
!
!*       Deposition of water vapor on r_i: RVDEPI
!        ----------------------------------------
!
!      ZZW(:) = 0.0
!      WHERE ( (PRIT(:)>XRTMIN(4)) .AND. (PCIT(:)>XCTMIN(4)) )
!         ZZW(:) = ( PSSI(:) / PAI(:) ) * PCIT(:) *        &
!              ( X0DEPI/PLBDI(:)+X2DEPI*PCJ(:)*PCJ(:)/PLBDI(:)**(XDI+2.0) )
!      END WHERE
!      P_RI_DEPI(:) = ZZW(:)
!
!*       Conversion of pristine ice to r_s: RICNVS
!        -----------------------------------------
!
!      ZZW(:) = 0.0
!      ZZW2(:) = 0.0
!      WHERE ( (PLBDI(:)<XLBDAICNVS_LIM) .AND. (PCIT(:)>XCTMIN(4)) &
!                                        .AND. (PSSI(:)>0.0)       )
!         ZZW(:) = (PLBDI(:)*XDICNVS_LIM)**(XALPHAI)
!         ZZX(:) = (PSSI(:)/PAI(:))*PCIT(:) * (ZZW(:)**XNUI) *EXP(-ZZW(:))
!
!         ZZW(:) = (XR0DEPIS + XR1DEPIS*PCJ(:))*ZZX(:)                             
!
!         ZZW2(:) = ZZW(:) * (XC0DEPIS+XC1DEPIS*PCJ(:)) / (XR0DEPIS+XR1DEPIS*PCJ(:))
!      END WHERE
!      P_RI_CNVS(:) = - ZZW(:)
!      P_CI_CNVS(:) = - ZZW2(:)
!   END WHERE
!--cb--
END IF
!
IF (NMOM_S.EQ.0) THEN
   P_RI_CNVS(:) = 0.
   P_CI_CNVS(:) = 0.
END IF
!
END SUBROUTINE LIMA_ICE_DEPOSITION
END MODULE MODE_LIMA_ICE_DEPOSITION
