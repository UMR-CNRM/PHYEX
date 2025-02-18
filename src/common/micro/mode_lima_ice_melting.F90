!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_LIMA_ICE_MELTING
  IMPLICIT NONE
CONTAINS
!     ########################################################################
  SUBROUTINE LIMA_ICE_MELTING (CST, LIMAP, KSIZE, PTSTEP, ODCOMPUTE,                 &
                               PEXNREF, PPABST,                          &
                               PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                               PCIT, PCIT_SHAPE, PINT,                               &
                               P_TH_IMLT, P_RC_IMLT, P_CC_IMLT, P_SHCI_IMLT,         &
                               PB_TH, PB_RC, PB_CC, PB_RI, PB_CI, PB_CI_SHAPE, PB_IFNN)
!     ########################################################################
!
!!    PURPOSE
!!    -------
!!      Melting of pristine ice crystals
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
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE MODD_CST, ONLY:CST_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_T),INTENT(IN)::CST
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
INTEGER,              INTENT(IN)    :: KSIZE
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRCT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRRT    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRST    ! 
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRGT    ! 
!
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PCIT    ! Rain water C. at t
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(IN)    :: PCIT_SHAPE    ! Pristine ice C. at t for different shapes
REAL, DIMENSION(KSIZE,LIMAP%NMOD_IFN), INTENT(IN)    :: PINT    ! Nucleated IFN C. at t
!
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_TH_IMLT
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_RC_IMLT
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: P_CC_IMLT
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: P_SHCI_IMLT
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_RC
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_CC
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_RI
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_CI
REAL, DIMENSION(KSIZE,LIMAP%NNB_CRYSTAL_SHAPE), INTENT(INOUT) :: PB_CI_SHAPE
REAL, DIMENSION(KSIZE,LIMAP%NMOD_IFN), INTENT(INOUT) :: PB_IFNN
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PTHT)) ::  &
     ZW,       &
     ZT,       &
     ZTCELSIUS,&
     ZLSFACT,  &
     ZLVFACT,  &
     ZMASK,    &
     ZCIT
!
INTEGER :: IMOD_IFN, ISH
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE_MELTING', 0, ZHOOK_HANDLE)
P_TH_IMLT(:) = 0.
P_RC_IMLT(:) = 0.
P_CC_IMLT(:) = 0.
!
! Temperature
ZT(:)  = PTHT(:) * ( PPABST(:)/CST%XP00 ) ** (CST%XRD/CST%XCPD)
ZTCELSIUS(:) = ZT(:)-CST%XTT
!
ZW(:) = PEXNREF(:)*( CST%XCPD+CST%XCPV*PRVT(:)+CST%XCL*(PRCT(:)+PRRT(:)) &
     +CST%XCI*(PRIT(:)+PRST(:)+PRGT(:)) )
ZLSFACT(:) = (CST%XLSTT+(CST%XCPV-CST%XCI)*ZTCELSIUS(:))/ZW(:)          ! L_s/(Pi_ref*C_ph)
ZLVFACT(:) = (CST%XLVTT+(CST%XCPV-CST%XCL)*ZTCELSIUS(:))/ZW(:)          ! L_v/(Pi_ref*C_ph)
!
ZW(:) = 0.0
!
ZMASK(:) = 0.
!
IF (LIMAP%LCRYSTAL_SHAPE) THEN
  ZCIT(:) = SUM(PCIT_SHAPE,DIM=2)
ELSE
  ZCIT(:) = PCIT(:)
END IF
WHERE( (ZT(:)>CST%XTT) .AND. (PRIT(:)>LIMAP%XRTMIN(4)) .AND. ODCOMPUTE(:) )
   P_TH_IMLT(:) = - PRIT(:)*(ZLSFACT(:)-ZLVFACT(:))
   P_RC_IMLT(:) = PRIT(:)
   P_CC_IMLT(:) = ZCIT(:)
   PB_TH(:) = PB_TH(:) + P_TH_IMLT(:)
   PB_RC(:) = PB_RC(:) + PRIT(:)
   PB_CC(:) = PB_CC(:) + ZCIT(:)
   PB_RI(:) = PB_RI(:) - PRIT(:)
   PB_CI(:) = PB_CI(:) - ZCIT(:)
   ZMASK(:) = 1.
ENDWHERE
!
IF (LIMAP%LCRYSTAL_SHAPE) THEN
  DO ISH = 1, LIMAP%NNB_CRYSTAL_SHAPE
    WHERE( (ZT(:)>CST%XTT) .AND. (PRIT(:)>LIMAP%XRTMIN(4)) .AND. ODCOMPUTE(:) )
      P_SHCI_IMLT(:,ISH) = -PCIT_SHAPE(:,ISH)
      PB_CI_SHAPE(:,ISH) = PB_CI_SHAPE(:,ISH) - PCIT_SHAPE(:,ISH)
    END WHERE
  END DO
END IF

DO IMOD_IFN = 1,LIMAP%NMOD_IFN
   PB_IFNN(:,IMOD_IFN) = PB_IFNN(:,IMOD_IFN) - PINT(:,IMOD_IFN)* ZMASK(:)
ENDDO
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_ICE_MELTING', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_ICE_MELTING
END MODULE MODE_LIMA_ICE_MELTING
