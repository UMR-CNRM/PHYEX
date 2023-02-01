!MNH_LIC Copyright 2013-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      #################################
       MODULE MODI_LIMA_ICE_MELTING
!      #################################
!
INTERFACE
      SUBROUTINE LIMA_ICE_MELTING (PTSTEP, LDCOMPUTE,                        &
                                   PEXNREF, PPABST,                          &
                                   PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                   PCIT, PINT,                               &
                                   P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,          &
                                   PB_TH, PB_RC, PB_CC, PB_RI, PB_CI, PB_IFNN)
!
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Rain water C. at t
REAL, DIMENSION(:,:), INTENT(IN)    :: PINT    ! Nucleated IFN C. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_IMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_IMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_IMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_CC
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_CI
REAL, DIMENSION(:,:), INTENT(INOUT) :: PB_IFNN
!
END SUBROUTINE LIMA_ICE_MELTING
END INTERFACE
END MODULE MODI_LIMA_ICE_MELTING
!
!     ########################################################################
      SUBROUTINE LIMA_ICE_MELTING (PTSTEP, LDCOMPUTE,                        &
                                   PEXNREF, PPABST,                          &
                                   PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                   PCIT, PINT,                               &
                                   P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,          &
                                   PB_TH, PB_RC, PB_CC, PB_RI, PB_CI, PB_IFNN)
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
USE MODD_CST,             ONLY : XP00, XRD, XCPD, XCPV, XCL, XCI, XTT, XLSTT, XLVTT
USE MODD_PARAM_LIMA,      ONLY : XRTMIN, NMOD_IFN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                 INTENT(IN)    :: PTSTEP 
LOGICAL, DIMENSION(:),INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:),   INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:),   INTENT(IN)    :: PRVT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRCT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRRT    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:),   INTENT(IN)    :: PRST    ! 
REAL, DIMENSION(:),   INTENT(IN)    :: PRGT    ! 
!
REAL, DIMENSION(:),   INTENT(IN)    :: PCIT    ! Rain water C. at t
REAL, DIMENSION(:,:), INTENT(IN)    :: PINT    ! Nucleated IFN C. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: P_TH_IMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_RC_IMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: P_CC_IMLT
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_TH
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_RC
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_CC
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_RI
REAL, DIMENSION(:),   INTENT(INOUT) :: PB_CI
REAL, DIMENSION(:,:), INTENT(INOUT) :: PB_IFNN
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PTHT)) ::  &
     ZW,       &
     ZT,       &
     ZTCELSIUS,&
     ZLSFACT,  &
     ZLVFACT,  &
     ZMASK
!
INTEGER :: JMOD_IFN
!
!
!
!-------------------------------------------------------------------------------
!
P_TH_IMLT(:) = 0.
P_RC_IMLT(:) = 0.
P_CC_IMLT(:) = 0.
!
! Temperature
ZT(:)  = PTHT(:) * ( PPABST(:)/XP00 ) ** (XRD/XCPD)
ZTCELSIUS(:) = ZT(:)-XTT
!
ZW(:) = PEXNREF(:)*( XCPD+XCPV*PRVT(:)+XCL*(PRCT(:)+PRRT(:)) &
     +XCI*(PRIT(:)+PRST(:)+PRGT(:)) )
ZLSFACT(:) = (XLSTT+(XCPV-XCI)*ZTCELSIUS(:))/ZW(:)          ! L_s/(Pi_ref*C_ph)
ZLVFACT(:) = (XLVTT+(XCPV-XCL)*ZTCELSIUS(:))/ZW(:)          ! L_v/(Pi_ref*C_ph)
!
ZW(:) = 0.0
!
ZMASK(:) = 0.
!
WHERE( (ZT(:)>XTT) .AND. (PRIT(:)>XRTMIN(4)) .AND. LDCOMPUTE(:) )
   P_TH_IMLT(:) = - PRIT(:)*(ZLSFACT(:)-ZLVFACT(:))
   P_RC_IMLT(:) = PRIT(:)
   P_CC_IMLT(:) = PCIT(:)
   PB_TH(:) = PB_TH(:) + P_TH_IMLT(:)
   PB_RC(:) = PB_RC(:) + PRIT(:)
   PB_CC(:) = PB_CC(:) + PCIT(:)
   PB_RI(:) = PB_RI(:) - PRIT(:)
   PB_CI(:) = PB_CI(:) - PCIT(:)
   ZMASK(:) = 1.
ENDWHERE
!
DO JMOD_IFN = 1,NMOD_IFN
   PB_IFNN(:,JMOD_IFN) = PB_IFNN(:,JMOD_IFN) - PINT(:,JMOD_IFN)* ZMASK(:)
ENDDO
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ICE_MELTING
