!      #################################
       MODULE MODI_LIMA_ICE_MELTING
!      #################################
!
INTERFACE
      SUBROUTINE LIMA_ICE_MELTING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,   &
                                   PEXNREF, PPABST,                          &
                                   PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                   PCIT, PINT,                               &
                                   P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,          &
                                   PB_TH, PB_RC, PB_CC, PB_RI, PB_CI, PB_IFNN)
!
REAL,                     INTENT(IN)    :: PTSTEP 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT
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
!     ######################################################################
      SUBROUTINE LIMA_ICE_MELTING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,   &
                                   PEXNREF, PPABST,                          &
                                   PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                   PCIT, PINT,                               &
                                   P_TH_IMLT, P_RC_IMLT, P_CC_IMLT,          &
                                   PB_TH, PB_RC, PB_CC, PB_RI, PB_CI, PB_IFNN)
!     ######################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the cold-phase homogeneous
!!    freezing of CCN, droplets and drops (T<-35°C)
!!
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy*   jan. 2014  add budgets
!!      B.Vie 10/2016 Bug zero division
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
REAL,                     INTENT(IN)    :: PTSTEP 
CHARACTER(LEN=*),         INTENT(IN)    :: HFMFILE
LOGICAL,                  INTENT(IN)    :: OCLOSE_OUT
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
WHERE( (ZT(:)<XTT-35.0) .AND. (PRIT(:)>XRTMIN(3)) .AND. LDCOMPUTE(:) )
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
! Correction BVIE aerosols not released but in droplets
!      ZIFS(:,JMOD_IFN) = ZIFS(:,JMOD_IFN) + ZINS(:,JMOD_IFN)*(1.-ZMASK(:)) 
   PB_IFNN(:,JMOD_IFN) = PB_IFNN(:,JMOD_IFN) - PINT(:,JMOD_IFN)* ZMASK(:)
ENDDO
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_ICE_MELTING
