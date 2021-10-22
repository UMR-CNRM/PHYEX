!      #################################
       MODULE MODI_LIMA_DROPS_HOM_FREEZING
!      #################################
!
INTERFACE
      SUBROUTINE LIMA_DROPS_HOM_FREEZING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE, &
                                        PEXNREF, PPABST,                          &
                                        PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                        PCRT,                                     &
                                        P_TH_HONR, P_RR_HONR, P_CR_HONR,          &
                                        PB_TH, PB_RR, PB_CR, PB_RG                )
!
REAL,                  INTENT(IN)    :: PTSTEP 
CHARACTER(LEN=*),      INTENT(IN)    :: HFMFILE
LOGICAL,               INTENT(IN)    :: OCLOSE_OUT
LOGICAL, DIMENSION(:), INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),    INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:),    INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:),    INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:),    INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:),    INTENT(IN)    :: PCRT    ! Rain water C. at t
!
REAL, DIMENSION(:),    INTENT(INOUT) :: P_TH_HONR
REAL, DIMENSION(:),    INTENT(INOUT) :: P_RR_HONR
REAL, DIMENSION(:),    INTENT(INOUT) :: P_CR_HONR
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_TH
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_RR
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_CR
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_RG
!
END SUBROUTINE LIMA_DROPS_HOM_FREEZING
END INTERFACE
END MODULE MODI_LIMA_DROPS_HOM_FREEZING
!
!     ######################################################################
      SUBROUTINE LIMA_DROPS_HOM_FREEZING (PTSTEP, HFMFILE, OCLOSE_OUT, LDCOMPUTE,   &
                                          PEXNREF, PPABST,                          &
                                          PTHT, PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                                          PCRT,                                     &
                                          P_TH_HONR, P_RR_HONR, P_CR_HONR,          &
                                          PB_TH, PB_RR, PB_CR, PB_RG                )
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
USE MODD_PARAM_LIMA,      ONLY : XRTMIN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                  INTENT(IN)    :: PTSTEP 
CHARACTER(LEN=*),      INTENT(IN)    :: HFMFILE
LOGICAL,               INTENT(IN)    :: OCLOSE_OUT
LOGICAL, DIMENSION(:), INTENT(IN)    :: LDCOMPUTE
!
REAL, DIMENSION(:),    INTENT(IN)    :: PEXNREF ! Reference Exner function
REAL, DIMENSION(:),    INTENT(IN)    :: PPABST  ! abs. pressure at time t
!
REAL, DIMENSION(:),    INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:),    INTENT(IN)    :: PRVT    ! Water vapor m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRCT    ! Cloud water m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRRT    ! Rain water m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRIT    ! Cloud ice m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRST    ! Snow/aggregate m.r. at t 
REAL, DIMENSION(:),    INTENT(IN)    :: PRGT    ! Graupel m.r. at t 
!
REAL, DIMENSION(:),    INTENT(IN)    :: PCRT    ! Rain water C. at t
!
REAL, DIMENSION(:),    INTENT(INOUT) :: P_TH_HONR
REAL, DIMENSION(:),    INTENT(INOUT) :: P_RR_HONR
REAL, DIMENSION(:),    INTENT(INOUT) :: P_CR_HONR
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_TH
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_RR
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_CR
REAL, DIMENSION(:),    INTENT(INOUT) :: PB_RG
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(PTHT)) ::  &
     ZW,       &
     ZT,       &
     ZLSFACT,  &
     ZLVFACT,  &
     ZTCELSIUS
!
!
!
!
!-------------------------------------------------------------------------------
!
P_TH_HONR(:) = 0.
P_RR_HONR(:) = 0.
P_CR_HONR(:) = 0.
!
! Temperature
ZT(:)  = PTHT(:) * ( PPABST(:)/XP00 ) ** (XRD/XCPD)
ZTCELSIUS(:) = ZT(:)-XTT                                    ! T [°C]
!
ZW(:) = PEXNREF(:)*( XCPD+XCPV*PRVT(:)+XCL*(PRCT(:)+PRRT(:)) &
     +XCI*(PRIT(:)+PRST(:)+PRGT(:)) )
ZLSFACT(:) = (XLSTT+(XCPV-XCI)*ZTCELSIUS(:))/ZW(:)          ! L_s/(Pi_ref*C_ph)
ZLVFACT(:) = (XLVTT+(XCPV-XCL)*ZTCELSIUS(:))/ZW(:)          ! L_v/(Pi_ref*C_ph)
!
ZW(:) = 0.0

WHERE( (ZT(:)<XTT-35.0) .AND. (PRRT(:)>XRTMIN(3)) .AND. LDCOMPUTE(:) )
   P_TH_HONR(:) = PRRT(:)*(ZLSFACT(:)-ZLVFACT(:))
   P_RR_HONR(:) = - PRRT(:)
   P_CR_HONR(:) = - PCRT(:)
   PB_TH(:) = PB_TH(:) + P_TH_HONR(:)
   PB_RR(:) = PB_RR(:) - PRRT(:)
   PB_CR(:) = PB_CR(:) - PCRT(:)
   PB_RG(:) = PB_RG(:) + PRRT(:)
ENDWHERE
!
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_DROPS_HOM_FREEZING
