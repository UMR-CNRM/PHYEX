!     ######spl
SUBROUTINE CH_QSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                   PSLOW, PFAST, KQSSAITER)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #######################################################################
!!
!!*** *CH_QSSA*
!!
!!    PURPOSE
!!    -------
!!    solve one time-step of the chemical differential equation d/dt C = P-L*C
!!
!!    METHOD
!!    ------
!!    Quasy steady state approximation (QSSA):
!!      C_n+1 = P_n / L_n                      if lifetime < Dt/10
!!      C_n+1 = C_n + ( P_n - L_n * C_n) Dt    if lifetime > 100 * Dt
!!      C_n+1 = P_n / L_n + ( C_n - P_n / L_n ) exp( - L_n Dt)  else
!!
!!    REFERENCE
!!    ---------
!!    O. Hertel et al ,
!!    Test of two Numerical schemes for use in atmospheric transport
!!    chemistry models,
!!    Atmospheric environment, 16, 2591--2611, 1993.
!!
!!    AUTHOR
!!    ------
!!    A. Fassi Fihri   *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 18/12/95
!!    10/04/96 (K. Suhre) take out all performance analysis and
!!                        internal time-stepping
!!    31/07/96 (K. Suhre) restructured
!!    01/08/01 (C. Mari)  add PSLOW, PFAST, KQSSAITER arguments
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_PRODLOSS
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI      ! model number
REAL,    INTENT(IN) :: PSLOW, PFAST
INTEGER, INTENT(IN) :: KQSSAITER
!
!*       0.2  declaration of local variables
!
REAL,    DIMENSION(KVECNPT,KEQ)  :: ZPROD, ZLOSS, ZCONC0, ZTMP
LOGICAL, DIMENSION(KVECNPT,KEQ)  :: GSHORT, GMED, GLONG
                                ! mask for short, medium and long lived species
INTEGER :: ICOUNT         ! counter for iteration number
!
!------------------------------------------------------------------------------
!
!*       1.   PREPARE VARIABLES
!        ----------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_QSSA',0,ZHOOK_HANDLE)
PNEWCONC(:,:) = PCONC(:,:)
!
!*       2.   ITERATE NQSSAITER TIMES
!        ----------------------------
!
DO ICOUNT = 1, KQSSAITER
!
  ZCONC0(:,:) = PNEWCONC(:,:)
!
!*       2.1  calculate production and loss terms
!
  CALL CH_PRODLOSS(PTSIMUL,ZCONC0,ZPROD,ZLOSS,KMI,KVECNPT,KEQ)
!
!*       2.2  decide which species are short, medium and long lived species
!
  GSHORT(:,:) = ZLOSS(:,:) .GT. (1./(PFAST*PDTACT))
  GLONG(:,:)  = ZLOSS(:,:) .LT. (1./(PSLOW*PDTACT))
  GMED(:,:)   = .NOT. (GSHORT(:,:) .OR. GLONG(:,:))
!
!*       2.3  analytical solution assuming PROD and LOSS constant throughout
!             a time-step (medium long lived species)
!
  WHERE (GMED)
    ZTMP(:,:) = ZPROD(:,:)/ZLOSS(:,:)
    PNEWCONC(:,:) = ZTMP(:,:) + (ZCONC0(:,:) - ZTMP(:,:))*EXP(-ZLOSS(:,:)*PDTACT)
  END WHERE
!
!*       2.4  euler explicit integration (long lived species)
!
  WHERE ( GLONG )
    PNEWCONC(:,:) = ZCONC0(:,:) + (ZPROD(:,:) - PDTACT*ZLOSS(:,:)*ZCONC0(:,:))
  ENDWHERE
!
!*       2.5  steady state approximation (short lived species)
!
  WHERE ( GSHORT )
    PNEWCONC(:,:) = ZPROD(:,:) / ZLOSS(:,:)
  END WHERE
!
END DO
!
IF (LHOOK) CALL DR_HOOK('CH_QSSA',1,ZHOOK_HANDLE)
END SUBROUTINE CH_QSSA
