!     ######spl
SUBROUTINE CH_LINSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                     KSSA, KSSAINDEX)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #########################################################################
!!
!!*** *CH_LINSSA*
!!
!!**  PURPOSE
!!    -------
!!    solve one time-step of the chemical differential equation d/dt C = f(C)
!!
!!    METHOD
!!    ------
!!    Linearised Steady (LINSSA):
!!    y^n+1 = y^n + dt/2 * (P-dt/2 J^n)^-1 * (P+I) * f^n
!!
!!    REFERENCE
!!    ---------
!!    K. Suhre and R. Rosset,
!!    Modification of a linearized semi-implicit scheme for chemical
!!    reactions using a steady-state-approximation,
!!    Ann. Geophysicae, 12, 359--361, 1994.
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 25/04/95
!!    Modification   01/12/03  (Gazen)   change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_FCN
USE MODI_CH_JAC
USE MODI_CH_GAUSS
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
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
INTEGER, INTENT(IN) :: KMI ! model number
INTEGER, INTENT(IN) :: KSSA
INTEGER, DIMENSION(1000)  ::  KSSAINDEX
!
!*       0.2  declaration of local variables
!
REAL,   DIMENSION(KVECNPT,KEQ)     :: ZDCDTVECT
REAL,   DIMENSION(KEQ)     :: ZDCDT
REAL,   DIMENSION(KVECNPT,KEQ,KEQ) :: ZWORKVECT
REAL,   DIMENSION(KEQ,KEQ) :: ZWORK,ZINV
INTEGER                    :: IFAIL, JI, JJ
REAL,   DIMENSION(:,:), SAVE, ALLOCATABLE :: ZPROJECTOR, ZFAC
LOGICAL,              SAVE              :: GSFIRSTCALL = .TRUE.
!
!------------------------------------------------------------------------------
!
!
!*       1.   ON FIRSTCALL, THE PROJECTOR WILL BE DEFINED
!        ------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_LINSSA',0,ZHOOK_HANDLE)
first_call : IF (GSFIRSTCALL) THEN
  GSFIRSTCALL = .FALSE.
  ALLOCATE(ZPROJECTOR(KVECNPT,KEQ))
  ALLOCATE(ZFAC(KVECNPT,KEQ))
  ZPROJECTOR(:,:) = 1.0
  DO JI = 1, KSSA
    ZPROJECTOR(:,KSSAINDEX(JI)) = 0.0
  ENDDO
!
  ZFAC(:,:) = 0.5*(1.0 + ZPROJECTOR(:,:))
!
ENDIF first_call
!
!*       2.   CALCULATE FIRST DERIVATIVE AND JACOBIAN
!        --------------------------------------------
!
CALL CH_FCN(PTSIMUL,PCONC,ZDCDTVECT, KMI,KVECNPT,KEQ)
!
!%kloos%
!ZJAC replaced by ZWORK
CALL CH_JAC(PTSIMUL,PCONC,ZWORKVECT, KMI,KVECNPT,KEQ)
!
!*       3.   CALCULATE (P-dt/2 J^n)^-1
!        ------------------------------
!
DO JJ=1,KVECNPT

ZWORK(:,:) = - (0.5*PDTACT) * ZWORKVECT(JJ,:,:)
ZDCDT(:)=ZDCDTVECT(JJ,:)
!%
DO JI = 1, KEQ
  ZWORK(JI,JI) = ZWORK(JI,JI) + ZPROJECTOR(JJ,JI)
ENDDO
!
!%kloos%
! inversion of matrix modified
!
!*       4.   CALCULATE LU FACTORIZATION FOR MATRIX (P-dt/2 J^n)^-1
!        ----------------------------------------------------------
!           (result is put in ZWORK)
!
IFAIL = 1
CALL CH_GAUSS(ZWORK,ZINV,KEQ,IFAIL)
IF (IFAIL.NE.0) THEN
  STOP 'CH_LinSSA ERROR: matrix cannot be inverted by CH_GAUSS'
ENDIF
!
!*       5.   CALCULATE 1/2 * (P+I) * f^n
!        --------------------------------
!
ZDCDT(:) = ZFAC(JJ,:) * ZDCDT(:)
!
!*       6.   CALCULATE (1-dt/2 J^n)^-1 * (P+I)/2 * f^n (result is put in ZDCDT)
!        -----------------------------------------------------------------------
!
ZDCDT(:) = MATMUL(ZINV,ZDCDT)
!
!*       7.   CALCULATE y^n+1 = y^n + dt/2 * (P-dt/2 J^n)^-1 * (P+I) * f^n
!        -----------------------------------------------------------------
!
PNEWCONC(JJ,:) = PCONC(JJ,:) + PDTACT * ZDCDT(:)
!%
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CH_LINSSA',1,ZHOOK_HANDLE)
END SUBROUTINE CH_LINSSA
