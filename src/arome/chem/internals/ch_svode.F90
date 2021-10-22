!     ######spl
SUBROUTINE CH_SVODE(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                    PRTOL, PATOL, KPED )
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ########################################################################
!!
!!****  *CH_SVODE*

!!    PURPOSE
!!    -------
!!    solve one time-step of the chemical differential equation d/dt C = f(C)

!!    METHOD
!!    ------
!!    the stiff-solver SVODE (written in FORTRAN77)
!!    will be used to solve the system
!!    NOTE: this subroutine is not 100% doctorized since we want
!!          to conserve quasi-standard variable names used by SVODE

!!    REFERENCE
!!    ---------
!!    MesoNH book 2

!!    AUTHOR
!!    ------
!!    K. Suhre

!!    MODIFICATIONS
!!    -------------
!!    Original 10/11/95
!!    01/08/01 (C. Mari)  add arguments
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
!!    SVODE: stiff solver from NETLIB
!!    CH_SVODE_FCN, CH_SVODE_JAC: argument list converters for CH_FCN, CH_JAC

!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
REAL,    INTENT(IN) :: PTSIMUL  ! time of simulation
REAL,    INTENT(IN) :: PDTACT   ! actual time-step
INTEGER, INTENT(IN) :: KEQ      ! dimension of the problem to solve
INTEGER, INTENT(IN) :: KVECNPT
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC
                                ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC
                                ! solution at PTSIMUL + PDTACT
INTEGER, INTENT(IN) :: KMI      ! model index
REAL,    INTENT(IN) :: PRTOL, PATOL
INTEGER, INTENT(IN) :: KPED

!!    DECLARATION OF LOCAL VARIABLES
!!    ------------------------------
REAL, DIMENSION(KEQ) :: ZCONC  ! concentration vector
REAL :: ZTBEGIN, ZTEND         ! begin and end of integration step

INTEGER :: ITOL    = 1      ! type of error, don't play with this parameter
INTEGER :: ITASK   = 1      ! integrate from PTSIMUL to PTSIMUL + PDTACT
INTEGER :: ISTATE           ! integer flag (input and output)
INTEGER :: IOPT    = 0      ! 0 to indicate no optional input used
INTEGER :: IMF              ! 21=use stiff full analytical Jacobian
                            ! 22=use stiff full numerical Jacobian
REAL, DIMENSION(KEQ) :: ZRTOL          ! relative tolerance (XRTOL)
REAL, DIMENSION(KEQ) :: ZATOL          ! absolute tolerance (XATOL)

! workspace declaration
INTEGER :: IRW ! will be set to (22 + 9*KEQ + 2*KEQ*KEQ) in the code
INTEGER :: IIW ! will be set to (30 + KEQ) in the code
REAL, DIMENSION(22 + 9*KEQ + 2*KEQ*KEQ) :: ZWORK ! workspace
INTEGER, DIMENSION(30 + KEQ)            :: IWORK ! workspace

! dummy parameters
INTEGER, DIMENSION(1) :: IPAR ! dummy parameter
REAL, DIMENSION(1)    :: ZPAR ! dummy parameter

INTEGER :: JI ! loop counter

EXTERNAL CH_SVODE_FCN, CH_SVODE_JAC

!*    EXECUTABLE STATEMENTS
!     ---------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_SVODE',0,ZHOOK_HANDLE)
ZTBEGIN = PTSIMUL
ZTEND = PTSIMUL + PDTACT

ZRTOL(:) = PRTOL
ZATOL(:) = PATOL

! set array dimensions
IRW = (22 + 9*KEQ + 2*KEQ*KEQ)
IIW = (30 + KEQ)

! choose calculation of Jacobian
IF (KPED .EQ. 0) THEN
  IMF = 22        ! numerical calculation of Jacobian
ELSE
  IMF = 21        ! analytical calculation of Jacobian using CH_JAC
ENDIF

! at each call to SVODE start a new iteration cycle
ISTATE = 1

! this solver is not vectorized, we loop over all elements
DO JI = 1, KVECNPT

  ZCONC(:) = PCONC(JI,:)

  ! call SVODE solver
  CALL SVODE (CH_SVODE_FCN, KEQ, ZCONC, ZTBEGIN, ZTEND, &
              ITOL, ZRTOL, ZATOL, ITASK, &
              ISTATE, IOPT, ZWORK, IRW, IWORK, IIW, CH_SVODE_JAC, IMF, &
              ZPAR, IPAR, KMI, JI)

  IF (ISTATE.LT.0) THEN
    PRINT *, "Problems !!! ISTATE = ", ISTATE
    PRINT *, "at vector element ", JI, " out of ", KVECNPT
    STOP "CH_SVODE: program stopped due to SVODE error!"
  ENDIF

  PNEWCONC(JI,:) = ZCONC(:)

END DO

IF (LHOOK) CALL DR_HOOK('CH_SVODE',1,ZHOOK_HANDLE)
END SUBROUTINE CH_SVODE
