!     ######spl
      MODULE MODD_CH_SOLVER_n
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #########################
!!
!!*** *MODD_CH_SOLVER$n*
!!
!!    PURPOSE
!!    -------
!     contains parameters for the stiff solvers
!!
!!**  AUTHOR
!!    ------
!!    K. Suhre     *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/03/95
!!    27/07/96 (K. Suhre) restructured
!!    01/08/01 (C. Mari)  change routine to $n
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    none
!!
!------------------------------------------------------------------------------
!
!*       0.  DECLARATIONS
!        ----------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_SOLVER_t
!
!*       0.1  choice of the stiff solver
!
  CHARACTER*32 :: CSOLVER = 'EXQSSA'  ! name of the solver to be used
!
!*       0.2  parameters for LinSSA solver
!
  INTEGER                  :: NSSA = 0  ! number of variables to be treated
                                        ! as "steady state"
!JUAN
  INTEGER, DIMENSION(:), POINTER  :: NSSAINDEX => NULL() ! index set of steady state variables
!JUAN
!
!*       0.3  tolerances (used by NAG and SVODE solvers)
!
  REAL :: XRTOL = 0.001 ! relative tolerance for SVODE
                        ! and D02EAF,D02EBF,D02NBF
                      !
  REAL :: XATOL = 0.1   ! absolute tolerance for SVODE
                        ! and D02NBF
!
!*       0.4  parameters for NAG's D02EBF solver
!
  INTEGER :: NRELAB = 2 ! choose relative error:
                      ! 1 for correct decimal places
                      ! 2 for correct significant digits
                      ! 0 for a mixture
!
!*       0.5  parameters for SVODE and NAG's D02EBF/D02NBF solvers
!
  INTEGER :: NPED   = 1 ! calculation of the Jacobian matric:
                      ! 0 for numerical Jacobian
                      ! 1 for analytical Jacobian (using subroutine CH_JAC)
!
!*       0.6  parameters for NAG's D02NBF solver
!
  INTEGER     :: NMAXORD = 5   ! maximum order for the BDF method (0<NMAXORD<=5)
  LOGICAL     :: LPETZLD = .TRUE. ! perform Petzold local error test (recommended)
  CHARACTER*1 :: CMETHOD = "N" ! method to use non-linear system
                             ! N or D for modified Newton iteration
                             ! F for functional iteration
  CHARACTER*1 :: CNORM = "A"   ! type of norm to be used
                             ! A or D for averaged L2 norm
                             ! M for maximum norm
  INTEGER     :: NTRACE = 0    ! level of output from D02NBF
                             ! -1 (no output) <= NTRACE <= 3 (maximum)
                             !  0 only warnings are printed
                             ! >0 details on Jacobian entries, nonlinear
                             !    iteration and time integration are given
!
!*       0.7  parameters for CRANCK solver
!
  REAL :: XALPHA = 0.5         ! the Cranck-Nicholson parameter (0,1)
!
!*       0.8  parameters for (EX)QSSA solvers
!
  REAL    :: XSLOW     = 100.0 ! slow species, lifetime > XSLOW * timestep
  REAL    :: XFAST     = 0.1   ! fast species, lifetime < XFAST * timestep
  INTEGER :: NQSSAITER = 1     ! number of iterations in QSSA
  REAL    :: XDTMIN    = 0.1   ! minimal allowed timestep for EXQSSA
  REAL    :: XDTMAX    = 10.   ! maximal allowed timestep for EXQSSA
  REAL    :: XDTFIRST  = 10.   ! timestep for first integration step of EXQSSA
!
END TYPE CH_SOLVER_t

TYPE(CH_SOLVER_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_SOLVER_MODEL
!JUAN
LOGICAL          , DIMENSION(JPMODELMAX),         SAVE :: CH_SOLVER_FIRST_CALL = .TRUE.
!JUAN

CHARACTER*32, POINTER :: CSOLVER=>NULL()
INTEGER, POINTER :: NSSA=>NULL()
INTEGER, DIMENSION(:), POINTER :: NSSAINDEX=>NULL()
REAL, POINTER :: XRTOL=>NULL()
REAL, POINTER :: XATOL=>NULL()
INTEGER, POINTER :: NRELAB=>NULL()
INTEGER, POINTER :: NPED=>NULL()
INTEGER, POINTER :: NMAXORD=>NULL()
LOGICAL, POINTER :: LPETZLD=>NULL()
CHARACTER*1, POINTER :: CMETHOD=>NULL()
CHARACTER*1, POINTER :: CNORM=>NULL()
INTEGER, POINTER :: NTRACE=>NULL()
REAL, POINTER :: XALPHA=>NULL()
REAL, POINTER :: XSLOW=>NULL()
REAL, POINTER :: XFAST=>NULL()
INTEGER, POINTER :: NQSSAITER=>NULL()
REAL, POINTER :: XDTMIN=>NULL()
REAL, POINTER :: XDTMAX=>NULL()
REAL, POINTER :: XDTFIRST=>NULL()

CONTAINS

SUBROUTINE CH_SOLVER_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
!JUAN
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODD_CH_SOLVER_N:CH_SOLVER_GOTO_MODEL',0,ZHOOK_HANDLE)
IF (CH_SOLVER_FIRST_CALL(KTO)) THEN
ALLOCATE (CH_SOLVER_MODEL(KTO)%NSSAINDEX(1000))
CH_SOLVER_FIRST_CALL(KTO) = .FALSE.
ENDIF
!JUAN

! Save current state for allocated arrays
!
! Current model is set to model KTO
CSOLVER=>CH_SOLVER_MODEL(KTO)%CSOLVER
NSSA=>CH_SOLVER_MODEL(KTO)%NSSA
NSSAINDEX=>CH_SOLVER_MODEL(KTO)%NSSAINDEX
XRTOL=>CH_SOLVER_MODEL(KTO)%XRTOL
XATOL=>CH_SOLVER_MODEL(KTO)%XATOL
NRELAB=>CH_SOLVER_MODEL(KTO)%NRELAB
NPED=>CH_SOLVER_MODEL(KTO)%NPED
NMAXORD=>CH_SOLVER_MODEL(KTO)%NMAXORD
LPETZLD=>CH_SOLVER_MODEL(KTO)%LPETZLD
CMETHOD=>CH_SOLVER_MODEL(KTO)%CMETHOD
CNORM=>CH_SOLVER_MODEL(KTO)%CNORM
NTRACE=>CH_SOLVER_MODEL(KTO)%NTRACE
XALPHA=>CH_SOLVER_MODEL(KTO)%XALPHA
XSLOW=>CH_SOLVER_MODEL(KTO)%XSLOW
XFAST=>CH_SOLVER_MODEL(KTO)%XFAST
NQSSAITER=>CH_SOLVER_MODEL(KTO)%NQSSAITER
XDTMIN=>CH_SOLVER_MODEL(KTO)%XDTMIN
XDTMAX=>CH_SOLVER_MODEL(KTO)%XDTMAX
XDTFIRST=>CH_SOLVER_MODEL(KTO)%XDTFIRST

IF (LHOOK) CALL DR_HOOK('MODD_CH_SOLVER_N:CH_SOLVER_GOTO_MODEL',1,ZHOOK_HANDLE)
END SUBROUTINE CH_SOLVER_GOTO_MODEL

END MODULE MODD_CH_SOLVER_n
