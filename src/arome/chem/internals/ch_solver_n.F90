!     ######spl
      SUBROUTINE CH_SOLVER_n(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ######################################################################### 
!!
!!*** *CH_SOLVER_n*
!!
!!    PURPOSE
!!    -------
!!      solution of one timestep of the chemical differential equation
!!
!!**  METHOD
!!    ------
!!      Calls the individual solvers and passes the corresponding parameters
!!    on. A fixed integration step is used, however, some solvers internally
!!    use variable time steps.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/04/95
!!    10/11/95 KS: integrate SVODE solver (AFF)
!!    10/04/96 KS: integrate QSSA and EXQSSA solver (AFF)
!!    31/07/96 KS: add TPK to parameterlist (some solver desactivated !@)
!!    01/08/01 (C. Mari)  change routine to _n
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!    01/06/07 (P. Tulet)  model number in argument (for AROME)
!!
!!    EXTERNAL
!!    --------
!!    calls the different solvers like SIS, LinSSA, NAG, etc 
USE MODI_CH_SIS
USE MODI_CH_LINSSA
USE MODI_CH_CRANCK
!USE MODI_CH_SVODE
USE MODI_CH_QSSA
USE MODI_CH_EXQSSA
!@USE MODI_CH_D02EAF
!@USE MODI_CH_D02EBF
!@USE MODI_CH_D02NBF
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_SOLVER_n
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
INTEGER, INTENT(IN) :: KMI      ! model number
REAL, INTENT(IN),  DIMENSION(KVECNPT,KEQ) :: PCONC    ! concentration vector at PTSIMUL
REAL, INTENT(OUT), DIMENSION(KVECNPT,KEQ) :: PNEWCONC ! solution at PTSIMUL + PDTACT
!
!------------------------------------------------------------------------------
!
!
!*       1.   CALL THE SOLVERS
!        ---------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_SOLVER_N',0,ZHOOK_HANDLE)
SELECT CASE (CSOLVER)
!
CASE ('SIS')
!
  ! call SIS
  CALL CH_SIS(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT,KMI)
  !
CASE ('LINSSA', 'LinSSA')
!
  ! call LinSSA
  CALL CH_LINSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                 NSSA, NSSAINDEX)
!
CASE ('CRANCK')
!
  ! call Cranck-Nicholson method
  CALL CH_CRANCK(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, & 
                 XALPHA)
!
CASE ('D02EAF')
!
  ! call NAG's stiff-solver D02EAF
  !@CALL CH_D02EAF(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KMI)
  STOP 'CH_SOLVER_n SORRY: requested solver currently not supported (CSOLVER)'
!
CASE ('D02EBF')
!
  ! call NAG's stiff-solver D02EBF
  !@CALL CH_D02EBF(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KMI)
  STOP 'CH_SOLVER_n SORRY: requested solver currently not supported (CSOLVER)'
!
CASE ('D02NBF')
!
  ! call NAG's stiff-solver D02NBF
  !@CALL CH_D02NBF(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KMI)
  STOP 'CH_SOLVER_n SORRY: requested solver currently not supported (CSOLVER)'
!
CASE ('SVODE')
!
  ! call SVODE solver

 ! CALL CH_SVODE(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
 !               XRTOL, XATOL, NPED)
  STOP 'CH_SOLVER_n SORRY: requested solver currently not supported (CSOLVER) until Masdev47'
!
CASE ('QSSA')
!
  ! call QSSA
  CALL CH_QSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
               XSLOW, XFAST, NQSSAITER)
  !
CASE ('EXQSSA')
!
  ! call EXQSSA
  CALL CH_EXQSSA(PTSIMUL, PDTACT, PCONC, PNEWCONC, KEQ, KVECNPT, KMI, &
                 XSLOW, XFAST, XDTMAX)
  !
CASE ('NONE')
!
  ! no integration at all (for debugging purposes)
  PNEWCONC(:,:) = PCONC(:,:)
!
CASE DEFAULT
  STOP 'CH_SOLVER_n ERROR: requested solver not supported (CSOLVER)'
END SELECT
!
IF (LHOOK) CALL DR_HOOK('CH_SOLVER_N',1,ZHOOK_HANDLE)
END SUBROUTINE CH_SOLVER_n
