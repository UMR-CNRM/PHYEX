!     ######spl
      SUBROUTINE CH_SHOW_CHEM(PCONC, HNAMES)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #############################
!!
!!****  *CH_SHOW_CHEM*
!!
!!    PURPOSE
!!    -------
!!    print a set of values to the screen
!!
!!    METHOD
!!    ------
!!    NEQ values will be written to the screen,
!!    the format will be as follows:
!!    line 1:     'name1'     value1
!!      ...
!!    line NEQ:   'nameNEQ'     valueNEQ
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 21/04/95
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
!!    PCONC: concentration vector to be read
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN) :: PCONC ! PCONC: concentration vector
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES ! name of each species
!
!!    DECLARATION OF LOCAL VARIABLES
!!    ------------------------------
INTEGER :: JI
!
!!    EXECUTABLE STATEMENTS
!!    ---------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_SHOW_CHEM',0,ZHOOK_HANDLE)
DO JI = 1, SIZE(PCONC)
  PRINT '(3A,E20.8)', "'", HNAMES(JI), "' ", PCONC(JI)
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CH_SHOW_CHEM',1,ZHOOK_HANDLE)
END SUBROUTINE CH_SHOW_CHEM
