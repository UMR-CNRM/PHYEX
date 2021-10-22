!     ######spl
      SUBROUTINE CH_READ_VECTOR(KEQ, HNAMES, PVAR, PDEFAULT, KIN, KOUT, KVERB)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #########################################################################
!!
!!*** *CH_READ_VECTOR*
!!
!!    PURPOSE
!!    -------
!       read a vector of chemical parameters
!!
!!**  METHOD
!!    ------
!!      This subroutine reads from an already opened file (channel KIO)
!!    pairs of character variables and values. It then associates the
!!    values to an array of variable names that has been passed as argument.
!!    In practice, this program will be used to read for example the initial
!!    conditions of a vector of chemical variables.
!!      The format of the data to be entered is as follows:
!!    line 1:   number of variables to be read
!!    line 2:   FORTRAN format         e.g.  (A10,E12.5)
!!    line 3++: species   value        e.g.  'O3' 5.5E10
!!
!!      If a species is not defined in the file it will be assigned the
!!    default value. If KVERB is greater or equal to 10, the full
!!    association will be printed to KOUT.
!!
!!    REFERENCE
!!    ---------
!!      none
!!
!!    AUTHOR
!!    ------
!!    K. Suhre    *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 01/04/1999
!!
!!    EXTERNAL
!!    --------
!!      none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      none
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,                       INTENT(IN) :: KEQ
                                    ! number of variables to be defined
CHARACTER*(*), DIMENSION(KEQ), INTENT(IN) :: HNAMES
                                    ! names of the variables to be defined
REAL,          DIMENSION(KEQ), INTENT(OUT):: PVAR
                                    ! value of the variable to be read
REAL,                          INTENT(IN) :: PDEFAULT
                                    ! default value
INTEGER,                       INTENT(IN) :: KIN
                                    ! I/O channel for file input
INTEGER,                       INTENT(IN) :: KOUT
                                    ! I/O channel for printing
INTEGER,                       INTENT(IN) :: KVERB
                                    ! verbosity level
!
!*      0.2    declarations of local variables
!
CHARACTER(LEN=40) :: YFORMAT    ! format for input
!
INTEGER :: IINPUT               ! number of surface values to be read
CHARACTER(LEN=40), DIMENSION(:), ALLOCATABLE :: YINPUTNAME ! species names
REAL             , DIMENSION(:), ALLOCATABLE :: ZINPUTVAL  ! surface val.
!
INTEGER :: JI, JN, JNREAL ! loop control variables
INTEGER :: INACT          ! array pointer
!
!-------------------------------------------------------------------------------
!
!*       1.    ASSOCIATE VALUES TO NAMES
!              -------------------------
!
!*       1.1   read input values
!
! read number of input lines IINPUT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_READ_VECTOR',0,ZHOOK_HANDLE)
READ(KIN, *) IINPUT
!
! read data input format
READ(KIN,"(A)") YFORMAT
!
! allocate input fields
ALLOCATE(YINPUTNAME(IINPUT))
ALLOCATE(ZINPUTVAL(IINPUT))
!
! read input values
IF (KVERB >= 10) THEN
  WRITE(KOUT,'(A)') '----------------------------------------------------'
  WRITE(KOUT,'(A)') 'CH_READ_VECTOR: reading ...'
END IF
DO JI = 1, IINPUT
  READ(KIN,YFORMAT) YINPUTNAME(JI), ZINPUTVAL(JI)
  IF (KVERB >= 10) WRITE(KOUT,YFORMAT) YINPUTNAME(JI), ZINPUTVAL(JI)
END DO
!
!*       1.2   map values on PVAR
!
IF (KVERB >= 10) THEN
  WRITE(KOUT,'(A)') '----------------------------------------------------'
  WRITE(KOUT,'(A)') 'CH_READ_VECTOR: associating ...'
  WRITE(KOUT,'(I4)') KEQ
  WRITE(KOUT,'(A)') YFORMAT
END IF
!
PVAR(:) = PDEFAULT
!
DO JNREAL = 1, KEQ
  INACT = 0
  search_loop1 : DO JN = 1, IINPUT
    IF (HNAMES(JNREAL) .EQ. YINPUTNAME(JN)) THEN
      INACT = JN
      EXIT search_loop1
    END IF
  END DO search_loop1
  IF (INACT .NE. 0) PVAR(JNREAL) = ZINPUTVAL(INACT)
  IF (KVERB >= 10) THEN
    WRITE(KOUT,YFORMAT) HNAMES(JNREAL), PVAR(JNREAL)
  END IF
END DO
!
IF (KVERB >= 10) &
  WRITE(KOUT,'(A)') '----------------------------------------------------'
!
DEALLOCATE(YINPUTNAME)
DEALLOCATE(ZINPUTVAL)
!
IF (LHOOK) CALL DR_HOOK('CH_READ_VECTOR',1,ZHOOK_HANDLE)
RETURN
!
END SUBROUTINE CH_READ_VECTOR
