!     ######spl
      SUBROUTINE CH_INIT_OUTPUT(TPM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ##############################
!!
!!*** *CH_INIT_OUTPUT*
!!
!!    PURPOSE
!!    -------
!!    prepare regular output of results
!!
!!**  METHOD
!!    ------
!!    the result file CRESULTFILE (default: "BOX.RESULT") will be opened
!!    and the headder will be written, containing general information
!!    on the saved variables, presently the following variables will be saved:
!!    - concentration of all prognostic variables
!!    - the reaction constants and photolysis rates
!!    - the meteo variables
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
!!    Original 25/04/95
!!    27/07/96 (K. Suhre) restructured
!!    01/12/03 (D. Gazen) change Chemical scheme interface
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D,  ONLY: CRESULTFILE, CRUNID, CRESULTFORMAT, &
                            NRESULTIO, XTSIMUL
USE MODD_CH_M9,       ONLY: NEQ, NREAC, NMETEOVARS, CNAMES, CREACS, METEOTRANSTYPE
USE MODD_CH_AEROSOL
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
TYPE(METEOTRANSTYPE), INTENT(IN) :: TPM  ! the meteo variables
!
!*       0.2  declaration of local variables
!
CHARACTER*8            :: YDATE  ! for retrieval of date and time
CHARACTER*10           :: YTIME  ! dito
INTEGER                :: JI     ! loop control
INTEGER                :: NAERO
!
!------------------------------------------------------------------------------
!
!*       1.   OPEN OUTPUT FILE
!        ----------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_INIT_OUTPUT',0,ZHOOK_HANDLE)
PRINT *, "CH_INIT_OUTPUT: opening unit ", NRESULTIO, " for file ", CRESULTFILE
OPEN(UNIT   = NRESULTIO,   &
     FILE   = CRESULTFILE, &
     FORM   = "FORMATTED", &
     STATUS = "UNKNOWN"    )
!
!*       2.   WRITE HEADDER
!        ------------------
!
!*       2.1  write two lines of comment (date, time, runid)
!
IF (LORILAM) THEN
NAERO=SIZE(CAERONAMES)
ELSE
NAERO=0
END IF

CALL DATE_AND_TIME(YDATE, YTIME)
WRITE(NRESULTIO,'(4A)')'MODEL0D RESULTS VERSION 1.1 AT ',YDATE,":",YTIME
WRITE(NRESULTIO,'(A)')  CRUNID
!
!*       2.2  write number of variables and their names
!
WRITE(NRESULTIO,'(4I6)') 1+NEQ+NAERO+NREAC+NMETEOVARS, NEQ, NREAC, NMETEOVARS
!
!*       2.3  write simulation time
!
WRITE(NRESULTIO, '(A)') "XTSIMUL"
!
!*       2.4  write chemical concentrations
!
DO JI = 1, NEQ
  WRITE(NRESULTIO, '(A)') CNAMES(JI)
ENDDO
DO JI = 1, NAERO
  WRITE(NRESULTIO, '(A)') CAERONAMES(JI)
ENDDO
!
!*       2.5  write reaction and photolysis rates
!
DO JI = 1, NREAC
  WRITE(NRESULTIO, '(A)') CREACS(JI)
ENDDO
!
!*       2.6  write meteo variables
!
DO JI = 1, NMETEOVARS
  WRITE(NRESULTIO, '(A)') TPM%CMETEOVAR(JI)
ENDDO
!
!*       2.7  write data output format
!
WRITE(NRESULTIO,'(A)') CRESULTFORMAT
!
IF (LHOOK) CALL DR_HOOK('CH_INIT_OUTPUT',1,ZHOOK_HANDLE)
END SUBROUTINE CH_INIT_OUTPUT
