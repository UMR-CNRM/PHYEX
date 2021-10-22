!     ######spl
      SUBROUTINE CH_READ_METEO(TPM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #############################
!!
!!***  *CH_READ_METEO*
!!
!!    PURPOSE
!!    -------
!!    Read a set of meteo variables
!!
!!**  METHOD
!!    ------
!!    read NMETEOVARS values and the time for the next update XTNEXTMETEO
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 21/04/95
!!    27/07/96 (K. Suhre) restructured
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_MODEL0D,  ONLY: NMETEOIO, XTNEXTMETEO, NVERB
USE MODD_CH_M9,       ONLY: NMETEOVARS, METEOTRANSTYPE
!!
!------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!        -----------------
IMPLICIT NONE
!
!*       0.1  declaration of arguments
!
TYPE(METEOTRANSTYPE), INTENT(INOUT) :: TPM  ! the meteo variables
!
!*       0.2  declaration of local variables
!     ----------------
INTEGER :: JI   ! loop control
!
!------------------------------------------------------------------------------
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
! read meteo variables and time of next update
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_READ_METEO',0,ZHOOK_HANDLE)
READ(NMETEOIO,*) (TPM%XMETEOVAR(JI), JI = 1, NMETEOVARS)
READ(NMETEOIO,*) XTNEXTMETEO
!
! print what has been read
IF (NVERB >= 7) THEN
  PRINT *, 'CH_READ_METEO: new set of meteo variables has been read:'
  DO JI = 1, NMETEOVARS
    PRINT *, TPM%CMETEOVAR(JI), ': ', TPM%XMETEOVAR(JI)
  ENDDO
END IF
IF (NVERB >= 5) THEN
  PRINT *, 'CH_READ_METEO: next update at XTNEXTMETEO = ', XTNEXTMETEO
END IF
!
IF (LHOOK) CALL DR_HOOK('CH_READ_METEO',1,ZHOOK_HANDLE)
END SUBROUTINE CH_READ_METEO
