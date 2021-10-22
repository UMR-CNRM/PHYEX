!     ######spl
      SUBROUTINE CH_UPDATE_METEO(TPM,PTIME)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #####################################
!!
!!***  *CH_UPDATE_METEO*
!!
!!    PURPOSE
!!    -------
!!    update a set of meteo variables
!!
!!**  METHOD
!!    ------
!!    interpolate NMETEOVARS variables in time
!!
!!    AUTHOR
!!    ------
!!    K. Suhre
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 21/04/95
!!    27/07/96 (K. Suhre) restructured
!!    20/04/99 (K. Suhre) meteo variables are now interpolated in time
!!    01/12/03 (D. Gazen)   change Chemical scheme interface!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_M9,       ONLY: NMETEOVARS, METEOTRANSTYPE
USE MODD_CH_METEO,   ONLY: NMETEORECS, XMETEOTIME, XMETEODATA, NMETEORECACT
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
REAL,                 INTENT(IN)    :: PTIME ! current simulation time
!
!*       0.2  declaration of local variables
!     ----------------
REAL :: ZALPHA ! interpolation weight
!
!------------------------------------------------------------------------------
!
!*    EXECUTABLE STATEMENTS
!     ---------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_UPDATE_METEO',0,ZHOOK_HANDLE)
IF (PTIME .LE. XMETEOTIME(1)) THEN
!
! take first record
  TPM%XMETEOVAR(:) = XMETEODATA(:,1)
!
ELSE IF (PTIME .GE. XMETEOTIME(NMETEORECS)) THEN
!
! take last record
  TPM%XMETEOVAR(1:NMETEOVARS) = XMETEODATA(1:NMETEOVARS,NMETEORECS)
!
ELSE
!
! interpolate meteo variables in time
  IF (PTIME .GE. XMETEOTIME(NMETEORECACT+1)) NMETEORECACT = NMETEORECACT+1
!
  ZALPHA = (PTIME                      - XMETEOTIME(NMETEORECACT)) &
         / (XMETEOTIME(NMETEORECACT+1) - XMETEOTIME(NMETEORECACT))
!
  TPM%XMETEOVAR(1:NMETEOVARS) = ZALPHA * XMETEODATA(1:NMETEOVARS,NMETEORECACT+1) &
                              + (1.-ZALPHA) * XMETEODATA(1:NMETEOVARS,NMETEORECACT)
!
END IF
!
IF (LHOOK) CALL DR_HOOK('CH_UPDATE_METEO',1,ZHOOK_HANDLE)
END SUBROUTINE CH_UPDATE_METEO
