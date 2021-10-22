!     ######spl
     SUBROUTINE DUST_FILTER(PSV, PRHODREF)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   #######################################
!!
!!   PURPOSE
!!   -------
!!
!!   REFERENCE
!!   ---------
!!   none
!!
!!   AUTHOR
!!    ------
!!    Pierre TULET (CNRM/GMEI) 
!!
!!   MODIFICATIONS
!!    -------------
!!   Original
!!
! Entry variables:
!
! PRSVS(INOUT)       -Array of moments included in PRSVS
!
!*************************************************************
! Exit variables:
!
!*************************************************************
! Variables used during the deposition velocity calculation
! 
! ZVGK       -Polydisperse settling velocity of the kth moment (m/s)
!************************************************************
!!
!!   IMPLICIT ARGUMENTS
!
USE MODD_DUST
USE MODD_CSTS_DUST
!!
IMPLICIT NONE

!! Declarations d'arguments

REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSV
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF

!! Declarations de variables internes

INTEGER :: JN
INTEGER :: JMODEIDX
REAL,    DIMENSION(NMODE_DST*3) :: ZPMIN
REAL,    DIMENSION(NMODE_DST)   :: ZINIRADIUS
REAL,    DIMENSION(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_DST*3)  :: ZM                  ! [aerosol units] local array which goes to output later
REAL,    DIMENSION(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3))  :: ZSIGMA !standard deviation

REAL :: ZRGMIN, ZSIGMIN
REAL :: ZPI, ZRHOP, ZFAC, ZMI
INTEGER,DIMENSION(NMODE_DST) :: NM0, NM3, NM6

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DUST_FILTER',0,ZHOOK_HANDLE)
PSV(:,:,:,:) =  MAX(PSV(:,:,:,:), 1.E-80)


IF (LHOOK) CALL DR_HOOK('DUST_FILTER',1,ZHOOK_HANDLE)
END SUBROUTINE DUST_FILTER
