!     ######spl
     SUBROUTINE SALT_FILTER(PSV, PRHODREF)
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
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!   IMPLICIT ARGUMENTS
!
USE MODD_SALT
USE MODD_CSTS_SALT
!
USE MODE_SALT_PSD
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,  DIMENSION(:,:,:,:),  INTENT(INOUT) :: PSV
REAL,  DIMENSION(:,:,:),    INTENT(IN)    :: PRHODREF
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JN
INTEGER :: IMODEIDX
REAL,    DIMENSION(NMODE_SLT*3) :: ZPMIN
REAL,    DIMENSION(NMODE_SLT)   :: ZINIRADIUS
REAL,    DIMENSION(SIZE(PSV,1), SIZE(PSV,2), SIZE(PSV,3), NMODE_SLT*3)  :: ZM                  ! [aerosol units] local array which goes to output later

REAL :: ZRGMIN, ZSIGMIN
REAL :: ZRHOP, ZMI
INTEGER,DIMENSION(NMODE_SLT) :: NM0, NM3, NM6
!
!*       0.3  initialize constant
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SALT_FILTER',0,ZHOOK_HANDLE)
ZRHOP = XDENSITY_SALT
ZMI   = XMOLARWEIGHT_SALT ! molecular mass in kg/mol
!
!-------------------------------------------------------------------------------

!
PSV(:,:,:,:) =  MAX(PSV(:,:,:,:), 1.E-80)
!

IF (LHOOK) CALL DR_HOOK('SALT_FILTER',1,ZHOOK_HANDLE)
END SUBROUTINE SALT_FILTER
