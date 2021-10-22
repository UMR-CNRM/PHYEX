!     ######spl
     SUBROUTINE CH_AER_SURF( PM, PRG, PSIG,  PSURF)
     USE PARKIND1, ONLY : JPRB
     USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!   ############################################################
!!
!!    PURPOSE
!!    -------
!!      Compute surface aerosol 
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    P. Tulet (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!    EXTERNAL
!!    --------
!!    None
!!
USE MODD_CH_AEROSOL, ONLY : XPI, JPMODE, NM0, NM3, NM6
!
IMPLICIT NONE
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
!*      0.1    declarations of arguments
!
REAL,   DIMENSION(:,:), INTENT(IN) :: PM      ! moments
REAL,   DIMENSION(:,:), INTENT(IN) :: PRG     ! radius
REAL,   DIMENSION(:,:), INTENT(IN) :: PSIG    ! dispersion
REAL,   DIMENSION(:,:), INTENT(OUT) :: PSURF  ! aerosol surface
!
!
!*      0.2    declarations local variables
!
REAL,   DIMENSION(SIZE(PM,1),JPMODE) :: ZSIG0, ZRG0, ZN0
REAL,   DIMENSION(SIZE(PM,1),SIZE(PM,2)) :: ZM
INTEGER :: JN
!
!-------------------------------------------------------------------------------
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_AER_SURF',0,ZHOOK_HANDLE)
DO JN=1, JPMODE
  ZM(:,NM0(JN)) = PM(:,NM0(JN)) ! 1/m3
  ZM(:,NM3(JN)) = PM(:,NM3(JN)) ! 
  ZM(:,NM6(JN)) = PM(:,NM6(JN)) 

  ZN0(:,JN) = ZM(:,NM0(JN)) * 1E-6 ! convert from 1/m3 to 1/cc
  ZRG0(:,JN)= PRG(:,JN) * 1E-6   ! convert from micrometers to meters 
  ZSIG0(:,JN)= PSIG(:,JN)
!
! Surface = Pi * M2 with M2 = N Rg**2 exp (2*(ln(sigma)**2))
!
!surface in m2 / cc (psurf will be convert into SI units in BASIC)
  PSURF(:,JN) =  XPI * ZN0(:,JN) * (ZRG0(:,JN))**2 * EXP(2.*ZSIG0(:,JN)**2)

ENDDO
!
!
IF (LHOOK) CALL DR_HOOK('CH_AER_SURF',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AER_SURF
