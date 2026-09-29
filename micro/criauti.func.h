ELEMENTAL SUBROUTINE CRIAUTI(PRCRIAUTI, PT0CRIAUTI, PBCRIAUTI, PACRIAUTI)
!!
!!**  PURPOSE
!!    -------
!!      Meaningful text
!!
!!    AUTHOR
!!    ------
!!      U. Andrae (2026)
!!
!!    MODIFICATIONS
!!    -------------
!!
!
!
!*      0. DECLARATIONS
!          ------------
!
IMPLICIT NONE
!
REAL,    INTENT(IN)    :: PRCRIAUTI   ! Constants for pristine ice autoconversion
REAL,    INTENT(IN)    :: PT0CRIAUTI  ! Temp degC at which cirrus law starts to be used
REAL,    INTENT(OUT)   :: PACRIAUTI   ! A Coef. for cirrus law 
REAL,    INTENT(OUT)   :: PBCRIAUTI   ! B Coef. for cirrus law 
!
REAL :: ZTCRI0, ZCRI0
!
! Second point to determine 10**(aT+b) law
!
ZTCRI0=-40.0
ZCRI0=1.25E-6
!
PBCRIAUTI=-( LOG10(PRCRIAUTI) - LOG10(ZCRI0)*PT0CRIAUTI/ZTCRI0 )&
            *ZTCRI0/(PT0CRIAUTI-ZTCRI0)
PACRIAUTI=(LOG10(ZCRI0)-PBCRIAUTI)/ZTCRI0
!
RETURN
!
END SUBROUTINE CRIAUTI
