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
REAL,    INTENT(IN)    :: PRCRIAUTI   ! Explanation
REAL,    INTENT(IN)    :: PT0CRIAUTI  ! Explanation
REAL,    INTENT(OUT)   :: PACRIAUTI   ! Explanation
REAL,    INTENT(OUT)   :: PBCRIAUTI   ! Explanation
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
