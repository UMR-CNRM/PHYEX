!     ######spl
       SUBROUTINE CH_INIT_SCHEME(KLUOUT)
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!      ###################################
!!
!!***  *CH_INIT_SCHEME*
!!
!!    PURPOSE
!!    -------
!      initialize module MODD_CH_M9 variables (now MesoNH variables) from internal 
!      scheme constants defined in module MODD_CH_M9_SCHEME from BASIC.f90 (scheme dependant).
!!
!!**  METHOD
!!    ------
!!    Names of variables are identical in both MODD_CH_M9 and MODD_CH_M9_SCHEME
!!    modules to minimize sources modifications and keep users habits.
!!    So we rename variables/constant from MODD_CH_M9_SCHEME by prefixing them with
!!    I_ (I means scheme Internal).
!!
!!==============================================================================
!!
!!    REFERENCES
!!    ----------
!!    MesoNH-chemistry book 3
!!
!!    AUTHOR
!!    ------
!!    D. Gazen
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 19/10/03
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!    None
!----------------------------------------------------------------------------
!
USE MODD_CH_M9 !! the mesonh interface for chemical variables
USE MODD_CH_M9_SCHEME, ONLY :  TACCS,  &
     & I_NEQ=>NEQ,                     &
     & I_NREAC=>NREAC,                 &
     & I_NMETEOVARS=>NMETEOVARS,       &
     & I_NNONZEROTERMS=>NNONZEROTERMS, &
     & I_CNAMES=>CNAMES,               &          
     & I_CREACS=>CREACS,               &
     & I_CFULLREACS=>CFULLREACS,       &
     & JP_CO
!
IMPLICIT NONE
!
!*       0.   DECLARATIONS
!        -----------------
!
INTEGER, INTENT(IN) :: KLUOUT
!
!----------------------------------------------------------------------------
!
!*       1.   INITIALISATION
!        -------------------
!Left member belongs to MODD_CH_M9 module variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_INIT_SCHEME',0,ZHOOK_HANDLE)
NEQ           = I_NEQ
NREAC         = I_NREAC
NMETEOVARS    = I_NMETEOVARS
NNONZEROTERMS = I_NNONZEROTERMS
! 
CNAMES=>I_CNAMES
CREACS=>I_CREACS
CFULLREACS=>I_CFULLREACS
!
!
WRITE(KLUOUT,*) 'CH_INIT_SCHEME done : NEQ, NREAC, NMETEOVARS, NNONZEROTERMS = ',&
     & NEQ, NREAC, NMETEOVARS, NNONZEROTERMS

IF (LHOOK) CALL DR_HOOK('CH_INIT_SCHEME',1,ZHOOK_HANDLE)
END SUBROUTINE CH_INIT_SCHEME
