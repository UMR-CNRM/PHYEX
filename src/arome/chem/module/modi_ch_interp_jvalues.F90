!     ######spl
       MODULE MODI_CH_INTERP_JVALUES
!      #############################
!
       INTERFACE 
!
       SUBROUTINE CH_INTERP_JVALUES(PJVALUES,PZZ, PZS, PALBUV, PZENITH,  &
                                      KLON, KLAT, KLEV,                  &
                                      NIB,NIE,NJB,NJE,NIU,NJU, KVERB, KLUOUT)
!
USE MODD_CH_INIT_JVALUES, ONLY : JPJVMAX
       IMPLICIT NONE
       INTEGER,                  INTENT(IN) :: KLUOUT  
       INTEGER,                  INTENT(IN)   :: KLON     ! dimension I
       INTEGER,                  INTENT(IN)   :: KLAT     ! dimension J
       INTEGER,                  INTENT(IN)   :: KLEV     ! Number of vertical levels
       REAL, DIMENSION(KLON,KLAT),       INTENT(IN) :: PZENITH, PALBUV, PZS
       REAL, DIMENSION(KLON,KLAT,KLEV),  INTENT(IN) :: PZZ
       REAL, DIMENSION(KLON,KLAT,KLEV,JPJVMAX), INTENT(INOUT) :: PJVALUES
       INTEGER,                  INTENT(IN)    :: NIB,NIE,NJB,NJE,NIU,NJU   !  domain dim
       INTEGER,                  INTENT(IN)    :: KVERB      ! verbosity level
       END SUBROUTINE CH_INTERP_JVALUES
       END INTERFACE 
!
       END MODULE MODI_CH_INTERP_JVALUES
