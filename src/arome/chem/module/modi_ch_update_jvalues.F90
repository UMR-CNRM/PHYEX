!     ######spl
      MODULE MODI_CH_UPDATE_JVALUES
!!    #############################

      INTERFACE
!
      SUBROUTINE CH_UPDATE_JVALUES(KLUOUT, PZENITH, PRT, &
       PALB_UV, PZS, PZZ, PLAT0, PLON0,      &
       KLON, KLAT, KLEV,KRR,                 &
       KDAY, KMONTH, KYEAR, PTIME,           &
       OCH_TUV_ONLINE,  HCH_TUV_CLOUDS,      &
       PALBNEW, PDOBNEW, PRHODREF, PJVALUES, &
       NIB,NIE,NJB,NJE,NIU,NJU, KVERB )
!
!!    EXTERNAL
!!    --------
USE MODD_CH_INIT_JVALUES, ONLY : JPJVMAX
!
IMPLICIT NONE
INTEGER,                  INTENT(IN)   :: KLON     ! dimension I
INTEGER,                  INTENT(IN)   :: KLAT     ! dimension J
INTEGER,                  INTENT(IN)   :: KLEV     ! Number of vertical levels
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
REAL,                     INTENT(IN)   :: PLAT0, PLON0
REAL, DIMENSION(KLON,KLAT),      INTENT(IN) :: PZENITH, PALB_UV, PZS
REAL, DIMENSION(KLON,KLAT,KLEV), INTENT(IN) :: PZZ
REAL, DIMENSION(KLON,KLAT,KLEV),      INTENT(IN) :: PRHODREF
REAL, DIMENSION(KLON,KLAT,KLEV,KRR),  INTENT(IN) :: PRT
INTEGER,                   INTENT(IN) :: KDAY, KMONTH, KYEAR ! current date
REAL,                      INTENT(IN) :: PTIME    ! current time (s)
INTEGER,                   INTENT(IN) :: KLUOUT             
LOGICAL,                   INTENT(IN) :: OCH_TUV_ONLINE ! online/lookup table
CHARACTER*4,               INTENT(IN) :: HCH_TUV_CLOUDS ! clouds and radiation
REAL,                      INTENT(IN) :: PALBNEW  ! surface albedo
REAL,                      INTENT(IN) :: PDOBNEW  ! ozone column dobson
REAL,DIMENSION(KLON,KLAT,KLEV,JPJVMAX), INTENT(INOUT) :: PJVALUES    ! Tuv coefficients
INTEGER,                   INTENT(IN)    :: KVERB      ! verbosity level
INTEGER,                   INTENT(IN)    :: NIB,NIE,NJB,NJE,NIU,NJU   !  domain dim
!!
!
END SUBROUTINE CH_UPDATE_JVALUES
!
END INTERFACE
!
END MODULE MODI_CH_UPDATE_JVALUES
