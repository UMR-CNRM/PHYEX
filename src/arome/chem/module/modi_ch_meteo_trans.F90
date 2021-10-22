!     ######spl
      MODULE MODI_CH_METEO_TRANS
!!    ############################# 
!!
!
INTERFACE
!!
SUBROUTINE CH_METEO_TRANS(KL, PRHODREF, PRT, PTHT, PABST,  &
                            KVECNPT, KVECMASK, TPM, KDAY,&
                            KMONTH, KYEAR, PLAT, PLON, PLAT0, PLON0,&
                            OUSERV, OUSERC, KLUOUT )

USE MODD_CH_M9, ONLY: METEOTRANSTYPE
IMPLICIT NONE
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF    ! air density
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT         ! moist variables at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PABST ! theta and pressure at t
INTEGER, DIMENSION(:,:), INTENT(IN) :: KVECMASK
TYPE(METEOTRANSTYPE), DIMENSION(:), INTENT(INOUT) :: TPM
                                                   ! meteo variable for CCS
INTEGER,              INTENT(IN)    :: KYEAR       ! Current Year
INTEGER,              INTENT(IN)    :: KMONTH      ! Current Month
INTEGER,              INTENT(IN)    :: KDAY        ! Current Day
INTEGER,              INTENT(IN)    :: KLUOUT     ! channel for output listing
INTEGER,              INTENT(IN)    :: KL, KVECNPT
REAL, DIMENSION(:,:), INTENT(IN)    :: PLAT, PLON
REAL,                 INTENT(IN)    :: PLAT0, PLON0
LOGICAL,              INTENT(IN)    :: OUSERV, OUSERC

END SUBROUTINE CH_METEO_TRANS
!!
END INTERFACE
!!
END MODULE MODI_CH_METEO_TRANS
