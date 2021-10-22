!     ######spl
      SUBROUTINE CH_METEO_TRANS(KL, PRHODREF, PRT, PTHT, PABST,  &
                                   KVECNPT, KVECMASK, TPM, KDAY,&
                                   KMONTH, KYEAR, PLAT, PLON, PLAT0, PLON0,&
                                   OUSERV, OUSERC, KLUOUT )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK

!!    ###############################################################
!!
!!*** *CH_METEO_TRANS*
!!
!!    PURPOSE
!!    -------
!     Transfer of meteorological data, such as temperature, pressure
!     and water vapor mixing ratio for one point into the variable TPM(JM+1)
!!
!!    METHOD
!!    ------
!!    For the given grid-point KI,KJ,KK, the meteorological parameters
!!    will be transfered for use by CH_SET_RATES and CH_SET_PHOTO_RATES.
!!    Presently, the variables altitude, air density, temperature,
!!    water vapor mixing ratio, cloud water, longitude, latitude and date
!!    will be transfered. In the chemical definition file (.chf)
!!    these variables have to be transfered into variables like O2, H2O etc.
!!    Also, consistency is checked between the number of
!!    variables expected by the CCS (as defined in the .chf file) and
!!    the number of variables to be transfered here. If you change
!!    the meaning of XMETEOVARS in your .chf file, make sure to modify
!!    this subroutine accordingly.
!!      If the model is run in 1D mode, the model level instead of altitude
!!    is passed. In 2D and 3D, altitude is passed with a negative sign
!!    so that the radiation scheme TUV can make the difference between
!!    model levels and altitude.
!!
!!    AUTHOR
!!    ------
!!    K. Suhre    *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 24/05/95
!!    04/08/96 (K. Suhre) restructured
!!    21/02/97 (K. Suhre) add XLAT0 and XLON0 for LCARTESIAN=T case
!!    27/08/98 (P. Tulet) add temperature at t for kinetic coefficient
!!    09/03/99 (V. Crassier & K. Suhre) vectorization
!!    09/03/99 (K. Suhre) modification for TUV
!!    09/03/99 (C. Mari & J. Escobar) Code optimization 
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!    01/12/03 (D. Gazen)   change Chemical scheme interface
!!    01/12/04 (P. Tulet)   update ch_meteo_transn.f90 for Arome
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
USE MODD_CH_M9,       ONLY: NMETEOVARS,    &! number of meteorological variables
                            METEOTRANSTYPE  !type for meteo . transfer
!!
USE MODD_CST,       ONLY: XP00,      &! Surface pressure
                          XRD,       &! R gas constant
                          XCPD        !specific heat for dry air
!!
USE MODD_CONF,      ONLY: LCARTESIAN  ! Logical for cartesian geometry
!!
!!
!!
!
 
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
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
REAL,                 INTENT(IN)    :: PLAT0, PLON0
LOGICAL,              INTENT(IN)    :: OUSERV, OUSERC
REAL, DIMENSION(:,:), INTENT(IN)    :: PLAT, PLON
 
 
!
!*      0.2    declarations of local variables
!
REAL,DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),2) :: ZRT
REAL,DIMENSION(SIZE(PRT,1),SIZE(PRT,2)) :: ZLAT, ZLON
LOGICAL, SAVE :: GSFIRSTCALL = .TRUE.
INTEGER :: JI,JJ,JK,JM
INTEGER :: IDTI,IDTJ,IDTK
!
!
!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZE METEO VARIABLE TRANSFER
!              ----------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_METEO_TRANS',0,ZHOOK_HANDLE)
firstcall : IF (GSFIRSTCALL) THEN
!
  GSFIRSTCALL = .FALSE.
!
!*       1.1   check if number of variables NMETEOVARS
!              corresponds to what the CCS expects
!
  IF (NMETEOVARS /= 10) THEN
    WRITE(KLUOUT,*) "CH_METEO_TRANS ERROR: number of meteovars to transfer"
    WRITE(KLUOUT,*) "does not correspond to the number expected by the CCS:"
    WRITE(KLUOUT,*) "     meteovars to transfer: ", 10
    WRITE(KLUOUT,*) "     NMETEOVARS expected:   ", NMETEOVARS
    WRITE(KLUOUT,*) "Check the definition of NMETEOVARS in your .chf file."
    WRITE(KLUOUT,*) "The program will be stopped now!"
    STOP 1
  END IF
!
!
!*       1.2   initialize names of meteo vars
!
  TPM(:)%CMETEOVAR(1) = "Model level"
  TPM(:)%CMETEOVAR(2) = "Air density  (kg/m3)"
  TPM(:)%CMETEOVAR(3) = "Temperature  (K)"
  TPM(:)%CMETEOVAR(4) = "Water vapor  (kg/kg)"
  TPM(:)%CMETEOVAR(5) = "Cloud water  (kg/kg)"
  TPM(:)%CMETEOVAR(6) = "Latitude     (rad)"
  TPM(:)%CMETEOVAR(7) = "Longitude    (rad)"
  TPM(:)%CMETEOVAR(8) = "Current date (year)"
  TPM(:)%CMETEOVAR(9) = "Current date (month)"
  TPM(:)%CMETEOVAR(10)= "Current date (day)"
!
ENDIF firstcall
!
! "Water vapor (kg/kg)"
  IF (OUSERV) THEN
    ZRT(:,:,:,1) = PRT(:,:,:, 1)
  ELSE
    ZRT(:,:,:,1) = 0.0
  ENDIF
!
! "Cloud water (kg/kg)"
  IF (OUSERC) THEN
    ZRT(:,:,:,2) = PRT(:,:,:, 2)
  ELSE
    ZRT(:,:,:,2) = 0.0
  ENDIF

  IF(LCARTESIAN) THEN
!  "Latitude (rad)"
    ZLAT(:,:) = PLAT0
!  "Longitude (rad)"
    ZLON(:,:) = PLON0
  ELSE
!  "Latitude (rad)"
    ZLAT(:,:) = PLAT(:,:)
!  "Longitude (rad)"
    ZLON(:,:) = PLON(:,:)
  END IF
!
!
!*       2.    TRANSFER METEO VARIABLES
!              ------------------------
!
IDTI=KVECMASK(2,KL)-KVECMASK(1,KL)+1
IDTJ=KVECMASK(4,KL)-KVECMASK(3,KL)+1
IDTK=KVECMASK(6,KL)-KVECMASK(5,KL)+1
!Vectorization:
!ocl novrec
!cdir nodep
DO JM=0,KVECNPT-1
  JI=JM-IDTI*(JM/IDTI)+KVECMASK(1,KL)
  JJ=JM/IDTI-IDTJ*(JM/(IDTI*IDTJ))+KVECMASK(3,KL)
  JK=JM/(IDTI*IDTJ)-IDTK*(JM/(IDTI*IDTJ*IDTK))+KVECMASK(5,KL)
!
!"Model Altitude"
  TPM(JM+1)%XMETEOVAR(1) = JK-1 ! assuming first model level is level 2
!   TPM(JM+1)%XMETEOVAR(1) = JK ! assuming first model level is level 1
!
! "Air density (kg/m3)"
  TPM(JM+1)%XMETEOVAR(2) = PRHODREF(JI, JJ, JK)
!
! "Temperature (K)"
  TPM(JM+1)%XMETEOVAR(3) = PTHT(JI,JJ,JK)*((PABST(JI,JJ,JK)/XP00)**(XRD/XCPD))
!
! "Water vapor (kg/kg)"
    TPM(JM+1)%XMETEOVAR(4) = ZRT(JI, JJ, JK, 1)
!
! "Cloud water (kg/kg)"
    TPM(JM+1)%XMETEOVAR(5) = ZRT(JI, JJ, JK, 2)
!
!  "Latitude (rad)"
    TPM(JM+1)%XMETEOVAR(6) = ZLAT(JI, JJ)
!  "Longitude (rad)"
    TPM(JM+1)%XMETEOVAR(7) = ZLON(JI, JJ)
!
!  "Current date"
  TPM(JM+1)%XMETEOVAR(8) = FLOAT(KYEAR)
  TPM(JM+1)%XMETEOVAR(9) = FLOAT(KMONTH)
  TPM(JM+1)%XMETEOVAR(10)= FLOAT(KDAY)
!
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CH_METEO_TRANS',1,ZHOOK_HANDLE)
END SUBROUTINE CH_METEO_TRANS
