!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!    ############################### 
      MODULE MODI_CH_METEO_TRANS_KESS
!!    ############################### 
!!
!
INTERFACE
!!
SUBROUTINE CH_METEO_TRANS_KESS(KL, PRHODJ, PRHODREF, PRTSM, PTHT, PABST, &
                               KVECNPT, KVECMASK, TPM, KDAY, KMONTH,     &
                               KYEAR, PLAT, PLON, PLAT0, PLON0, OUSERV,  &
                               OUSERC, OUSERR, KLUOUT, HCLOUD, PTSTEP    )
!
USE MODD_CH_M9_n, ONLY: METEOTRANSTYPE
!
IMPLICIT NONE
REAL,                     INTENT(IN), OPTIONAL    :: PTSTEP  !timestep
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD  ! Cloud parameterization
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF    ! air density
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRTSM  ! moist variables at t or t-dt or water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PABST ! theta and pressure at t
INTEGER, DIMENSION(:,:),  INTENT(IN)    :: KVECMASK
!
TYPE(METEOTRANSTYPE), DIMENSION(:), INTENT(INOUT) :: TPM
                                                   ! meteo variable for CCS
INTEGER,                  INTENT(IN)    :: KYEAR       ! Current Year
INTEGER,                  INTENT(IN)    :: KMONTH      ! Current Month
INTEGER,                  INTENT(IN)    :: KDAY        ! Current Day
INTEGER,                  INTENT(IN)    :: KLUOUT     ! channel for output listing
INTEGER,                  INTENT(IN)    :: KL, KVECNPT
REAL, DIMENSION(:,:),     INTENT(IN)    :: PLAT, PLON
REAL,                     INTENT(IN)    :: PLAT0, PLON0
LOGICAL,                  INTENT(IN)    :: OUSERV, OUSERC, OUSERR
END SUBROUTINE CH_METEO_TRANS_KESS
!!
END INTERFACE
!!
END MODULE MODI_CH_METEO_TRANS_KESS
!!
!!    ####################################################################
SUBROUTINE CH_METEO_TRANS_KESS(KL, PRHODJ, PRHODREF, PRTSM, PTHT, PABST, &
                               KVECNPT, KVECMASK, TPM, KDAY, KMONTH,     &
                               KYEAR, PLAT, PLON, PLAT0, PLON0, OUSERV,  &
                               OUSERC, OUSERR, KLUOUT, HCLOUD, PTSTEP    )
!!    ####################################################################
!!
!!*** *CH_METEO_TRANS_KESS*
!!
!!    PURPOSE
!!    -------
!     Transfer of meteorological data, such as temperature, pressure
!     and water vapor mixing ratio for one point into the variable TPM(JM+1)
!     here LWC, LWR and mean radius computed from Kessler or ICEx schemes
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
!!    01/12/07 (M. Leriche) include rain
!!    14/05/08 (M. Leriche) include raindrops and cloud droplets mean radius
!!    05/06/08 (M. Leriche) calculate LWC and LWR in coherence with time spliting scheme
!!    05/11/08 (M. Leriche) split in two routines for 1-moment and 2-moment cloud schemes
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!!
!!    EXTERNAL
!!    --------
!!    GAMMA    :  gamma function
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
USE MODD_CH_M9_n,     ONLY: NMETEOVARS,    & ! number of meteorological variables
                            METEOTRANSTYPE   !type for meteo . transfer
!!
USE MODD_CST,         ONLY: XP00,          & ! Surface pressure
                          XRD,             & ! R gas constant
                          XCPD,            & !specific heat for dry air
                          XPI,             & !pie
                          XRHOLW             !density of water
!!
USE MODD_CONF,        ONLY: LCARTESIAN     ! Logical for cartesian geometry
!!
USE MODD_RAIN_ICE_DESCR_n, ONLY: XNUC, XALPHAC,  & !Cloud droplets distrib. param.
                                 XRTMIN,         & ! min values of the water m. r.
                                 XLBC,  XLBEXC,  & !shape param. of the cloud droplets
                                 XLBR,  XLBEXR,  & !shape param. of the raindrops
                                 XCONC_LAND
!!
use mode_msg

USE MODI_GAMMA
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
REAL,                     INTENT(IN), OPTIONAL    :: PTSTEP  ! Double timestep
CHARACTER (LEN=4),        INTENT(IN)    :: HCLOUD  ! Cloud parameterization
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF    ! air density
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRTSM  ! moist variables at t or t-dt or water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT, PABST ! theta and pressure at t
INTEGER, DIMENSION(:,:),  INTENT(IN)    :: KVECMASK
!
TYPE(METEOTRANSTYPE), DIMENSION(:), INTENT(INOUT) :: TPM
                                                   ! meteo variable for CCS
INTEGER,                  INTENT(IN)    :: KYEAR       ! Current Year
INTEGER,                  INTENT(IN)    :: KMONTH      ! Current Month
INTEGER,                  INTENT(IN)    :: KDAY        ! Current Day
INTEGER,                  INTENT(IN)    :: KLUOUT     ! channel for output listing
INTEGER,                  INTENT(IN)    :: KL, KVECNPT
REAL, DIMENSION(:,:),     INTENT(IN)    :: PLAT, PLON
REAL,                     INTENT(IN)    :: PLAT0, PLON0
LOGICAL,                  INTENT(IN)    :: OUSERV, OUSERC, OUSERR
!
!*      0.2    declarations of local variables
!
REAL,DIMENSION(SIZE(PRTSM,1),SIZE(PRTSM,2),SIZE(PRTSM,3),3) :: ZRTSM
REAL,DIMENSION(SIZE(PRTSM,1),SIZE(PRTSM,2))                 :: ZLAT, ZLON
REAL,DIMENSION(SIZE(PRTSM,1),SIZE(PRTSM,2),SIZE(PRTSM,3))   :: ZRAYC, ZWLBDC, &
                                                               ZWLBDC3, ZCONC
REAL,DIMENSION(SIZE(PRTSM,1),SIZE(PRTSM,2),SIZE(PRTSM,3))   :: ZRAYR, ZWLBDR, ZWLBDR3
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
firstcall : IF (GSFIRSTCALL) THEN
!
  GSFIRSTCALL = .FALSE.
!
!*       1.1   check if number of variables NMETEOVARS
!              corresponds to what the CCS expects
!
  IF (NMETEOVARS /= 13) THEN
    WRITE(KLUOUT,*) "CH_METEO_TRANS ERROR: number of meteovars to transfer"
    WRITE(KLUOUT,*) "does not correspond to the number expected by the CCS:"
    WRITE(KLUOUT,*) "     meteovars to transfer: ", 13
    WRITE(KLUOUT,*) "     NMETEOVARS expected:   ", NMETEOVARS
    WRITE(KLUOUT,*) "Check the definition of NMETEOVARS in your .chf file."
    WRITE(KLUOUT,*) "The program will be stopped now!"
    call Print_msg( NVERB_FATAL, 'GEN', 'CH_METEO_TRANS_KESS', &
                    'number of meteovars to transfer does not correspond to the expected number.' )
  END IF
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
  TPM(:)%CMETEOVAR(11)= "Rain water   (kg/kg)"
  TPM(:)%CMETEOVAR(12)= "Mean cloud droplets radius (m)"
  TPM(:)%CMETEOVAR(13)= "Mean raindrops radius      (m)"
!
ENDIF firstcall
!
! "Water vapor (kg/kg)"
!
IF (OUSERV) THEN
! if split option, use tendency
  IF (PRESENT(PTSTEP)) THEN
    ZRTSM(:,:,:,1) = (PRTSM(:,:,:, 1)/ PRHODJ(:,:,:))*PTSTEP
  ELSE
    ZRTSM(:,:,:,1) = PRTSM(:,:,:, 1)
  ENDIF
ELSE
  ZRTSM(:,:,:,1) = 0.0
ENDIF
!
! "Cloud water (kg/kg)" and  "Mean cloud droplets radius (m)"
!
IF (OUSERC) THEN
  IF (PRESENT(PTSTEP)) THEN    
    ZRTSM(:,:,:,2) = (PRTSM(:,:,:, 2)/ PRHODJ(:,:,:))*PTSTEP
  ELSE
    ZRTSM(:,:,:,2) = PRTSM(:,:,:, 2)
  ENDIF
  ZRAYC(:,:,:) = 10.e-6 ! avoid division by zero
  SELECT CASE (HCLOUD)
  CASE ('KESS')
    WHERE (ZRTSM(:,:,:, 2)>1.e-20) !default value for Kessler
      ZRAYC(:,:,:) = 10.e-6 ! assume a cloud droplet radius of 10 µm
    ENDWHERE
  CASE ('ICE3','ICE4')
    WHERE (ZRTSM(:,:,:, 2)>XRTMIN(2))
      ZCONC(:,:,:) = XCONC_LAND
      ZWLBDC3(:,:,:) = XLBC(1) * ZCONC(:,:,:) / (PRHODREF(:,:,:) * ZRTSM(:,:,:, 2))
      ZWLBDC(:,:,:)  = ZWLBDC3(:,:,:)**XLBEXC
      ZRAYC(:,:,:) = 0.5*GAMMA(XNUC+1./XALPHAC)/(GAMMA(XNUC)*ZWLBDC(:,:,:))
!      ZRAYC(:,:,:) = 10.e-6 ! assume a cloud droplet radius of 10 µm
    ENDWHERE
  END SELECT
ELSE
  ZRTSM(:,:,:,2) = 0.0
  ZRAYC(:,:,:) = 10.e-6 ! avoid division by zero
ENDIF
!
! "Rain water (kg/kg)" and "Mean raindrops radius (m)"
!
IF (OUSERR) THEN
  IF (PRESENT(PTSTEP)) THEN
    ZRTSM(:,:,:,3) = (PRTSM(:,:,:, 3)/ PRHODJ(:,:,:))*PTSTEP
  ELSE
    ZRTSM(:,:,:,3) = PRTSM(:,:,:, 3)
  ENDIF
  ZRAYR(:,:,:) = 500.e-6 ! avoid division by zero
  SELECT CASE (HCLOUD)
  CASE ('KESS')
    WHERE (ZRTSM(:,:,:, 3)>1.e-20) !default value for Kessler
      ZRAYR(:,:,:) = 0.5*((XPI*XRHOLW*1.E7)/ &
                     (PRHODREF(:,:,:)*ZRTSM(:,:,:,3)))**(-1./4.)
    ENDWHERE
  CASE ('ICE3','ICE4')
    WHERE (ZRTSM(:,:,:, 3)>XRTMIN(3))
      ZRAYR(:,:,:) = 0.5*(1./(XLBR*(PRHODREF(:,:,:)*ZRTSM(:,:,:,3))**XLBEXR))
    ENDWHERE
  END SELECT
ELSE
  ZRTSM(:,:,:,3) = 0.0
  ZRAYR(:,:,:) = 500.e-6 ! avoid division by zero
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
!!
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
!
  TPM(JM+1)%XMETEOVAR(1) = JK-1 ! assuming first model level is level 2
!   TPM(JM+1)%XMETEOVAR(1) = JK ! assuming first model level is level 1
!
! "Air density (kg/m3)"
!
  TPM(JM+1)%XMETEOVAR(2) = PRHODREF(JI, JJ, JK)
!
! "Temperature (K)"
!
  TPM(JM+1)%XMETEOVAR(3) = PTHT(JI,JJ,JK)*((PABST(JI,JJ,JK)/XP00)**(XRD/XCPD))
!
! "Water vapor (kg/kg)"
!
  TPM(JM+1)%XMETEOVAR(4) = ZRTSM(JI, JJ, JK, 1)
!
! "Cloud water (kg/kg)"
!
  TPM(JM+1)%XMETEOVAR(5) = ZRTSM(JI, JJ, JK, 2)
!
!  "Latitude (rad)"
!
  TPM(JM+1)%XMETEOVAR(6) = ZLAT(JI, JJ)
!
!  "Longitude (rad)"
!
  TPM(JM+1)%XMETEOVAR(7) = ZLON(JI, JJ)
!
!  "Current date"
!
  TPM(JM+1)%XMETEOVAR(8) = REAL(KYEAR)
  TPM(JM+1)%XMETEOVAR(9) = REAL(KMONTH)
  TPM(JM+1)%XMETEOVAR(10)= REAL(KDAY)
!
! "Rain water (kg/kg)"
!
  TPM(JM+1)%XMETEOVAR(11) = ZRTSM(JI, JJ, JK, 3)
!
! "Mean cloud droplets radius (m)"
!
  TPM(JM+1)%XMETEOVAR(12) = ZRAYC(JI, JJ, JK)
!
! "Mean raindrops radius (m)"
!
  TPM(JM+1)%XMETEOVAR(13) = ZRAYR(JI, JJ, JK)
!
ENDDO
!
END SUBROUTINE CH_METEO_TRANS_KESS
