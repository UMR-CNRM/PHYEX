!MNH_LIC Copyright 2000-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Author: Valery Masson (Meteo-France *)
!     Original 15/05/2000
! Modifications:
!  G. Jaubert  19/04/2001: add CVBALL type
!  P. Lacarrere   03/2008: add 3D fluxes
!  M. Leriche  12/12/2008: move ZTDIST out from if.not.(tpflyer%fly)
!  V. Masson   15/12/2008: correct do while aircraft move
!  O. Caumont     03/2013: add radar reflectivities
!  C. Lac         04/2014: allow RARE calculation only if CCLOUD=ICE3
!  O. Caumont     05/2014: modify RARE for hydrometeors containing ice + add bright band calculation for RARE
!  C. Lac 02/2015: correction to prevent aircraft crash
!  O. Nuissier/F. Duffourg 07/2015: add microphysics diagnostic for aircraft, ballon and profiler
!  G. Delautier   10/2016: LIMA
!  P. Wautelet 28/03/2018: replace TEMPORAL_DIST by DATETIME_DISTANCE
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  P. Wautelet 01/10/2020: bugfix: initialize GSTORE
!  P. Wautelet 14/01/2021: bugfixes: -ZXCOEF and ZYCOEF were not computed if CVBALL
!                                    -PCIT was used if CCLOUD/=ICEx (not allocated)
!                                    -PSEA was always used even if not allocated (CSURF/=EXTE)
!                                    -do not use PMAP if cartesian domain
!  P. Wautelet    06/2022: reorganize flyers
!  P. Wautelet 01/06/2023: deduplicate code => moved to modd/mode_sensors.f90
!-----------------------------------------------------------------
!      ##########################
MODULE MODE_AIRCRAFT_BALLOON_EVOL
!      ##########################

USE MODE_MSG

IMPLICIT NONE

PRIVATE

PUBLIC :: AIRCRAFT_BALLOON_EVOL

PUBLIC :: AIRCRAFT_COMPUTE_POSITION

PUBLIC :: FLYER_GET_RANK_MODEL_ISCRASHED

CONTAINS
!     ########################################################
      SUBROUTINE AIRCRAFT_BALLOON_EVOL(PTSTEP,               &
                       PZ, PMAP, PLONOR, PLATOR,             &
                       PU, PV, PW, PP, PTH, PR, PSV, PTKE,   &
                       PTS, PRHODREF, PCIT, TPFLYER,         &
                       KRANK_CUR, KRANK_NXT, PSEA            )
!     ########################################################
!
!
!!****  *AIRCRAFT_BALLOON_EVOL* - (advects and) stores
!!                                balloons/aircrafts in the model
!!
!!    PURPOSE
!!    -------
!
!
!!**  METHOD
!!    ------
!!
!! 1) All the balloons are tested. If the current balloon is
!!     a) in the current model
!!     b) not crashed
!!   the following computations are done.
!!
!! 2) The balloon position is computed.
!!       Interpolations at balloon positions are performed according to mass
!! points (because density is computed here for iso-density balloons).
!! Therefore, all model variables are used at mass points. Shuman averaging
!! are performed on X, Y, Z, U, V, W.
!!
!! 3) Storage of balloon data
!!       If storage is asked for this time-step, the data are recorded in the
!! balloon time-series.
!!
!! 4) Balloon advection
!!       If the balloon is launched, it is advected according its type
!!    a) iso-density balloons are advected following horizontal wind.
!!          the slope of the iso-density surfaces is neglected.
!!    b) radio-sounding balloons are advected according to all wind velocities.
!!          the vertical ascent speed is added to the vertical wind speed.
!!    c) Constant Volume balloons are advected according to all wind velocities.
!!          the vertical ascent speed is computed using the balloon equation
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_AIRCRAFT_BALLOON
USE MODD_CST,              ONLY: XCPD, XLVTT
USE MODD_IO,               ONLY: ISP
USE MODD_TIME_n,           ONLY: TDTCUR
USE MODD_TURB_FLUX_AIRCRAFT_BALLOON, ONLY: XRCW_FLUX, XSVW_FLUX, XTHW_FLUX
!
USE MODE_DATETIME
USE MODE_NEST_ll,          ONLY: GET_MODEL_NUMBER_ll
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
REAL,                     INTENT(IN)     :: PTSTEP ! time step
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PZ     ! z array
REAL, DIMENSION(:,:),     INTENT(IN)     :: PMAP   ! map factor
REAL,                     INTENT(IN)     :: PLONOR ! origine longitude
REAL,                     INTENT(IN)     :: PLATOR ! origine latitude
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PU     ! horizontal wind X component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PV     ! horizontal wind Y component
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PW     ! vertical wind
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PP     ! pressure
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTH    ! potential temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PR     ! water mixing ratios
REAL, DIMENSION(:,:,:,:), INTENT(IN)     :: PSV    ! Scalar variables
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PTKE   ! turbulent kinetic energy
REAL, DIMENSION(:,:),     INTENT(IN)     :: PTS    ! surface temperature
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PRHODREF ! dry air density of the reference state
REAL, DIMENSION(:,:,:),   INTENT(IN)     :: PCIT     ! pristine ice concentration
!
CLASS(TFLYERDATA),        INTENT(INOUT)  :: TPFLYER ! balloon/aircraft
INTEGER,                  INTENT(IN)     :: KRANK_CUR
INTEGER,                  INTENT(OUT)    :: KRANK_NXT
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER :: IMI        ! model index
INTEGER :: IKB        ! vertical domain sizes
INTEGER :: IKE
INTEGER :: IKU
!
REAL, DIMENSION(2,2,SIZE(PZ,3))     :: ZZM    ! mass point coordinates
REAL, DIMENSION(2,2,SIZE(PZ,3))     :: ZZU    ! U points z coordinates
REAL, DIMENSION(2,2,SIZE(PZ,3))     :: ZZV    ! V points z coordinates
REAL, DIMENSION(2,2,SIZE(PZ,3))     :: ZWM    ! mass point wind
!
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZEXN   ! Exner function
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZTH_EXN ! potential temperature multiplied by Exner function
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZRHO   ! air density
REAL                                :: ZFLYER_EXN ! balloon/aircraft Exner func.
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZTHW_FLUX  !
REAL, DIMENSION(2,2,SIZE(PTH,3))    :: ZRCW_FLUX  !
REAL, DIMENSION(2,2,SIZE(PSV,3),SIZE(PSV,4)) :: ZSVW_FLUX
!
LOGICAL :: GLAUNCH  ! launch/takeoff is effective at this time-step (if true)
LOGICAL :: GOWNER_CUR ! The process is the current owner of the flyer
!
INTEGER :: II_M     ! mass balloon position (x index)
INTEGER :: IJ_M     ! mass balloon position (y index)
INTEGER :: II_U     ! U flux point balloon position (x index)
INTEGER :: IJ_V     ! V flux point balloon position (y index)
!
INTEGER :: ISTORE          ! time index for storage
!
REAL    :: ZTSTEP
TYPE(DATE_TIME) :: TZNEXT ! Time for next position
!----------------------------------------------------------------------------
IKU = SIZE(PZ,3)

CALL GET_MODEL_NUMBER_ll(IMI)

! Set initial value for KRANK_NXT
! It needs to be 0 on all processes except the one where it is when this subroutine is called
! If the flyer flies to an other process, KRANK_NXT will be set accordingly by the current owner
IF ( TPFLYER%NRANK_CUR == ISP ) THEN
  GOWNER_CUR = .TRUE. ! This variable is set and used because NRANK_CUR could change in this subroutine
  KRANK_NXT = ISP
ELSE
  GOWNER_CUR = .FALSE.
  KRANK_NXT = 0
END IF

SELECT TYPE ( TPFLYER )
  CLASS IS ( TAIRCRAFTDATA)
    ! Take-off?
    TAKEOFF: IF ( .NOT. TPFLYER%LTOOKOFF ) THEN
      ! Do the take-off positioning only once
      ! (on model 1 for 'MOB', if aircraft is on an other model, data will be available on the right one anyway)
      IF (      ( TPFLYER%CMODEL == 'MOB' .AND. IMI == 1              ) &
           .OR. ( TPFLYER%CMODEL == 'FIX' .AND. IMI == TPFLYER%NMODEL ) ) THEN
        ! Is the aircraft in flight ?
        IF ( TDTCUR >= TPFLYER%TLAUNCH .AND. TDTCUR <= TPFLYER%TLAND ) THEN
          TPFLYER%LFLY     = .TRUE.
          TPFLYER%LTOOKOFF = .TRUE.
        END IF
      END IF
    END IF TAKEOFF

    !Do we have to store aircraft data?
    IF ( IMI == TPFLYER%NMODEL ) THEN
      TPFLYER%LSTORE = TPFLYER%TFLYER_TIME%STORESTEP_CHECK_AND_SET( ISTORE )
      IF ( TPFLYER%LSTORE ) TPFLYER%NSTORE_CUR = ISTORE
    END IF


    ! For aircrafts, data has only to be computed at store moments
    IF ( IMI == TPFLYER%NMODEL .AND. TPFLYER%LFLY .AND. TPFLYER%LSTORE ) THEN
      ! Check if it is the right moment to store data
      IF ( ABS( TDTCUR - TPFLYER%TFLYER_TIME%TPDATES(ISTORE) ) < 1e-10 ) THEN
        ISOWNERAIR: IF ( TPFLYER%NRANK_CUR == ISP ) THEN
          CALL FLYER_INTERP_TO_MASSPOINTS()

          ZEXN(:,:,:) = FLYER_COMPUTE_EXNER( )
          ZRHO(:,:,:) = FLYER_COMPUTE_RHO( )

          ZTHW_FLUX(:,:,:) = ZRHO(:,:,:)*XCPD *XTHW_FLUX(II_M:II_M+1,IJ_M:IJ_M+1,:)
          ZRCW_FLUX(:,:,:) = ZRHO(:,:,:)*XLVTT*XRCW_FLUX(II_M:II_M+1,IJ_M:IJ_M+1,:)
          ZSVW_FLUX(:,:,:,:) = XSVW_FLUX(II_M:II_M+1,IJ_M:IJ_M+1,:,:)

          ! Compute coefficents for horizontal interpolations
          CALL FLYER_COMPUTE_INTERP_COEFF_HOR_STAGE1( )
          ! Compute coefficents for vertical interpolations
          CALL FLYER_COMPUTE_INTERP_COEFF_VER( )
          ! Compute coefficents for horizontal interpolations
          CALL FLYER_COMPUTE_INTERP_COEFF_HOR_STAGE2( )

          CALL FLYER_RECORD_DATA( )
        END IF ISOWNERAIR

        ! Store has been done
        TPFLYER%LSTORE = .FALSE.
      END IF
    END IF

    ! Compute next position if the previous store has just been done (right moment on right model)
    IF ( IMI == TPFLYER%NMODEL .AND. ISTORE > 0 ) THEN
      ! This condition may only be tested if ISTORE > 0
      IF (ABS( TDTCUR - TPFLYER%TFLYER_TIME%TPDATES(ISTORE) ) < 1e-10 ) THEN
        ! Next store moment
        TZNEXT = TDTCUR + TPFLYER%TFLYER_TIME%XTSTEP

        ! Is the aircraft in flight ?
        IF ( TZNEXT >= TPFLYER%TLAUNCH .AND. TZNEXT <= TPFLYER%TLAND ) THEN
          TPFLYER%LFLY = .TRUE.
          ! Force LTOOKOFF to prevent to do it again (at a next timestep)
          TPFLYER%LTOOKOFF = .TRUE.

          ! Compute next position
          CALL AIRCRAFT_COMPUTE_POSITION( TZNEXT, TPFLYER )

          ! Get rank of the process where the aircraft is and the model number
          CALL FLYER_GET_RANK_MODEL_ISCRASHED( TPFLYER )
        ELSE
          TPFLYER%LFLY = .FALSE.
        END IF
      END IF
    END IF

    IF ( GOWNER_CUR ) KRANK_NXT = TPFLYER%NRANK_CUR

  CLASS IS ( TBALLOONDATA)
    GLAUNCH = .FALSE. !Set to true only at the launch instant (set to false in flight after launch)

    ! Launch?
    LAUNCH: IF ( .NOT. TPFLYER%LFLY .AND. .NOT. TPFLYER%LCRASH .AND. TPFLYER%NMODEL == IMI ) THEN
      ! Check if it is launchtime
      LAUNCHTIME: IF ( ( TDTCUR - TPFLYER%TLAUNCH ) >= -1.e-10 ) THEN
        TPFLYER%LFLY = .TRUE.
        GLAUNCH = .TRUE.

        TPFLYER%XX_CUR = TPFLYER%XXLAUNCH
        TPFLYER%XY_CUR = TPFLYER%XYLAUNCH
        TPFLYER%TPOS_CUR = TDTCUR
      END IF LAUNCHTIME
    END IF LAUNCH

    ! Check if it is time to store data. This has also to be checked if the balloon
    ! is not yet launched or is crashed (data is also written in these cases, but with default values)
    IF ( TPFLYER%NMODEL == IMI .AND. &
          ( .NOT. TPFLYER%LFLY .OR. TPFLYER%LCRASH .OR. ABS( TPFLYER%TPOS_CUR - TDTCUR ) < 1.e-8 ) ) THEN
      !Do we have to store balloon data?
      TPFLYER%LSTORE = TPFLYER%TFLYER_TIME%STORESTEP_CHECK_AND_SET( ISTORE )
      IF ( TPFLYER%LSTORE ) TPFLYER%NSTORE_CUR = ISTORE
    END IF

    ! In flight
    INFLIGHTONMODEL: IF ( TPFLYER%LFLY .AND. .NOT. TPFLYER%LCRASH .AND. TPFLYER%NMODEL == IMI &
                          .AND. ABS( TPFLYER%TPOS_CUR - TDTCUR ) < 1.e-8 ) THEN
      ISOWNERBAL: IF ( TPFLYER%NRANK_CUR == ISP ) THEN
        CALL FLYER_INTERP_TO_MASSPOINTS()

        ZEXN(:,:,:) = FLYER_COMPUTE_EXNER( )
        ZRHO(:,:,:) = FLYER_COMPUTE_RHO( )

        ZTHW_FLUX(:,:,:) = ZRHO(:,:,:)*XCPD *XTHW_FLUX(II_M:II_M+1,IJ_M:IJ_M+1,:)
        ZRCW_FLUX(:,:,:) = ZRHO(:,:,:)*XLVTT*XRCW_FLUX(II_M:II_M+1,IJ_M:IJ_M+1,:)
        ZSVW_FLUX(:,:,:,:) = XSVW_FLUX(II_M:II_M+1,IJ_M:IJ_M+1,:,:)

        ! Compute coefficents for horizontal interpolations
        CALL FLYER_COMPUTE_INTERP_COEFF_HOR_STAGE1( )

        IF ( GLAUNCH ) CALL BALLOON_COMPUTE_INITIAL_VERTICAL_POSITION( TPFLYER )

        ! Compute coefficents for vertical interpolations
        CALL FLYER_COMPUTE_INTERP_COEFF_VER( )

        CRASH_VERT: IF ( TPFLYER%LCRASH ) THEN
          TPFLYER%LFLY = .FALSE.
          WRITE( CMNHMSG(1), "( 'Balloon ', A, ' crashed the ', I2, '/', I2, '/', I4, ' at ', F18.12, &
                                's (too low or too high)' )" )                                        &
                 TRIM( TPFLYER%CNAME ), TDTCUR%NDAY, TDTCUR%NMONTH, TDTCUR%NYEAR, TDTCUR%XTIME
          CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'AIRCRAFT_BALLOON_EVOL', OLOCAL = .TRUE. )
        ELSE CRASH_VERT
          !No vertical crash

          ! Compute coefficents for horizontal interpolations
          CALL FLYER_COMPUTE_INTERP_COEFF_HOR_STAGE2( )

          ! Check if it is the right moment to store data
          IF ( TPFLYER%LSTORE ) THEN
            ISTORE = TPFLYER%TFLYER_TIME%N_CUR
            IF ( ABS( TDTCUR - TPFLYER%TFLYER_TIME%TPDATES(ISTORE) ) < 1e-10 ) THEN
              CALL FLYER_RECORD_DATA( )
            END IF
          END IF

          ! Compute next horizontal position (balloon advection)
          CALL BALLOON_ADVECTION_HOR( TPFLYER )

          ! Compute next vertical position (balloon advection)
          CALL BALLOON_ADVECTION_VER( TPFLYER )

          TPFLYER%TPOS_CUR = TDTCUR + ZTSTEP
        END IF CRASH_VERT !end of no vertical crash branch
      END IF ISOWNERBAL
    END IF INFLIGHTONMODEL

    IF ( GOWNER_CUR ) KRANK_NXT = TPFLYER%NRANK_CUR
END SELECT

CONTAINS

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE BALLOON_COMPUTE_INITIAL_VERTICAL_POSITION( TPBALLOON )

USE MODD_CST, ONLY: XCPD, XP00, XRD

IMPLICIT NONE

CLASS(TBALLOONDATA), INTENT(INOUT) :: TPBALLOON

LOGICAL :: GLOW, GHIGH

SELECT CASE ( TPBALLOON%CTYPE )
  !
  ! Iso-density balloon
  !
  CASE ( 'ISODEN' )
    IF ( TPBALLOON%XALTLAUNCH /= XNEGUNDEF ) THEN
      CALL TPBALLOON%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPBALLOON%XALTLAUNCH, ZZM, GLOW, GHIGH )
      TPBALLOON%XRHO = TPBALLOON%INTERP_FROM_MASSPOINT( ZRHO )
    ELSE IF ( TPBALLOON%XPRES /= XNEGUNDEF ) THEN
      ZFLYER_EXN = (TPBALLOON%XPRES/XP00)**(XRD/XCPD)
      CALL TPBALLOON%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', ZFLYER_EXN, ZEXN, GLOW, GHIGH )
      TPBALLOON%XRHO = TPBALLOON%INTERP_FROM_MASSPOINT( ZRHO )
    ELSE
      CMNHMSG(1) = 'Error in balloon initial position (balloon ' // TRIM(TPBALLOON%CNAME) // ' )'
      CMNHMSG(2) = 'neither initial ALTITUDE or PRESsure is given'
      CMNHMSG(3) = 'Check your INI_BALLOON routine'
      CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'AIRCRAFT_BALLOON_EVOL', OLOCAL = .TRUE. )
    END IF
  !
  ! Radiosounding balloon
  !
  CASE ( 'RADIOS' )
    TPBALLOON%XZ_CUR = TPBALLOON%XALTLAUNCH
    TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(1,1,IKB) )
    TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(2,1,IKB) )
    TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(1,2,IKB) )
    TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(2,2,IKB) )
    IF ( TPBALLOON%XZ_CUR > TPBALLOON%XALTLAUNCH ) THEN
      WRITE( CMNHMSG(1), '(A)' ) 'initial vertical position of ' // TRIM( TPBALLOON%CNAME ) // ' was too low'
      WRITE( CMNHMSG(2), '( "forced to ", EN12.3, " (instead of ", EN12.3, ")" )' ) TPBALLOON%XZ_CUR, TPBALLOON%XALTLAUNCH
      CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'BALLOON_COMPUTE_INITIAL_VERTICAL_POSITION', OLOCAL = .TRUE. )
    END IF
  !
  ! Constant Volume Balloon
  !
  CASE ( 'CVBALL' )
    IF ( TPBALLOON%XALTLAUNCH /= XNEGUNDEF ) THEN
      CALL TPBALLOON%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPBALLOON%XALTLAUNCH, ZZM, GLOW, GHIGH )
      IF ( GLOW ) THEN
        TPBALLOON%XZ_CUR = TPBALLOON%XALTLAUNCH
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(1,1,IKB) )
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(2,1,IKB) )
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(1,2,IKB) )
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(2,2,IKB) )

        WRITE( CMNHMSG(1), '(A)' ) 'initial vertical position of ' // TRIM( TPBALLOON%CNAME ) // ' was too low'
        WRITE( CMNHMSG(2), '( "forced to ", EN12.3, " (instead of ", EN12.3, ")" )' ) TPBALLOON%XZ_CUR, TPBALLOON%XALTLAUNCH
        CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'BALLOON_COMPUTE_INITIAL_VERTICAL_POSITION', OLOCAL = .TRUE. )

        !Recompute the vertical interpolation coefficients at the corrected vertical position
        CALL TPBALLOON%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPBALLOON%XALTLAUNCH, ZZM, GLOW, GHIGH )
      ELSE
        TPBALLOON%XZ_CUR = TPBALLOON%INTERP_FROM_MASSPOINT( ZZM )
      END IF
      TPBALLOON%XRHO     = TPBALLOON%INTERP_FROM_MASSPOINT( ZRHO )
    ELSE IF ( TPBALLOON%XPRES /= XNEGUNDEF ) THEN
      ZFLYER_EXN = (TPBALLOON%XPRES/XP00)**(XRD/XCPD)
      CALL TPBALLOON%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', ZFLYER_EXN, ZEXN, GLOW, GHIGH )
      IF ( GLOW ) THEN
        TPBALLOON%XZ_CUR = ZZM(1,1,IKB)
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(2,1,IKB) )
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(1,2,IKB) )
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(2,2,IKB) )

        WRITE( CMNHMSG(1), '(A)' ) 'initial vertical position of ' // TRIM( TPBALLOON%CNAME ) // ' was too low'
        WRITE( CMNHMSG(2), '( "forced to ", EN12.3 )' ) TPBALLOON%XZ_CUR
        CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'BALLOON_COMPUTE_INITIAL_VERTICAL_POSITION', OLOCAL = .TRUE. )

        !Recompute the vertical interpolation coefficients at the corrected vertical position
        CALL TPBALLOON%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPBALLOON%XZ_CUR, ZZM, GLOW, GHIGH )
      ELSE
        TPBALLOON%XZ_CUR = TPBALLOON%INTERP_FROM_MASSPOINT( ZZM )
      END IF
      TPBALLOON%XRHO     = TPBALLOON%INTERP_FROM_MASSPOINT( ZRHO )
    ELSE
      TPBALLOON%XRHO = TPBALLOON%XMASS / TPBALLOON%XVOLUME
      CALL TPBALLOON%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPBALLOON%XRHO, ZRHO, GLOW, GHIGH )
      IF ( GLOW ) THEN
        TPBALLOON%XZ_CUR = ZZM(1,1,IKB)
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(2,1,IKB) )
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(1,2,IKB) )
        TPBALLOON%XZ_CUR = MAX ( TPBALLOON%XZ_CUR , ZZM(2,2,IKB) )

        WRITE( CMNHMSG(1), '(A)' ) 'initial vertical position of ' // TRIM( TPBALLOON%CNAME ) // ' was too low'
        WRITE( CMNHMSG(2), '( "forced to ", EN12.3 )' ) TPBALLOON%XZ_CUR
        CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'BALLOON_COMPUTE_INITIAL_VERTICAL_POSITION', OLOCAL = .TRUE. )

        !Recompute the vertical interpolation coefficients at the corrected vertical position
        CALL TPBALLOON%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPBALLOON%XZ_CUR, ZZM, GLOW, GHIGH )
      ELSE
        TPBALLOON%XZ_CUR = TPBALLOON%INTERP_FROM_MASSPOINT( ZZM )
      END IF
    END IF
END SELECT

END SUBROUTINE BALLOON_COMPUTE_INITIAL_VERTICAL_POSITION
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE BALLOON_ADVECTION_HOR( TPBALLOON )

USE MODD_AIRCRAFT_BALLOON, ONLY: TBALLOONDATA
USE MODD_CONF,             ONLY: LCARTESIAN
USE MODD_NESTING,          ONLY: NDAD, NDTRATIO
USE MODD_TIME,             only: TDTSEG
USE MODD_TIME_n,           ONLY: TDTCUR

IMPLICIT NONE

CLASS(TBALLOONDATA), INTENT(INOUT) :: TPBALLOON

INTEGER :: IMODEL
INTEGER :: IMODEL_OLD
REAL    :: ZX_OLD, ZY_OLD
REAL    :: ZDELTATIME
REAL    :: ZDIVTMP
REAL    :: ZMAP     ! map factor at balloon location
REAL    :: ZU_BAL   ! horizontal wind speed at balloon location (along x)
REAL    :: ZV_BAL   ! horizontal wind speed at balloon location (along y)

ZTSTEP = PTSTEP

ZU_BAL = TPBALLOON%INTERP_FROM_UPOINT( PU )
ZV_BAL = TPBALLOON%INTERP_FROM_VPOINT( PV )
if ( .not. lcartesian ) then
  ZMAP = TPBALLOON%INTERP_HOR_FROM_MASSPOINT( PMAP )
else
  ZMAP = 1.
end if
!
ZX_OLD = TPBALLOON%XX_CUR
ZY_OLD = TPBALLOON%XY_CUR

TPBALLOON%XX_CUR = TPBALLOON%XX_CUR + ZU_BAL * ZTSTEP * ZMAP
TPBALLOON%XY_CUR = TPBALLOON%XY_CUR + ZV_BAL * ZTSTEP * ZMAP

! Compute rank and model for next position
! This is done here because we need to check if there is a change of model (for 'MOB' balloons)
! because position has to be adapted to the timestep of a coarser model (if necessary)
IMODEL_OLD = TPBALLOON%NMODEL

! Get rank of the process where the balloon is and the model number
CALL FLYER_GET_RANK_MODEL_ISCRASHED( TPBALLOON )

IF ( TPBALLOON%LCRASH ) THEN
  WRITE( CMNHMSG(1), "( 'Balloon ', A, ' crashed the ', I2, '/', I2, '/', I4, ' at ', F18.12, &
                      's (out of the horizontal boundaries)' )" ) &
    TRIM( TPBALLOON%CNAME ), TDTCUR%NDAY, TDTCUR%NMONTH, TDTCUR%NYEAR, TDTCUR%XTIME
  CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'AIRCRAFT_BALLOON_EVOL', OLOCAL = .TRUE. )
END IF

IF ( TPBALLOON%NMODEL /= IMODEL_OLD .AND. .NOT. TPBALLOON%LCRASH ) THEN
  ! Balloon has changed of model
  IF ( NDAD(TPBALLOON%NMODEL ) == IMODEL_OLD ) THEN
    ! Nothing special to do when going to child model
  ELSE IF ( TPBALLOON%NMODEL == NDAD(IMODEL_OLD) ) THEN
    ! Balloon go to parent model
    ! Recompute position to be compatible with parent timestep
    ! Parent timestep could be bigger (factor NDTRATIO) and therefore next position is not the one computed just before

    ! Determine step compatible with parent model at next parent timestep
    ZDELTATIME = TDTCUR - TDTSEG
    ZDIVTMP = ZDELTATIME / ( PTSTEP * NDTRATIO(IMODEL_OLD) )
    IF ( ABS( ZDIVTMP - NINT( ZDIVTMP ) ) < 1E-6 * PTSTEP * NDTRATIO(IMODEL_OLD) ) THEN
      ! Current time is a multiple of parent timestep => next position is parent timestep
      ZTSTEP = ZTSTEP * NDTRATIO(IMODEL_OLD)
    ELSE
      ! Current time is not a multiple of parent timestep
      ! Next position must be a multiple of parent timestep
      ! NINT( NDTRATIO(IMODEL_OLD) * ( 1 - ( ZDIVTMP - INT( ZDIVTMP ) ) ) ) corresponds to the number
      ! of child timesteps to go to the next parent timestep
      ! We skip one timestep (+NDTRATIO(IMODEL_OLD)) because it has already been computed for the parent model
      ZTSTEP = ZTSTEP * ( NINT( NDTRATIO(IMODEL_OLD) * ( 1 - ( ZDIVTMP - INT( ZDIVTMP ) ) ) ) + NDTRATIO(IMODEL_OLD) )

      ! Detect if we need to skip a store (if time of next position is after time of next store)
      ! This can happen when a ballon goes to its parent model
      IF ( TDTCUR + ZTSTEP > TPBALLOON%TFLYER_TIME%TPDATES(TPBALLOON%TFLYER_TIME%N_CUR) + TPBALLOON%TFLYER_TIME%XTSTEP + 1e-6 ) THEN
        !Force a dummy store (nothing is computed, therefore default/initial values will be stored)
        TPBALLOON%LSTORE = .TRUE.

        TPBALLOON%TFLYER_TIME%N_CUR = TPBALLOON%TFLYER_TIME%N_CUR + 1
        ISTORE = TPBALLOON%TFLYER_TIME%N_CUR

        !Remark: by construction here, ISTORE is always > 1 => no risk with ISTORE-1 value
        TPBALLOON%TFLYER_TIME%TPDATES(ISTORE) = TPBALLOON%TFLYER_TIME%TPDATES(ISTORE-1) + TPBALLOON%TFLYER_TIME%XTSTEP

        WRITE( CMNHMSG(1), "( 'Balloon ', A, ': store skipped at ', I2, '/', I2, '/', I4, ' at ', F18.12, 's' )" ) &
               TRIM( TPBALLOON%CNAME ),                                                                            &
               TPBALLOON%TFLYER_TIME%TPDATES(ISTORE)%NDAY,  TPBALLOON%TFLYER_TIME%TPDATES(ISTORE)%NMONTH,          &
               TPBALLOON%TFLYER_TIME%TPDATES(ISTORE)%NYEAR, TPBALLOON%TFLYER_TIME%TPDATES(ISTORE)%XTIME
        CMNHMSG(2) = 'due to change of model (child to its parent)'
        CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'AIRCRAFT_BALLOON_EVOL', OLOCAL = .TRUE. )
      END IF
    END IF

    ! Compute new horizontal position
    TPBALLOON%XX_CUR = TPBALLOON%XX_CUR + ZU_BAL * ZTSTEP * ZMAP
    TPBALLOON%XY_CUR = TPBALLOON%XY_CUR + ZV_BAL * ZTSTEP * ZMAP

    ! Get rank of the process where the balloon is and the model number
    ! Model number is now imposed
    IMODEL = TPBALLOON%NMODEL
    CALL FLYER_GET_RANK_MODEL_ISCRASHED( TPBALLOON, KMODEL = IMODEL )
    IF ( TPBALLOON%LCRASH ) THEN
      WRITE( CMNHMSG(1), "( 'Balloon ', A, ' crashed the ', I2, '/', I2, '/', I4, ' at ', F18.12, &
                         's (out of the horizontal boundaries)' )" ) &
        TRIM( TPBALLOON%CNAME ), TDTCUR%NDAY, TDTCUR%NMONTH, TDTCUR%NYEAR, TDTCUR%XTIME
      CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'AIRCRAFT_BALLOON_EVOL', OLOCAL = .TRUE. )
    END IF
  ELSE
    ! Special case not-managed (different dads, change of several models in 1 step (going to grand parent/grand children)...)
    ! This situation should be very infrequent => reasonable risk, error on the trajectory should be relatively small in most cases
    CMNHMSG(1) = 'unmanaged change of model for ballon ' // TPBALLOON%CNAME
    CMNHMSG(2) = 'its trajectory might be wrong'
    CALL PRINT_MSG( NVERB_WARNING, 'GEN', 'AIRCRAFT_BALLOON_EVOL', OLOCAL = .TRUE. )
  END IF
END IF

END SUBROUTINE BALLOON_ADVECTION_HOR
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE BALLOON_ADVECTION_VER( TPBALLOON )

USE MODD_AIRCRAFT_BALLOON, ONLY: TBALLOONDATA
USE MODD_CST,              ONLY: XG

IMPLICIT NONE

CLASS(TBALLOONDATA), INTENT(INOUT) :: TPBALLOON

INTEGER :: JK ! loop index
REAL    :: ZRO_BAL  ! air density at balloon location
REAL    :: ZW_BAL   ! vertical   wind speed at balloon location (along z)

IF ( TPBALLOON%CTYPE == 'RADIOS' ) THEN
  ZW_BAL = TPBALLOON%INTERP_FROM_MASSPOINT( ZWM )
  TPBALLOON%XZ_CUR = TPBALLOON%XZ_CUR + ( ZW_BAL + TPBALLOON%XWASCENT ) * ZTSTEP
END IF

IF ( TPBALLOON%CTYPE == 'CVBALL' ) THEN
  ZW_BAL  = TPBALLOON%INTERP_FROM_MASSPOINT( ZWM )
  ZRO_BAL = TPBALLOON%INTERP_FROM_MASSPOINT( ZRHO )
  ! calculation with a time step of 1 second or less
  IF (INT(ZTSTEP) .GT. 1 ) THEN
    DO JK=1,INT(ZTSTEP)
      TPBALLOON%XWASCENT = TPBALLOON%XWASCENT &
        -  ( 1. / (1. + TPBALLOON%XINDDRAG ) ) * 1. * &
           ( XG * ( ( TPBALLOON%XMASS / TPBALLOON%XVOLUME ) - ZRO_BAL ) / ( TPBALLOON%XMASS / TPBALLOON%XVOLUME ) &
              + TPBALLOON%XWASCENT * ABS ( TPBALLOON%XWASCENT ) * &
                TPBALLOON%XDIAMETER * TPBALLOON%XAERODRAG / ( 2. * TPBALLOON%XVOLUME ) &
            )
      TPBALLOON%XZ_CUR = TPBALLOON%XZ_CUR + ( ZW_BAL + TPBALLOON%XWASCENT ) * 1.
    END DO
  END IF
  IF (ZTSTEP .GT. INT(ZTSTEP)) THEN
    TPBALLOON%XWASCENT = TPBALLOON%XWASCENT &
      -  ( 1. / (1. + TPBALLOON%XINDDRAG ) ) * (ZTSTEP-INT(ZTSTEP)) * &
         ( XG * ( ( TPBALLOON%XMASS / TPBALLOON%XVOLUME ) - ZRO_BAL ) / ( TPBALLOON%XMASS / TPBALLOON%XVOLUME ) &
            + TPBALLOON%XWASCENT * ABS ( TPBALLOON%XWASCENT ) * &
              TPBALLOON%XDIAMETER * TPBALLOON%XAERODRAG / ( 2. * TPBALLOON%XVOLUME ) &
          )
    TPBALLOON%XZ_CUR = TPBALLOON%XZ_CUR + ( ZW_BAL + TPBALLOON%XWASCENT ) * (ZTSTEP-INT(ZTSTEP))
  END IF
END IF

END SUBROUTINE BALLOON_ADVECTION_VER
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE FLYER_INTERP_TO_MASSPOINTS()

USE MODD_GRID_n,     ONLY: XXHAT, XXHATM, XYHAT, XYHATM
USE MODD_PARAMETERS, ONLY: JPVEXT

IMPLICIT NONE

INTEGER :: IDU      ! difference between II_U and II_M
INTEGER :: IDV      ! difference between IJ_V and IJ_M

! Indices
IKB =   1   + JPVEXT
IKE = SIZE(PZ,3) - JPVEXT

! Interpolations of model variables to mass points
! ------------------------------------------------

! X position
TPFLYER%NI_U = COUNT( XXHAT (:) <= TPFLYER%XX_CUR )
TPFLYER%NI_M = COUNT( XXHATM(:) <= TPFLYER%XX_CUR )
II_U = TPFLYER%NI_U
II_M = TPFLYER%NI_M

! Y position
TPFLYER%NJ_V = COUNT( XYHAT (:)<=TPFLYER%XY_CUR )
TPFLYER%NJ_M = COUNT( XYHATM(:)<=TPFLYER%XY_CUR )
IJ_V = TPFLYER%NJ_V
IJ_M = TPFLYER%NJ_M

ZZM(:,:,1:IKU-1)=0.5 *PZ(II_M  :II_M+1,IJ_M  :IJ_M+1,1:IKU-1)+0.5 *PZ(II_M  :II_M+1,IJ_M  :IJ_M+1,2:IKU  )
ZZM(:,:,  IKU  )=1.5 *PZ(II_M  :II_M+1,IJ_M  :IJ_M+1,  IKU-1)-0.5 *PZ(II_M  :II_M+1,IJ_M  :IJ_M+1,  IKU-2)

IDU = II_U - II_M
ZZU(:,:,1:IKU-1)=0.25*PZ(IDU+II_M-1:IDU+II_M,  IJ_M  :IJ_M+1,1:IKU-1)+0.25*PZ(IDU+II_M-1:IDU+II_M  ,IJ_M  :IJ_M+1,2:IKU  ) &
                +0.25*PZ(IDU+II_M  :IDU+II_M+1,IJ_M  :IJ_M+1,1:IKU-1)+0.25*PZ(IDU+II_M  :IDU+II_M+1,IJ_M  :IJ_M+1,2:IKU  )
ZZU(:,:,  IKU  )=0.75*PZ(IDU+II_M-1:IDU+II_M  ,IJ_M  :IJ_M+1,  IKU-1)-0.25*PZ(IDU+II_M-1:IDU+II_M  ,IJ_M  :IJ_M+1,  IKU-2) &
                +0.75*PZ(IDU+II_M  :IDU+II_M+1,IJ_M  :IJ_M+1,  IKU-1)-0.25*PZ(IDU+II_M  :IDU+II_M+1,IJ_M  :IJ_M+1,  IKU-2)

IDV = IJ_V - IJ_M
ZZV(:,:,1:IKU-1)=0.25*PZ(II_M  :II_M+1,IDV+IJ_M-1:IDV+IJ_M  ,1:IKU-1)+0.25*PZ(II_M  :II_M+1,IDV+IJ_M-1:IDV+IJ_M  ,2:IKU  ) &
                +0.25*PZ(II_M  :II_M+1,IDV+IJ_M  :IDV+IJ_M+1,1:IKU-1)+0.25*PZ(II_M  :II_M+1,IDV+IJ_M  :IDV+IJ_M+1,2:IKU  )
ZZV(:,:,  IKU  )=0.75*PZ(II_M  :II_M+1,IDV+IJ_M-1:IDV+IJ_M  ,  IKU-1)-0.25*PZ(II_M  :II_M+1,IDV+IJ_M-1:IDV+IJ_M  ,  IKU-2) &
                +0.75*PZ(II_M  :II_M+1,IDV+IJ_M  :IDV+IJ_M+1,  IKU-1)-0.25*PZ(II_M  :II_M+1,IDV+IJ_M  :IDV+IJ_M+1,  IKU-2)

ZWM(:,:,1:IKU-1)=0.5*PW(II_M:II_M+1,IJ_M:IJ_M+1,1:IKU-1)+0.5*PW(II_M:II_M+1,IJ_M:IJ_M+1,2:IKU  )
ZWM(:,:,  IKU  )=1.5*PW(II_M:II_M+1,IJ_M:IJ_M+1,  IKU-1)-0.5*PW(II_M:II_M+1,IJ_M:IJ_M+1,  IKU-2)

END SUBROUTINE FLYER_INTERP_TO_MASSPOINTS
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
PURE FUNCTION FLYER_COMPUTE_EXNER( ) RESULT( PEXN )

USE MODD_CST, ONLY: XCPD, XP00, XRD

IMPLICIT NONE

REAL, DIMENSION(2,2,SIZE(PTH,3)) :: PEXN

INTEGER :: JK

PEXN(:,:,:) = ( PP(II_M:II_M+1, IJ_M:IJ_M+1, :) / XP00) ** ( XRD / XCPD )
DO JK = IKB-1, 1, -1
  PEXN(:,:,JK) = 1.5 * PEXN(:,:,JK+1) - 0.5 * PEXN(:,:,JK+2)
END DO
DO JK = IKE+1, IKU
  PEXN(:,:,JK) = 1.5 * PEXN(:,:,JK-1) - 0.5 * PEXN(:,:,JK-2)
END DO

END FUNCTION FLYER_COMPUTE_EXNER
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
PURE FUNCTION FLYER_COMPUTE_RHO( ) RESULT( PRHO )

USE MODD_CST, ONLY: XRD, XRV

USE MODI_WATER_SUM

IMPLICIT NONE

REAL, DIMENSION(2,2,SIZE(PTH,3)) :: PRHO

INTEGER                          :: JK
REAL, DIMENSION(2,2,SIZE(PTH,3)) :: ZTHV ! virtual potential temperature

ZTHV(:,:,:) = PTH(II_M:II_M+1, IJ_M:IJ_M+1, :)
IF ( SIZE( PR, 4 ) > 0 )                                                              &
  ZTHV(:,:,:) = ZTHV(:,:,:) * ( 1. + XRV / XRD * PR(II_M:II_M+1, IJ_M:IJ_M+1, :, 1) ) &
                            / ( 1. + WATER_SUM( PR(II_M:II_M+1, IJ_M:IJ_M+1, :, :)) )
!
PRHO(:,:,:) = PP(II_M:II_M+1, IJ_M:IJ_M+1, :) / ( XRD * ZTHV(:,:,:) * ZEXN(:,:,:) )
DO JK = IKB-1, 1, -1
  PRHO(:,:,JK) = 1.5 * PRHO(:,:,JK+1) - 0.5 * PRHO(:,:,JK+2)
END DO
DO JK = IKE+1, IKU
  PRHO(:,:,JK) = 1.5 * PRHO(:,:,JK-1) - 0.5 * PRHO(:,:,JK-2)
END DO

END FUNCTION FLYER_COMPUTE_RHO
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE FLYER_COMPUTE_INTERP_COEFF_HOR_STAGE1( )
! Compute coefficents for horizontal interpolations (1st stage)

USE MODD_GRID_n, ONLY: XXHAT, XXHATM, XYHAT, XYHATM

IMPLICIT NONE

! Interpolation coefficient for X
TPFLYER%XXMCOEF = ( TPFLYER%XX_CUR - XXHATM(II_M) ) / ( XXHATM(II_M+1) - XXHATM(II_M) )
TPFLYER%XXMCOEF = MAX( 0., MIN( TPFLYER%XXMCOEF, 1. ) )

! Interpolation coefficient for y
TPFLYER%XYMCOEF = ( TPFLYER%XY_CUR - XYHATM(IJ_M) ) / ( XYHATM(IJ_M+1) - XYHATM(IJ_M) )
TPFLYER%XYMCOEF = MAX( 0., MIN( TPFLYER%XYMCOEF, 1. ) )

! Interpolation coefficient for X (for U)
TPFLYER%XXUCOEF = ( TPFLYER%XX_CUR - XXHAT(II_U) ) / ( XXHAT(II_U+1) - XXHAT(II_U) )
TPFLYER%XXUCOEF = MAX( 0., MIN( TPFLYER%XXUCOEF, 1. ) )

! Interpolation coefficient for y (for V)
TPFLYER%XYVCOEF = ( TPFLYER%XY_CUR - XYHAT(IJ_V) ) / ( XYHAT(IJ_V+1) - XYHAT(IJ_V) )
TPFLYER%XYVCOEF = MAX( 0., MIN( TPFLYER%XYVCOEF, 1. ) )

END SUBROUTINE FLYER_COMPUTE_INTERP_COEFF_HOR_STAGE1
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE FLYER_COMPUTE_INTERP_COEFF_VER( )
! Compute coefficent for vertical interpolations

USE MODD_CST,    ONLY: XCPD, XP00, XRD
USE MODD_TIME_n, ONLY: TDTCUR

IMPLICIT NONE

LOGICAL :: GLOW, GHIGH

! Find indices surrounding the vertical box where the flyer is
SELECT TYPE ( TPFLYER )
  CLASS IS ( TAIRCRAFTDATA)
    IF ( TPFLYER%LALTDEF ) THEN
      ZFLYER_EXN = (TPFLYER%XP_CUR/XP00)**(XRD/XCPD)
      CALL TPFLYER%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', ZFLYER_EXN,     ZEXN, GLOW, GHIGH, ODONOLOWCRASH = .TRUE. )
    ELSE
      CALL TPFLYER%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPFLYER%XZ_CUR, ZZM,  GLOW, GHIGH, ODONOLOWCRASH = .TRUE. )
    END IF

  CLASS IS ( TBALLOONDATA)
    IF ( TPFLYER%CTYPE == 'ISODEN' ) THEN
      CALL TPFLYER%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPFLYER%XRHO,   ZRHO, GLOW, GHIGH, ODONOLOWCRASH = .TRUE. )
    ELSE IF ( TPFLYER%CTYPE == 'RADIOS' .OR. TPFLYER%CTYPE == 'CVBALL' ) THEN
      CALL TPFLYER%COMPUTE_VERTICAL_INTERP_COEFF( 'MASS', TPFLYER%XZ_CUR, ZZM,  GLOW, GHIGH, ODONOLOWCRASH = .TRUE. )
    END IF

END SELECT

! Check if the flyer crashed vertically (higher bound)
IF ( GHIGH ) THEN
  TPFLYER%LCRASH = .TRUE.
  TPFLYER%NCRASH = NCRASH_OUT_HIGH
END IF

SELECT TYPE ( TPFLYER )
  CLASS IS ( TAIRCRAFTDATA)
    IF ( TPFLYER%LALTDEF ) THEN
      TPFLYER%XZ_CUR = TPFLYER%INTERP_FROM_MASSPOINT( ZZM )
    ELSE
      TPFLYER%XP_CUR = TPFLYER%INTERP_FROM_MASSPOINT( PP )
    END IF

  CLASS IS ( TBALLOONDATA)
    IF ( TPFLYER%CTYPE == 'ISODEN' ) THEN
      TPFLYER%XZ_CUR = TPFLYER%INTERP_FROM_MASSPOINT( ZZM )
    ELSE IF ( TPFLYER%CTYPE == 'RADIOS' .OR. TPFLYER%CTYPE == 'CVBALL' ) THEN
      !Nothing to do
    END IF

END SELECT

END SUBROUTINE FLYER_COMPUTE_INTERP_COEFF_VER
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE FLYER_COMPUTE_INTERP_COEFF_HOR_STAGE2( )
! Compute coefficents for horizontal interpolations (2nd stage)
! This stage must be done after FLYER_COMPUTE_INTERP_COEFF_VER because we should need XZ_CUR computed in it

IMPLICIT NONE

LOGICAL :: GLOW, GHIGH

! Interpolation coefficients for the 4 surroundings verticals (for U)
! ODONOLOWCRASH = .TRUE. because check for low crash has already been done
CALL TPFLYER%COMPUTE_VERTICAL_INTERP_COEFF( 'U', TPFLYER%XZ_CUR, ZZU, GLOW, GHIGH, ODONOLOWCRASH = .TRUE. )

! Interpolation coefficients for the 4 suroundings verticals (for V)
CALL TPFLYER%COMPUTE_VERTICAL_INTERP_COEFF( 'V', TPFLYER%XZ_CUR, ZZV, GLOW, GHIGH, ODONOLOWCRASH = .TRUE. )

END SUBROUTINE FLYER_COMPUTE_INTERP_COEFF_HOR_STAGE2
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE FLYER_RECORD_DATA( )

USE MODD_CST,              ONLY: XP00, XPI, XRD
USE MODD_DIAG_IN_RUN,      ONLY: XCURRENT_TKE_DISS
USE MODD_GRID,             ONLY: XBETA, XLON0, XRPK
USE MODD_NSV,              ONLY: NSV_LIMA_NC, NSV_LIMA_NR, NSV_LIMA_NI
USE MODD_PARAMETERS,       ONLY: JPVEXT
USE MODD_PARAM_n,          ONLY: CCLOUD, CRAD

USE MODE_GRIDPROJ,         ONLY: SM_LATLON
USE MODE_SENSOR,           ONLY: Sensor_rare_compute, Sensor_wc_compute

IMPLICIT NONE

INTEGER                        :: JLOOP    ! loop counter
REAL                           :: ZGAM     ! rotation between meso-nh base and spherical lat-lon base.
REAL                           :: ZU_BAL   ! horizontal wind speed at balloon location (along x)
REAL                           :: ZV_BAL   ! horizontal wind speed at balloon location (along y)
REAL, DIMENSION(SIZE(PZ,3))    :: ZZ       ! altitude of model levels at station location
REAL, DIMENSION(SIZE(PR,1),SIZE(PR,2),SIZE(PR,3))    :: ZR

TPFLYER%NMODELHIST(ISTORE) = TPFLYER%NMODEL

TPFLYER%XX(ISTORE) = TPFLYER%XX_CUR
TPFLYER%XY(ISTORE) = TPFLYER%XY_CUR
TPFLYER%XZ(ISTORE) = TPFLYER%XZ_CUR
!
CALL SM_LATLON( PLATOR, PLONOR,                    &
                TPFLYER%XX_CUR,   TPFLYER%XY_CUR,  &
                TPFLYER%XLAT_CUR, TPFLYER%XLON_CUR )
TPFLYER%XLAT(ISTORE) = TPFLYER%XLAT_CUR
TPFLYER%XLON(ISTORE) = TPFLYER%XLON_CUR
!
ZU_BAL = TPFLYER%INTERP_FROM_UPOINT( PU )
ZV_BAL = TPFLYER%INTERP_FROM_VPOINT( PV )
ZGAM   = (XRPK * (TPFLYER%XLON_CUR - XLON0) - XBETA)*(XPI/180.)
TPFLYER%XZON (1,ISTORE) = ZU_BAL * COS(ZGAM) + ZV_BAL * SIN(ZGAM)
TPFLYER%XMER (1,ISTORE) = - ZU_BAL * SIN(ZGAM) + ZV_BAL * COS(ZGAM)
!
TPFLYER%XW   (1,ISTORE) = TPFLYER%INTERP_FROM_MASSPOINT( ZWM )
TPFLYER%XTH  (1,ISTORE) = TPFLYER%INTERP_FROM_MASSPOINT( PTH )
!
ZFLYER_EXN = TPFLYER%INTERP_FROM_MASSPOINT( ZEXN )
TPFLYER%XP   (1,ISTORE) = XP00 * ZFLYER_EXN**(XCPD/XRD)

ZR(:,:,:) = 0.
DO JLOOP=1,SIZE(PR,4)
  TPFLYER%XR   (1,ISTORE,JLOOP) = TPFLYER%INTERP_FROM_MASSPOINT( PR(:,:,:,JLOOP) )
  IF (JLOOP>=2) ZR(:,:,:) = ZR(:,:,:) + PR(:,:,:,JLOOP)
END DO
DO JLOOP=1,SIZE(PSV,4)
  TPFLYER%XSV  (1,ISTORE,JLOOP) = TPFLYER%INTERP_FROM_MASSPOINT( PSV(:,:,:,JLOOP) )
END DO
TPFLYER%XRTZ  (:,ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( ZR(:,:,:) )
DO JLOOP=1,SIZE(PR,4)
  TPFLYER%XRZ  (:,ISTORE,JLOOP) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( PR(:,:,:,JLOOP) )
END DO

TPFLYER%XFFZ  (:,ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( SQRT(PU**2+PV**2) )

TPFLYER%XRHOD (:,ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( PRHODREF )

IF (CCLOUD=="LIMA") THEN
  TPFLYER%XCIZ  (:,ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( PSV(:,:,:,NSV_LIMA_NI) )
  TPFLYER%XCCZ  (:,ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( PSV(:,:,:,NSV_LIMA_NC) )
  TPFLYER%XCRZ  (:,ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( PSV(:,:,:,NSV_LIMA_NR) )
ELSE IF ( CCLOUD=="ICE3" .OR. CCLOUD=="ICE4" ) THEN
  TPFLYER%XCIZ  (:,ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( PCIT(:,:,:) )
END IF

ZTH_EXN(:,:,:) = PTH(TPFLYER%NI_M:TPFLYER%NI_M+1, TPFLYER%NJ_M:TPFLYER%NJ_M+1, :) * ZEXN(:,:,:)
ZZ(:) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( ZZM(:,:,:) )
TPFLYER%XZZ(:,ISTORE) = ZZ(:)

CALL Sensor_wc_compute(   TPFLYER, ISTORE, PR, PRHODREF )
CALL Sensor_rare_compute( TPFLYER, ISTORE, PR, PSV, PRHODREF, PCIT, ZTH_EXN, ZZ, PSEA )

! vertical wind
TPFLYER%XWZ  (:,ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT( ZWM(:,:,:) )

! Dry air density at flyer position
TPFLYER%XRHOD_SENSOR(ISTORE) = TPFLYER%INTERP_FROM_MASSPOINT( PRHODREF )

IF (SIZE(PTKE)>0) TPFLYER%XTKE  (1,ISTORE) = TPFLYER%INTERP_FROM_MASSPOINT( PTKE )
IF ( CRAD /= 'NONE' ) TPFLYER%XTSRAD(ISTORE) = TPFLYER%INTERP_HOR_FROM_MASSPOINT(PTS )
TPFLYER%XTKE_DISS(ISTORE) = TPFLYER%INTERP_FROM_MASSPOINT( XCURRENT_TKE_DISS )
TPFLYER%XZS(ISTORE)       = TPFLYER%INTERP_HOR_FROM_MASSPOINT( PZ(:,:,1+JPVEXT) )
TPFLYER%XTHW_FLUX(ISTORE) = TPFLYER%INTERP_FROM_MASSPOINT( ZTHW_FLUX )
TPFLYER%XRCW_FLUX(ISTORE) = TPFLYER%INTERP_FROM_MASSPOINT( ZRCW_FLUX )
DO JLOOP=1,SIZE(PSV,4)
TPFLYER%XSVW_FLUX(ISTORE,JLOOP) = TPFLYER%INTERP_FROM_MASSPOINT( ZSVW_FLUX(:,:,:,JLOOP) )
END DO

END SUBROUTINE FLYER_RECORD_DATA
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
END SUBROUTINE AIRCRAFT_BALLOON_EVOL
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE AIRCRAFT_COMPUTE_POSITION( TPDATE, TPAIRCRAFT )

USE MODD_AIRCRAFT_BALLOON, ONLY: TAIRCRAFTDATA
USE MODD_TYPE_DATE,        ONLY: DATE_TIME

USE MODE_DATETIME
USE MODE_POSITION_TOOLS,   ONLY: FIND_PROCESS_AND_MODEL_FROM_XY_POS

IMPLICIT NONE

TYPE(DATE_TIME),      INTENT(IN)    :: TPDATE
CLASS(TAIRCRAFTDATA), INTENT(INOUT) :: TPAIRCRAFT !aircraft

INTEGER :: IL        ! flight segment index
REAL    :: ZTDIST    ! time since launch (sec)
REAL    :: ZSEG_FRAC ! fraction of flight in the current segment

! Find the flight segment
ZTDIST = TPDATE - TPAIRCRAFT%TLAUNCH
IL = TPAIRCRAFT%NPOSCUR
DO WHILE ( ZTDIST > TPAIRCRAFT%XPOSTIME(IL+1) )
  IL = IL + 1
  IF ( IL > TPAIRCRAFT%NPOS-1 ) THEN
    !Security (should not happen)
    IL = TPAIRCRAFT%NPOS-1
    EXIT
  END IF
END DO
TPAIRCRAFT%NPOSCUR = IL

! Compute the current position
ZSEG_FRAC = ( ZTDIST - TPAIRCRAFT%XPOSTIME(IL) ) / ( TPAIRCRAFT%XPOSTIME(IL+1) - TPAIRCRAFT%XPOSTIME(IL) )

TPAIRCRAFT%XX_CUR = (1.-ZSEG_FRAC) * TPAIRCRAFT%XPOSX(IL  ) &
               +     ZSEG_FRAC  * TPAIRCRAFT%XPOSX(IL+1)
TPAIRCRAFT%XY_CUR = (1.-ZSEG_FRAC) * TPAIRCRAFT%XPOSY(IL  ) &
               +     ZSEG_FRAC  * TPAIRCRAFT%XPOSY(IL+1)

IF (TPAIRCRAFT%LALTDEF) THEN
  TPAIRCRAFT%XP_CUR = (1.-ZSEG_FRAC) * TPAIRCRAFT%XPOSP(IL  ) &
                 +     ZSEG_FRAC  * TPAIRCRAFT%XPOSP(IL+1)
ELSE
  TPAIRCRAFT%XZ_CUR = (1.-ZSEG_FRAC) * TPAIRCRAFT%XPOSZ(IL ) &
                 +     ZSEG_FRAC  * TPAIRCRAFT%XPOSZ(IL +1)
END IF

END SUBROUTINE AIRCRAFT_COMPUTE_POSITION
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
SUBROUTINE FLYER_GET_RANK_MODEL_ISCRASHED( TPFLYER, PX, PY, KMODEL )

USE MODD_AIRCRAFT_BALLOON, ONLY: NCRASH_NO, NCRASH_OUT_HORIZ, TFLYERDATA

USE MODE_POSITION_TOOLS,   ONLY: FIND_PROCESS_AND_MODEL_FROM_XY_POS

IMPLICIT NONE

CLASS(TFLYERDATA), INTENT(INOUT) :: TPFLYER ! balloon/aircraft
REAL,    OPTIONAL, INTENT(IN)    :: PX      ! X position (if not provided, takes current flyer position)
REAL,    OPTIONAL, INTENT(IN)    :: PY      ! Y position (if not provided, takes current flyer position)
INTEGER, OPTIONAL, INTENT(IN)    :: KMODEL  ! if provided, model number is imposed (if not 0)

INTEGER :: IMODEL
INTEGER :: IRANK
REAL    :: ZX, ZY

IF ( PRESENT( KMODEL ) ) THEN
  IMODEL = KMODEL
ELSE
  IF ( TPFLYER%CMODEL == 'FIX' ) THEN
    IMODEL = TPFLYER%NMODEL
  ELSE
    IMODEL = 0
  END IF
END IF

IF ( PRESENT( PX ) ) THEN
  ZX = PX
ELSE
  ZX = TPFLYER%XX_CUR
END IF

IF ( PRESENT( PY ) ) THEN
  ZY = PY
ELSE
  ZY = TPFLYER%XY_CUR
END IF

CALL FIND_PROCESS_AND_MODEL_FROM_XY_POS( ZX, ZY, IRANK, IMODEL )

IF ( IRANK < 1 ) THEN
  ! Flyer is outside of horizontal domain
  ! TPFLYER%NMODEL !Do not change to keep a valid value
  TPFLYER%LCRASH    = .TRUE.
  TPFLYER%NCRASH    = NCRASH_OUT_HORIZ
  TPFLYER%LFLY      = .FALSE.
ELSE
  TPFLYER%NMODEL    = IMODEL
  TPFLYER%LCRASH    = .FALSE.
  TPFLYER%NCRASH    = NCRASH_NO
  !TPFLYER%LFLY      = !Do not touch LFLY (flyer could be in flight or not)
  TPFLYER%NRANK_CUR = IRANK
END IF

END SUBROUTINE FLYER_GET_RANK_MODEL_ISCRASHED
!----------------------------------------------------------------------------

END MODULE MODE_AIRCRAFT_BALLOON_EVOL
