!MNH_LIC Copyright 2010-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_SERIES_CLOUD_ELEC
!     #############################
!
INTERFACE
  SUBROUTINE SERIES_CLOUD_ELEC (KTCOUNT, PTSTEP,                &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, &
                                PRT, PRS, PSVT,                 &
                                PTHT, PWT, PPABST, PCIT,        &
                                TPFILE_SERIES_CLOUD_ELEC,       &
                                PINPRR                          )
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
!
REAL,                     INTENT(IN)    :: PTSTEP   ! Double time step except for cold start
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! Moist variables at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRS     ! Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT    ! Scalar variable at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT   ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PWT    ! Vertical velocity at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST ! ab. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT   ! Pristine ice number
                                                  ! concentration at time t
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE_SERIES_CLOUD_ELEC
REAL, DIMENSION(:,:),     INTENT(IN)    :: PINPRR  ! Rain instant precip
!
END SUBROUTINE SERIES_CLOUD_ELEC
END INTERFACE
END MODULE MODI_SERIES_CLOUD_ELEC
!
!      
! ###############################################################
  SUBROUTINE SERIES_CLOUD_ELEC (KTCOUNT, PTSTEP,                &
                                PZZ, PRHODJ, PRHODREF, PEXNREF, &
                                PRT, PRS, PSVT,                 &
                                PTHT, PWT, PPABST, PCIT,        &
                                TPFILE_SERIES_CLOUD_ELEC,       &
                                PINPRR                          )
! ###############################################################
!
!!****  * -
!!
!!    PURPOSE
!!    -------
!!
!!    METHOD
!!    ------
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
!!    AUTHOR
!!    ------
!!      C. Bovalo   * LA * 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original : Avril 2010
!!      Modifications:
!!      C. Barthe  * LACy *  Dec. 2010    add some parameters
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      Philippe Wautelet: 10/01/2019: use NEWUNIT argument of OPEN
!!      Philippe Wautelet: 22/01/2019: use standard FLUSH statement instead of non standard intrinsics
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  C. Barthe   20/03/2023: PRINPRR passed as input argument only
!
!-------------------------------------------------------------------------------
!
! 0. DECLARATIONS
! ------------
!
USE MODD_CONF,            ONLY: CEXP
USE MODD_CST
USE MODD_DYN_n,           ONLY: XDXHATM, XDYHATM
USE MODD_ELEC_DESCR
USE MODD_ELEC_PARAM
USE MODD_IO,              ONLY: TFILEDATA
USE MODD_NSV,             ONLY: NSV_ELECBEG, NSV_ELECEND
USE MODD_PARAMETERS
USE MODD_RAIN_ICE_DESCR_n
USE MODD_RAIN_ICE_PARAM_n
USE MODD_REF

USE MODI_MOMG
USE MODI_RADAR_RAIN_ICE

USE MODE_ELEC_ll
USE MODE_ll
use mode_tools,           only: Countjv

IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
!
REAL,                     INTENT(IN)    :: PTSTEP   ! Double time step except for cold start
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ     ! Height (z)
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PEXNREF ! Reference Exner function
!
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT     ! Moist variables at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRS     ! Moist  variable sources
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSVT    ! Scalar variable at time t
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT   ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PWT    ! Vertical velocity at t-dt
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST ! ab. pressure at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT   ! Pristine ice number
                                                  ! concentration at time t
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE_SERIES_CLOUD_ELEC
REAL, DIMENSION(:,:),     INTENT(IN)    :: PINPRR  ! Rain instant precip
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: II, IJ, IK
INTEGER :: IIB,IIE    ! Indices for the first and last inner mass point along x
INTEGER :: IJB,IJE    ! Indices for the first and last inner mass point along y
INTEGER :: IKB,IKE    ! Indices for the first and last inner mass point along z
INTEGER :: JCOUNT_STOP
INTEGER :: ICOUNT     ! counter for iwp computation
INTEGER :: IPROC      ! my proc number
INTEGER :: IPROC_MAX  ! proc that contains max value
INTEGER :: IINFO_ll   ! return code of parallel routine
INTEGER :: ILU        ! unit number for IO
!
INTEGER, SAVE :: JCOUNT
!
INTEGER, DIMENSION(SIZE(PRT,1),SIZE(PRT,2)) :: IFLAG
!
REAL :: ZRHO00     ! Surface reference air density
REAL :: ZMASS_SP   ! Precipitation snow mass (kg)
REAL :: ZMASS_GP   ! Precipitation graupel mass (kg)
REAL :: ZFLUX_I    ! Ice crystal mass flux (kg m/s)
REAL :: ZFLUX_SP   ! Precipitation snow mass flux (kg m/s)
REAL :: ZFLUX_SNP  ! Non precipitation snow mass flux (kg m/s)
REAL :: ZFLUX_G    ! Graupel mass flux (kg m/s)
REAL :: ZCLD_TOP_REF ! Cloud top height (m) from radar refl.
REAL :: ZCLD_TOP_MR  ! Cloud top height (m) from mixing ratio
REAL :: ZICE_MASS  ! Ice mass (kg)
!
REAL, SAVE :: ZMASS_C      ! Cloud water mass (kg)
REAL, SAVE :: ZMASS_R      ! Rain water mass (kg)
REAL, SAVE :: ZMASS_I      ! Ice crystal mass (kg)
REAL, SAVE :: ZMASS_S      ! Snow mass (kg)
REAL, SAVE :: ZMASS_G      ! Graupel mass (kg)
REAL, SAVE :: ZMASS_ICE_P  ! Precipitation ice mass (kg)
REAL, SAVE :: ZFLUX_PROD   ! Ice mass flux product (kg^2 m^2/s^2)
REAL, SAVE :: ZFLUX_PRECIP ! Precipitation ice mass flux (kg m/s)
REAL, SAVE :: ZFLUX_NPRECIP ! Non-precipitation ice mass flux (kg m/s)
REAL, SAVE :: ZVOL_UP5     ! Updraft volume for W > 5 m/s (m^3)
REAL, SAVE :: ZVOL_UP10    ! Updraft volume for W > 10 m/s (m^3)
REAL, SAVE :: ZWMAX        ! Maximum vertical velocity (m/s)
REAL, SAVE :: ZVOL_G       ! Graupel volume (m^3)
REAL, SAVE :: ZIWP         ! Ice water path (kg/m^2)
REAL, SAVE :: ZCTH_MR      ! Cloud top height / m.r. > 1.e-4 kg/kg (m)
REAL, SAVE :: ZCTH_REF     ! Cloud top height / Z > 20 dBZ (m)
REAL, SAVE :: ZCLD_VOL     ! Cloud volume (m^3)
REAL, SAVE :: ZDBZMAX      ! Max radar reflectivity (dBZ)
REAL, SAVE :: ZINPRR       ! Rain instant precip. (mm/H)
REAL, SAVE :: ZMAX_INPRR   ! Maximum rain instant. precip. (mm/H)
!
REAL, DIMENSION(SIZE(XRTMIN)) :: ZRTMIN
! XRTMIN = Minimum value for the mixing ratio
! ZRTMIN = Minimum value for the source (tendency)
!
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: &
         ZTCT  ! Temperature in Degrees Celsius
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: &
         ZWORK31, ZWORK32, ZWORK33, ZWORK34
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCLOUD
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLAMBDAS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLAMBDAG
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZVTS
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZVTG
!
LOGICAL, SAVE :: GFIRSTCALL = .TRUE.
!
!
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS AND SOME PARAMETERS
!   	        -------------------------------------------
!
JCOUNT_STOP = INT(NTSAVE_SERIES/PTSTEP)
!
!*       1.1     beginning and end indexes of the physical subdomain
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PZZ,3) - JPVEXT
!
!
!*       1.2     compute some parameters
!
! temperature : K -> C
ZTCT(:,:,:) = (PTHT(:,:,:) * (PPABST(:,:,:) / XP00)**(XRD/XCPD)) - XTT
!
! total mixing ratio
ALLOCATE(ZCLOUD(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
ZCLOUD(:,:,:) = 0.
ZCLOUD(IIB:IIE,IJB:IJE,IKB:IKE) = PRT(IIB:IIE,IJB:IJE,IKB:IKE,2) + &
                                  PRT(IIB:IIE,IJB:IJE,IKB:IKE,3) + &
                                  PRT(IIB:IIE,IJB:IJE,IKB:IKE,4) + &
                                  PRT(IIB:IIE,IJB:IJE,IKB:IKE,5) + &
                                  PRT(IIB:IIE,IJB:IJE,IKB:IKE,6)
!
!
!*       1.3     compute the terminal fall speed
!
! the mean terminal fall speed is computed following:
! V_mean = Int(v(D) n(D) dD) / Int(n(D) dD)
!
ALLOCATE(ZLAMBDAS(SIZE(PRT,1), SIZE(PRT,2), SIZE(PRT,3)))
ALLOCATE(ZLAMBDAG(SIZE(PRT,1), SIZE(PRT,2), SIZE(PRT,3)))
ALLOCATE(ZVTS(SIZE(PRT,1), SIZE(PRT,2), SIZE(PRT,3)))
ALLOCATE(ZVTG(SIZE(PRT,1), SIZE(PRT,2), SIZE(PRT,3)))
!
ZLAMBDAS(:,:,:) = 0.
ZLAMBDAG(:,:,:) = 0.
ZVTS(:,:,:) = 0.
ZVTG(:,:,:) = 0.
!
! Surface reference air density
ZRHO00 = XP00 / (XRD * XTHVREFZ(IKB))
!
! for snow
WHERE (PRT(:,:,:,5) .GT. 1.E-12)
  ZLAMBDAS(:,:,:) = MIN(XLBDAS_MAX, &
                        XLBS * (PRHODREF(:,:,:) * &
                        MAX(PRT(:,:,:,5), XRTMIN(5)))**XLBEXS)
  ZVTS(:,:,:) = XCS * MOMG(XALPHAS, XNUS, XBS+XDS) * ZLAMBDAS(:,:,:)**(-XDS) * &
               (ZRHO00 / PRHODREF(:,:,:))**XCEXVT / MOMG(XALPHAS, XNUS, XBS)
ELSEWHERE
  ZLAMBDAS(:,:,:) = 0.
  ZVTS(:,:,:) = 0.
END WHERE
!
! for graupel
WHERE(PRT(:,:,:,6) .GT. 1.E-12)
  ZLAMBDAG(:,:,:) = XLBG * (PRHODREF(:,:,:) * &
                  MAX(PRT(:,:,:,6), XRTMIN(6)))**XLBEXG
  ZVTG(:,:,:) = XCG * MOMG(XALPHAG, XNUG, XBG+XDG) * ZLAMBDAG(:,:,:)**(-XDG) * &
               (ZRHO00 / PRHODREF(:,:,:))**XCEXVT / MOMG(XALPHAG, XNUG, XBG)
ELSEWHERE
  ZLAMBDAG(:,:,:) = 0.
  ZVTG(:,:,:) = 0.
END WHERE
!
DEALLOCATE(ZLAMBDAS)
DEALLOCATE(ZLAMBDAG)
!
!
!-------------------------------------------------------------------------------
!
!*       2.     INITIALIZE THE VARIABLES
!   	        ------------------------
!
IF (GFIRSTCALL) THEN
  GFIRSTCALL = .FALSE.
!
  JCOUNT = 0
  ZMASS_C     = 0.
  ZMASS_R     = 0.
  ZMASS_I     = 0.
  ZMASS_S     = 0.
  ZMASS_G     = 0.
  ZMASS_ICE_P = 0.
  ZFLUX_PROD  = 0.
  ZFLUX_PRECIP  = 0.
  ZFLUX_NPRECIP = 0.
  ZVOL_UP5    = 0.
  ZVOL_UP10   = 0.
  ZVOL_G      = 0.
  ZWMAX       = 0.
  ZDBZMAX     = 0.
  ZCTH_REF    = 0.
  ZCTH_MR     = 0.
  ZCLD_VOL    = 0.
  ZINPRR      = 0.
  ZMAX_INPRR  = 0.
END IF
!
ZICE_MASS    = 0.
ZMASS_SP     = 0.  
ZMASS_GP     = 0.
ZFLUX_I      = 0.
ZFLUX_SP     = 0.
ZFLUX_SNP    = 0.
ZFLUX_G      = 0.
ZCLD_TOP_REF = 0.
ZCLD_TOP_MR  = 0.
!
!-------------------------------------------------------------------------------
!
!*       3.     COMPUTE THE DYNAMICAL AND MICROPHYSICAL PARAMETERS
!               --------------------------------------------------
!
JCOUNT = JCOUNT + 1
!
!*       3.1   compute the maximum vertical velocity
!
ZWMAX = ZWMAX + MAXVAL(PWT(IIB:IIE,IJB:IJE,IKB:IKE))
!
!
!*       3.2   compute the maximum radar reflectivity
!
CALL RADAR_RAIN_ICE (PRT, PCIT, PRHODREF, ZTCT, &
                     ZWORK31, ZWORK32, ZWORK33, ZWORK34)
!
ZDBZMAX = ZDBZMAX + MAXVAL(ZWORK31(IIB:IIE,IJB:IJE,IKB:IKE))
!  
!
!*       3.3    compute the mass of the different microphysical species
!
ZMASS_C = ZMASS_C + SUM(PRT(IIB:IIE,IJB:IJE,IKB:IKE,2) * &
                          PRHODJ(IIB:IIE,IJB:IJE,IKB:IKE))
!
ZMASS_R = ZMASS_R + SUM(PRT(IIB:IIE,IJB:IJE,IKB:IKE,3) * &
                          PRHODJ(IIB:IIE,IJB:IJE,IKB:IKE))
!
ZMASS_I = ZMASS_I + SUM(PRT(IIB:IIE,IJB:IJE,IKB:IKE,4) * &
                          PRHODJ(IIB:IIE,IJB:IJE,IKB:IKE))
!
ZMASS_S = ZMASS_S + SUM(PRT(IIB:IIE,IJB:IJE,IKB:IKE,5) * &
                          PRHODJ(IIB:IIE,IJB:IJE,IKB:IKE))
!
ZMASS_G = ZMASS_G + SUM(PRT(IIB:IIE,IJB:IJE,IKB:IKE,6) * &
                          PRHODJ(IIB:IIE,IJB:IJE,IKB:IKE))
!
!
!*       3.4    compute the ice mass fluxes
!
!*       3.4.1  non-precipitation ice mass flux
!
IFLAG(:,:) = 0
ICOUNT = 0
!
DO II = IIB, IIE
  DO IJ = IJB, IJE
    DO IK = IKB, IKE
!
!*       3.4.1  non-precipitation ice crystal mass flux
!
      IF (ZTCT(II,IJ,IK) .LT. 0. .AND. PWT(II,IJ,IK) .GT. 0.) THEN
        ZFLUX_I = ZFLUX_I + &
                  PWT(II,IJ,IK) * PRT(II,IJ,IK,4) * PRHODJ(II,IJ,IK)
      END IF
!
!*       3.4.2  non-precipitation snow mass flux
!
      IF (ZTCT(II,IJ,IK) .LT. 0. .AND. PWT(II,IJ,IK) .GT. ZVTS(II,IJ,IK)) THEN
        ZFLUX_SNP = ZFLUX_SNP + &
                    (PWT(II,IJ,IK) - ZVTS(II,IJ,IK)) * PRT(II,IJ,IK,5) * &
                    PRHODJ(II,IJ,IK)
      END IF
!
!*       3.4.3  precipitation snow mass flux
!
      IF (ZTCT(II,IJ,IK) .LT. 0. .AND. PWT(II,IJ,IK) .LT. ZVTS(II,IJ,IK)) THEN
        ZMASS_SP = ZMASS_SP + PRT(II,IJ,IK,5) * PRHODJ(II,IJ,IK)
        ZFLUX_SP = ZFLUX_SP + &
                  (PWT(II,IJ,IK) - ZVTS(II,IJ,IK)) * PRT(II,IJ,IK,5) * &
                   PRHODJ(II,IJ,IK)
      END IF
!
!*       3.4.4  precipitation graupel mass flux
!
      IF (ZTCT(II,IJ,IK) .LT. 0. .AND. PWT(II,IJ,IK) .LT. ZVTG(II,IJ,IK)) THEN
        ZMASS_GP = ZMASS_GP + PRT(II,IJ,IK,6) * PRHODJ(II,IJ,IK)
        ZFLUX_G = ZFLUX_G + &
                 (PWT(II,IJ,IK) - ZVTG(II,IJ,IK)) * PRT(II,IJ,IK,6) * &
                  PRHODJ(II,IJ,IK)
      END IF
!
!
!*       3.5   compute the updraft volume
!
! Updraft volume for W > 5 m/s
      IF (ZTCT(II,IJ,IK) .LT. -5. .AND. PWT(II,IJ,IK) .GT. 5.) THEN
        ZVOL_UP5 = ZVOL_UP5 + XDXHATM * XDYHATM * &
                              (PZZ(II,IJ,IK+1) - PZZ(II,IJ,IK-1)) / 2. 
      END IF
!
! Updraft volume for W > 10 m/s
      IF (ZTCT(II,IJ,IK) .LT. -5. .AND. PWT(II,IJ,IK) .GT. 10.) THEN
        ZVOL_UP10 = ZVOL_UP10 + XDXHATM * XDYHATM * &
                              (PZZ(II,IJ,IK+1) - PZZ(II,IJ,IK-1)) / 2.
      END IF
!
!
!*       3.6  total ice mass
!
      IF (ZTCT(II,IJ,IK) .LT. -10. .AND. ZWORK31(II,IJ,IK) .GT. 18.) THEN
        ZICE_MASS = ZICE_MASS + (PRT(II,IJ,IK,4) + PRT(II,IJ,IK,5) + PRT(II,IJ,IK,6)) * &
                    PRHODJ(II,IJ,IK)
        IFLAG(II,IJ) = IFLAG(II,IJ) + 1
      END IF
    END DO  ! end loop ik
!
    IF (IFLAG(II,IJ) .GE. 1) THEN
      ICOUNT = ICOUNT + 1
    END IF
  END DO   ! end loop ij
END DO   ! end loop ii
!
DEALLOCATE(ZVTS)
DEALLOCATE(ZVTG)
!
!
!*       3.7   precipitation and non precipitation ice mass flux product
!  
IF (ZFLUX_G .LT. 0. .AND. ZFLUX_I .GT. 0.) THEN
  ZFLUX_PROD = ZFLUX_PROD - (ZFLUX_I + ZFLUX_SNP) * (ZFLUX_G + ZFLUX_SP)
END IF
!
! precipitation ice mass flux
IF ((ZFLUX_G+ZFLUX_SP) .LT. 0.) THEN
  ZFLUX_PRECIP = ZFLUX_PRECIP - (ZFLUX_G + ZFLUX_SP)
END IF
!
! non-precipitation ice mass flux
IF ((ZFLUX_I+ZFLUX_SNP) .GT. 0.) THEN
  ZFLUX_NPRECIP = ZFLUX_NPRECIP + (ZFLUX_I + ZFLUX_SNP)
END IF
!
!
!*       3.8   compute the precipitation ice mass
!
IF ((ZMASS_GP .GT. 0.) .OR. (ZMASS_SP .GT. 0.)) THEN
  ZMASS_ICE_P = ZMASS_ICE_P + ZMASS_GP + ZMASS_SP
END IF
!
!
!*       3.9   compute the ice water path
!
CALL SUM_ELEC_ll(ZICE_MASS)
CALL SUM_ELEC_ll(ICOUNT)
!
IF (ICOUNT .GT. 0) THEN
  ZIWP = ZIWP + ZICE_MASS / (REAL(ICOUNT) * XDXHATM * XDYHATM)
END IF
!
!
!*       3.10  compute the cloud top height
!
DO II = IIB, IIE
  DO IJ = IJB, IJE
    DO IK = IKB, IKE
! maximum height of the 20 dBZ echo
      IF (ZWORK31(II,IJ,IK) .GT. 20. .AND. PZZ(II,IJ,IK) .GT. ZCLD_TOP_REF) THEN
        ZCLD_TOP_REF = PZZ(II,IJ,IK)
      END IF
!
! maximum height with mixing ratio > 1.e-4
      IF (ZCLOUD(II,IJ,IK) .GT. 1.E-4 .AND. PZZ(II,IJ,IK) .GT. ZCLD_TOP_REF) THEN
        ZCLD_TOP_MR = PZZ(II,IJ,IK)
      END IF
!
!
!*       3.11  compute the cloud volume
!
      IF (ZCLOUD(II,IJ,IK) .GT. 1.E-4) THEN
        ZCLD_VOL = ZCLD_VOL + XDXHATM * XDYHATM * &
                              (PZZ(II,IJ,IK+1) - PZZ(II,IJ,IK-1)) / 2.
      END IF
!
    END DO
  END DO
END DO
!
DEALLOCATE(ZCLOUD)
!
ZCTH_MR  = ZCTH_MR  + ZCLD_TOP_MR
ZCTH_REF = ZCTH_REF + ZCLD_TOP_REF
!
!
!*       3.12  compute the instantaneous precipitation rate
!
ZMAX_INPRR = ZMAX_INPRR + MAXVAL(PINPRR(IIB:IIE,IJB:IJE))
ZINPRR     = ZINPRR + SUM(PINPRR(IIB:IIE,IJB:IJE))
!
!-------------------------------------------------------------------------------
!
!*       4.     FROM LOCAL TO GLOBAL VARIABLES
!               ------------------------------
!
CALL MAX_ELEC_ll (ZCTH_REF, IPROC_MAX)
CALL MAX_ELEC_ll (ZCTH_MR,  IPROC_MAX)
CALL MAX_ELEC_ll (ZDBZMAX,  IPROC_MAX)
CALL MAX_ELEC_ll (ZMAX_INPRR,IPROC_MAX)
CALL MAX_ELEC_ll (ZWMAX,    IPROC_MAX)
!
!
!-------------------------------------------------------------------------------
!
!*       5.     SAVE THE DATA IN AN ASCII FILE
!               ------------------------------
!
CALL MYPROC_ELEC_ll(IPROC)
!
IF (JCOUNT == JCOUNT_STOP) THEN
!
  ZINPRR     = ZINPRR * 3.6E6          ! m/s --> mm/H
  ZMAX_INPRR = ZMAX_INPRR * 3.6E6      ! m/s --> mm/H
!
  CALL REDUCESUM_ll (ZVOL_UP5,      IINFO_ll)
  CALL REDUCESUM_ll (ZVOL_UP10,     IINFO_ll)
  CALL REDUCESUM_ll (ZMASS_C,       IINFO_ll)
  CALL REDUCESUM_ll (ZMASS_R,       IINFO_ll)
  CALL REDUCESUM_ll (ZMASS_I,       IINFO_ll)
  CALL REDUCESUM_ll (ZMASS_S,       IINFO_ll)
  CALL REDUCESUM_ll (ZMASS_G,       IINFO_ll)
  CALL REDUCESUM_ll (ZMASS_ICE_P,   IINFO_ll)
  CALL REDUCESUM_ll (ZFLUX_PROD,    IINFO_ll)
  CALL REDUCESUM_ll (ZFLUX_PRECIP,  IINFO_ll)
  CALL REDUCESUM_ll (ZFLUX_NPRECIP, IINFO_ll)
  CALL REDUCESUM_ll (ZCLD_VOL,      IINFO_ll)
  CALL REDUCESUM_ll (ZINPRR,        IINFO_ll)
!
  IF (IPROC == 0) THEN
    ILU = TPFILE_SERIES_CLOUD_ELEC%NLU
    WRITE (ILU, FMT='(I6,19(E12.4))') &
             INT(KTCOUNT*PTSTEP),         & ! time
             ZCTH_REF/REAL(JCOUNT),      & ! cloud top height from Z
             ZCTH_MR/REAL(JCOUNT),       & ! cloud top height from m.r.
             ZDBZMAX/REAL(JCOUNT),       & ! maximum radar reflectivity
             ZWMAX/REAL(JCOUNT),         & ! maximum vertical velocity
             ZVOL_UP5/REAL(JCOUNT),      & ! updraft volume for W > 5 m/s
             ZVOL_UP10/REAL(JCOUNT),     & ! updraft volume for W > 10 m/s
             ZMASS_C/REAL(JCOUNT),       & ! cloud droplets mass
             ZMASS_R/REAL(JCOUNT),       & ! rain mass
             ZMASS_I/REAL(JCOUNT),       & ! ice crystal mass
             ZMASS_S/REAL(JCOUNT),       & ! snow mass
             ZMASS_G/REAL(JCOUNT),       & ! graupel mass
             ZMASS_ICE_P/REAL(JCOUNT),   & ! precipitation ice mass
             ZFLUX_PROD/REAL(JCOUNT),    & ! ice mass flux product
             ZFLUX_PRECIP/REAL(JCOUNT),  & ! precipitation ice mass flux
             ZFLUX_NPRECIP/REAL(JCOUNT), & ! non-precipitation ice mass flux
             ZIWP/REAL(JCOUNT),          & ! ice water path
             ZCLD_VOL/REAL(JCOUNT),      & ! cloud volume
             ZINPRR/REAL(JCOUNT),        & ! Rain instant precip
             ZMAX_INPRR/REAL(JCOUNT)       ! maximum rain instant. precip.
    FLUSH(UNIT=ILU)
  END IF
!
  JCOUNT = 0
  ZMASS_C     = 0.
  ZMASS_R     = 0.
  ZMASS_I     = 0.
  ZMASS_S     = 0.
  ZMASS_G     = 0.
  ZMASS_ICE_P = 0.
  ZFLUX_PROD  = 0.
  ZFLUX_PRECIP  = 0.
  ZFLUX_NPRECIP = 0.
  ZVOL_UP5    = 0.
  ZVOL_UP10   = 0.
  ZWMAX       = 0.
  ZDBZMAX     = 0.
  ZCTH_REF    = 0.
  ZCTH_MR     = 0.
  ZIWP        = 0.
  ZCLD_VOL    = 0.
  ZINPRR      = 0.
  ZMAX_INPRR  = 0.
END IF
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
!     ##############################################
      FUNCTION MOMG0D(PALPHA, PNU, PP)  RESULT(PMOMG)
!     ##############################################
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PALPHA, PNU
REAL, INTENT(IN) :: PP
REAL             :: PMOMG
!
!
PMOMG = GAMMA(PNU+PP/PALPHA) / GAMMA(PNU)
!
END FUNCTION MOMG0D
!
!-------------------------------------------------------------------------------

!
END SUBROUTINE SERIES_CLOUD_ELEC
