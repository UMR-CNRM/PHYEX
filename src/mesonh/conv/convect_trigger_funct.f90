!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 conv 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #################
      MODULE MODI_CONVECT_TRIGGER_FUNCT
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_TRIGGER_FUNCT( KLON, KLEV,                           &
                                        PPRES, PTH, PTHV, PTHES,              &
                                        PRV, PW, PZ, PDXDY,                   &
                                        PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL,  &
                                        PTHVELCL, KLCL, KDPL, KPBL, OTRIG,    &
                                        PCAPE )
!
INTEGER, INTENT(IN)                   :: KLON      ! horizontal loop index
INTEGER, INTENT(IN)                   :: KLEV      ! vertical loop index
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY     ! grid area
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH, PTHV ! theta, theta_v
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHES     ! envir. satur. theta_e
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRV       ! vapor mixing ratio 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PPRES     ! pressure
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PZ        ! height of grid point (m)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PW        ! vertical velocity
!
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHLCL    ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PTLCL     ! temp. at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PWLCL     ! parcel velocity at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PZLCL     ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(OUT):: OTRIG     ! logical mask for convection 
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KLCL    ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KDPL    ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KPBL    ! contains index of source layer top
REAL, DIMENSION(KLON),     INTENT(OUT):: PCAPE     ! CAPE (J/kg) for diagnostics
!
END SUBROUTINE CONVECT_TRIGGER_FUNCT
!
END INTERFACE
!
END MODULE MODI_CONVECT_TRIGGER_FUNCT
!     #########################################################################
      SUBROUTINE CONVECT_TRIGGER_FUNCT( KLON, KLEV,                           &
                                        PPRES, PTH, PTHV, PTHES,              &
                                        PRV, PW, PZ, PDXDY,                   &
                                        PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL,  &
                                        PTHVELCL, KLCL, KDPL, KPBL, OTRIG,    &
                                        PCAPE )
!     #########################################################################
!
!!**** Determine convective columns as well as the cloudy values of theta,
!!     and qv at the lifting condensation level (LCL) 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine convective columns
!!   
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      What we look for is the undermost unstable level at each grid point.
!!      
!!     
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_SATMIXRATIO
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! Reference pressure
!!          XRD, XRV           ! Gaz  constants for dry air and water vapor
!!          XCPD               ! Cpd (dry air)
!!          XTT                ! triple point temperature
!!          XBETAW, XGAMW      ! constants for vapor saturation pressure
!!
!!      Module MODD_CONVPAR
!!          XA25               ! reference grid area
!!          XZLCL              ! maximum height difference between
!!                             ! the surface and the DPL
!!          XZPBL              ! minimum mixed layer depth to sustain convection
!!          XWTRIG             ! constant in vertical velocity trigger
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book2 of documentation ( routine TRIGGER_FUNCT)
!!      Fritsch and Chappell (1980), J. Atm. Sci., Vol. 37, 1722-1761.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95 
!!   Last modified  20/03/97  Select first departure level
!!                            that produces a cloud thicker than XCDEPTH
!!   Last modified  12/12/97  add small perturbation
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAR
USE MODD_CONVPAREXT
USE MODI_CONVECT_SATMIXRATIO
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN)                   :: KLON      ! horizontal loop index
INTEGER, INTENT(IN)                   :: KLEV      ! vertical loop index
REAL, DIMENSION(KLON),     INTENT(IN) :: PDXDY     ! grid area
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTH, PTHV ! theta, theta_v
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PTHES     ! envir. satur. theta_e
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PRV       ! vapor mixing ratio 
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PPRES     ! pressure
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PZ        ! height of grid point (m)
REAL, DIMENSION(KLON,KLEV),INTENT(IN) :: PW        ! vertical velocity
!
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHLCL    ! theta at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PTLCL     ! temp. at LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PWLCL     ! parcel velocity at  LCL
REAL, DIMENSION(KLON),     INTENT(OUT):: PZLCL     ! height at LCL (m)
REAL, DIMENSION(KLON),     INTENT(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(KLON),  INTENT(OUT):: OTRIG     ! logical mask for convection 
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KLCL    ! contains vert. index of LCL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KDPL    ! contains vert. index of DPL
INTEGER, DIMENSION(KLON),  INTENT(INOUT):: KPBL    ! contains index of source layer top
REAL, DIMENSION(KLON),     INTENT(OUT):: PCAPE     ! CAPE (J/kg) for diagnostics
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JKK, JK, JKP, JKM, JKDL, JL, JKT, JT! vertical loop index
INTEGER :: JI                                  ! horizontal loop index 
INTEGER :: IIE, IKB, IKE                       ! horizontal + vertical loop bounds
REAL    :: ZEPS, ZEPSA                         ! R_d / R_v, R_v / R_d 
REAL    :: ZCPORD, ZRDOCP                      ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(KLON) :: ZTHLCL, ZTLCL, ZRVLCL, & ! locals for PTHLCL,PTLCL
                               ZWLCL,  ZZLCL, ZTHVELCL  ! PRVLCL, ....
INTEGER, DIMENSION(KLON) :: IDPL, IPBL, ILCL      ! locals for KDPL, ...
REAL, DIMENSION(KLON) :: ZPLCL    ! pressure at LCL
REAL, DIMENSION(KLON) :: ZZDPL    ! height of DPL 
REAL, DIMENSION(KLON) :: ZTHVLCL  ! theta_v at LCL = mixed layer value
REAL, DIMENSION(KLON) :: ZTMIX    ! mixed layer temperature
REAL, DIMENSION(KLON) :: ZEVMIX   ! mixed layer water vapor pressure 
REAL, DIMENSION(KLON) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL, DIMENSION(KLON) :: ZCAPE    ! convective available energy (m^2/s^2/g)
REAL, DIMENSION(KLON) :: ZTHEUL   ! updraft equiv. pot. temperature (K)
REAL, DIMENSION(KLON) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL, DIMENSION(KLON) :: ZDP      ! pressure between LCL and model layer
REAL, DIMENSION(KLON) :: ZTOP     ! estimated cloud top (m)
REAL, DIMENSION(KLON,KLEV):: ZCAP ! CAPE at every level for diagnostics
REAL, DIMENSION(KLON) :: ZWORK1, ZWORK2, ZWORK3    ! work arrays
LOGICAL, DIMENSION(KLON) :: GTRIG, GTRIG2          ! local arrays for OTRIG
LOGICAL, DIMENSION(KLON) :: GWORK1                 ! work array
!
!
!-------------------------------------------------------------------------------
!
!*       0.3    Compute array bounds
!               --------------------
!
IIE = KLON
IKB = 1 + JCVEXB 
IKE = KLEV - JCVEXT 
!
!
!*       1.     Initialize local variables
!               --------------------------
!
ZEPS       = XRD / XRV
ZEPSA      = XRV / XRD 
ZCPORD     = XCPD / XRD
ZRDOCP     = XRD / XCPD
OTRIG(:)   = .FALSE.
IDPL(:)    = KDPL(:)
IPBL(:)    = KPBL(:)
ILCL(:)    = KLCL(:)
PWLCL(:)   = 0.
ZWLCL(:)   = 0.
PTHLCL(:)  = 1.
PTHVELCL(:)= 1.
PTLCL(:)   = 1.
PRVLCL(:)  = 0.
PWLCL(:)   = 0.
PZLCL(:)   = PZ(:,IKB)
ZZDPL(:)   = PZ(:,IKB)
GTRIG2(:)  = .TRUE.
ZCAP(:,:)  = 0.
!
!
!
!       1.     Determine highest necessary loop test layer
!              -------------------------------------------
!
JT = IKE - 2
DO JK = IKB + 1, IKE - 2
  IF ( PZ(1,JK) - PZ(1,IKB) < 12.E3 ) JT = JK 
END DO
!
!
!*       2.     Enter loop for convection test
!               ------------------------------
!
JKP = MINVAL( IDPL(:) ) + 1
JKT = JT
DO JKK = JKP, JKT
!
     GWORK1(:) = ZZDPL(:) - PZ(:,IKB) < XZLCL
          ! we exit the trigger test when the center of the mixed layer is more
          ! than 3500 m  above soil level.
     WHERE ( GWORK1(:) )
        ZDPTHMIX(:) = 0.
        ZPRESMIX(:) = 0.
        ZTHLCL(:)   = 0.
        ZRVLCL(:)   = 0.
        ZZDPL(:)    = PZ(:,JKK)
        IDPL(:)     = JKK
     END WHERE
!
!
!*       3.     Construct a mixed layer of at least 60 hPa (XZPBL)
!               ------------------------------------------
!
     DO JK = JKK, IKE - 1
       JKM = JK + 1
       DO JI = 1, IIE     
         IF ( GWORK1(JI) .AND. ZDPTHMIX(JI) < XZPBL ) THEN
            IPBL(JI)     = JK
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            ZTHLCL(JI)   = ZTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            ZRVLCL(JI)   = ZRVLCL(JI)   + PRV(JI,JK)   * ZWORK1(JI)
         END IF
       END DO
        IF ( MINVAL ( ZDPTHMIX(:) ) >= XZPBL ) EXIT
     END DO
!
!
     WHERE ( GWORK1(:) )
!
        ZPRESMIX(:) = ZPRESMIX(:) / ZDPTHMIX(:)
        ZTHLCL(:)   = ZTHLCL(:)   / ZDPTHMIX(:) +.3   ! add small perturbation
        ZRVLCL(:)   = ZRVLCL(:)   / ZDPTHMIX(:) +1.e-4
        ZTHVLCL(:)  = ZTHLCL(:) * ( 1. + ZEPSA * ZRVLCL(:) )                 &
				/ ( 1. + ZRVLCL(:) )
!
!*       4.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL. 
!               Nota: the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               ----------------------------------------------------
!
! 
        ZTMIX(:)  = ZTHLCL(:) * ( ZPRESMIX(:) / XP00 ) ** ZRDOCP 
        ZEVMIX(:) = ZRVLCL(:) * ZPRESMIX(:) / ( ZRVLCL(:) + ZEPS )
        ZEVMIX(:) = MAX( 1.E-8, ZEVMIX(:) )
        ZWORK1(:) = LOG( ZEVMIX(:) / 613.3 )
              ! dewpoint temperature
        ZWORK1(:) = ( 4780.8 - 32.19 * ZWORK1(:) ) / ( 17.502 - ZWORK1(:) ) 
              ! adiabatic saturation temperature
        ZTLCL(:)  = ZWORK1(:) - ( .212 + 1.571E-3 * ( ZWORK1(:) - XTT )      &
                   - 4.36E-4 * ( ZTMIX(:) - XTT ) ) * ( ZTMIX(:) - ZWORK1(:) )
        ZTLCL(:)  = MIN( ZTLCL(:), ZTMIX(:) )
        ZPLCL(:)  = XP00 * ( ZTLCL(:) / ZTHLCL(:) ) ** ZCPORD
!
     END WHERE
!
!
!*       4.2    Correct ZTLCL in order to be completely consistent
!               with MNH saturation formula
!               ---------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( KLON, ZPLCL, ZTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( GWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTLCL(:) * ( XBETAW / ZTLCL(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
                        ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) 
        ZTLCL(:)  = ZTLCL(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
!
     END WHERE
!
!
!*       4.3    If ZRVLCL = PRVMIX is oversaturated set humidity 
!               and temperature to saturation values. 
!               ---------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( KLON, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     WHERE( GWORK1(:) .AND. ZRVLCL(:) > ZWORK1(:) )
        ZWORK2(:) = ZWORK1(:) / ZTMIX(:) * ( XBETAW / ZTMIX(:) - XGAMW ) ! dr_sat/dT
        ZWORK2(:) = ( ZWORK1(:) - ZRVLCL(:) ) /                              &
                       ( 1. + ZLV(:) / ZCPH(:) * ZWORK2(:) ) 
        ZTLCL(:)  = ZTMIX(:) - ZLV(:) / ZCPH(:) * ZWORK2(:)
        ZRVLCL(:) = ZRVLCL(:) - ZWORK2(:)
        ZPLCL(:)  = ZPRESMIX(:)
        ZTHLCL(:) = ZTLCL(:) * ( XP00 / ZPLCL(:) ) ** ZRDOCP
        ZTHVLCL(:)= ZTHLCL(:) * ( 1. + ZEPSA * ZRVLCL(:) )                   &
                              / ( 1. + ZRVLCL(:) )
     END WHERE
!
!
!*        5.1   Determine  vertical loop index at the LCL and DPL
!               --------------------------------------------------
!
    DO JK = JKK, IKE - 1
       DO JI = 1, IIE
         IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. GWORK1(JI) ) ILCL(JI) = JK + 1
       END DO
    END DO
!
!
!*        5.2   Estimate height and environm. theta_v at LCL
!               --------------------------------------------------
!
    DO JI = 1, IIE
        JK   = ILCL(JI)
        JKM  = JK - 1
        ZDP(JI)    = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                     LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI) = PTHV(JI,JKM) + ( PTHV(JI,JK) - PTHV(JI,JKM) ) * ZDP(JI) 
           ! we compute the precise value of the LCL
           ! The precise height is between the levels ILCL and ILCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    END DO
    WHERE( GWORK1(:) )
        ZTHVELCL(:) = ZWORK1(:) 
        ZZLCL(:)    = ZWORK2(:)
    END WHERE
!        
!
!*       6.     Check to see if cloud is bouyant 
!               --------------------------------
!
!*      6.1    Compute grid scale vertical velocity perturbation term ZWORK1
!               -------------------------------------------------------------
! 
             !  normalize w grid scale to a 25 km refer. grid
     DO JI = 1, IIE
        JK  = ILCL(JI)
        JKM = JK - 1 
        JKDL= IDPL(JI)
       !ZWORK1(JI) =  ( PW(JI,JKM)  + ( PW(JI,JK) - PW(JI,JKM) ) * ZDP(JI) )  &
        ZWORK1(JI) =  ( PW(JI,JK)  +  PW(JI,JKDL)*ZZLCL(JI)/PZ(JI,JKDL) ) * .5 &   
                           * SQRT( PDXDY(JI) / XA25 )
!                         - 0.02 * ZZLCL(JI) / XZLCL ! avoid spurious convection
     END DO
             ! compute sign of normalized grid scale w
        ZWORK2(:) = SIGN( 1., ZWORK1(:) ) 
        ZWORK1(:) = XWTRIG * ZWORK2(:) * ABS( ZWORK1(:) ) ** 0.333       &
                           * ( XP00 / ZPLCL(:) ) ** ZRDOCP
!
!*       6.2    Compute parcel vertical velocity at LCL
!               ---------------------------------------
!                   
     DO JI = 1, IIE
        JKDL = IDPL(JI)
        ZWORK3(JI) = XG * ZWORK1(JI) * ( ZZLCL(JI) - PZ(JI,JKDL) )       &
                       / ( PTHV(JI,JKDL) + ZTHVELCL(JI) )
     END DO
     WHERE( GWORK1(:) )
       ZWLCL(:)  = 1. + .5 * ZWORK2(:) * SQRT( ABS( ZWORK3(:) ) ) 
       GTRIG(:)  = ZTHVLCL(:) - ZTHVELCL(:) + ZWORK1(:) > 0. .AND.       &
                   ZWLCL(:) > 0. 
     END WHERE
!
!
!*       6.3    Look for parcel that produces sufficient cloud depth.
!               The cloud top is estimated as the level where the CAPE 
!               is smaller  than a given value (based on vertical velocity eq.)
!               --------------------------------------------------------------
!
     ZTHEUL(:) = ZTLCL(:) * ( ZTHLCL(:) / ZTLCL(:) )                       &
                                             ** ( 1. - 0.28 * ZRVLCL(:) )  &
                          * EXP( ( 3374.6525 / ZTLCL(:) - 2.5403 ) *       &
                               ZRVLCL(:) * ( 1. + 0.81 * ZRVLCL(:) ) )
!
     ZCAPE(:) = 0.
     ZTOP(:)  = 0.
     ZWORK3(:)= 0.
     JKM = MINVAL( ILCL(:) )
     DO JL = JKM, JT
        JK = JL + 1
        DO JI = 1, IIE
           ZWORK1(JI) = ( 2. * ZTHEUL(JI) /                                &
            ( PTHES(JI,JK) + PTHES(JI,JL) ) - 1. ) * ( PZ(JI,JK) - PZ(JI,JL) )
           IF ( JL < ILCL(JI) ) ZWORK1(JI) = 0.
        !  IF ( JL <= ILCL(JI) ) ZWORK1(JI) = 0.
           ZCAPE(JI)  = ZCAPE(JI) + ZWORK1(JI)
           ZCAP(JI,JKK) = ZCAP(JI,JKK) + XG * MAX( 0., ZWORK1(JI) ) ! actual CAPE
           ZWORK2(JI) = XNHGAM * XG * ZCAPE(JI) + 1.05 * ZWLCL(JI) * ZWLCL(JI)
               ! the factor 1.05 takes entrainment into account
           ZWORK2(JI) = SIGN( 1., ZWORK2(JI) )
           ZWORK3(JI) = ZWORK3(JI) + MIN(0., ZWORK2(JI) )
           ZWORK3(JI) = MAX( -1., ZWORK3(JI) )
               ! Nota, the factors ZWORK2 and ZWORK3 are only used to avoid
               ! if and goto statements, the difficulty is to extract only
               ! the level where the criterium is first fullfilled
           ZTOP(JI)   = PZ(JI,JL) * .5 * ( 1. + ZWORK2(JI) ) * ( 1. + ZWORK3(JI) ) + &
                        ZTOP(JI) * .5 * ( 1. - ZWORK2(JI) )
         END DO
     END DO
!
!
     WHERE( ZTOP(:) - ZZLCL(:) .GE. XCDEPTH  .AND. GTRIG(:) .AND. GTRIG2(:) )
        GTRIG2(:)   = .FALSE.
        OTRIG(:)    = GTRIG(:)     ! we  select the first departure level
        PTHLCL(:)   = ZTHLCL(:)    ! that gives sufficient cloud depth
        PRVLCL(:)   = ZRVLCL(:)
        PTLCL(:)    = ZTLCL(:)
        PWLCL(:)    = ZWLCL(:)
        PZLCL(:)    = ZZLCL(:)
        PTHVELCL(:) = ZTHVELCL(:)
        KDPL(:)     = IDPL(:)
        KPBL(:)     = IPBL(:)
        KLCL(:)     = ILCL(:)
     END WHERE
!
END DO
!
     DO JI = 1, IIE
       PCAPE(JI) = MAXVAL( ZCAP(JI,:) ) ! maximum CAPE for diagnostics
     END DO
!
!
END SUBROUTINE CONVECT_TRIGGER_FUNCT
