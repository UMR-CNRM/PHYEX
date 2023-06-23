!     ######spl
      SUBROUTINE CONVECT_TRIGGER_SHAL(  CVP_SHAL, CVPEXT, CST, D,  &
                                        PPRES, PTH, PTHV, PTHES,             &
                                        PRV, PW, PZ, PTKECLS,                &
                                        PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, &
                                        PTHVELCL, KLCL, KDPL, KPBL, OTRIG)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ########################################################################
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
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XCDEPTH_D          ! maximum allowed cloud depth
!!          XDTPERT            ! add small Temp peturbation
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!          XAW, XBW, XATPERT, XBTPERT
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
!!   F. Bouyssel    05/11/08  Modifications for reproductibility
!!   E. Bazile      05/05/09  Modifications for using really W and the tempe.
!!                            perturbation function of the TKE.
!!   F. Bouyssel    08/11/13  Modifications for reproductibility

!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : CST_T
USE MODD_CONVPAR_SHAL, ONLY : CONVPAR_SHAL
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CONVPAR_SHAL)                         ,INTENT(IN)     :: CVP_SHAL
TYPE(CONVPAREXT)                           ,INTENT(IN)     :: CVPEXT
TYPE(CST_T)                                ,INTENT(IN)     :: CST
TYPE(DIMPHYEX_T)                           ,INTENT(IN)     :: D
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PPRES     ! pressure
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PTH,PTHV ! theta, theta_v
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PTHES     ! envir. satur. theta_e
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PRV       ! vapor mixing ratio
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PW        ! vertical velocity
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PZ        ! height of grid point (m)
REAL               ,DIMENSION(D%NIT)       ,INTENT(IN)     :: PTKECLS   ! TKE CLS
!
REAL               ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: PTHLCL    ! theta at LCL
REAL               ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: PTLCL     ! temp. at LCL
REAL               ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: PRVLCL    ! vapor mixing ratio at  LCL
REAL               ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: PWLCL     ! parcel velocity at  LCL
REAL               ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: PZLCL     ! height at LCL (m)
REAL               ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: PTHVELCL  ! environm. theta_v at LCL (K)
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(INOUT)  :: KLCL    ! contains vert. index of LCL
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(INOUT)  :: KDPL    ! contains vert. index of DPL
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(INOUT)  :: KPBL    ! contains index of source layer top
LOGICAL            ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: OTRIG     ! logical mask for convection
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JKK, JK, JKM, JL, JT! vertical loop index
INTEGER :: JI                  ! horizontal loop index
INTEGER :: IKB, IKE            ! horizontal + vertical loop bounds
REAL    :: ZEPS, ZEPSA         ! R_d / R_v, R_v / R_d
REAL    :: ZCPORD, ZRDOCP      ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(D%NIT) :: ZTHLCL, ZTLCL, ZRVLCL, & ! locals for PTHLCL,PTLCL
                               ZWLCL,  ZZLCL, ZTHVELCL  ! PRVLCL, ....
INTEGER, DIMENSION(D%NIT) :: IDPL, IPBL, ILCL      ! locals for KDPL, ...
REAL, DIMENSION(D%NIT) :: ZPLCL    ! pressure at LCL
REAL, DIMENSION(D%NIT) :: ZZDPL    ! height of DPL
REAL, DIMENSION(D%NIT) :: ZTHVLCL  ! theta_v at LCL = mixed layer value
REAL, DIMENSION(D%NIT) :: ZTMIX    ! mixed layer temperature
REAL, DIMENSION(D%NIT) :: ZEVMIX   ! mixed layer water vapor pressure
REAL, DIMENSION(D%NIT) :: ZDPTHMIX, ZPRESMIX ! mixed layer depth and pressure
REAL, DIMENSION(D%NIT) :: ZCAPE    ! convective available energy (m^2/s^2/g)
REAL, DIMENSION(D%NIT) :: ZCAP     ! pseudo fro CAPE
REAL, DIMENSION(D%NIT) :: ZTHEUL   ! updraft equiv. pot. temperature (K)
REAL, DIMENSION(D%NIT) :: ZLV, ZCPH! specific heats of vaporisation, dry air
REAL, DIMENSION(D%NIT) :: ZDP      ! pressure between LCL and model layer
REAL, DIMENSION(D%NIT) :: ZTOP,ZTOPP     ! estimated cloud top (m)
REAL, DIMENSION(D%NIT) :: ZWORK1, ZWORK2, ZWORK3    ! work arrays
LOGICAL, DIMENSION(D%NIT) :: GTRIG2          ! local arrays for OTRIG
LOGICAL, DIMENSION(D%NIT) :: GWORK1                 ! work array
!
!
!-------------------------------------------------------------------------------
!
!*       0.3    Compute array bounds
!               --------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "convect_satmixratio.h"

IF (LHOOK) CALL DR_HOOK('CONVECT_TRIGGER_SHAL',0,ZHOOK_HANDLE)
IKB = 1 + CVPEXT%JCVEXB
IKE = D%NKT - CVPEXT%JCVEXT
!
!
!*       1.     Initialize local variables
!               --------------------------
!
ZEPS       = CST%XRD / CST%XRV
ZEPSA      = CST%XRV / CST%XRD
ZCPORD     = CST%XCPD / CST%XRD
ZRDOCP     = CST%XRD / CST%XCPD
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
!
!
!
!       1.     Determine highest necessary loop test layer
!              -------------------------------------------
!
JT = IKE - 2
!
!*       2.     Enter loop for convection test
!               ------------------------------
!
DO JKK = IKB + 1, IKE - 2
!
     GWORK1(D%NIB:D%NIE) = ZZDPL(D%NIB:D%NIE) - PZ(D%NIB:D%NIE,IKB) < CVP_SHAL%XZLCL
          ! we exit the trigger test when the center of the mixed layer is more
          ! than 1500 m  above soil level.
     DO JI=D%NIB, D%NIE
       IF ( GWORK1(JI) ) THEN
          ZDPTHMIX(JI) = 0.
          ZPRESMIX(JI) = 0.
          ZTHLCL(JI)   = 0.
          ZRVLCL(JI)   = 0.
          ZZDPL(JI)    = PZ(JI,JKK)
          IDPL(JI)     = JKK
       END IF
     ENDDO
!
!
!*       3.     Construct a mixed layer of at least 50 hPa (XZPBL)
!               ------------------------------------------
!
     DO JK = JKK, IKE - 1
       JKM = JK + 1
       DO JI = D%NIB, D%NIE
         IF ( GWORK1(JI) .AND. ZDPTHMIX(JI) < CVP_SHAL%XZPBL ) THEN
            IPBL(JI)     = JK
            ZWORK1(JI)   = PPRES(JI,JK) - PPRES(JI,JKM)
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZWORK1(JI)
            ZPRESMIX(JI) = ZPRESMIX(JI) + PPRES(JI,JK) * ZWORK1(JI)
            ZTHLCL(JI)   = ZTHLCL(JI)   + PTH(JI,JK)   * ZWORK1(JI)
            ZRVLCL(JI)   = ZRVLCL(JI)   + MAX(0., PRV(JI,JK))   * ZWORK1(JI)
         END IF
       END DO
     END DO
!
!
     DO JI=D%NIB, D%NIE
     IF ( GWORK1(JI) ) THEN
!
        ZPRESMIX(JI) = ZPRESMIX(JI) / ZDPTHMIX(JI)
        ZTHLCL(JI)   = ZTHLCL(JI)   / ZDPTHMIX(JI) + &
      & (CVP_SHAL%XATPERT * MIN(3.,PTKECLS(JI))/CST%XCPD +CVP_SHAL%XBTPERT) * CVP_SHAL%XDTPERT ! add small Temp Perturb.
        ZRVLCL(JI)   = ZRVLCL(JI)   / ZDPTHMIX(JI)
        ZTHVLCL(JI)  = ZTHLCL(JI) * ( 1. + ZEPSA * ZRVLCL(JI) )                 &
                    / ( 1. + ZRVLCL(JI) )
!
!*       4.1    Use an empirical direct solution ( Bolton formula )
!               to determine temperature and pressure at LCL.
!               NotaJI the adiabatic saturation temperature is not
!                     equal to the dewpoint temperature
!               ----------------------------------------------------
!
!
        ZTMIX(JI)  = ZTHLCL(JI) * ( ZPRESMIX(JI) / CST%XP00 ) ** ZRDOCP
        ZEVMIX(JI) = ZRVLCL(JI) * ZPRESMIX(JI) / ( ZRVLCL(JI) + ZEPS )
        ZEVMIX(JI) = MAX( 1.E-8, ZEVMIX(JI) )
        ZWORK1(JI) = LOG( ZEVMIX(JI) / 613.3 )
              ! dewpoint temperature
        ZWORK1(JI) = ( 4780.8 - 32.19 * ZWORK1(JI) ) / ( 17.502 - ZWORK1(JI) )
              ! adiabatic saturation temperature
        ZTLCL(JI)  = ZWORK1(JI) - ( .212 + 1.571E-3 * ( ZWORK1(JI) - CST%XTT )      &
                   - 4.36E-4 * ( ZTMIX(JI) - CST%XTT ) ) * ( ZTMIX(JI) - ZWORK1(JI) )
        ZTLCL(JI)  = MIN( ZTLCL(JI), ZTMIX(JI) )
        ZPLCL(JI)  = CST%XP00 * ( ZTLCL(JI) / ZTHLCL(JI) ) ** ZCPORD
!
     END IF
     ENDDO
!
!
!*       4.2    Correct ZTLCL in order to be completely consistent
!               with MNH saturation formula
!               ---------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( CST, D, ZPLCL, ZTLCL, ZWORK1, ZLV, ZWORK2, ZCPH )
     DO JI=D%NIB, D%NIE
     IF( GWORK1(JI) ) THEN
        ZWORK2(JI) = ZWORK1(JI) / ZTLCL(JI) * ( CST%XBETAW / ZTLCL(JI) - CST%XGAMW ) ! dr_sat/dT
        ZWORK2(JI) = ( ZWORK1(JI) - ZRVLCL(JI) ) /                              &
                        ( 1. + ZLV(JI) / ZCPH(JI) * ZWORK2(JI) )
        ZTLCL(JI)  = ZTLCL(JI) - ZLV(JI) / ZCPH(JI) * ZWORK2(JI)
!
     END IF
     ENDDO
!
!
!*       4.3    If ZRVLCL = PRVMIX is oversaturated set humidity
!               and temperature to saturation values.
!               ---------------------------------------------
!
     CALL CONVECT_SATMIXRATIO( CST, D, ZPRESMIX, ZTMIX, ZWORK1, ZLV, ZWORK2, ZCPH )
     DO JI=D%NIB, D%NIE
     IF( GWORK1(JI) .AND. ZRVLCL(JI) > ZWORK1(JI) ) THEN
        ZWORK2(JI) = ZWORK1(JI) / ZTMIX(JI) * ( CST%XBETAW / ZTMIX(JI) - CST%XGAMW ) ! dr_sat/dT
        ZWORK2(JI) = ( ZWORK1(JI) - ZRVLCL(JI) ) /                              &
                       ( 1. + ZLV(JI) / ZCPH(JI) * ZWORK2(JI) )
        ZTLCL(JI)  = ZTMIX(JI) - ZLV(JI) / ZCPH(JI) * ZWORK2(JI)
        ZRVLCL(JI) = ZRVLCL(JI) - ZWORK2(JI)
        ZPLCL(JI)  = ZPRESMIX(JI)
        ZTHLCL(JI) = ZTLCL(JI) * ( CST%XP00 / ZPLCL(JI) ) ** ZRDOCP
        ZTHVLCL(JI)= ZTHLCL(JI) * ( 1. + ZEPSA * ZRVLCL(JI) )                   &
                              / ( 1. + ZRVLCL(JI) )
     END IF
     ENDDO
!
!
!*        5.1   Determine  vertical loop index at the LCL and DPL
!               --------------------------------------------------
!
    DO JK = JKK, IKE - 1
       DO JI = D%NIB, D%NIE
         IF ( ZPLCL(JI) <= PPRES(JI,JK) .AND. GWORK1(JI) ) ILCL(JI) = JK + 1
       END DO
    END DO
!
!
!*        5.2   Estimate height and environm. theta_v at LCL
!               --------------------------------------------------
!
    DO JI = D%NIB, D%NIE
        JK   = ILCL(JI)
        JKM  = JK - 1
        ZDP(JI)    = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                     LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZWORK1(JI) = PTHV(JI,JKM) + ( PTHV(JI,JK) - PTHV(JI,JKM) ) * ZDP(JI)
           ! we compute the precise value of the LCL
           ! The precise height is between the levels ILCL and ILCL-1.
        ZWORK2(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
    END DO
    DO JI = D%NIB, D%NIE
    IF( GWORK1(JI) ) THEN
        ZTHVELCL(JI) = ZWORK1(JI)
        ZZLCL(JI)    = ZWORK2(JI)
    END IF
    END DO
!
!
!*       6.     Check to see if cloud is bouyant
!               --------------------------------
!
!*      6.1    Compute grid scale vertical velocity perturbation term ZWORK1
!               -------------------------------------------------------------
!
!            !  normalize w grid scale to a 25 km refer. grid
!    DO JI = 1, D%NIT
!       JK  = ILCL(JI)
!       JKM = JK - 1
!       ZWORK1(JI) =  ( PW(JI,JKM)  + ( PW(JI,JK) - PW(JI,JKM) ) * ZDP(JI) )  &
!                          * SQRT( PDXDY(JI) / XA25 )
!                         - 0.02 * ZZLCL(JI) / XZLCL ! avoid spurious convection
!    END DO
!            ! compute sign of normalized grid scale w
!       ZWORK2(:) = SIGN( 1., ZWORK1(:) )
!       ZWORK1(:) = XWTRIG * ZWORK2(:) * ABS( ZWORK1(:) ) ** 0.333       &
!                          * ( XP00 / ZPLCL(:) ) ** ZRDOCP
!
!*       6.2    Compute parcel vertical velocity at LCL
!               ---------------------------------------
!
!    DO JI = 1, D%NIT
!       JKDL = IDPL(JI)
!       ZWORK3(JI) = XG * ZWORK1(JI) * ( ZZLCL(JI) - PZ(JI,JKDL) )       &
!                      / ( PTHV(JI,JKDL) + ZTHVELCL(JI) )
!    END DO
!    WHERE( GWORK1(:) )
!      ZWLCL(:)  = 1. + .5 * ZWORK2(:) * SQRT( ABS( ZWORK3(:) ) )
!      GTRIG(:)  = ZTHVLCL(:) - ZTHVELCL(:) + ZWORK1(:) > 0. .AND.       &
!                  ZWLCL(:) > 0.
!    END WHERE
     ZWLCL(D%NIB:D%NIE) = CVP_SHAL%XAW * MAX(0.,PW(D%NIB:D%NIE,IKB)) + CVP_SHAL%XBW
!
!
!*       6.3    Look for parcel that produces sufficient cloud depth.
!               The cloud top is estimated as the level where the CAPE
!               is smaller  than a given value (based on vertical velocity eq.)
!               --------------------------------------------------------------
!
     ZTHEUL(D%NIB:D%NIE) = ZTLCL(D%NIB:D%NIE) * ( ZTHLCL(D%NIB:D%NIE) / ZTLCL(D%NIB:D%NIE) )                       &
                                             ** ( 1. - 0.28 * ZRVLCL(D%NIB:D%NIE) )  &
                          * EXP( ( 3374.6525 / ZTLCL(D%NIB:D%NIE) - 2.5403 ) *       &
                               ZRVLCL(D%NIB:D%NIE) * ( 1. + 0.81 * ZRVLCL(D%NIB:D%NIE) ) )
!
     ZCAPE(D%NIB:D%NIE) = 0.
     ZCAP(D%NIB:D%NIE)  = 0.
     ZTOP(D%NIB:D%NIE)  = 0.
     ZTOPP(D%NIB:D%NIE)  = 0.
     ZWORK3(D%NIB:D%NIE)= 0.
     JKM = IKB
     DO JL = JKM, JT
        JK = JL + 1
        DO JI = D%NIB, D%NIE
           ZWORK1(JI) = ( 2. * ZTHEUL(JI) /                                &
            ( PTHES(JI,JK) + PTHES(JI,JL) ) - 1. ) * ( PZ(JI,JK) - PZ(JI,JL) )
           IF ( JL < ILCL(JI) ) ZWORK1(JI) = 0.
           ZCAPE(JI)  = ZCAPE(JI) + CST%XG * MAX( 1., ZWORK1(JI) )
           ZCAP(JI)   = ZCAP(JI) + ZWORK1(JI)
           ZWORK2(JI) = CVP_SHAL%XNHGAM * CST%XG * ZCAP(JI) + 1.05 * ZWLCL(JI) * ZWLCL(JI)
               ! the factor 1.05 takes entrainment into account
           ZWORK2(JI) = SIGN( 1., ZWORK2(JI) )
           ZWORK3(JI) = ZWORK3(JI) + MIN(0., ZWORK2(JI) )
           ZWORK3(JI) = MAX( -1., ZWORK3(JI) )
               ! Nota, the factors ZWORK2 and ZWORK3 are only used to avoid
               ! if and goto statements, the difficulty is to extract only
               ! the level where the criterium is first fullfilled
           ZTOPP(JI)=ZTOP(JI)
           ZTOP(JI)   = PZ(JI,JL) * .5 * ( 1. + ZWORK2(JI) ) * ( 1. + ZWORK3(JI) ) + &
                        ZTOP(JI) * .5 * ( 1. - ZWORK2(JI) )
           ZTOP(JI)=MAX(ZTOP(JI),ZTOPP(JI))
           ZTOPP(JI)=ZTOP(JI)
         END DO
     END DO
!
!
     ZWORK2(D%NIB:D%NIE) = ZTOP(D%NIB:D%NIE) - ZZLCL(D%NIB:D%NIE)
   ! WHERE( ZWORK2(:)  .GE. XCDEPTH  .AND. ZWORK2(:) < XCDEPTH_D .AND. GTRIG2(:) &
     DO JI=D%NIB, D%NIE
     IF( ZWORK2(JI) .GE. CVP_SHAL%XCDEPTH .AND. GTRIG2(JI) .AND. ZCAPE(JI) > 10. )THEN
        GTRIG2(JI)   = .FALSE.
        OTRIG(JI)    = .TRUE.
      ! OTRIG(JI)    = GTRIG(JI)     ! we  select the first departure level
        PTHLCL(JI)   = ZTHLCL(JI)    ! that gives sufficient cloud depth
        PRVLCL(JI)   = ZRVLCL(JI)
        PTLCL(JI)    = ZTLCL(JI)
        PWLCL(JI)    = ZWLCL(JI)
        PZLCL(JI)    = ZZLCL(JI)
        PTHVELCL(JI) = ZTHVELCL(JI)
        KDPL(JI)     = IDPL(JI)
        KPBL(JI)     = IPBL(JI)
        KLCL(JI)     = ILCL(JI)
     END IF
     ENDDO
!
END DO
!
!
IF (LHOOK) CALL DR_HOOK('CONVECT_TRIGGER_SHAL',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_TRIGGER_SHAL

