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
TYPE(CONVPAR_SHAL),        INTENT(IN) :: CVP_SHAL
TYPE(CONVPAREXT),          INTENT(IN) :: CVPEXT
TYPE(CST_T),               INTENT(IN) :: CST
TYPE(DIMPHYEX_T),          INTENT(IN) :: D
REAL, DIMENSION(D%NIT),     INTENT(IN) :: PTKECLS   ! TKE CLS
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PTH, PTHV ! theta, theta_v
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PTHES     ! envir. satur. theta_e
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PRV       ! vapor mixing ratio
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PPRES     ! pressure
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PZ        ! height of grid point (m)
REAL, DIMENSION(D%NIT,D%NKT),INTENT(IN) :: PW        ! vertical velocity
!
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PTHLCL    ! theta at LCL
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PTLCL     ! temp. at LCL
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PRVLCL    ! vapor mixing ratio at  LCL
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PWLCL     ! parcel velocity at  LCL
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PZLCL     ! height at LCL (m)
REAL, DIMENSION(D%NIT),     INTENT(OUT):: PTHVELCL  ! environm. theta_v at LCL (K)
LOGICAL, DIMENSION(D%NIT),  INTENT(OUT):: OTRIG     ! logical mask for convection
INTEGER, DIMENSION(D%NIT),  INTENT(INOUT):: KLCL    ! contains vert. index of LCL
INTEGER, DIMENSION(D%NIT),  INTENT(INOUT):: KDPL    ! contains vert. index of DPL
INTEGER, DIMENSION(D%NIT),  INTENT(INOUT):: KPBL    ! contains index of source layer top
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JKK, JK, JKM, JL, JT! vertical loop index
INTEGER :: JI                  ! horizontal loop index
INTEGER :: IKB, IKE            ! horizontal + vertical loop bounds
REAL    :: ZEPS, ZEPSA         ! R_d / R_v, R_v / R_d
REAL    :: ZCPORD, ZRDOCP      ! C_pd / R_d, R_d / C_pd
REAL    :: ZX1, ZX2
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
REAL, DIMENSION(D%NIT) :: ZLVA, ZCPHA! specific heats of vaporisation, dry air
REAL, DIMENSION(D%NIT) :: ZLVB, ZCPHB! specific heats of vaporisation, dry air
REAL, DIMENSION(D%NIT) :: ZEWA, ZEWB ! vapor saturation mixing ratios
REAL, DIMENSION(D%NIT) :: ZLSA, ZLSB  ! latent heat L_s
REAL, DIMENSION(D%NIT) :: ZDP      ! pressure between LCL and model layer
REAL, DIMENSION(D%NIT) :: ZTOP,ZTOPP     ! estimated cloud top (m)
REAL, DIMENSION(D%NIT) :: ZWORK1, ZWORK2, ZWORK3    ! work arrays
LOGICAL, DIMENSION(D%NIT) :: GTRIG2          ! local arrays for OTRIG
LOGICAL, DIMENSION(D%NIT) :: GWORK1                 ! work array
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

REAL, DIMENSION(D%NIT,D%NKT) :: ZZZX1, ZZPPRES,ZZPTH,ZZPRV,ZZPTHESINV,ZZPZ
REAL, DIMENSION(D%NIT)       :: ZWLCLSQRENT
INTEGER  :: JLSTEP,JLSIZE,JLSTART
INTEGER  :: JLCLMIN   !!MIN value of LCL on all grid points, used to remove
                      !!unnecessary computation for smaller vertical levels.
INTEGER  :: IJIMIN,IJIMAX,IJIMIN2,IJIMAX2
LOGICAL  :: LLCOMPUTE
!
!-------------------------------------------------------------------------------
!
!*       0.3    Compute array bounds
!               --------------------
!

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
!       0.     Auxiliary arrays 
!              ----------------
!
!!!Auxiliary arrays, computed once here and used for each JKK in 3.
DO JK=IKB+1,IKE-1
  DO JI=D%NIB,D%NIE
    ZZZX1(JI,JK)=PPRES(JI,JK)-PPRES(JI,JK+1)
    ZZPPRES(JI,JK)=PPRES(JI,JK)*ZZZX1(JI,JK)
    ZZPTH(JI,JK)=PTH(JI,JK)*ZZZX1(JI,JK)
    ZZPRV(JI,JK)=MAX(0.,PRV(JI,JK))*ZZZX1(JI,JK)
  ENDDO
ENDDO
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
     DO JI=D%NIB, D%NIE
       GWORK1(JI) = ZZDPL(JI) - PZ(JI,IKB) < CVP_SHAL%XZLCL
     ENDDO
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
            ZDPTHMIX(JI) = ZDPTHMIX(JI) + ZZZX1(JI,JK)        !!!uses the auxiliary arrays of 0. 
            ZPRESMIX(JI) = ZPRESMIX(JI) + ZZPPRES(JI,JK) 
            ZTHLCL(JI)   = ZTHLCL(JI)   + ZZPTH(JI,JK)   
            ZRVLCL(JI)   = ZRVLCL(JI)   + ZZPRV(JI,JK) 
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
        ZX1 = LOG( ZEVMIX(JI) / 613.3 )
              ! dewpoint temperature
        ZX1 = ( 4780.8 - 32.19 * ZX1 ) / ( 17.502 - ZX1 )
              ! adiabatic saturation temperature
        ZTLCL(JI)  = ZX1 - ( .212 + 1.571E-3 * ( ZX1 - CST%XTT )      &
                   - 4.36E-4 * ( ZTMIX(JI) - CST%XTT ) ) * ( ZTMIX(JI) - ZX1 )
        ZTLCL(JI)  = MIN( ZTLCL(JI), ZTMIX(JI) )
        ZPLCL(JI)  = CST%XP00 * ( ZTLCL(JI) / ZTHLCL(JI) ) ** ZCPORD
!
     END IF
     ENDDO
!
     DO JI=D%NIB, D%NIE
     !arrays ZEWA, ZLVA, ZLSA, ZCPHA, 
     !       ZEWB, ZLVB, ZLSB, ZCPHB are only used if GWORK1(JI)=.TRUE.
     !in the rest of the code.
     IF ( GWORK1(JI) ) THEN
       CALL CONVECT_SATMIXRATIO( ZPLCL(JI), ZTLCL(JI), ZEPS, ZEWA(JI), ZLVA(JI), ZLSA(JI), ZCPHA(JI) )
       CALL CONVECT_SATMIXRATIO( ZPRESMIX(JI), ZTMIX(JI), ZEPS, ZEWB(JI), ZLVB(JI), ZLSB(JI), ZCPHB(JI) )
     ENDIF
     ENDDO
!
!*       4.2    Correct ZTLCL in order to be completely consistent
!               with MNH saturation formula
!               ---------------------------------------------
!
     DO JI=D%NIB, D%NIE
     IF( GWORK1(JI) ) THEN
        ZLSA(JI) = ZEWA(JI) / ZTLCL(JI) * ( CST%XBETAW / ZTLCL(JI) - CST%XGAMW ) ! dr_sat/dT
        ZLSA(JI) = ( ZEWA(JI) - ZRVLCL(JI) ) /                              &
                        ( 1. + ZLVA(JI) / ZCPHA(JI) * ZLSA(JI) )
        ZTLCL(JI)  = ZTLCL(JI) - ZLVA(JI) / ZCPHA(JI) * ZLSA(JI)
!
     END IF
     ENDDO
!
!
!*       4.3    If ZRVLCL = PRVMIX is oversaturated set humidity
!               and temperature to saturation values.
!               ---------------------------------------------
!
     DO JI=D%NIB, D%NIE
     IF( GWORK1(JI) .AND. ZRVLCL(JI) > ZEWB(JI) ) THEN
        ZLSB(JI) = ZEWB(JI) / ZTMIX(JI) * ( CST%XBETAW / ZTMIX(JI) - CST%XGAMW ) ! dr_sat/dT
        ZLSB(JI) = ( ZEWB(JI) - ZRVLCL(JI) ) /                              &
                       ( 1. + ZLVB(JI) / ZCPHB(JI) * ZLSB(JI) )
        ZTLCL(JI)  = ZTMIX(JI) - ZLVB(JI) / ZCPHB(JI) * ZLSB(JI)
        ZRVLCL(JI) = ZRVLCL(JI) - ZLSB(JI)
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
    !!!Every 4 vertical level, the first loop checks if computations need 
    !!!to be done in the second loop for a given JI.
    !!!PPRES comes from array PREHYDF, which is decreasing with
    !!!increasing vertical level.
    JLSTEP=4
    DO JLSTART=JKK,IKE-1,JLSTEP
      JLSIZE=MIN(JLSTEP,JT-JLSTART+1) 
      IF (JLSTART .EQ. JKK) THEN
        IJIMIN=D%NIB
        IJIMAX=D%NIE
        LLCOMPUTE=.TRUE.
      ELSE
        IF (LLCOMPUTE) THEN
          LLCOMPUTE=.FALSE.
          IJIMIN2=IJIMIN
          IJIMAX2=IJIMAX
          DO JI=D%NIB,D%NIE
            IF ((JI .GE. IJIMIN2) .AND. (JI .LE. IJIMAX2) .AND. GWORK1(JI) .AND. (ZPLCL(JI) <= PPRES(JI,JLSTART))) THEN
              IF (.NOT. LLCOMPUTE) IJIMIN=JI      
              IJIMAX=JI 
              LLCOMPUTE=.TRUE.           
            ENDIF
          ENDDO
        ENDIF
      ENDIF

      IF (LLCOMPUTE) THEN
        DO JK = JLSTART, JLSTART+JLSIZE-1
           DO JI = D%NIB,D%NIE
             IF ( (JI .GE. IJIMIN) .AND. (JI .LE. IJIMAX) .AND. ZPLCL(JI) <= PPRES(JI,JK) .AND. GWORK1(JI) ) ILCL(JI) = JK + 1
           END DO
        END DO
      ENDIF
    ENDDO

!
!
!*        5.2   Estimate height and environm. theta_v at LCL
!               --------------------------------------------------
!
    JLCLMIN=ILCL(D%NIB)
    DO JI = D%NIB, D%NIE
      IF ( GWORK1(JI) ) THEN
        JK   = ILCL(JI)
        JKM  = JK - 1
        ZDP(JI)    = LOG( ZPLCL(JI) / PPRES(JI,JKM) ) /                     &
                     LOG( PPRES(JI,JK) / PPRES(JI,JKM) )
        ZTHVELCL(JI) = PTHV(JI,JKM) + ( PTHV(JI,JK) - PTHV(JI,JKM) ) * ZDP(JI)
           ! we compute the precise value of the LCL
           ! The precise height is between the levels ILCL and ILCL-1.
        ZZLCL(JI) = PZ(JI,JKM) + ( PZ(JI,JK) - PZ(JI,JKM) ) * ZDP(JI)
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
    DO JI = D%NIB, D%NIE
      JLCLMIN=MIN(JLCLMIN,ILCL(JI))  
      ZWLCL(JI) = CVP_SHAL%XAW * MAX(0.,PW(JI,IKB)) + CVP_SHAL%XBW
                 ! the factor 1.05 takes entrainment into account
      ZWLCLSQRENT(JI) = 1.05*ZWLCL(JI)*ZWLCL(JI)
!
!
!*       6.3    Look for parcel that produces sufficient cloud depth.
!               The cloud top is estimated as the level where the CAPE
!               is smaller  than a given value (based on vertical velocity eq.)
!               --------------------------------------------------------------
!
      ZTHEUL(JI) = ZTLCL(JI) * ( ZTHLCL(JI) / ZTLCL(JI) )                       &
                                             ** ( 1. - 0.28 * ZRVLCL(JI) )  &
                          * EXP( ( 3374.6525 / ZTLCL(JI) - 2.5403 ) *       &
                               ZRVLCL(JI) * ( 1. + 0.81 * ZRVLCL(JI) ) )
    END DO

     ZCAPE(D%NIB:D%NIE)  = MAX(JLCLMIN-1-IKB,0)*CST%XG 
     ZCAP(D%NIB:D%NIE) = 0.
     ZTOP(D%NIB:D%NIE)  = 0.
     ZWORK3(D%NIB:D%NIE)= 0.
     JKM=MAX(IKB,JLCLMIN-1) 

!REK : should try if's instead of sign & max

     !!!Every 4 vertical level, the first loop checks if computations need to be done in the second loop.
     !!!If ZWORK3(JI)==-1 and ZCAPE(JI) >10., no more computation are necessary for JI :
     !!!if ZWORK3(JI)==-1, it remains equal to -1, and ZTOP(JI) does not change.
     !!!if ZCAPE(JI)>10., it remains >10., because MAX(1.,ZX1) is a positive
     !!!number ; the last loop of CONVECT_TRIGGER_SHAL uses ZCAPE(JI) as a
     !!!criterion. 
     JLSTEP=4
     DO JLSTART=JKM,JT,JLSTEP
       JLSIZE=MIN(JLSTEP,JT-JLSTART+1)

       IF (JLSTART .EQ. JKM) THEN
         LLCOMPUTE=.TRUE.
       ELSE
         IF (LLCOMPUTE) THEN
           LLCOMPUTE=.FALSE.
           DO JI=D%NIB,D%NIE
             IF (GTRIG2(JI) .AND. ((ZWORK3(JI) .GT. -0.5) .OR. (ZCAPE(JI) .LE. 10.))) LLCOMPUTE=.TRUE.
           ENDDO
         ENDIF
      ENDIF
      IF (LLCOMPUTE) THEN       

         DO JL = JLSTART,JLSTART+JLSIZE-1
           JK = JL + 1
  
           DO JI = D%NIB, D%NIE
             !!!it was not possible to use a mask IF (JI>=IJIMIN .AND. JI
             !!!<=IJIMAX) here, as the loop did not vectorize.
             ZX1 = ( 2. * ZTHEUL(JI) / ( PTHES(JI,JK) + PTHES(JI,JL) ) - 1. ) *&
                    & (PZ(JI,JK) - PZ(JI,JL))
             IF  ((JL < ILCL(JI))) ZX1 = 0.
             ZCAPE(JI)  = ZCAPE(JI) + CST%XG * MAX( 1., ZX1 )
             ZCAP(JI)   = ZCAP(JI) + ZX1
                 ! the factor 1.05 in ZWLCLSQRENT takes entrainment into account
             ZX2 = SIGN( 1., ( CVP_SHAL%XNHGAM * CST%XG * ZCAP(JI) + ZWLCLSQRENT(JI) ) )
             ZWORK3(JI) = MAX( -1., (ZWORK3(JI) + MIN(0.,ZX2)) )
                 ! Nota, the factors ZX2 and ZWORK3 are only used to avoid
                 ! if and goto statements, the difficulty is to extract only
                 ! the level where the criterium is first fullfilled
             ZTOPP(JI)=ZTOP(JI)
             ZTOP(JI)   = PZ(JI,JL) * .5 * ( 1. + ZX2 ) * ( 1. + ZWORK3(JI) ) + &
                        & ZTOP(JI) * .5 * ( 1. - ZX2 )
             ZTOP(JI)=MAX(ZTOP(JI),ZTOPP(JI))
           ENDDO

         ENDDO

      ENDIF !!LLCOMPUTE

     END DO
!
!

     DO JI=D%NIB, D%NIE
     IF( (ZTOP(JI)-ZZLCL(JI)) .GE. CVP_SHAL%XCDEPTH .AND. GTRIG2(JI) .AND. ZCAPE(JI) > 10. )THEN
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
CONTAINS
INCLUDE "convect_satmixratio.h"
!
END SUBROUTINE CONVECT_TRIGGER_SHAL
