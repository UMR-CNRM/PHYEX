!     ######spl
    SUBROUTINE CONVECT_UPDRAFT_SHAL( CVP_SHAL, CVPEXT, CST, D, CONVPAR,              &
                                     KICE, PPRES, PDPRES, PZ, PTHL, PTHV, PTHES, PRW,&
                                     PTHLCL, PTLCL, PRVLCL, PWLCL, PZLCL, PTHVELCL,  &
                                     PMFLCL, OTRIG, KLCL, KDPL, KPBL,                &
                                     PUMF, PUER, PUDR, PUTHL, PUTHV, PURW,           &
                                     PURC, PURI, PCAPE, KCTL, KETL,GTRIG1 )
    USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!    ###############################################################################
!
!!**** Compute updraft properties from DPL to CTL.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine updraft properties
!!      ( mass flux, thermodynamics, precipitation )
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!
!!
!!
!!    EXTERNAL
!!    --------
!!     Routine CONVECT_MIXING_FUNCT
!!     Routine CONVECT_CONDENS
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                 ! gravity constant
!!          XP00               ! reference pressure
!!          XRD, XRV           ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV, XCL    ! Cp of dry air, water vapor and liquid water
!!          XTT                ! triple point temperature
!!          XLVTT              ! vaporisation heat at XTT
!!
!!
!!      Module MODD_CONVPAR_SHAL
!!          XA25               ! reference grid area
!!          XCRAD              ! cloud radius
!!          XCDEPTH            ! minimum necessary cloud depth
!!          XENTR              ! entrainment constant
!!          XNHGAM             ! coefficient for buoyancy term in w eq.
!!                             ! accounting for nh-pressure
!!          XTFRZ1             ! begin of freezing interval
!!          XTFRZ2             ! begin of freezing interval
!!
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_UPDRAFT)
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/11/95
!!   Last modified  10/12/97
!!   F. Bouyssel    05/11/08  Modifications for reproductibility
!!   F. Bouyssel    08/11/13  Modifications for reproductibility
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : CST_T
USE MODD_CONVPAR, ONLY : CONVPAR_T
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
TYPE(CONVPAR_T)                            ,INTENT(IN)     :: CONVPAR
INTEGER                                    ,INTENT(IN)     :: KICE  ! flag for ice ( 1 = yes,
                                                !                0 = no ice )
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PPRES ! pressure (P)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PDPRES! pressure difference between
                                                ! bottom and top of layer (Pa)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PZ    ! height of model layer (m)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PTHL  ! grid scale enthalpy (J/kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PTHV  ! grid scale theta_v
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PTHES ! grid scale saturated theta_e
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(IN)     :: PRW   ! grid scale total water
                                                ! mixing ratio
REAL               ,DIMENSION(D%NIT)       ,INTENT(IN)     :: PTHLCL ! theta at LCL
REAL               ,DIMENSION(D%NIT)       ,INTENT(IN)     :: PTLCL  ! temp. at LCL
REAL               ,DIMENSION(D%NIT)       ,INTENT(IN)     :: PRVLCL ! vapor mixing ratio at  LCL
REAL               ,DIMENSION(D%NIT)       ,INTENT(IN)     :: PWLCL  ! parcel velocity at LCL (m/s)
REAL               ,DIMENSION(D%NIT)       ,INTENT(IN)     :: PZLCL  ! height at LCL (m)
REAL               ,DIMENSION(D%NIT)       ,INTENT(IN)     :: PTHVELCL  ! environm. theta_v at LCL (K)
REAL                                       ,INTENT(IN)     :: PMFLCL ! cloud  base unit mass flux
                                                ! (kg/s)
LOGICAL            ,DIMENSION(D%NIT)       ,INTENT(INOUT)  :: OTRIG! logical mask for convection
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: KLCL   ! contains vert. index of LCL
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: KDPL   ! contains vert. index of DPL
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: KPBL   !  " vert. index of source layertop
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PUMF  ! updraft mass flux (kg/s)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PUER  ! updraft entrainment (kg/s)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PUDR  ! updraft detrainment (kg/s)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PUTHL ! updraft enthalpy (J/kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PUTHV ! updraft theta_v (K)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PURW  ! updraft total water (kg/kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PURC  ! updraft cloud water (kg/kg)
REAL               ,DIMENSION(D%NIT,D%NKT) ,INTENT(OUT)    :: PURI  ! updraft cloud ice   (kg/kg)
REAL               ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: PCAPE  ! available potent. energy
!
!
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: KCTL   ! contains vert. index of CTL
INTEGER            ,DIMENSION(D%NIT)       ,INTENT(OUT)    :: KETL   ! contains vert. index of        &
                                                !equilibrium (zero buoyancy) level
LOGICAL            ,DIMENSION(D%NIT)       ,INTENT(IN)     :: GTRIG1! logical mask for convection
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB, IKE  ! horizontal and vertical loop bounds
INTEGER :: JI             ! horizontal loop index
INTEGER :: JK, JKP, JKM, JK1, JK2   ! vertical loop index
REAL    :: ZEPSA          ! R_v / R_d, C_pv / C_pd
REAL    :: ZRDOCP         ! C_pd / R_d, R_d / C_pd
!
REAL, DIMENSION(D%NIT)    :: ZUT             ! updraft temperature (K)
REAL, DIMENSION(D%NIT)    :: ZUW1, ZUW2      ! square of updraft vert.
                                            ! velocity at levels k and k+1
REAL, DIMENSION(D%NIT)    :: ZE1,ZE2,ZD1,ZD2 ! fractional entrainm./detrain
                                            ! rates at levels k and k+1
REAL, DIMENSION(D%NIT)    :: ZMIXF           ! critical mixed fraction
REAL, DIMENSION(D%NIT)    :: ZCPH            ! specific heat C_ph
REAL, DIMENSION(D%NIT)    :: ZLV, ZLS        ! latent heat of vaporis., sublim.
REAL, DIMENSION(D%NIT)    :: ZURV            ! updraft water vapor at level k+1
REAL, DIMENSION(D%NIT)    :: ZPI             ! Pi=(P0/P)**(Rd/Cpd)
REAL, DIMENSION(D%NIT)    :: ZTHEUL          ! theta_e for undilute ascent
REAL, DIMENSION(D%NIT)    :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5,   &
                            ZWORK6          ! work arrays
INTEGER, DIMENSION(D%NIT) :: IWORK           ! wok array
LOGICAL, DIMENSION(D%NIT) :: GWORK1, GWORK2, GWORK4
                                            ! work arrays
!
!
!-------------------------------------------------------------------------------
!
!        0.3   Set loop bounds
!              ---------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "convect_condens.h"
#include "convect_mixing_funct.h"

IF (LHOOK) CALL DR_HOOK('CONVECT_UPDRAFT_SHAL',0,ZHOOK_HANDLE)
IKB = 1 + CVPEXT%JCVEXB
IKE = D%NKT - CVPEXT%JCVEXT
!
!
!*       1.     Initialize updraft properties and local variables
!               -------------------------------------------------
!
ZEPSA      = CST%XRV / CST%XRD
ZRDOCP     = CST%XRD / CST%XCPD
!
PUMF(:,:)  = 0.
PUER(:,:)  = 0.
PUDR(:,:)  = 0.
PUTHL(:,:) = 0.
PUTHV(:,:) = 0.
PURW(:,:)  = 0.
PURC(:,:)  = 0.
PURI(:,:)  = 0.
ZUW1(:)    = PWLCL(:) * PWLCL(:)
ZUW2(:)    = 0.
ZE1(:)     = 0.
ZD1(:)     = 0.
PCAPE(:)   = 0.
KCTL(:)    = IKB
KETL(:)    = KLCL(:)
GWORK2(:)  = .TRUE.
ZPI(:)     = 1.
ZWORK3(:)  = 0.
ZWORK4(:)  = 0.
ZWORK5(:)  = 0.
ZWORK6(:)  = 0.
GWORK1(:)  = .FALSE.
GWORK4(:)  = .FALSE.
!
!
!*       1.1    Compute undilute updraft theta_e for CAPE computations
!               Bolton (1980) formula.
!               Define accurate enthalpy for updraft
!               -----------------------------------------------------
!
ZTHEUL(D%NIB:D%NIE) = PTLCL(D%NIB:D%NIE) * ( PTHLCL(D%NIB:D%NIE) / PTLCL(D%NIB:D%NIE) ) ** ( 1. - 0.28 * PRVLCL(D%NIB:D%NIE) )  &
            * EXP( ( 3374.6525 / PTLCL(D%NIB:D%NIE) - 2.5403 ) *                        &
                                   PRVLCL(D%NIB:D%NIE) * ( 1. + 0.81 * PRVLCL(D%NIB:D%NIE) ) )
!
!
ZWORK1(D%NIB:D%NIE) = ( CST%XCPD + PRVLCL(D%NIB:D%NIE) * CST%XCPV ) * PTLCL(D%NIB:D%NIE)                            &
            + ( 1. + PRVLCL(D%NIB:D%NIE) ) * CST%XG * PZLCL(D%NIB:D%NIE)
!
!
!*       2.     Set updraft properties between DPL and LCL
!               ------------------------------------------
!
JKP=IKE
JKM=IKB

DO JK = JKM, JKP
   DO JI = D%NIB, D%NIE
   IF ( JK >= KDPL(JI) .AND. JK < KLCL(JI) ) THEN
        PUMF(JI,JK)  = PMFLCL
        PUTHL(JI,JK) = ZWORK1(JI)
        PUTHV(JI,JK) = PTHLCL(JI) * ( 1. + ZEPSA * PRVLCL(JI) ) /             &
                                  ( 1. + PRVLCL(JI) )
        PURW(JI,JK)  = PRVLCL(JI)
   END IF
   END DO
END DO
!
!
!*       3.     Enter loop for updraft computations
!               ------------------------------------
!
DO JK = IKB + 1, IKE - 1
  ZWORK6(:) = 1.
  JKP = JK + 1
!
  GWORK4(D%NIB:D%NIE) = JK >= KLCL(D%NIB:D%NIE) - 1
  GWORK1(D%NIB:D%NIE) = GWORK4(D%NIB:D%NIE) .AND. GWORK2(D%NIB:D%NIE) ! this mask is used to confine
                           ! updraft computations between the LCL and the CTL
!
  DO JI=D%NIB, D%NIE
    IF( JK == KLCL(JI) - 1 ) ZWORK6(JI) = 0. ! factor that is used in buoyancy
  ENDDO                                      ! computation at first level above LCL
!
!
!*       4.     Estimate condensate, L_v L_i, Cph and theta_v at level k+1
!               ----------------------------------------------------------
!
    ZWORK1(D%NIB:D%NIE) = PURC(D%NIB:D%NIE,JK)
    ZWORK2(D%NIB:D%NIE) = PURI(D%NIB:D%NIE,JK)
    CALL CONVECT_CONDENS(CST, D, CONVPAR, KICE, PPRES(D%NIB:D%NIE,JKP),&
                         PUTHL(D%NIB:D%NIE,JK), PURW(D%NIB:D%NIE,JK),  &
                         ZWORK1, ZWORK2, PZ(D%NIB:D%NIE,JKP), ZUT,ZURV,&
                         PURC(D%NIB:D%NIE,JKP), PURI(D%NIB:D%NIE,JKP), &
                         ZLV, ZLS, ZCPH )
!
!
  ZPI(D%NIB:D%NIE) = ( CST%XP00 / PPRES(D%NIB:D%NIE,JKP) ) ** ZRDOCP
  DO JI=D%NIB, D%NIE
    IF ( GWORK1(JI) ) THEN
!
      PUTHV(JI,JKP) = ZPI(JI) * ZUT(JI) * ( 1. + ZEPSA * ZURV(JI) )           &
                           / ( 1. + PURW(JI,JK) )
!
!
!*         5.     Compute square of vertical velocity using entrainment
!                 at level k
!                 -----------------------------------------------------
!
      ZWORK3(JI) = PZ(JI,JKP) - PZ(JI,JK) * ZWORK6(JI) -         &
                       ( 1. - ZWORK6(JI) ) * PZLCL(JI)          ! level thickness
      ZWORK4(JI) = PTHV(JI,JK) * ZWORK6(JI) +                   &
                   ( 1. - ZWORK6(JI) ) * PTHVELCL(JI)
      ZWORK5(JI) = 2. * ZUW1(JI) * PUER(JI,JK) / MAX( .1, PUMF(JI,JK) )
      ZUW2(JI)   = ZUW1(JI) + ZWORK3(JI) * CVP_SHAL%XNHGAM * CST%XG *        &
                    ( ( PUTHV(JI,JK) + PUTHV(JI,JKP) ) /       &
                    ( ZWORK4(JI) + PTHV(JI,JKP) ) - 1. )       & ! buoyancy term
                  - ZWORK5(JI)                                  ! entrainment term
!
!
!*         6.     Update total precipitationJI dr_r=(r_c+r_i)*exp(-rate*dz)
!                 --------------------------------------------------------
!
!                      compute level mean vertical velocity
      ZWORK2(JI)   = 0.5 *                                                    &
                         ( SQRT( MAX( 1.E-2, ZUW2(JI) ) ) +                   &
                           SQRT( MAX( 1.E-2, ZUW1(JI) ) ) )
!
!
!*         7.     Update r_c, r_i, enthalpy, r_w  for precipitation
!                 -------------------------------------------------------
!
      PURW(JI,JKP)  = PURW(JI,JK)
      PURC(JI,JKP)  = PURC(JI,JKP)
      PURI(JI,JKP)  = PURI(JI,JKP)
      PUTHL(JI,JKP) = PUTHL(JI,JK)
!
      ZUW1(JI)      = ZUW2(JI)
!
    END IF
  ENDDO
!
!
!*       8.     Compute entrainment and detrainment using conservative
!               variables adjusted for precipitation ( not for entrainment)
!               -----------------------------------------------------------
!
!*       8.1    Compute critical mixed fraction by estimating unknown
!               T^mix r_c^mix and r_i^mix from enthalpy^mix and r_w^mix
!               We determine the zero crossing of the linear curve
!               evaluating the derivative using ZMIXF=0.1.
!               -----------------------------------------------------
!
    ZMIXF(D%NIB:D%NIE)  = 0.1   ! starting value for critical mixed fraction
    ZWORK1(D%NIB:D%NIE) = ZMIXF(D%NIB:D%NIE) * PTHL(D%NIB:D%NIE,JKP)                                     &
                     + ( 1. - ZMIXF(D%NIB:D%NIE) ) * PUTHL(D%NIB:D%NIE,JKP) ! mixed enthalpy
    ZWORK2(D%NIB:D%NIE) = ZMIXF(D%NIB:D%NIE) * PRW(D%NIB:D%NIE,JKP)                                      &
                     + ( 1. - ZMIXF(D%NIB:D%NIE) ) * PURW(D%NIB:D%NIE,JKP)  ! mixed r_w
!
    CALL CONVECT_CONDENS(CST, D, CONVPAR, KICE, PPRES(D%NIB:D%NIE,JKP),&
                         ZWORK1, ZWORK2, PURC(D%NIB:D%NIE,JKP),        &
                         PURI(D%NIB:D%NIE,JKP), PZ(D%NIB:D%NIE,JKP),   &
                         ZUT, ZWORK3, ZWORK4, ZWORK5, ZLV, ZLS, ZCPH)
!        put in enthalpy and r_w and get T r_c, r_i (ZUT, ZWORK4-5)
!
     ! compute theta_v of mixture
    ZWORK3(D%NIB:D%NIE) = ZUT(D%NIB:D%NIE) * ZPI(D%NIB:D%NIE) * ( 1. + ZEPSA * (                         &
                ZWORK2(D%NIB:D%NIE) - ZWORK4(D%NIB:D%NIE) - ZWORK5(D%NIB:D%NIE) ) ) / ( 1. + ZWORK2(D%NIB:D%NIE) )
     ! compute final value of critical mixed fraction using theta_v
     ! of mixture, grid-scale and updraft
    ZMIXF(D%NIB:D%NIE) = MAX( 0., PUTHV(D%NIB:D%NIE,JKP) - PTHV(D%NIB:D%NIE,JKP) ) * ZMIXF(D%NIB:D%NIE) /          &
                              ( PUTHV(D%NIB:D%NIE,JKP) - ZWORK3(D%NIB:D%NIE) + 1.E-10 )
    ZMIXF(D%NIB:D%NIE) = MAX( 0., MIN( 1., ZMIXF(D%NIB:D%NIE) ) )
!
!
!*       8.2     Compute final midlevel values for entr. and detrainment
!                after call of distribution function
!                -------------------------------------------------------
!
!
    CALL CONVECT_MIXING_FUNCT ( D, ZMIXF, 1, ZE2, ZD2 )
!       NoteD%NIB:D%NIE routine MIXING_FUNCT returns fractional entrainm/detrainm. rates
!
  ZE2=MIN(ZD2,MAX(.3,ZE2))
!
! ZWORK1(D%NIB:D%NIE) = XENTR * PMFLCL * PDPRES(D%NIB:D%NIE,JKP) / XCRAD ! rate of env. inflow
!*MOD
  zwork1(D%NIB:D%NIE) = CVP_SHAL%xentr * CST%xg / CVP_SHAL%xcrad * pumf(D%NIB:D%NIE,jk) * ( pz(D%NIB:D%NIE,jkp) - pz(D%NIB:D%NIE,jk) )
! ZWORK1(D%NIB:D%NIE) = XENTR * pumf(D%NIB:D%NIE,jk) * PDPRES(D%NIB:D%NIE,JKP) / XCRAD ! rate of env. inflow
!*MOD
  ZWORK2(:) = 0.
  DO JI=D%NIB, D%NIE
    IF( GWORK1(JI) ) ZWORK2(JI) = 1.
  ENDDO
  DO JI=D%NIB, D%NIE
  IF ( PUTHV(JI,JKP) > PTHV(JI,JKP) ) THEN
    PUER(JI,JKP) = 0.5 * ZWORK1(JI) * ( ZE1(JI) + ZE2(JI) ) * ZWORK2(JI)
    PUDR(JI,JKP) = 0.5 * ZWORK1(JI) * ( ZD1(JI) + ZD2(JI) ) * ZWORK2(JI)
  ELSE
    PUER(JI,JKP) = 0.
    PUDR(JI,JKP) = ZWORK1(JI) * ZWORK2(JI)
  END IF
  ENDDO
!
!*       8.3     Determine equilibrium temperature level
!                --------------------------------------
!
  DO JI=D%NIB, D%NIE
    IF ( PUTHV(JI,JKP) > PTHV(JI,JKP) .AND. JK > KLCL(JI) + 1 .AND. GWORK1(JI) )THEN
         KETL(JI) = JKP            ! equilibrium temperature level
    END IF
  ENDDO
!
!*       8.4     If the calculated detrained mass flux is greater than
!                the total updraft mass flux, or vertical velocity is
!                negative, all cloud mass detrains at previous model level,
!                exit updraft calculations - CTL is attained
!                -------------------------------------------------------
!
  DO JI=D%NIB, D%NIE
    IF( GWORK1(JI) ) THEN
      GWORK2(JI) = PUMF(JI,JK) - PUDR(JI,JKP) > 10. .AND. ZUW2(JI) > 0.
    ENDIF
  ENDDO
  DO JI=D%NIB, D%NIE
    IF ( GWORK2(JI) ) KCTL(JI) = JKP   ! cloud top level
  ENDDO
  GWORK1(D%NIB:D%NIE) = GWORK2(D%NIB:D%NIE) .AND. GWORK4(D%NIB:D%NIE)
!
  !IF ( COUNT( GWORK2(:) ) == 0 ) EXIT        
!
!
!*       9.   Compute CAPE for undilute ascent using theta_e and
!             theta_es instead of theta_v. This estimation produces
!             a significantly larger value for CAPE than the actual one.
!             ----------------------------------------------------------
!
  DO JI=D%NIB, D%NIE
    IF ( GWORK1(JI) )THEN
!
      ZWORK3(JI)   = PZ(JI,JKP) - PZ(JI,JK) * ZWORK6(JI) -                      &
                    ( 1. - ZWORK6(JI) ) *  PZLCL(JI)              ! level thickness
      ZWORK2(JI)   = PTHES(JI,JK) + ( 1. - ZWORK6(JI) ) *                      &
       ( PTHES(JI,JKP) - PTHES(JI,JK) ) / ( PZ(JI,JKP) - PZ(JI,JK) ) *          &
       ( PZLCL(JI) - PZ(JI,JK) ) ! linear interpolation for theta_es at LCL
                              ! ( this is only done for model level just above LCL
!
      ZWORK1(JI) = ( 2. * ZTHEUL(JI) ) / ( ZWORK2(JI) + PTHES(JI,JKP) ) - 1.
      PCAPE(JI)  = PCAPE(JI) + CST%XG * ZWORK3(JI) * MAX( 0., ZWORK1(JI) )
!
!
!*         10.   Compute final values of updraft mass flux, enthalpy, r_w
!                at level k+1
!                --------------------------------------------------------
!
      PUMF(JI,JKP)  = PUMF(JI,JK) - PUDR(JI,JKP) + PUER(JI,JKP)
      PUMF(JI,JKP)  = MAX( PUMF(JI,JKP), 0.1 )
      PUTHL(JI,JKP) = ( PUMF(JI,JK) * PUTHL(JI,JK) +                              &
                       PUER(JI,JKP) * PTHL(JI,JK) - PUDR(JI,JKP) * PUTHL(JI,JK) )  &
                      / PUMF(JI,JKP)
      PURW(JI,JKP)  = ( PUMF(JI,JK) * PURW(JI,JK) +                               &
                       PUER(JI,JKP) * PRW(JI,JK) - PUDR(JI,JKP) * PURW(JI,JK) )    &
                      / PUMF(JI,JKP)
!
!
      ZE1(JI) = ZE2(JI) ! update fractional entrainment/detrainment
      ZD1(JI) = ZD2(JI)
!
    END IF
  ENDDO
!
END DO
!
!*       12.1    Set OTRIG to False if cloud thickness < 0.5km
!                or > 3km (deep convection) or CAPE < 1
!                ------------------------------------------------
!
    DO JI = D%NIB, D%NIE
          JK  = KCTL(JI)
          ZWORK1(JI) = PZ(JI,JK) - PZLCL(JI)
          OTRIG(JI) = ZWORK1(JI) >= CVP_SHAL%XCDEPTH  .AND. ZWORK1(JI) < CVP_SHAL%XCDEPTH_D        &
                     .AND. PCAPE(JI) > 1.
    END DO
    DO JI = D%NIB, D%NIE
    IF( .NOT. OTRIG(JI) ) KCTL(JI) = IKB
    ENDDO
KETL(D%NIB:D%NIE) = MAX( KETL(D%NIB:D%NIE), KLCL(D%NIB:D%NIE) + 2 )
KETL(D%NIB:D%NIE) = MIN( KETL(D%NIB:D%NIE), KCTL(D%NIB:D%NIE) )
!
!
!*       12.2    If the ETL and CTL are the same detrain updraft mass
!                flux at this level
!                -------------------------------------------------------
!
ZWORK1(:) = 0.
DO JI=D%NIB, D%NIE
  IF ( KETL(JI) == KCTL(JI) ) ZWORK1(JI) = 1.
ENDDO
!
DO JI = D%NIB, D%NIE
    JK = KETL(JI)
    PUDR(JI,JK)   = PUDR(JI,JK) +                                    &
                          ( PUMF(JI,JK) - PUER(JI,JK) )  * ZWORK1(JI)
    PUER(JI,JK)   = PUER(JI,JK) * ( 1. - ZWORK1(JI) )
    PUMF(JI,JK)   = PUMF(JI,JK) * ( 1. - ZWORK1(JI) )
    JKP = KCTL(JI) + 1
    PUER(JI,JKP)  = 0. ! entrainm/detr rates have been already computed
    PUDR(JI,JKP)  = 0. ! at level KCTL+1, set them to zero
    PURW(JI,JKP)  = 0.
    PURC(JI,JKP)  = 0.
    PURI(JI,JKP)  = 0.
    PUTHL(JI,JKP) = 0.
    PURC(JI,JKP+1)= 0.
    PURI(JI,JKP+1)= 0.
END DO
!
!*       12.3    Adjust mass flux profiles, detrainment rates, and
!                precipitation fallout rates to reflect linear decrease
!                in mass flux between the ETL and CTL
!                -------------------------------------------------------
!
ZWORK1(:) = 0.
JK1 = IKB
JK2 = IKE

DO JK = JK1, JK2
    DO JI = D%NIB, D%NIE
    IF( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
        ZWORK1(JI) = ZWORK1(JI) + PDPRES(JI,JK)
    END IF
    END DO
END DO
!
DO JI = D%NIB, D%NIE
    JK = KETL(JI)
    ZWORK1(JI) = PUMF(JI,JK) / MAX( 1., ZWORK1(JI) )
END DO
!
DO JK = JK1 + 1, JK2
    JKP = JK - 1
    DO JI = D%NIB, D%NIE
    IF ( JK > KETL(JI) .AND. JK <= KCTL(JI) ) THEN
        PUDR(JI,JK)  = PDPRES(JI,JK) * ZWORK1(JI)
        PUMF(JI,JK)  = PUMF(JI,JKP) - PUDR(JI,JK)
    END IF
    END DO
END DO
!
!         12.4   Set mass flux and entrainment in the source layer.
!                Linear increase throughout the source layer.
!                -------------------------------------------------------
!
!IWORK(:) = MIN( KPBL(:), KLCL(:) - 1 )
IWORK(D%NIB:D%NIE) = KPBL(D%NIB:D%NIE)
DO JI = D%NIB, D%NIE
     JK  = KDPL(JI)
     JKP = IWORK(JI)
!          mixed layer depth
     ZWORK2(JI) = PPRES(JI,JK) - PPRES(JI,JKP) + PDPRES(JI,JK)
END DO
!
JKP=IKE
DO JK = JKM, JKP
   DO JI = D%NIB, D%NIE
   IF ( JK >= KDPL(JI)  .AND. JK <= IWORK(JI) .AND. GTRIG1(JI)) THEN
       PUER(JI,JK) = PUER(JI,JK) + PMFLCL * PDPRES(JI,JK) / ( ZWORK2(JI) + 0.1 )
       PUMF(JI,JK) = PUMF(JI,JK-1) + PUER(JI,JK)
   END IF
   END DO
END DO
!
!
!*       13.   If cloud thickness is smaller than  .5 km or > 3 km
!              no shallow convection is allowed
!              Nota: For technical reasons, we stop the convection
!                    computations in this case and do not go back to
!                    TRIGGER_FUNCT to look for the next unstable LCL
!                    which could produce a thicker cloud.
!              ---------------------------------------------------
!
DO JK=1,D%NKT
    DO JI=D%NIB, D%NIE
        IF(.NOT. OTRIG(JI))THEN
            PUMF(JI,JK)  = 0.
            PUDR(JI,JK)  = 0.
            PUER(JI,JK)  = 0.
            PUTHL(JI,JK) = PTHL(JI,JK)
            PURW(JI,JK)  = PRW(JI,JK)
            PURC(JI,JK)  = 0.
            PURI(JI,JK)  = 0.
        ENDIF
    ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('CONVECT_UPDRAFT_SHAL',1,ZHOOK_HANDLE)
END SUBROUTINE CONVECT_UPDRAFT_SHAL

