!     ######spl
    SUBROUTINE SHALLOW_CONVECTION_PART1&
                                 (CVPEXT, CVP_SHAL, CST, D, NSV, CONVPAR, KBDIA, KTDIA, &
                                  KICE, OSETTADJ, PTADJS, PPABST, PZZ, &
                                  PTKECLS, PTT, PRVT, PRCT, PRIT, PWT, &
                                  PTTEN, PRVTEN, PRCTEN, PRITEN,       &
                                  KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
                                  PCH1, PCH1TEN, PTHT, PSTHV, PSTHES,  &
                                  KSDPL, KSPBL, KSLCL, PSTHLCL, PSTLCL,&
                                  PSRVLCL, PSWLCL, PSZLCL, PSTHVELCL, OTRIG1)
    USE PARKIND1, ONLY : JPRB
    USE YOMHOOK , ONLY : LHOOK, JPHOOK, DR_HOOK
!   ###############################################################################
!
!!**** Monitor routine to compute all convective tendencies by calls
!!     of several subroutines.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to determine the convective
!!      tendencies. The routine first prepares all necessary grid-scale
!!      variables. The final convective tendencies are then computed by
!!      calls of different subroutines.
!!
!!
!!**  METHOD
!!    ------
!!      We start by selecting convective columns in the model domain through
!!      the call of routine TRIGGER_FUNCT. Then, we allocate memory for the
!!      convection updraft and downdraft variables and gather the grid scale
!!      variables in convective arrays.
!!      The updraft and downdraft computations are done level by level starting
!!      at the  bottom and top of the domain, respectively.
!!      All computations are done on MNH thermodynamic levels. The depth
!!      of the current model layer k is defined by DP(k)=P(k-1)-P(k)
!!
!!
!!
!!    EXTERNAL
!!    --------
!!    CONVECT_TRIGGER_SHAL
!!    CONVECT_SATMIXRATIO
!!    CONVECT_UPDRAFT_SHAL
!!        CONVECT_CONDENS
!!        CONVECT_MIXING_FUNCT
!!    CONVECT_CLOSURE_SHAL
!!        CONVECT_CLOSURE_THRVLCL
!!        CONVECT_CLOSURE_ADJUST_SHAL
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!          XG                   ! gravity constant
!!          XPI                  ! number Pi
!!          XP00                 ! reference pressure
!!          XRD, XRV             ! gaz  constants for dry air and water vapor
!!          XCPD, XCPV           ! specific heat for dry air and water vapor
!!          XRHOLW               ! density of liquid water
!!          XALPW, XBETAW, XGAMW ! constants for water saturation pressure
!!          XTT                  ! triple point temperature
!!          XLVTT, XLSTT         ! vaporization, sublimation heat constant
!!          XCL, XCI             ! specific heat for liquid water and ice
!!
!!      Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT       ! extra levels on the vertical boundaries
!!
!!      Module MODD_CONVPAR
!!          XA25                 ! reference grid area
!!          XCRAD                ! cloud radius
!!
!!
!!    REFERENCE
!!    ---------
!!
!!      Bechtold, 1997 : Meso-NH scientific  documentation (31 pp)
!!      Fritsch and Chappell, 1980, J. Atmos. Sci., Vol. 37, 1722-1761.
!!      Kain and Fritsch, 1990, J. Atmos. Sci., Vol. 47, 2784-2801.
!!      Kain and Fritsch, 1993, Meteor. Monographs, Vol. 24, 165-170.
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96
!!   Peter Bechtold 15/11/96 replace theta_il by enthalpy
!!         "        10/12/98 changes for ARPEGE
!!         "        01/01/02 Apply conservation correction
!!   F Bouyssel     05/11/08 Modifications for reproductibility
!!   E. Bazile      20/07/09 Input of TKECLS.
!!   F. Bouyssel    08/11/13 Modifications for reproductibility
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST, ONLY : CST_T
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_CONVPAR, ONLY: CONVPAR_T
USE MODD_CONVPAR_SHAL, ONLY: CONVPAR_SHAL
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_NSV, ONLY: NSV_T
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
TYPE(CONVPAREXT)   ,INTENT(IN)     :: CVPEXT
TYPE(CONVPAR_SHAL) ,INTENT(IN)     :: CVP_SHAL
TYPE(CST_T)        ,INTENT(IN)     :: CST
TYPE(DIMPHYEX_T)   ,INTENT(IN)     :: D
TYPE(NSV_T)        ,INTENT(IN)     :: NSV
TYPE(CONVPAR_T)    ,INTENT(IN)     :: CONVPAR
INTEGER            ,INTENT(IN)     :: KBDIA    ! vertical  computations start at
!                                                    ! KBDIA that is at least 1
INTEGER            ,INTENT(IN)     :: KTDIA    ! vertical computations can be
                                                     ! limited to D%NKT + 1 - KTDIA
                                                     ! default=1
                                                     ! scheme
INTEGER            ,INTENT(IN)     :: KICE     ! flag for ice ( 1 = yes,
                                                     !                0 = no ice )
LOGICAL            ,INTENT(IN)     :: OSETTADJ ! logical to set convective
                                                     ! adjustment time by user
REAL               ,INTENT(IN)     :: PTADJS   ! user defined adjustment time
REAL               ,INTENT(IN)     :: PPABST(D%NIT,D%NKT)   ! grid scale pressure at t
REAL               ,INTENT(IN)     :: PZZ(D%NIT,D%NKT)      ! height of model layer (m)
REAL               ,INTENT(IN)     :: PTKECLS(D%NIT)  ! TKE in the CLS  (m2/s2)
REAL               ,INTENT(IN)     :: PTT(D%NIT,D%NKT)      ! grid scale temperature at t
REAL               ,INTENT(IN)     :: PRVT(D%NIT,D%NKT)     ! grid scale water vapor "
REAL               ,INTENT(IN)     :: PRCT(D%NIT,D%NKT)     ! grid scale r_c  "
REAL               ,INTENT(IN)     :: PRIT(D%NIT,D%NKT)     ! grid scale r_i "
REAL               ,INTENT(IN)     :: PWT(D%NIT,D%NKT)      ! grid scale vertical
                                                                           ! velocity (m/s)
!
REAL               ,INTENT(INOUT)  :: PTTEN(D%NIT,D%NKT)  ! convective temperature
                                                                           ! tendency (K/s)
REAL               ,INTENT(INOUT)  :: PRVTEN(D%NIT,D%NKT) ! convective r_v tendency (1/s)
REAL               ,INTENT(INOUT)  :: PRCTEN(D%NIT,D%NKT) ! convective r_c tendency (1/s)
REAL               ,INTENT(INOUT)  :: PRITEN(D%NIT,D%NKT) ! convective r_i tendency (1/s)
INTEGER            ,INTENT(INOUT)  :: KCLTOP(D%NIT) ! cloud top level
INTEGER            ,INTENT(INOUT)  :: KCLBAS(D%NIT) ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL               ,INTENT(INOUT)  :: PUMF(D%NIT,D%NKT)   ! updraft mass flux (kg/s m2)
LOGICAL            ,INTENT(IN)     :: OCH1CONV ! include tracer transport
INTEGER            ,INTENT(IN)     :: KCH1     ! number of species
REAL               ,INTENT(IN)     :: PCH1(D%NIT,D%NKT,KCH1)! grid scale chemical species
REAL               ,INTENT(INOUT)  :: PCH1TEN(D%NIT,D%NKT,KCH1)! species conv. tendency (1/s)
REAL               ,INTENT(OUT)    :: PTHT(D%NIT,D%NKT) 
REAL               ,INTENT(OUT)    :: PSTHV(D%NIT,D%NKT) 
REAL               ,INTENT(OUT)    :: PSTHES(D%NIT,D%NKT)  ! grid scale theta, theta_v
INTEGER            ,INTENT(OUT)    :: KSDPL(D%NIT)   ! index for parcel departure level
INTEGER            ,INTENT(OUT)    :: KSPBL(D%NIT)   ! index for source layer top
INTEGER            ,INTENT(OUT)    :: KSLCL(D%NIT)   ! index for lifting condensation level
REAL               ,INTENT(OUT)    :: PSTHLCL(D%NIT) ! updraft theta at LCL/L
REAL               ,INTENT(OUT)    :: PSTLCL(D%NIT)  ! updraft temp. at LCL
REAL               ,INTENT(OUT)    :: PSRVLCL(D%NIT) ! updraft rv at LCL
REAL               ,INTENT(OUT)    :: PSWLCL(D%NIT)  ! updraft w at LCL
REAL               ,INTENT(OUT)    :: PSZLCL(D%NIT)  ! LCL height
REAL               ,INTENT(OUT)    :: PSTHVELCL(D%NIT)! envir. theta_v at LCL
LOGICAL            ,INTENT(OUT)    :: OTRIG1(D%NIT)  ! logical mask for convection
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: IKB, IKE                ! vertical loop bounds
INTEGER  :: JI                      ! horizontal loop index
INTEGER  :: JK                      ! vertical loop index
INTEGER  :: ICONV
REAL     :: ZEPS, ZEPSA             ! R_d / R_v, R_v / R_d
REAL     :: ZRDOCP                  ! R_d/C_p
!
REAL                        :: ZES     ! saturation vapor mixng ratio
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

#include "convect_trigger_shal.h"

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_PART1',0,ZHOOK_HANDLE)

IKB = 1 + CVPEXT%JCVEXB
IKE = D%NKT - CVPEXT%JCVEXT
!
!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------
!
PTTEN(:,:)  = 0.
PRVTEN(:,:) = 0.
PRCTEN(:,:) = 0.
PRITEN(:,:) = 0.
PUMF(:,:)   = 0.
KCLTOP(:)  = 0
KCLBAS(:)  = 0

IF ( OCH1CONV ) THEN
  PCH1TEN(:,:,:) = 0.
END IF
!
!
!*       1.     Initialize  local variables
!               ----------------------------
!
ZEPS   = CST%XRD / CST%XRV
ZEPSA  = CST%XRV / CST%XRD
ZRDOCP = CST%XRD / CST%XCPD
!
!-------------------------------------------------------------------------------
!
!*       1.1    Set up grid scale theta, theta_v, theta_es
!               ------------------------------------------
!
PTHT(:,:) = 300.
PSTHV(:,:)= 300.
PSTHES(:,:)= 400.
DO JK = IKB, IKE
DO JI = D%NIB, D%NIE
  IF ( PPABST(JI,JK) > 40.E2 ) THEN
    PTHT(JI,JK)  = PTT(JI,JK) * ( CST%XP00 / PPABST(JI,JK) ) ** ZRDOCP
    PSTHV(JI,JK) = PTHT(JI,JK) * ( 1. + ZEPSA * PRVT(JI,JK) ) /              &
                   ( 1. + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )
!
        ! use conservative Bolton (1980) formula for theta_e
        ! it is used to compute CAPE for undilute parcel ascent
        ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here
!
    ZES = EXP( CST%XALPW - CST%XBETAW / PTT(JI,JK) - CST%XGAMW * LOG( PTT(JI,JK) ) )
    ZES = MIN( 1., ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
    PSTHES(JI,JK) = PTT(JI,JK) * ( PTHT(JI,JK) / PTT(JI,JK) ) **             &
              ( 1. - 0.28 * ZES ) * EXP( ( 3374.6525 / PTT(JI,JK) - 2.5403 ) &
                                        * ZES * ( 1. + 0.81 * ZES ) )
  END IF
END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*       2.     Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------
!
!*       2.3    Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------
!
KSLCL(:) = MAX( IKB, 2 )   ! initialize DPL PBL and LCL
KSDPL(:) = IKB
KSPBL(:) = IKB
!
CALL CONVECT_TRIGGER_SHAL(CVP_SHAL, CVPEXT, CST, D, PPABST, PTHT,      &
                          PSTHV, PSTHES, PRVT, PWT, PZZ, PTKECLS,      &
                          PSTHLCL, PSTLCL, PSRVLCL, PSWLCL, PSZLCL,    &
                          PSTHVELCL, KSLCL, KSDPL, KSPBL, OTRIG1)

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION_PART1',1,ZHOOK_HANDLE)

END SUBROUTINE


