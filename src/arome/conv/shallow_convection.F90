!     ######spl
    SUBROUTINE SHALLOW_CONVECTION(CVP_SHAL, CST, D, NSV, CONVPAR, KBDIA, KTDIA, &
                                  KICE, OSETTADJ, PTADJS, PPABST, PZZ, &
                                  PTKECLS, PTT, PRVT, PRCT, PRIT, PWT, &
                                  PTTEN, PRVTEN, PRCTEN, PRITEN,       &
                                  KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
                                  PCH1, PCH1TEN)
    USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODI_SHALLOW_CONVECTION_PART1
USE MODI_SHALLOW_CONVECTION_PART2
USE MODI_SHALLOW_CONVECTION_PART2_SELECT
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
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
!
LOGICAL            ,INTENT(IN)     :: OCH1CONV ! include tracer transport
INTEGER            ,INTENT(IN)     :: KCH1     ! number of species
REAL               ,INTENT(IN)     :: PCH1(D%NIT,D%NKT,KCH1)! grid scale chemical species
REAL               ,INTENT(INOUT)  :: PCH1TEN(D%NIT,D%NKT,KCH1)! species conv. tendency (1/s)
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: ICONV
REAL     :: ZRDOCP                  ! R_d/C_p
!
REAL         :: ZTHT(D%NIT,D%NKT), ZSTHV(D%NIT,D%NKT), ZSTHES(D%NIT,D%NKT)  ! grid scale theta, theta_v
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER  :: ISDPL(D%NIT)   ! index for parcel departure level
INTEGER  :: ISPBL(D%NIT)   ! index for source layer top
INTEGER  :: ISLCL(D%NIT)   ! index for lifting condensation level
REAL     :: ZSTHLCL(D%NIT) ! updraft theta at LCL/L
REAL     :: ZSTLCL(D%NIT)  ! updraft temp. at LCL
REAL     :: ZSRVLCL(D%NIT) ! updraft rv at LCL
REAL     :: ZSWLCL(D%NIT)  ! updraft w at LCL
REAL     :: ZSZLCL(D%NIT)  ! LCL height
REAL     :: ZSTHVELCL(D%NIT)! envir. theta_v at LCL
!
LOGICAL    :: GTRIG1(D%NIT)  ! logical mask for convection
!
TYPE(CONVPAREXT) :: CVPEXT
!
!-------------------------------------------------------------------------------
!
!
!*       0.3    Compute loop bounds
!               -------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION',0,ZHOOK_HANDLE)

CVPEXT%JCVEXB = MAX( 0, KBDIA - 1 )
CVPEXT%JCVEXT = MAX( 0, KTDIA - 1)

ZRDOCP = CST%XRD / CST%XCPD

CALL SHALLOW_CONVECTION_PART1&
   (CVPEXT, CVP_SHAL, CST, D, NSV, CONVPAR, KBDIA, KTDIA, &
    KICE, OSETTADJ, PTADJS, PPABST, PZZ, &
    PTKECLS, PTT, PRVT, PRCT, PRIT, PWT, &
    PTTEN, PRVTEN, PRCTEN, PRITEN,       &
    KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
    PCH1, PCH1TEN, ZTHT, ZSTHV, ZSTHES,  &
    ISDPL, ISPBL, ISLCL, ZSTHLCL, ZSTLCL,&
    ZSRVLCL, ZSWLCL, ZSZLCL, ZSTHVELCL, GTRIG1)

ICONV = COUNT(GTRIG1(D%NIB:D%NIE))

IF(ICONV==0)THEN
  ! Do nothing if there are no selected columns
ELSE IF (ICONV < D%NIT*9/10) THEN
  CALL SHALLOW_CONVECTION_PART2_SELECT &
                             & (CVP_SHAL, CVPEXT, CST, D, NSV, CONVPAR, KICE, &
                                OSETTADJ, PTADJS, PPABST, PZZ, PTT,  &
                                PRVT, PRCT, PRIT, OCH1CONV, KCH1,    &
                                PCH1, ZRDOCP, ZTHT, ZSTHV, ZSTHES,   &
                                ISDPL, ISPBL, ISLCL, ZSTHLCL, ZSTLCL,&
                                ZSRVLCL, ZSWLCL, ZSZLCL, ZSTHVELCL,  &
                                GTRIG1, PUMF, PTTEN, PRVTEN, PRCTEN, &
                                PRITEN, KCLTOP, KCLBAS, PCH1TEN, ICONV)
ELSE
  CALL SHALLOW_CONVECTION_PART2 &
                             & (CVP_SHAL, CVPEXT, CST, D, NSV, CONVPAR, KICE, &
                                OSETTADJ, PTADJS, PPABST, PZZ, PTT,  &
                                PRVT, PRCT, PRIT, OCH1CONV, KCH1,    &
                                PCH1, ZRDOCP, ZTHT, ZSTHV, ZSTHES,   &
                                ISDPL, ISPBL, ISLCL, ZSTHLCL, ZSTLCL,&
                                ZSRVLCL, ZSWLCL, ZSZLCL, ZSTHVELCL,  &
                                GTRIG1, PUMF, PTTEN, PRVTEN, PRCTEN, &
                                PRITEN, KCLTOP, KCLBAS, PCH1TEN)
ENDIF

IF (LHOOK) CALL DR_HOOK('SHALLOW_CONVECTION',1,ZHOOK_HANDLE)
END SUBROUTINE SHALLOW_CONVECTION

