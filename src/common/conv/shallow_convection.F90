!MNH_LIC Copyright 1996-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!   ###############################################################################
    SUBROUTINE SHALLOW_CONVECTION(CVP_SHAL, CST, D, NSV, CONVPAR, KBDIA, KTDIA, &
                                  KICE, OSETTADJ, PTADJS, PPABST, PZZ, &
                                  PTKECLS, PTT, PRVT, PRCT, PRIT, PWT, &
                                  PTTEN, PRVTEN, PRCTEN, PRITEN,       &
                                  KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
                                  PCH1, PCH1TEN)
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
!  P. Wautelet 03/06/2019: simplify code (remove always true masks) + replace PACK intrinsics
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_CST, ONLY: CST_t
USE MODD_CONVPAREXT, ONLY : CONVPAREXT
USE MODD_CONVPAR, ONLY: CONVPAR_T
USE MODD_CONVPAR_SHAL, ONLY: CONVPAR_SHAL
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_NSV, ONLY: NSV_T
USE MODI_SHALLOW_CONVECTION_PART1, ONLY: SHALLOW_CONVECTION_PART1
USE MODI_SHALLOW_CONVECTION_PART2, ONLY: SHALLOW_CONVECTION_PART2
USE MODI_SHALLOW_CONVECTION_PART2_SELECT, ONLY: SHALLOW_CONVECTION_PART2_SELECT
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
INTEGER,                    INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                  ! KBDIA that is at least 1
INTEGER,                    INTENT(IN) :: KTDIA    ! vertical computations can be
                                                   ! limited to D%NKT + 1 - KTDIA
                                                   ! default=1
INTEGER,                    INTENT(IN) :: KICE     ! flag for ice ( 1 = yes,
                                                   !                0 = no ice )
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
                                                   ! adjustment time by user
REAL,                       INTENT(IN) :: PTADJS   ! user defined adjustment time
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PZZ      ! height of model layer (m)
REAL, DIMENSION(D%NIJT),       INTENT(IN) :: PTKECLS  ! TKE in the CLS  (m2/s2)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PWT      ! grid scale vertical
                                                      ! velocity (m/s)
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT):: PTTEN  ! convective temperature tendency (K/s)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
INTEGER, DIMENSION(D%NIJT),    INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(D%NIJT),    INTENT(INOUT):: KCLBAS ! cloud base level
                                                      ! they are given a value of
                                                      ! 0 if no convection
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(D%NIJT,D%NKT,KCH1), INTENT(IN)   :: PCH1   ! grid scale chemical species
REAL, DIMENSION(D%NIJT,D%NKT,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: ICONV                   ! number of convective columns
REAL     :: ZRDOCP                  ! R_d/C_p
!
REAL, DIMENSION(D%NIJT,D%NKT)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, theta_v
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER  :: ISDPL(D%NIJT)   ! index for parcel departure level
INTEGER  :: ISPBL(D%NIJT)   ! index for source layer top
INTEGER  :: ISLCL(D%NIJT)   ! index for lifting condensation level
REAL     :: ZSTHLCL(D%NIJT) ! updraft theta at LCL/L
REAL     :: ZSTLCL(D%NIJT)  ! updraft temp. at LCL
REAL     :: ZSRVLCL(D%NIJT) ! updraft rv at LCL
REAL     :: ZSWLCL(D%NIJT)  ! updraft w at LCL
REAL     :: ZSZLCL(D%NIJT)  ! LCL height
REAL     :: ZSTHVELCL(D%NIJT)! envir. theta_v at LCL
!
LOGICAL    :: GTRIG1(D%NIJT)  ! logical mask for convection
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

CALL SHALLOW_CONVECTION_PART1(CVPEXT, CVP_SHAL, CST, D, NSV, CONVPAR, KBDIA, KTDIA, &
                              KICE, OSETTADJ, PTADJS, PPABST, PZZ, &
                              PTKECLS, PTT, PRVT, PRCT, PRIT, PWT, &
                              PTTEN, PRVTEN, PRCTEN, PRITEN,       &
                              KCLTOP, KCLBAS, PUMF, OCH1CONV, KCH1,&
                              PCH1, PCH1TEN, ZTHT, ZSTHV, ZSTHES,  &
                              ISDPL, ISPBL, ISLCL, ZSTHLCL, ZSTLCL,&
                              ZSRVLCL, ZSWLCL, ZSZLCL, ZSTHVELCL, GTRIG1)

ICONV = COUNT(GTRIG1(D%NIJB:D%NIJE))

IF(ICONV==0)THEN
  ! Do nothing if there are no selected columns
ELSE IF (ICONV < D%NIJT*9/10) THEN
  CALL SHALLOW_CONVECTION_PART2_SELECT(CVP_SHAL, CVPEXT, CST, D, NSV, CONVPAR, KICE, &
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
