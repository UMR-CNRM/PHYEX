!MNH_LIC Copyright 1996-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_SHALLOW_CONVECTION
!     ######################
!
INTERFACE
!
    SUBROUTINE SHALLOW_CONVECTION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                                   PDTCONV, KICE, OSETTADJ, PTADJS,               &
                                   PPABST, PZZ, PTKECLS,                          &
                                   PTT, PRVT, PRCT, PRIT, PWT,                    &
                                   PTTEN, PRVTEN, PRCTEN, PRITEN,                 &
                                   KCLTOP, KCLBAS, PUMF,                          &
                                   OCH1CONV, KCH1, PCH1, PCH1TEN                  )
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                    INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                    INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                    INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                  ! KBDIA that is at least 1
INTEGER,                    INTENT(IN) :: KTDIA    ! vertical computations can be
                                                   ! limited to KLEV + 1 - KTDIA
                                                   ! default=1
REAL,                       INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                   ! calls of the deep convection
                                                   ! scheme
INTEGER,                    INTENT(IN) :: KICE     ! flag for ice ( 1 = yes, 
                                                   !                0 = no ice )
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
                                                   ! adjustment time by user
REAL,                       INTENT(IN) :: PTADJS   ! user defined adjustment time 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical 
                                                   ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PTKECLS  ! TKE in the CLS  (m2/s2) 
!   
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature 
                                                   ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
!
END SUBROUTINE SHALLOW_CONVECTION
!
END INTERFACE
!
END MODULE MODI_SHALLOW_CONVECTION
!   ###############################################################################
    SUBROUTINE SHALLOW_CONVECTION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                                   PDTCONV, KICE, OSETTADJ, PTADJS,               &
                                   PPABST, PZZ, PTKECLS,                          &
                                   PTT, PRVT, PRCT, PRIT, PWT,                    &
                                   PTTEN, PRVTEN, PRCTEN, PRITEN,                 &
                                   KCLTOP, KCLBAS, PUMF,                          &
                                   OCH1CONV, KCH1, PCH1, PCH1TEN                  )
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
!  P. Wautelet 03/06/2019: simplify code (remove always true masks) + replace PACK intrinsics
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAREXT
USE MODD_CONVPAR_SHAL
USE MODD_NSV,       ONLY : NSV_LGBEG,NSV_LGEND
!
USE MODI_CONVECT_TRIGGER_SHAL
USE MODI_CONVECT_UPDRAFT_SHAL
USE MODI_CONVECT_CLOSURE_SHAL
USE MODI_CONVECT_CHEM_TRANSPORT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                    INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                    INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                    INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                  ! KBDIA that is at least 1
INTEGER,                    INTENT(IN) :: KTDIA    ! vertical computations can be
                                                   ! limited to KLEV + 1 - KTDIA
                                                   ! default=1
REAL,                       INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                   ! calls of the deep convection
                                                   ! scheme
INTEGER,                    INTENT(IN) :: KICE     ! flag for ice ( 1 = yes, 
                                                   !                0 = no ice )
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
                                                   ! adjustment time by user
REAL,                       INTENT(IN) :: PTADJS   ! user defined adjustment time 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical 
                                                   ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PTKECLS  ! TKE in the CLS  (m2/s2) 
!   
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature 
                                                   ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: ICONV                   ! number of convective columns
INTEGER  :: IIB, IIE                ! horizontal loop bounds
INTEGER  :: IKB, IKE                ! vertical loop bounds
INTEGER  :: IKS                     ! vertical dimension
INTEGER  :: JI, JL                  ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKM, JKP            ! vertical loop index
INTEGER  :: IFTSTEPS                ! only used for chemical tracers
REAL     :: ZEPS, ZEPSA             ! R_d / R_v, R_v / R_d
REAL     :: ZRDOCP                  ! R_d/C_p
!
REAL, DIMENSION(KLON,KLEV)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, theta_v
REAL, DIMENSION(KLON)              :: ZWORK2, ZWORK2B ! work array 
REAL                               :: ZW1     ! work variable
!
!
!*       0.2   Declarations of local allocatable  variables :
!
INTEGER, DIMENSION(:),ALLOCATABLE  :: IDPL    ! index for parcel departure level
INTEGER, DIMENSION(:),ALLOCATABLE  :: IPBL    ! index for source layer top
INTEGER, DIMENSION(:),ALLOCATABLE  :: ILCL    ! index for lifting condensation level 
INTEGER, DIMENSION(:),ALLOCATABLE  :: IETL    ! index for zero buoyancy level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ICTL    ! index for cloud top level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ILFS    ! index for level of free sink
!
INTEGER, DIMENSION(:), ALLOCATABLE :: ISDPL   ! index for parcel departure level
INTEGER, DIMENSION(:),ALLOCATABLE  :: ISPBL   ! index for source layer top
INTEGER, DIMENSION(:), ALLOCATABLE :: ISLCL   ! index for lifting condensation level 
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTHLCL ! updraft theta at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTLCL  ! updraft temp. at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSRVLCL ! updraft rv at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSWLCL  ! updraft w at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSZLCL  ! LCL height
REAL, DIMENSION(:), ALLOCATABLE    :: ZSTHVELCL! envir. theta_v at LCL
REAL, DIMENSION(:), ALLOCATABLE    :: ZSDXDY  ! grid area (m^2)
!
! grid scale variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZZ      ! height of model layer (m) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPRES   ! grid scale pressure
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDPRES  ! pressure difference between 
                                              ! bottom and top of layer (Pa)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZW      ! grid scale vertical velocity on theta grid
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTT     ! temperature
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTH     ! grid scale theta     
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHV    ! grid scale theta_v     
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHL    ! grid scale enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHES, ZTHEST ! grid scale saturated theta_e
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRW     ! grid scale total water (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRV     ! grid scale water vapor (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRC     ! grid scale cloud water (kg/kg) 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRI     ! grid scale cloud ice (kg/kg) 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZDXDY   ! grid area (m^2)
!
! updraft variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUMF    ! updraft mass flux (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUER    ! updraft entrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUDR    ! updraft detrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTHL   ! updraft enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTHV   ! updraft theta_v (K)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURW    ! updraft total water (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURC    ! updraft cloud water (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURI    ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZMFLCL  ! cloud base unit mass flux(kg/s) 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZCAPE   ! available potent. energy     
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTHLCL  ! updraft theta at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTLCL   ! updraft temp. at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZRVLCL  ! updraft rv at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZWLCL   ! updraft w at LCL
REAL, DIMENSION(:),   ALLOCATABLE  :: ZZLCL   ! LCL height
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTHVELCL! envir. theta_v at LCL
!
! downdraft variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDMF    ! downdraft mass flux (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDER    ! downdraft entrainment (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDDR    ! downdraft detrainment (kg/s)
!
! closure variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZLMASS  ! mass of model layer (kg)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTIMEC  ! advective time period
!
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRVC    ! conv. adj. grid scale r_w 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRCC    ! conv. adj. grid scale r_c 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRIC    ! conv. adj. grid scale r_i 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZWSUB   ! envir. compensating subsidence (Pa/s)
!
LOGICAL, DIMENSION(:),ALLOCATABLE  :: GTRIG1  ! logical mask for convection
INTEGER, DIMENSION(:),ALLOCATABLE  :: IJINDEX ! hor.index
REAL, DIMENSION(:),   ALLOCATABLE  :: ZCPH    ! specific heat C_ph 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZLV, ZLS! latent heat of vaporis., sublim.
REAL                               :: ZES     ! saturation vapor mixng ratio
!
! Chemical Tracers:
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1    ! grid scale chemical specy (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1C   ! conv. adjust. chemical specy 1
REAL, DIMENSION(:,:),   ALLOCATABLE:: ZWORK3  ! conv. adjust. chemical specy 1
!
!-------------------------------------------------------------------------------
!
!
!*       0.3    Compute loop bounds
!               -------------------
!
IIB = KIDIA
IIE = KFDIA
JCVEXB = MAX( 0, KBDIA - 1 )
IKB = 1 + JCVEXB 
IKS = KLEV
JCVEXT = MAX( 0, KTDIA - 1)
IKE = IKS - JCVEXT 
!
!
!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------
!
PTTEN(:,:)  = 0.
PRVTEN(:,:) = 0.
PRCTEN(:,:) = 0.
PRITEN(:,:) = 0.
! PUTEN(:,:)  = 0.
! PVTEN(:,:)  = 0.
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
ZEPS   = XRD / XRV
ZEPSA  = XRV / XRD 
ZRDOCP = XRD / XCPD
!
!-------------------------------------------------------------------------------
!
!*       1.1    Set up grid scale theta, theta_v, theta_es 
!               ------------------------------------------
!
ZTHT(:,:) = 300.
ZSTHV(:,:)= 300.
ZSTHES(:,:)= 400.
DO JK = IKB, IKE
DO JI = IIB, IIE
  IF ( PPABST(JI,JK) > 40.E2 ) THEN
    ZTHT(JI,JK)  = PTT(JI,JK) * ( XP00 / PPABST(JI,JK) ) ** ZRDOCP
    ZSTHV(JI,JK) = ZTHT(JI,JK) * ( 1. + ZEPSA * PRVT(JI,JK) ) /              &
                   ( 1. + PRVT(JI,JK) + PRCT(JI,JK) + PRIT(JI,JK) )
!
        ! use conservative Bolton (1980) formula for theta_e
        ! it is used to compute CAPE for undilute parcel ascent
        ! For economical reasons we do not use routine CONVECT_SATMIXRATIO here
!
    ZES = EXP( XALPW - XBETAW / PTT(JI,JK) - XGAMW * LOG( PTT(JI,JK) ) )
    ZES = MIN( 1., ZEPS * ZES / ( PPABST(JI,JK) - ZES ) )
    ZSTHES(JI,JK) = PTT(JI,JK) * ( ZTHT(JI,JK) / PTT(JI,JK) ) **             &
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
!*       2.1    Allocate arrays depending on number of model columns that need
!               to be tested for convection (i.e. where no convection is present
!               at the moment.
!               --------------------------------------------------------------
!
ALLOCATE( ZPRES(KLON,IKS) )
ALLOCATE( ZZ(KLON,IKS) )
ALLOCATE( ZW(KLON,IKS) )
ALLOCATE( ZTH(KLON,IKS) )
ALLOCATE( ZTHV(KLON,IKS) )
ALLOCATE( ZTHEST(KLON,IKS) )
ALLOCATE( ZRV(KLON,IKS) )
ALLOCATE( ZSTHLCL(KLON) )
ALLOCATE( ZSTLCL(KLON) )
ALLOCATE( ZSRVLCL(KLON) )
ALLOCATE( ZSWLCL(KLON) )
ALLOCATE( ZSZLCL(KLON) )
ALLOCATE( ZSTHVELCL(KLON) )
ALLOCATE( ISDPL(KLON) )
ALLOCATE( ISPBL(KLON) )
ALLOCATE( ISLCL(KLON) )
ALLOCATE( ZSDXDY(KLON) )
ALLOCATE( GTRIG1(KLON) )
!
DO JK = IKB, IKE
DO JI = 1, KLON
  JL = JI
  ZPRES(JI,JK)  = PPABST(JI,JK)
  ZZ(JI,JK)     = PZZ(JI,JK)
  ZTH(JI,JK)    = ZTHT(JI,JK)
  ZTHV(JI,JK)   = ZSTHV(JI,JK)
  ZTHEST(JI,JK) = ZSTHES(JI,JK)
  ZRV(JI,JK)    = MAX( 0., PRVT(JI,JK) )
  ZW(JI,JK)     = PWT(JI,JK)
END DO
END DO
ZSDXDY(:)    = XA25
!
!*       2.2    Compute environm. enthalpy and total water = r_v + r_i + r_c
!               and envir. saturation theta_e
!               ------------------------------------------------------------
!
!
!*       2.3    Test for convective columns and determine properties at the LCL
!               --------------------------------------------------------------
!
ISLCL(:) = MAX( IKB, 2 )   ! initialize DPL PBL and LCL 
ISDPL(:) = IKB
ISPBL(:) = IKB
!
CALL CONVECT_TRIGGER_SHAL(  KLON, KLEV,                              &
                            ZPRES, ZTH, ZTHV, ZTHEST,                 &
                            ZRV, ZW, ZZ, ZSDXDY, PTKECLS,             &
                            ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
                            ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1    )
!
DEALLOCATE( ZPRES )
DEALLOCATE( ZZ )
DEALLOCATE( ZTH )
DEALLOCATE( ZTHV )
DEALLOCATE( ZTHEST )
DEALLOCATE( ZRV )
DEALLOCATE( ZW )
!
!-------------------------------------------------------------------------------
!
!*       3.     After the call of TRIGGER_FUNCT we allocate all the dynamic
!               arrays used in the convection scheme using the mask GTRIG, i.e.
!               we do calculus only in convective columns. This corresponds to
!               a GATHER operation.
!               --------------------------------------------------------------
!
ICONV = COUNT( GTRIG1(:) )
IF ( ICONV == 0 )  THEN 
  DEALLOCATE( ZSTHLCL )
  DEALLOCATE( ZSTLCL )
  DEALLOCATE( ZSRVLCL )
  DEALLOCATE( ZSWLCL )
  DEALLOCATE( ZSZLCL )
  DEALLOCATE( ZSTHVELCL )
  DEALLOCATE( ZSDXDY )
  DEALLOCATE( ISLCL )
  DEALLOCATE( ISDPL )
  DEALLOCATE( ISPBL )
  DEALLOCATE( GTRIG1 )
  RETURN   ! no convective column has been found, exit DEEP_CONVECTION
ENDIF
!
     ! vertical index variables
!
ALLOCATE( IDPL(ICONV) )
ALLOCATE( IPBL(ICONV) )
ALLOCATE( ILCL(ICONV) )
ALLOCATE( ICTL(ICONV) )
ALLOCATE( IETL(ICONV) )
!
	 ! grid scale variables
!
ALLOCATE( ZZ(ICONV,IKS) ) ;   ZZ  = 0.0
ALLOCATE( ZPRES(ICONV,IKS) );   ZPRES = 0.0
ALLOCATE( ZDPRES(ICONV,IKS) ) ;   ZDPRES = 0.0
ALLOCATE( ZTT(ICONV, IKS) ) ;   ZTT    = 0.0
ALLOCATE( ZTH(ICONV,IKS) ) ;   ZTH    = 0.0
ALLOCATE( ZTHV(ICONV,IKS) )   ;   ZTHV   = 0.0
ALLOCATE( ZTHL(ICONV,IKS) )  ; ZTHL  = 0.0
ALLOCATE( ZTHES(ICONV,IKS) )  ; ZTHES = 0.0
ALLOCATE( ZRV(ICONV,IKS) ) ; ZRV   = 0.0
ALLOCATE( ZRC(ICONV,IKS) )   ; ZRC   = 0.0
ALLOCATE( ZRI(ICONV,IKS) )   ; ZRI   = 0.0
ALLOCATE( ZRW(ICONV,IKS) )   ; ZRW   = 0.0
ALLOCATE( ZDXDY(ICONV) )  ; ZDXDY = 0.0
!
         ! updraft variables
!
ALLOCATE( ZUMF(ICONV,IKS) )
ALLOCATE( ZUER(ICONV,IKS) )
ALLOCATE( ZUDR(ICONV,IKS) )
ALLOCATE( ZUTHL(ICONV,IKS) )
ALLOCATE( ZUTHV(ICONV,IKS) )
ALLOCATE( ZURW(ICONV,IKS) )
ALLOCATE( ZURC(ICONV,IKS) )
ALLOCATE( ZURI(ICONV,IKS) )
ALLOCATE( ZTHLCL(ICONV) )
ALLOCATE( ZTLCL(ICONV) )
ALLOCATE( ZRVLCL(ICONV) )
ALLOCATE( ZWLCL(ICONV) )
ALLOCATE( ZMFLCL(ICONV) )
ALLOCATE( ZZLCL(ICONV) )
ALLOCATE( ZTHVELCL(ICONV) )
ALLOCATE( ZCAPE(ICONV) )
!
         ! work variables
!
ALLOCATE( IJINDEX(ICONV) )
ALLOCATE( ZCPH(ICONV) )
ALLOCATE( ZLV(ICONV) )
ALLOCATE( ZLS(ICONV) )
!
!
!*           3.1    Gather grid scale and updraft base variables in
!                   arrays using mask GTRIG
!                   ---------------------------------------------------
!
JL = 1
DO JI = 1, KLON
  IF ( GTRIG1(JI) ) THEN
    IJINDEX(JL) = JI
    JL = JL +1
  END IF
END DO
!
DO JK = IKB, IKE
DO JI = 1, ICONV
  JL = IJINDEX(JI)
  ZZ(JI,JK)     = PZZ(JL,JK)
  ZPRES(JI,JK)  = PPABST(JL,JK)
  ZTT(JI,JK)    = PTT(JL,JK)
  ZTH(JI,JK)    = ZTHT(JL,JK)
  ZTHES(JI,JK)  = ZSTHES(JL,JK)
  ZRV(JI,JK)    = MAX( 0., PRVT(JL,JK) )
  ZRC(JI,JK)    = MAX( 0., PRCT(JL,JK) )
  ZRI(JI,JK)    = MAX( 0., PRIT(JL,JK) )
  ZTHV(JI,JK)   = ZSTHV(JL,JK)
END DO
END DO
!
DO JI = 1, ICONV
  JL = IJINDEX(JI)
  IDPL(JI)      = ISDPL(JL)
  IPBL(JI)      = ISPBL(JL)
  ILCL(JI)      = ISLCL(JL)
  ZTHLCL(JI)    = ZSTHLCL(JL)
  ZTLCL(JI)     = ZSTLCL(JL)
  ZRVLCL(JI)    = ZSRVLCL(JL)
  ZWLCL(JI)     = ZSWLCL(JL)
  ZZLCL(JI)     = ZSZLCL(JL)
  ZTHVELCL(JI)  = ZSTHVELCL(JL)
  ZDXDY(JI)     = ZSDXDY(JL)
END DO

DEALLOCATE( GTRIG1 )
ALLOCATE( GTRIG1(ICONV) )
GTRIG1(:) = .true.

DEALLOCATE( ISDPL )
DEALLOCATE( ISPBL )
DEALLOCATE( ISLCL )
DEALLOCATE( ZSTHLCL )
DEALLOCATE( ZSTLCL )
DEALLOCATE( ZSRVLCL )
DEALLOCATE( ZSWLCL )
DEALLOCATE( ZSZLCL )
DEALLOCATE( ZSTHVELCL )
DEALLOCATE( ZSDXDY )
!
!
!*           3.2    Compute pressure difference
!                   ---------------------------------------------------
!
ZDPRES(:,IKB) = 0.
DO JK = IKB + 1, IKE
  ZDPRES(:,JK)  = ZPRES(:,JK-1) - ZPRES(:,JK)
END DO
!
!*           3.3   Compute environm. enthalpy and total water = r_v + r_i + r_c 
!                  ----------------------------------------------------------
!
DO JK = IKB, IKE, 1
  ZRW(:,JK)  = ZRV(:,JK) + ZRC(:,JK) + ZRI(:,JK)
  ZCPH(:)    = XCPD + XCPV * ZRW(:,JK)
  ZLV(:)     = XLVTT + ( XCPV - XCL ) * ( ZTT(:,JK) - XTT ) ! compute L_v
  ZLS(:)     = XLSTT + ( XCPV - XCI ) * ( ZTT(:,JK) - XTT ) ! compute L_i
  ZTHL(:,JK) = ZCPH(:) * ZTT(:,JK) + ( 1. + ZRW(:,JK) ) * XG * ZZ(:,JK) &
               - ZLV(:) * ZRC(:,JK) - ZLS(:) * ZRI(:,JK)
END DO
!
DEALLOCATE( ZCPH )
DEALLOCATE( ZLV )
DEALLOCATE( ZLS )
!
!-------------------------------------------------------------------------------
!
!*           4.     Compute updraft properties 
!                   ----------------------------
!
!*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s ) 
!                   -------------------------------------------------------------
!
ZDXDY(:)  = XA25
ZMFLCL(:) = XA25 * 1.E-3
!
!
!
CALL CONVECT_UPDRAFT_SHAL( ICONV, KLEV,                                     &
                           KICE, ZPRES, ZDPRES, ZZ, ZTHL, ZTHV, ZTHES, ZRW, &
                           ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL,   & 
                           ZMFLCL, GTRIG1, ILCL, IDPL, IPBL,                &
                           ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZURW,            &
                           ZURC, ZURI, ZCAPE, ICTL, IETL                    )
!
!
!
!*           4.2    In routine UPDRAFT GTRIG1 has been set to false when cloud 
!                   thickness is smaller than 3 km
!                   -----------------------------------------------------------
!
!
!
!*       4.3    Allocate memory for downdraft variables
!               ---------------------------------------
!
! downdraft variables
!
  ALLOCATE( ZDMF(ICONV,IKS) )
  ALLOCATE( ZDER(ICONV,IKS) )
  ALLOCATE( ZDDR(ICONV,IKS) )
  ALLOCATE( ILFS(ICONV) )
  ALLOCATE( ZLMASS(ICONV,IKS) )
  ZDMF(:,:) = 0.
  ZDER(:,:) = 0.
  ZDDR(:,:) = 0.
  ILFS(:)   = IKB
  DO JK = IKB, IKE
    ZLMASS(:,JK)  = ZDXDY(:) * ZDPRES(:,JK) / XG  ! mass of model layer
  END DO
  ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
!
! closure variables
!
  ALLOCATE( ZTIMEC(ICONV) )
  ALLOCATE( ZTHC(ICONV,IKS) )
  ALLOCATE( ZRVC(ICONV,IKS) )
  ALLOCATE( ZRCC(ICONV,IKS) )
  ALLOCATE( ZRIC(ICONV,IKS) )
  ALLOCATE( ZWSUB(ICONV,IKS) )
!
!-------------------------------------------------------------------------------
!
!*           5.     Compute downdraft properties 
!                   ----------------------------
!
  ZTIMEC(:) = XCTIME_SHAL
  IF ( OSETTADJ ) ZTIMEC(:) = PTADJS
!
!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------
!
  CALL CONVECT_CLOSURE_SHAL( ICONV, KLEV,                         &
                             ZPRES, ZDPRES, ZZ, ZDXDY, ZLMASS,    &
                             ZTHL, ZTH, ZRW, ZRC, ZRI, GTRIG1,    &
                             ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB,       &
                             ILCL, IDPL, IPBL, ICTL,              &
                             ZUMF, ZUER, ZUDR, ZUTHL, ZURW,       &
                             ZURC, ZURI, ZCAPE, ZTIMEC, IFTSTEPS  )
!
!-------------------------------------------------------------------------------
!
!*           8.     Determine the final grid-scale (environmental) convective 
!                   tendencies and set convective counter
!                   --------------------------------------------------------
!
!
!*           8.1    Grid scale tendencies
!                   ---------------------
!
          ! in order to save memory, the tendencies are temporarily stored
          ! in the tables for the adjusted grid-scale values
!
  DO JK = IKB, IKE
     ZTHC(:,JK) = ( ZTHC(:,JK) - ZTH(:,JK) ) / ZTIMEC(:)             &
       * ( ZPRES(:,JK) / XP00 ) ** ZRDOCP ! change theta in temperature
     ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + ZRC(:,JK) + ZRI(:,JK) ) &
                                          / ZTIMEC(:) 

     ZRCC(:,JK) = ( ZRCC(:,JK) - ZRC(:,JK) ) / ZTIMEC(:)
     ZRIC(:,JK) = ( ZRIC(:,JK) - ZRI(:,JK) ) / ZTIMEC(:)
!
  END DO
!
!
!*           8.2    Apply conservation correction
!                   -----------------------------
!
          ! adjustment at cloud top to smooth possible discontinuous profiles at PBL inversions
          ! (+ - - tendencies for moisture )
!
!
IF (LLSMOOTH) THEN
  DO JI = 1, ICONV
     JK = ICTL(JI)
     JKM= MAX(2,ICTL(JI)-1)
     JKP= MAX(2,ICTL(JI)-2)
     ZRVC(JI,JKM) = ZRVC(JI,JKM) + .5 * ZRVC(JI,JK)
     ZRCC(JI,JKM) = ZRCC(JI,JKM) + .5 * ZRCC(JI,JK)
     ZRIC(JI,JKM) = ZRIC(JI,JKM) + .5 * ZRIC(JI,JK)
     ZTHC(JI,JKM) = ZTHC(JI,JKM) + .5 * ZTHC(JI,JK)
     ZRVC(JI,JKP) = ZRVC(JI,JKP) + .3 * ZRVC(JI,JK)
     ZRCC(JI,JKP) = ZRCC(JI,JKP) + .3 * ZRCC(JI,JK)
     ZRIC(JI,JKP) = ZRIC(JI,JKP) + .3 * ZRIC(JI,JK)
     ZTHC(JI,JKP) = ZTHC(JI,JKP) + .3 * ZTHC(JI,JK)
     ZRVC(JI,JK)  = .2 * ZRVC(JI,JK)
     ZRCC(JI,JK)  = .2 * ZRCC(JI,JK)
     ZRIC(JI,JK)  = .2 * ZRIC(JI,JK)
     ZTHC(JI,JK)  = .2 * ZTHC(JI,JK)
  END DO
ENDIF
!
!
          ! Compute vertical integrals - Fluxes
!
  JKM = MAXVAL( ICTL(:) )
  ZWORK2(:) = 0.
  ZWORK2B(:) = 0.
  DO JK = IKB+1, JKM
    JKP = JK + 1
    DO JI = 1, ICONV
      IF ( JK <= ICTL(JI) ) THEN
      ZW1 =  ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK)
      ZWORK2(JI) = ZWORK2(JI) +  ZW1 *          & ! moisture
                                  .5 * (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / XG
      ZW1 = ( XCPD + XCPV * ZRW(JI,JK) )* ZTHC(JI,JK)   - &
            ( XLVTT + ( XCPV - XCL ) * ( ZTT(JI,JK) - XTT ) ) * ZRCC(JI,JK) - &
            ( XLSTT + ( XCPV - XCL ) * ( ZTT(JI,JK) - XTT ) ) * ZRIC(JI,JK)
      ZWORK2B(JI) = ZWORK2B(JI) + ZW1 *         & ! energy
                                  .5 * (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / XG
      END IF
    END DO
  END DO
!
          ! Budget error (integral must be zero)
!
  DO JI = 1, ICONV
    IF ( ICTL(JI) > IKB+1 ) THEN
      JKP = ICTL(JI)
      ZW1 = XG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - &
                .5 * (ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1)) )
      ZWORK2(JI) =  ZWORK2(JI) * ZW1
      ZWORK2B(JI) = ZWORK2B(JI)* ZW1
    END IF
  END DO
!
          ! Apply uniform correction
!
  DO JK = JKM, IKB+1, -1
  DO JI = 1, ICONV
    IF ( ICTL(JI) > IKB+1 .AND. JK <= ICTL(JI) ) THEN
      ! ZW1 = ABS(ZRVC(JI,JK)) +  ABS(ZRCC(JI,JK)) +  ABS(ZRIC(JI,JK)) + 1.E-12
      ! ZRVC(JI,JK) = ZRVC(JI,JK) - ABS(ZRVC(JI,JK))/ZW1*ZWORK2(JI)           ! moisture
      ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)                                ! moisture
      ! ZRCC(JI,JK) = ZRCC(JI,JK) - ABS(ZRCC(JI,JK))/ZW1*ZWORK2(JI)
      ! ZRIC(JI,JK) = ZRIC(JI,JK) - ABS(ZRIC(JI,JK))/ZW1*ZWORK2(JI)
      ZTHC(JI,JK) = ZTHC(JI,JK) - ZWORK2B(JI) /  XCPD                       ! enthalpy
    END IF
  END DO
  END DO
!
	      ! execute a "scatter"= pack command to store the tendencies in
	      ! the final 2D tables
!
  DO JK = IKB, IKE
  DO JI = 1, ICONV
    JL = IJINDEX(JI)
    PTTEN(JL,JK)   = ZTHC(JI,JK)
    PRVTEN(JL,JK)  = ZRVC(JI,JK)
    PRCTEN(JL,JK)  = ZRCC(JI,JK)
    PRITEN(JL,JK)  = ZRIC(JI,JK)
  END DO
  END DO
!
!
!                   Cloud base and top levels
!                   -------------------------
!
  ILCL(:) = MIN( ILCL(:), ICTL(:) )
  DO JI = 1, ICONV
    JL = IJINDEX(JI)
    KCLTOP(JL) = ICTL(JI)
    KCLBAS(JL) = ILCL(JI)
  END DO
!
!
!*           8.7    Compute convective tendencies for Tracers
!                   ------------------------------------------
!
  IF ( OCH1CONV ) THEN
!
    ALLOCATE( ZCH1(ICONV,IKS,KCH1) )
    ALLOCATE( ZCH1C(ICONV,IKS,KCH1) )
    ALLOCATE( ZWORK3(ICONV,KCH1) )
!
    DO JK = IKB, IKE
    DO JI = 1, ICONV
      JL = IJINDEX(JI)
      ZCH1(JI,JK,:) = PCH1(JL,JK,:)
    END DO
    END DO
!
    CALL CONVECT_CHEM_TRANSPORT( ICONV, KLEV, KCH1, ZCH1, ZCH1C,          &
                                 IDPL, IPBL, ILCL, ICTL, ILFS, ILFS,      &
                                 ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,      &
                                 ZTIMEC, ZDXDY, ZDMF(:,1), ZLMASS, ZWSUB, &
                                 IFTSTEPS )
!
!
!*           8.8    Apply conservation correction
!                   -----------------------------
!
          ! Compute vertical integrals
!
    JKM = MAXVAL( ICTL(:) )
    DO JN = 1, KCH1
      IF(JN < NSV_LGBEG .OR. JN>NSV_LGEND-1) THEN ! no correction for xy lagrangian variables
        ZWORK3(:,JN) = 0.
        ZWORK2(:)    = 0.
        DO JK = IKB+1, JKM
          JKP = JK + 1
          DO JI = 1, ICONV
            ZW1 = .5 * (ZPRES(JI,JK-1) - ZPRES(JI,JKP))
            ZWORK3(JI,JN) = ZWORK3(JI,JN) + (ZCH1C(JI,JK,JN)-ZCH1(JI,JK,JN)) * ZW1
            ZWORK2(JI)    = ZWORK2(JI)    + ABS(ZCH1C(JI,JK,JN)) * ZW1
          END DO
        END DO
!
             ! Apply concentration weighted correction
!
        DO JK = JKM, IKB+1, -1
          DO JI = 1, ICONV
            IF ( ICTL(JI) > IKB+1 .AND. JK <= ICTL(JI) ) THEN
              ZCH1C(JI,JK,JN) = ZCH1C(JI,JK,JN) -   &
                                ZWORK3(JI,JN)*ABS(ZCH1C(JI,JK,JN))/MAX(1.E-30,ZWORK2(JI))
            END IF
          END DO
        END DO
      END IF
!
      DO JK = IKB, IKE
        DO JI = 1, ICONV
          JL = IJINDEX(JI)
          PCH1TEN(JL,JK,JN) = (ZCH1C(JI,JK,JN)-ZCH1(JI,JK,JN) ) / ZTIMEC(JI)
        END DO
      END DO
    END DO
  END IF
!
!-------------------------------------------------------------------------------
!
!*           9.     Write up- and downdraft mass fluxes 
!                   ------------------------------------
!
  DO JK = IKB, IKE
    ZUMF(:,JK)  = ZUMF(:,JK) / ZDXDY(:) ! Mass flux per unit area
  END DO
  ZWORK2(:) = 1.
  DO JK = IKB, IKE
  DO JI = 1, ICONV
    JL = IJINDEX(JI)
    IF ( KCLTOP(JL) <= IKB+1 ) ZWORK2(JL) = 0.
    PUMF(JL,JK) = ZUMF(JI,JK) * ZWORK2(JL)
  END DO
  END DO
!
!-------------------------------------------------------------------------------
!
!*           10.    Deallocate all local arrays
!                   ---------------------------
!
! downdraft variables
!
  DEALLOCATE( ZDMF )
  DEALLOCATE( ZDER )
  DEALLOCATE( ZDDR )
  DEALLOCATE( ILFS )
  DEALLOCATE( ZLMASS )
!
!   closure variables
!
  DEALLOCATE( ZTIMEC )
  DEALLOCATE( ZTHC )
  DEALLOCATE( ZRVC )
  DEALLOCATE( ZRCC )
  DEALLOCATE( ZRIC )
  DEALLOCATE( ZWSUB )
!
  IF ( OCH1CONV ) THEN
    DEALLOCATE( ZCH1 )
    DEALLOCATE( ZCH1C )
    DEALLOCATE( ZWORK3 )
  END IF
!
!    vertical index
!
DEALLOCATE( IDPL )
DEALLOCATE( IPBL )
DEALLOCATE( ILCL )
DEALLOCATE( ICTL )
DEALLOCATE( IETL )
!
! grid scale variables
!
DEALLOCATE( ZZ )
DEALLOCATE( ZPRES )
DEALLOCATE( ZDPRES )
DEALLOCATE( ZTT )
DEALLOCATE( ZTH )
DEALLOCATE( ZTHV )
DEALLOCATE( ZTHL )
DEALLOCATE( ZTHES )
DEALLOCATE( ZRW )
DEALLOCATE( ZRV )
DEALLOCATE( ZRC )
DEALLOCATE( ZRI )
DEALLOCATE( ZDXDY )
!
! updraft variables
!
DEALLOCATE( ZUMF )
DEALLOCATE( ZUER )
DEALLOCATE( ZUDR )
DEALLOCATE( ZUTHL )
DEALLOCATE( ZUTHV )
DEALLOCATE( ZURW )
DEALLOCATE( ZURC )
DEALLOCATE( ZURI )
DEALLOCATE( ZTHLCL )
DEALLOCATE( ZTLCL )
DEALLOCATE( ZRVLCL )
DEALLOCATE( ZWLCL )
DEALLOCATE( ZZLCL )
DEALLOCATE( ZTHVELCL )
DEALLOCATE( ZMFLCL )
DEALLOCATE( ZCAPE )
!
! work arrays
!
DEALLOCATE( IJINDEX )
DEALLOCATE( GTRIG1 )
!
!
END SUBROUTINE SHALLOW_CONVECTION
