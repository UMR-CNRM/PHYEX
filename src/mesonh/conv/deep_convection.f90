!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!     ######################
      MODULE MODI_DEEP_CONVECTION
!     ######################
!
INTERFACE
!
    SUBROUTINE DEEP_CONVECTION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                                PDTCONV, KICE, OREFRESH, ODOWN, OSETTADJ,      &
                                PPABST, PZZ, PDXDY, PTIMEC,                    &
                                PTT, PRVT, PRCT, PRIT, PUT, PVT, PWT,          &
                                KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,         &
                                PPRLTEN, PPRSTEN,                              &
                                KCLTOP, KCLBAS, PPRLFLX, PPRSFLX,              &
                                PUMF, PDMF, PCAPE,                             &
                                OCH1CONV, KCH1, PCH1, PCH1TEN,                 &
                                OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX,       &
                                ODUST, OSALT, PRHODREF, PIC_RATE, PCG_RATE     )

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
LOGICAL,                    INTENT(IN) :: OREFRESH ! refresh or not tendencies
                                                   ! at every call
LOGICAL,                    INTENT(IN) :: ODOWN    ! take or not convective
                                                   ! downdrafts into account
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
                                                   ! adjustment time by user 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUT      ! grid scale horiz. wind u "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PVT      ! grid scale horiz. wind v "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical 
                                                   ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PDXDY    ! horizontal grid area (m-a2)
REAL, DIMENSION(KLON),      INTENT(IN) :: PTIMEC   ! value of convective adjustment
                                                   ! time if OSETTADJ=.TRUE.
!   
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCOUNT ! convective counter (recompute
                                                   ! tendency or keep it)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature
                                                   ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRLTEN! liquid surf. precipitation
                                                   ! tendency (m/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRSTEN! solid surf. precipitation
                                                   ! tendency (m/s)
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRLFLX! liquid precip flux (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRSFLX! solid  precip flux (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF   ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PCAPE  ! maximum CAPE (J/kg)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
LOGICAL,                    INTENT(IN) :: OUSECHEM      ! flag for chemistry 
LOGICAL,                    INTENT(IN) :: OCH_CONV_SCAV !  & scavenging
LOGICAL,                    INTENT(IN) :: OCH_CONV_LINOX ! & LiNOx
LOGICAL,                    INTENT(IN) :: ODUST         ! flag for dust
LOGICAL,                    INTENT(IN) :: OSALT         ! flag for sea salt
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRHODREF      ! grid scale density
REAL, DIMENSION(KLON), INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(KLON), INTENT(INOUT) :: PCG_RATE ! CG lightning frequency

!
END SUBROUTINE DEEP_CONVECTION
!
END INTERFACE
!
END MODULE MODI_DEEP_CONVECTION
!
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 conv 2006/09/21 10:55:01
!-----------------------------------------------------------------
!   ############################################################################
    SUBROUTINE DEEP_CONVECTION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                                PDTCONV, KICE, OREFRESH, ODOWN, OSETTADJ,      &
                                PPABST,  PZZ, PDXDY, PTIMEC,                   &
                                PTT, PRVT, PRCT, PRIT, PUT, PVT, PWT,          &
                                KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,         &
                                PPRLTEN, PPRSTEN,                              &
                                KCLTOP, KCLBAS, PPRLFLX, PPRSFLX,              &
                                PUMF, PDMF, PCAPE,                             &
                                OCH1CONV, KCH1, PCH1, PCH1TEN,                 &
                                OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX,       &
                                ODUST, OSALT, PRHODREF, PIC_RATE, PCG_RATE     )
!   ############################################################################
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
!!    CONVECT_TRIGGER_FUNCT 
!!    CONVECT_SATMIXRATIO
!!    CONVECT_UPDRAFT
!!        CONVECT_CONDENS
!!        CONVECT_MIXING_FUNCT
!!    CONVECT_TSTEP_PREF
!!    CONVECT_DOWNDRAFT
!!    CONVECT_PRECIP_ADJUST
!!    CONVECT_CLOSURE
!!        CONVECT_CLOSURE_THRVLCL
!!        CONVECT_CLOSURE_ADJUST
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
!!      Bechtold et al., 2001, Quart. J. Roy. Met. Soc.
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
!!   Peter Bechtold 04/10/97 replace theta_il by enthalpy
!!         "        10/12/98 changes for ARPEGE
!!         "        12/12/00 add conservation correction
!!      C. Mari     13/02/01 add scavenging of chemical species in updraft
!!     P. Jabouille 02/07/01 case of lagragian variables
!!     P. Tulet     02/03/05 update for dust
!!     C.Lac        27/09/10 modification loop index for reproducibility
!!    Juan 24/09/2012: for BUG Pgi rewrite PACK function on mode_pack_pgi
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_CONVPAREXT
USE MODD_CONVPAR
USE MODD_NSV,       ONLY : NSV_LGBEG,NSV_LGEND, &
                           NSV_CHEMBEG,NSV_CHEMEND, &
                           NSV_LNOXBEG
USE MODD_CH_M9_n,   ONLY : CNAMES
!
USE MODI_CH_CONVECT_LINOX
USE MODI_CONVECT_TRIGGER_FUNCT
USE MODI_CONVECT_UPDRAFT
USE MODI_CONVECT_TSTEP_PREF
USE MODI_CONVECT_DOWNDRAFT
USE MODI_CONVECT_PRECIP_ADJUST
USE MODI_CONVECT_CLOSURE
USE MODI_CH_CONVECT_SCAVENGING
USE MODI_CONVECT_CHEM_TRANSPORT
!
#ifdef MNH_PGI
USE MODE_PACK_PGI
#endif
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
LOGICAL,                    INTENT(IN) :: OREFRESH ! refresh or not tendencies
                                                   ! at every call
LOGICAL,                    INTENT(IN) :: ODOWN    ! take or not convective
                                                   ! downdrafts into account
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
                                                   ! adjustment time by user 
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUT      ! grid scale horiz. wind u "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PVT      ! grid scale horiz. wind v "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical 
                                                   ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m) 
REAL, DIMENSION(KLON),      INTENT(IN) :: PDXDY    ! horizontal grid area (m-a2)
REAL, DIMENSION(KLON),      INTENT(IN) :: PTIMEC   ! value of convective adjustment
                                                   ! time if OSETTADJ=.TRUE.
!   
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCOUNT ! convective counter (recompute
                                                   ! tendency or keep it)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature
                                                   ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRLTEN! liquid surf. precipitation
                                                   ! tendency (m/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRSTEN! solid surf. precipitation
                                                   ! tendency (m/s)
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRLFLX! liquid precip flux (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRSFLX! solid  precip flux (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF   ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PCAPE  ! maximum CAPE (J/kg)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
LOGICAL,                    INTENT(IN) :: OUSECHEM      ! flag for chemistry 
LOGICAL,                    INTENT(IN) :: OCH_CONV_SCAV !  & scavenging
LOGICAL,                    INTENT(IN) :: OCH_CONV_LINOX ! & LiNOx
LOGICAL,                    INTENT(IN) :: ODUST         ! flag for dust
LOGICAL,                    INTENT(IN) :: OSALT         ! flag for sea salt
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRHODREF      ! grid scale density
REAL, DIMENSION(KLON), INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(KLON), INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
!
!
!*       0.2   Declarations of local fixed memory variables :
!
INTEGER  :: ITEST, ICONV, ICONV1    ! number of convective columns
INTEGER  :: IIB, IIE                ! horizontal loop bounds
INTEGER  :: IKB, IKE                ! vertical loop bounds
INTEGER  :: IKS                     ! vertical dimension
INTEGER  :: JI, JL, JJ              ! horizontal loop index
INTEGER  :: JN                      ! number of tracers
INTEGER  :: JK, JKP, JKM            ! vertical loop index
INTEGER  :: IFTSTEPS                ! only used for chemical tracers
REAL     :: ZEPS, ZEPSA             ! R_d / R_v, R_v / R_d
REAL     :: ZRDOCP                  ! R_d/C_p
!
LOGICAL, DIMENSION(KLON, KLEV)     :: GTRIG3 ! 3D logical mask for convection 
LOGICAL, DIMENSION(KLON)           :: GTRIG  ! 2D logical mask for trigger test
REAL, DIMENSION(KLON,KLEV)         :: ZTHT, ZSTHV, ZSTHES  ! grid scale theta, 
                                                           ! theta_v, theta_es
REAL, DIMENSION(KLON)              :: ZTIME  ! convective time period
REAL, DIMENSION(KLON)              :: ZWORK2, ZWORK2B ! work array 
REAL                               :: ZW1    ! work variable
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
INTEGER, DIMENSION(:),ALLOCATABLE  :: IDBL    ! index for downdraft base level  
INTEGER, DIMENSION(:),ALLOCATABLE  :: IML     ! melting level  
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
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZU      ! grid scale horiz. u component on theta grid
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZV      ! grid scale horiz. v component on theta grid
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
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUPR    ! updraft precipitation in
                                              ! flux units (kg water / s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTHL   ! updraft enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTHV   ! updraft theta_v (K)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZUTT    ! updraft temperature (K)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURW    ! updraft total water (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURC    ! updraft cloud water (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURI    ! updraft cloud ice   (kg/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURR    ! liquid precipit. (kg/kg)
                                              ! produced in  model layer
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZURS    ! solid precipit. (kg/kg)
                                              ! produced in  model layer
REAL, DIMENSION(:),   ALLOCATABLE  :: ZUTPR   ! total updraft precipitation (kg/s)
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
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDTHL   ! downdraft enthalpy (J/kg)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDRW    ! downdraft total water (kg/kg)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZMIXF   ! mixed fraction at LFS        
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTPR    ! total surf precipitation (kg/s)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZSPR    ! solid surf precipitation (kg/s)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZDTEVR  ! donwndraft evapor. (kg/s)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZPREF   ! precipitation efficiency
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZDTEVRF ! donwndraft evapor. (kg/s)
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPRLFLX ! liquid precip flux
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZPRSFLX ! solid precip flux
!
! closure variables
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZLMASS  ! mass of model layer (kg)
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTIMEA  ! advective time period
REAL, DIMENSION(:),   ALLOCATABLE  :: ZTIMEC, ZTIMED! time during which convection is
                                              ! active at grid point (as ZTIME)
!
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZTHC    ! conv. adj. grid scale theta
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRVC    ! conv. adj. grid scale r_w 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRCC    ! conv. adj. grid scale r_c 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZRIC    ! conv. adj. grid scale r_i 
REAL, DIMENSION(:,:), ALLOCATABLE  :: ZWSUB   ! envir. compensating subsidence (Pa/s)
!
LOGICAL, DIMENSION(:),ALLOCATABLE  :: GTRIG1  ! logical mask for convection    
LOGICAL, DIMENSION(:),ALLOCATABLE  :: GWORK   ! logical work array
INTEGER, DIMENSION(:),ALLOCATABLE  :: IINDEX, IJINDEX, IJSINDEX, IJPINDEX!hor.index
REAL, DIMENSION(:),   ALLOCATABLE  :: ZCPH    ! specific heat C_ph 
REAL, DIMENSION(:),   ALLOCATABLE  :: ZLV, ZLS! latent heat of vaporis., sublim.
REAL                               :: ZES     ! saturation vapor mixng ratio
!
! Chemical Tracers:
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1    ! grid scale chemical specy (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE:: ZCH1C   ! conv. adjust. chemical specy 1
REAL, DIMENSION(:,:),   ALLOCATABLE:: ZWORK3  ! work array 
LOGICAL, DIMENSION(:,:,:),ALLOCATABLE::GTRIG4 ! logical mask
INTEGER                            :: JN_NO   ! index of NO compound in PCH1
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZWORK4, ZWORK4C
                                         ! LiNOx conc. and tendency
REAL, DIMENSION(:,:),ALLOCATABLE   :: ZZZ, ZRHODREF
REAL, DIMENSION(:),ALLOCATABLE     :: ZIC_RATE,ZCG_RATE
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
JCVEXT = MAX( 0, KTDIA - 1 )
IKE = IKS - JCVEXT
!
!
!*       0.5    Update convective counter ( where KCOUNT > 0 
!               convection is still active ).
!               ---------------------------------------------
!
KCOUNT(IIB:IIE) = KCOUNT(IIB:IIE) - 1 
!
IF ( OREFRESH ) THEN
  KCOUNT(:) = 1
  KCOUNT(IIB:IIE) = 0 ! refresh or not at every call
END IF
!
GTRIG(:)  = KCOUNT(:) <= 0
ITEST = COUNT( GTRIG(:) )
IF ( ITEST == 0 )  THEN   ! if convection is already active at every grid point
  RETURN 
ENDIF
                          ! exit DEEP_CONVECTION
!
!
!*       0.7    Reset convective tendencies to zero if convective
!               counter becomes negative
!               -------------------------------------------------
!
DO JJ=1,KLEV ; DO JI=1,KLON
 GTRIG3(JI,JJ)=GTRIG(JI)
ENDDO ; ENDDO 
WHERE ( GTRIG3(:,:) ) 
  PTTEN(:,:)  = 0.
  PRVTEN(:,:) = 0.
  PRCTEN(:,:) = 0.
  PRITEN(:,:) = 0.
  PPRLFLX(:,:)= 0.
  PPRSFLX(:,:)= 0.
! PUTEN(:,:)  = 0.
! PVTEN(:,:)  = 0.
  PUMF(:,:)   = 0.
  PDMF(:,:)   = 0.
END WHERE
WHERE ( GTRIG(:) ) 
  PPRLTEN(:) = 0.
  PPRSTEN(:) = 0.
  KCLTOP(:)  = 0
  KCLBAS(:)  = 0
  PCAPE(:)   = 0.
END WHERE
ALLOCATE( GTRIG4(KLON,KLEV,KCH1) )
DO JK=1,KCH1; DO JJ=1,KLEV ; DO JI=1,KLON
!GTRIG4(:,:,:) = SPREAD( GTRIG3(:,:), DIM=3, NCOPIES=KCH1 )
  GTRIG4(JI,JJ,JK) = GTRIG3(JI,JJ)
ENDDO ; ENDDO ; ENDDO
WHERE( GTRIG4(:,:,:) ) PCH1TEN(:,:,:) = 0.
DEALLOCATE( GTRIG4 )
!
!-------------------------------------------------------------------------------
!
!*       1.     Initialize  local variables
!               ----------------------------
!
ZEPS   = XRD / XRV
ZEPSA  = XRV / XRD 
ZRDOCP = XRD / XCPD
!
!
!*       1.1    Set up grid scale theta, theta_v, theta_es 
!               ------------------------------------------
!
ZTHT(:,:) = 300.
ZSTHV(:,:)= 300.
ZSTHES(:,:)=400.
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
    ZES = ZEPS * ZES / ( PPABST(JI,JK) - ZES )
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
ALLOCATE( ZPRES(ITEST,IKS) )
ALLOCATE( ZZ(ITEST,IKS) )
ALLOCATE( ZW(ITEST,IKS) )
ALLOCATE( ZTH(ITEST,IKS) )
ALLOCATE( ZTHV(ITEST,IKS) )
ALLOCATE( ZTHEST(ITEST,IKS) )
ALLOCATE( ZRV(ITEST,IKS) )
ALLOCATE( ZSTHLCL(ITEST) )
ALLOCATE( ZSTLCL(ITEST) )
ALLOCATE( ZSRVLCL(ITEST) )
ALLOCATE( ZSWLCL(ITEST) )
ALLOCATE( ZSZLCL(ITEST) )
ALLOCATE( ZSTHVELCL(ITEST) )
ALLOCATE( ISDPL(ITEST) )
ALLOCATE( ISPBL(ITEST) )
ALLOCATE( ISLCL(ITEST) )
ALLOCATE( ZSDXDY(ITEST) )
ALLOCATE( GTRIG1(ITEST) )
ALLOCATE( ZCAPE(ITEST) )
ALLOCATE( IINDEX(KLON) )
ALLOCATE( IJSINDEX(ITEST) )
DO JI = 1, KLON
  IINDEX(JI) = JI
END DO
IJSINDEX(:) = PACK( IINDEX(:), MASK=GTRIG(:) )
!
ZPRES = 0.     
ZZ    = 0.
ZTH   = 0.
ZTHV  = 0.
ZTHEST = 0.
ZRV   = 0.
ZW    = 0.
!
DO JK = IKB, IKE
DO JI = 1, ITEST
  JL = IJSINDEX(JI)
  ZPRES(JI,JK)  = PPABST(JL,JK)
  ZZ(JI,JK)     = PZZ(JL,JK)
  ZTH(JI,JK)    = ZTHT(JL,JK)
  ZTHV(JI,JK)   = ZSTHV(JL,JK)
  ZTHEST(JI,JK) = ZSTHES(JL,JK)
  ZRV(JI,JK)    = MAX( 0., PRVT(JL,JK) )
  ZW(JI,JK)     = PWT(JL,JK)
END DO
END DO
DO JI = 1, ITEST
  JL = IJSINDEX(JI)
  ZSDXDY(JI)    = PDXDY(JL)
END DO
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
!
CALL CONVECT_TRIGGER_FUNCT( ITEST, KLEV,                              &
                            ZPRES, ZTH, ZTHV, ZTHEST,                 &
                            ZRV, ZW, ZZ, ZSDXDY,                      &
                            ZSTHLCL, ZSTLCL, ZSRVLCL, ZSWLCL, ZSZLCL, &
                            ZSTHVELCL, ISLCL, ISDPL, ISPBL, GTRIG1,   &
                            ZCAPE )
!
DO JI = 1, ITEST
  JL = IJSINDEX(JI)
  PCAPE(JL) = ZCAPE(JI)
END DO
!
DEALLOCATE( ZPRES )
DEALLOCATE( ZZ )
DEALLOCATE( ZTH )
DEALLOCATE( ZTHV )
DEALLOCATE( ZTHEST )
DEALLOCATE( ZRV )
DEALLOCATE( ZW )
DEALLOCATE( ZCAPE )
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
  DEALLOCATE( IINDEX )
  DEALLOCATE( IJSINDEX )
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
ALLOCATE( ZZ(ICONV,IKS) )     ;   ZZ  = 0.0
ALLOCATE( ZPRES(ICONV,IKS) )  ;   ZPRES = 0.0
ALLOCATE( ZDPRES(ICONV,IKS) ) ;   ZDPRES = 0.0
ALLOCATE( ZU(ICONV,IKS) )     ;   ZU     = 0.0
ALLOCATE( ZV(ICONV,IKS) )     ;   ZV     = 0.0
ALLOCATE( ZTT(ICONV, IKS) )   ;   ZTT    = 0.0
ALLOCATE( ZTH(ICONV,IKS) )    ;   ZTH    = 0.0
ALLOCATE( ZTHV(ICONV,IKS) )   ;   ZTHV   = 0.0
ALLOCATE( ZTHL(ICONV,IKS) )    ; ZTHL  = 0.0
ALLOCATE( ZTHES(ICONV,IKS) )   ; ZTHES = 0.0
ALLOCATE( ZRV(ICONV,IKS) )     ; ZRV   = 0.0
ALLOCATE( ZRC(ICONV,IKS) )     ; ZRC   = 0.0
ALLOCATE( ZRI(ICONV,IKS) )     ; ZRI   = 0.0
ALLOCATE( ZRW(ICONV,IKS) )     ; ZRW   = 0.0
ALLOCATE( ZDXDY(ICONV) )       ; ZDXDY = 0.0
!
         ! updraft variables
!
ALLOCATE( ZUMF(ICONV,IKS) )
ALLOCATE( ZUER(ICONV,IKS) )
ALLOCATE( ZUDR(ICONV,IKS) )
ALLOCATE( ZUPR(ICONV,IKS) )
ALLOCATE( ZUTHL(ICONV,IKS) )
ALLOCATE( ZUTHV(ICONV,IKS) )
ALLOCATE( ZUTT(ICONV,IKS) )
ALLOCATE( ZURW(ICONV,IKS) )
ALLOCATE( ZURC(ICONV,IKS) )
ALLOCATE( ZURI(ICONV,IKS) )
ALLOCATE( ZURR(ICONV,IKS) )
ALLOCATE( ZURS(ICONV,IKS) )
ALLOCATE( ZUTPR(ICONV) )
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
ALLOCATE( IJPINDEX(ICONV) )
ALLOCATE( ZCPH(ICONV) )
ALLOCATE( ZLV(ICONV) )
ALLOCATE( ZLS(ICONV) )
!
!
!*           3.1    Gather grid scale and updraft base variables in
!                   arrays using mask GTRIG
!                   ---------------------------------------------------
!
GTRIG(:)      = UNPACK( GTRIG1(:), MASK=GTRIG, FIELD=.FALSE. )  
IJINDEX(:)    = PACK( IINDEX(:), MASK=GTRIG(:) )
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
  ZU(JI,JK)     = PUT(JL,JK)
  ZV(JI,JK)     = PVT(JL,JK)
END DO
END DO
IF ( OSETTADJ ) THEN
  ALLOCATE( ZTIMED(ICONV) )
  DO JI = 1, ICONV
    JL = IJINDEX(JI)
    ZTIMED(JI) = PTIMEC(JL)
   END DO
END IF
!
DO JI = 1, ITEST
  IJSINDEX(JI) = JI
END DO
IJPINDEX(:) = PACK( IJSINDEX(:), MASK=GTRIG1(:) )
DO JI = 1, ICONV
  JL = IJPINDEX(JI)
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
ALLOCATE( GWORK(ICONV) )
GWORK(:)      = PACK( GTRIG1(:),  MASK=GTRIG1(:) ) 
DEALLOCATE( GTRIG1 )
ALLOCATE( GTRIG1(ICONV) )
GTRIG1(:)     = GWORK(:)
!                 
DEALLOCATE( GWORK )
DEALLOCATE( IJPINDEX )
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
!            ---------------------------------------------------
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
!-------------------------------------------------------------------------------
!
!*           4.     Compute updraft properties 
!                   ----------------------------
!
!*           4.1    Set mass flux at LCL ( here a unit mass flux with w = 1 m/s ) 
!                   -------------------------------------------------------------
!
DO JI = 1, ICONV
  JK = ILCL(JI) - 1
  ZMFLCL(JI) = ZPRES(JI,JK) / ( XRD * ZTT(JI,JK) *                &
            ( 1. + ZEPS * ZRVLCL(JI) ) ) * XPI * XCRAD * XCRAD    &
              * MAX ( 1., ZDXDY(JI)/XA25 )
END DO
!
DEALLOCATE( ZCPH )
DEALLOCATE( ZLV )
DEALLOCATE( ZLS )
!
!
CALL CONVECT_UPDRAFT( ICONV, KLEV,                                     &
                      KICE, ZPRES, ZDPRES, ZZ, ZTHL, ZTHV, ZTHES, ZRW, &
                      ZTHLCL, ZTLCL, ZRVLCL, ZWLCL, ZZLCL, ZTHVELCL,   & 
                      ZMFLCL, GTRIG1, ILCL, IDPL, IPBL,                &
                      ZUMF, ZUER, ZUDR, ZUTHL, ZUTHV, ZURW,            &
                      ZURC, ZURI, ZURR, ZURS, ZUPR,                    &
                      ZUTPR, ZCAPE, ICTL, IETL, ZUTT                    )
!
!
!
!*           4.2    In routine UPDRAFT GTRIG1 has been set to false when cloud 
!                   thickness is smaller than 3 km
!                   -----------------------------------------------------------
!
!
ICONV1 = COUNT(GTRIG1) 
!
IF ( ICONV1 > 0 )  THEN
!
!*       4.3    Allocate memory for downdraft variables
!               ---------------------------------------
!
! downdraft variables
!
  ALLOCATE( ILFS(ICONV) )
  ALLOCATE( IDBL(ICONV) )
  ALLOCATE( IML(ICONV) )
  ALLOCATE( ZDMF(ICONV,IKS) )
  ALLOCATE( ZDER(ICONV,IKS) )
  ALLOCATE( ZDDR(ICONV,IKS) )
  ALLOCATE( ZDTHL(ICONV,IKS) )
  ALLOCATE( ZDRW(ICONV,IKS) )
  ALLOCATE( ZLMASS(ICONV,IKS) ) ;  ZLMASS = 0.0
  DO JK = IKB, IKE
    ZLMASS(:,JK)  = ZDXDY(:) * ZDPRES(:,JK) / XG  ! mass of model layer
  END DO
  ZLMASS(:,IKB) = ZLMASS(:,IKB+1)
  ALLOCATE( ZMIXF(ICONV) )
  ALLOCATE( ZTPR(ICONV) )
  ALLOCATE( ZSPR(ICONV) )
  ALLOCATE( ZDTEVR(ICONV) )
  ALLOCATE( ZPREF(ICONV) )
  ALLOCATE( ZDTEVRF(ICONV,IKS) )
  ALLOCATE( ZPRLFLX(ICONV,IKS) )
  ALLOCATE( ZPRSFLX(ICONV,IKS) )
!
! closure variables
!
  ALLOCATE( ZTIMEA(ICONV) )
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
!*           5.1    Compute advective time period and precipitation 
!                   efficiency as a function of mean ambient wind (shear) 
!                   --------------------------------------------------------
!
  CALL CONVECT_TSTEP_PREF( ICONV, KLEV,                          &
                           ZU, ZV, ZPRES, ZZ, ZDXDY, ILCL, ICTL, &
                           ZTIMEA, ZPREF )
!
          ! exclude convective downdrafts if desired
  IF ( .NOT. ODOWN ) ZPREF(:) = 1.
!
! Compute the period during which convection is active
  ZTIMEC(:) = MAX( 1800., MIN( 3600., ZTIMEA(:) ) )
  ZTIMEC(:) = REAL( INT( ZTIMEC(:) / PDTCONV ) ) * PDTCONV
  ZTIMEC(:) = MAX( PDTCONV, ZTIMEC(:) ) ! necessary if PDTCONV > 1800
  IF ( OSETTADJ ) THEN
     ZTIMEC(:) = MAX( PDTCONV, ZTIMED(:) )
  END IF
!
!
!*           5.2    Compute melting level
!                   ----------------------
!
  IML(:) = IKB
  DO JK = IKE, IKB, -1
    WHERE( ZTT(:,JK) <= XTT )  IML(:) = JK
  END DO
!
  CALL CONVECT_DOWNDRAFT( ICONV, KLEV,                               &
                          KICE, ZPRES, ZDPRES, ZZ, ZTH, ZTHES,       & 
                          ZRW, ZRC, ZRI,                             &
                          ZPREF, ILCL, ICTL, IETL,                   &
                          ZUTHL, ZURW, ZURC, ZURI,                   &
                          ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,             &
                          ZMIXF, ZDTEVR, ILFS, IDBL, IML,            &
                          ZDTEVRF                                    )
!
!-------------------------------------------------------------------------------
!
!*           6.     Adjust up and downdraft mass flux to be consistent
!                   with precipitation efficiency relation.
!                   --------------------------------------------------- 
!
  CALL CONVECT_PRECIP_ADJUST( ICONV, KLEV,                              &
                              ZPRES,ZUMF, ZUER, ZUDR, ZUPR, ZUTPR, ZURW,&
                              ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,            &
                              ZPREF, ZTPR, ZMIXF, ZDTEVR,               &
                              ILFS, IDBL, ILCL, ICTL, IETL,             &
                              ZDTEVRF                                   )
!
!-------------------------------------------------------------------------------
!
!*           7.     Determine adjusted environmental values assuming
!                   that all available buoyant energy must be removed
!                   within an advective time step ZTIMEC.
!                   ---------------------------------------------------
!
  CALL CONVECT_CLOSURE( ICONV, KLEV,                                &
                        ZPRES, ZDPRES, ZZ, ZDXDY, ZLMASS,           &
                        ZTHL, ZTH, ZRW, ZRC, ZRI, GTRIG1,           &
                        ZTHC, ZRVC, ZRCC, ZRIC, ZWSUB,              &
                        ILCL, IDPL, IPBL, ILFS, ICTL, IML,          &
                        ZUMF, ZUER, ZUDR, ZUTHL, ZURW,              &
                        ZURC, ZURI, ZUPR,                           &
                        ZDMF, ZDER, ZDDR, ZDTHL, ZDRW,              &
                        ZTPR, ZSPR, ZDTEVR,                         &
                        ZCAPE, ZTIMEC,                              &
                        IFTSTEPS,                                   &
                        ZDTEVRF, ZPRLFLX, ZPRSFLX )
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
!     in order to save memory, the tendencies are temporarily stored
!     in the tables for the adjusted grid-scale values
!
  DO JK = IKB, IKE
    ZTHC(:,JK) = ( ZTHC(:,JK) - ZTH(:,JK) ) / ZTIMEC(:)             &
      * ( ZPRES(:,JK) / XP00 ) ** ZRDOCP ! change theta in temperature
    ZRVC(:,JK) = ( ZRVC(:,JK) - ZRW(:,JK) + ZRC(:,JK) + ZRI(:,JK) )/ ZTIMEC(:)  
    ZRCC(:,JK) = ( ZRCC(:,JK) - ZRC(:,JK) ) / ZTIMEC(:)
    ZRIC(:,JK) = ( ZRIC(:,JK) - ZRI(:,JK) ) / ZTIMEC(:) 
!
    ZPRLFLX(:,JK) = ZPRLFLX(:,JK) / ( XRHOLW * ZDXDY(:) )
    ZPRSFLX(:,JK) = ZPRSFLX(:,JK) / ( XRHOLW * ZDXDY(:) )
!
  END DO
!
  ZPRLFLX(:,IKB) = ZPRLFLX(:,IKB+1)
  ZPRSFLX(:,IKB) = ZPRSFLX(:,IKB+1)
!
!
!*           8.2    Apply conservation correction
!                   -----------------------------
!
          ! Compute vertical integrals
!
! Reproducibility
! JKM = MAXVAL( ICTL(:) )
  JKM = IKE - 1
  ZWORK2(:) = 0.
  ZWORK2B(:) = 0.
  DO JK = IKB+1, JKM
    JKP = JK + 1
    DO JI = 1, ICONV
      ZW1 = .5 * (ZPRES(JI,JK-1) - ZPRES(JI,JKP)) / XG
      ZWORK2(JI) = ZWORK2(JI) + ( ZRVC(JI,JK) + ZRCC(JI,JK) + ZRIC(JI,JK) ) * ZW1 ! moisture
      ZWORK2B(JI) = ZWORK2B(JI) + (    ( XCPD + XCPV * ZRW(JI,JK) )* ZTHC(JI,JK)   - &
                 ( XLVTT + ( XCPV - XCL ) * ( ZTT(JI,JK) - XTT ) ) * ZRCC(JI,JK)   - &
                 ( XLSTT + ( XCPV - XCL ) * ( ZTT(JI,JK) - XTT ) ) * ZRIC(JI,JK) ) * &
                                        ZW1                                       ! enthalpy
    END DO
  END DO
!
          ! Budget error (compare integral to surface precip.)
!
  DO JI = 1, ICONV
    IF ( ZTPR(JI) > 0.) THEN
      ZW1 = XG / ( ZPRES(JI,IKB) - ZPRES(JI,JKP) - .5 * ( &
                 ZDPRES(JI,IKB+1) - ZDPRES(JI,JKP+1) ) )
      ZWORK2(JI) = ( ZTPR(JI) / ZDXDY(JI) + ZWORK2(JI) ) * ZW1
      ZWORK2B(JI) = ( ZTPR(JI) / ZDXDY(JI) *                                &
         ( XLVTT + ( XCPV - XCL ) * ( ZTT(JI,IKB) - XTT ) ) - ZWORK2B(JI) ) &
                                       * ZW1
    END IF
  END DO
!
          ! Apply uniform correction
!
  DO JK = JKM, IKB+1, -1
  DO JI = 1, ICONV
    IF ( ZTPR(JI) > 0. .AND. JK <= ICTL(JI) ) THEN
      ! ZW1 = ABS(ZRVC(JI,JK)) +  ABS(ZRCC(JI,JK)) +  ABS(ZRIC(JI,JK)) + 1.E-12
      ! ZRVC(JI,JK) = ZRVC(JI,JK) - ABS(ZRVC(JI,JK))/ZW1*ZWORK2(JI)           ! moisture
      ZRVC(JI,JK) = ZRVC(JI,JK) - ZWORK2(JI)                                ! moisture
      ! ZRCC(JI,JK) = ZRCC(JI,JK) - ABS(ZRCC(JI,JK))/ZW1*ZWORK2(JI)
      ! ZRIC(JI,JK) = ZRIC(JI,JK) - ABS(ZRIC(JI,JK))/ZW1*ZWORK2(JI)
      ZTHC(JI,JK) = ZTHC(JI,JK) + ZWORK2B(JI) / ( XCPD + XCPV * ZRW(JI,JK) )! energy
    END IF
  END DO
  END DO
!
!
!     execute a "scatter"= pack command to store the tendencies in
!     the final 2D tables
!
  DO JK = IKB, IKE
  DO JI = 1, ICONV
    JL = IJINDEX(JI)
    PTTEN(JL,JK)   = ZTHC(JI,JK)
    PRVTEN(JL,JK)  = ZRVC(JI,JK)
    PRCTEN(JL,JK)  = ZRCC(JI,JK)
    PRITEN(JL,JK)  = ZRIC(JI,JK)
!
    PPRLFLX(JL,JK) = ZPRLFLX(JI,JK)
    PPRSFLX(JL,JK) = ZPRSFLX(JI,JK)
  END DO
  END DO
!
!
!*           8.3    Convective rainfall tendency
!                   ----------------------------
!
         ! liquid and solid surface rainfall tendency in m/s
  ZTPR(:)   = ZTPR(:) / ( XRHOLW * ZDXDY(:) ) ! total surf precip
  ZSPR(:)   = ZSPR(:) / ( XRHOLW * ZDXDY(:) ) ! solid surf precip
  ZTPR(:)   = ZTPR(:) - ZSPR(:) ! compute liquid part
!
  DO JI = 1, ICONV
    JL = IJINDEX(JI)
    PPRLTEN(JL) = ZTPR(JI)
    PPRSTEN(JL) = ZSPR(JI)
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
!*           8.4    Set convective counter
!                   ----------------------
!
         ! compute convective counter for just activated convective
         ! grid points
         ! If the advective time period is less than specified
         ! minimum for convective period, allow feedback to occur only
         ! during advective time
!
  ZTIME(:) = 1.
  ZWORK2(:) = 0.
  DO JI = 1, ICONV
    JL = IJINDEX(JI)
    ZTIME(JL)  =  ZTIMEC(JI)
    ZWORK2(JL) =  ZTIMEA(JI)
    ZWORK2(JL) =  MIN( ZWORK2(JL), ZTIME(JL) )
    ZWORK2(JL) =  MAX( ZWORK2(JL), PDTCONV )
    IF ( GTRIG(JL) )  KCOUNT(JL) = INT( ZWORK2(JL) / PDTCONV )
    IF ( GTRIG(JL) .AND. PPRLTEN(JL)<1.E-14 ) KCOUNT(JL) = 0
  END DO
!
!
!*           8.7    Compute convective tendencies for Tracers
!                   ------------------------------------------
!
  IF ( OCH1CONV ) THEN
!
    ALLOCATE( ZCH1(ICONV,IKS,KCH1) )    ; ZCH1 = 0.0
    ALLOCATE( ZCH1C(ICONV,IKS,KCH1) )   ; ZCH1C = 0.0
    ALLOCATE( ZWORK3(ICONV,KCH1) )
!
    ALLOCATE( ZRHODREF(ICONV,IKS) )
    ZRHODREF=0.
    IF ( OCH_CONV_LINOX ) THEN
      ALLOCATE( ZZZ(ICONV,IKS) )
      ALLOCATE( ZIC_RATE(ICONV) )
      ALLOCATE( ZCG_RATE(ICONV) )
      ALLOCATE( ZWORK4(ICONV,IKS) )
      ALLOCATE( ZWORK4C(ICONV,IKS) )
      ZZZ=0.
      ZIC_RATE=0.
      ZCG_RATE=0.
      ZWORK4=0.
      ZWORK4C=0.
    END IF
!
    DO JI = 1, ICONV
     DO JK = IKB, IKE
      JL = IJINDEX(JI)
      ZCH1(JI,JK,:) = PCH1(JL,JK,:)
      ZRHODREF(JI,JK)=PRHODREF(JL,JK)
     END DO
     ZRHODREF(JI,1) = PRHODREF(JL,IKB) 
     ZRHODREF(JI,IKS) = PRHODREF(JL,IKE) 
    END DO
    ZCH1(:,1,:) = ZCH1(:,IKB,:) 
    ZCH1(:,IKS,:) = ZCH1(:,IKE,:) 
!
    JN_NO = 0
    IF ( OCH_CONV_LINOX ) THEN
      DO JK = IKB, IKE
      DO JI = 1, ICONV
        JL = IJINDEX(JI)
        ZZZ(JI,JK)=PZZ(JL,JK)
        ZIC_RATE(JI)=PIC_RATE(JL)
        ZCG_RATE(JI)=PCG_RATE(JL)
      END DO
      END DO
      IF (OUSECHEM) THEN
        DO JN = NSV_CHEMBEG,NSV_CHEMEND
          IF (CNAMES(JN-NSV_CHEMBEG+1)=='NO') JN_NO = JN
        END DO
      ELSE
        JN_NO = NSV_LNOXBEG
      ENDIF
      ZWORK4(:,:) = ZCH1(:,:,JN_NO)
      CALL CH_CONVECT_LINOX( ICONV, KLEV, ZWORK4, ZWORK4C,        &
                             IDPL, IPBL, ILCL, ICTL, ILFS, IDBL,  &
                             ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,  &
                             ZTIMEC, ZDXDY, ZMIXF, ZLMASS, ZWSUB, &
                             IFTSTEPS, ZUTT, ZRHODREF,            &
                             OUSECHEM, ZZZ, ZIC_RATE, ZCG_RATE    )
      DO JI = 1, ICONV
        JL = IJINDEX(JI)
        PIC_RATE(JL)=ZIC_RATE(JI)
        PCG_RATE(JL)=ZCG_RATE(JI)
      ENDDO
    ENDIF
!
    IF ((OUSECHEM .AND. OCH_CONV_SCAV).OR.(ODUST .AND.  OCH_CONV_SCAV).OR.&
        (OSALT .AND.  OCH_CONV_SCAV)  ) THEN
! 
      CALL CH_CONVECT_SCAVENGING( ICONV, KLEV, KCH1, ZCH1, ZCH1C,      &
                                  IDPL, IPBL, ILCL, ICTL, ILFS, IDBL,  &
                                  ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,  &
                                  ZTIMEC, ZDXDY, ZMIXF, ZLMASS, ZWSUB, &
                                  IFTSTEPS,                            &
                                  ZURC, ZURR, ZURI, ZURS, ZUTT, ZPRES, &
                                  ZRHODREF, PPABST, ZTHT               )
!
      IF (OCH_CONV_LINOX) THEN
        ZCH1C(:,:,JN_NO) = ZWORK4C(:,:)
      ENDIF
!    no conservation correction for scavenging
     DO JI = 1, ICONV
        JL = IJINDEX(JI)
        IF ( ZTPR(JI) > 0. ) THEN
          DO JK = IKB, IKE
            PCH1TEN(JL,JK,:) = (ZCH1C(JI,JK,:)- ZCH1(JI,JK,:)) /ZTIMEC(JI)
          END DO
        ELSE
          DO JK = IKB, IKE
            PCH1TEN(JL,JK,:) = 0.
          END DO
        ENDIF
      END DO

!
    ELSE
!
      CALL CONVECT_CHEM_TRANSPORT( ICONV, KLEV, KCH1, ZCH1, ZCH1C,      &
                                   IDPL, IPBL, ILCL, ICTL, ILFS, IDBL,  &
                                   ZUMF, ZUER, ZUDR, ZDMF, ZDER, ZDDR,  &
                                   ZTIMEC, ZDXDY, ZMIXF, ZLMASS, ZWSUB, &
                                   IFTSTEPS )
!
      IF (OCH_CONV_LINOX) THEN
        ZCH1C(:,:,JN_NO) = ZWORK4C(:,:)
      ENDIF
!
!*           8.8    Apply conservation correction
!                   -----------------------------
!
          ! Compute vertical integrals
!
! Reproducibility
!     JKM = MAXVAL( ICTL(:) )
      JKM = IKE - 1
      DO JN = 1, KCH1
        IF((JN < NSV_LGBEG .OR. JN>NSV_LGEND-1) .AND. JN .NE. JN_NO ) THEN 
          ! no correction for Lagrangian and LiNOx variables
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
            IF ( ZTPR(JI) > 0. .AND. JK <= ICTL(JI) ) THEN
              ZCH1C(JI,JK,JN) = ZCH1C(JI,JK,JN) -       &
              ZWORK3(JI,JN)*ABS(ZCH1C(JI,JK,JN))/MAX(1.E-30,ZWORK2(JI))
                ! ZCH1C(JI,JK,JN) = MAX( ZCH1C(JI,JK,JN), -ZCH1(JI,JK,JN)/ZTIMEC(JI) )
            END IF
          END DO
          END DO
        END IF
!
        DO JI = 1, ICONV
          JL = IJINDEX(JI)
          IF ( ZTPR(JI) > 0. ) THEN
            DO JK = IKB, IKE
              PCH1TEN(JL,JK,JN) = (ZCH1C(JI,JK,JN)-ZCH1(JI,JK,JN) ) /ZTIMEC(JI)
            END DO
          ELSE
            DO JK = IKB, IKE
              PCH1TEN(JL,JK,JN) = 0.
            END DO
          ENDIF
        END DO
      END DO
    END IF
  END IF
!
!-------------------------------------------------------------------------------
!
!*           9.     Write up- and downdraft mass fluxes 
!                   ------------------------------------
!
  DO JK = IKB, IKE
    ZUMF(:,JK)  = ZUMF(:,JK) / ZDXDY(:) ! Mass flux per unit area
    ZDMF(:,JK)  = ZDMF(:,JK) / ZDXDY(:)
  END DO
  ZWORK2(:) = 1.
  WHERE ( PPRLTEN(:)<1.E-14 ) ZWORK2(:) = 0.
  DO JK = IKB, IKE
  DO JI = 1, ICONV
    JL = IJINDEX(JI)
    PUMF(JL,JK) = ZUMF(JI,JK) * ZWORK2(JL)
    PDMF(JL,JK) = ZDMF(JI,JK) * ZWORK2(JL)
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
  DEALLOCATE( ZDTHL )
  DEALLOCATE( ZDRW )
  DEALLOCATE( ZLMASS )
  DEALLOCATE( ZMIXF )
  DEALLOCATE( ZTPR )
  DEALLOCATE( ZSPR )
  DEALLOCATE( ZDTEVR )
  DEALLOCATE( ZPREF )
  DEALLOCATE( IML )
  DEALLOCATE( ILFS )
  DEALLOCATE( IDBL )
  DEALLOCATE( ZDTEVRF )
  DEALLOCATE( ZPRLFLX )
  DEALLOCATE( ZPRSFLX )
!
!   closure variables
!
  DEALLOCATE( ZTIMEA )
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
    DEALLOCATE( ZRHODREF )
    IF ( OCH_CONV_LINOX ) THEN
      DEALLOCATE( ZZZ )
      DEALLOCATE( ZIC_RATE )
      DEALLOCATE( ZCG_RATE )
      DEALLOCATE( ZWORK4 )
      DEALLOCATE( ZWORK4C )
    END IF
  END IF
!
ENDIF
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
DEALLOCATE( ZU )
DEALLOCATE( ZV )
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
DEALLOCATE( ZUTT )
DEALLOCATE( ZURW )
DEALLOCATE( ZURC )
DEALLOCATE( ZURI )
DEALLOCATE( ZURR )
DEALLOCATE( ZURS )
DEALLOCATE( ZUPR )
DEALLOCATE( ZUTPR )
DEALLOCATE( ZTHLCL )
DEALLOCATE( ZTLCL )
DEALLOCATE( ZRVLCL )
DEALLOCATE( ZWLCL )
DEALLOCATE( ZZLCL )
DEALLOCATE( ZTHVELCL )
DEALLOCATE( ZMFLCL )
DEALLOCATE( ZCAPE )
IF ( OSETTADJ ) DEALLOCATE( ZTIMED )
!
! work arrays
!
DEALLOCATE( IINDEX )
DEALLOCATE( IJINDEX )
DEALLOCATE( IJSINDEX )
DEALLOCATE( GTRIG1 )
!
!
END SUBROUTINE DEEP_CONVECTION
