!MNH_LIC Copyright 2008-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ####################################
       MODULE MODI_CH_AQUEOUS_TMICICE
!      ####################################
!
INTERFACE
      SUBROUTINE CH_AQUEOUS_TMICICE( PTSTEP, PRHODREF, PRHODJ, PTHT, PPABST,      &
                                     PRTMIN_AQ, OUSECHIC, OCH_RET_ICE, HNAMES,    &
                                     HICNAMES, KEQ, KEQAQ, PRVT, PRCT, PRRT, PRIT,&
                                     PRST, PRGT, PCIT, PRCS, PRRS, PRIS, PRSS,    &
                                     PRGS, PGSVT, PGRSVS, PCSVT, PCRSVS, PRSVT,   &
                                     PRRSVS, PSGSVT, PSGRSVS                      )   
!
REAL,                     INTENT(IN)    :: PTSTEP    ! Time step          
REAL,                     INTENT(IN)    :: PRTMIN_AQ ! LWC threshold liq. chem.
INTEGER,                  INTENT(IN)    :: KEQ   ! Number of chem. spec.
INTEGER,                  INTENT(IN)    :: KEQAQ   ! Number of liq. chem. spec.
LOGICAL,                  INTENT(IN)    :: OUSECHIC ! flag for ice chem.
LOGICAL,                  INTENT(IN)    :: OCH_RET_ICE ! flag for retention in ice
!
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES ! name of chem. species
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HICNAMES ! name of ice chem. species
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rainwater m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Pristine conc. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCS    ! cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rainwater m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIS    ! Pristine m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRSS    ! Snow m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGS    ! graupel m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PGSVT   ! gas species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PGRSVS  ! gas species source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PCSVT   ! cloud water aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PCRSVS  ! cloud water aq. species source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRSVT   ! Rainwater aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRSVS  ! Rainwater aq. species source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSGSVT  ! ice species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSGRSVS ! ice species source
!
END SUBROUTINE CH_AQUEOUS_TMICICE
END INTERFACE
END MODULE MODI_CH_AQUEOUS_TMICICE
!
!     ################################################################################
      SUBROUTINE CH_AQUEOUS_TMICICE( PTSTEP, PRHODREF, PRHODJ, PTHT, PPABST,      &
                                     PRTMIN_AQ, OUSECHIC, OCH_RET_ICE, HNAMES,    &
                                     HICNAMES, KEQ, KEQAQ, PRVT, PRCT, PRRT, PRIT,&
                                     PRST, PRGT, PCIT, PRCS, PRRS, PRIS, PRSS,    &
                                     PRGS, PGSVT, PGRSVS, PCSVT, PCRSVS, PRSVT,   &
                                     PRRSVS, PSGSVT, PSGRSVS                      )   
!     ################################################################################
!
!!****  * -  compute the explicit microphysical sources 
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the microphysical sources
!!    corresponding to collision/coalescence processes (autoconversion + accretion)
!!    and to the freezing, rimin and melting processes for snow and graupel
!!    for the ICE3(4) cloud microphysics parameterization (see rain_ice)
!!
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      None
!!     
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS
!!          JPHEXT       : Horizontal external points number
!!          JPVEXT       : Vertical external points number
!!
!!    REFERENCE
!!    ---------
!!      Book1 of the documentation ( routine CH_AQUEOUS_TMICICE )
!!
!!    AUTHOR
!!    ------
!!      C. Mari J.P. Pinty M. Leriche      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/08
!!    M. Leriche 19/07/2010 add riming, freezing and melting for ice phase(ICE3)
!!    M. Leriche 17/09/2010 add OUSECHIC flag
!!    Juan 24/09/2012: for BUG Pgi rewrite PACK function on mode_pack_pgi
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!    M.Leriche 2015 correction bug
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS,      ONLY : JPHEXT,    &! number of horizontal External points
                                 JPVEXT      ! number of vertical External points
USE MODD_CST,             ONLY : XP00, XRD, XRV, XCPD, XTT, XLMTT, XLVTT, XCPV, &
                                 XCL, XCI, XESTT, XMV, XMD
USE MODD_RAIN_ICE_DESCR_n,  ONLY : XLBR, XLBEXR, XCEXVT, XLBDAS_MAX, XLBS, XLBEXS, &
                                 XLBG, XLBEXG, XCXS, XCXG, XDG, XBS
USE MODD_RAIN_ICE_PARAM_n,  ONLY : XTIMAUTC, XCRIAUTC, XFCACCR, XEXCACCR, &
                                 XRIMINTP1, XRIMINTP2, XCRIMSS, XCRIMSG,&
                                 XEXCRIMSS, XEXCRIMSG, NGAMINC, XGAMINC_RIM1, &
                                 XFRACCSS, XLBRACCS1, XLBRACCS2, XLBRACCS3,    &
                                 XACCINTP1S, XACCINTP2S, XACCINTP1R, XACCINTP2R, &
                                 NACCLBDAS, NACCLBDAR, XKER_RACCSS, XKER_RACCS, &
                                 XEXRCFRI, XRCFRI, X0DEPG, XEX0DEPG, X1DEPG,    &
                                 XEX1DEPG, XSCFAC, XFCDRYG, XFIDRYG, XCOLEXIG,  &
                                 XCOLEXSG, XFSDRYG, NDRYLBDAG, XDRYINTP1G,      &
                                 XDRYINTP2G, NDRYLBDAS, XDRYINTP1S, XDRYINTP2S, &
                                 XKER_SDRYG, XLBSDRYG1, XLBSDRYG2, XLBSDRYG3,   &
                                 XFRDRYG, NDRYLBDAR, XDRYINTP1R, XDRYINTP2R,    &
                                 XKER_RDRYG, XLBRDRYG1, XLBRDRYG2, XLBRDRYG3,   &
                                 XCOLIG, XCOLEXIG, XCOLSG, XCOLEXSG
USE MODD_CH_ICE                              ! value of retention coefficient
USE MODD_CH_ICE_n                            ! index for ice phase chemistry with IC3/4
!
#ifdef MNH_PGI
USE MODE_PACK_PGI
#endif
use mode_tools,           only: Countjv
use mode_tools_ll,        only: GET_INDICE_ll
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,                     INTENT(IN)    :: PTSTEP    ! Time step          
REAL,                     INTENT(IN)    :: PRTMIN_AQ ! LWC threshold liq. chem.
INTEGER,                  INTENT(IN)    :: KEQ   ! Number of chem. spec.
INTEGER,                  INTENT(IN)    :: KEQAQ   ! Number of liq. chem. spec.
LOGICAL,                  INTENT(IN)    :: OUSECHIC ! flag for ice chem.
LOGICAL,                  INTENT(IN)    :: OCH_RET_ICE ! flag for retention in ice
!
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HNAMES ! name of chem. species
CHARACTER(LEN=32), DIMENSION(:), INTENT(IN) :: HICNAMES ! name of ice chem. species
!
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT    ! Theta at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRVT    ! Vapor m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCT    ! cloud water m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRT    ! Rainwater m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIT    ! Pristine m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRST    ! Snow m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGT    ! Graupel m.r. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PCIT    ! Pristine conc. at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRCS    ! cloud water m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRRS    ! Rainwater m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRIS    ! Pristine m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRSS    ! Snow m.r. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRGS    ! graupel m.r. source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PGSVT   ! gas species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PGRSVS  ! gas species source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PCSVT   ! cloud water aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PCRSVS  ! cloud water aq. species source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRSVT   ! Rainwater aq. species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRSVS  ! Rainwater aq. species source
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PSGSVT  ! ice species at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PSGRSVS ! ice species source
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JLC, JLR, JLI, JLG, JLW ! Loop index for cloud water, rainwater and ice species
INTEGER :: JJ                 ! Loop index 
INTEGER :: IIB                !  Define the domain where is 
INTEGER :: IIE                !  the microphysical sources have to be computed
INTEGER :: IJB            
INTEGER :: IJE           
INTEGER :: IKB           
INTEGER :: IKE           
!
INTEGER :: IMICRO        ! case number of r_x>0 locations
LOGICAL, DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: GMICRO   ! where to compute mic. processes
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZT     ! Temperature
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZRCS   ! Cloud water m.r. source phys.tendency
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZRRS   ! Rain water m.r. source phys. tendency
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZRIS   ! Pristine m.r. source phys. tendency
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZRSS   ! Snow m.r. source phys. tendency
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZRGS   ! Graupel m.r. source phys. tendency
REAL,    DIMENSION(SIZE(PGRSVS,1),SIZE(PGRSVS,2),SIZE(PGRSVS,3),SIZE(PGRSVS,4))   &
                                :: ZZGRSVS   ! Gas species source
REAL,    DIMENSION(SIZE(PCRSVS,1),SIZE(PCRSVS,2),SIZE(PCRSVS,3),SIZE(PCRSVS,4))   &
                                :: ZZCRSVS   ! Cloud water aq. species source
REAL,    DIMENSION(SIZE(PRRSVS,1),SIZE(PRRSVS,2),SIZE(PRRSVS,3),SIZE(PRRSVS,4))   &
                                :: ZZRRSVS   ! Rain water aq. species source
REAL,    DIMENSION(SIZE(PSGRSVS,1),SIZE(PSGRSVS,2),SIZE(PSGRSVS,3),SIZE(PSGRSVS,4))   &
                                :: ZZSGRSVS   !  Ice (snow+graupel) species source
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZCW      ! work array
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZRW      ! work array
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZSGW      ! work array
REAL,    DIMENSION(SIZE(PRCT,1),SIZE(PRCT,2),SIZE(PRCT,3))   &
                                :: ZGW      ! work array
REAL, DIMENSION(:),   ALLOCATABLE :: ZZT    ! Temperature
REAL, DIMENSION(:),   ALLOCATABLE :: ZPRES  ! Pressure
REAL, DIMENSION(:),   ALLOCATABLE :: ZRVT   ! Vapor m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRCT   ! Cloud water m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRRT   ! Rain water m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRIT   ! Pristine m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRST   ! Snow m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZRGT   ! Graupel m.r. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZCIT   ! Pristine conc. at t
REAL, DIMENSION(:),   ALLOCATABLE :: ZZRCS  ! Cloud water m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZZRRS  ! Rain water m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZZRIS  ! Pristine m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZZRSS  ! snow m.r. source
REAL, DIMENSION(:),   ALLOCATABLE :: ZZRGS  ! graupel m.r. source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCSVT  ! Cloud water aq. species at t
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRSVT  ! Rain water aq. species at t
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSGSVT ! Ice (snow + graupel) species at t
REAL, DIMENSION(:,:), ALLOCATABLE :: ZGRSVS ! Gas species source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCRSVS ! Cloud water aq. species source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZRRSVS ! Rain water aq. species source
REAL, DIMENSION(:,:), ALLOCATABLE :: ZSGRSVS! Ice (snow+graupel) species source
REAL, DIMENSION(:),   ALLOCATABLE :: ZCJ    ! Function to compute the ventilation coefficient
REAL, DIMENSION(:),   ALLOCATABLE :: ZKA    ! Thermal conductivity of the air
REAL, DIMENSION(:),   ALLOCATABLE :: ZDV    ! Diffusivity of water vapor in the air
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZRHODREF, & ! RHO Dry REFerence
                                     ZZW,      & ! Work array
                                     ZLBDAR,   & ! Slope parameter of the raindrop distribution
                                     ZLBDAS,   & ! Slope parameter of the snow distribution
                                     ZLBDAG,   & ! Slope parameter of the graupel distribution
                                     ZRDRYG,   & ! Dry growth rate of the graupel
                                     ZRWETG      ! Wet growth rate of the graupel
!
INTEGER :: IGRIM, IGACC                    ! Case number of riming, accretion
INTEGER :: IGDRY
!, IGWET                    ! dry growth and wet growth locations for graupels
LOGICAL, DIMENSION(:), ALLOCATABLE :: GRIM ! Test where to compute riming
LOGICAL, DIMENSION(:), ALLOCATABLE :: GACC ! Test where to compute accretion
LOGICAL, DIMENSION(:), ALLOCATABLE :: GDRY ! Test where to compute dry growth
!LOGICAL, DIMENSION(:), ALLOCATABLE :: GWET  ! Test where to compute wet growt
INTEGER, DIMENSION(:), ALLOCATABLE :: IVEC1,IVEC2       ! Vectors of indices for
                                                        ! interpolations
REAL,    DIMENSION(:), ALLOCATABLE :: ZVEC1,ZVEC2,ZVEC3 ! Work vectors for
                                                        ! interpolations
REAL, DIMENSION(:,:), ALLOCATABLE :: ZZW1, ZZW2, ZZW3, ZZW4   ! Work arrays
!
INTEGER , DIMENSION(SIZE(GMICRO)) :: I1,I2,I3 ! Used to replace the COUNT
INTEGER                           :: JL       ! and PACK intrinsics
!
!
!  compute the temperature
!
ZT(:,:,:) = PTHT(:,:,:) * ( PPABST(:,:,:) / XP00 ) ** (XRD/XCPD)
!  
!-------------------------------------------------------------------------------
!
!*       1.     COMPUTE THE LOOP BOUNDS
!   	        -----------------------
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=SIZE(PRCT,3) - JPVEXT
!
!-------------------------------------------------------------------------------
!
!!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!                ---------------------------------------
!
ZRCS(:,:,:)  = PRCS(:,:,:) / PRHODJ(:,:,:)
ZRRS(:,:,:)  = PRRS(:,:,:) / PRHODJ(:,:,:)
ZRSS(:,:,:)  = PRSS(:,:,:) / PRHODJ(:,:,:)
ZRIS(:,:,:)  = PRIS(:,:,:) / PRHODJ(:,:,:)
ZRGS(:,:,:)  = PRGS(:,:,:) / PRHODJ(:,:,:)
!
DO JLC= 1, SIZE(PCRSVS,4)
  ZZCRSVS(:,:,:,JLC) = PCRSVS(:,:,:,JLC) / PRHODJ(:,:,:)
ENDDO
DO JLR= 1, SIZE(PRRSVS,4)
  ZZRRSVS(:,:,:,JLR) = PRRSVS(:,:,:,JLR) / PRHODJ(:,:,:)
ENDDO
IF (OUSECHIC) THEN
  DO JLG= 1, SIZE(PGRSVS,4)
    ZZGRSVS(:,:,:,JLG) = PGRSVS(:,:,:,JLG) / PRHODJ(:,:,:)
  ENDDO
  DO JLI= 1, SIZE(PSGRSVS,4)
    ZZSGRSVS(:,:,:,JLI) = PSGRSVS(:,:,:,JLI) / PRHODJ(:,:,:)
  ENDDO
ELSE
  IF (.NOT.(OCH_RET_ICE)) THEN
    DO JLG= 1, SIZE(PGRSVS,4)
      ZZGRSVS(:,:,:,JLG) = PGRSVS(:,:,:,JLG) / PRHODJ(:,:,:)
    ENDDO
  ENDIF
ENDIF
!
!-------------------------------------------------------------------------------
!
!*       3.     OPTIMIZATION: looking for locations where m.r. hydro. > min value
!   	        -----------------------------------------------------------------
!
GMICRO(:,:,:) = .FALSE.
GMICRO(IIB:IIE,IJB:IJE,IKB:IKE) =                                                           &
      (PRCT(IIB:IIE,IJB:IJE,IKB:IKE)>PRTMIN_AQ*1.e3/PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE)) .OR. &
      (PRRT(IIB:IIE,IJB:IJE,IKB:IKE)>PRTMIN_AQ*1.e3/PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE)) .OR. &
      (PRST(IIB:IIE,IJB:IJE,IKB:IKE)>PRTMIN_AQ*1.e3/PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE)) .OR. &
      (PRGT(IIB:IIE,IJB:IJE,IKB:IKE)>PRTMIN_AQ*1.e3/PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE))
!
IMICRO = COUNTJV( GMICRO(:,:,:),I1(:),I2(:),I3(:))
IF( IMICRO >= 1 ) THEN
  ALLOCATE(ZZT(IMICRO))
  ALLOCATE(ZPRES(IMICRO))
  ALLOCATE(ZRVT(IMICRO))
  ALLOCATE(ZRCT(IMICRO))
  ALLOCATE(ZRRT(IMICRO))  
  ALLOCATE(ZRIT(IMICRO))  
  ALLOCATE(ZRST(IMICRO))  
  ALLOCATE(ZRGT(IMICRO))  
  ALLOCATE(ZCIT(IMICRO))  
  ALLOCATE(ZCSVT(IMICRO,SIZE(PCSVT,4)))
  ALLOCATE(ZRSVT(IMICRO,SIZE(PRSVT,4)))
  ALLOCATE(ZZRCS(IMICRO)) 
  ALLOCATE(ZZRRS(IMICRO)) 
  ALLOCATE(ZZRIS(IMICRO)) 
  ALLOCATE(ZZRSS(IMICRO)) 
  ALLOCATE(ZZRGS(IMICRO)) 
  ALLOCATE(ZCRSVS(IMICRO,SIZE(PCRSVS,4)))
  ALLOCATE(ZRRSVS(IMICRO,SIZE(PRRSVS,4)))
  ALLOCATE(ZRHODREF(IMICRO)) 
  ALLOCATE(ZZW(IMICRO))
  ALLOCATE(ZZW2(IMICRO,SIZE(PCSVT,4)))
  ALLOCATE(ZZW4(IMICRO,SIZE(PCSVT,4)))
  ALLOCATE(ZZW1(IMICRO,6))
  ALLOCATE(ZLBDAR(IMICRO))
  ALLOCATE(ZLBDAS(IMICRO))
  ALLOCATE(ZLBDAG(IMICRO))
  ALLOCATE(ZRDRYG(IMICRO))
  ALLOCATE(ZRWETG(IMICRO))
  ALLOCATE(ZKA(IMICRO))
  ALLOCATE(ZDV(IMICRO))
  ALLOCATE(ZCJ(IMICRO))
  DO JL=1,IMICRO   
    ZCSVT(JL,:) = PCSVT(I1(JL),I2(JL),I3(JL),:)
    ZCRSVS(JL,:) = ZZCRSVS(I1(JL),I2(JL),I3(JL),:)
    ZRSVT(JL,:) = PRSVT(I1(JL),I2(JL),I3(JL),:)
    ZRRSVS(JL,:) = ZZRRSVS(I1(JL),I2(JL),I3(JL),:)
!
    ZRVT(JL) = PRVT(I1(JL),I2(JL),I3(JL))
    ZRCT(JL) = PRCT(I1(JL),I2(JL),I3(JL))
    ZRRT(JL) = PRRT(I1(JL),I2(JL),I3(JL))
    ZRIT(JL) = PRIT(I1(JL),I2(JL),I3(JL))
    ZRST(JL) = PRST(I1(JL),I2(JL),I3(JL))
    ZRGT(JL) = PRGT(I1(JL),I2(JL),I3(JL))
    ZCIT(JL) = PCIT(I1(JL),I2(JL),I3(JL))
!
    ZZRCS(JL) = ZRCS(I1(JL),I2(JL),I3(JL))
    ZZRRS(JL) = ZRRS(I1(JL),I2(JL),I3(JL))
    ZZRIS(JL) = ZRIS(I1(JL),I2(JL),I3(JL))
    ZZRSS(JL) = ZRSS(I1(JL),I2(JL),I3(JL))
    ZZRGS(JL) = ZRGS(I1(JL),I2(JL),I3(JL))
!
    ZZT(JL) = ZT(I1(JL),I2(JL),I3(JL))
    ZRHODREF(JL) = PRHODREF(I1(JL),I2(JL),I3(JL))
    ZPRES(JL) = PPABST(I1(JL),I2(JL),I3(JL))
  ENDDO
  IF (OUSECHIC) THEN
    ALLOCATE(ZSGSVT(IMICRO,SIZE(PSGSVT,4)))
    ALLOCATE(ZGRSVS(IMICRO,SIZE(PGRSVS,4)))
    ALLOCATE(ZSGRSVS(IMICRO,SIZE(PSGRSVS,4)))
    ALLOCATE(ZZW3(IMICRO,SIZE(PSGSVT,4)))
    DO JL=1,IMICRO   
      ZGRSVS(JL,:) = ZZGRSVS(I1(JL),I2(JL),I3(JL),:)
      ZSGSVT(JL,:) = PSGSVT(I1(JL),I2(JL),I3(JL),:)
      ZSGRSVS(JL,:) = ZZSGRSVS(I1(JL),I2(JL),I3(JL),:)
    ENDDO
  ELSE
    IF (.NOT.(OCH_RET_ICE)) THEN
      ALLOCATE(ZGRSVS(IMICRO,SIZE(PGRSVS,4)))
      DO JL=1,IMICRO   
        ZGRSVS(JL,:) = ZZGRSVS(I1(JL),I2(JL),I3(JL),:)
      ENDDO
    ENDIF
  ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*       4.     COMPUTES THE SLOW WARM PROCESS SOURCES
!   	        --------------------------------------
!
!*       4.1    compute the slope parameter Lbda_r
!
  WHERE( ZRRT(:)>0.0 )
    ZLBDAR(:)  = XLBR*( ZRHODREF(:)*MAX( ZRRT(:),PRTMIN_AQ*1.e3/ZRHODREF(:)) )**XLBEXR
  END WHERE
!
!*       4.2    compute the autoconversion of r_c for r_r production: RCAUTR
!
    ZZW(:) = 0.0
    ZZW2(:,:) = 0.0
!
  DO JL=1,IMICRO
    IF ( (ZRCT(JL)>0.0) .AND. (ZZRCS(JL)>0.0) ) THEN
      ZZW(JL) = MIN( ZZRCS(JL),XTIMAUTC*MAX( ZRCT(JL)-XCRIAUTC/ZRHODREF(JL),0.0))
!
      ZZW2(JL,:) = ZZW(JL) * ZCSVT(JL,:)/ZRCT(JL)
      ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZCSVT(JL,:)/PTSTEP)),0.0)
      ZCRSVS(JL,:)  = ZCRSVS(JL,:) - ZZW2(JL,:) 
      ZRRSVS(JL,:)  = ZRRSVS(JL,:) + ZZW2(JL,:) 
    END IF
  END DO 
!
!*       4.3    compute the accretion of r_c for r_r production: RCACCR
!
    ZZW(:) = 0.0
    ZZW2(:,:) = 0.0
!    
  DO JL = 1,IMICRO
    IF( (ZRCT(JL)>0.0) .AND. (ZRRT(JL)>0.0) .AND. (ZZRCS(JL)>0.0) ) THEN
      ZZW(JL) = MIN( ZZRCS(JL),XFCACCR * ZRCT(JL)                &
                                    * ZLBDAR(JL)**XEXCACCR      &
                                    * ZRHODREF(JL)**(-XCEXVT) )
!
      ZZW2(JL,:) = ZZW(JL) * ZCSVT(JL,:)/ZRCT(JL)
      ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZCSVT(JL,:)/PTSTEP)),0.0)
      ZCRSVS(JL,:)  = ZCRSVS(JL,:) - ZZW2(JL,:) 
      ZRRSVS(JL,:)  = ZRRSVS(JL,:) + ZZW2(JL,:)
    END IF
  END DO
!
!
!*       4.4    compute the evaporation of r_r: RREVAV
!
! calculated by the kinetic mass transfer equation (BASIC.f90)
!
!
!-------------------------------------------------------------------------------
!
!*       5.     COMPUTES THE SLOW COLD PROCESS SOURCES
!               --------------------------------------
!
!*       5.1    compute the spontaneous freezing source: RRHONG
!
 ZZW(:) = 0.0
 ZZW2(:,:) = 0.0
! 
 DO JL = 1,IMICRO
   IF( (ZZT(JL)<XTT-35.0) .AND. (ZRRT(JL)>0.) .AND. (ZZRRS(JL)>0.) ) THEN
     ZZW(JL) = MIN( ZZRRS(JL),ZRRT(JL)/PTSTEP )
     ZZW2(JL,:) = ZZW(JL) * ZRSVT(JL,:)/ZRRT(JL)
     ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZRSVT(JL,:)/PTSTEP)),0.0)
     ZRRSVS(JL,:) = ZRRSVS(JL,:) - ZZW2(JL,:)
     IF (OUSECHIC) THEN
     DO JLI = 1, SIZE(PSGRSVS,4)
       IF (TRIM(HICNAMES(JLI)) == 'IC_HNO3' .OR. TRIM(HICNAMES(JLI)) == 'IC_SULF' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_H2SO4' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_NH3' .OR. TRIM(HICNAMES(JLI)) == 'IC_HCL' &
          .OR. HICNAMES(JLI)(1:4) == 'IC_A' .OR. HICNAMES(JLI)(1:4) == 'IC_B' &
          .OR. NINDEXGI(JLI).EQ.0) THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETNA * ZZW2(JL,NINDEXWI(JLI))           
       ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_H2O2' .OR. TRIM(HICNAMES(JLI)) == 'IC_HO2' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_HONO' .OR. TRIM(HICNAMES(JLI)) == 'IC_HNO4'&
          .OR. TRIM(HICNAMES(JLI)) == 'IC_HCHO' .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA1'&
          .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA2') THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETHP * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETHP) * ZZW2(JL,NINDEXWI(JLI))  
       ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_SO2' .OR. TRIM(HICNAMES(JLI)) == 'IC_OH' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_MO2' .OR. &
               TRIM(HICNAMES(JLI)) == 'IC_OP1') THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETSU * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETSU) * ZZW2(JL,NINDEXWI(JLI))
       ELSE
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETDF * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETDF) * ZZW2(JL,NINDEXWI(JLI))
       ENDIF
     ENDDO
     ELSE
       IF (.NOT.(OCH_RET_ICE)) THEN
         DO JLW = 1, SIZE(PRRSVS,4)
            IF (.NOT.(NINDEXWG(JLW).EQ.0)) THEN
              ZGRSVS(JL,NINDEXWG(JLW)) = ZGRSVS(JL,NINDEXWG(JLW)) + ZZW2(JL,JLW)
            ENDIF
         ENDDO
       ENDIF
     ENDIF
   ENDIF  
 ENDDO
!
!
!-------------------------------------------------------------------------------
!
!*       6.     COMPUTES THE FAST COLD PROCESS SOURCES
!               --------------------------------------
!
!*       6.1    compute the slope parameter Lbda_s and Lbda_g
!
 WHERE ( ZRST(:)>0.0 )
   ZLBDAS(:)  = MIN( XLBDAS_MAX,                                           &
                     XLBS*( ZRHODREF(:)*MAX( ZRST(:),PRTMIN_AQ*1.e3/ZRHODREF(:)) )**XLBEXS )
 END WHERE
!
 WHERE ( ZRGT(:)>0.0 )
   ZLBDAG(:)  = XLBG*( ZRHODREF(:)*MAX( ZRGT(:),PRTMIN_AQ*1.e3/ZRHODREF(:)))**XLBEXG
 END WHERE
!
!*       6.2    cloud droplet riming of the aggregates
!
 ZZW1(:,:) = 0.0
 ZZW(:) = 0.0

 ALLOCATE(GRIM(IMICRO))
 GRIM(:) = (ZRCT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. & 
           (ZRST(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. &
           (ZZRCS(:)>0.0) .AND. (ZZT(:)<XTT)
 IGRIM = COUNT( GRIM(:) )
!
 IF( IGRIM>0 ) THEN
!
!        6.2.0  allocations
!
   ALLOCATE(ZVEC1(IGRIM))
   ALLOCATE(ZVEC2(IGRIM))
   ALLOCATE(IVEC1(IGRIM))
   ALLOCATE(IVEC2(IGRIM))
!
!        6.2.1  select the ZLBDAS
!
   ZVEC1(:) = PACK( ZLBDAS(:),MASK=GRIM(:) )
!
!        6.2.2  find the next lower indice for the ZLBDAS in the geometrical
!               set of Lbda_s used to tabulate some moments of the incomplete
!               gamma function
!
   ZVEC2(1:IGRIM) = MAX( 1.00001, MIN( REAL(NGAMINC)-0.00001,           &
                         XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + XRIMINTP2 ) )
   IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
   ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
!
!        6.2.3  perform the linear interpolation of the normalized
!               "2+XDS"-moment of the incomplete gamma function
!
   ZVEC1(1:IGRIM) =   XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
   ZZW(:) = UNPACK( VECTOR=ZVEC1(:),MASK=GRIM,FIELD=0.0 )
!
!        6.2.4  riming of the small sized aggregates
!
   ZZW2(:,:) = 0.0
   DO JL = 1,IMICRO
     IF ( GRIM(JL) ) THEN
       ZZW1(JL,1) = MIN( ZZRCS(JL), XCRIMSS * ZZW(JL) * ZRCT(JL) * ZRST(JL)  & ! RCRIMSS
                         * ZLBDAS(JL)**(XBS+XEXCRIMSS) * ZRHODREF(JL)**(-XCEXVT+1) )
       ZZW2(JL,:) = ZZW1(JL,1) * ZCSVT(JL,:)/ZRCT(JL)
       ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZCSVT(JL,:)/PTSTEP)),0.0)
       ZCRSVS(JL,:) = ZCRSVS(JL,:) - ZZW2(JL,:)
       IF (OUSECHIC) THEN
       DO JLI = 1, SIZE(PSGRSVS,4)
         IF (TRIM(HICNAMES(JLI)) == 'IC_HNO3' .OR. TRIM(HICNAMES(JLI)) == 'IC_SULF' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_H2SO4' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_NH3' .OR. TRIM(HICNAMES(JLI)) == 'IC_HCL' &
            .OR. HICNAMES(JLI)(1:4) == 'IC_A' .OR. HICNAMES(JLI)(1:4) == 'IC_B' &
            .OR. NINDEXGI(JLI).EQ.0) THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETNA * ZZW2(JL,NINDEXWI(JLI))
         ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_H2O2' .OR. TRIM(HICNAMES(JLI)) == 'IC_HO2' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_HONO' .OR. TRIM(HICNAMES(JLI)) == 'IC_HNO4'&
            .OR. TRIM(HICNAMES(JLI)) == 'IC_HCHO' .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA1'&
            .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA2') THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETHP * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETHP) * ZZW2(JL,NINDEXWI(JLI))
         ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_SO2' .OR. TRIM(HICNAMES(JLI)) == 'IC_OH' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_MO2' .OR. &
                 TRIM(HICNAMES(JLI)) == 'IC_OP1') THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETSU * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETSU) * ZZW2(JL,NINDEXWI(JLI))
         ELSE
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETDF * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETDF) * ZZW2(JL,NINDEXWI(JLI))
         ENDIF
       ENDDO
       ELSE
         IF (.NOT.(OCH_RET_ICE)) THEN
           DO JLW = 1, SIZE(PCRSVS,4)
             IF (.NOT.(NINDEXWG(JLW).EQ.0)) THEN
               ZGRSVS(JL,NINDEXWG(JLW)) = ZGRSVS(JL,NINDEXWG(JLW)) + ZZW2(JL,JLW)
             ENDIF
           ENDDO
         ENDIF
       ENDIF
     ENDIF
   ENDDO
!
!        6.2.5  riming-conversion of the large sized aggregates into graupel
!
   ZZW2(:,:) = 0.0
   DO JL = 1,IMICRO
     IF ( GRIM(JL) .AND. (ZZRSS(JL)>0.0) ) THEN
       ZZW1(JL,2) = MIN( ZZRCS(JL), XCRIMSG * ZRCT(JL) * ZRST(JL) * ZLBDAS(JL)**(XBS+XEXCRIMSG)  & ! RCRIMSG
                        * ZRHODREF(JL)**(-XCEXVT+1) - ZZW1(JL,1) )
       ZZW2(JL,:) = ZZW1(JL,2) * ZCSVT(JL,:)/ZRCT(JL)
       ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZCSVT(JL,:)/PTSTEP)),0.0)
       ZCRSVS(JL,:) = ZCRSVS(JL,:) - ZZW2(JL,:)
       IF (OUSECHIC) THEN
       DO JLI = 1, SIZE(PSGRSVS,4)
         IF (TRIM(HICNAMES(JLI)) == 'IC_HNO3' .OR. TRIM(HICNAMES(JLI)) == 'IC_SULF' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_H2SO4' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_NH3' .OR. TRIM(HICNAMES(JLI)) == 'IC_HCL' &
            .OR. HICNAMES(JLI)(1:4) == 'IC_A' .OR. HICNAMES(JLI)(1:4) == 'IC_B' &
            .OR. NINDEXGI(JLI).EQ.0) THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETNA * ZZW2(JL,NINDEXWI(JLI))
         ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_H2O2' .OR. TRIM(HICNAMES(JLI)) == 'IC_HO2' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_HONO' .OR. TRIM(HICNAMES(JLI)) == 'IC_HNO4'&
            .OR. TRIM(HICNAMES(JLI)) == 'IC_HCHO' .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA1'&
            .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA2') THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETHP * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETHP) * ZZW2(JL,NINDEXWI(JLI))
         ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_SO2' .OR. TRIM(HICNAMES(JLI)) == 'IC_OH' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_MO2' .OR. &
                 TRIM(HICNAMES(JLI)) == 'IC_OP1') THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETSU * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETSU) * ZZW2(JL,NINDEXWI(JLI))
         ELSE
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETDF * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETDF) * ZZW2(JL,NINDEXWI(JLI))
         ENDIF
       ENDDO
       ELSE
         IF (.NOT.(OCH_RET_ICE)) THEN
           DO JLW = 1, SIZE(PCRSVS,4)
             IF (.NOT.(NINDEXWG(JLW).EQ.0)) THEN
               ZGRSVS(JL,NINDEXWG(JLW)) = ZGRSVS(JL,NINDEXWG(JLW)) + ZZW2(JL,JLW)
             ENDIF
           ENDDO
         ENDIF
       ENDIF
     ENDIF
   ENDDO

   DEALLOCATE(IVEC2)
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)
 END IF
 DEALLOCATE(GRIM)
!
!*       6.3    rain accretion onto the aggregates
!
 ZZW(:) = 0.0
 ZZW1(:,2:3) = 0.0
 ALLOCATE(GACC(IMICRO))
 GACC(:) = (ZRRT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. &
           (ZRST(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. &
           (ZZRRS(:)>0.0) .AND. (ZZT(:)<XTT)
 IGACC = COUNT( GACC(:) )
!
 IF( IGACC>0 ) THEN
!
!        6.3.0  allocations
!
   ALLOCATE(ZVEC1(IGACC))
   ALLOCATE(ZVEC2(IGACC))
   ALLOCATE(ZVEC3(IGACC))
   ALLOCATE(IVEC1(IGACC))
   ALLOCATE(IVEC2(IGACC))
!
!        6.3.1  select the (ZLBDAS,ZLBDAR) couplet
!
   ZVEC1(:) = PACK( ZLBDAS(:),MASK=GACC(:) )
   ZVEC2(:) = PACK( ZLBDAR(:),MASK=GACC(:) )
!
!        6.3.2  find the next lower indice for the ZLBDAS and for the ZLBDAR
!               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
!               tabulate the RACCSS-kernel
!
   ZVEC1(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAS)-0.00001,           &
                         XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + XACCINTP2S ) )
   IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
   ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - REAL( IVEC1(1:IGACC) )
!
   ZVEC2(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAR)-0.00001,           &
                         XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + XACCINTP2R ) )
   IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
   ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - REAL( IVEC2(1:IGACC) )
!
!        6.3.3  perform the bilinear interpolation of the normalized
!               RACCSS-kernel
!
   DO JJ = 1,IGACC
     ZVEC3(JJ) =  (  XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                - (  XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                            * (ZVEC1(JJ) - 1.0)
   END DO
   ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GACC,FIELD=0.0 )
!
!        6.3.4  raindrop accretion on the small sized aggregates
!
   ZZW2(:,:) = 0.0
   DO JL = 1,IMICRO
     IF ( GACC(JL) ) THEN
       ZZW1(JL,2) =                                            & !! coef of RRACCS
              XFRACCSS*( ZRST(JL)*ZLBDAS(JL)**XBS )*( ZRHODREF(JL)**(-XCEXVT) ) &
         *( XLBRACCS1/((ZLBDAS(JL)**2)               ) +                  &
            XLBRACCS2/( ZLBDAS(JL)    * ZLBDAR(JL)    ) +                  &
            XLBRACCS3/(               (ZLBDAR(JL)**2)) )/ZLBDAR(JL)**4
       ZZW1(JL,4) = MIN( ZZRRS(JL),ZZW1(JL,2)*ZZW(JL) )           ! RRACCSS
       ZZW2(JL,:) = ZZW1(JL,4) * ZRSVT(JL,:)/ZRRT(JL)
       ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZRSVT(JL,:)/PTSTEP)),0.0)
       ZRRSVS(JL,:) = ZRRSVS(JL,:) - ZZW2(JL,:) 
       IF (OUSECHIC) THEN
       DO JLI = 1, SIZE(PSGRSVS,4)
         IF (TRIM(HICNAMES(JLI)) == 'IC_HNO3' .OR. TRIM(HICNAMES(JLI)) == 'IC_SULF' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_H2SO4' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_NH3' .OR. TRIM(HICNAMES(JLI)) == 'IC_HCL' &
            .OR. HICNAMES(JLI)(1:4) == 'IC_A' .OR. HICNAMES(JLI)(1:4) == 'IC_B' &
            .OR. NINDEXGI(JLI).EQ.0) THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETNA * ZZW2(JL,NINDEXWI(JLI))
         ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_H2O2' .OR. TRIM(HICNAMES(JLI)) == 'IC_HO2' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_HONO' .OR. TRIM(HICNAMES(JLI)) == 'IC_HNO4'&
            .OR. TRIM(HICNAMES(JLI)) == 'IC_HCHO' .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA1'&
            .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA2') THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETHP * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETHP) * ZZW2(JL,NINDEXWI(JLI))
         ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_SO2' .OR. TRIM(HICNAMES(JLI)) == 'IC_OH' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_MO2' .OR. &
                 TRIM(HICNAMES(JLI)) == 'IC_OP1') THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETSU * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETSU) * ZZW2(JL,NINDEXWI(JLI))
         ELSE
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETDF * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETDF) * ZZW2(JL,NINDEXWI(JLI))
         ENDIF
       ENDDO
       ELSE
         IF (.NOT.(OCH_RET_ICE)) THEN
           DO JLW = 1, SIZE(PRRSVS,4)
             IF (.NOT.(NINDEXWG(JLW).EQ.0)) THEN
               ZGRSVS(JL,NINDEXWG(JLW)) = ZGRSVS(JL,NINDEXWG(JLW)) + ZZW2(JL,JLW)
             ENDIF
           ENDDO
         ENDIF
       ENDIF
     ENDIF
   ENDDO
!
!        6.3.4b perform the bilinear interpolation of the normalized
!               RACCS-kernel
!
   DO JJ = 1,IGACC
     ZVEC3(JJ) =  (   XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                   -  XKER_RACCS(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                                  * ZVEC2(JJ) &
                - (   XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                   -  XKER_RACCS(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                          * (ZVEC2(JJ) - 1.0)
   END DO
   ZZW1(:,2) = ZZW1(:,2)*UNPACK( VECTOR=ZVEC3(:),MASK=GACC(:),FIELD=0.0 )
!
!        6.3.5  raindrop accretion-conversion of the large sized aggregates
!               into graupeln
!
   ZZW2(:,:) = 0.0
   WHERE ( GACC(:) .AND. (ZZRSS(:)>0.0) )
     ZZW1(:,2) = MAX( MIN( ZZRRS(:),ZZW1(:,2)-ZZW1(:,4) ),0.0 )       ! RRACCSG
   END WHERE
   DO JL = 1,IMICRO
     IF ( GACC(JL) .AND. (ZZRSS(JL)>0.0) .AND. ZZW1(JL,2)>0.0 ) THEN
       ZZW2(JL,:) = ZZW1(JL,2) * ZRSVT(JL,:)/ZRRT(JL)
       ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZRSVT(JL,:)/PTSTEP)),0.0)
       ZRRSVS(JL,:) = ZRRSVS(JL,:) - ZZW2(JL,:) 
       IF (OUSECHIC) THEN
       DO JLI = 1, SIZE(PSGRSVS,4)
         IF (TRIM(HICNAMES(JLI)) == 'IC_HNO3' .OR. TRIM(HICNAMES(JLI)) == 'IC_SULF' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_H2SO4' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_NH3' .OR. TRIM(HICNAMES(JLI)) == 'IC_HCL' &
            .OR. HICNAMES(JLI)(1:4) == 'IC_A' .OR. HICNAMES(JLI)(1:4) == 'IC_B' &
            .OR. NINDEXGI(JLI).EQ.0) THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETNA * ZZW2(JL,NINDEXWI(JLI))
         ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_H2O2' .OR. TRIM(HICNAMES(JLI)) == 'IC_HO2' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_HONO' .OR. TRIM(HICNAMES(JLI)) == 'IC_HNO4'&
            .OR. TRIM(HICNAMES(JLI)) == 'IC_HCHO' .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA1'&
            .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA2') THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETHP * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETHP) * ZZW2(JL,NINDEXWI(JLI))
         ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_SO2' .OR. TRIM(HICNAMES(JLI)) == 'IC_OH' &
            .OR. TRIM(HICNAMES(JLI)) == 'IC_MO2' .OR. &
                 TRIM(HICNAMES(JLI)) == 'IC_OP1') THEN
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETSU * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETSU) * ZZW2(JL,NINDEXWI(JLI))
         ELSE
           ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETDF * ZZW2(JL,NINDEXWI(JLI))
           ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                     (1. - XRETDF) * ZZW2(JL,NINDEXWI(JLI))
         ENDIF
       ENDDO
       ELSE
         IF (.NOT.(OCH_RET_ICE)) THEN
           DO JLW = 1, SIZE(PRRSVS,4)
             IF (.NOT.(NINDEXWG(JLW).EQ.0)) THEN
               ZGRSVS(JL,NINDEXWG(JLW)) = ZGRSVS(JL,NINDEXWG(JLW)) + ZZW2(JL,JLW)
             ENDIF
           ENDDO
         ENDIF
       ENDIF
     ENDIF
   ENDDO
!                                                           
   DEALLOCATE(IVEC2)
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC3)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)
 END IF
 DEALLOCATE(GACC)
!
!*       6.4    rain contact freezing
!
 ZZW1(:,4) = 0.0
 ZZW2(:,:) = 0.0
 DO JL = 1,IMICRO
   IF ( (ZRIT(JL)>PRTMIN_AQ*1.e3/ZRHODREF(JL)) .AND.  &
        (ZRRT(JL)>PRTMIN_AQ*1.e3/ZRHODREF(JL)) .AND.  &
        (ZZRIS(JL)>0.0) .AND. (ZZRRS(JL)>0.0) ) THEN
     ZZW1(JL,4) = MIN( ZZRRS(JL), XRCFRI * ZCIT(JL)       & ! RRCFRIG
                                 * ZLBDAR(JL)**XEXRCFRI    &
                                 * ZRHODREF(JL)**(-XCEXVT-1.) )
     ZZW2(JL,:) = ZZW1(JL,4) * ZRSVT(JL,:)/ZRRT(JL)
     ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZRSVT(JL,:)/PTSTEP)),0.0)
     ZRRSVS(JL,:) = ZRRSVS(JL,:) - ZZW2(JL,:) 
     IF (OUSECHIC) THEN
     DO JLI = 1, SIZE(PSGRSVS,4)
       IF (TRIM(HICNAMES(JLI)) == 'IC_HNO3' .OR. TRIM(HICNAMES(JLI)) == 'IC_SULF' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_H2SO4' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_NH3' .OR. TRIM(HICNAMES(JLI)) == 'IC_HCL' &
          .OR. HICNAMES(JLI)(1:4) == 'IC_A' .OR. HICNAMES(JLI)(1:4) == 'IC_B' &
          .OR. NINDEXGI(JLI).EQ.0) THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETNA * ZZW2(JL,NINDEXWI(JLI))
       ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_H2O2' .OR. TRIM(HICNAMES(JLI)) == 'IC_HO2' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_HONO' .OR. TRIM(HICNAMES(JLI)) == 'IC_HNO4'&
          .OR. TRIM(HICNAMES(JLI)) == 'IC_HCHO' .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA1'&
          .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA2') THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETHP * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETHP) * ZZW2(JL,NINDEXWI(JLI))
       ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_SO2' .OR. TRIM(HICNAMES(JLI)) == 'IC_OH' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_MO2' .OR. &
               TRIM(HICNAMES(JLI)) == 'IC_OP1') THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETSU * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETSU) * ZZW2(JL,NINDEXWI(JLI))
       ELSE
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETDF * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETDF) * ZZW2(JL,NINDEXWI(JLI))
       ENDIF
     ENDDO
     ELSE
       IF (.NOT.(OCH_RET_ICE)) THEN
         DO JLW = 1, SIZE(PRRSVS,4)
           IF (.NOT.(NINDEXWG(JLW).EQ.0)) THEN
             ZGRSVS(JL,NINDEXWG(JLW)) = ZGRSVS(JL,NINDEXWG(JLW)) + ZZW2(JL,JLW)
           ENDIF
         ENDDO
       ENDIF
     ENDIF
   ENDIF
 ENDDO
!
!*       6.5    compute the Dry growth case of graupel
!
 ZZW(:) = 0.0
 ZZW1(:,:) = 0.0
 WHERE( (ZRGT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND.   &
        ((ZRCT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:) .AND. ZZRCS(:)>0.0)) )
   ZZW(:) = ZLBDAG(:)**(XCXG-XDG-2.0) * ZRHODREF(:)**(-XCEXVT)
   ZZW1(:,1) = MIN( ZZRCS(:),XFCDRYG * ZRCT(:) * ZZW(:) )             ! RCDRYG
 END WHERE
 WHERE( (ZRGT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. &
        ((ZRIT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:) .AND. ZZRIS(:)>0.0)) )
   ZZW(:) = ZLBDAG(:)**(XCXG-XDG-2.0) * ZRHODREF(:)**(-XCEXVT)
   ZZW1(:,2) = MIN( ZZRIS(:),XFIDRYG * EXP( XCOLEXIG*(ZZT(:)-XTT) ) &
                                    * ZRIT(:) * ZZW(:) )             ! RIDRYG
 END WHERE
!
!        6.5.1  accretion of aggregates on the graupeln
!
 ALLOCATE(GDRY(IMICRO))
 GDRY(:) = (ZRST(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. &
           (ZRGT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. (ZZRSS(:)>0.0)
 IGDRY = COUNT( GDRY(:) )
!
 IF( IGDRY>0 ) THEN
!
!        6.5.2  allocations
!
   ALLOCATE(ZVEC1(IGDRY))
   ALLOCATE(ZVEC2(IGDRY))
   ALLOCATE(ZVEC3(IGDRY))
   ALLOCATE(IVEC1(IGDRY))
   ALLOCATE(IVEC2(IGDRY))
!
!        6.5.3  select the (ZLBDAG,ZLBDAS) couplet
!
   ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
   ZVEC2(:) = PACK( ZLBDAS(:),MASK=GDRY(:) )
!
!        6.5.4  find the next lower indice for the ZLBDAG and for the ZLBDAS
!               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
!               tabulate the SDRYG-kernel
!
   ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAG)-0.00001,           &
                         XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + XDRYINTP2G ) )
   IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
   ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
   ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAS)-0.00001,           &
                         XDRYINTP1S * LOG( ZVEC2(1:IGDRY) ) + XDRYINTP2S ) )
   IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
   ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
!
!        6.5.5  perform the bilinear interpolation of the normalized
!               SDRYG-kernel
!
   DO JJ = 1,IGDRY
     ZVEC3(JJ) =  (  XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                - (  XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
   END DO
   ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
   WHERE( GDRY(:) )
     ZZW1(:,3) = MIN( ZZRSS(:),XFSDRYG*ZZW(:)                         & ! RSDRYG
                                     * EXP( XCOLEXSG*(ZZT(:)-XTT) )  &
                   *ZRST(:)*( ZLBDAG(:)**XCXG )    &
                   *( ZRHODREF(:)**(-XCEXVT) )                    &
                        *( XLBSDRYG1/( ZLBDAG(:)**2              ) + &
                           XLBSDRYG2/( ZLBDAG(:)   * ZLBDAS(:)   ) + &
                           XLBSDRYG3/(               ZLBDAS(:)**2) ) )
   END WHERE
   DEALLOCATE(IVEC2)
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC3)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)
 END IF
!
!        6.5.6  accretion of raindrops on the graupeln
!
 GDRY(:) = (ZRRT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. &
           (ZRGT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:)) .AND. (ZZRRS(:)>0.0)
 IGDRY = COUNT( GDRY(:) )
!
 IF( IGDRY>0 ) THEN
!
!        6.5.7  allocations
!
   ALLOCATE(ZVEC1(IGDRY))
   ALLOCATE(ZVEC2(IGDRY))
   ALLOCATE(ZVEC3(IGDRY))
   ALLOCATE(IVEC1(IGDRY))
   ALLOCATE(IVEC2(IGDRY))
!
!        6.5.8  select the (ZLBDAG,ZLBDAR) couplet
!
   ZVEC1(:) = PACK( ZLBDAG(:),MASK=GDRY(:) )
   ZVEC2(:) = PACK( ZLBDAR(:),MASK=GDRY(:) )
!
!        6.5.9  find the next lower indice for the ZLBDAG and for the ZLBDAR
!               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
!               tabulate the RDRYG-kernel
!
   ZVEC1(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAG)-0.00001,           &
                         XDRYINTP1G * LOG( ZVEC1(1:IGDRY) ) + XDRYINTP2G ) )
   IVEC1(1:IGDRY) = INT( ZVEC1(1:IGDRY) )
   ZVEC1(1:IGDRY) = ZVEC1(1:IGDRY) - REAL( IVEC1(1:IGDRY) )
!
   ZVEC2(1:IGDRY) = MAX( 1.00001, MIN( REAL(NDRYLBDAR)-0.00001,           &
                         XDRYINTP1R * LOG( ZVEC2(1:IGDRY) ) + XDRYINTP2R ) )
   IVEC2(1:IGDRY) = INT( ZVEC2(1:IGDRY) )
   ZVEC2(1:IGDRY) = ZVEC2(1:IGDRY) - REAL( IVEC2(1:IGDRY) )
!
!        6.5.10 perform the bilinear interpolation of the normalized
!               RDRYG-kernel
!
   DO JJ = 1,IGDRY
     ZVEC3(JJ) =  (  XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                - (  XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                            * (ZVEC1(JJ) - 1.0)
   END DO
   ZZW(:) = UNPACK( VECTOR=ZVEC3(:),MASK=GDRY,FIELD=0.0 )
!
   WHERE( GDRY(:) )
     ZZW1(:,4) = MIN( ZZRRS(:), XFRDRYG*ZZW(:)                  & ! RRDRYG
                       *( ZLBDAR(:)**(-4) )*( ZLBDAG(:)**XCXG ) &
                               *( ZRHODREF(:)**(-XCEXVT-1.) )   &
                   *( XLBRDRYG1/( ZLBDAG(:)**2              ) + &
                      XLBRDRYG2/( ZLBDAG(:)   * ZLBDAR(:)   ) + &
                      XLBRDRYG3/(               ZLBDAR(:)**2) ) )
   END WHERE
   DEALLOCATE(IVEC2)
   DEALLOCATE(IVEC1)
   DEALLOCATE(ZVEC3)
   DEALLOCATE(ZVEC2)
   DEALLOCATE(ZVEC1)
 END IF
!
 ZRDRYG(:) = ZZW1(:,1) + ZZW1(:,2) + ZZW1(:,3) + ZZW1(:,4)
 DEALLOCATE(GDRY)
!
!*       6.6    compute the Wet growth case of the graupel
!
 ZZW(:) = 0.0
 ZRWETG(:) = 0.0
!
 ZKA(:) = 2.38E-2 + 0.0071E-2 * ( ZZT(:) - XTT )          ! k_a
 ZDV(:) = 0.211E-4 * (ZZT(:)/XTT)**1.94 * (XP00/ZPRES(:)) ! D_v
 ZCJ(:) = XSCFAC * ZRHODREF(:)**0.3 / SQRT( 1.718E-5+0.0049E-5*(ZZT(:)-XTT) )
                                         !c^prime_j (in the ventilation factor)
 WHERE( ZRGT(:)>PRTMIN_AQ*1.e3/ZRHODREF(:) )
   ZZW1(:,5) = MIN( ZZRIS(:),                                   &
               ZZW1(:,2) / (XCOLIG*EXP(XCOLEXIG*(ZZT(:)-XTT)) ) ) ! RIWETG
   ZZW1(:,6) = MIN( ZZRSS(:),                                    &
               ZZW1(:,3) / (XCOLSG*EXP(XCOLEXSG*(ZZT(:)-XTT)) ) ) ! RSWETG
!
   ZZW(:) = ZRVT(:)*ZPRES(:)/((XMV/XMD)+ZRVT(:)) ! Vapor pressure
   ZZW(:) =   ZKA(:)*(XTT-ZZT(:)) +                              &
            ( ZDV(:)*(XLVTT + ( XCPV - XCL ) * ( ZZT(:) - XTT )) &
                          *(XESTT-ZZW(:))/(XRV*ZZT(:))           )
! compute RWETG
!
   ZRWETG(:)=MAX( 0.0,                                               &
                ( ZZW(:) * ( X0DEPG*       ZLBDAG(:)**XEX0DEPG +     &
                             X1DEPG*ZCJ(:)*ZLBDAG(:)**XEX1DEPG ) +   &
                ( ZZW1(:,5)+ZZW1(:,6) ) *                            &
                ( ZRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-ZZT(:)))   ) ) / &
                                ( ZRHODREF(:)*(XLMTT-XCL*(XTT-ZZT(:))) )   )
 END WHERE
!
!*       6.7    Select Wet or Dry case for the growth of the graupel
!
 ZZW(:) = 0.0
 ZZW2(:,:) = 0.0
 ZZW4(:,:) = 0.0
 DO JL = 1,IMICRO
   IF ( (ZRGT(JL)>PRTMIN_AQ*1.e3/ZRHODREF(JL)) .AND.     &    ! wet case
         ZZT(JL)<XTT .AND. ZRDRYG(JL)>=ZRWETG(JL) .AND.  &
         ZRWETG(JL)>0.0 .AND. ZRCT(JL)>0.0 .AND. ZRRT(JL)>0.0)  THEN 
     ZZW(JL)  = ZRWETG(JL)
     ZZW2(JL,:) = ZZW(JL) * ZRSVT(JL,:)/ZRRT(JL)
     ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZRSVT(JL,:)/PTSTEP)),0.0)
     ZRRSVS(JL,:) = ZRRSVS(JL,:) - ZZW2(JL,:)  ! rain -> graupel
     IF (OUSECHIC) THEN
     ZZW3(:,:) = 0.0
     DO JLI = 1, SIZE(PSGRSVS,4)
       IF (TRIM(HICNAMES(JLI)) == 'IC_HNO3' .OR. TRIM(HICNAMES(JLI)) == 'IC_SULF' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_H2SO4' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_NH3' .OR. TRIM(HICNAMES(JLI)) == 'IC_HCL' &
          .OR. HICNAMES(JLI)(1:4) == 'IC_A' .OR. HICNAMES(JLI)(1:4) == 'IC_B' &
          .OR. NINDEXGI(JLI).EQ.0) THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETNA * ZZW2(JL,NINDEXWI(JLI))
       ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_H2O2' .OR. TRIM(HICNAMES(JLI)) == 'IC_HO2' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_HONO' .OR. TRIM(HICNAMES(JLI)) == 'IC_HNO4'&
          .OR. TRIM(HICNAMES(JLI)) == 'IC_HCHO' .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA1'&
          .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA2') THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETHP * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETHP) * ZZW2(JL,NINDEXWI(JLI))
       ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_SO2' .OR. TRIM(HICNAMES(JLI)) == 'IC_OH' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_MO2' .OR. &
               TRIM(HICNAMES(JLI)) == 'IC_OP1') THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETSU * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETSU) * ZZW2(JL,NINDEXWI(JLI))
       ELSE
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETDF * ZZW2(JL,NINDEXWI(JLI))
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) +                        &
                                   (1. - XRETDF) * ZZW2(JL,NINDEXWI(JLI))
       ENDIF
     ENDDO
     IF (ZRST(JL)>0.0) THEN
       ZZW3(JL,:) = ZZW1(JL,6) * ZSGSVT(JL,:)/ZRST(JL)
       ZZW3(JL,:) = MAX(MIN(ZZW3(JL,:),(ZSGSVT(JL,:)/PTSTEP)),0.0)
       ZSGRSVS(JL,:) = ZSGRSVS(JL,:) - ZZW3(JL,:) !snow->rain
       DO JLI = 1, SIZE(PSGRSVS,4)
         ZRRSVS(JL,NINDEXWI(JLI)) = ZRRSVS(JL,NINDEXWI(JLI)) + ZZW3(JL,JLI)
       ENDDO
     ENDIF
     ELSE
       IF (.NOT.(OCH_RET_ICE)) THEN
         DO JLW = 1, SIZE(PRRSVS,4)
           IF (.NOT.(NINDEXWG(JLW).EQ.0)) THEN
             ZGRSVS(JL,NINDEXWG(JLW)) = ZGRSVS(JL,NINDEXWG(JLW)) + ZZW2(JL,JLW)
           ENDIF
         ENDDO
       ENDIF
     ENDIF
     ZZW4(JL,:) = ZZW1(JL,1) * ZCSVT(JL,:)/ZRCT(JL)
     ZZW4(JL,:) = MAX(MIN(ZZW4(JL,:),(ZCSVT(JL,:)/PTSTEP)),0.0)
     ZCRSVS(JL,:) = ZCRSVS(JL,:) - ZZW4(JL,:) !cloud->rain 
     ZRRSVS(JL,:) = ZRRSVS(JL,:) + ZZW4(JL,:) 
   ELSE IF ( (ZRGT(JL)>PRTMIN_AQ*1.e3/ZRHODREF(JL)) .AND.     &        ! dry case
              ZZT(JL)<XTT .AND. ZRDRYG(JL)<ZRWETG(JL) .AND.  &
              ZRDRYG(JL)>0.0 .AND. ZRCT(JL)>0.0 .AND. ZRRT(JL)>0.0) THEN
     ZZW2(JL,:) = ZZW1(JL,1) * ZCSVT(JL,:)/ZRCT(JL)
     ZZW2(JL,:) = MAX(MIN(ZZW2(JL,:),(ZCSVT(JL,:)/PTSTEP)),0.0)
     ZZW4(JL,:) = ZZW1(JL,4) * ZRSVT(JL,:)/ZRRT(JL) 
     ZZW4(JL,:) = MAX(MIN(ZZW4(JL,:),(ZRSVT(JL,:)/PTSTEP)),0.0)
     ZCRSVS(JL,:) = ZCRSVS(JL,:) - ZZW2(JL,:) 
     ZRRSVS(JL,:) = ZRRSVS(JL,:) - ZZW4(JL,:)
     IF (OUSECHIC) THEN
     DO JLI = 1, SIZE(PSGRSVS,4)
       IF (TRIM(HICNAMES(JLI)) == 'IC_HNO3' .OR. TRIM(HICNAMES(JLI)) == 'IC_SULF' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_H2SO4' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_NH3' .OR. TRIM(HICNAMES(JLI)) == 'IC_HCL' &
          .OR. HICNAMES(JLI)(1:4) == 'IC_A' .OR. HICNAMES(JLI)(1:4) == 'IC_B' &
          .OR. NINDEXGI(JLI).EQ.0) THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETNA * (        &
                           ZZW2(JL,NINDEXWI(JLI)) +  ZZW4(JL,NINDEXWI(JLI)) )
       ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_H2O2' .OR. TRIM(HICNAMES(JLI)) == 'IC_HO2' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_HONO' .OR. TRIM(HICNAMES(JLI)) == 'IC_HNO4'&
          .OR. TRIM(HICNAMES(JLI)) == 'IC_HCHO' .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA1'&
          .OR. TRIM(HICNAMES(JLI)) == 'IC_ORA2') THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETHP *  (       &
                           ZZW2(JL,NINDEXWI(JLI)) +  ZZW4(JL,NINDEXWI(JLI)) )
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) + (1. - XRETHP) * (  &
                           ZZW2(JL,NINDEXWI(JLI)) +  ZZW4(JL,NINDEXWI(JLI)) )
       ELSE IF (TRIM(HICNAMES(JLI)) == 'IC_SO2' .OR. TRIM(HICNAMES(JLI)) == 'IC_OH' &
          .OR. TRIM(HICNAMES(JLI)) == 'IC_MO2' .OR. &
               TRIM(HICNAMES(JLI)) == 'IC_OP1') THEN
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETSU *  (       &
                           ZZW2(JL,NINDEXWI(JLI)) +  ZZW4(JL,NINDEXWI(JLI)) )
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) + (1. - XRETSU) * (  &
                           ZZW2(JL,NINDEXWI(JLI)) +  ZZW4(JL,NINDEXWI(JLI)) )
       ELSE
         ZSGRSVS(JL,JLI) = ZSGRSVS(JL,JLI) + XRETDF *  (       &
                           ZZW2(JL,NINDEXWI(JLI)) +  ZZW4(JL,NINDEXWI(JLI)) )
         ZGRSVS(JL,NINDEXGI(JLI)) = ZGRSVS(JL,NINDEXGI(JLI)) + (1. - XRETDF) * (  &
                           ZZW2(JL,NINDEXWI(JLI)) +  ZZW4(JL,NINDEXWI(JLI)) )
       ENDIF
     ENDDO
     ELSE
       IF (.NOT.(OCH_RET_ICE)) THEN
         DO JLW = 1, SIZE(PRRSVS,4)
           IF (.NOT.(NINDEXWG(JLW).EQ.0)) THEN
             ZGRSVS(JL,NINDEXWG(JLW)) = ZGRSVS(JL,NINDEXWG(JLW)) + ZZW2(JL,JLW) &
                                                              + ZZW4(JL,JLW)
           ENDIF
         ENDDO
       ENDIF
     ENDIF
   ENDIF
 ENDDO
!
!*       6.8    Melting of the graupel
!
 IF (OUSECHIC) THEN
 ZZW(:) = 0.0
 ZZW3(:,:) = 0.0
 DO JL = 1,IMICRO
   IF ( (ZRGT(JL)>PRTMIN_AQ*1.e3/ZRHODREF(JL)) .AND.  &
        (ZZRGS(JL)>0.0) .AND. (ZZT(JL)>XTT) ) THEN
     ZZW(JL) = ZRVT(JL)*ZPRES(JL)/((XMV/XMD)+ZRVT(JL)) ! Vapor pressure
     ZZW(JL) = ZKA(JL)*(XTT-ZZT(JL)) +                                 &
             ( ZDV(JL)*(XLVTT + ( XCPV - XCL ) * ( ZZT(JL) - XTT )) &
                           *(XESTT-ZZW(JL))/(XRV*ZZT(JL))             )
! compute RGMLTR
     ZZW(JL)  = MIN( ZZRGS(JL), MAX( 0.0,( -ZZW(JL) *                     &
                           ( X0DEPG*       ZLBDAG(JL)**XEX0DEPG +     &
                             X1DEPG*ZCJ(JL)*ZLBDAG(JL)**XEX1DEPG ) -   &
                                     ( ZZW1(JL,1)+ZZW1(JL,4) ) *       &
                            ( ZRHODREF(JL)*XCL*(XTT-ZZT(JL))) ) /    &
                                             ( ZRHODREF(JL)*XLMTT ) ) )
     ZZW3(JL,:) = ZZW(JL) * ZSGSVT(JL,:)/ZRGT(JL)
     ZZW3(JL,:) = MAX(MIN(ZZW3(JL,:),(ZSGSVT(JL,:)/PTSTEP)),0.0)
     ZSGRSVS(JL,:) = ZSGRSVS(JL,:) - ZZW3(JL,:) !graupel->rain
     DO JLI = 1, SIZE(PSGRSVS,4)
       ZRRSVS(JL,NINDEXWI(JLI)) = ZRRSVS(JL,NINDEXWI(JLI)) + ZZW3(JL,JLI)
     ENDDO
   ENDIF
 ENDDO
 ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*       7.     UNPACK RESULTS AND DEALLOCATE ARRAYS
!               ------------------------------------


 DO JLC= 1, SIZE(PCRSVS,4)
   ZCW(:,:,:) = ZZCRSVS(:,:,:,JLC)
   ZZCRSVS(:,:,:,JLC) = UNPACK(ZCRSVS(:,JLC), MASK=GMICRO(:,:,:), FIELD=ZCW(:,:,:))
   PCRSVS(:,:,:,JLC) = ZZCRSVS(:,:,:,JLC) * PRHODJ(:,:,:)
 END DO
 DO JLR= 1, SIZE(PRRSVS,4)
   ZRW(:,:,:) = ZZRRSVS(:,:,:,JLR)
   ZZRRSVS(:,:,:,JLR) = UNPACK(ZRRSVS(:,JLR), MASK=GMICRO(:,:,:), FIELD=ZRW(:,:,:))
   PRRSVS(:,:,:,JLR) = ZZRRSVS(:,:,:,JLR) * PRHODJ(:,:,:)
 END DO
 IF (OUSECHIC) THEN
   DO JLI= 1, SIZE(PSGRSVS,4)
     ZSGW(:,:,:) = ZZSGRSVS(:,:,:,JLI)
     ZZSGRSVS(:,:,:,JLI) = UNPACK(ZSGRSVS(:,JLI), MASK=GMICRO(:,:,:), FIELD=ZSGW(:,:,:))
     PSGRSVS(:,:,:,JLI) = ZZSGRSVS(:,:,:,JLI) * PRHODJ(:,:,:)
   END DO
   DO JLG= 1, SIZE(PGRSVS,4)
     ZGW(:,:,:) = ZZGRSVS(:,:,:,JLG)
     ZZGRSVS(:,:,:,JLG) = UNPACK(ZGRSVS(:,JLG), MASK=GMICRO(:,:,:), FIELD=ZGW(:,:,:))
     PGRSVS(:,:,:,JLG) = ZZGRSVS(:,:,:,JLG) * PRHODJ(:,:,:)
   END DO
   DEALLOCATE(ZGRSVS)
   DEALLOCATE(ZSGRSVS)
   DEALLOCATE(ZSGSVT)
   DEALLOCATE(ZZW3)
 ELSE
   IF (.NOT.(OCH_RET_ICE)) THEN
     DO JLG= 1, SIZE(PGRSVS,4)
       ZGW(:,:,:) = ZZGRSVS(:,:,:,JLG)
       ZZGRSVS(:,:,:,JLG) = UNPACK(ZGRSVS(:,JLG), MASK=GMICRO(:,:,:), FIELD=ZGW(:,:,:))
       PGRSVS(:,:,:,JLG) = ZZGRSVS(:,:,:,JLG) * PRHODJ(:,:,:)
     END DO
     DEALLOCATE(ZGRSVS)
   ENDIF
 ENDIF

  DEALLOCATE(ZRHODREF) 
  DEALLOCATE(ZZT) 
  DEALLOCATE(ZPRES) 
  DEALLOCATE(ZKA) 
  DEALLOCATE(ZDV) 
  DEALLOCATE(ZCJ) 
  DEALLOCATE(ZZW) 
  DEALLOCATE(ZZW1) 
  DEALLOCATE(ZZW2) 
  DEALLOCATE(ZZW4)
  DEALLOCATE(ZZRCS)
  DEALLOCATE(ZZRRS)
  DEALLOCATE(ZZRIS)
  DEALLOCATE(ZZRSS)
  DEALLOCATE(ZZRGS)
  DEALLOCATE(ZCRSVS)
  DEALLOCATE(ZRRSVS)
  DEALLOCATE(ZRVT)
  DEALLOCATE(ZRCT)
  DEALLOCATE(ZRRT)
  DEALLOCATE(ZRIT)
  DEALLOCATE(ZRST)
  DEALLOCATE(ZRGT)
  DEALLOCATE(ZCIT)
  DEALLOCATE(ZCSVT)
  DEALLOCATE(ZRSVT)
  DEALLOCATE(ZLBDAR)
  DEALLOCATE(ZLBDAS)
  DEALLOCATE(ZLBDAG)
  DEALLOCATE(ZRDRYG)
!
END IF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE CH_AQUEOUS_TMICICE
