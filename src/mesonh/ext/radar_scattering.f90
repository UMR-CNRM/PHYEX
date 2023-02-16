!MNH_LIC Copyright 2004-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
       MODULE MODI_RADAR_SCATTERING 
!      #############################
!
INTERFACE
    SUBROUTINE RADAR_SCATTERING(PT_RAY,PRHODREF_RAY,PR_RAY,PI_RAY,PCIT_RAY,PS_RAY,PG_RAY,PVDOP_RAY, &
   PELEV,PX_H,PX_V,PW_H,PW_V,PZE,PBU_MASK_RAY,PCR_RAY,PH_RAY)
REAL, DIMENSION(:,:,:,:,:,:),INTENT(IN)  :: PT_RAY ! temperature interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),INTENT(IN)  :: PRHODREF_RAY ! 
REAL, DIMENSION(:,:,:,:,:,:),INTENT(IN)  :: PR_RAY  ! rainwater mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),INTENT(IN)  :: PI_RAY  ! pristine ice mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PCIT_RAY  ! pristine ice concentration interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PS_RAY !aggregates mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PG_RAY  ! graupel         mixing ratio interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PVDOP_RAY !Doppler radial velocity interpolated along the rays
REAL, DIMENSION(:,:,:,:),     INTENT(IN) :: PELEV ! elevation
REAL, DIMENSION(:),           INTENT(IN) :: PX_H ! Gaussian horizontal nodes
REAL, DIMENSION(:),           INTENT(IN) :: PX_V ! Gaussian vertical nodes
REAL, DIMENSION(:),           INTENT(IN) :: PW_H ! Gaussian horizontal weights
REAL, DIMENSION(:),           INTENT(IN) :: PW_V ! Gaussian vertical weights
REAL,DIMENSION(:,:,:,:,:),    INTENT(INOUT) :: PZE ! 5D matrix (iradar, ielev, iaz, irangestep, ivar) containing the radar variables that will be calculated
!in polar or cartesian projection (same projection as the observation grid)
! convective/stratiform
REAL, DIMENSION(:,:,:,:,:,:),INTENT(INOUT) :: PBU_MASK_RAY
REAL, DIMENSION(:,:,:,:,:,:),OPTIONAL,INTENT(IN)  :: PCR_RAY  ! rainwater concentration interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),OPTIONAL,INTENT(IN)  :: PH_RAY   ! hail mixing ratio interpolated along the rays
    END SUBROUTINE RADAR_SCATTERING
END INTERFACE
END MODULE MODI_RADAR_SCATTERING
!
!     ######spl
       SUBROUTINE RADAR_SCATTERING(PT_RAY,PRHODREF_RAY,PR_RAY,PI_RAY,PCIT_RAY, &
            PS_RAY,PG_RAY,PVDOP_RAY,PELEV,PX_H,PX_V,PW_H,PW_V,PZE,PBU_MASK_RAY,PCR_RAY,PH_RAY)
!     ##############################
!
!!****  *RADAR_SCATTERING* - computes radar reflectivities.
!!
!!    PURPOSE
!!    -------
!!      Compute equivalent reflectivities of a mixed phase cloud.
!!
!!**  METHOD
!!    ------
!!      The reflectivities are computed using the n(D) * sigma(D) formula. The 
!!    equivalent reflectiviy is the sum of the reflectivity produced by the
!!    the raindrops and the equivalent reflectivities of the ice crystals.
!!    The latter are computed using the mass-equivalent diameter.
!!    Four types of diffusion are possible : Rayleigh, Mie, T-matrix, and
!!    Rayleigh-Gans (Kerker, 1969, Chap. 10; Battan, 1973, Sec. 5.4; van de
!!    Hulst, 1981, Sec. 6.32; Doviak and Zrnic, 1993, p. 249; Bringi and 
!!    Chandrasekar, 2001, Chap. 2).
!!    The integration over diameters for Mie and T-matrix methods is done by
!!    using Gauss-Laguerre quadrature (Press et al. 1986). Attenuation is taken
!!    into account by computing the extinction efficiency and correcting 
!!    reflectivities along the beam path.
!!    Gaussian quadrature methods are used to model the beam broadening (Gauss-
!!    Hermite or Gauss-Legendre, see Press et al. 1986).
!!      
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XLIGHTSPEED
!!        XPI
!!      Module MODD_ARF
!!
!!    REFERENCE
!!    ---------
!!      Press, W. H., B. P. Flannery, S. A. Teukolsky et W. T. Vetterling, 1986: 
!!    Numerical Recipes: The Art of Scientific Computing. Cambridge University 
!!    Press, 818 pp.
!!      Probert-Jones, J. R., 1962 : The radar equation in meteorology. Quart. 
!!    J. Roy. Meteor. Soc., 88, 485-495.
!!
!!    AUTHOR
!!    ------
!!      O. Caumont & V. Ducrocq      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   26/03/2004 
!!      O. Caumont 09/09/2009 minor changes to compute radial velocities when no
!!                              hydrometeors so as to emulate wind lidar
!!      O. Caumont 21/12/2009 correction of bugs to compute KDP.
!!      O. Caumont 11/02/2010 thresholding and conversion from linear to 
!!          log values after interpolation instead of before.
!!      G.Tanguy 25/03/2010 Introduction of MODD_TMAT and ALLOCATE/DEALLOCATE
!!      C.Augros  2014 New simulator for T matrice
!!      G.Delautier 10/2014 : Mise a jour simulateur T-matrice pour LIMA
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_IO,           ONLY: TFILEDATA
USE MODD_LUNIT
USE MODD_PARAMETERS
USE MODD_PARAM_ICE,      ONLY: LSNOW_T_I=>LSNOW_T
USE MODD_RAIN_ICE_DESCR, ONLY: XALPHAR_I=>XALPHAR,XNUR_I=>XNUR,XDR_I=>XDR,XLBEXR_I=>XLBEXR,&
                               XLBR_I=>XLBR,XCCR_I=>XCCR,XBR_I=>XBR,XCR_I=>XCR,&
                               XALPHAS_I=>XALPHAS,XNUS_I=>XNUS,XDS_I=>XDS,XLBEXS_I=>XLBEXS,&
                               XLBS_I=>XLBS,XCCS_I=>XCCS,XNS_I=>XNS,XAS_I=>XAS,XBS_I=>XBS,XCXS_I=>XCXS,XCS_I=>XCS,&
                               XALPHAG_I=>XALPHAG,XNUG_I=>XNUG,XDG_I=>XDG,XLBEXG_I=>XLBEXG,&
                               XLBG_I=>XLBG,XCCG_I=>XCCG,XAG_I=>XAG,XBG_I=>XBG,XCXG_I=>XCXG,XCG_I=>XCG,&
                               XALPHAH_I=>XALPHAH,XNUH_I=>XNUH,XDH_I=>XDH,XLBEXH_I=>XLBEXH,&
                               XLBH_I=>XLBH,XCCH_I=>XCCH,XAH_I=>XAH,XBH_I=>XBH,XCXH_I=>XCXH,XCH_I=>XCH,&
                               XALPHAI_I=>XALPHAI,XNUI_I=>XNUI,XDI_I=>XDI,XLBEXI_I=>XLBEXI,&
                               XLBI_I=>XLBI,XAI_I=>XAI,XBI_I=>XBI,XC_I_I=>XC_I,&
                               XRTMIN_I=>XRTMIN
!!LIMA         
USE MODD_PARAM_LIMA_WARM, ONLY: XDR_L=>XDR,XLBEXR_L=>XLBEXR,XLBR_L=>XLBR,XBR_L=>XBR,XCR_L=>XCR
USE MODD_PARAM_LIMA_COLD, ONLY: XDI_L=>XDI,XLBEXI_L=>XLBEXI,XLBI_L=>XLBI,XAI_L=>XAI,XBI_L=>XBI,XC_I_L=>XC_I,&
                                XDS_L=>XDS,XLBEXS_L=>XLBEXS,XLBS_L=>XLBS,XCCS_L=>XCCS,XNS_L=>XNS,XAS_L=>XAS,XBS_L=>XBS,&
                                XCXS_L=>XCXS,XCS_L=>XCS,XLBDAS_MIN,XLBDAS_MAX

USE MODD_PARAM_LIMA_MIXED, ONLY:XDG_L=>XDG,XLBEXG_L=>XLBEXG,XLBG_L=>XLBG,XCCG_L=>XCCG,XAG_L=>XAG,XBG_L=>XBG,XCXG_L=>XCXG,XCG_L=>XCG
USE MODD_PARAM_LIMA, ONLY: XALPHAR_L=>XALPHAR,XNUR_L=>XNUR,XALPHAS_L=>XALPHAS,XNUS_L=>XNUS,&
                           XALPHAG_L=>XALPHAG,XNUG_L=>XNUG, XALPHAI_L=>XALPHAI,XNUI_L=>XNUI,&
                           XRTMIN_L=>XRTMIN, LSNOW_T_L=>LSNOW_T
!!LIMA
USE MODD_RADAR, ONLY:XLAM_RAD,XSTEP_RAD,NBELEV,NDIFF,LATT,NPTS_GAULAG,LQUAD,XVALGROUND,NDGS, &
     LFALL,LWBSCS,LWREFL,XREFLVDOPMIN,XREFLMIN,LSNRT,XSNRMIN
USE MODD_TMAT
! 
USE MODE_ARF
USE MODE_FSCATTER
USE MODE_READTMAT
USE MODE_FGAU , ONLY:GAULAG
USE MODI_GAMMA, ONLY:GAMMA
!
USE MODD_LUNIT
USE MODE_IO_FILE,          ONLY: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_File_add2list
USE MODE_MSG

!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PT_RAY ! temperature interpolated along the rays
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PRHODREF_RAY ! 
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PR_RAY   ! rainwater mixing ratio interpolated along the rays
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PI_RAY   ! pristine ice mixing ratio interpolated along the rays
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PCIT_RAY !pristine ice concentration interpolated along the rays
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PS_RAY !aggregates mixing ratio interpolated along the rays
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PG_RAY   ! graupel mixing ratio interpolated along the rays
REAL,DIMENSION(:,:,:,:,:,:), INTENT(IN) :: PVDOP_RAY !Doppler radial velocity interpolated along the rays
REAL,DIMENSION(:,:,:,:),     INTENT(IN) :: PELEV ! elevation
REAL,DIMENSION(:),           INTENT(IN) :: PX_H ! Gaussian horizontal nodes
REAL,DIMENSION(:),           INTENT(IN) :: PX_V ! Gaussian vertical nodes
REAL,DIMENSION(:),           INTENT(IN) :: PW_H ! Gaussian horizontal weights
REAL,DIMENSION(:),           INTENT(IN) :: PW_V ! Gaussian vertical weights
REAL,DIMENSION(:,:,:,:,:),   INTENT(INOUT) :: PZE ! gate equivalent reflectivity factor (horizontal & vertical)
! convective/stratiform
REAL,DIMENSION(:,:,:,:,:,:),INTENT(INOUT) :: PBU_MASK_RAY
! /convective/stratiform
REAL, DIMENSION(:,:,:,:,:,:),OPTIONAL,INTENT(IN)  :: PCR_RAY  ! rainwater concentration interpolated along the rays
REAL, DIMENSION(:,:,:,:,:,:),OPTIONAL,INTENT(IN)  :: PH_RAY   ! hail mixing ratio interpolated along the rays
!
!*       0.2   Declarations of local variables :
!
REAL,   DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: ZREFL
!1: ZHH (dBZ), 2: ZDR, 3: KDP, 4: CSR (0 pr air clair, 1 pour stratiforme, 2 pour convectif)
!5-8: ZER, ZEI, ZES,ZEG
!9 : VRU (vitesse radiale)
!10-13 : AER, AEI, AES, AEG
!14-17: ATR, ATI, ATS, ATG
!18-20: RhoHV, PhiDP, DeltaHV

REAL,   DIMENSION(:,:,:,:,:,:,:),ALLOCATABLE :: ZAELOC ! local attenuation
REAL,   DIMENSION(:,:,:),ALLOCATABLE :: ZAETOT ! 1: total attenuation, 2: // vertical
REAL :: ZAERINT,ZAEIINT,ZAESINT,ZAEGINT,ZAEHINT ! total attenuation horizontal
REAL :: ZAVRINT,ZAVSINT,ZAVGINT,ZAVHINT ! total attenuation vertical 
!
REAL,DIMENSION(:),ALLOCATABLE :: ZX,ZW ! Gauss-Laguerre points and weights
!
REAL,DIMENSION(4) :: ZREFLOC
REAL,DIMENSION(2) :: ZAETMP
REAL,DIMENSION(:),ALLOCATABLE :: ZVTEMP ! temp var for Gaussian quadrature 8 : r_r, 9 : r_i, 10 : r_s , 11 : r_g
REAL :: ZCXR=-1.0   ! for rain N ~ 1/N_0 (in Kessler parameterization)
REAL :: ZDMELT_FACT ! factor used to compute the equivalent melted diameter
REAL :: ZEQICE=0.224! factor used to convert the ice crystals reflectivity into an equivalent  liquid water reflectivity (from Smith, JCAM 84)
REAL :: ZEXP        ! anciliary parameter
REAL :: ZLBDA   ! slope distribution parameter
REAL :: ZN      ! Number concentration
REAL :: ZFRAC_ICE,ZD,ZDE ! auxiliary variables
REAL :: ZQSCA
REAL,DIMENSION(2) :: ZQEXT
REAL,DIMENSION(3) :: ZQBACK ! Q_b(HH),Q_b(VV) (backscattering efficiencies at horizontal and vertical polarizations, resp.)
!REAL :: P=DACOS(-1D0)
REAL :: ZRHOI ! pristine ice density (from m=a*D**b), 
REAL :: ZRHOPI=916. !pure ice density (kg/m3)
COMPLEX :: ZNUM, ZDEN !for calculation of ice dielectri cconstant
COMPLEX  :: ZQM,ZQMW,ZQMI,ZQK,ZQB, ZEPSI ! dielectric parameters
REAL :: ZS11_CARRE_R,ZS22_CARRE_R,ZRE_S22S11_R,ZIM_S22S11_R
REAL :: ZS11_CARRE_I,ZS22_CARRE_I,ZRE_S22S11_I,ZIM_S22S11_I
REAL :: ZS11_CARRE_S,ZS22_CARRE_S,ZRE_S22S11_S,ZIM_S22S11_S
REAL :: ZS11_CARRE_G,ZS22_CARRE_G,ZRE_S22S11_G,ZIM_S22S11_G
REAL :: ZS11_CARRE_H,ZS22_CARRE_H,ZRE_S22S11_H,ZIM_S22S11_H
REAL :: ZS11_CARRE_T,ZS22_CARRE_T,ZRE_S22S11_T,ZIM_S22S11_T
REAL :: ZRE_S22FMS11F,ZIM_S22FT,ZIM_S11FT 

REAL :: ZM
!
INTEGER  :: INBRAD,IIELV,INBAZIM,INBSTEPMAX,INPTS_H,INPTS_V ! sizes of the arrays
INTEGER  :: IEL
INTEGER  :: JI,JL,JEL,JAZ,JH,JV,JJ,JT ! Loop variables of control
REAL :: ZLB ! depolarization factor along the spheroid symmetry axis
REAL :: ZCXI=0. ! should be defined with other parameters of microphysical scheme
REAL :: ZCR,ZCI,ZCS,ZCG,ZCH ! coefficients to take into account fall speeds when simulating Doppler winds
REAL, DIMENSION(:,:,:,:),ALLOCATABLE :: ZCONC_BIN
INTEGER :: IMAX
LOGICAL :: LPART_MASK ! indicates a partial mask along the beam

!
INTEGER :: IZER,IZEI,IZES,IZEG
INTEGER :: IVDOP,IRHV,IPDP,IDHV
INTEGER :: IAER,IAEI,IAES,IAEG
INTEGER :: IAVR,IAVI,IAVS,IAVG
INTEGER :: IATR,IATI,IATS,IATG
INTEGER :: IRHR, IRHS, IRHG, IZDA, IZDS, IZDG, IKDR, IKDS, IKDG
INTEGER :: IZEH, IRHH,IKDH,IZDH ! hail
INTEGER :: IAEH,IAVH,IATH
!
!for ZSNR threshold
REAL ::ZDISTRAD,ZSNR,ZSNR_R,ZSNR_S,ZSNR_I,ZSNR_G,ZSNR_H,ZZHH,ZZE_R,ZZE_I,ZZE_S,ZZE_G,ZZE_H
LOGICAL :: GTHRESHOLD_V, GTHRESHOLD_Z,GTHRESHOLD_ZR,GTHRESHOLD_ZI,GTHRESHOLD_ZS,GTHRESHOLD_ZG,GTHRESHOLD_ZH

!--------- TO READ T-MATRIX TABLE --------
CHARACTER(LEN=6) :: YBAND
CHARACTER(LEN=1) ::YTYPE
CHARACTER(LEN=1),DIMENSION(5) :: YTAB_TYPE
CHARACTER(LEN=25),DIMENSION(5) :: YFILE_COEFINT

REAL,DIMENSION(5) :: ZELEV_MIN,ZELEV_MAX,ZELEV_STEP,&
ZTC_MIN,ZTC_MAX,ZTC_STEP,ZFW_MIN,ZFW_MAX,ZFW_STEP
INTEGER :: IRESP,ILINE,INB_M
INTEGER,DIMENSION(5) :: INB_ELEV,INB_TC,INB_FW,INB_LINE

REAL, DIMENSION(:),ALLOCATABLE :: ZTC_T_R,        ZTC_T_S,        ZTC_T_G,        ZTC_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZELEV_T_R,      ZELEV_T_S,      ZELEV_T_G,      ZELEV_T_W
REAL, DIMENSION(:),ALLOCATABLE ::                 ZFW_T_S,        ZFW_T_G,        ZFW_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZM_T_R,         ZM_T_S,         ZM_T_G,         ZM_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZS11_CARRE_T_R, ZS11_CARRE_T_S, ZS11_CARRE_T_G, ZS11_CARRE_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZS22_CARRE_T_R, ZS22_CARRE_T_S, ZS22_CARRE_T_G, ZS22_CARRE_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZRE_S22S11_T_R, ZRE_S22S11_T_S, ZRE_S22S11_T_G, ZRE_S22S11_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZIM_S22S11_T_R, ZIM_S22S11_T_S, ZIM_S22S11_T_G, ZIM_S22S11_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZIM_S22FT_T_R,  ZIM_S22FT_T_S,  ZIM_S22FT_T_G,  ZIM_S22FT_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZIM_S11FT_T_R,  ZIM_S11FT_T_S,  ZIM_S11FT_T_G,  ZIM_S11FT_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZRE_S22FMS11FT_T_R, ZRE_S22FMS11FT_T_S, ZRE_S22FMS11FT_T_G, ZRE_S22FMS11FT_T_W
REAL, DIMENSION(:),ALLOCATABLE :: ZTC_T_H ,ZELEV_T_H ,ZFW_T_H,ZM_T_H,ZS11_CARRE_T_H,ZS22_CARRE_T_H,ZRE_S22S11_T_H
REAL, DIMENSION(:),ALLOCATABLE :: ZIM_S22S11_T_H,ZIM_S22FT_T_H,ZIM_S11FT_T_H,ZRE_S22FMS11FT_T_H
INTEGER,DIMENSION(16):: ITMAT
REAL:: ZELEV_RED,ZTC_RED,ZM_RED,ZFW_RED
INTEGER :: JIND
REAL,DIMENSION(7,16) :: KMAT_COEF !matrice contenant tous les coef interpolés
                                !pour chaque val inf et sup de ELEV_t
REAL :: ZEXPM_MIN, ZEXPM_STEP, ZEXPM_MAX,ZM_MIN
REAL :: ZFW !water fraction inside melting graupel (ZFW=0 for rain, snow and dry graupel). used only with NDIFF=7: Tmatrix
INTEGER :: ILUOUT0,IUNIT
!
! MODIF GAELLE POUR LIMA
!
LOGICAL :: GLIMA,GHAIL
REAL,DIMENSION(5) :: ZCC_MIN,ZCC_MAX, ZCC_STEP
INTEGER,DIMENSION(5):: INB_CC
REAL, DIMENSION(:),ALLOCATABLE :: ZCC_T_R
REAL :: ZCC_RED
LOGICAL :: GCALC
REAL :: ZCC
REAL,   DIMENSION(:,:,:,:,:,:),ALLOCATABLE :: ZM_6D,ZCC_6D
REAL :: ZC
!
REAL ::  ZCCR,ZLBR,ZLBEXR,ZDR,ZALPHAR,ZNUR,ZBR
REAL ::  ZCCS,ZLBS,ZLBEXS,ZDS,ZALPHAS,ZNUS,ZAS,ZBS,ZCXS,ZNS
REAL ::  ZCCG,ZLBG,ZLBEXG,ZDG,ZALPHAG,ZNUG,ZAG,ZBG,ZCXG
REAL ::  ZCCH,ZLBH,ZLBEXH,ZDH,ZALPHAH,ZNUH,ZAH,ZBH,ZCXH
REAL ::       ZLBI,ZLBEXI,ZDI,ZALPHAI,ZNUI,ZAI,ZBI
REAL,DIMENSION(:),ALLOCATABLE :: ZRTMIN
CHARACTER(LEN=100) :: YMSG
TYPE(TFILEDATA),POINTER :: TZFILE
!
!*       1.     INITIALISATION 
!--------------
ILUOUT0 = TLUOUT0%NLU
TZFILE => NULL()
! 
IF (PRESENT(PCR_RAY)) THEN
  GLIMA=.TRUE.
ELSE
  GLIMA=.FALSE.
ENDIF
IF (PRESENT(PH_RAY)) THEN
  GHAIL=.TRUE.
ELSE
  GHAIL=.FALSE.
ENDIF
!
!
!
  ZS11_CARRE_R=0
  ZS22_CARRE_R=0
  ZRE_S22S11_R=0
  ZIM_S22S11_R=0
  ZS11_CARRE_I=0
  ZS22_CARRE_I=0
  ZRE_S22S11_I=0
  ZIM_S22S11_I=0
  ZS11_CARRE_S=0
  ZS22_CARRE_S=0
  ZRE_S22S11_S=0
  ZIM_S22S11_S=0
  ZS11_CARRE_G=0
  ZS22_CARRE_G=0
  ZRE_S22S11_G=0
  ZIM_S22S11_G=0
  ZS11_CARRE_H=0
  ZS22_CARRE_H=0
  ZRE_S22S11_H=0
  ZIM_S22S11_H=0  
! Initialisation varibales microphysiques
IF (GLIMA) THEN ! LIMA
  ZLBR=XLBR_L
  ZLBEXR=XLBEXR_L
  ZDR=XDR_L
  ZALPHAR=XALPHAR_L
  ZNUR=XNUR_L
  ZBR=XBR_L
  ZCCS=XCCS_L
  ZCXS=XCXS_L
  ZLBS=XLBS_L
  ZLBEXS=XLBEXS_L
  ZNS=XNS_L
  ZDS=XDS_L
  ZALPHAS=XALPHAS_L
  ZNUS=XNUS_L
  ZAS=XAS_L
  ZBS=XBS_L
  ZCCG=XCCG_L
  ZCXG=XCXG_L
  ZLBG=XLBG_L
  ZLBEXG=XLBEXG_L
  ZDG=XDG_L
  ZALPHAG=XALPHAG_L
  ZNUG=XNUG_L
  ZAG=XAG_L
  ZBG=XBG_L
  ZLBI=XLBI_L
  ZLBEXI=XLBEXI_L
  ZDI=XDI_L
  ZALPHAI=XALPHAI_L
  ZNUI=XNUI_L
  ZAI=XAI_L
  ZBI=XBI_L
  ALLOCATE(ZRTMIN(SIZE(XRTMIN_L)))
  ZRTMIN=XRTMIN_L
ELSE ! ICE3
  ZCCR=XCCR_I
  ZLBR=XLBR_I
  ZLBEXR=XLBEXR_I
  ZDR=XDR_I
  ZALPHAR=XALPHAR_I
  ZNUR=XNUR_I
  ZBR=XBR_I
  ZCCS=XCCS_I
  ZCXS=XCXS_I
  ZLBS=XLBS_I
  ZLBEXS=XLBEXS_I
  ZNS=XNS_I
  ZDS=XDS_I
  ZALPHAS=XALPHAS_I
  ZNUS=XNUS_I
  ZAS=XAS_I
  ZBS=XBS_I
  ZCCG=XCCG_I
  ZCXG=XCXG_I
  ZLBG=XLBG_I
  ZLBEXG=XLBEXG_I
  ZDG=XDG_I
  ZALPHAG=XALPHAG_I
  ZNUG=XNUG_I
  ZAG=XAG_I
  ZBG=XBG_I
  ZLBI=XLBI_I
  ZLBEXI=XLBEXI_I
  ZDI=XDI_I
  ZALPHAI=XALPHAI_I
  ZNUI=XNUI_I
  ZAI=XAI_I
  ZBI=XBI_I
  ALLOCATE(ZRTMIN(SIZE(XRTMIN_I)))
  ZRTMIN=XRTMIN_I
  IF (GHAIL) THEN
    ZCCH=XCCH_I
    ZCXH=XCXH_I
    ZLBH=XLBH_I
    ZLBEXH=XLBEXH_I
    ZDH=XDH_I
    ZALPHAH=XALPHAH_I
    ZNUH=XNUH_I
    ZAH=XAH_I
    ZBH=XBH_I
  ENDIF
ENDIF
!
! initialisation of refractivity indices
! 1 : ZHH
! 2 : ZDR
! 3 : KDP
! 4 : CSR
IZER=5 ! ZER
IZEI=IZER+1 ! ZEI
IZES=IZEI+1 ! ZES
IZEG=IZES+1 ! ZEG
IF (GHAIL) THEN
  IZEH=IZEG+1 !ZEH
  IVDOP=IZEH+1 !VRU
ELSE
  IVDOP=IZEG+1 !VRU
END IF
IF (LATT) THEN
  IF (GHAIL) THEN
    IAER=IVDOP+1
    IAEI=IAER+1
    IAES=IAEI+1
    IAEG=IAES+1
    IAEH=IAEG+1
    IAVR=IAEH+1
    IAVI=IAVR+1
    IAVS=IAVI+1
    IAVG=IAVS+1
    IAVH=IAVG+1
    IATR=IAVH+1
    IATI=IATR+1
    IATS=IATI+1
    IATG=IATS+1
    IATH=IATG+1
    IRHV=IATH+1
  ELSE
    IAER=IVDOP+1
    IAEI=IAER+1
    IAES=IAEI+1
    IAEG=IAES+1
    IAVR=IAEG+1
    IAVI=IAVR+1
    IAVS=IAVI+1
    IAVG=IAVS+1
    IATR=IAVG+1
    IATI=IATR+1
    IATS=IATI+1
    IATG=IATS+1
    IRHV=IATG+1
  ENDIF
ELSE
    IRHV=IVDOP+1        
ENDIF
IPDP=IRHV+1
IDHV=IPDP+1
IRHR=IDHV+1
IRHS=IRHR+1
IRHG=IRHS+1
IF (GHAIL) THEN
  IRHH=IRHG+1
  IZDA=IRHH+1
ELSE
  IZDA=IRHG+1
ENDIF
IZDS=IZDA+1
IZDG=IZDS+1
IF (GHAIL) THEN
  IZDH=IZDG+1
  IKDR=IZDH+1
ELSE
  IKDR=IZDG+1
ENDIF
IKDS=IKDR+1
IKDG=IKDS+1
IF (GHAIL) THEN
  IKDH=IKDG+1
ENDIF
!
!
!
INBRAD=SIZE(PT_RAY,1)
IIELV=SIZE(PT_RAY,2)
INBAZIM=SIZE(PT_RAY,3)
INBSTEPMAX=SIZE(PT_RAY,4)
INPTS_H=SIZE(PT_RAY,5)
INPTS_V=SIZE(PT_RAY,6)
!
! Initialisation for radial winds
IF(LFALL) THEN
  IF (GLIMA) THEN
    ZCR=XCR_L
    ZCI=XC_I_L
    ZCS=XCS_L
    ZCG=XCG_L
  ELSE
    ZCR=XCR_I 
    ZCI=XC_I_I
    ZCS=XCS_I
    ZCG=XCG_I
    IF (GHAIL) ZCH=XCH_I
  ENDIF
ELSE
  ZCR=0.
  ZCI=0.
  ZCS=0.
  ZCG=0.
  IF (GHAIL)  ZCH=0.
END IF

! Calculation of nodes and weights for the Gauss-Laguerre quadrature
! for Mie and T-matrix and RG
IF(NDIFF/=0) THEN
  ALLOCATE(ZX(NPTS_GAULAG),ZW(NPTS_GAULAG)) !NPTS_GAULAG : number of points for the quadrature
  CALL GAULAG(NPTS_GAULAG,ZX,ZW)
END IF
!
!
IMAX=SIZE(PZE,5)
WRITE(ILUOUT0,*) "-----------------"
WRITE(ILUOUT0,*) "Radar scattering"
WRITE(ILUOUT0,*) "-----------------"
WRITE(ILUOUT0,*) 'Nombre de variables dans PZE: ',IMAX

IF(.NOT.LWREFL) IMAX=IMAX+1

ALLOCATE(ZREFL(INBRAD,IIELV,INBAZIM,INBSTEPMAX,INPTS_H,INPTS_V,IMAX))
ZREFL(:,:,:,:,:,:,:)=0.
IF(LATT) THEN
  ZREFL(:,:,:,:,:,:,IATR:IATG)=1.
  IF (GHAIL)  ZREFL(:,:,:,:,:,:,IATH)=1.
END IF
PZE(:,:,:,:,:)=0.
IF (LATT)THEN
  ALLOCATE(ZAELOC(INBRAD,IIELV,INBAZIM,INBSTEPMAX,INPTS_H,INPTS_V,2))
  ALLOCATE(ZAETOT(INPTS_H,INPTS_V,2))
  ZAELOC(:,:,:,:,:,:,:)=0. ! initialization of attenuation stuff (alpha_e for first gate)
  ZAETOT(:,:,:)=1. ! initialization of attenuation stuff (total attenuation)
END IF
WRITE(ILUOUT0,*) 'BEFORE LOOP DIFFUSION'

IF(LWBSCS) THEN
  ALLOCATE(ZCONC_BIN(INBRAD,IIELV,INBAZIM,INBSTEPMAX))
  ZCONC_BIN(:,:,:,:)=0.
END IF

WRITE(ILUOUT0,*) "XCCR:",ZCCR
WRITE(ILUOUT0,*) "XLBR:",ZLBR
WRITE(ILUOUT0,*) "XLBEXR:",ZLBEXR

WRITE(ILUOUT0,*) "XCCS:",ZCCS
WRITE(ILUOUT0,*) "XLBS:",ZLBS
WRITE(ILUOUT0,*) "XLBEXS:",ZLBEXS

WRITE(ILUOUT0,*) "XCCG:",ZCCG
WRITE(ILUOUT0,*) "XLBG:",ZLBG
WRITE(ILUOUT0,*) "XLBEXG:",ZLBEXG

IF (GHAIL)  THEN
  WRITE(ILUOUT0,*) "XCCH:",ZCCH
  WRITE(ILUOUT0,*) "XLBH:",ZLBH
  WRITE(ILUOUT0,*) "XLBEXH:",ZLBEXH
ENDIF
!
IF (GLIMA .AND. NDIFF==7) THEN
  IF (ZALPHAR/=1 .AND. ZNUR /=2.) THEN
    WRITE(ILUOUT0,*) " ERROR : TMATRICE TABLE ARE MADE WITH XALPHAR=1 XNUR=2"
    WRITE(ILUOUT0,*) " FOR CCLOUD=LIMA. PLEASE CHANGE THIS VALUES OR PROVIDE "
    WRITE(ILUOUT0,*) " NEW TMATRICE TABLES "
    CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SCATTERING','')
  ENDIF
ELSE
  IF (ZALPHAR/=1 .AND. ZNUR /=1.) THEN
    WRITE(ILUOUT0,*) " ERROR : TMATRICE TABLE ARE MADE WITH XALPHAR=1 XNUR=1"
    WRITE(ILUOUT0,*) " FOR CCLOUD=ICE3. PLEASE CHANGE THIS VALUEs OR PROVIDE "
    WRITE(ILUOUT0,*) " NEW TMATRICE TABLES "
    CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SCATTERING','')
  ENDIF
ENDIF

!---------------------------------------------
! LOOP OVER EVERYTHING
!--------------------------------------------
IF(NDIFF==7) THEN
  YTAB_TYPE(1)='r'
  YTAB_TYPE(2)='s'
  YTAB_TYPE(3)='g'
  YTAB_TYPE(4)='w'
  YTAB_TYPE(5)='h'
  ! definition des paramètres de lecture de la table T-matrice
  ! all mixing ratio
  ZEXPM_MIN=-7. 
  ZEXPM_STEP=0.01
  ZEXPM_MAX=-2. 
  ZM_MIN=10**ZEXPM_MIN 
  ! rain
  ZELEV_MIN(1)=0.0
  ZELEV_STEP(1)=4.0
  ZELEV_MAX(1)=12.0
  ZTC_MIN(1)=-20.0
  ZTC_STEP(1)=1.0
  ZTC_MAX(1)=40.0
  ZFW_MIN(1)=0.0
  ZFW_STEP(1)=0.1
  ZFW_MAX(1)=0.0
  IF (GLIMA) THEN
    ZCC_MIN(1)=1.8 
    ZCC_STEP(1)=0.02
    ZCC_MAX(1)=6
  ELSE
    ZCC_MIN(1)=1. 
    ZCC_STEP(1)=1.
    ZCC_MAX(1)=1.
  ENDIF
  ! snow + graupel
  ZELEV_MIN(2:3)=0.0
  ZELEV_STEP(2:3)=12.0
  ZELEV_MAX(2:3)=12.0
  ZTC_MIN(2:3)=-70.0
  ZTC_STEP(2:3)=1.0
  ZTC_MAX(2:3)=10.0
  ZFW_MIN(2:3)=0.0
  ZFW_STEP(2:3)=0.1
  ZFW_MAX(2:3)=0.0
  ZCC_MIN(2:3)=1.
  ZCC_STEP(2:3)=1.
  ZCC_MAX(2:3)=1.
  ! wet graupel
  ZELEV_MIN(4)=0.0
  ZELEV_STEP(4)=4.0
  ZELEV_MAX(4)=12.0
  ZTC_MIN(4)=-10.0
  ZTC_STEP(4)=1.0
  ZTC_MAX(4)=10.0
  ZFW_MIN(4)=0.0
  ZFW_STEP(4)=0.1
  ZFW_MAX(4)=1.0
  ZCC_MIN(4)=1.
  ZCC_STEP(4)=1.
  ZCC_MAX(4)=1. 
  ! hail
  ZELEV_MIN(5)=0.0
  ZELEV_STEP(5)=4.0
  ZELEV_MAX(5)=12.0
  ZTC_MIN(5)=-20.0
  ZTC_STEP(5)=1.0
  ZTC_MAX(5)=30.0
  ZFW_MIN(5)=0.
  ZFW_STEP(5)=0.1
  ZFW_MAX(5)=0.0
  ZCC_MIN(5)=1.
  ZCC_STEP(5)=1.
  ZCC_MAX(5)=1. 
  DO JT=1,5
    INB_ELEV(JT)=NINT((ZELEV_MAX(JT)-ZELEV_MIN(JT))/ZELEV_STEP(JT))+1
    INB_TC(JT)=NINT((ZTC_MAX(JT)-ZTC_MIN(JT))/ZTC_STEP(JT))+1
    INB_FW(JT)=NINT((ZFW_MAX(JT)-ZFW_MIN(JT))/ZFW_STEP(JT))+1
    INB_M=NINT((ZEXPM_MAX-ZEXPM_MIN)/ZEXPM_STEP)+1
    INB_CC(JT)=NINT((ZCC_MAX(JT)-ZCC_MIN(JT))/ZCC_STEP(JT))+1          
    INB_LINE(JT)=INB_ELEV(JT)*INB_TC(JT)*INB_FW(JT)*INB_M*INB_CC(JT)
  ENDDO
ENDIF                                 

!---------------------------------------------
! LOOP OVER EVERYTHING
!--------------------------------------------
 !==============         loop over radars        ================= 
WRITE(ILUOUT0,*) "INBRAD",INBRAD
DO JI=1,INBRAD
  WRITE(ILUOUT0,*) "JI",JI
  WRITE(ILUOUT0,*) "XLAM_RAD(JI):",XLAM_RAD(JI)

  IF(NDIFF==7) THEN ! If T-MATRIX
  !---------------------------------------------------------------------------------------------
  !	  0. LECTURE DES TABLES TMAT POUR PLUIE, NEIGE, GRAUPEL
  !          en fonction de la bande frequence
  !---------------------------------------------------------------------------------------------
    IF ( XLAM_RAD(JI)==0.1062) THEN
      YBAND='S106.2'     
    ELSEIF (XLAM_RAD(JI) ==0.0532 ) THEN
      YBAND='C053.2'
    ELSEIF (XLAM_RAD(JI)==0.0319 ) THEN
      YBAND='X031.9'
    ELSE
      WRITE(ILUOUT0,*) "ERROR RADAR_SCATTERING"
      WRITE(ILUOUT0,*) "Tmatrice tables are only available for XLAM_RAD=0.1062"
      WRITE(ILUOUT0,*) "or XLAM_RAD=0.0532 or XLAM_RAD=0.0319"
      WRITE(ILUOUT0,*) "change XLAM_RAD in namelist or compute new tmatrice table"
    CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SCATTERING','')
    ENDIF

    !************ fichiers Min Max Pas et Coef Tmat ***********
    DO JT=1,5  !types (r, s, g, w, h)
  
      YTYPE=YTAB_TYPE(JT)
      IF (JT .EQ. 1) THEN
        IF (GLIMA) THEN
          YFILE_COEFINT(JT)='TmatCoefInt_LIMA_'//YBAND//YTYPE
        ELSE
          YFILE_COEFINT(JT)='TmatCoefInt_ICE3_'//YBAND//YTYPE
        ENDIF
      ELSE
       YFILE_COEFINT(JT)='TmatCoefInt_'//YBAND//YTYPE
      ENDIF
      YFILE_COEFINT(JT)=TRIM(ADJUSTL(YFILE_COEFINT(JT)))  
    ENDDO
    !lookup tables for rain
    ALLOCATE (ZTC_T_R(INB_LINE(1)),ZELEV_T_R(INB_LINE(1)),ZCC_T_R(INB_LINE(1)),ZM_T_R(INB_LINE(1)),&
    ZS11_CARRE_T_R(INB_LINE(1)),ZS22_CARRE_T_R(INB_LINE(1)), ZRE_S22S11_T_R(INB_LINE(1)),ZIM_S22S11_T_R(INB_LINE(1)),&
    ZRE_S22FMS11FT_T_R(INB_LINE(1)),ZIM_S22FT_T_R(INB_LINE(1)),ZIM_S11FT_T_R(INB_LINE(1)))
   
    !lookup tables for snow
    ALLOCATE (ZTC_T_S(INB_LINE(2)),ZELEV_T_S(INB_LINE(2)),ZFW_T_S(INB_LINE(2)),ZM_T_S(INB_LINE(2)),&
    ZS11_CARRE_T_S(INB_LINE(2)),ZS22_CARRE_T_S(INB_LINE(2)),ZRE_S22S11_T_S(INB_LINE(2)),ZIM_S22S11_T_S(INB_LINE(2)),&
    ZRE_S22FMS11FT_T_S(INB_LINE(2)),ZIM_S22FT_T_S(INB_LINE(2)),ZIM_S11FT_T_S(INB_LINE(2)))

    !lookup tables for graupel
    ALLOCATE (ZTC_T_G(INB_LINE(3)),ZELEV_T_G(INB_LINE(3)),ZFW_T_G(INB_LINE(3)),ZM_T_G(INB_LINE(3)),&
    ZS11_CARRE_T_G(INB_LINE(3)),ZS22_CARRE_T_G(INB_LINE(3)), ZRE_S22S11_T_G(INB_LINE(3)),ZIM_S22S11_T_G(INB_LINE(3)),&
    ZRE_S22FMS11FT_T_G(INB_LINE(3)),ZIM_S22FT_T_G(INB_LINE(3)),ZIM_S11FT_T_G(INB_LINE(3)))

    !lookup tables for wet graupel
    ALLOCATE (ZTC_T_W(INB_LINE(4)),ZELEV_T_W(INB_LINE(4)),ZFW_T_W(INB_LINE(4)),ZM_T_W(INB_LINE(4)),&
    ZS11_CARRE_T_W(INB_LINE(4)),ZS22_CARRE_T_W(INB_LINE(4)), ZRE_S22S11_T_W(INB_LINE(4)),ZIM_S22S11_T_W(INB_LINE(4)),&
    ZRE_S22FMS11FT_T_W(INB_LINE(4)),ZIM_S22FT_T_W(INB_LINE(4)),ZIM_S11FT_T_W(INB_LINE(4)))

    IF (GHAIL) THEN
    !lookup tables for hail
    ALLOCATE (ZTC_T_H(INB_LINE(5)),ZELEV_T_H(INB_LINE(5)),ZFW_T_H(INB_LINE(5)),ZM_T_H(INB_LINE(5)),&
    ZS11_CARRE_T_H(INB_LINE(5)),ZS22_CARRE_T_H(INB_LINE(5)), ZRE_S22S11_T_H(INB_LINE(5)),ZIM_S22S11_T_H(INB_LINE(5)),&
    ZRE_S22FMS11FT_T_H(INB_LINE(5)),ZIM_S22FT_T_H(INB_LINE(5)),ZIM_S11FT_T_H(INB_LINE(5)))
    ENDIF
    !=====  Lecture des tables  ===========  
    
    6003 FORMAT (E11.4,2X,E9.3,2X,E10.4,2X,E10.4,2X,E12.5,2X,E12.5,2X,&
                 E12.5,2X,E12.5,2X,E12.5,2X,E12.5,2X,E12.5)

    !rain
    CALL IO_File_add2list(TZFILE,YFILE_COEFINT(1),'TXT','READ')
    CALL IO_File_open(TZFILE,KRESP=IRESP)
    IUNIT = TZFILE%NLU
    IF ( IRESP /= 0 ) THEN       
      WRITE(YMSG,*) "problem opening file ",TRIM(YFILE_COEFINT(1))
      CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SCATTERING',YMSG)
    ENDIF
    ILINE=1
    DO WHILE (ILINE .LE. INB_LINE(1))
      READ( UNIT=IUNIT,FMT=6003, IOSTAT=IRESP ) ZTC_T_R(ILINE),ZELEV_T_R(ILINE),&
      ZCC_T_R(ILINE),ZM_T_R(ILINE),ZS11_CARRE_T_R(ILINE),ZS22_CARRE_T_R(ILINE),ZRE_S22S11_T_R(ILINE),&
      ZIM_S22S11_T_R(ILINE),ZRE_S22FMS11FT_T_R(ILINE),ZIM_S22FT_T_R(ILINE),ZIM_S11FT_T_R(ILINE) 
      ILINE=ILINE+1
    ENDDO
    CALL IO_File_close(TZFILE)
    TZFILE => NULL()    
    WRITE(ILUOUT0,*) "NLIGNE rain",ILINE      
    ILINE=2
    WRITE(ILUOUT0,*) "ILINE=",ILINE
    WRITE(ILUOUT0,*) "ZTC_T_R(ILINE),ZELEV_T_R(ILINE),ZCC_T_R(ILINE)",&
    ZTC_T_R(ILINE),ZELEV_T_R(ILINE),ZCC_T_R(ILINE)
    WRITE(ILUOUT0,*) "ZM_T_R(ILINE),ZS11_CARRE_T_R(ILINE),ZS22_CARRE_T_R(ILINE),ZRE_S22S11_T_R(ILINE)",&
    ZM_T_R(ILINE),ZS11_CARRE_T_R(ILINE),ZS22_CARRE_T_R(ILINE),ZRE_S22S11_T_R(ILINE)
    WRITE(ILUOUT0,*) "ZIM_S22S11_T_R(ILINE),ZRE_S22FMS11FT_T_R(ILINE),ZIM_S22FT_T_R(ILINE),ZIM_S11FT_T_R(ILINE)",&
    ZIM_S22S11_T_R(ILINE),ZRE_S22FMS11FT_T_R(ILINE),ZIM_S22FT_T_R(ILINE),ZIM_S11FT_T_R(ILINE)
   
    !snow
    CALL IO_File_add2list(TZFILE,YFILE_COEFINT(2),'TXT','READ')
    CALL IO_File_open(TZFILE,KRESP=IRESP)
    IUNIT = TZFILE%NLU
    IF ( IRESP /= 0 ) THEN       
      WRITE(YMSG,*) "problem opening file ",TRIM(YFILE_COEFINT(2))
      CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SCATTERING',YMSG)
    ENDIF
    ILINE=1
    DO WHILE (ILINE .LE. INB_LINE(2))
      READ( UNIT=IUNIT,FMT=6003, IOSTAT=IRESP ) ZTC_T_S(ILINE),ZELEV_T_S(ILINE),&
      ZFW_T_S(ILINE),ZM_T_S(ILINE),ZS11_CARRE_T_S(ILINE),ZS22_CARRE_T_S(ILINE),ZRE_S22S11_T_S(ILINE),&
      ZIM_S22S11_T_S(ILINE),ZRE_S22FMS11FT_T_S(ILINE),ZIM_S22FT_T_S(ILINE),ZIM_S11FT_T_S(ILINE)
      ILINE=ILINE+1
    ENDDO
    CALL IO_File_close(TZFILE)
    TZFILE => NULL()
    WRITE(ILUOUT0,*) "NLIGNE snow",ILINE    
    ILINE=2
    WRITE(ILUOUT0,*) "ILINE=",ILINE
    WRITE(ILUOUT0,*) "ZTC_T_S(ILINE),ZELEV_T_S(ILINE),ZFW_T_S(ILINE)",&
    ZTC_T_S(ILINE),ZELEV_T_S(ILINE),ZFW_T_S(ILINE)
    WRITE(ILUOUT0,*) "ZM_T_S(ILINE),ZS11_CARRE_T_S(ILINE),ZS22_CARRE_T_S(ILINE),ZRE_S22S11_T_S(ILINE)",&
    ZM_T_S(ILINE),ZS11_CARRE_T_S(ILINE),ZS22_CARRE_T_S(ILINE),ZRE_S22S11_T_S(ILINE)
    WRITE(ILUOUT0,*) "ZIM_S22S11_T_S(ILINE),ZRE_S22FMS11FT_T_S(ILINE),ZIM_S22FT_T_S(ILINE),ZIM_S11FT_T_S(ILINE)",&
    ZIM_S22S11_T_S(ILINE),ZRE_S22FMS11FT_T_S(ILINE),ZIM_S22FT_T_S(ILINE),ZIM_S11FT_T_S(ILINE)
 
    !graupel
    CALL IO_File_add2list(TZFILE,YFILE_COEFINT(3),'TXT','READ')
    CALL IO_File_open(TZFILE,KRESP=IRESP)
    IUNIT = TZFILE%NLU
    IF ( IRESP /= 0 ) THEN       
      WRITE(YMSG,*) "problem opening file ",TRIM(YFILE_COEFINT(3))
      CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SCATTERING',YMSG)
    ENDIF
    ILINE=1
    DO WHILE (ILINE .LE. INB_LINE(3))
      READ( UNIT=IUNIT, FMT=6003,IOSTAT=IRESP ) ZTC_T_G(ILINE),ZELEV_T_G(ILINE),&
      ZFW_T_G(ILINE),ZM_T_G(ILINE),ZS11_CARRE_T_G(ILINE),ZS22_CARRE_T_G(ILINE),ZRE_S22S11_T_G(ILINE),&
      ZIM_S22S11_T_G(ILINE),ZRE_S22FMS11FT_T_G(ILINE),ZIM_S22FT_T_G(ILINE),ZIM_S11FT_T_G(ILINE)
      ILINE=ILINE+1
    ENDDO
    CALL IO_File_close(TZFILE)
    TZFILE => NULL()
    WRITE(ILUOUT0,*) "NLIGNE graupel",ILINE    
    ILINE=2
    WRITE(ILUOUT0,*) "ILINE=",ILINE
    WRITE(ILUOUT0,*) "ZTC_T_G(ILINE),ZELEV_T_G(ILINE)",&
    ZTC_T_G(ILINE),ZELEV_T_G(ILINE)
    WRITE(ILUOUT0,*) "ZM_T_G(ILINE),ZS11_CARRE_T_G(ILINE),ZS22_CARRE_T_G(ILINE),ZRE_S22S11_T_G(ILINE)",&
    ZM_T_G(ILINE),ZS11_CARRE_T_G(ILINE),ZS22_CARRE_T_G(ILINE),ZRE_S22S11_T_G(ILINE)
    WRITE(ILUOUT0,*) "ZIM_S22S11_T_G(ILINE),ZRE_S22FMS11FT_T_G(ILINE),ZIM_S22FT_T_G(ILINE),ZIM_S11FT_T_G(ILINE)",&
    ZIM_S22S11_T_G(ILINE),ZRE_S22FMS11FT_T_G(ILINE),ZIM_S22FT_T_G(ILINE),ZIM_S11FT_T_G(ILINE)

    !wet graupel
    CALL IO_File_add2list(TZFILE,YFILE_COEFINT(4),'TXT','READ')
    CALL IO_File_open(TZFILE,KRESP=IRESP)
    IUNIT = TZFILE%NLU
    IF ( IRESP /= 0 ) THEN       
      WRITE(YMSG,*) "problem opening file ",TRIM(YFILE_COEFINT(4))
      CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SCATTERING',YMSG)
    ENDIF
    ILINE=1
    DO WHILE (ILINE .LE. INB_LINE(4))
      READ( UNIT=IUNIT, FMT=6003,IOSTAT=IRESP ) ZTC_T_W(ILINE),ZELEV_T_W(ILINE),&
      ZFW_T_W(ILINE),ZM_T_W(ILINE),ZS11_CARRE_T_W(ILINE),ZS22_CARRE_T_W(ILINE),ZRE_S22S11_T_W(ILINE),&
      ZIM_S22S11_T_W(ILINE),ZRE_S22FMS11FT_T_W(ILINE),ZIM_S22FT_T_W(ILINE),ZIM_S11FT_T_W(ILINE)
      ILINE=ILINE+1
    ENDDO
    CALL IO_File_close(TZFILE)
    TZFILE => NULL()
    WRITE(ILUOUT0,*) "NLIGNE wet graupel",ILINE    
    ILINE=2
    WRITE(ILUOUT0,*) "ILINE=",ILINE
    WRITE(ILUOUT0,*) "ZTC_T_W(ILINE),ZELEV_T_W(ILINE)", ZTC_T_W(ILINE),ZELEV_T_W(ILINE)
    WRITE(ILUOUT0,*) "ZM_T_W(ILINE),ZS11_CARRE_T_W(ILINE),ZS22_CARRE_T_W(ILINE),ZRE_S22S11_T_W(ILINE)",&
    ZM_T_W(ILINE),ZS11_CARRE_T_W(ILINE),ZS22_CARRE_T_W(ILINE),ZRE_S22S11_T_W(ILINE)
    WRITE(ILUOUT0,*) "ZIM_S22S11_T_W(ILINE),ZRE_S22FMS11FT_T_W(ILINE),ZIM_S22FT_T_W(ILINE),ZIM_S11FT_T_W(ILINE)",&
    ZIM_S22S11_T_W(ILINE),ZRE_S22FMS11FT_T_W(ILINE),ZIM_S22FT_T_W(ILINE),ZIM_S11FT_T_W(ILINE)

    !hail
    IF (GHAIL) THEN
    CALL IO_File_add2list(TZFILE,YFILE_COEFINT(5),'TXT','READ')
    CALL IO_File_open(TZFILE,KRESP=IRESP)
    IUNIT = TZFILE%NLU
    IF ( IRESP /= 0 ) THEN       
      WRITE(YMSG,*) "problem opening file ",TRIM(YFILE_COEFINT(5))
      CALL PRINT_MSG(NVERB_FATAL,'GEN','RADAR_SCATTERING',YMSG)
    ENDIF
      ILINE=1
      DO WHILE (ILINE .LE. INB_LINE(5))
        READ( UNIT=IUNIT, FMT=6003,IOSTAT=IRESP ) ZTC_T_H(ILINE),ZELEV_T_H(ILINE),&
        ZFW_T_H(ILINE),ZM_T_H(ILINE),ZS11_CARRE_T_H(ILINE),ZS22_CARRE_T_H(ILINE),ZRE_S22S11_T_H(ILINE),&
        ZIM_S22S11_T_H(ILINE),ZRE_S22FMS11FT_T_H(ILINE),ZIM_S22FT_T_H(ILINE),ZIM_S11FT_T_H(ILINE)
        ILINE=ILINE+1
      ENDDO
      CALL IO_File_close(TZFILE)
      TZFILE => NULL()
      WRITE(ILUOUT0,*) "NLIGNE hail",ILINE    
      ILINE=2
      WRITE(ILUOUT0,*) "ILINE=",ILINE
      WRITE(ILUOUT0,*) "ZTC_T_H(ILINE),ZELEV_T_H(ILINE)", ZTC_T_H(ILINE),ZELEV_T_H(ILINE)
      WRITE(ILUOUT0,*) "ZM_T_H(ILINE),ZS11_CARRE_T_H(ILINE),ZS22_CARRE_T_H(ILINE),ZRE_S22S11_T_H(ILINE)",&
      ZM_T_W(ILINE),ZS11_CARRE_T_H(ILINE),ZS22_CARRE_T_H(ILINE),ZRE_S22S11_T_H(ILINE)
      WRITE(ILUOUT0,*) "ZIM_S22S11_T_H(ILINE),ZRE_S22FMS11FT_T_H(ILINE),ZIM_S22FT_T_H(ILINE),ZIM_S11FT_T_H(ILINE)",&
      ZIM_S22S11_T_H(ILINE),ZRE_S22FMS11FT_T_H(ILINE),ZIM_S22FT_T_H(ILINE),ZIM_S11FT_T_H(ILINE)
    ENDIF
  ENDIF !END IF T-MATRIX => END OF LOOKUP TABLE READING

 !==============         loop over elevations      ================= 
  IEL=NBELEV(JI)
  WRITE(ILUOUT0,*) "NBELEV(JI)",NBELEV(JI)
  WRITE(ILUOUT0,*) "INPTS_V",INPTS_V
  DO JEL=1,IEL  
    WRITE(ILUOUT0,*) "JEL",JEL
    JL=1
    JV=1
    WRITE(ILUOUT0,*) "JL,JV",JL,JV
    WRITE(ILUOUT0,*) "PELEV(JI,JEL,JL,JV)*180./XPI",PELEV(JI,JEL,JL,JV)*180./XPI
    JL=INBSTEPMAX
    JV=INPTS_V
    WRITE(ILUOUT0,*) "JL,JV",JL,JV
    WRITE(ILUOUT0,*) "PELEV(JI,JEL,JL,JV)*180./XPI",PELEV(JI,JEL,JL,JV)*180./XPI
    !==============         loop over azimuths     ================= 
    DO JAZ=1,INBAZIM   
      DO JH=1,INPTS_H !horizontal discretization of the beam
        DO JV=1,INPTS_V ! vertical discretization (we go down to check partial masks)
          IF(LATT) THEN
            ZAERINT=1.
            ZAVRINT=1.
            ZAEIINT=1.
            ZAESINT=1.
            ZAVSINT=1.
            ZAEGINT=1.
            ZAVGINT=1.
            ZAEHINT=1.
            ZAVHINT=1.
          END IF
          !Loop over the ranges for one azimuth. If the range is masked, the reflectivity for all the consecutive ranges is set to 0
          LPART_MASK=.FALSE.
          LOOPJL: DO JL=1,INBSTEPMAX
            IF(LPART_MASK) THEN ! THIS RAY IS MASKED
              ZREFL(JI,JEL,JAZ,JL:INBSTEPMAX,JH,JV,1)=0. 
              EXIT LOOPJL
            ELSE
              ! if not underground or outside of the MESO-NH domain (PT_RAY : temperature interpolated along the rays)             
              IF(PT_RAY(JI,JEL,JAZ,JL,JH,JV) /= -XUNDEF) THEN
                !
                !---------------------------------------------------------------------------------------------------
                !*       2.    RAINDROPS
                !              ---------
                !
                IF(SIZE(PR_RAY,1) > 0) THEN ! if PR_RAY is available for at least one radar 
                  !contenu en hydrometeore 
                  ZM=PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PR_RAY(JI,JEL,JAZ,JL,JH,JV)
                  IF (GLIMA) ZCC=PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PCR_RAY(JI,JEL,JAZ,JL,JH,JV)
                  !ZM_MIN : min value for rain content (10**-7 <=> Z=-26 dBZ)mixing ratio
                  IF (GLIMA) THEN
                    GCALC=((ZM > ZM_MIN).AND.(ZCC > 10**ZCC_MIN(1)))  
                  ELSE
                    GCALC=(ZM > ZM_MIN)
                  ENDIF
                  IF(GCALC ) THEN 
                    !calculation of the dielectrique constant (permittitivité relative)
                    ! for liquid water from function QEPSW
                    !(defined in mode_fscatter.f90 => equation 3.6 p 64)
                    YTYPE='r'  
                    ZQMW=SQRT(QEPSW(PT_RAY(JI,JEL,JAZ,JL,JH,JV),XLIGHTSPEED/XLAM_RAD(JI)))
                    !ZLBDA : slope distribution parameter (equation 2.6 p 23)  
                     IF (GLIMA) THEN
                       ZLBDA=( ZLBR*ZCC / ZM )**ZLBEXR
                     ELSE
                       ZLBDA=ZLBR*(ZM)**ZLBEXR
                     ENDIF
                    ZQK=(ZQMW**2-1.)/(ZQMW**2+2.) !dielectric factor (3.43 p 56)
                    ZFW=0 !Liquid water fraction (only for melting graupel => 0 for rain) 
        
                    !compteur=compteur+1
                    !---------------------------------------------------  
                    ! ------------    DIFFUSION      --------------
                    !---------------------------------------------------  
                    !******************************* NDIFF=0 or 4 *********************************
                    IF(NDIFF==0.OR.NDIFF==4) THEN ! Rayleigh
                      !ZREFLOC(1:2) : Zh et Zv = int(sigma(D)*N(D)) (eq 1.6 p 16)
                      !with N(D) formulation (eq 2.2 p 23) and sigma Rayleigh (3.41 p 55)
                      !MOMG : gamma function defined in mong.f90
                      !XCCR = 1.E7; XLBEXR = -0.25! Marshall-Palmer law (radar_rain_ice.f90)
                      !ZCXR : -1 (Xi coeff in equation 2.3 p 23)
                      ZREFLOC(1:2)=1.E18*ZCCR*ZLBDA**(ZCXR-6.)*MOMG(ZALPHAR,ZNUR,6.)
                      IF(LWREFL) THEN ! weighting by reflectivities
                        !ZREFL(...,IVDOP)=radial velocity (IVDOP=9), weighted by reflectivity and
                        !taking into account raindrops fall velocity (ZCR = 842, XDR = 0.8 -> 2.8 p23 et 2.1 p24)
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=-ZCR*SIN(PELEV(JI,JEL,JL,JV)) &
                              *1.E18*ZCCR*ZLBDA**(ZCXR-6.-ZDR)*MOMG(ZALPHAR,ZNUR,6.+ZDR)
                      ELSE 
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)=ZCCR*ZLBDA**ZCXR ! N0j of equation 2.3 p23 (density of particules)
                        !projection of fall velocity only
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=-ZCR*SIN(PELEV(JI,JEL,JL,JV)) &
                                                         *ZCCR*ZLBDA**(ZCXR-ZDR)*MOMG(ZALPHAR,ZNUR,ZDR)
                      END IF ! end weighting by reflectivities
                      IF(LATT) THEN ! Calculation of Extinction coefficient
                        IF(NDIFF==0) THEN ! Rayleigh 3rd order : calculation from equations
                          ! 3.39 p55 : extinction coeff = int(extinction_section(D) * N(D))
                          ! 2.2 and 2.3 p23: simplification of int(D**p * N(D)) and N0j
                          ! 3.42 p57 : extinction_section(D)
                          ZAETMP(:)=ZCCR*ZLBDA**ZCXR*(XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)&
                                    *MOMG(ZALPHAR,ZNUR,ZBR)/ZLBDA**ZBR)
                        ELSE  ! Rayleigh 6th order ! eq 3.52 p 58 for extinction coefficient
                          ZAETMP(:)=ZCCR*ZLBDA**ZCXR*(XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)&
                                    *MOMG(ZALPHAR,ZNUR,ZBR)/ZLBDA**ZBR             &
                                    +XPI**4/15./XLAM_RAD(JI)**3*AIMAG(ZQK**2*(ZQMW**4+27.*ZQMW**2+38.) &
                                    /(2.*ZQMW**2+3.))*MOMG(ZALPHAR,ZNUR,5.*ZBR/3.)/ZLBDA**(5.*ZBR/3.)&
                                    +2.*XPI**5/3. /XLAM_RAD(JI)**4*REAL(ZQK**2)                      &
                                    *MOMG(ZALPHAR,ZNUR,2.*ZBR)   /ZLBDA**(2.*ZBR))
                        END IF
                      END IF ! end IF(LATT)
                      ZRE_S22S11_R=0
                      ZIM_S22S11_R=0
                      ZS22_CARRE_R=0
                      ZS11_CARRE_R=0       
                    !******************************* NDIFF==7 ************************************                 
                    ELSE IF(NDIFF==7) THEN !T-matrix 
                      ZREFLOC(:)=0
                      IF(LATT) ZAETMP(:)=0       
                      IF (GLIMA) THEN
                        CALL CALC_KTMAT_LIMA(PELEV(JI,JEL,JL,JV),&
                                        PT_RAY(JI,JEL,JAZ,JL,JH,JV),ZCC,ZM,&
                                        ZELEV_MIN(1),ZELEV_MAX(1),ZELEV_STEP(1),&
                                        ZTC_MIN(1),ZTC_MAX(1),ZTC_STEP(1),&
                                        ZCC_MIN(1),ZCC_MAX(1),ZCC_STEP(1),&
                                        ZEXPM_MIN,ZEXPM_MAX,ZEXPM_STEP,&
                                        ITMAT,ZELEV_RED,ZTC_RED,ZCC_RED,ZM_RED)
                      ELSE
                        CALL CALC_KTMAT(PELEV(JI,JEL,JL,JV),&
                                        PT_RAY(JI,JEL,JAZ,JL,JH,JV),ZFW,ZM,&
                                        ZELEV_MIN(1),ZELEV_MAX(1),ZELEV_STEP(1),&
                                        ZTC_MIN(1),ZTC_MAX(1),ZTC_STEP(1),&
                                        ZFW_MIN(1),ZFW_MAX(1),ZFW_STEP(1),&
                                        ZEXPM_MIN,ZEXPM_MAX,ZEXPM_STEP,&
                                        ITMAT,ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED)
                      ENDIF
                      IF (ITMAT(1) .NE. -NUNDEF) THEN 
                        DO JIND=1,SIZE(KMAT_COEF,2),1
                          KMAT_COEF(1,JIND)=ZS11_CARRE_T_R(ITMAT(JIND))
                          KMAT_COEF(2,JIND)=ZS22_CARRE_T_R(ITMAT(JIND))
                          KMAT_COEF(3,JIND)=ZRE_S22S11_T_R(ITMAT(JIND))
                          KMAT_COEF(4,JIND)=ZIM_S22S11_T_R(ITMAT(JIND))
                          KMAT_COEF(5,JIND)=ZRE_S22FMS11FT_T_R(ITMAT(JIND))
                          KMAT_COEF(6,JIND)=ZIM_S22FT_T_R(ITMAT(JIND))
                          KMAT_COEF(7,JIND)=ZIM_S11FT_T_R(ITMAT(JIND))
                        ENDDO
                        IF (GLIMA) THEN
                          CALL INTERPOL(ZELEV_RED,ZTC_RED,ZCC_RED,ZM_RED,KMAT_COEF,ZS11_CARRE_R,ZS22_CARRE_R,&
                                       ZRE_S22S11_R,ZIM_S22S11_R,ZRE_S22FMS11F,ZIM_S22FT,ZIM_S11FT)  
                        ELSE
                          CALL INTERPOL(ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED,KMAT_COEF,ZS11_CARRE_R,ZS22_CARRE_R,&
                                       ZRE_S22S11_R,ZIM_S22S11_R,ZRE_S22FMS11F,ZIM_S22FT,ZIM_S11FT)
                        ENDIF
                      ELSE
                        ZS11_CARRE_R=0
                        ZS22_CARRE_R=0
                        ZRE_S22S11_R=0
                        ZIM_S22S11_R=0
                        ZRE_S22FMS11F=0
                        ZIM_S22FT=0
                        ZIM_S11FT=0
                      END IF
                      ZREFLOC(1)=1.E18*(XLAM_RAD(JI))**4/(XPI**5*.93)*4*XPI*ZS22_CARRE_R
                      ZREFLOC(2)=1.E18*(XLAM_RAD(JI))**4/(XPI**5*.93)*4*XPI*ZS11_CARRE_R
                      ZREFLOC(3)=180.E3/XPI*XLAM_RAD(JI)*ZRE_S22FMS11F
                      IF (GLIMA) THEN
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1) &
                                                       -ZCR*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(1) &
                                                       *1.E18*(XLAM_RAD(JI)/XPI)**4/.93*ZCC/4./ZLBDA**(2+ZDR)
                      ELSE
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1) &
                                                       -ZCR*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(1) &
                                                       *1.E18*(XLAM_RAD(JI)/XPI)**4/.93*ZCCR/4./ZLBDA**(3+ZDR)
                      ENDIF
                      IF(LATT) THEN
                        ZAETMP(1)=ZIM_S22FT*XLAM_RAD(JI)*2
                        ZAETMP(2)=ZIM_S11FT*XLAM_RAD(JI)*2
                      END IF   
                    !******************************* NDIFF=1 or 3 *********************************
                    !       Gauss Laguerre integration
                    ELSE ! MIE OR T-MATRIX OR RAYLEIGH FOR ELLIPSOIDES
                      ZREFLOC(:)=0.
                      IF(LATT) ZAETMP(:)=0.
                      DO JJ=1,NPTS_GAULAG ! ****** Gauss-Laguerre quadrature
                        SELECT CASE(NDIFF)
                        CASE(1) ! *************** NDIFF=1 MIE *****************
                          ! subroutine BHMIE defined in mode_fscatter.f90
                          ! calculate extinction coefficient ZQEXT(1),scattering : ZQSCA
                          ! and backscattering ZQBACK(1) on the horizontal plan (spheroid)
                          CALL BHMIE(XPI/XLAM_RAD(JI)*ZX(JJ)/ZLBDA,ZQMW,ZQEXT(1),ZQSCA,ZQBACK(1))
                          ZQBACK(2)=ZQBACK(1) !=> same because sphere
                          ZQEXT(2)=ZQEXT(1) ! modif Clotilde 23/04/2012
                          ZQBACK(3)=0. !=> 0 because sphere
                        CASE(3) !****************** NDIFF==3 RG RAYLEIGH FOR ELLIPSOIDES ***********************
                          IF(ARF(ZX(JJ)/ZLBDA)==1.) THEN
                            ZLB=1./3.
                          ELSE
                            ZLB=1./(ARF(ZX(JJ)/ZLBDA))**2-1. ! f**2
                            ZLB=(1.+ZLB)/ZLB*(1.-ATAN(SQRT(ZLB))/SQRT(ZLB)) ! lambda_b
                            IF(ZX(JJ)/ZLBDA>16.61E-3) PRINT*, 'Negative axis ratio; reduce NPTS_GAULAG.'
                          END IF
                          ! equation 3.44 p 56 (ZX**4 instead of ZX**6 but ZQBACK is multiplied after by ZX**2)
                          ZQBACK(1)=4.*(XPI/XLAM_RAD(JI)*ZX(JJ)/ZLBDA)**4&
                                    *ABS((ZQMW**2-1.)/3./(1.+.5*(1.-ZLB)*(ZQMW**2-1.)))**2
                          ! equation 3.45 p 56
                          ZQBACK(2)=4.*(XPI/XLAM_RAD(JI)*ZX(JJ)/ZLBDA)**4*ABS((ZQMW**2-1.)/3.*&
                                    (SIN(PELEV(JI,JEL,JL,JV))**2/(1.+.5*(1.-ZLB)*(ZQMW**2-1.))+& ! PELEV=PI+THETA_I
                                    COS(PELEV(JI,JEL,JL,JV))**2/(1.+ZLB*(ZQMW**2-1.))) )**2 !
                          ! KDP from equation 3.49
                          ZQBACK(3)=ZX(JJ)/ZLBDA**3*REAL((ZQMW**2-1.)**2*(3.*ZLB-1.)/(2.+(ZQMW**2-1.)*(ZLB+1.) &
                          +ZLB*(1.-ZLB)*(ZQMW**2-1.)**2))
                          IF(LATT) THEN
                            ! equations 3.48 and 3.49 p57
                            ZQEXT(1)=4.*(XPI/XLAM_RAD(JI)*ZX(JJ)/ZLBDA)*AIMAG((ZQMW**2-1.)/3./(1.+.5*(1.-ZLB)*(ZQMW**2-1.)))
                            ZQEXT(2)=4.*(XPI/XLAM_RAD(JI)*ZX(JJ)/ZLBDA)*AIMAG((ZQMW**2-1.)/3.*&
                                     (SIN(PELEV(JI,JEL,JL,JV))**2/(1.+.5*(1.-ZLB)*(ZQMW**2-1.))+& ! PELEV=PI+THETA_I
                                     COS(PELEV(JI,JEL,JL,JV))**2/(1.+ZLB*(ZQMW**2-1.))))
                          END IF
                        END SELECT !end SELECT NDIFF
                        !incrementation of the reflectivity and Kdp(1,2,3,4 for Zh, Zv, )
                        !with the backscattering coefficients for each point of the GAULAG distribution
                        ! or each diameter D
                        ZREFLOC(1:3)=ZREFLOC(1:3)+ZQBACK(1:3)*ZX(JJ)**2*ZW(JJ)
                        ZREFLOC(4)=ZREFLOC(4)+ZQBACK(1)*ZX(JJ)**(2+ZDR)*ZW(JJ)
                          !same for attenuation with extinction coefficient						
                        IF(LATT) ZAETMP(:)=ZAETMP(:)+ZQEXT(:)*ZX(JJ)**2*ZW(JJ)
                      END DO ! ****** end loop Gauss-Laguerre quadrature
               
                      ZREFLOC(1:2)=1.E18*ZREFLOC(1:2)*(XLAM_RAD(JI)/XPI)**4/.93*ZCCR/4./ZLBDA**3
                      ZREFLOC(3)=ZREFLOC(3)*XPI**2/6./XLAM_RAD(JI)*ZCCR/ZLBDA &
                                 *180.E3/XPI ! (in deg/km)		
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1) &
                                                       -ZCR*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(4) &
                                                       *1.E18*(XLAM_RAD(JI)/XPI)**4/.93*ZCCR/4./ZLBDA**(3+ZDR)
          
                      !********* for all cases with Gauss-Laguerre integration
                      ZRE_S22S11_R=0
                      ZIM_S22S11_R=0
                      ZS22_CARRE_R=0
                      ZS11_CARRE_R=0
                      IF(LATT) ZAETMP(:)=ZAETMP(:)*XPI*ZCCR*ZLBDA**(ZCXR-2.*ZBR/3.)/(4.*GAMMA(ZNUR))
                    END IF ! ****************** End if for each type of diffusion ************************
                    !incrementation of ZHH, ZDR and KDP
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)=ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)+ZREFLOC(1:3)
                    ! ZER (Z due to raindrops)
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZER)=ZREFLOC(1)
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDA)=ZREFLOC(2) !Zvv for ZDR due to rain
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IKDR)=ZREFLOC(3) !Zvv for ZDR due to rain  

                    ! RhoHV due to rain
                    IF (ZS22_CARRE_R*ZS11_CARRE_R .GT. 0) THEN
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHR)=SQRT(ZRE_S22S11_R**2+ZIM_S22S11_R**2)/SQRT(ZS22_CARRE_R*ZS11_CARRE_R)
                    ELSE
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHR)=1
                    END IF
                    IF(LATT) THEN
                      ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)=ZAETMP(:) ! specific attenuation due to rain
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAER)=ZAETMP(1)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAVR)=ZAETMP(2)
                      ! for ranges over 1, correction of attenuation on reflectivity due to rain
                      IF(JL>1) THEN
                        ZAERINT=ZAERINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAER)*XSTEP_RAD)
                        ZAVRINT=ZAVRINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAVR)*XSTEP_RAD)
                      END IF
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZER)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZER)*ZAERINT ! Z_r attenuated
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDA)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDA)*ZAVRINT ! ZVr attenuated		
                    END IF !end IF(LATT)
                  END IF 
                        ! mimimum rainwater mixing ratio
                        ! Total attenuation even if no hydrometeors (equation 1.7 p 17)
                  IF(LATT.AND.JL>1) ZREFL(JI,JEL,JAZ,JL,JH,JV,IATR)=ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IATR) &
                                          *EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAER)*XSTEP_RAD)
                END IF ! **************** end RAIN (end IF SIZE(PR_RAY,1) > 0)
                !
                !---------------------------------------------------------------------------------------------------
                !*       3.    PRISTINE ICE
                !              ---------
                !
                IF (SIZE(PI_RAY,1)>0)  THEN
                  ZM=PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PI_RAY(JI,JEL,JAZ,JL,JH,JV) !ice content
                  IF (PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)==-XUNDEF .OR. PI_RAY(JI,JEL,JAZ,JL,JH,JV)==-XUNDEF) ZM=-XUNDEF
                  IF (GLIMA) THEN
                    ZC=PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PCIT_RAY(JI,JEL,JAZ,JL,JH,JV)
                    IF (PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)==-XUNDEF .OR. PCIT_RAY(JI,JEL,JAZ,JL,JH,JV)==-XUNDEF) ZC=-XUNDEF
                  ELSE
                    ZC=PCIT_RAY(JI,JEL,JAZ,JL,JH,JV)
                  ENDIF
                  IF(ZM>ZM_MIN .AND. ZC> 527.82) THEN  
                    ! cit > 527.82 otherwise pbs due to interpolation
                    !ice dielectric constant (QPESI defined in mode_fscatter, equation 3.65 p 65)
                    ZEPSI=QEPSI(PT_RAY(JI,JEL,JAZ,JL,JH,JV),XLIGHTSPEED/XLAM_RAD(JI))
                    ZQMI=SQRT(ZEPSI)
                    ZQK=(ZQMI**2-1.)/(ZQMI**2+2.)
                    !see 3.77 p68 : to replace Dg by an equivalent diameter De of pure ice, a multiplicative
                    !melting factor has to be added                      
                    ZDMELT_FACT=(6.*ZAI)/(XPI*.92*XRHOLW)
                    ZEXP=2.*ZBI !XBI = 2.5 (Plates) in ini_radar.f90 (bj tab 2.1 p24)
                    !ZLBDA : slope distribution parameter (equation 2.6 p 23)
                    IF (GLIMA) THEN                   
                      ZLBDA=(ZLBI*ZC/ZM)**ZLBEXI
                    ELSE
                      ZLBDA=ZLBI*(ZM/ZC)**ZLBEXI
                    ENDIF
                    ! Rayleigh or Rayleigh-Gans (=> Rayleigh) or Rayleigh with 6th order for attenuation
                    ! (pristine ice = sphere),
                    IF(NDIFF==0.OR.NDIFF==3.OR.NDIFF==4) THEN
                      !ZREFLOC(1:2) : Zh et Zv from equation 2.2 p23 and Cristals parameters
                      !ZEQICE=0.224 (radar_rain_ice.f90) factor used to convert the ice crystals
                      !reflectivity into an equivalent liquid water reflectivity (from Smith, JCAM 84)
                      ZREFLOC(1:2)=ZEQICE*.92**2*ZDMELT_FACT**2*1.E18*ZC &
                                   *ZLBDA**(ZCXI-ZEXP)*MOMG(ZALPHAI,ZNUI,ZEXP)
                      ZREFLOC(3)=0.
                      IF(LWREFL) THEN ! weighting by reflectivities
                        !calculation of radial velocity
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                               -ZCI*SIN(PELEV(JI,JEL,JL,JV))*ZEQICE*.92**2*ZDMELT_FACT**2&
                               *1.E18*ZC*ZLBDA**(ZCXI-ZEXP-ZDI)&
                               *MOMG(ZALPHAI,ZNUI,ZEXP+ZDI)
                      ELSE
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)&
                                              +ZC*ZLBDA**ZCXI
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                                                            -ZCI*SIN(PELEV(JI,JEL,JL,JV))&
                                                           *ZC&
                                                           *ZLBDA**(ZCXI-ZDI)*MOMG(ZALPHAI,ZNUI,ZDI)
                      END IF
                      IF(LATT) THEN ! Calculation of Extinction coefficient
                        ! Rayleigh 3rd order
                        IF(NDIFF==0.OR.NDIFF==3) THEN
                          ZAETMP(:)=ZC*ZLBDA**ZCXI&
                                   *(ZDMELT_FACT*XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)&
                                   *MOMG(ZALPHAI,ZNUI,ZBI)/ZLBDA**ZBI)
                          ! Rayleigh 6th order
                        ELSE
                          ZAETMP(:)=ZC*ZLBDA**ZCXI*(&
                                    ZDMELT_FACT*XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)&
                                    *MOMG(ZALPHAI,ZNUI,ZBI)/ZLBDA**ZBI&
                                    +ZDMELT_FACT**(5./3.)*XPI**4/15./XLAM_RAD(JI)**3&
                                    *AIMAG(ZQK**2*(ZQMI**4+27.*ZQMI**2+38.)&
                                    /(2.*ZQMI**2+3.))*MOMG(ZALPHAI,ZNUI,5.*ZBI/3.)/ZLBDA**(5.*ZBI/3.) &
                                    +ZDMELT_FACT**2*2.*XPI**5/3. /XLAM_RAD(JI)**4*REAL(ZQK**2)&
                                    *MOMG(ZALPHAI,ZNUI,2.*ZBI)/ZLBDA**(2.*ZBI))
                        END IF
                      END IF
                    ELSE ! (if NDIFF=1 or NDIFF=7) => MIE (if choice=T-Matrix => Mie)
                      ZREFLOC(:)=0.
                      IF(LATT) ZAETMP(:)=0.
                      DO JJ=1,NPTS_GAULAG ! ****** Gauss-Laguerre quadrature
                        ZD=ZX(JJ)**(1./ZALPHAI)/ZLBDA !equivaut au ZDELTA_EQUIV olivier
                        ZRHOI=6*ZAI*ZD**(ZBI-3.)/XPI !pristine ice density
                        ZNUM=1.+2.*ZRHOI*(ZEPSI-1.)/(ZRHOPI*(ZEPSI+2.))
                        ZDEN=1.-ZRHOI*(ZEPSI-1.)/(ZRHOPI*(ZEPSI+2.))
                        ZQM=sqrt(ZNUM/ZDEN)
                        CALL BHMIE(XPI/XLAM_RAD(JI)*ZD,ZQM,ZQEXT(1),ZQSCA,ZQBACK(1))
                        ZQBACK(2)=ZQBACK(1)
                        ZQEXT(2)=ZQEXT(1) ! modif Clotilde 23/04/2012
                        ZQBACK(3)=0.
                        ZREFLOC(1:3)=ZREFLOC(1:3)+ZQBACK(1:3)*ZX(JJ)**(ZNUI-1.)*ZD**2*ZW(JJ)
                        ZREFLOC(4)=ZREFLOC(4)+ZQBACK(1)*ZX(JJ)**(ZNUI-1.+ZDI/ZALPHAI)*ZD**2*ZW(JJ)
                        IF(LATT) ZAETMP(:)=ZAETMP(:)+ZQEXT(:)*ZX(JJ)**(ZNUI-1.)*ZD**2*ZW(JJ)                  
                      END DO ! **************** end loop Gauss-Laguerre quadrature

                      ZREFLOC(1:2)=ZREFLOC(1:2)*1.E18*(XLAM_RAD(JI)/XPI)**4/.93*ZC &
                                   *ZLBDA**(ZCXI)/(4.*GAMMA(ZNUI))

                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                                           +PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1) &
                                           -ZCI*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(4) &
                                           *1.E18*(XLAM_RAD(JI)/XPI)**4*ZC &
                                           *ZLBDA**(ZCXI-ZDI)/(4.*GAMMA(ZNUI)*.93)
                      IF(LATT) ZAETMP(:)=ZAETMP(:)*XPI*ZC*ZLBDA**(ZCXI)/(4.*GAMMA(ZNUI))
                    END IF !**************** end loop for each type of diffusion
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)=ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)+ZREFLOC(1:3)
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEI)=ZREFLOC(1) ! z_e due to pristine ice
                    IF(LATT) THEN
                      ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)=ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)+ZAETMP(:)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAEI)=ZAETMP(1)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAVI)=ZAETMP(2)
                      IF(JL>1) ZAEIINT=ZAEIINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAEI)*XSTEP_RAD)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEI)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEI)*ZAEIINT ! Z_i attenuated
                    END IF            
                  END IF !********************* end IF (SIZE(PI_RAY,1)>0)
         
                  ! Total attenuation even if no hydrometeors
                  IF(LATT.AND.JL>1) ZREFL(JI,JEL,JAZ,JL,JH,JV,IATI)=ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IATI) &
                                                *EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAEI)*XSTEP_RAD)
                  ZRE_S22S11_I=0
                  ZIM_S22S11_I=0
                  ZS22_CARRE_I=0
                  ZS11_CARRE_I=0              
                END IF !******************** end IF (SIZE(PI_RAY,1)>0)
                !---------------------------------------------------------------------------------------------------
                !*       4.    SNOW 
                !              -----
                IF (SIZE(PS_RAY,1)>0)  THEN
                  ZM=PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PS_RAY(JI,JEL,JAZ,JL,JH,JV) !snow content
                  IF(ZM > ZM_MIN) THEN
                    YTYPE='s' 
                    !ZQMI: same formulation than for ice because snow is simulated only
                    !above melting leyer (3.5.4 p 67)
                    ZFW=0 
                    ZQMI=SQRT(QEPSI(PT_RAY(JI,JEL,JAZ,JL,JH,JV),XLIGHTSPEED/XLAM_RAD(JI)))
                    ZQK=(ZQMI**2-1.)/(ZQMI**2+2.) !ajout de Clotilde 23/04/2012
                    ZDMELT_FACT=6.*ZAS/(XPI*.92*XRHOLW)
                    ZEXP=2.*ZBS !XBS = 1.9 in ini_radar.f90 (bj tab 2.1 p24)
                    !dans ini_rain_ice.f90 :
                    IF ( (GLIMA .AND. LSNOW_T_L) .OR. (.NOT.GLIMA .AND. LSNOW_T_I) ) THEN
                       IF (PT_RAY(JI,JEL,JAZ,JL,JH,JV)>-10.) THEN
                          ZLBDA = MAX(MIN(XLBDAS_MAX, 10**(14.554-0.0423*(PT_RAY(JI,JEL,JAZ,JL,JH,JV)+273.15))),XLBDAS_MIN)
                       ELSE
                          ZLBDA = MAX(MIN(XLBDAS_MAX, 10**(6.226-0.0106*(PT_RAY(JI,JEL,JAZ,JL,JH,JV)+273.15))),XLBDAS_MIN)
                       END IF
                       ZN=ZNS*ZM*ZLBDA**ZBS
                    ELSE
                       ZLBDA= ZLBS*(ZM)**ZLBEXS
                       ZN=ZCCS*ZLBDA**ZCXS
                    END IF
                    ! Rayleigh or Rayleigh-Gans or Rayleigh with 6th order for attenuation
                    IF(NDIFF==0.OR.NDIFF==3.OR.NDIFF==4) THEN
                      ZREFLOC(1:2)=ZEQICE*.92**2*ZDMELT_FACT**2*1.E18*ZN*ZLBDA**(ZEXP)*MOMG(ZALPHAS,ZNUS,ZEXP)
                      ZREFLOC(3)=0.
                      IF(LWREFL) THEN ! weighting by reflectivities
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                               -ZCS*SIN(PELEV(JI,JEL,JL,JV))*ZEQICE*.92**2*ZDMELT_FACT**2&
                               *1.E18*ZN*ZLBDA**(ZEXP-ZDS)*MOMG(ZALPHAS,ZNUS,ZEXP+ZDS)
                      ELSE
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)+ZN
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                                                            -ZCS*SIN(PELEV(JI,JEL,JL,JV))&
                                           *ZN*ZLBDA**(ZDS)*MOMG(ZALPHAS,ZNUS,ZDS)
                      END IF
                      IF(LATT) THEN
                        IF(NDIFF==0.OR.NDIFF==3) THEN
                          ZAETMP(:)=ZN*(ZDMELT_FACT*XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)&
                                    *MOMG(ZALPHAS,ZNUS,ZBS)/ZLBDA**ZBS)
                        ELSE
                          ZAETMP(:)=ZN*(ZDMELT_FACT*XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)      &
                                    *MOMG(ZALPHAS,ZNUS,ZBS)/ZLBDA**ZBS                               &
                                    +ZDMELT_FACT**(5./3.)*XPI**4/15./XLAM_RAD(JI)**3                 &
                                    *AIMAG(ZQK**2*(ZQMI**4+27.*ZQMI**2+38.)                             &
                                    /(2.*ZQMI**2+3.))*MOMG(ZALPHAS,ZNUS,5.*ZBS/3.)/ZLBDA**(5.*ZBS/3.) &
                                    +ZDMELT_FACT**2   *2.*XPI**5/3. /XLAM_RAD(JI)**4*REAL(ZQK**2)     &
                                    *MOMG(ZALPHAS,ZNUS,2.*ZBS)/ZLBDA**(2.*ZBS))
                        END IF
                      END IF
                      ZRE_S22S11_S=0
                      ZIM_S22S11_S=0
                      ZS22_CARRE_S=0
                      ZS11_CARRE_S=0                        
                    !******************************* NDIFF==7 ************************************ 
                    ELSE IF(NDIFF==7) THEN

                      ZREFLOC(:)=0
                      IF(LATT) ZAETMP(:)=0      
                      CALL CALC_KTMAT(PELEV(JI,JEL,JL,JV), PT_RAY(JI,JEL,JAZ,JL,JH,JV),&
                                      ZFW,ZM,&
                                      ZELEV_MIN(2),ZELEV_MAX(2),ZELEV_STEP(2),&
                                      ZTC_MIN(2),ZTC_MAX(2),ZTC_STEP(2),&
                                      ZFW_MIN(2),ZFW_MAX(2),ZFW_STEP(2),&
                                      ZEXPM_MIN,ZEXPM_MAX,ZEXPM_STEP,&
                                      ITMAT,ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED)

                      IF (ITMAT(1) .NE. -NUNDEF) THEN 
                        DO JIND=1,SIZE(KMAT_COEF,2),1
                          KMAT_COEF(1,JIND)=ZS11_CARRE_T_S(ITMAT(JIND))
                          KMAT_COEF(2,JIND)=ZS22_CARRE_T_S(ITMAT(JIND))
                          KMAT_COEF(3,JIND)=ZRE_S22S11_T_S(ITMAT(JIND))
                          KMAT_COEF(4,JIND)=ZIM_S22S11_T_S(ITMAT(JIND))
                          KMAT_COEF(5,JIND)=ZRE_S22FMS11FT_T_S(ITMAT(JIND))
                          KMAT_COEF(6,JIND)=ZIM_S22FT_T_S(ITMAT(JIND))
                          KMAT_COEF(7,JIND)=ZIM_S11FT_T_S(ITMAT(JIND))
                        ENDDO 
                          CALL  INTERPOL(ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED,KMAT_COEF,ZS11_CARRE_S,ZS22_CARRE_S,&
                                       ZRE_S22S11_S,ZIM_S22S11_S,ZRE_S22FMS11F,ZIM_S22FT,ZIM_S11FT)   
                      ELSE
                        ZS11_CARRE_S=0
                        ZS22_CARRE_S=0
                        ZRE_S22S11_S=0
                        ZIM_S22S11_S=0
                        ZRE_S22FMS11F=0
                        ZIM_S22FT=0
                        ZIM_S11FT=0
                      END IF
                      ZREFLOC(1)=1.E18*(XLAM_RAD(JI))**4/(XPI**5*.93)*4*XPI*ZS22_CARRE_S
                      ZREFLOC(2)=1.E18*(XLAM_RAD(JI))**4/(XPI**5*.93)*4*XPI*ZS11_CARRE_S
                      ZREFLOC(3)=180.E3/XPI*XLAM_RAD(JI)*ZRE_S22FMS11F
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1) &
                                                        -ZCS*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(1) &
                                         *1.E18*(XLAM_RAD(JI)/XPI)**4/.93*(ZN*ZLBDA**(-ZCXS))/4./ZLBDA**(3+ZDS)
                      IF(LATT) THEN
                        ZAETMP(1)=ZIM_S22FT*XLAM_RAD(JI)*2
                        ZAETMP(2)=ZIM_S11FT*XLAM_RAD(JI)*2
                      END IF     
                    ELSE ! MIE 
                      ZREFLOC(:)=0.
                      IF(LATT) ZAETMP(:)=0.
                      DO JJ=1,NPTS_GAULAG ! ****** Gauss-Laguerre quadrature
                        ZD=ZX(JJ)**(1./ZALPHAS)/ZLBDA
                        ZDE=ZDMELT_FACT**(1./3.)*ZD**(ZBS/3.)
                        CALL BHMIE(XPI/XLAM_RAD(JI)*ZDE,ZQMI,ZQEXT(1),ZQSCA,ZQBACK(1))
                        ZQBACK(2)=ZQBACK(1)
                        ZQEXT(2)=ZQEXT(1) ! modif Clotilde 23/04/2012
                        ZQBACK(3)=0.
                        ZREFLOC(1:3)=ZREFLOC(1:3)+ZQBACK(1:3)*ZX(JJ)**(ZNUS-1.+2.*ZBS/3./ZALPHAS)*ZW(JJ)
                        ZREFLOC(4)=ZREFLOC(4)+ZQBACK(1)*ZX(JJ)**(ZNUS-1.+2.*ZBS/3./ZALPHAS+ZDS/ZALPHAS)*ZW(JJ)
                        IF(LATT) ZAETMP(:)=ZAETMP(:)+ZQEXT(:)*ZX(JJ)**(ZNUS-1.+2.*ZBS/3./ZALPHAS)*ZW(JJ)
                      END DO ! ****** end loop Gauss-Laguerre quadrature
                      ZREFLOC(1:2)=1.E18*(XLAM_RAD(JI)/XPI)**4*ZN*ZLBDA**(-2.*ZBS/3.)/&
                                  (4.*GAMMA(ZNUS)*.93)*ZDMELT_FACT**(2./3.)*ZREFLOC(1:2)               
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                                            +PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1) &
                                            -ZCS*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(4)   &
                                            *1.E18*(XLAM_RAD(JI)/XPI)**4*ZN          &
                                            *ZLBDA**(2.*ZBS/3.-ZDS)/              &
                                            (4.*GAMMA(ZNUS)*.93)*ZDMELT_FACT**(2./3.)                              
                      IF(LATT) ZAETMP(:)=ZAETMP(:)*XPI*ZN*ZLBDA**(-2.*ZBS/3.)/(4.*GAMMA(ZNUS))&
                                         *ZDMELT_FACT**(2./3.)            
                      ZRE_S22S11_S=0
                      ZIM_S22S11_S=0
                      ZS22_CARRE_S=0
                      ZS11_CARRE_S=0 
                    END IF !**************** end loop for each type of diffusion
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)=ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)+ZREFLOC(1:3)
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZES)=ZREFLOC(1) ! Z_e due to snow
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDS)=ZREFLOC(2) !Zvv for ZDR due to snow
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IKDS)=ZREFLOC(3) !Zvv for ZDR due to snow 
                    IF (ZS22_CARRE_S*ZS11_CARRE_S .GT. 0) THEN
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHS)=SQRT(ZRE_S22S11_S**2+ZIM_S22S11_S**2)/SQRT(ZS22_CARRE_S*ZS11_CARRE_S)
                    ELSE
                     ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHS)=1
                    END IF
                    IF(LATT) THEN
                     ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)=ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)+ZAETMP(:)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAES)=ZAETMP(1)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAVS)=ZAETMP(2) 
                      IF(JL>1) THEN
                        ZAESINT=ZAESINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAES)*XSTEP_RAD)
                        ZAVSINT=ZAVSINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAVS)*XSTEP_RAD)
                      ENDIF
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZES)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZES)*ZAESINT ! Z_s attenuated
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDS)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDS)*ZAVSINT ! ZVs attenuated
                    END IF !end IF(LATT)
                  END IF !end IF(PS_RAY(JI,JEL,JAZ,JL,JH,JV) > ...)


                 ! Total attenuation even if no hydrometeors
                  IF(LATT.AND.JL>1) ZREFL(JI,JEL,JAZ,JL,JH,JV,IATS)=ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IATS) &
                                                *EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAES)*XSTEP_RAD)
                END IF !END IF (SIZE(PS_RAY,1)>0)
                !---------------------------------------------------------------------------------------------------
                !*       5.    GRAUPEL
                !              -------
                !
                !ZDG=.5 ! from Bringi & Chandrasekar 2001, p. 433
                IF (SIZE(PG_RAY,1)>0)  THEN
                  ZM=PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PG_RAY(JI,JEL,JAZ,JL,JH,JV) !graupel content
                  IF(ZM > ZM_MIN) THEN
                    YTYPE='g'
                    ZQMI=SQRT(QEPSI(MIN(PT_RAY(JI,JEL,JAZ,JL,JH,JV),XTT),XLIGHTSPEED/XLAM_RAD(JI)))
                    ZQMW=SQRT(QEPSW(MAX(PT_RAY(JI,JEL,JAZ,JL,JH,JV),XTT),XLIGHTSPEED/XLAM_RAD(JI)))
                    !ini_radar.f90 : ZCXG = -0.5 XBG = 2.8 ( Xj et bj tab 2.1 p 24)
                    !ini_rain_ice.f90 : XLBEXG = 1.0/(XCXG-XBG) XAG = 19.6 (aj tab 2.1 p 24)
                    !XLBG   = ( XAG*XCCG*MOMG(XALPHAG,XNUG,XBG) )**(-XLBEXG) (eq 2.6 p 23)
                    IF (PR_RAY(JI,JEL,JAZ,JL,JH,JV) > ZRTMIN(3) ) THEN
                      ZFW=PR_RAY(JI,JEL,JAZ,JL,JH,JV)/(PR_RAY(JI,JEL,JAZ,JL,JH,JV)+PG_RAY(JI,JEL,JAZ,JL,JH,JV))
                    ELSE
                      ZFW=0.
                    END IF
                    ZLBDA=ZLBG*(PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PG_RAY(JI,JEL,JAZ,JL,JH,JV))**ZLBEXG
                    !XTT : température du point triple de l'eau (273.16 K <=> 0.1 °C)
                    IF(PT_RAY(JI,JEL,JAZ,JL,JH,JV) > XTT) THEN ! mixture of ice and water
                      ZFRAC_ICE = .85 !(see p 68)
                    ELSE ! only ice
                      ZFRAC_ICE=1.
                    END IF
                    ! from eq 3.77 p 68
                    !XRHOLW=1000 (initialized in ini_cst.f90)
                    ZDMELT_FACT=6.*ZAG/(XPI*XRHOLW*((1.-ZFRAC_ICE)+ZFRAC_ICE*0.92))
                    ZEXP=2.*ZBG
                    !Calculation of the refractive index from Bohren and Battan (3.72 p66)
                    ZQB=2.*ZQMW**2*(2.*ZQMI**2*LOG(ZQMI/ZQMW)/(ZQMI**2-ZQMW**2)-1.)/(ZQMI**2-ZQMW**2) !Beta (3.73 p66)
                    ZQM=SQRT(((1.-ZFRAC_ICE)*ZQMW**2+ZFRAC_ICE*ZQB*ZQMI**2)/(1.-ZFRAC_ICE+ZFRAC_ICE*ZQB)) ! Bohren & Battan (1982) 3.72 p66
                    ZQK=(ZQM**2-1.)/(ZQM**2+2.)                                       
                    !Rayleigh, Rayleigh for ellipsoides or Rayleigh 6th order
                    IF(NDIFF==0.OR.NDIFF==3.OR.NDIFF==4) THEN
                      ZREFLOC(1:2)=ABS(ZQK)**2/.93*ZDMELT_FACT**2*1.E18*ZCCG*ZLBDA**(ZCXG-ZEXP)*MOMG(ZALPHAG,ZNUG,ZEXP)
                      ZREFLOC(3)=0.
                      IF(LWREFL) THEN ! weighting by reflectivities
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                              -ZCG*SIN(PELEV(JI,JEL,JL,JV))*ABS(ZQK)**2/.93*ZDMELT_FACT**2&
                              *1.E18*ZCCG*ZLBDA**(ZCXG-ZEXP-ZDG)*MOMG(ZALPHAG,ZNUG,ZEXP+ZDG)
                      ELSE
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)+ZCCG*ZLBDA**ZCXG
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                                                        -ZCG*SIN(PELEV(JI,JEL,JL,JV))&
                                                        *ZCCG*ZLBDA**(ZCXG-ZDG)*MOMG(ZALPHAG,ZNUG,ZDG)
                      END IF !end IF(LWREFL)
                      IF(LATT) THEN
                        IF(NDIFF==0.OR.NDIFF==3) THEN
                          ZAETMP(:)=ZCCG*ZLBDA**ZCXG*(ZDMELT_FACT*XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)   &
                                    *MOMG(ZALPHAG,ZNUG,ZBG)/ZLBDA**ZBG)
                        ELSE
                          ZAETMP(:)=ZCCG*ZLBDA**ZCXG*(ZDMELT_FACT*XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)   &
                                   *MOMG(ZALPHAG,ZNUG,ZBG)/ZLBDA**ZBG&
                                   +ZDMELT_FACT**(5./3.)*XPI**4/15./XLAM_RAD(JI)**3               &
                                   *AIMAG(ZQK**2*(ZQM**4+27.*ZQM**2+38.)                             &
                                   /(2.*ZQM**2+3.))*MOMG(ZALPHAG,ZNUG,5.*ZBG/3.)/ZLBDA**(5.*ZBG/3.)&
                                   +ZDMELT_FACT**2   *2.*XPI**5/3. /XLAM_RAD(JI)**4*REAL(ZQK**2)   &
                                   *MOMG(ZALPHAG,ZNUG,2.*ZBG)   /ZLBDA**(2.*ZBG))
                        END IF ! end IF(NDIFF==0.OR.NDIFF==3)
                      END IF ! end IF(LATT)
                      ZRE_S22S11_G=0
                      ZIM_S22S11_G=0
                      ZS22_CARRE_G=0
                      ZS11_CARRE_G=0                    
                      !******************************* NDIFF==7 TmatInt ************************************ 
                    ELSE IF(NDIFF==7) THEN
                      ZREFLOC(:)=0
                      IF(LATT) ZAETMP(:)=0                
                      IF (ZFW < 0.01) THEN !******** DRY GRAUPEL 
                         CALL CALC_KTMAT(PELEV(JI,JEL,JL,JV), PT_RAY(JI,JEL,JAZ,JL,JH,JV),&
                                        ZFW,ZM,&
                                        ZELEV_MIN(3),ZELEV_MAX(3),ZELEV_STEP(3),&
                                        ZTC_MIN(3),ZTC_MAX(3),ZTC_STEP(3),&
                                        ZFW_MIN(3),ZFW_MAX(3),ZFW_STEP(3),&
                                        ZEXPM_MIN,ZEXPM_MAX,ZEXPM_STEP,&
                                        ITMAT,ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED)
                        IF (ITMAT(1) .NE. -NUNDEF) THEN 
                          DO JIND=1,SIZE(KMAT_COEF,2),1
                            KMAT_COEF(1,JIND)=ZS11_CARRE_T_G(ITMAT(JIND))
                            KMAT_COEF(2,JIND)=ZS22_CARRE_T_G(ITMAT(JIND))
                            KMAT_COEF(3,JIND)=ZRE_S22S11_T_G(ITMAT(JIND))
                            KMAT_COEF(4,JIND)=ZIM_S22S11_T_G(ITMAT(JIND))
                            KMAT_COEF(5,JIND)=ZRE_S22FMS11FT_T_G(ITMAT(JIND))
                            KMAT_COEF(6,JIND)=ZIM_S22FT_T_G(ITMAT(JIND))
                            KMAT_COEF(7,JIND)=ZIM_S11FT_T_G(ITMAT(JIND))
                          ENDDO
                         CALL  INTERPOL(ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED,KMAT_COEF,ZS11_CARRE_G,ZS22_CARRE_G,&
                                         ZRE_S22S11_G,ZIM_S22S11_G,ZRE_S22FMS11F,ZIM_S22FT,ZIM_S11FT)
  ELSE
                          ZS11_CARRE_G=0
                          ZS22_CARRE_G=0
                          ZRE_S22S11_G=0
                          ZIM_S22S11_G=0
                          ZRE_S22FMS11F=0
                          ZIM_S22FT=0
                          ZIM_S11FT=0
                        END IF                    
                      ELSE !ZFW >= 0.01 ************** WET GRAUPEL
                        CALL CALC_KTMAT(PELEV(JI,JEL,JL,JV),PT_RAY(JI,JEL,JAZ,JL,JH,JV),&
                                        ZFW,ZM,&
                                        ZELEV_MIN(4),ZELEV_MAX(4),ZELEV_STEP(4),&
                                        ZTC_MIN(4),ZTC_MAX(4),ZTC_STEP(4),&
                                        ZFW_MIN(4),ZFW_MAX(4),ZFW_STEP(4),&
                                        ZEXPM_MIN,ZEXPM_MAX,ZEXPM_STEP,&
                                        ITMAT,ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED)
                        IF (ITMAT(1) .NE. -NUNDEF) THEN  
                          DO JIND=1,SIZE(KMAT_COEF,2),1
                            KMAT_COEF(1,JIND)=ZS11_CARRE_T_W(ITMAT(JIND))
                            KMAT_COEF(2,JIND)=ZS22_CARRE_T_W(ITMAT(JIND))
                            KMAT_COEF(3,JIND)=ZRE_S22S11_T_W(ITMAT(JIND))
                            KMAT_COEF(4,JIND)=ZIM_S22S11_T_W(ITMAT(JIND))
                            KMAT_COEF(5,JIND)=ZRE_S22FMS11FT_T_W(ITMAT(JIND))
                            KMAT_COEF(6,JIND)=ZIM_S22FT_T_W(ITMAT(JIND))
                            KMAT_COEF(7,JIND)=ZIM_S11FT_T_W(ITMAT(JIND))
                          ENDDO
                          CALL  INTERPOL(ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED,KMAT_COEF,ZS11_CARRE_G,ZS22_CARRE_G,&
                                         ZRE_S22S11_G,ZIM_S22S11_G,ZRE_S22FMS11F,ZIM_S22FT,ZIM_S11FT) 
                        ELSE
                          ZS11_CARRE_G=0
                          ZS22_CARRE_G=0
                          ZRE_S22S11_G=0
                          ZIM_S22S11_G=0
                          ZRE_S22FMS11F=0
                          ZIM_S22FT=0
                          ZIM_S11FT=0
                        END IF                    
                      END IF!END IF (ZFW<0.01)
                      ZREFLOC(1)=1.E18*(XLAM_RAD(JI))**4/(XPI**5*.93)*4*XPI*ZS22_CARRE_G
                      ZREFLOC(2)=1.E18*(XLAM_RAD(JI))**4/(XPI**5*.93)*4*XPI*ZS11_CARRE_G
                      ZREFLOC(3)=180.E3/XPI*XLAM_RAD(JI)*ZRE_S22FMS11F
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1) &
                                                        -ZCG*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(1) &
                                         *1.E18*(XLAM_RAD(JI)/XPI)**4/.93*ZCCG/4./ZLBDA**(3+ZDG)
                      IF(LATT) THEN
                        ZAETMP(1)=ZIM_S22FT*XLAM_RAD(JI)*2
                        ZAETMP(2)=ZIM_S11FT*XLAM_RAD(JI)*2
                      END IF            
                    ELSE ! Mie (NDIFF=1)
                      ZREFLOC(:)=0.
                      IF(LATT) ZAETMP(:)=0.
                      DO JJ=1,NPTS_GAULAG ! ****** Gauss-Laguerre quadrature
                        ZD=ZX(JJ)**(1./ZALPHAG)/ZLBDA
                        ZDE=ZDMELT_FACT**(1./3.)*ZD**(ZBG/3.)
                        CALL BHMIE(XPI/XLAM_RAD(JI)*ZDE,ZQM,ZQEXT(1),ZQSCA,ZQBACK(1))
                        ZQBACK(2)=ZQBACK(1)
                        ZQEXT(2)=ZQEXT(1) ! modif Clotilde 23/04/2012
                        ZQBACK(3)=0.
                        ZREFLOC(1:3)=ZREFLOC(1:3)+ZQBACK(1:3)*ZX(JJ)**(ZNUG-1.+2.*ZBG/3./ZALPHAG)*ZW(JJ)
                        ZREFLOC(4)=ZREFLOC(4)+ZQBACK(1)*ZX(JJ)**(ZNUG-1.+2.*ZBG/3./ZALPHAG+ZDG/ZALPHAG)*ZW(JJ)
                        IF(LATT) ZAETMP(:)=ZAETMP(:)+ZQEXT(:)*ZX(JJ)**(ZNUG-1.+2.*ZBG/3./ZALPHAG)*ZW(JJ)
                      END DO ! ****** end loop on diameter (Gauss-Laguerre)                     
                      ZREFLOC(1:2)=ZREFLOC(1:2)*1.E18*(XLAM_RAD(JI)/XPI)**4*ZCCG                      &
                                    *ZLBDA**(ZCXG-2.*ZBG/3.)/(4.*GAMMA(ZNUG)*.93)*ZDMELT_FACT**(2./3.)                    
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)               &
                                    +PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1)                        &
                                    -ZCG*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(4)                          &
                                    *1.E18*(XLAM_RAD(JI)/XPI)**4*ZCCG                                 &
                                    *ZLBDA**(ZCXG-2.*ZBG/3.-ZDG)/(4.*GAMMA(ZNUG)*.93)*ZDMELT_FACT**(2./3.)               
                      IF(LATT) ZAETMP(:)=ZAETMP(:)*XPI*ZCCG*ZLBDA**(ZCXG-2.*ZBG/3.)/(4.*GAMMA(ZNUG))  &
                                    *ZDMELT_FACT**(2./3.)
                      ZRE_S22S11_G=0
                      ZIM_S22S11_G=0
                      ZS22_CARRE_G=0
                      ZS11_CARRE_G=0  !0 in case of Mie
                    END IF !**************** end loop for each type of diffusion : IF(NDIFF==0.OR.NDIFF==3.OR.NDIFF==4)          
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)=ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)+ZREFLOC(1:3)
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEG)=ZREFLOC(1) ! z_e due to graupel
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDG)=ZREFLOC(2) !Zvv for ZDR due to graupel
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IKDG)=ZREFLOC(3) !Zvv for ZDR due to graupel  

                    IF (ZS22_CARRE_G*ZS11_CARRE_G .GT. 0) THEN
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHG)=SQRT(ZRE_S22S11_G**2+ZIM_S22S11_G**2)/SQRT(ZS22_CARRE_G*ZS11_CARRE_G)
                    ELSE
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHG)=1
                    END IF
                    IF(LATT) THEN
                      ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)=ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)+ZAETMP(:)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAEG)=ZAETMP(1)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAVG)=ZAETMP(2)
                      IF(JL>1) THEN
                        ZAEGINT=ZAEGINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAEG)*XSTEP_RAD)
                        ZAVGINT=ZAVGINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAVG)*XSTEP_RAD)
                      END IF
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEG)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEG)*ZAEGINT ! Z_g attenuated
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDG)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDG)*ZAVGINT ! Z_g attenuated
                    END IF !end IF(LATT)            
                  END IF !**************** IF(PG_RAY(JI,JEL,JAZ,JL,JH,JV) > XRTMIN(6))
         
                  ! Total attenuation even if no hydrometeors
                  IF(LATT.AND.JL>1) ZREFL(JI,JEL,JAZ,JL,JH,JV,IATG)=ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IATG) &
                                                *EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAEG)*XSTEP_RAD)
              
                END IF ! **************** end GRAUPEL (end IF SIZE(PG_RAY,1) > 0)
                !-----------------------------------------------------------------------------------------------
                !-----------------------------------------------------------------------------------------------
!**********************************
!**********************************
!**********************************
!**********************************


!---------------------------------------------------------------------------------------------------
                !*       6. HAIL
                !              -------
                !
                !
                IF (GHAIL)  THEN
                  ZM=PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PH_RAY(JI,JEL,JAZ,JL,JH,JV) !graupel content
                  IF(ZM > ZM_MIN) THEN
                    YTYPE='h'
                    ZQMI=SQRT(QEPSI(MIN(PT_RAY(JI,JEL,JAZ,JL,JH,JV),XTT),XLIGHTSPEED/XLAM_RAD(JI)))
                    ZQMW=SQRT(QEPSW(MAX(PT_RAY(JI,JEL,JAZ,JL,JH,JV),XTT),XLIGHTSPEED/XLAM_RAD(JI)))
                    !ini_radar.f90 : ZCXG = -0.5 XBG = 2.8 ( Xj et bj tab 2.1 p 24)
                    !ini_rain_ice.f90 : XLBEXG = 1.0/(XCXG-XBG) XAG = 19.6 (aj tab 2.1 p 24)
                    !XLBG   = ( XAG*XCCG*MOMG(XALPHAG,XNUG,XBG) )**(-XLBEXG) (eq 2.6 p 23)
ZFW=0    !????????
                    ZLBDA=ZLBH*(PRHODREF_RAY(JI,JEL,JAZ,JL,JH,JV)*PH_RAY(JI,JEL,JAZ,JL,JH,JV))**ZLBEXH
                    !XTT : température du point triple de l'eau (273.16 K <=> 0.1 °C)
                    IF(PT_RAY(JI,JEL,JAZ,JL,JH,JV) > XTT) THEN ! mixture of ice and water
                      ZFRAC_ICE = .85 !(see p 68)
                    ELSE ! only ice
                      ZFRAC_ICE=1.
                    END IF
                    ! from eq 3.77 p 68
                    !XRHOLW=1000 (initialized in ini_cst.f90)
                    ZDMELT_FACT=6.*ZAG/(XPI*XRHOLW*((1.-ZFRAC_ICE)+ZFRAC_ICE*0.92))
                    ZEXP=2.*ZBH
                    !Calculation of the refractive index from Bohren and Battan (3.72 p66)
                    ZQB=2.*ZQMW**2*(2.*ZQMI**2*LOG(ZQMI/ZQMW)/(ZQMI**2-ZQMW**2)-1.)/(ZQMI**2-ZQMW**2) !Beta (3.73 p66)
                    ZQM=SQRT(((1.-ZFRAC_ICE)*ZQMW**2+ZFRAC_ICE*ZQB*ZQMI**2)/(1.-ZFRAC_ICE+ZFRAC_ICE*ZQB)) ! Bohren & Battan (1982) 3.72 p66
                    ZQK=(ZQM**2-1.)/(ZQM**2+2.)                                       
                    !Rayleigh, Rayleigh for ellipsoides or Rayleigh 6th order
                    IF(NDIFF==0.OR.NDIFF==3.OR.NDIFF==4) THEN
                      ZREFLOC(1:2)=ABS(ZQK)**2/.93*ZDMELT_FACT**2*1.E18*ZCCH*ZLBDA**(ZCXH-ZEXP)*MOMG(ZALPHAH,ZNUH,ZEXP)
                      ZREFLOC(3)=0.
                      IF(LWREFL) THEN ! weighting by reflectivities
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                              -ZCH*SIN(PELEV(JI,JEL,JL,JV))*ABS(ZQK)**2/.93*ZDMELT_FACT**2&
                              *1.E18*ZCCH*ZLBDA**(ZCXH-ZEXP-ZDH)*MOMG(ZALPHAH,ZNUH,ZEXP+ZDH)
                      ELSE
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)+ZCCH*ZLBDA**ZCXH                       
                        ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                                                        -ZCH*SIN(PELEV(JI,JEL,JL,JV))&
                                                        *ZCCH*ZLBDA**(ZCXH-ZDH)*MOMG(ZALPHAH,ZNUH,ZDH)
                      END IF !end IF(LWREFL)
                      IF(LATT) THEN
                        IF(NDIFF==0.OR.NDIFF==3) THEN
                          ZAETMP(:)=ZCCH*ZLBDA**ZCXH*(ZDMELT_FACT*XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)   &
                                    *MOMG(ZALPHAH,ZNUH,ZBH)/ZLBDA**ZBH)
                        ELSE
                          ZAETMP(:)=ZCCH*ZLBDA**ZCXH*(ZDMELT_FACT*XPI**2/XLAM_RAD(JI)*AIMAG(ZQK)   &
                                   *MOMG(ZALPHAH,ZNUH,ZBH)/ZLBDA**ZBH&
                                   +ZDMELT_FACT**(5./3.)*XPI**4/15./XLAM_RAD(JI)**3               &
                                   *AIMAG(ZQK**2*(ZQM**4+27.*ZQM**2+38.)                             &
                                   /(2.*ZQM**2+3.))*MOMG(ZALPHAH,ZNUH,5.*ZBH/3.)/ZLBDA**(5.*ZBH/3.)&
                                   +ZDMELT_FACT**2   *2.*XPI**5/3. /XLAM_RAD(JI)**4*REAL(ZQK**2)   &
                                   *MOMG(ZALPHAH,ZNUH,2.*ZBH)   /ZLBDA**(2.*ZBH))
                        END IF ! end IF(NDIFF==0.OR.NDIFF==3)
                      END IF ! end IF(LATT)
                      ZRE_S22S11_H=0
                      ZIM_S22S11_H=0
                      ZS22_CARRE_H=0
                      ZS11_CARRE_H=0                    
                      !******************************* NDIFF==7 TmatInt ************************************ 
                    ELSE IF(NDIFF==7) THEN
                      ZREFLOC(:)=0
                      IF(LATT) ZAETMP(:)=0                
                         CALL CALC_KTMAT(PELEV(JI,JEL,JL,JV), PT_RAY(JI,JEL,JAZ,JL,JH,JV),&
                                        ZFW,ZM,&
                                        ZELEV_MIN(3),ZELEV_MAX(3),ZELEV_STEP(3),&
                                        ZTC_MIN(3),ZTC_MAX(3),ZTC_STEP(3),&
                                        ZFW_MIN(3),ZFW_MAX(3),ZFW_STEP(3),&
                                        ZEXPM_MIN,ZEXPM_MAX,ZEXPM_STEP,&
                                        ITMAT,ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED)
                        IF (ITMAT(1) .NE. -NUNDEF) THEN 
                          DO JIND=1,SIZE(KMAT_COEF,2),1
                            KMAT_COEF(1,JIND)=ZS11_CARRE_T_H(ITMAT(JIND))
                            KMAT_COEF(2,JIND)=ZS22_CARRE_T_H(ITMAT(JIND))
                            KMAT_COEF(3,JIND)=ZRE_S22S11_T_H(ITMAT(JIND))
                            KMAT_COEF(4,JIND)=ZIM_S22S11_T_H(ITMAT(JIND))
                            KMAT_COEF(5,JIND)=ZRE_S22FMS11FT_T_H(ITMAT(JIND))
                            KMAT_COEF(6,JIND)=ZIM_S22FT_T_H(ITMAT(JIND))
                            KMAT_COEF(7,JIND)=ZIM_S11FT_T_H(ITMAT(JIND))
                          ENDDO
                         CALL  INTERPOL(ZELEV_RED,ZTC_RED,ZFW_RED,ZM_RED,KMAT_COEF,ZS11_CARRE_H,ZS22_CARRE_H,&
                                         ZRE_S22S11_H,ZIM_S22S11_H,ZRE_S22FMS11F,ZIM_S22FT,ZIM_S11FT)
                        ELSE
                          ZS11_CARRE_H=0
                          ZS22_CARRE_H=0
                          ZRE_S22S11_H=0
                          ZIM_S22S11_H=0
                          ZRE_S22FMS11F=0
                          ZIM_S22FT=0
                          ZIM_S11FT=0
                        END IF                    
                      ZREFLOC(1)=1.E18*(XLAM_RAD(JI))**4/(XPI**5*.93)*4*XPI*ZS22_CARRE_H
                      ZREFLOC(2)=1.E18*(XLAM_RAD(JI))**4/(XPI**5*.93)*4*XPI*ZS11_CARRE_H
                      ZREFLOC(3)=180.E3/XPI*XLAM_RAD(JI)*ZRE_S22FMS11F
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1) &
                                                        -ZCH*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(1) &
                                         *1.E18*(XLAM_RAD(JI)/XPI)**4/.93*ZCCH/4./ZLBDA**(3+ZDH)
                      IF(LATT) THEN
                        ZAETMP(1)=ZIM_S22FT*XLAM_RAD(JI)*2
                        ZAETMP(2)=ZIM_S11FT*XLAM_RAD(JI)*2
                      END IF            
                    ELSE ! Mie (NDIFF=1)
                      ZREFLOC(:)=0.
                      IF(LATT) ZAETMP(:)=0.
                      DO JJ=1,NPTS_GAULAG ! ****** Gauss-Laguerre quadrature
                        ZD=ZX(JJ)**(1./ZALPHAH)/ZLBDA
                        ZDE=ZDMELT_FACT**(1./3.)*ZD**(ZBH/3.)
                        CALL BHMIE(XPI/XLAM_RAD(JI)*ZDE,ZQM,ZQEXT(1),ZQSCA,ZQBACK(1))
                        ZQBACK(2)=ZQBACK(1)
                        ZQEXT(2)=ZQEXT(1) ! modif Clotilde 23/04/2012
                        ZQBACK(3)=0.
                        ZREFLOC(1:3)=ZREFLOC(1:3)+ZQBACK(1:3)*ZX(JJ)**(ZNUH-1.+2.*ZBH/3./ZALPHAH)*ZW(JJ)
                        ZREFLOC(4)=ZREFLOC(4)+ZQBACK(1)*ZX(JJ)**(ZNUH-1.+2.*ZBH/3./ZALPHAH+ZDH/ZALPHAH)*ZW(JJ)
                        IF(LATT) ZAETMP(:)=ZAETMP(:)+ZQEXT(:)*ZX(JJ)**(ZNUH-1.+2.*ZBH/3./ZALPHAH)*ZW(JJ)
                      END DO ! ****** end loop on diameter (Gauss-Laguerre)                     
                      ZREFLOC(1:2)=ZREFLOC(1:2)*1.E18*(XLAM_RAD(JI)/XPI)**4*ZCCH                      &
                                    *ZLBDA**(ZCXH-2.*ZBH/3.)/(4.*GAMMA(ZNUH)*.93)*ZDMELT_FACT**(2./3.)                    
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)               &
                                    +PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFLOC(1)                        &
                                    -ZCH*SIN(PELEV(JI,JEL,JL,JV))*ZREFLOC(4)                          &
                                    *1.E18*(XLAM_RAD(JI)/XPI)**4*ZCCH                                 &
                                    *ZLBDA**(ZCXH-2.*ZBH/3.-ZDH)/(4.*GAMMA(ZNUH)*.93)*ZDMELT_FACT**(2./3.)               
                      IF(LATT) ZAETMP(:)=ZAETMP(:)*XPI*ZCCH*ZLBDA**(ZCXH-2.*ZBH/3.)/(4.*GAMMA(ZNUH))  &
                                    *ZDMELT_FACT**(2./3.)
                      ZRE_S22S11_H=0
                      ZIM_S22S11_H=0
                      ZS22_CARRE_H=0
                      ZS11_CARRE_H=0  !0 in case of Mie
                    END IF !**************** end loop for each type of diffusion : IF(NDIFF==0.OR.NDIFF==3.OR.NDIFF==4)          
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)=ZREFL(JI,JEL,JAZ,JL,JH,JV,1:3)+ZREFLOC(1:3)
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEH)=ZREFLOC(1) ! z_e due to graupel
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDH)=ZREFLOC(2) !Zvv for ZDR due to graupel
                    ZREFL(JI,JEL,JAZ,JL,JH,JV,IKDH)=ZREFLOC(3) !Zvv for ZDR due to graupel  

                    IF (ZS22_CARRE_H*ZS11_CARRE_H .GT. 0) THEN
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHH)=SQRT(ZRE_S22S11_H**2+ZIM_S22S11_H**2)/SQRT(ZS22_CARRE_H*ZS11_CARRE_H)
                    ELSE
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHH)=1
                    END IF
                    IF(LATT) THEN
                      ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)=ZAELOC(JI,JEL,JAZ,JL,JH,JV,:)+ZAETMP(:)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAEH)=ZAETMP(1)
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IAVH)=ZAETMP(2)
                      IF(JL>1) THEN
                        ZAEHINT=ZAEHINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAEH)*XSTEP_RAD)
                        ZAVHINT=ZAVHINT*EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAVH)*XSTEP_RAD)
                      END IF
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEH)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZEH)*ZAEHINT ! Z_g attenuated
                      ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDH)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IZDH)*ZAVHINT ! Z_g attenuated
                    END IF !end IF(LATT)            
                  END IF !**************** IF(PH_RAY(JI,JEL,JAZ,JL,JH,JV) > XRTMIN(6))
         
                  ! Total attenuation even if no hydrometeors
                  IF(LATT.AND.JL>1) ZREFL(JI,JEL,JAZ,JL,JH,JV,IATH)=ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IATH) &
                                                *EXP(-2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IAEH)*XSTEP_RAD)
              
                END IF ! **************** end HAIL (end IF SIZE(PH_RAY,1) > 0)
                !-----------------------------------------------------------------------------------------------
                !-----------------------------------------------------------------------------------------------
!**********************************
!**********************************
!**********************************
!**********************************

                IF(LWREFL) THEN ! weighting by reflectivities
                  ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                                   +PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFL(JI,JEL,JAZ,JL,JH,JV,1)             
                ELSE IF(LWBSCS) THEN ! weighting by hydrometeor concentrations
                  ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)&
                                   +PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)*ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)
                ELSE IF(ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)/=0.) THEN ! no weighting
                  ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)/ZREFL(JI,JEL,JAZ,JL,JH,JV,IMAX)&
                                   +PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)
                END IF
                !Calculation of  Phidp  (ZREFL(JI,JEL,JAZ,JL,JH,JV,IPDP) is initialized to 0 before the loop
                IF (JL>1) ZREFL(JI,JEL,JAZ,JL,JH,JV,IPDP)=ZREFL(JI,JEL,JAZ,JL-1,JH,JV,IPDP)+ &
                                   2.*ZREFL(JI,JEL,JAZ,JL-1,JH,JV,3)*XSTEP_RAD*1D-3

                !Calculation of RhoHV and DeltaHV
                ZRE_S22S11_T=ZRE_S22S11_R+ZRE_S22S11_I+ZRE_S22S11_S+ZRE_S22S11_G+ZRE_S22S11_H
                ZIM_S22S11_T=ZIM_S22S11_R+ZIM_S22S11_I+ZIM_S22S11_S+ZIM_S22S11_G+ZIM_S22S11_H
                ZS22_CARRE_T=ZS22_CARRE_R+ZS22_CARRE_I+ZS22_CARRE_S+ZS22_CARRE_G+ZS22_CARRE_H
                ZS11_CARRE_T=ZS11_CARRE_R+ZS11_CARRE_I+ZS11_CARRE_S+ZS11_CARRE_G+ZS11_CARRE_H
                !RhoHV 
                IF ((ZS22_CARRE_T*ZS11_CARRE_T)>0.) THEN
                  ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHV)=SQRT(ZRE_S22S11_T**2+ZIM_S22S11_T**2)/SQRT(ZS22_CARRE_T*ZS11_CARRE_T)   
                ELSE
                  ZREFL(JI,JEL,JAZ,JL,JH,JV,IRHV)=-XUNDEF
                END IF 
                !DeltaHV
                IF (ZRE_S22S11_T/=0) THEN
                  ZREFL(JI,JEL,JAZ,JL,JH,JV,IDHV)=180/XPI*ATAN(ZIM_S22S11_T/ZRE_S22S11_T)
                ELSE
                  ZREFL(JI,JEL,JAZ,JL,JH,JV,IDHV)=0
                END IF
              ELSE !if temperature is not defined
                ZREFL(JI,JEL,JAZ,JL,JH,JV,1:2)=XVALGROUND
                ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=XVALGROUND
                LPART_MASK=.TRUE.
              END IF !end condition : IF(PT_RAY(JI,JEL,JAZ,JL,JH,JV) /= -XUNDEF) => if temperature is defined
            END IF !end condition : IF(LPART_MASK) => if pixel is not masked
          END DO LOOPJL
        END DO !JV
      END DO !JH
    END DO !JAZ
  END DO !JEL   
  !
  !lookup tables for rain
  DEALLOCATE (ZTC_T_R,ZELEV_T_R,ZM_T_R,ZS11_CARRE_T_R,ZS22_CARRE_T_R,&
  ZRE_S22S11_T_R,ZIM_S22S11_T_R,ZRE_S22FMS11FT_T_R,ZIM_S22FT_T_R,ZIM_S11FT_T_R)
  !lookup tables for snow
  DEALLOCATE (ZTC_T_S,ZELEV_T_S,ZM_T_S,ZS11_CARRE_T_S,ZS22_CARRE_T_S,&
  ZRE_S22S11_T_S,ZIM_S22S11_T_S,ZRE_S22FMS11FT_T_S,ZIM_S22FT_T_S,ZIM_S11FT_T_S)
  !lookup tables for graupel
  DEALLOCATE (ZTC_T_G,ZELEV_T_G,ZM_T_G,ZS11_CARRE_T_G,ZS22_CARRE_T_G,&
  ZRE_S22S11_T_G,ZIM_S22S11_T_G,ZRE_S22FMS11FT_T_G,ZIM_S22FT_T_G,ZIM_S11FT_T_G)
  !lookup tables for wet graupel
  DEALLOCATE (ZTC_T_W,ZELEV_T_W,ZM_T_W,ZS11_CARRE_T_W,ZS22_CARRE_T_W,&
  ZRE_S22S11_T_W,ZIM_S22S11_T_W,ZRE_S22FMS11FT_T_W,ZIM_S22FT_T_W,ZIM_S11FT_T_W)   
  IF (GHAIL) THEN
    !lookup tables for hail
    DEALLOCATE (ZTC_T_H,ZELEV_T_H,ZM_T_H,ZS11_CARRE_T_H,ZS22_CARRE_T_H,&
    ZRE_S22S11_T_H,ZIM_S22S11_T_H,ZRE_S22FMS11FT_T_H,ZIM_S22FT_T_H,ZIM_S11FT_T_H)  
  ENDIF
END DO !JI
!
! attenuation in dB/km
IF(LATT) ZREFL(:,:,:,:,:,:,IAER:IAEH)=4343.*2.*ZREFL(:,:,:,:,:,:,IAER:IAEH) ! horizontal specific attenuation
IF(LATT) ZREFL(:,:,:,:,:,:,IAVR:IAVH)=4343.*2.*ZREFL(:,:,:,:,:,:,IAVR:IAVH) ! vertical specific attenuation
! convective/stratiform
ZREFL(:,:,:,:,:,:,4)=PBU_MASK_RAY(:,:,:,:,:,:) ! CSR
! /convective/stratiform

WRITE(ILUOUT0,*) 'NB ZREFL  MIN MAX :', MINVAL(ZREFL(:,:,:,:,:,:,:)),MAXVAL(ZREFL(:,:,:,:,:,:,:))
WRITE(ILUOUT0,*) 'NB ZREFL VALGROUND :', COUNT(ZREFL(:,:,:,:,:,:,:) ==XVALGROUND)
WRITE(ILUOUT0,*) 'NB ZREFL -XUNDEF :', COUNT(ZREFL(:,:,:,:,:,:,:) ==-XUNDEF)
WRITE(ILUOUT0,*) 'NB ZREFL > 0 :', COUNT(ZREFL(:,:,:,:,:,:,:)>0.)
WRITE(ILUOUT0,*) 'NB ZREFL = 0 :', COUNT(ZREFL(:,:,:,:,:,:,:)==0.)
WRITE(ILUOUT0,*) 'NB ZREFL < 0 :', COUNT(ZREFL(:,:,:,:,:,:,:) < 0.)-COUNT( ZREFL(:,:,:,:,:,:,:)==XVALGROUND)
!---------------------------------------------------------------------------------------------------
!*       6.    FINAL STEP : TOTAL ATTENUATION AND EQUIVALENT REFLECTIVITY FACTOR
!              ---------------------------------------------------------------
!
ALLOCATE(ZVTEMP(IMAX))
DO JI=1,INBRAD  
  IEL=NBELEV(JI)
  DO JEL=1,IEL  
    DO JAZ=1,INBAZIM 
      IF (LATT) ZAETOT(:,:,1:2)=1.
      PZE(JI,JEL,JAZ,1,IPDP)=0
      DO JL=1,INBSTEPMAX
        ! if no undef point in gate JL and at least one point where T is defined
        IF(COUNT(ZREFL(JI,JEL,JAZ,JL,:,:,1)==-XUNDEF)==0.AND. &
           COUNT(ZREFL(JI,JEL,JAZ,JL,:,:,1)==XVALGROUND)==0.AND. &
           COUNT(PT_RAY(JI,JEL,JAZ,JL,:,:)/=-XUNDEF)/=0) THEN 
          DO JH=1,INPTS_H
            ZVTEMP(:)=0.
            DO JV=1,INPTS_V  ! Loop on Jv
              !if range is over 1, attenuation is added 
              IF (JL > 1) THEN      
                IF(LATT) THEN ! we use ZALPHAE0=alpha_0 from last gate
                  !Total attenuation
                  ZAETOT(JH,JV,1:2)=ZAETOT(JH,JV,1:2)*EXP(-2.*ZAELOC(JI,JEL,JAZ,JL-1,JH,JV,:)*XSTEP_RAD)
                  !Zhh, Zvv
                  ZREFL(JI,JEL,JAZ,JL,JH,JV,1:2)=ZREFL(JI,JEL,JAZ,JL,JH,JV,1:2)*ZAETOT(JH,JV,1:2)!attenuated reflectivity
                  !Z for Radial velocity                         
                  IF(LWREFL) ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)*ZAETOT(JH,JV,1)
                END IF !end IF(LATT)
              END IF !end IF (JL > 1)
              IF(.NOT.(LWREFL.AND.LWBSCS)) THEN
                ZREFL(JI,JEL,JAZ,JL,JH,JV,IVDOP)=PVDOP_RAY(JI,JEL,JAZ,JL,JH,JV)
              END IF
              ! Quadrature on vertical reflectivities +VDOP
              IF(LQUAD) THEN
                ZVTEMP(:)=ZVTEMP(:)+ZREFL(JI,JEL,JAZ,JL,JH,JV,:)*PW_V(ABS((2*JV-INPTS_V-1)/2)+1) &
                          *EXP(-2.*LOG(2.)*PX_V(ABS((2*JV-INPTS_V-1)/2)+1)**2)
              ELSE
                ZVTEMP(:)=ZVTEMP(:)+ZREFL(JI,JEL,JAZ,JL,JH,JV,:)*PW_V(ABS((2*JV-INPTS_V-1)/2)+1)
              END IF
            END DO ! End loop on JV                  		
!
            IF(LQUAD) THEN
              PZE(JI,JEL,JAZ,JL,:)=PZE(JI,JEL,JAZ,JL,:)+ZVTEMP(1:SIZE(PZE,5))*PW_H(ABS((2*JH-INPTS_H-1)/2)+1) &
                                   *EXP(-2.*LOG(2.)*PX_H(ABS((2*JH-INPTS_H-1)/2)+1)**2)
              IF(LWBSCS) ZCONC_BIN(JI,JEL,JAZ,JL)=ZCONC_BIN(JI,JEL,JAZ,JL)+ZVTEMP(IMAX)* &
              PW_H(ABS((2*JH-INPTS_H-1)/2)+1)*EXP(-2.*LOG(2.)*PX_H(ABS((2*JH-INPTS_H-1)/2)+1)**2)
            ELSE
              PZE(JI,JEL,JAZ,JL,:)=PZE(JI,JEL,JAZ,JL,:)+ZVTEMP(1:SIZE(PZE,5))*PW_H(ABS((2*JH-INPTS_H-1)/2)+1)
              IF(LWBSCS) ZCONC_BIN(JI,JEL,JAZ,JL)=ZCONC_BIN(JI,JEL,JAZ,JL)+ZVTEMP(IMAX)* &
                                                   PW_H(ABS((2*JH-INPTS_H-1)/2)+1)
            END IF !end IF(LQUAD)
          END DO ! End loop on JH 

          IF(LQUAD) THEN
            PZE(JI,JEL,JAZ,JL,:)=PZE(JI,JEL,JAZ,JL,:)*2.*LOG(2.)/XPI
            IF(LWBSCS) ZCONC_BIN(JI,JEL,JAZ,JL)=ZCONC_BIN(JI,JEL,JAZ,JL)*2.*LOG(2.)/XPI
          ELSE
            PZE(JI,JEL,JAZ,JL,:)=PZE(JI,JEL,JAZ,JL,:)/XPI
            IF(LWBSCS) ZCONC_BIN(JI,JEL,JAZ,JL)=ZCONC_BIN(JI,JEL,JAZ,JL)/XPI
          END IF !end IF(LQUAD)
! 
          !**** Thresholding: with ZSNR, or with XREFLVDOPMIN and XREFLMIN
          ZSNR=-XUNDEF
          ZSNR_R=-XUNDEF
          ZSNR_I=-XUNDEF
          ZSNR_S=-XUNDEF
          ZSNR_G=-XUNDEF
          ZSNR_H=-XUNDEF
          ZZHH=PZE(JI,JEL,JAZ,JL,1)
          ZZE_R=PZE(JI,JEL,JAZ,JL,IZER)
          ZZE_I=PZE(JI,JEL,JAZ,JL,IZEI)
          ZZE_S=PZE(JI,JEL,JAZ,JL,IZES)
          ZZE_G=PZE(JI,JEL,JAZ,JL,IZEG)
          IF (GHAIL) ZZE_H=PZE(JI,JEL,JAZ,JL,IZEH)
          ZDISTRAD=JL*XSTEP_RAD !radar distance in meters
          IF (LSNRT) THEN 
            IF (ZZHH/=XVALGROUND .AND. ZZHH/=-XUNDEF.AND.ZZHH/=0) THEN
              ZSNR=10*LOG10(ZZHH)-20*LOG10(ZDISTRAD/(100*10**3))
            END IF
            IF (ZZE_R/=XVALGROUND .AND. ZZE_R/=-XUNDEF.AND.ZZE_R/=0) THEN
              ZSNR_R=10*LOG10(ZZE_R)-20*LOG10(ZDISTRAD/(100*10**3))
            END IF
            IF (ZZE_I/=XVALGROUND .AND. ZZE_I/=-XUNDEF.AND.ZZE_I/=0) THEN
              ZSNR_I=10*LOG10(ZZE_I)-20*LOG10(ZDISTRAD/(100*10**3))
            END IF
            IF (ZZE_S/=XVALGROUND .AND. ZZE_S/=-XUNDEF.AND.ZZE_S/=0) THEN
              ZSNR_S=10*LOG10(ZZE_S)-20*LOG10(ZDISTRAD/(100*10**3))
            END IF
            IF (ZZE_G/=XVALGROUND .AND. ZZE_G/=-XUNDEF.AND.ZZE_G/=0) THEN
              ZSNR_G=10*LOG10(ZZE_G)-20*LOG10(ZDISTRAD/(100*10**3))
            END IF
            IF (GHAIL) THEN
              IF (ZZE_H/=XVALGROUND .AND. ZZE_H/=-XUNDEF.AND.ZZE_H/=0) THEN
                ZSNR_H=10*LOG10(ZZE_H)-20*LOG10(ZDISTRAD/(100*10**3))
              END IF
            END IF
            GTHRESHOLD_V=(ZSNR>=XSNRMIN)
            GTHRESHOLD_Z=GTHRESHOLD_V
            GTHRESHOLD_ZR=(ZSNR_R>=XSNRMIN)
            GTHRESHOLD_ZI=(ZSNR_I>=XSNRMIN)
            GTHRESHOLD_ZS=(ZSNR_S>=XSNRMIN)
            GTHRESHOLD_ZG=(ZSNR_G>=XSNRMIN)
            IF (GHAIL) GTHRESHOLD_ZH=(ZSNR_H>=XSNRMIN)
          ELSE
            GTHRESHOLD_V=(ZZHH>=10**(XREFLVDOPMIN/10.))
            GTHRESHOLD_Z=(ZZHH>=10**(XREFLMIN/10.))
            GTHRESHOLD_ZR=(ZZE_R>=10**(XREFLMIN/10.))
            GTHRESHOLD_ZI=(ZZE_I>=10**(XREFLMIN/10.))
            GTHRESHOLD_ZS=(ZZE_S>=10**(XREFLMIN/10.))
            GTHRESHOLD_ZG=(ZZE_G>=10**(XREFLMIN/10.))
            IF (GHAIL) GTHRESHOLD_ZH=(ZZE_H>=10**(XREFLMIN/10.))
          END IF          
          !--- Doppler velocities 
          IF(GTHRESHOLD_V) THEN      
            IF(LWREFL) THEN
              !change Clotilde 27/04/2012 to avoid division by zero and floating point exception
              IF (PZE(JI,JEL,JAZ,JL,1)/=0) THEN
                PZE(JI,JEL,JAZ,JL,IVDOP)=PZE(JI,JEL,JAZ,JL,IVDOP)/PZE(JI,JEL,JAZ,JL,1) 
              END IF
            ELSE IF(LWBSCS) THEN
              IF(ZCONC_BIN(JI,JEL,JAZ,JL)>0.) THEN
                PZE(JI,JEL,JAZ,JL,IVDOP)=PZE(JI,JEL,JAZ,JL,IVDOP)/ZCONC_BIN(JI,JEL,JAZ,JL)
              ELSE
                PZE(JI,JEL,JAZ,JL,IVDOP)=-XUNDEF
              END IF !end IF(ZCONC_BIN(JI,JEL,JAZ,JL)>0.)
            END IF !end IF(LWREFL)
          ELSE
            PZE(JI,JEL,JAZ,JL,IVDOP)=-XUNDEF        
          END IF !end IF(GTHRESHOLD_V)   
      
          !--- Zhh, Zvv et variables globales
          IF(GTHRESHOLD_Z .EQV. .FALSE.) THEN
            PZE(JI,JEL,JAZ,JL,1:4)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IRHV:IDHV)=-XUNDEF
          END IF
          !--- ZER, ZDA, KDR, RHR
          IF(GTHRESHOLD_ZR .EQV. .FALSE.) THEN
            PZE(JI,JEL,JAZ,JL,IZER)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IZDA)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IKDR)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IRHR)=-XUNDEF
          END IF
          !--- ZES, ZDS, KDS, RHS
          IF(GTHRESHOLD_ZS .EQV. .FALSE.) THEN
            PZE(JI,JEL,JAZ,JL,IZES)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IZDS)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IKDS)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IRHS)=-XUNDEF
          END IF

          !--- ZEG, ZDG, KDG, RHG
          IF(GTHRESHOLD_ZG .EQV. .FALSE.) THEN
            PZE(JI,JEL,JAZ,JL,IZEG)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IZDG)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IKDG)=-XUNDEF
            PZE(JI,JEL,JAZ,JL,IRHG)=-XUNDEF
          END IF
          !--- ZEH, ZDH, KDH, RHH
          IF (GHAIL) THEN
            IF(GTHRESHOLD_ZH .EQV. .FALSE.) THEN
              PZE(JI,JEL,JAZ,JL,IZEH)=-XUNDEF
              PZE(JI,JEL,JAZ,JL,IZDH)=-XUNDEF
              PZE(JI,JEL,JAZ,JL,IKDH)=-XUNDEF
              PZE(JI,JEL,JAZ,JL,IRHH)=-XUNDEF
            END IF
          END IF
          !--- ZEI
          IF(GTHRESHOLD_ZI .EQV. .FALSE.) THEN
            PZE(JI,JEL,JAZ,JL,IZEI)=-XUNDEF
          END IF     
        ELSE 
          ! ground clutter or outside Meso-NH domain 
          !(IF T not defined or if one undef point at least in gate)
          PZE(JI,JEL,JAZ,JL,:)=XVALGROUND
        END IF 
        IF(PZE(JI,JEL,JAZ,JL,1) < 0. .AND. PZE(JI,JEL,JAZ,JL,1)/=-XUNDEF) THEN  ! flag bin when underground => xvalground si < 0?                 
          PZE(JI,JEL,JAZ,JL,:)=XVALGROUND
        END IF ! end IF(PZE(JI,JEL,JAZ,JL,1) < 0.)
      END DO ! end DO JL=1,INBSTEPMAX
    END DO !end DO JAZ=1,INBAZIM
  END DO !end DO JEL=1,IEL
END DO !end DO JI=1,INBRAD
DEALLOCATE(ZREFL,ZVTEMP,ZRTMIN)
WRITE(ILUOUT0,*) '*****************FIN RADAR_SCATTERING ***********************' 
WRITE(ILUOUT0,*) 'NB PZE MIN MAX :', MINVAL(PZE(:,:,:,:,IZEI)),MAXVAL(PZE(:,:,:,:,IZEI))
WRITE(ILUOUT0,*) 'NB PZE VALGROUND :', COUNT(PZE(:,:,:,:,IZEI) ==XVALGROUND)
WRITE(ILUOUT0,*) 'NB PZE -XUNDEF :', COUNT(PZE(:,:,:,:,IZEI) ==-XUNDEF)
WRITE(ILUOUT0,*) 'NB PZE > 0 :', COUNT(PZE(:,:,:,:,IZEI)>0.)
WRITE(ILUOUT0,*) 'NB PZE = 0 :', COUNT(PZE(:,:,:,:,IZEI)==0.)
WRITE(ILUOUT0,*) 'NB PZE < 0 :', COUNT(PZE(:,:,:,:,IZEI) < 0.)-COUNT(PZE(:,:,:,:,IZEI) ==XVALGROUND)
IF(NDIFF/=0) DEALLOCATE(ZX,ZW)
IF (LATT) DEALLOCATE(ZAELOC,ZAETOT)
WRITE(ILUOUT0,*) 'END OF RADAR SCATTERING'
END SUBROUTINE RADAR_SCATTERING

