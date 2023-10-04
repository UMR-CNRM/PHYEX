!MNH_LIC Copyright 1995-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!    ########################
     MODULE MODI_PHYS_PARAM_n
!    ########################
!
!
INTERFACE
!
      SUBROUTINE PHYS_PARAM_n( KTCOUNT, TPFILE,                                              &
                               PRAD, PSHADOWS, PKAFR, PGROUND, PMAFL, PDRAG,PEOL, PTURB,     &
                               PTRACER, PTIME_BU, PWETDEPAER, OMASKkids, OCLOUD_ONLY         )
!
USE MODD_IO,        ONLY: TFILEDATA
use modd_precision, only: MNHTIME
!
INTEGER,           INTENT(IN)     :: KTCOUNT   ! temporal iteration count
TYPE(TFILEDATA),   INTENT(IN)     :: TPFILE    ! Synchronous output file
! advection schemes
REAL(kind=MNHTIME), DIMENSION(2), INTENT(INOUT) :: PRAD,PSHADOWS,PKAFR,PGROUND,PTURB,PMAFL,PDRAG,PTRACER,PEOL ! to store CPU
                                                                                                         ! time for computing time
REAL(kind=MNHTIME), DIMENSION(2), INTENT(INOUT) :: PTIME_BU  ! time used in budget&LES budgets statistics
REAL, DIMENSION(:,:,:,:), INTENT(INOUT)  :: PWETDEPAER
LOGICAL, DIMENSION(:,:), INTENT(IN) :: OMASKkids ! kids domains mask
LOGICAL, INTENT(OUT) :: OCLOUD_ONLY ! conditionnal radiation computations for
                                !      the only cloudy columns
                                !
END SUBROUTINE PHYS_PARAM_n
!
END INTERFACE
!
END MODULE MODI_PHYS_PARAM_n
!
!     ########################################################################################
      SUBROUTINE PHYS_PARAM_n( KTCOUNT, TPFILE,                                              &
                               PRAD, PSHADOWS, PKAFR, PGROUND, PMAFL, PEOL, PDRAG, PTURB,    &
                               PTRACER, PTIME_BU, PWETDEPAER, OMASKkids, OCLOUD_ONLY         )
!     ########################################################################################
!
!!****  *PHYS_PARAM_n * -monitor of the parameterizations used by model _n
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to update the sources by adding the
!     parameterized terms. This is realized by sequentially calling the
!     specialized routines.
!    
!!**  METHOD
!!    ------
!!      The first parametrization is the radiation scheme:
!!                                       ----------------
!!     *  CRAD = 'FIXE'
!!     In this case, a temporal interpolation  is performed for the downward
!!     surface fluxes XFLALWD and XFLASWD.
!!     *  CRAD = 'ECMWF'
!!     Several tests are performed before calling the radiation computations
!!     interface with the ECMWF radiation scheme code. A control is made to
!!     ensure that:
!!         - the full radiation code is called at the first model timestep
!!         - there is a priority for calling the full radiation instead of the
!!           cloud-only approximation if both must be called at the current
!!           timestep
!!         - the cloud-only option (approximation) is coherent with the
!!           occurence of one cloudy vertical column at least
!!      If all the above conditions are fulfilled (GRAD is .TRUE.) then the
!!     position of the sun is computed in routine SUNPOS_n and the interfacing
!!     routine RADIATIONS is called to update the radiative tendency XDTHRAD
!!     and the downward surface fluxes XFLALWD and XFLASWD. Finally, the
!!     radiative tendency is integrated as a source term in the THETA prognostic
!!     equation.
!!
!!      The second parameterization is the soil scheme:
!!                                         -----------
!!
!!     externalized surface
!!
!!       The third parameterization is the turbulence scheme:
!!                                         -----------------
!!     * CTURB='NONE'
!!     no turbulent mixing is taken into account
!!     * CTURB='TKEL'
!!     The turbulent fluxes are computed according to a one and half order
!!     closure of the hydrodynamical equations. This scheme is based on a
!!     prognostic for the turbulent kinetic energy and a mixing length
!!     computation ( the mesh size or a physically based length). Other
!!     turbulent moments are diagnosed according to a stationarization of the
!!     second order turbulent moments. This turbulent scheme forecasts
!!     either a purely vertical turbulent mixing or 3-dimensional mixing
!!     according to its internal degrees of freedom.
!!
!!
!!       The LAST parameterization is the chemistry scheme:
!!                                        -----------------
!!     The chemistry part of MesoNH has two namelists, NAM_SOLVER for the
!!     parameters concerning the stiff solver, and NAM_MNHCn concerning the
!!     configuration and options of the chemistry module itself.
!!     The switch LUSECHEM in NAM_CONF acitvates or deactivates the chemistry.
!!     The only variables of MesoNH that are modified by chemistry are the
!!     scalar variables. If calculation of chemical surface fluxes is
!!     requested, those fluxes are calculated before
!!     entering the turbulence scheme, since those fluxes are taken into
!!     account by TURB as surface boundary conditions.
!!     CAUTION: chemistry has allways to be called AFTER ALL OTHER TERMS
!!     that affect the scalar variables (dynamical terms, forcing,
!!     parameterizations (like TURB, CONVECTION), since it uses the variables
!!     XRSVS as input in case of the time-split option.
!!
!!    EXTERNAL
!!    --------
!!      Subroutine SUNPOS_n     : computes the position of the sun
!!      Subroutine RADIATIONS   : computes the radiative tendency and fluxes
!!      Subroutine TSZ0         : computes the surface from temporally
!!                                interpolated Ts and given z0
!!      Subroutine ISBA         : computes the surface fluxes from a soil scheme
!!      Subroutine TURB         : computes the turbulence source terms
!!      Subroutine CONVECTION   : computes the convection source term
!!      Subroutine CH_SURFACE_FLUX_n: computes the surface flux for chemical
!!                                species
!!      Subroutine CH_MONITOR_n : computes the chemistry source terms
!!                                that are applied to the scalar variables
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      USE MODD_DYN
!!      USE MODD_CONF
!!      USE MODD_CONF_n
!!      USE MODD_CURVCOR_n
!!      USE MODD_DYN_n
!!      USE MODD_FIELD_n
!!      USE MODD_GR_FIELD_n
!!      USE MODD_LSFIELD_n
!!      USE MODD_GRID_n
!!      USE MODD_LBC_n
!!      USE MODD_PARAM_RAD_n
!!      USE MODD_RADIATIONS_n
!!      USE MODD_REF_n
!!      USE MODD_LUNIT_n
!!      USE MODD_TIME_n
!!      USE MODD_CH_MNHC_n
!!
!!    REFERENCE
!!    ---------
!!      None
!!
!!    AUTHOR
!!    ------
!!      J. Stein           * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/01/95
!!      Modifications  Feb 14, 1995 (J.Cuxart)  add the I/O arguments,
!!             the director cosinus and change the names of the surface fluxes
!!      Modifications March 21, 1995 (J.M.Carriere) take into account liquid
!!                                             water
!!                    June 30,1995  (J.Stein)  initialize at 0 the surf. fluxes
!!      Modifications Sept. 1, 1995 (S.Belair) ISBA scheme
!!      Modifications Sept.25, 1995 (J.Stein)  switch on the radiation scheme
!!      Modifications Sept. 11, 1995 (J.-P. Pinty) radiation scheme
!!                    Nov.  15, 1995 (J.Stein) cleaning + change the temporal
!!                                   algorithm for the soil scheme-turbulence
!!                    Jan.  23, 1996 (J.Stein) add a new option for the surface
!!                                   fluxes where Ts and z0 are given
!!                    March 18, 1996 (J.Stein) add the cloud fraction
!!                    March 28, 1996 (J.Stein) the soil scheme gives energy
!!                                             fluxes + cleaning
!!                    June  17, 1996 (Lafore)  statistics of computing time
!!                    August 4, 1996 (K. Suhre) add chemistry
!!                    Oct.  12, 1996 (J.Stein) use XSRCM in the turbulence
!!                                             scheme
!!                    Nov.  18, 1996 (J.-P. Pinty) add domain translation
!!                                                 change arg. in radiations
!!                    Fev.   4, 1997 (J.Viviand) change isba's calling for ice
!!                    Jun.  22, 1997 (J.Stein) change the equation system and use
!!                                             the absolute pressure
!!                    Jul.  09, 1997 (V.Masson) add directional z0
!!                    Jan.  24, 1998 (P.Bechtold) add convective transport for tracers
!!                    Jan.  24, 1998 (J.-P. Pinty) split SW and LW part for radiation
!!                    Mai.  10, 1999 (P.Bechtold) shallow convection
!!                    Oct.  20, 1999 (P.Jabouille) domain translation for turbulence
!!                    Jan.  04, 2000 (V.Masson) removes TSZ0 case
!!                    Jan.  04, 2000 (V.Masson) modifies albedo computation
!                     Jul   02, 2000 (F.Solmon/V.Masson) adaptation for patch approach
!!                    Nov.  15, 2000 (V.Masson) LES routines
!!                    Nov.  15, 2000 (V.Masson) effect of slopes on surface fluxes
!!                    Feb.  02, 2001 (P.Tulet) add friction velocities and aerodynamical
!!                                             resistance (patch approach)
!!                    Jan.  04, 2000 (V.Masson) modify surf_rad_modif computation
!!                    Mar.  04, 2002 (F.Solmon) new interface for radiation call
!!                    Nov.  06, 2002 (V.Masson) LES budgets & budget time counters
!!                    Jan. 2004      (V.Masson) surface externalization
!!                    Jan.  13, 2004 (J.Escobar) bug correction : compute "GRAD" in parallel
!!                    Jan.  20, 2005 (P. Tulet)  add dust sedimentation 
!!                    Jan.  20, 2005 (P. Tulet)  climatologic SSA
!!                    Jan.  20, 2005 (P. Tulet)  add aerosol / dust scavenging
!!                    Jul. 2005       (N. Asencio) use the two-way result-fields
!!                                  before ground_param call
!!                    May 2006        Remove EPS
!!                    Oct. 2007      (J.Pergaud) Add shallow_MF
!!                    Oct. 2009     (C.Lac) Introduction of different PTSTEP according to the
!!                              advection schemes
!!                    Oct. 2009     (V. MAsson) optimization of Pergaud et al massflux scheme
!!                    Aug. 2010     (V.Masson, C.Lac) Exchange of SBL_DEPTH for
!!                                  reproducibility
!!                    Oct. 2010   (J.Escobar) init  ZTIME_LES_MF ( pb detected with g95 )
!!                    Feb. 2011 (V.Masson, C.Lac) SBL_DEPTH values on outer pts
!!                               for RMC01
!!                    Sept.2011 (J.Escobar) init YINST_SFU ='M'
!!
!!                        Specific for 2D modeling : 
!! 
!!                    06/2010    (P.Peyrille)  add Call to aerozon.f90 if LAERO_FT=T
!!                                to update 
!!                                aerosols and ozone climatology at each call to
!!                                phys_param otherwise it is constant to monthly average
!!                    03/2013  (C.Lac) FIT temporal scheme
!!                    01/2014 (C.Lac) correction for the nesting of 2D surface
!!                           fields if the number of the son model does not
!!                           follow the number of the dad model
!!      J.Escobar 21/03/2013: for HALOK comment all NHALO=1 test
!!                       2014  (M.Faivre)
!!  06/2016     (G.Delautier) phasage surfex 8
!!  2016 B.VIE LIMA
!!      M. Leriche 02/2017 Avoid negative fluxes if sv=0 outside the physics domain
!!      C.Lac  10/2017 : ch_monitor and aer_monitor extracted from phys_param
!!                       to be called directly by modeln as the last process 
!!                   02/2018 Q.Libois ECRAD
!  P. Wautelet 28/03/2018: replace TEMPORAL_DIST by DATETIME_DISTANCE
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 28/03/2019: use MNHTIME for time measurement variables
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  P. Wautelet 21/11/2019: ZRG_HOUR and ZRAT_HOUR are now parameter arrays
!  C. Lac         11/2019: correction in the drag formula and application to building in addition to tree
!  F. Auguste     02/2021: add IBM
!  JL Redelsperger 03/2021: add the SW flux penetration for Ocean model case
!  R. Schoetter    12/2021: multi-level coupling between MesoNH and SURFEX  
!  P. Wautelet 30/11/2022: compute XTHW_FLUX, XRCW_FLUX and XSVW_FLUX only when needed
!  A. Costes      12/2021: add Blaze fire model
!  Q. Rodier      2022   : integration with PHYEX
!!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_ADV_n,       ONLY : XRTKEMS
USE MODD_AIRCRAFT_BALLOON, ONLY: LFLYER
USE MODD_ARGSLIST_ll, ONLY : LIST_ll
USE MODD_BLOWSNOW,    ONLY : LBLOWSNOW,XRSNOW
USE MODD_BUDGET,      ONLY: NBUDGET_TH, NBUDGET_RV, NBUDGET_RC, NBUDGET_RI, NBUDGET_SV1, &
                            TBUDGETS, xtime_bu_process, TBUCONF
USE MODD_CH_AEROSOL
USE MODD_CH_MNHC_n, ONLY : LUSECHEM,         &! indicates if chemistry is used
                           LCH_CONV_SCAV,    &
                           LCH_CONV_LINOX
USE MODD_CLOUD_MF_n
USE MODD_CONDSAMP
USE MODD_CONF
USE MODD_CONF_n
USE MODD_CST, ONLY : CST
USE MODD_CTURB, ONLY : CSTURB
USE MODD_CURVCOR_n
USE MODD_DEEP_CONVECTION_n
USE MODD_DEF_EDDY_FLUX_n           ! Ajout PP
USE MODD_DEF_EDDYUV_FLUX_n         ! Ajout PP
USE MODD_DIAG_IN_RUN, ONLY: LDIAG_IN_RUN, XCURRENT_TKE_DISS
USE MODD_DIM_n, ONLY: NIMAX_ll, NJMAX_ll
USE MODD_DRAGBLDG_n
USE MODD_DRAGTREE_n
USE MODD_DUST
USE MODD_DYN
USE MODD_DYN_n
USE MODD_EOL_MAIN, ONLY: LMAIN_EOL, CMETH_EOL, NMODEL_EOL
USE MODD_FIELD_n
USE MODD_FRC
USE MODD_FRC_n
USE MODD_GRID
USE MODD_GRID_n
USE MODD_IBM_PARAM_n,      ONLY: LIBM, XIBM_EPSI, XIBM_LS, XIBM_XMUT
USE MODD_ICE_C1R3_DESCR,  ONLY : XRTMIN_C1R3=>XRTMIN
USE MODD_IO, ONLY: TFILEDATA
USE MODD_LATZ_EDFLX
USE MODD_LBC_n
USE MODD_LES
USE MODD_LES_n, ONLY: NLES_TIMES
USE MODD_LES_BUDGET
USE MODD_LSFIELD_n
USE MODD_LUNIT_n
USE MODD_METRICS_n
USE MODD_MNH_SURFEX_n
USE MODD_NESTING, ONLY : XWAY,NDAD, NDXRATIO_ALL, NDYRATIO_ALL
USE MODD_NSV, ONLY : NSV, NSV_LGBEG, NSV_LGEND, &
                     NSV_SLTBEG,NSV_SLTEND,NSV_SLT,&
                     NSV_AERBEG,NSV_AEREND, &
                     NSV_DSTBEG,NSV_DSTEND, NSV_DST,&
                     NSV_LIMA_NR,NSV_LIMA_NS,NSV_LIMA_NG,NSV_LIMA_NH
USE MODD_OCEANH
USE MODD_OUT_n
USE MODD_PARAM_C2R2,       ONLY : LSEDC
USE MODD_PARAMETERS
USE MODD_PARAM_ICE_n,        ONLY : LSEDIC
USE MODD_PARAM_KAFR_n
USE MODD_PARAM_LIMA,       ONLY : MSEDC => LSEDC, XRTMIN_LIMA=>XRTMIN
USE MODD_PARAM_MFSHALL_n,  ONLY: CMF_CLOUD
USE MODD_PARAM_n
USE MODD_PARAM_RAD_n
USE MODD_PASPOL
USE MODD_PASPOL_n
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_PRECIP_n
use modd_precision,        only: MNHTIME
USE MODD_RADIATIONS_n
USE MODD_RAIN_ICE_DESCR_n,   ONLY: XRTMIN
USE MODD_REF,              ONLY: LCOUPLES
USE MODD_REF_n
USE MODD_SALT
USE MODD_SHADOWS_n
USE MODD_SUB_PHYS_PARAM_n
USE MODD_TIME_n
USE MODD_TIME_n
USE MODD_TIME, ONLY : TDTEXP  ! Ajout PP
USE MODD_TURB_FLUX_AIRCRAFT_BALLOON, ONLY : XTHW_FLUX, XRCW_FLUX, XSVW_FLUX
USE MODD_TURB_n
USE MODD_NEB_n, ONLY: NEBN

USE MODE_AERO_PSD
use mode_budget,            only: Budget_store_end, Budget_store_init
USE MODE_DATETIME
USE MODE_DUST_PSD
USE MODE_ll
USE MODE_GATHER_ll
USE MODE_MNH_TIMING
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
USE MODE_SALT_PSD

USE MODI_AEROZON          ! Ajout PP
USE MODI_CONDSAMP
USE MODI_CONVECTION
USE MODI_DRAG_BLD
USE MODI_DRAG_VEG
USE MODI_DUST_FILTER
USE MODI_EDDY_FLUX_n               ! Ajout PP
USE MODI_EDDY_FLUX_ONE_WAY_n       ! Ajout PP
USE MODI_EDDYUV_FLUX_n             ! Ajout PP
USE MODI_EDDYUV_FLUX_ONE_WAY_n     ! Ajout PP
USE MODI_EOL_MAIN
USE MODI_GROUND_PARAM_n
USE MODI_GRADIENT_M
USE MODI_GRADIENT_W
USE MODI_PASPOL
USE MODI_RADIATIONS
USE MODI_SALT_FILTER
USE MODI_SEDIM_DUST
USE MODI_SEDIM_SALT
USE MODI_SHALLOW_MF_PACK
USE MODI_SUNPOS_n
USE MODI_SURF_RAD_MODIF
USE MODI_SWITCH_SBG_LES_N
USE MODI_TURB

IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
INTEGER,           INTENT(IN)     :: KTCOUNT   ! temporal iteration count
TYPE(TFILEDATA),   INTENT(IN)     :: TPFILE    ! Synchronous output file
! advection schemes
REAL(kind=MNHTIME), DIMENSION(2), INTENT(INOUT) :: PRAD,PSHADOWS,PKAFR,PGROUND,PTURB,PMAFL,PDRAG,PTRACER,PEOL ! to store CPU
                                                                                                         ! time for computing time
REAL(kind=MNHTIME), DIMENSION(2), INTENT(INOUT) :: PTIME_BU  ! time used in budget&LES budgets statistics
REAL, DIMENSION(:,:,:,:), INTENT(INOUT)  :: PWETDEPAER
LOGICAL, DIMENSION(:,:), INTENT(IN) :: OMASKkids ! kids domains mask
LOGICAL, INTENT(OUT) :: OCLOUD_ONLY ! conditionnal radiation computations for
                                !      the only cloudy columns
                                !
!
!*      0.2    declarations of local variables
!
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFU  ! surface flux of x and
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFV  ! y component of wind
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFTH ! surface flux of theta
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFRV ! surface flux of vapor
REAL, DIMENSION(:,:,:), ALLOCATABLE   :: ZSFSV ! surface flux of scalars
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFCO2! surface flux of CO2
!
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFTH_WALL
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFTH_ROOF
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZCD_ROOF
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFRV_WALL 
REAL, DIMENSION(:,:), ALLOCATABLE     :: ZSFRV_ROOF
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDIR_ALB ! direct albedo
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZSCA_ALB ! diffuse albedo
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZEMIS    ! emissivity
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZTSRAD   ! surface temperature
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZRGDST,ZSIGDST,ZNDST,ZSVDST
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZRGSLT,ZSIGSLT,ZNSLT,ZSVSLT
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZRGAER,ZSIGAER,ZNAER,ZSVAER
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZSVT
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZEXN   ! Atmospheric density and Exner
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZSIGMF   ! MF contribution to XSIGS
!
REAL, DIMENSION(0:24), parameter :: ZRG_HOUR =  (/ 0., 0., 0., 0., 0., 32.04, 114.19, &
                                                   228.01, 351.25, 465.49, 557.24,    &
                                                   616.82, 638.33, 619.43, 566.56,    &
                                                   474.71, 359.20, 230.87, 115.72,    &
                                                   32.48, 0., 0., 0., 0., 0. /)
!
REAL, DIMENSION(0:24), parameter :: ZRAT_HOUR = (/ 326.00, 325.93, 325.12, 324.41, &
                                                   323.16, 321.95, 322.51, 325.16, &
                                                   328.01, 331.46, 335.58, 340.00, &
                                                   345.20, 350.32, 354.20, 356.58, &
                                                   356.56, 355.33, 352.79, 351.34, &
                                                   347.00, 342.00, 337.00, 332.00, &
                                                   326.00     /)
!
!
character(len=6) :: ynum
INTEGER  :: IHOUR               ! parameters necessary for the temporal
REAL     :: ZTIME, ZDT          ! interpolation
REAL     :: ZTEMP_DIST          ! time between 2 instants (in seconds)
!
LOGICAL :: GRAD                 ! conditionnal call for the full radiation
                                !         computations
REAL    :: ZRAD_GLOB_ll         ! 'real' global parallel mask of 'GRAD'
INTEGER :: INFO_ll              ! error report of parallel routines
                                !      the only cloudy columns
!
REAL(kind=MNHTIME), DIMENSION(2) :: ZTIME1, ZTIME2, ZTIME3, ZTIME4 ! for computing time analysis
REAL(kind=MNHTIME), DIMENSION(2) :: ZTIME_LES_MF                   ! time spent in LES computation in shallow conv.
LOGICAL :: GDCONV               ! conditionnal call for the deep convection
                                !         computations
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZRC, ZRI, ZWT ! additional dummies
REAL, DIMENSION(:,:),   ALLOCATABLE  :: ZDXDY         ! grid area
                    ! for rc, ri, w required if main variables not allocated
!
INTEGER :: IIU, IJU, IKU                              ! dimensional indexes
!
INTEGER     :: JSV              ! Loop index for Scalar Variables
INTEGER     :: JSWB             ! loop on SW spectral bands
INTEGER     :: IIB,IIE,IJB,IJE, IKB, IKE, JI,JJ
INTEGER     :: IMODEIDX
              ! index values for the Beginning or the End of the physical
              ! domain in x and y directions
TYPE(LIST_ll), POINTER :: TZFIELDS_ll   ! list of fields to exchange
INTEGER                :: IINFO_ll       ! return code of parallel routine
!
!* variables for writing in a fm file
!
INTEGER           :: IRESP          ! IRESP  : return-code if a problem appears
                                    !in LFI subroutines at the open of the file
INTEGER           :: ILUOUT         ! logical unit numbers of output-listing
INTEGER           :: IMI            ! model index
INTEGER           :: JKID           ! loop index to look for the KID models
REAL              :: ZINIRADIUSI, ZINIRADIUSJ ! ORILAM initial radius
REAL, DIMENSION(NMODE_DST)    :: ZINIRADIUS  ! DUST initial radius
REAL, DIMENSION(NMODE_SLT)    :: ZINIRADIUS_SLT  ! Sea Salt initial radius
REAL, DIMENSION(SIZE(XRSVS,1), SIZE(XRSVS,2), SIZE(XRSVS,3), SIZE(XRSVS,4))  :: ZRSVS
LOGICAL :: GCLD                     ! conditionnal call for dust wet deposition
! * arrays to store the surface fields before radiation and convection scheme
!  calls
INTEGER           :: IMODSON        ! Number of son models of IMI with XWAY=2
INTEGER           :: IKIDM          ! index loop                                 
INTEGER           :: IGRADIENTS     ! Number of horizontal gradients in turb
REAL, DIMENSION(:,:,:),   ALLOCATABLE  :: ZSAVE_INPRR,ZSAVE_INPRS,ZSAVE_INPRG,ZSAVE_INPRH
REAL, DIMENSION(:,:,:),   ALLOCATABLE  :: ZSAVE_INPRC,ZSAVE_PRCONV,ZSAVE_PRSCONV
REAL, DIMENSION(:,:,:,:), ALLOCATABLE  :: ZSAVE_DIRFLASWD, ZSAVE_SCAFLASWD,ZSAVE_DIRSRFSWD
! for ocean model
INTEGER           :: JKM , JSW         ! vertical index loop                                 
REAL :: ZSWA,TINTSW     ! index for SW interpolation and int time betwenn forcings (ocean model)
REAL, DIMENSION(:), ALLOCATABLE :: ZIZOCE(:) ! Solar flux penetrating in ocean
REAL, DIMENSION(:), ALLOCATABLE :: ZPROSOL1(:),ZPROSOL2(:) ! Funtions for penetrating solar flux
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLENGTHM, ZLENGTHH, ZMFMOIST !OHARAT turb option from AROME (not allocated in MNH)
                                                                    ! to be moved as optional args for turb
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTDIFF, ZTDISS
REAL, DIMENSION(:),ALLOCATABLE  :: ZXHAT_ll,ZYHAT_ll  !  Position x/y in the conformal
                                                 ! plane (array on the complete domain)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZDIST ! distance from the center of the cooling 
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZHGRAD ! horizontal gradient used in turb
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
LOGICAL :: GCOMPUTE_SRC ! flag to define dimensions of SIGS and SRCT variables 
!-----------------------------------------------------------------------------

NULLIFY(TZFIELDS_ll)
IMI=GET_CURRENT_MODEL_INDEX()
!
ILUOUT = TLUOUT%NLU
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
IKU=SIZE(XTHT,3)
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL FILL_DIMPHYEX( YLDIMPHYEX, SIZE(XTHT,1), SIZE(XTHT,2), SIZE(XTHT,3), LTURB=.TRUE., KLES_TIMES=NLES_TIMES, KLES_K=NLES_K )
!
ZTIME1 = 0.0_MNHTIME
ZTIME2 = 0.0_MNHTIME
ZTIME3 = 0.0_MNHTIME
ZTIME4 = 0.0_MNHTIME
PTIME_BU = 0._MNHTIME
ZTIME_LES_MF = 0.0_MNHTIME
PWETDEPAER(:,:,:,:) = 0.
!
!* allocation of variables used in more than one parameterization
!
ALLOCATE(ZSFU  (IIU,IJU))         ! surface schemes + turbulence
ALLOCATE(ZSFV  (IIU,IJU))
ALLOCATE(ZSFTH (IIU,IJU))
ALLOCATE(ZSFRV (IIU,IJU))
ALLOCATE(ZSFSV (IIU,IJU,NSV))
ALLOCATE(ZSFCO2(IIU,IJU))
!
ALLOCATE(ZSFTH_WALL (IIU,IJU))
ALLOCATE(ZSFTH_ROOF (IIU,IJU))
ALLOCATE(ZCD_ROOF   (IIU,IJU))
ALLOCATE(ZSFRV_WALL (IIU,IJU))
ALLOCATE(ZSFRV_ROOF (IIU,IJU))
!
!* if XWAY(son)=2 save surface fields before radiation or convective scheme
!  calls
!
IMODSON = 0
DO JKID = IMI+1,NMODEL  ! min value of the possible kids
 IF (IMI == NDAD(JKID) .AND. XWAY(JKID) == 2. .AND. CPROGRAM=='MESONH' &
  .AND. (CCONF == 'RESTA' .OR. (CCONF == 'START' .AND. KTCOUNT /= 1))) THEN
  IMODSON = IMODSON + 1
 END IF
END DO
!
 IF (IMODSON /= 0 ) THEN
   IF (LUSERC .AND. (                                               &
       (LSEDIC .AND. CCLOUD(1:3) == 'ICE')                     .OR. &
       (LSEDC  .AND. (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO')) .OR. &
       (MSEDC  .AND. CCLOUD=='LIMA')                                &
      )) THEN
     ALLOCATE( ZSAVE_INPRC(SIZE(XINPRC,1),SIZE(XINPRC,2),IMODSON))
   ELSE
     ALLOCATE( ZSAVE_INPRC(0,0,0))
   END IF
   IF (LUSERR) THEN
     ALLOCATE( ZSAVE_INPRR(SIZE(XINPRR,1),SIZE(XINPRR,2),IMODSON))
   ELSE
     ALLOCATE( ZSAVE_INPRR(0,0,0))
   END IF
   IF (LUSERS) THEN
     ALLOCATE( ZSAVE_INPRS(SIZE(XINPRS,1),SIZE(XINPRS,2),IMODSON))
   ELSE
     ALLOCATE( ZSAVE_INPRS(0,0,0))                              
   END IF
   IF (LUSERG) THEN
     ALLOCATE( ZSAVE_INPRG(SIZE(XINPRG,1),SIZE(XINPRG,2),IMODSON))
   ELSE
     ALLOCATE( ZSAVE_INPRG(0,0,0))                                
   END IF
   IF (LUSERH) THEN
     ALLOCATE( ZSAVE_INPRH(SIZE(XINPRH,1),SIZE(XINPRH,2),IMODSON))
   ELSE
     ALLOCATE( ZSAVE_INPRH(0,0,0))                               
   END IF
   IF (CDCONV /= 'NONE') THEN
     ALLOCATE( ZSAVE_PRCONV(SIZE(XPRCONV,1),SIZE(XPRCONV,2),IMODSON))
     ALLOCATE( ZSAVE_PRSCONV(SIZE(XPRSCONV,1),SIZE(XPRSCONV,2),IMODSON))
   ELSE
     ALLOCATE( ZSAVE_PRCONV(0,0,0))                                
     ALLOCATE( ZSAVE_PRSCONV(0,0,0))                                    
   END IF
   IF (CRAD /= 'NONE') THEN
     ALLOCATE( ZSAVE_DIRFLASWD(SIZE(XDIRFLASWD,1),SIZE(XDIRFLASWD,2),SIZE(XDIRFLASWD,3),IMODSON))
     ALLOCATE( ZSAVE_SCAFLASWD(SIZE(XSCAFLASWD,1),SIZE(XSCAFLASWD,2),SIZE(XSCAFLASWD,3),IMODSON))
     ALLOCATE( ZSAVE_DIRSRFSWD(SIZE(XDIRSRFSWD,1),SIZE(XDIRSRFSWD,2),SIZE(XDIRSRFSWD,3),IMODSON))
   ELSE
     ALLOCATE( ZSAVE_DIRFLASWD(0,0,0,0))
     ALLOCATE( ZSAVE_SCAFLASWD(0,0,0,0))
     ALLOCATE( ZSAVE_DIRSRFSWD(0,0,0,0)) 
   END IF
 ENDIF
!
IKIDM=0
DO JKID = IMI+1,NMODEL  ! min value of the possible kids
 IF (IMI == NDAD(JKID) .AND. XWAY(JKID) == 2. .AND. CPROGRAM=='MESONH' &
  .AND. (CCONF == 'RESTA' .OR. (CCONF == 'START' .AND. KTCOUNT /= 1))) THEN
! BUG if number of the son does not follow the number of the dad
! IKIDM = JKID-IMI
  IKIDM = IKIDM + 1
   IF (LUSERC .AND. (                                               &
       (LSEDIC .AND. CCLOUD(1:3) == 'ICE')                     .OR. &
       (LSEDC  .AND. (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO')) .OR. &
       (MSEDC  .AND. CCLOUD=='LIMA')                                &
      )) THEN
     ZSAVE_INPRC(:,:,IKIDM) = XINPRC(:,:)
   END IF
   IF (LUSERR) THEN
     ZSAVE_INPRR(:,:,IKIDM) = XINPRR(:,:)
   END IF
   IF (LUSERS) THEN
     ZSAVE_INPRS(:,:,IKIDM) = XINPRS(:,:)
   END IF
   IF (LUSERG) THEN
     ZSAVE_INPRG(:,:,IKIDM) = XINPRG(:,:)
   END IF
   IF (LUSERH) THEN
     ZSAVE_INPRH(:,:,IKIDM) = XINPRH(:,:)
   END IF
   IF (CDCONV /= 'NONE') THEN
     ZSAVE_PRCONV(:,:,IKIDM) = XPRCONV(:,:)
     ZSAVE_PRSCONV(:,:,IKIDM) = XPRSCONV(:,:)
   END IF
   IF (CRAD /= 'NONE') THEN
     ZSAVE_DIRFLASWD(:,:,:,IKIDM) = XDIRFLASWD(:,:,:)
     ZSAVE_SCAFLASWD(:,:,:,IKIDM) = XSCAFLASWD(:,:,:)
     ZSAVE_DIRSRFSWD(:,:,:,IKIDM) = XDIRSRFSWD(:,:,:)
   END IF
 ENDIF
END DO
!
!-----------------------------------------------------------------------------
!
!*        1.    RADIATION SCHEME
!               ----------------
!
!
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
CALL SECOND_MNH2(ZTIME1)
!
!
!*        1.1   Tests to control how the radiation package should be called (at the current timestep)
!               -----------------------------------------------------------
!
!
GRAD = .FALSE.
OCLOUD_ONLY = .FALSE.
!
IF (CRAD /='NONE') THEN
!
!  test to see if the partial radiations for cloudy must be called
!
  IF (CRAD =='ECMW' .OR. CRAD =='ECRA') THEN
    CALL DATETIME_DISTANCE(TDTRAD_CLONLY,TDTCUR,ZTEMP_DIST)
    IF( MOD(NINT(ZTEMP_DIST/XTSTEP),NINT(XDTRAD_CLONLY/XTSTEP))==0 ) THEN
      TDTRAD_CLONLY = TDTCUR
      GRAD = .TRUE.
      OCLOUD_ONLY = .TRUE.
    END IF
  END IF
!   
! test to see if the full radiations must be called
!   
  CALL DATETIME_DISTANCE(TDTCUR,TDTRAD_FULL,ZTEMP_DIST)
  IF( MOD(NINT(ZTEMP_DIST/XTSTEP),NINT(XDTRAD/XTSTEP))==0 ) THEN
    TDTRAD_FULL = TDTCUR
    GRAD = .TRUE.
    OCLOUD_ONLY = .FALSE.
  END IF
!
! tests to see if any cloud exists
!   
  IF (CRAD =='ECMW' .OR. CRAD =='ECRA') THEN
    IF (GRAD .AND. NRR.LE.3 ) THEN 
      IF( MAX(MAXVAL(XCLDFR(:,:,:)),MAXVAL(XICEFR(:,:,:))).LE. 1.E-10 .AND. OCLOUD_ONLY ) THEN
          GRAD = .FALSE.                ! only the cloudy verticals would be 
                                        ! refreshed but there is no clouds 
      END IF
    END IF
!
    IF (GRAD .AND. NRR.GE.4 ) THEN 
      IF( CCLOUD(1:3)=='ICE' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN(4) .AND. OCLOUD_ONLY ) THEN
            GRAD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF
      IF( CCLOUD=='C3R5' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN_C1R3(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN_C1R3(4) .AND. OCLOUD_ONLY ) THEN
            GRAD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF
      IF( CCLOUD=='LIMA' )THEN
        IF( MAXVAL(XRT(:,:,:,2)).LE.XRTMIN_LIMA(2) .AND.             &
            MAXVAL(XRT(:,:,:,4)).LE.XRTMIN_LIMA(4) .AND. OCLOUD_ONLY ) THEN
            GRAD = .FALSE.            ! only the cloudy verticals would be 
                                      ! refreshed but there is no cloudwater and ice
        END IF
      END IF      
    END IF
  END IF
!
END IF
!
! global parallel mask for 'GRAD'
ZRAD_GLOB_ll = 0.0
IF (GRAD) ZRAD_GLOB_ll = 1.0
CALL REDUCESUM_ll(ZRAD_GLOB_ll,INFO_ll)
if (ZRAD_GLOB_ll .NE. 0.0 ) GRAD = .TRUE.
!
!
IF( GRAD ) THEN                                 
  ALLOCATE(ZCOSZEN(IIU,IJU))
  ALLOCATE(ZSINZEN(IIU,IJU))
  ALLOCATE(ZAZIMSOL(IIU,IJU))
!
!
!*        1.2.  Astronomical computations
!               -------------------------
!
! Ajout PP
IF (.NOT. OCLOUD_ONLY .AND. KTCOUNT /= 1)  THEN 
 IF (LAERO_FT) THEN 
  CALL AEROZON (XPABST,XTHT,XTSRAD,XLAT,XLON,TDTCUR,TDTEXP,   &
         NDLON,NFLEV,CAER,NAER,NSTATM,                             &
         XSINDEL,XCOSDEL,XTSIDER,XCORSOL,                          &
         XSTATM,XOZON, XAER)
  XAER_CLIM = XAER
 END IF
END IF
!
CALL SUNPOS_n   ( XZENITH, ZCOSZEN, ZSINZEN, ZAZIMSOL )
!
!*        1.3   Call to radiation scheme
!               ------------------------
!
  SELECT CASE ( CRAD )
!
!*        1.3.1 TOP of Atmposphere radiation
!               ----------------------------
    CASE('TOPA')
!
      XFLALWD   (:,:)   = 300.
      DO JSWB=1,NSWB_MNH
        XDIRFLASWD(:,:,JSWB) = CST%XI0 * MAX(COS(XZENITH(:,:)),0.)/REAL(NSWB_MNH)
        XSCAFLASWD(:,:,JSWB) = 0.
      END DO
      XDTHRAD(:,:,:) = 0.
     
!
!*        1.3.1 FIXEd radiative surface fluxes
!               ------------------------------
!
    CASE('FIXE')
      ZTIME = MOD(TDTCUR%xtime +XLON0*240., CST%XDAY)
      IHOUR = INT( ZTIME/3600. )
      IF (IHOUR < 0) IHOUR=IHOUR + 24
      ZDT = ZTIME/3600. - REAL(IHOUR)
      XDIRFLASWD(:,:,:) =(( ZRG_HOUR(IHOUR+1)-ZRG_HOUR(IHOUR) )*ZDT + ZRG_HOUR(IHOUR)) / REAL(NSWB_MNH)
      XFLALWD   (:,:)   = (ZRAT_HOUR(IHOUR+1)-ZRAT_HOUR(IHOUR))*ZDT + ZRAT_HOUR(IHOUR)
      DO JSWB=1,NSWB_MNH
        WHERE(ZCOSZEN(:,:)<0.) XDIRFLASWD(:,:,JSWB) = 0.
      END DO

      XSCAFLASWD(:,:,:) = XDIRFLASWD(:,:,:) * 0.2
      XDIRFLASWD(:,:,:) = XDIRFLASWD(:,:,:) * 0.8
      XDTHRAD(:,:,:) = 0.
      !
!
!*        1.3.2 ECMWF or ECRAD radiative surface and atmospheric fluxes
!               ----------------------------------------------
!
    CASE('ECMW' , 'ECRA')
      IF (LLES_MEAN) OCLOUD_ONLY=.FALSE.
      XRADEFF(:,:,:)=0.0
      XSWU(:,:,:)=0.0
      XSWD(:,:,:)=0.0
      XLWU(:,:,:)=0.0
      XLWD(:,:,:)=0.0
      XDTHRADSW(:,:,:)=0.0
      XDTHRADLW(:,:,:)=0.0
      CALL RADIATIONS( TPFILE,                                                                   &
                       LCLEAR_SKY, OCLOUD_ONLY, NCLEARCOL_TM1, CEFRADL, CEFRADI, COPWSW, COPISW, &
                       COPWLW, COPILW, XFUDG,                                                    &
                       NDLON, NFLEV, NRAD_DIAG, NFLUX, NRAD, NAER, NSWB_OLD, NSWB_MNH, NLWB_MNH, &
                       NSTATM, NRAD_COLNBR, ZCOSZEN, XSEA, XCORSOL,                              &
                       XDIR_ALB, XSCA_ALB, XEMIS, MAX(XCLDFR,XICEFR), XCCO2, XTSRAD, XSTATM, XTHT, XRT,      &
                       XPABST, XOZON, XAER,XDST_WL, XAER_CLIM, XSVT,                             &
                       XDTHRAD, XFLALWD, XDIRFLASWD, XSCAFLASWD, XRHODREF, XZZ ,                 &
                       XRADEFF, XSWU, XSWD, XLWU, XLWD, XDTHRADSW, XDTHRADLW                     )
!

      WRITE(UNIT=ILUOUT,FMT='("  RADIATIONS called for KTCOUNT=",I6,       &
         &  "with the CLOUD_ONLY option set ",L2)')   KTCOUNT,OCLOUD_ONLY
!
      !
      WHERE (XDIRFLASWD.LT.0.0)
         XDIRFLASWD=0.0
      ENDWHERE
      !
      WHERE (XDIRFLASWD.GT.1500.0)
         XDIRFLASWD=1500.0
      ENDWHERE
      !
      WHERE (XSCAFLASWD.LT.0.0) 
         XSCAFLASWD=0.0
      ENDWHERE
      !
      WHERE (XSCAFLASWD.GT.1500.0) 
         XSCAFLASWD=1500.0
      ENDWHERE
      !
      WHERE( XDIRFLASWD(:,:,1) + XSCAFLASWD(:,:,1) >0. )
        XALBUV(:,:) = (  XDIR_ALB(:,:,1) * XDIRFLASWD(:,:,1)   &
                       + XSCA_ALB(:,:,1) * XSCAFLASWD(:,:,1) ) &
                    / (XDIRFLASWD(:,:,1) + XSCAFLASWD(:,:,1) )
      ELSEWHERE
        XALBUV(:,:) = XDIR_ALB(:,:,1)
      END WHERE
!
  END SELECT
!
  CALL SECOND_MNH2(ZTIME2)
!
  PRAD = PRAD + ZTIME2 - ZTIME1
!
  ZTIME1 = ZTIME2
!
  CALL SURF_RAD_MODIF (XMAP, XDXHAT, XDYHAT, XXHAT, XYHAT, &
                  ZCOSZEN, ZSINZEN, ZAZIMSOL, XZS, XZS_XY,   &
                  XDIRFLASWD, XDIRSRFSWD                     )
!
!* Azimuthal angle to be sent later to surface processes
!  Defined in radian, clockwise, from North
!
  XAZIM = ZAZIMSOL
!
  CALL SECOND_MNH2(ZTIME2)
!
  PSHADOWS = PSHADOWS + ZTIME2 - ZTIME1
!
  ZTIME1 = ZTIME2
!
  DEALLOCATE(ZCOSZEN)
  DEALLOCATE(ZSINZEN)
  DEALLOCATE(ZAZIMSOL)
!
END IF
!
!
!*        1.4   control prints
!               --------------
!
!*        1.5   Radiative tendency integration
!               ------------------------------
!
IF (CRAD /='NONE') THEN
  if ( TBUCONF%LBUDGET_th ) call Budget_store_init( TBUDGETS(NBUDGET_TH), 'RAD', xrths(:, :, :) )
  XRTHS(:,:,:) = XRTHS(:,:,:) + XRHODJ(:,:,:)*XDTHRAD(:,:,:)
  if ( TBUCONF%LBUDGET_th ) call Budget_store_end ( TBUDGETS(NBUDGET_TH), 'RAD', xrths(:, :, :) )
END IF
!
!
!*        1.6   Ocean case:
! Sfc turbulent fluxes & Radiative tendency due to SW penetrating ocean
! 
IF (LCOUPLES) THEN
ZSFU(:,:)= XSSUFL_C(:,:,1)
ZSFV(:,:)= XSSVFL_C(:,:,1)
ZSFTH(:,:)= XSSTFL_C(:,:,1)
ZSFRV(:,:)=XSSRFL_C(:,:,1)
ELSE 
IF (LOCEAN) THEN
!
  ALLOCATE( ZIZOCE(IKU)); ZIZOCE(:)=0. 
  ALLOCATE( ZPROSOL1(IKU))
  ALLOCATE( ZPROSOL2(IKU))
  ALLOCATE(XSSOLA(IIU,IJU))
  ! Time interpolation
  JSW     = INT(TDTCUR%xtime/REAL(NINFRT))
  ZSWA    = TDTCUR%xtime/REAL(NINFRT)-REAL(JSW)
  ZSFRV = 0.
  ZSFTH  = (XSSTFL_T(JSW+1)*(1.-ZSWA)+XSSTFL_T(JSW+2)*ZSWA) 
  ZSFU = (XSSUFL_T(JSW+1)*(1.-ZSWA)+XSSUFL_T(JSW+2)*ZSWA)
  ZSFV = (XSSVFL_T(JSW+1)*(1.-ZSWA)+XSSVFL_T(JSW+2)*ZSWA)
!
  ZIZOCE(IKU)   = XSSOLA_T(JSW+1)*(1.-ZSWA)+XSSOLA_T(JSW+2)*ZSWA
  ZPROSOL1(IKU) = CST%XROC*ZIZOCE(IKU)
  ZPROSOL2(IKU) = (1.-CST%XROC)*ZIZOCE(IKU)
  if ( TBUCONF%LBUDGET_th ) call Budget_store_init( TBUDGETS(NBUDGET_TH), 'OCEAN', xrths(:, :, :) ) 
  DO JKM=IKU-1,2,-1
    ZPROSOL1(JKM) = ZPROSOL1(JKM+1)* exp(-XDZZ(2,2,JKM)/CST%XD1)
    ZPROSOL2(JKM) = ZPROSOL2(JKM+1)* exp(-XDZZ(2,2,JKM)/CST%XD2)
    ZIZOCE(JKM)   = (ZPROSOL1(JKM+1)-ZPROSOL1(JKM) + ZPROSOL2(JKM+1)-ZPROSOL2(JKM))/XDZZ(2,2,JKM)
    ! Adding to temperature tendency, the solar radiation penetrating in ocean
    XRTHS(:,:,JKM) = XRTHS(:,:,JKM) + XRHODJ(:,:,JKM)*ZIZOCE(JKM)
  END DO
  if ( TBUCONF%LBUDGET_th ) call Budget_store_end ( TBUDGETS(NBUDGET_TH), 'OCEAN', xrths(:, :, :) )
  DEALLOCATE (XSSOLA)
  DEALLOCATE( ZIZOCE) 
  DEALLOCATE (ZPROSOL1)
  DEALLOCATE (ZPROSOL2)
END IF! LOCEAN NO LCOUPLES
END IF!NO LCOUPLES
!
!
CALL SECOND_MNH2(ZTIME2)
!
PRAD = PRAD + ZTIME2 - ZTIME1 &
     - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
PTIME_BU = PTIME_BU + XTIME_LES_BU_PROCESS + XTIME_BU_PROCESS
!
!
!-----------------------------------------------------------------------------
!
!*        2.    DEEP CONVECTION SCHEME
!               ----------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
CALL SECOND_MNH2(ZTIME1)
!
IF( CDCONV == 'KAFR' .OR. CSCONV == 'KAFR' ) THEN

  if (  TBUCONF%LBUDGET_th ) call Budget_store_init( TBUDGETS(NBUDGET_TH), 'DCONV', xrths(:, :, :)    )
  if (  TBUCONF%LBUDGET_rv ) call Budget_store_init( TBUDGETS(NBUDGET_RV), 'DCONV', xrrs (:, :, :, 1) )
  if (  TBUCONF%LBUDGET_rc ) call Budget_store_init( TBUDGETS(NBUDGET_RC), 'DCONV', xrrs (:, :, :, 2) )
  if (  TBUCONF%LBUDGET_ri ) call Budget_store_init( TBUDGETS(NBUDGET_RI), 'DCONV', xrrs (:, :, :, 4) )
  if (  TBUCONF%LBUDGET_sv .and. lchtrans ) then
    do jsv = 1, size( xrsvs, 4 )
      call Budget_store_init( TBUDGETS(NBUDGET_SV1 - 1 + jsv), 'DCONV', xrsvs (:, :, :, jsv) )
    end do
  end if
!
! test to see if the deep convection scheme should be called
!
  GDCONV = .FALSE.
!
  CALL DATETIME_DISTANCE(TDTDCONV,TDTCUR,ZTEMP_DIST)
  IF( MOD(NINT(ZTEMP_DIST/XTSTEP),NINT(XDTCONV/XTSTEP))==0 ) THEN
    TDTDCONV = TDTCUR
    GDCONV   = .TRUE.
  END IF
!
  IF( GDCONV ) THEN
    IF (CDCONV == 'KAFR' .OR. CSCONV == 'KAFR' ) THEN
        ALLOCATE( ZRC(IIU,IJU,IKU) )
        ALLOCATE( ZRI(IIU,IJU,IKU) )
        ALLOCATE( ZWT(IIU,IJU,IKU) )
        ALLOCATE( ZDXDY(IIU,IJU) )
        ! Compute grid area
        ZDXDY(:,:) = SPREAD(XDXHAT(1:IIU),2,IJU) * SPREAD(XDYHAT(1:IJU),1,IIU)
        !
        IF( LUSERC .AND. LUSERI ) THEN
          ZRC(:,:,:) = XRT(:,:,:,2)
          ZRI(:,:,:) = XRT(:,:,:,4)
        ELSE IF( LUSERC .AND. (.NOT. LUSERI) ) THEN
          ZRC(:,:,:) = XRT(:,:,:,2)
          ZRI(:,:,:) = 0.0
        ELSE
          ZRC(:,:,:) = 0.0
          ZRI(:,:,:) = 0.0
        END IF
        WRITE(UNIT=ILUOUT,FMT='("  CONVECTION called for KTCOUNT=",I6)')  &
                                              KTCOUNT
        IF ( LFORCING .AND. L1D ) THEN
          ZWT(:,:,:) = XWTFRC(:,:,:)
        ELSE
          ZWT(:,:,:) = XWT(:,:,:)
        ENDIF
        IF (LDUST) CALL DUST_FILTER(XSVT(:,:,:,NSV_DSTBEG:NSV_DSTEND), XRHODREF(:,:,:))
        IF (LSALT) CALL SALT_FILTER(XSVT(:,:,:,NSV_SLTBEG:NSV_SLTEND), XRHODREF(:,:,:))
        IF (LCH_CONV_LINOX) THEN
          CALL CONVECTION( XDTCONV, CDCONV, CSCONV, LREFRESH_ALL, LDOWN, NICE, &
                         LSETTADJ, XTADJD, XTADJS, LDIAGCONV, NENSM,           &
                         XPABST, XZZ, ZDXDY,                                   &
                         XTHT, XRT(:,:,:,1), ZRC, ZRI, XUT, XVT,               &
                         ZWT,XTKET(:,:,IKB),                                   &
                         NCOUNTCONV, XDTHCONV, XDRVCONV, XDRCCONV, XDRICONV,   &
                         XPRCONV, XPRSCONV,                                    &
                         XUMFCONV,XDMFCONV,XMFCONV,XPRLFLXCONV,XPRSFLXCONV,    &
                         XCAPE, NCLTOPCONV, NCLBASCONV,                        &
                         LCHTRANS, XSVT, XDSVCONV,                             &
                         LUSECHEM, LCH_CONV_SCAV, LCH_CONV_LINOX,              &
                         LDUST, LSALT,                                         &
                         XRHODREF, XIC_RATE, XCG_RATE                          )
        ELSE
          CALL CONVECTION( XDTCONV, CDCONV, CSCONV, LREFRESH_ALL, LDOWN, NICE, &
                         LSETTADJ, XTADJD, XTADJS, LDIAGCONV, NENSM,           &
                         XPABST, XZZ, ZDXDY,                                   &
                         XTHT, XRT(:,:,:,1), ZRC, ZRI, XUT, XVT,               &
                         ZWT,XTKET(:,:,IKB),                                   &
                         NCOUNTCONV, XDTHCONV, XDRVCONV, XDRCCONV, XDRICONV,   &
                         XPRCONV, XPRSCONV,                                    &
                         XUMFCONV,XDMFCONV,XMFCONV,XPRLFLXCONV,XPRSFLXCONV,    &
                         XCAPE, NCLTOPCONV, NCLBASCONV,                        &
                         LCHTRANS, XSVT, XDSVCONV,                             &
                         LUSECHEM, LCH_CONV_SCAV, LCH_CONV_LINOX,              &
                         LDUST, LSALT,                                         &
                         XRHODREF )
        END IF
!
        DEALLOCATE( ZRC )
        DEALLOCATE( ZRI )
        DEALLOCATE( ZWT )
        DEALLOCATE( ZDXDY )
    END IF    
  END IF
!
!  Deep convection tendency integration
!
  XRTHS(:,:,:)  = XRTHS(:,:,:)  + XRHODJ(:,:,:) * XDTHCONV(:,:,:)
  XRRS(:,:,:,1) = XRRS(:,:,:,1) + XRHODJ(:,:,:) * XDRVCONV(:,:,:)
!
!
! Aerosols size distribution
! Compute Rg and sigma before tracers convection tendency (for orilam, dust and sea
! salt)
!

  IF ( LCHTRANS ) THEN  ! update tracers for chemical transport
    IF (LORILAM) ZRSVS(:,:,:,:) = XRSVS(:,:,:,:)    !
    IF ((LDUST)) THEN ! dust convective balance
      ALLOCATE(ZSIGDST(IIU,IJU,IKU,NMODE_DST))
      ALLOCATE(ZRGDST(IIU,IJU,IKU,NMODE_DST))
      ALLOCATE(ZNDST(IIU,IJU,IKU,NMODE_DST))
      ALLOCATE(ZSVDST(IIU,IJU,IKU,NSV_DST))
      !
      DO JSV=1,NMODE_DST
        IMODEIDX = JPDUSTORDER(JSV)
        IF (CRGUNITD=="MASS") THEN
          ZINIRADIUS(JSV) = XINIRADIUS(IMODEIDX) * EXP(-3.*(LOG(XINISIG(IMODEIDX)))**2)
        ELSE
          ZINIRADIUS(JSV) = XINIRADIUS(IMODEIDX)
        END IF
        ZSIGDST(:,:,:,JSV) = XINISIG(IMODEIDX)
        ZRGDST(:,:,:,JSV)  = ZINIRADIUS(JSV)
        ZNDST(:,:,:,JSV)   = XN0MIN(IMODEIDX)
      ENDDO
      !
    IF (CPROGRAM == "MESONH") THEN
      DO JSV=NSV_DSTBEG,NSV_DSTEND
        ZSVDST(:,:,:,JSV-NSV_DSTBEG+1) = XRSVS(:,:,:,JSV) * XTSTEP / XRHODJ(:,:,:) 
      ENDDO
    ELSE
      DO JSV=NSV_DSTBEG,NSV_DSTEND
        ZSVDST(:,:,:,JSV-NSV_DSTBEG+1) = XSVT(:,:,:,JSV)
      ENDDO
    ENDIF
      CALL PPP2DUST(ZSVDST(IIB:IIE,IJB:IJE,IKB:IKE,:), XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),&
              PSIG3D=ZSIGDST(IIB:IIE,IJB:IJE,IKB:IKE,:), PRG3D=ZRGDST(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
              PN3D=ZNDST(IIB:IIE,IJB:IJE,IKB:IKE,:))
    END IF
    !
    IF ((LSALT)) THEN ! sea salt convective balance
      ALLOCATE(ZSIGSLT(IIU,IJU,IKU,NMODE_SLT))
      ALLOCATE(ZRGSLT(IIU,IJU,IKU,NMODE_SLT))
      ALLOCATE(ZNSLT(IIU,IJU,IKU,NMODE_SLT))
      ALLOCATE(ZSVSLT(IIU,IJU,IKU,NSV_SLT))
      !
      DO JSV=1,NMODE_SLT
        IMODEIDX = JPSALTORDER(JSV)
        IF (CRGUNITS=="MASS") THEN
          ZINIRADIUS_SLT(JSV) = XINIRADIUS_SLT(IMODEIDX) * &
                            EXP(-3.*(LOG(XINISIG_SLT(IMODEIDX)))**2)
        ELSE
          ZINIRADIUS_SLT(JSV) = XINIRADIUS_SLT(IMODEIDX)
        END IF
        ZSIGSLT(:,:,:,JSV) = XINISIG_SLT(IMODEIDX)
        ZRGSLT(:,:,:,JSV)  = ZINIRADIUS_SLT(JSV)
        ZNSLT(:,:,:,JSV)   = XN0MIN_SLT(IMODEIDX)
      ENDDO
      !
    IF (CPROGRAM == "MESONH") THEN
      DO JSV=NSV_SLTBEG,NSV_SLTEND
        ZSVSLT(:,:,:,JSV-NSV_SLTBEG+1) = XRSVS(:,:,:,JSV) * XTSTEP / XRHODJ(:,:,:) 
      ENDDO
    ELSE
      DO JSV=NSV_SLTBEG,NSV_SLTEND
        ZSVSLT(:,:,:,JSV-NSV_SLTBEG+1) = XSVT(:,:,:,JSV)
      ENDDO
    END IF
      CALL PPP2SALT(ZSVSLT(IIB:IIE,IJB:IJE,IKB:IKE,:), XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE),&
              PSIG3D=ZSIGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:), PRG3D=ZRGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),   &
              PN3D=ZNSLT(IIB:IIE,IJB:IJE,IKB:IKE,:))
    END IF
    !
!
! Compute convective tendency for all tracers
!
  IF (LCHTRANS) THEN
    DO JSV = 1, SIZE(XRSVS,4)
      XRSVS(:,:,:,JSV) = XRSVS(:,:,:,JSV) + XRHODJ(:,:,:) * XDSVCONV(:,:,:,JSV)
    END DO
    IF (LORILAM) THEN
      DO JSV = NSV_AERBEG,NSV_AEREND
        PWETDEPAER(:,:,:,JSV-NSV_AERBEG+1) = XDSVCONV(:,:,:,JSV) * XRHODJ(:,:,:)
        XRSVS(:,:,:,JSV) = ZRSVS(:,:,:,JSV) 
      END DO
    END IF  
  END IF
!
  IF ((LDUST).AND.(LCHTRANS)) THEN ! dust convective balance
    IF (CPROGRAM == "MESONH") THEN
      DO JSV=NSV_DSTBEG,NSV_DSTEND
          ZSVDST(:,:,:,JSV-NSV_DSTBEG+1) = XRSVS(:,:,:,JSV) * XTSTEP / XRHODJ(:,:,:) 
      ENDDO
    ELSE
      DO JSV=NSV_DSTBEG,NSV_DSTEND
        ZSVDST(:,:,:,JSV-NSV_DSTBEG+1) = XSVT(:,:,:,JSV)
      ENDDO
    ENDIF
    CALL DUST2PPP(ZSVDST(IIB:IIE,IJB:IJE,IKB:IKE,:), &
                    XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE), ZSIGDST(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                    ZRGDST(IIB:IIE,IJB:IJE,IKB:IKE,:))
    DO JSV=NSV_DSTBEG,NSV_DSTEND
      XRSVS(:,:,:,JSV) =  ZSVDST(:,:,:,JSV-NSV_DSTBEG+1) * XRHODJ(:,:,:) / XTSTEP
    ENDDO
    !
    DEALLOCATE(ZSVDST)
    DEALLOCATE(ZNDST)
    DEALLOCATE(ZRGDST)
    DEALLOCATE(ZSIGDST)
  END IF
    !
  IF ((LSALT).AND.(LCHTRANS)) THEN ! sea salt convective balance
    IF (CPROGRAM == "MESONH") THEN
      DO JSV=NSV_SLTBEG,NSV_SLTEND
        ZSVSLT(:,:,:,JSV-NSV_SLTBEG+1) = XRSVS(:,:,:,JSV) * XTSTEP / XRHODJ(:,:,:) 
      ENDDO
    ELSE
      DO JSV=NSV_SLTBEG,NSV_SLTEND
        ZSVSLT(:,:,:,JSV-NSV_SLTBEG+1) = XSVT(:,:,:,JSV)
      ENDDO
    END IF
    CALL SALT2PPP(ZSVSLT(IIB:IIE,IJB:IJE,IKB:IKE,:), &
                  XRHODREF(IIB:IIE,IJB:IJE,IKB:IKE), ZSIGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:),&
                  ZRGSLT(IIB:IIE,IJB:IJE,IKB:IKE,:))
    DO JSV=NSV_SLTBEG,NSV_SLTEND
      XRSVS(:,:,:,JSV) =  ZSVSLT(:,:,:,JSV-NSV_SLTBEG+1) * XRHODJ(:,:,:) / XTSTEP
    ENDDO
    !
    DEALLOCATE(ZSVSLT)
    DEALLOCATE(ZNSLT)
    DEALLOCATE(ZRGSLT)
    DEALLOCATE(ZSIGSLT)
  END IF
  !
END IF
!
  IF( LUSERC .AND. LUSERI ) THEN
    XRRS(:,:,:,2) = XRRS(:,:,:,2) + XRHODJ(:,:,:) * XDRCCONV(:,:,:)
    XRRS(:,:,:,4) = XRRS(:,:,:,4) + XRHODJ(:,:,:) * XDRICONV(:,:,:)
!
  ELSE IF ( LUSERC .AND. (.NOT. LUSERI) ) THEN
!
!  If only cloud water but no cloud ice is used, the convective tendency
!     for cloud ice is added to the tendency for cloud water
!
      XRRS(:,:,:,2) = XRRS(:,:,:,2) + XRHODJ(:,:,:) * (XDRCCONV(:,:,:) + &
                                                       XDRICONV(:,:,:)   )
!     and cloud ice is melted
!
      XRTHS(:,:,:) = XRTHS(:,:,:) - XRHODJ(:,:,:) *                      &
         ( XP00/XPABST(:,:,:) )**(XRD/XCPD) * CST%XLMTT / XCPD * XDRICONV(:,:,:)
!
  ELSE IF ( (.NOT. LUSERC) .AND. (.NOT. LUSERI) ) THEN
!
!  If no cloud water and no cloud ice are used the convective tendencies for these
!     variables are added to the water vapor tendency
!
      XRRS(:,:,:,1) = XRRS(:,:,:,1) + XRHODJ(:,:,:) * (XDRCCONV(:,:,:) + &
                                                       XDRICONV(:,:,:)   )
!     and all cloud condensate is evaporated
!
      XRTHS(:,:,:) = XRTHS(:,:,:) - XRHODJ(:,:,:) / XCPD * (              &
                     CST%XLVTT * XDRCCONV(:,:,:) + CST%XLSTT * XDRICONV(:,:,:) ) *&
                    ( XP00 / XPABST(:,:,:) ) ** ( XRD / XCPD )
  END IF

  if (  TBUCONF%LBUDGET_th ) call Budget_store_end( TBUDGETS(NBUDGET_TH), 'DCONV', xrths(:, :, :)    )
  if (  TBUCONF%LBUDGET_rv ) call Budget_store_end( TBUDGETS(NBUDGET_RV), 'DCONV', xrrs (:, :, :, 1) )
  if (  TBUCONF%LBUDGET_rc ) call Budget_store_end( TBUDGETS(NBUDGET_RC), 'DCONV', xrrs (:, :, :, 2) )
  if (  TBUCONF%LBUDGET_ri ) call Budget_store_end( TBUDGETS(NBUDGET_RI), 'DCONV', xrrs (:, :, :, 4) )
  if (  TBUCONF%LBUDGET_sv .and. lchtrans ) then
    do jsv = 1, size( xrsvs, 4 )
      call Budget_store_end( TBUDGETS(NBUDGET_SV1 - 1 + jsv), 'DCONV', xrsvs (:, :, :, jsv) )
    end do
  end if
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
PKAFR = PKAFR + ZTIME2 - ZTIME1 &
       - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
PTIME_BU = PTIME_BU + XTIME_LES_BU_PROCESS + XTIME_BU_PROCESS
!
!-----------------------------------------------------------------------------
!
!*        3.    TURBULENT SURFACE FLUXES
!               ------------------------
!
ZTIME1 = ZTIME2
!
IF (CSURF=='EXTE') THEN
  CALL GOTO_SURFEX(IMI)
!
  IF( LTRANS ) THEN
    XUT(:,:,1+JPVEXT) = XUT(:,:,1+JPVEXT) + XUTRANS
    XVT(:,:,1+JPVEXT) = XVT(:,:,1+JPVEXT) + XVTRANS
  END IF
  !
  ALLOCATE(ZDIR_ALB(IIU,IJU,NSWB_MNH))
  ALLOCATE(ZSCA_ALB(IIU,IJU,NSWB_MNH))
  ALLOCATE(ZEMIS  (IIU,IJU,NLWB_MNH))
  ALLOCATE(ZTSRAD (IIU,IJU))
  !  
  IKIDM=0
  DO JKID = IMI+1,NMODEL  ! min value of the possible kids
    IF (IMI == NDAD(JKID) .AND. XWAY(JKID) == 2. .AND. &
     CPROGRAM=='MESONH' .AND. &
     (CCONF == 'RESTA' .OR. (CCONF == 'START' .AND. KTCOUNT /= 1))) THEN
    !  where kids exist, use the two-way output fields (i.e. OMASKkids true)
    !  rather than the farther calculations in radiation and convection schemes
! BUG if number of the son does not follow the number of the dad
!    IKIDM = JKID-IMI
      IKIDM = IKIDM + 1
     IF (LUSERC .AND. (                                               &
         (LSEDIC .AND. CCLOUD(1:3) == 'ICE')                     .OR. &
         (LSEDC  .AND. (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO')) .OR. &
         (MSEDC  .AND. CCLOUD=='LIMA')                                &
        )) THEN
         WHERE (OMASKkids(:,:) )
            XINPRC(:,:) = ZSAVE_INPRC(:,:,IKIDM)
          ENDWHERE
      END IF
      IF (LUSERR) THEN
        WHERE (OMASKkids(:,:) )
          XINPRR(:,:) = ZSAVE_INPRR(:,:,IKIDM)
        ENDWHERE
      END IF
      IF (LUSERS) THEN
        WHERE (OMASKkids(:,:) )
          XINPRS(:,:) = ZSAVE_INPRS(:,:,IKIDM)
       ENDWHERE
      END IF
      IF (LUSERG) THEN
        WHERE (OMASKkids(:,:) )
          XINPRG(:,:) = ZSAVE_INPRG(:,:,IKIDM)
        ENDWHERE
      END IF
      IF (LUSERH) THEN
        WHERE (OMASKkids(:,:) )
          XINPRH(:,:) = ZSAVE_INPRH(:,:,IKIDM)
        ENDWHERE
      END IF
      IF (CDCONV /= 'NONE') THEN
        WHERE (OMASKkids(:,:) )
          XPRCONV(:,:) = ZSAVE_PRCONV(:,:,IKIDM)
          XPRSCONV(:,:) = ZSAVE_PRSCONV(:,:,IKIDM)
        ENDWHERE
      END IF
      IF (CRAD /= 'NONE') THEN
        DO JSWB=1,NSWB_MNH
          WHERE (OMASKkids(:,:) ) 
            XDIRFLASWD(:,:,JSWB) = ZSAVE_DIRFLASWD(:,:,JSWB,IKIDM)
            XSCAFLASWD(:,:,JSWB) = ZSAVE_SCAFLASWD(:,:,JSWB,IKIDM)
            XDIRSRFSWD(:,:,JSWB) = ZSAVE_DIRSRFSWD(:,:,JSWB,IKIDM)
          ENDWHERE
        ENDDO
      END IF
    ENDIF
  END DO
  !
 IF (IMODSON /= 0 ) THEN
    DEALLOCATE( ZSAVE_INPRR,ZSAVE_INPRS,ZSAVE_INPRG,ZSAVE_INPRH)
    DEALLOCATE( ZSAVE_INPRC,ZSAVE_PRCONV,ZSAVE_PRSCONV)
    DEALLOCATE( ZSAVE_DIRFLASWD,ZSAVE_SCAFLASWD,ZSAVE_DIRSRFSWD)
 END IF
   CALL GROUND_PARAM_n(YLDIMPHYEX,ZSFTH, ZSFTH_WALL, ZSFTH_ROOF, ZCD_ROOF, ZSFRV, ZSFRV_WALL, ZSFRV_ROOF, &
                      ZSFSV, ZSFCO2, ZSFU, ZSFV, ZDIR_ALB, ZSCA_ALB, ZEMIS, ZTSRAD, KTCOUNT, TPFILE  )
  !
  IF (LIBM) THEN
    WHERE(XIBM_LS(:,:,IKB,1).GT.-XIBM_EPSI)
      ZSFTH(:,:)=0.
      ZSFRV(:,:)=0. 
      ZSFU (:,:)=0. 
      ZSFV (:,:)=0.
    ENDWHERE
    IF (NSV>0) THEN
      DO JSV = 1 , NSV
         WHERE(XIBM_LS(:,:,IKB,1).GT.-XIBM_EPSI) ZSFSV(:,:,JSV)=0.
      ENDDO
    ENDIF 
  ENDIF
  !
  IF (SIZE(XEMIS)>0) THEN
    XDIR_ALB = ZDIR_ALB
    XSCA_ALB = ZSCA_ALB
    XEMIS    = ZEMIS
    XTSRAD   = ZTSRAD
  END IF
  !
  DEALLOCATE(ZDIR_ALB)
  DEALLOCATE(ZSCA_ALB)
  DEALLOCATE(ZEMIS   )
  DEALLOCATE(ZTSRAD  )
  !
  !
  IF( LTRANS ) THEN
    XUT(:,:,1+JPVEXT) = XUT(:,:,1+JPVEXT) - XUTRANS
    XVT(:,:,1+JPVEXT) = XVT(:,:,1+JPVEXT) - XVTRANS
  END IF
!
ELSE ! case no SURFEX (CSURF logical)
  ZSFSV    = 0.
  ZSFCO2   = 0.
  ZSFTH_WALL = 0.
  ZSFTH_ROOF = 0.
  ZCD_ROOF   = 0.
  ZSFRV_WALL = 0.
  ZSFRV_ROOF = 0.
  IF (.NOT.LOCEAN) THEN
    ZSFTH    = 0.
    ZSFRV    = 0.
    ZSFSV    = 0.
    ZSFCO2   = 0.
    ZSFU     = 0.
    ZSFV     = 0.
  END IF
END IF !CSURF
!
CALL SECOND_MNH2(ZTIME2)
!
PGROUND = PGROUND + ZTIME2 - ZTIME1
!
!-----------------------------------------------------------------------------
!
!*        3.1    EDDY FLUXES PARAMETRIZATION
!               ------------------
!
IF (IMI==1) THEN  ! On calcule les flus turb. comme preconise par PP

   ! Heat eddy fluxes
   IF ( LTH_FLX ) CALL EDDY_FLUX_n(IMI,KTCOUNT,XVT,XTHT,XRHODJ,XRTHS,XVTH_FLUX_M,XWTH_FLUX_M)
   !
   ! Momentum eddy fluxes
   IF ( LUV_FLX ) CALL EDDYUV_FLUX_n(IMI,KTCOUNT,XVT,XTHT,XRHODJ,XRHODREF,XPABSM,XRVS,XVU_FLUX_M)

ELSE
   ! TEST pour maille infèrieure à 20km ? 
   !      car pb d'instabilités ?
   !      Pour le modèle fils, on spawne les flux du modèle père
   ! Heat eddy fluxes
   IF ( LTH_FLX ) CALL EDDY_FLUX_ONE_WAY_n (IMI,KTCOUNT,NDXRATIO_ALL(IMI),NDYRATIO_ALL(IMI),CLBCX,CLBCY)
   !
   ! Momentum eddy fluxes
   IF ( LUV_FLX ) CALL EDDYUV_FLUX_ONE_WAY_n (IMI,KTCOUNT,NDXRATIO_ALL(IMI),NDYRATIO_ALL(IMI),CLBCX,CLBCY)
   !
END IF
!-----------------------------------------------------------------------------
!
!*        4.    PASSIVE POLLUTANTS
!               ------------------
!
ZTIME1 = ZTIME2
!
IF (LPASPOL) CALL PASPOL(XTSTEP, ZSFSV, ILUOUT, NVERB, TPFILE)
!
!
!*        4b.  PASSIVE POLLUTANTS FOR MASS-FLUX SCHEME DIAGNOSTICS
!              ---------------------------------------------------
!
IF (LCONDSAMP) CALL CONDSAMP(XTSTEP, ZSFSV, ILUOUT, NVERB)
!
CALL SECOND_MNH2(ZTIME2)
!
PTRACER = PTRACER + ZTIME2 - ZTIME1
!-----------------------------------------------------------------------------
!
!*        5a.    Drag force 
!               ----------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZSFTH_WALL, 'PHYS_PARAM_n::ZSFTH_WALL')
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZSFTH_ROOF, 'PHYS_PARAM_n::ZSFTH_ROOF')
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZCD_ROOF, 'PHYS_PARAM_n::ZCD_ROOF')
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZSFRV_WALL, 'PHYS_PARAM_n::ZSFRV_WALL')
CALL ADD2DFIELD_ll(TZFIELDS_ll,ZSFRV_ROOF, 'PHYS_PARAM_n::ZSFRV_ROOF')
!
IF ( CLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
   ZSFTH_WALL(IIB-1,:)=ZSFTH_WALL(IIB,:)
   ZSFTH_ROOF(IIB-1,:)=ZSFTH_ROOF(IIB,:)
   ZCD_ROOF  (IIB-1,:)=ZCD_ROOF(IIB,:)
   ZSFRV_WALL(IIB-1,:)=ZSFRV_WALL(IIB,:)
   ZSFRV_ROOF(IIB-1,:)=ZSFRV_ROOF(IIB,:)
ENDIF
!
IF ( CLBCX(2) /= "CYCL" .AND. LEAST_ll()) THEN
   ZSFTH_WALL(IIE+1,:)=ZSFTH_WALL(IIE,:)
   ZSFTH_ROOF(IIE+1,:)=ZSFTH_ROOF(IIE,:)
   ZCD_ROOF(IIE+1,:)  =ZCD_ROOF(IIE,:)
   ZSFRV_WALL(IIE+1,:)=ZSFRV_WALL(IIE,:)
   ZSFRV_ROOF(IIE+1,:)=ZSFRV_ROOF(IIE,:)
ENDIF
!
IF ( CLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
   ZSFTH_WALL(:,IJB-1)=ZSFTH_WALL(:,IJB)
   ZSFTH_ROOF(:,IJB-1)=ZSFTH_ROOF(:,IJB)
   ZCD_ROOF(:,IJB-1)  =ZCD_ROOF(:,IJB)
   ZSFRV_WALL(:,IJB-1)=ZSFRV_WALL(:,IJB)
   ZSFRV_ROOF(:,IJB-1)=ZSFRV_ROOF(:,IJB)
ENDIF
!
IF ( CLBCY(2) /= "CYCL" .AND. LNORTH_ll()) THEN
   ZSFTH_WALL(:,IJE+1)=ZSFTH_WALL(:,IJE)
   ZSFTH_ROOF(:,IJE+1)=ZSFTH_ROOF(:,IJE)
   ZCD_ROOF(:,IJE+1)=ZCD_ROOF(:,IJE)
   ZSFRV_WALL(:,IJE+1)=ZSFRV_WALL(:,IJE)
   ZSFRV_ROOF(:,IJE+1)=ZSFRV_ROOF(:,IJE)
ENDIF
!
!
IF (LDRAGTREE) CALL DRAG_VEG( XTSTEP, XUT, XVT, XTKET, LDEPOTREE, XVDEPOTREE, &
                              CCLOUD, XPABST, XTHT, XRT, XSVT, XRHODJ, XZZ,   &
                              XRUS, XRVS, XRTKES, XRRS, XRSVS )
!
IF (LDRAGBLDG) CALL DRAG_BLD ( XTSTEP, XUT, XVT, XTKET, XPABST, XTHT, XRT, XSVT, &
                               XRHODJ, XZZ, XRUS, XRVS, XRTKES, XRTHS, XRRS,     &
                               ZSFTH_WALL, ZSFTH_ROOF, ZCD_ROOF, ZSFRV_WALL,     &
                               ZSFRV_ROOF       )
!
CALL SECOND_MNH2(ZTIME2)
!
PDRAG = PDRAG + ZTIME2 - ZTIME1 &
             - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
PTIME_BU = PTIME_BU + XTIME_LES_BU_PROCESS + XTIME_BU_PROCESS
!
!*        5b.   Drag force from wind turbines 
!               -----------------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
IF (LMAIN_EOL .AND. IMI == NMODEL_EOL) THEN
 CALL EOL_MAIN(KTCOUNT,XTSTEP,     &
               XDXX,XDYY,XDZZ,     &
               XRHODJ,             &
               XUT,XVT,XWT,        &
               XRUS, XRVS, XRWS    )
END IF
!
CALL SECOND_MNH2(ZTIME2)
!
PEOL = PEOL + ZTIME2 - ZTIME1 &
             - XTIME_LES_BU_PROCESS - XTIME_BU_PROCESS
!
PTIME_BU = PTIME_BU + XTIME_LES_BU_PROCESS + XTIME_BU_PROCESS
!
!*        
!-----------------------------------------------------------------------------
!
!*        6.    TURBULENCE SCHEME
!               -----------------
!
ZTIME1 = ZTIME2
XTIME_BU_PROCESS = 0.
XTIME_LES_BU_PROCESS = 0.
!
ZSFTH(:,:)  = ZSFTH(:,:) * XDIRCOSZW(:,:)
ZSFRV(:,:)  = ZSFRV(:,:) * XDIRCOSZW(:,:)
DO JSV=1,NSV
  ZSFSV(:,:,JSV)  = ZSFSV(:,:,JSV) * XDIRCOSZW(:,:)
END DO
!
IF (LLES_CALL) CALL SWITCH_SBG_LES_n
!
!
IF ( CTURB == 'TKEL' ) THEN
!

!*        6.1 complete surface flux fields on the border
!
!!$  IF(NHALO == 1) THEN
    CALL ADD2DFIELD_ll( TZFIELDS_ll, ZSFTH, 'PHYS_PARAM_n::ZSFTH' )
    CALL ADD2DFIELD_ll( TZFIELDS_ll, ZSFRV, 'PHYS_PARAM_n::ZSFRV' )
    CALL ADD2DFIELD_ll( TZFIELDS_ll, ZSFU,  'PHYS_PARAM_n::ZSFU' )
    CALL ADD2DFIELD_ll( TZFIELDS_ll, ZSFV,  'PHYS_PARAM_n::ZSFV' )
    IF(NSV >0)THEN
      DO JSV=1,NSV
        write ( ynum, '( I6 ) ' ) jsv
        CALL ADD2DFIELD_ll( TZFIELDS_ll, ZSFSV(:,:,JSV), 'PHYS_PARAM_n::ZSFSV:'//trim( adjustl( ynum ) ) )
      END DO
    END IF
    CALL ADD2DFIELD_ll( TZFIELDS_ll, ZSFCO2, 'PHYS_PARAM_n::ZSFCO2' )
    CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
    CALL CLEANLIST_ll(TZFIELDS_ll)
!!$  END IF
!
  CALL MPPDB_CHECK2D(ZSFU,"phys_param::ZSFU",PRECISION)
  !
  IF ( CLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
    ZSFTH(IIB-1,:)=ZSFTH(IIB,:)
    ZSFRV(IIB-1,:)=ZSFRV(IIB,:)
    ZSFU(IIB-1,:)=ZSFU(IIB,:)
    ZSFV(IIB-1,:)=ZSFV(IIB,:)
    IF (NSV>0)  THEN
      ZSFSV(IIB-1,:,:)=ZSFSV(IIB,:,:)
      WHERE ((ZSFSV(IIB-1,:,:).LT.0.).AND.(XSVT(IIB-1,:,IKB,:).EQ.0.))
          ZSFSV(IIB-1,:,:) = 0.
      END WHERE
    ENDIF
    ZSFCO2(IIB-1,:)=ZSFCO2(IIB,:)
  END IF
  !
  IF ( CLBCX(2) /= "CYCL" .AND. LEAST_ll()) THEN
    ZSFTH(IIE+1,:)=ZSFTH(IIE,:)
    ZSFRV(IIE+1,:)=ZSFRV(IIE,:)
    ZSFU(IIE+1,:)=ZSFU(IIE,:)
    ZSFV(IIE+1,:)=ZSFV(IIE,:)
    IF (NSV>0) THEN
      ZSFSV(IIE+1,:,:)=ZSFSV(IIE,:,:)
      WHERE ((ZSFSV(IIE+1,:,:).LT.0.).AND.(XSVT(IIE+1,:,IKB,:).EQ.0.))
          ZSFSV(IIE+1,:,:) = 0.
      END WHERE
    ENDIF
    ZSFCO2(IIE+1,:)=ZSFCO2(IIE,:)
  END IF
  !
  IF ( CLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
    ZSFTH(:,IJB-1)=ZSFTH(:,IJB)
    ZSFRV(:,IJB-1)=ZSFRV(:,IJB)
    ZSFU(:,IJB-1)=ZSFU(:,IJB)
    ZSFV(:,IJB-1)=ZSFV(:,IJB)
    IF (NSV>0) THEN
      ZSFSV(:,IJB-1,:)=ZSFSV(:,IJB,:)
      WHERE ((ZSFSV(:,IJB-1,:).LT.0.).AND.(XSVT(:,IJB-1,IKB,:).EQ.0.))
          ZSFSV(:,IJB-1,:) = 0.
      END WHERE
    ENDIF
    ZSFCO2(:,IJB-1)=ZSFCO2(:,IJB)
  END IF
  !
  IF ( CLBCY(2) /= "CYCL" .AND. LNORTH_ll()) THEN
    ZSFTH(:,IJE+1)=ZSFTH(:,IJE)
    ZSFRV(:,IJE+1)=ZSFRV(:,IJE)
    ZSFU(:,IJE+1)=ZSFU(:,IJE)
    ZSFV(:,IJE+1)=ZSFV(:,IJE)
    IF (NSV>0) THEN
      ZSFSV(:,IJE+1,:)=ZSFSV(:,IJE,:)
      WHERE ((ZSFSV(:,IJE+1,:).LT.0.).AND.(XSVT(:,IJE+1,IKB,:).EQ.0.))
          ZSFSV(:,IJE+1,:) = 0.
      END WHERE
    ENDIF
    ZSFCO2(:,IJE+1)=ZSFCO2(:,IJE)
  END IF
!
  IF( LTRANS ) THEN
    XUT(:,:,:) = XUT(:,:,:) + XUTRANS
    XVT(:,:,:) = XVT(:,:,:) + XVTRANS
  END IF
!
!
IF ( ALLOCATED( XTHW_FLUX ) ) DEALLOCATE( XTHW_FLUX )
IF ( LFLYER ) THEN
  ALLOCATE( XTHW_FLUX(SIZE( XTHT, 1 ), SIZE( XTHT, 2 ), SIZE( XTHT, 3 )) )
ELSE
  ALLOCATE( XTHW_FLUX(0, 0, 0) )
END IF

IF ( ALLOCATED( XRCW_FLUX ) ) DEALLOCATE( XRCW_FLUX )
IF ( LFLYER ) THEN
  ALLOCATE( XRCW_FLUX(SIZE( XTHT, 1 ), SIZE( XTHT, 2 ), SIZE( XTHT, 3 )) )
ELSE
  ALLOCATE( XRCW_FLUX(0, 0, 0) )
END IF

IF ( ALLOCATED( XSVW_FLUX ) ) DEALLOCATE( XSVW_FLUX )
IF ( LFLYER ) THEN
  ALLOCATE( XSVW_FLUX(SIZE( XSVT, 1 ), SIZE( XSVT, 2 ), SIZE( XSVT, 3 ), SIZE( XSVT, 4 )) )
ELSE
  ALLOCATE( XSVW_FLUX(0, 0, 0, 0) )
END IF
!
GCOMPUTE_SRC=SIZE(XSIGS, 3)/=0
!
ALLOCATE(ZTDIFF(IIU,IJU,IKU))
ALLOCATE(ZTDISS(IIU,IJU,IKU))
!
!! Compute Shape of sfc flux for Oceanic Deep Conv Case
!
IF (LOCEAN .AND. LDEEPOC) THEN
  ALLOCATE(ZDIST(IIU,IJU))
  !*       COMPUTES THE PHYSICAL SUBDOMAIN BOUNDS
  ALLOCATE(ZXHAT_ll(NIMAX_ll+2*JPHEXT),ZYHAT_ll(NJMAX_ll+2*JPHEXT))
  !compute ZXHAT_ll = position in the (0:Lx) domain 1 (Lx=Size of domain1 )
  !compute XXHAT_ll = position in the (L0_subproc,Lx_subproc) domain for the current subproc
  !                                     L0_subproc as referenced in the full domain 1
  CALL GATHERALL_FIELD_ll('XX',XXHAT,ZXHAT_ll,IRESP)
  CALL GATHERALL_FIELD_ll('YY',XYHAT,ZYHAT_ll,IRESP)
  CALL GET_DIM_EXT_ll('B',IIU,IJU)
  DO JJ = IJB,IJE
    DO JI = IIB,IIE
      ZDIST(JI,JJ) = SQRT(                         &
      (( (XXHAT(JI)+XXHAT(JI+1))*0.5 - XCENTX_OC ) / XRADX_OC)**2 + &
      (( (XYHAT(JJ)+XYHAT(JJ+1))*0.5 - XCENTY_OC ) / XRADY_OC)**2   &
                                )
    END DO
  END DO
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IF ( ZDIST(JI,JJ) > 1.) ZSFTH(JI,JJ)=0.
    END DO
  END DO
END IF !END DEEP OCEAN CONV CASE
!
IF(LLEONARD) THEN
  IGRADIENTS=6
  ALLOCATE(ZHGRAD(IIU,IJU,IKU,IGRADIENTS))
  ZHGRAD(:,:,:,1) = GX_W_UW(XWT(:,:,:), XDXX,XDZZ,XDZX,1,IKU,1)
  ZHGRAD(:,:,:,2) = GY_W_VW(XWT(:,:,:), XDXX,XDZZ,XDZX,1,IKU,1)
  ZHGRAD(:,:,:,3) = GX_M_M(XTHT(:,:,:), XDXX,XDZZ,XDZX,1,IKU,1)
  ZHGRAD(:,:,:,4) = GY_M_M(XTHT(:,:,:), XDXX,XDZZ,XDZX,1,IKU,1)
  ZHGRAD(:,:,:,5) = GX_M_M(XRT(:,:,:,1), XDXX,XDZZ,XDZX,1,IKU,1)
  ZHGRAD(:,:,:,6) = GY_M_M(XRT(:,:,:,1), XDXX,XDZZ,XDZX,1,IKU,1)
END IF
   CALL TURB( CST,CSTURB, TBUCONF, TURBN, NEBN, YLDIMPHYEX,TLES, &
              NRR, NRRL, NRRI, CLBCX, CLBCY, IGRADIENTS, NHALO, NTURBSPLIT,          &
              LCLOUDMODIFLM, NSV, NSV_LGBEG, NSV_LGEND,                              &
              NSV_LIMA_NR, NSV_LIMA_NS, NSV_LIMA_NG, NSV_LIMA_NH,                    &
              L2D, LNOMIXLG,LFLAT,                                                   &
              LCOUPLES, LBLOWSNOW, LIBM,LFLYER,                                      &
              GCOMPUTE_SRC, XRSNOW,                                                  &
              LOCEAN, LDEEPOC, LDIAG_IN_RUN,                                         &
              CTURBLEN_CLOUD, CCLOUD,                                                &
              XTSTEP, TPFILE,                                                        &
              XDXX, XDYY, XDZZ, XDZX, XDZY, XZZ,                                     &
              XDIRCOSXW, XDIRCOSYW, XDIRCOSZW, XCOSSLOPE, XSINSLOPE,                 &
              XRHODJ, XTHVREF, ZHGRAD, XZS,                                          &
              ZSFTH, ZSFRV, ZSFSV, ZSFU, ZSFV,                                       &
              XPABST, XUT, XVT, XWT, XTKET, XSVT, XSRCT,                             &
              ZLENGTHM, ZLENGTHH, ZMFMOIST,                                          &
              XBL_DEPTH, XSBL_DEPTH,                                                 &
              XCEI, XCEI_MIN, XCEI_MAX, XCOEF_AMPL_SAT,                              &
              XTHT, XRT,                                                             &
              XRUS, XRVS, XRWS, XRTHS, XRRS, XRSVS, XRTKES, XSIGS, XWTHVMF,          &
              XTHW_FLUX, XRCW_FLUX, XSVW_FLUX,XDYP, XTHP, ZTDIFF, ZTDISS,            &
              TBUDGETS, KBUDGETS=SIZE(TBUDGETS),PLEM=XLEM,PRTKEMS=XRTKEMS,           &
              PTR=XTR, PDISS=XDISS, PCURRENT_TKE_DISS=XCURRENT_TKE_DISS,             &
              PIBM_LS=XIBM_LS(:,:,:,1), PIBM_XMUT=XIBM_XMUT,                         & 
              PSSTFL=XSSTFL, PSSTFL_C=XSSTFL_C, PSSRFL_C=XSSRFL_C,                   &
              PSSUFL_C=XSSUFL_C, PSSVFL_C=XSSVFL_C, PSSUFL=XSSUFL, PSSVFL=XSSVFL     )
!
DEALLOCATE(ZTDIFF)
DEALLOCATE(ZTDISS)
IF(LLEONARD) DEALLOCATE(ZHGRAD)
!
IF (LRMC01) THEN
  CALL ADD2DFIELD_ll( TZFIELDS_ll, XSBL_DEPTH, 'PHYS_PARAM_n::XSBL_DEPTH' )
  CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
  CALL CLEANLIST_ll(TZFIELDS_ll)
  IF ( CLBCX(1) /= "CYCL" .AND. LWEST_ll()) THEN
    XSBL_DEPTH(IIB-1,:)=XSBL_DEPTH(IIB,:)
  END IF
  IF ( CLBCX(2) /= "CYCL" .AND. LEAST_ll()) THEN
    XSBL_DEPTH(IIE+1,:)=XSBL_DEPTH(IIE,:)
  END IF
  IF ( CLBCY(1) /= "CYCL" .AND. LSOUTH_ll()) THEN
    XSBL_DEPTH(:,IJB-1)=XSBL_DEPTH(:,IJB)
  END IF
  IF ( CLBCY(2) /= "CYCL" .AND. LNORTH_ll()) THEN
    XSBL_DEPTH(:,IJE+1)=XSBL_DEPTH(:,IJE)
  END IF
END IF
!
CALL SECOND_MNH2(ZTIME3)
!
!-----------------------------------------------------------------------------
!
!*        7.    EDMF SCHEME
!               -----------
!
IF (CSCONV == 'EDKF') THEN
     ALLOCATE(ZEXN (IIU,IJU,IKU))
     ALLOCATE(ZSIGMF (IIU,IJU,IKU))
     ZSIGMF(:,:,:)=0.    
     ZEXN(:,:,:)=(XPABST(:,:,:)/XP00)**(XRD/XCPD)  
     !$20131113 check3d on ZEXN
     CALL MPPDB_CHECK3D(ZEXN,"physparan.7::ZEXN",PRECISION)
     CALL ADD3DFIELD_ll( TZFIELDS_ll, ZEXN, 'PHYS_PARAM_n::ZEXN' )
     !$20131113 add update_halo_ll
     CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
     CALL CLEANLIST_ll(TZFIELDS_ll)
     CALL MPPDB_CHECK3D(ZEXN,"physparam.7::ZEXN",PRECISION)
 !    
     CALL SHALLOW_MF_PACK(NRR,NRRL,NRRI,                                  &
                   TPFILE,ZTIME_LES_MF,                                   &
                   XTSTEP,                                                &
                   XDZZ, XZZ,XDXHAT(1),XDYHAT(1),                         &
                   XRHODJ, XRHODREF, XPABST, ZEXN, ZSFTH, ZSFRV,          &
                   XTHT,XRT,XUT,XVT,XTKET,XSVT,                           &
                   XRTHS,XRRS,XRUS,XRVS,XRSVS,                            &
                   ZSIGMF,XRC_MF, XRI_MF, XCF_MF, XWTHVMF)
!
ELSE
    XWTHVMF(:,:,:)=0.
    XRC_MF(:,:,:)=0.
    XRI_MF(:,:,:)=0.
    XCF_MF(:,:,:)=0.
ENDIF   
!
CALL SECOND_MNH2(ZTIME4)

  IF( LTRANS ) THEN
    XUT(:,:,:) = XUT(:,:,:) - XUTRANS
    XVT(:,:,:) = XVT(:,:,:) - XVTRANS
  END IF
  IF (CMF_CLOUD == 'STAT') THEN
    XSIGS =SQRT( XSIGS**2 + ZSIGMF**2 )
  ENDIF
  IF (CSCONV == 'EDKF') THEN
    DEALLOCATE(ZSIGMF)
    DEALLOCATE(ZEXN)
  ENDIF
END IF
!
IF (LLES_CALL) CALL SWITCH_SBG_LES_n
!
CALL SECOND_MNH2(ZTIME2)
!
PTURB = PTURB + ZTIME2 - ZTIME1 - (XTIME_LES-ZTIME_LES_MF) - XTIME_LES_BU_PROCESS &
      - XTIME_BU_PROCESS - (ZTIME4 - ZTIME3)
!
PMAFL = PMAFL + ZTIME4 - ZTIME3 - ZTIME_LES_MF
!
PTIME_BU = PTIME_BU + XTIME_LES_BU_PROCESS + XTIME_BU_PROCESS
!
!
!-------------------------------------------------------------------------------
!
!* deallocation of variables used in more than one parameterization
!
DEALLOCATE(ZSFU  )         ! surface schemes + turbulence
DEALLOCATE(ZSFV  )
DEALLOCATE(ZSFTH )
DEALLOCATE(ZSFRV )
DEALLOCATE(ZSFSV )
DEALLOCATE(ZSFCO2)
!
DEALLOCATE(ZSFTH_WALL )
DEALLOCATE(ZSFTH_ROOF )
DEALLOCATE(ZCD_ROOF )
DEALLOCATE(ZSFRV_WALL )
DEALLOCATE(ZSFRV_ROOF )
!-------------------------------------------------------------------------------
!
END SUBROUTINE PHYS_PARAM_n

