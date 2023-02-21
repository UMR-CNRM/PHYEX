!MNH_LIC Copyright 1995-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
! Modifications:
!  P. Wautelet 17/08/2020: add Budget_preallocate subroutine
!-----------------------------------------------------------------
module mode_ini_budget

  use mode_msg

  implicit none

  private

  public :: Budget_preallocate, Ini_budget

  integer, parameter :: NSOURCESMAX = 60 !Maximum number of sources in a budget

contains

subroutine Budget_preallocate()

use modd_budget, only: nbudgets, tbudgets,                                         &
                       NBUDGET_U, NBUDGET_V, NBUDGET_W, NBUDGET_TH, NBUDGET_TKE,   &
                       NBUDGET_RV, NBUDGET_RC, NBUDGET_RR, NBUDGET_RI, NBUDGET_RS, &
                       NBUDGET_RG, NBUDGET_RH, NBUDGET_SV1
use modd_nsv,    only: nsv, tsvlist

integer          :: ibudget
integer          :: jsv

call Print_msg( NVERB_DEBUG, 'BUD', 'Budget_preallocate', 'called' )

if ( allocated( tbudgets ) ) then
  call Print_msg( NVERB_WARNING, 'BUD', 'Budget_preallocate', 'tbudgets already allocated' )
  return
end if

nbudgets = NBUDGET_SV1 - 1 + nsv
allocate( tbudgets( nbudgets ) )

tbudgets(NBUDGET_U)%cname    = "UU"
tbudgets(NBUDGET_U)%ccomment = "Budget for U"
tbudgets(NBUDGET_U)%nid      = NBUDGET_U

tbudgets(NBUDGET_V)%cname    = "VV"
tbudgets(NBUDGET_V)%ccomment = "Budget for V"
tbudgets(NBUDGET_V)%nid      = NBUDGET_V

tbudgets(NBUDGET_W)%cname    = "WW"
tbudgets(NBUDGET_W)%ccomment = "Budget for W"
tbudgets(NBUDGET_W)%nid      = NBUDGET_W

tbudgets(NBUDGET_TH)%cname    = "TH"
tbudgets(NBUDGET_TH)%ccomment = "Budget for potential temperature"
tbudgets(NBUDGET_TH)%nid      = NBUDGET_TH

tbudgets(NBUDGET_TKE)%cname    = "TK"
tbudgets(NBUDGET_TKE)%ccomment = "Budget for turbulent kinetic energy"
tbudgets(NBUDGET_TKE)%nid      = NBUDGET_TKE

tbudgets(NBUDGET_RV)%cname    = "RV"
tbudgets(NBUDGET_RV)%ccomment = "Budget for water vapor mixing ratio"
tbudgets(NBUDGET_RV)%nid      = NBUDGET_RV

tbudgets(NBUDGET_RC)%cname    = "RC"
tbudgets(NBUDGET_RC)%ccomment = "Budget for cloud water mixing ratio"
tbudgets(NBUDGET_RC)%nid      = NBUDGET_RC

tbudgets(NBUDGET_RR)%cname    = "RR"
tbudgets(NBUDGET_RR)%ccomment = "Budget for rain water mixing ratio"
tbudgets(NBUDGET_RR)%nid      = NBUDGET_RR

tbudgets(NBUDGET_RI)%cname    = "RI"
tbudgets(NBUDGET_RI)%ccomment = "Budget for cloud ice mixing ratio"
tbudgets(NBUDGET_RI)%nid      = NBUDGET_RI

tbudgets(NBUDGET_RS)%cname    = "RS"
tbudgets(NBUDGET_RS)%ccomment = "Budget for snow/aggregate mixing ratio"
tbudgets(NBUDGET_RS)%nid      = NBUDGET_RS

tbudgets(NBUDGET_RG)%cname    = "RG"
tbudgets(NBUDGET_RG)%ccomment = "Budget for graupel mixing ratio"
tbudgets(NBUDGET_RG)%nid      = NBUDGET_RG

tbudgets(NBUDGET_RH)%cname    = "RH"
tbudgets(NBUDGET_RH)%ccomment = "Budget for hail mixing ratio"
tbudgets(NBUDGET_RH)%nid      = NBUDGET_RH

do jsv = 1, nsv
  ibudget = NBUDGET_SV1 - 1 + jsv
  tbudgets(ibudget)%cname    = Trim( tsvlist(jsv)%cmnhname )
  tbudgets(ibudget)%ccomment = 'Budget for scalar variable ' // Trim( tsvlist(jsv)%cmnhname )
  tbudgets(ibudget)%nid      = ibudget
end do


end subroutine Budget_preallocate


!     #################################################################
      SUBROUTINE Ini_budget(KLUOUT,PTSTEP,KSV,KRR,                    &
      ONUMDIFU,ONUMDIFTH,ONUMDIFSV,                                   &
      OHORELAX_UVWTH,OHORELAX_RV,OHORELAX_RC,OHORELAX_RR,             &
      OHORELAX_RI,OHORELAX_RS, OHORELAX_RG, OHORELAX_RH,OHORELAX_TKE, &
      OHORELAX_SV, OVE_RELAX, ove_relax_grd, OCHTRANS,                &
      ONUDGING,ODRAGTREE,ODEPOTREE, OAERO_EOL,                        &
      HRAD,HDCONV,HSCONV,HTURB,HTURBDIM,HCLOUD                        )
!     #################################################################
!
!!****  *INI_BUDGET* - routine to initialize the parameters for the budgets
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to set or compute the parameters used
!     by the MESONH budgets. Names of files for budget recording are processed 
!     and storage arrays are initialized.               
!
!!**  METHOD
!!    ------
!!      The essential of information is passed by modules. The choice of budgets 
!!    and processes set by the user as integers is converted in "actions" 
!!    readable  by the subroutine BUDGET under the form of string characters. 
!!    For each complete process composed of several elementary processes, names 
!!    of elementary processes are concatenated in order to have an explicit name
!!    in the comment of the recording file for budget. 
!!
!!      
!!    EXTERNAL
!!    --------   
!!      None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Modules MODD_*
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INI_BUDGET)
!!      
!!
!!    AUTHOR
!!    ------
!!	P. Hereil      * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        01/03/95 
!!      J. Stein        25/06/95  put the sources in phase with the code
!!      J. Stein        20/07/95  reset to FALSE of all the switches when
!!                                CBUTYPE /= MASK or CART
!!      J. Stein        26/06/96  add the new sources + add the increment between
!!                                2 active processes
!!      J.-P. Pinty     13/12/96  Allowance of multiple SVs
!!      J.-P. Pinty     11/01/97  Includes deep convection ice and forcing processes
!!      J.-P. Lafore    10/02/98  Allocation of the RHODJs for budget
!!      V. Ducrocq      04/06/99  //  
!!      N. Asencio      18/06/99  // MASK case : delete KIMAX and KJMAX arguments,
!!                                GET_DIM_EXT_ll initializes the dimensions of the
!!                                extended local domain.
!!                                LBU_MASK and NBUSURF are allocated on the extended
!!                                local domain.
!!                                add 3 local variables IBUDIM1,IBUDIM2,IBUDIM3
!!                                to define the dimensions of the budget arrays
!!                                in the different cases CART and MASK
!!      J.-P. Pinty     23/09/00  add budget for C2R2
!!      V. Masson       18/11/02  add budget for 2way nesting
!!      O.Geoffroy      03/2006   Add KHKO scheme
!!      J.-P. Pinty     22/04/97  add the explicit hail processes
!!      C.Lac           10/08/07  Add ADV for PPM without contribution
!!                                of each direction
!!      C. Barthe       19/11/09  Add atmospheric electricity
!!      C.Lac           01/07/11  Add vegetation drag        
!!      P. Peyrille, M. Tomasini : include in the forcing term the 2D forcing
!!                                terms in term 2DFRC search for modif PP . but Not very clean! 
!!      C .Lac          27/05/14    add negativity corrections for chemical species
!!      C.Lac           29/01/15  Correction for NSV_USER
!!      J.Escobar       02/10/2015 modif for JPHEXT(JPVEXT) variable  
!!      C.Lac           04/12/15  Correction for LSUPSAT 
!  C. Lac         04/2016: negative contribution to the budget split between advection, turbulence and microphysics for KHKO/C2R2
!  C. Barthe      01/2016: add budget for LIMA
!  C. Lac         10/2016: add budget for droplet deposition
!  S. Riette      11/2016: new budgets for ICE3/ICE4
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 15/11/2019: remove unused CBURECORD variable
!  P. Wautelet 24/02/2020: bugfix: corrected condition for budget NCDEPITH
!  P. Wautelet 26/02/2020: bugfix: rename CEVA->REVA for budget for raindrop evaporation in C2R2 (necessary after commit 4ed805fc)
!  P. Wautelet 26/02/2020: bugfix: add missing condition on OCOLD for NSEDIRH budget in LIMA case
!  P. Wautelet 02-03/2020: use the new data structures and subroutines for budgets
!  B. Vie      02/03/2020: LIMA negativity checks after turbulence, advection and microphysics budgets
!  P .Wautelet 09/03/2020: add missing budgets for electricity
!  P. Wautelet 25/03/2020: add missing ove_relax_grd
!  P. Wautelet 23/04/2020: add nid in tbudgetdata datatype
!  P. Wautelet + Benoit Vié 11/06/2020: improve removal of negative scalar variables + adapt the corresponding budgets
!  P. Wautelet 30/06/2020: use NADVSV when possible
!  P. Wautelet 30/06/2020: add NNETURSV, NNEADVSV and NNECONSV variables
!  P. Wautelet 06/07/2020: bugfix: add condition on HTURB for NETUR sources for SV budgets
!  P. Wautelet 08/12/2020: add nbusubwrite and nbutotwrite
!  P. Wautelet 11/01/2021: ignore xbuwri for cartesian boxes (write at every xbulen interval)
!  P. Wautelet 01/02/2021: bugfix: add missing CEDS source terms for SV budgets
!  P. Wautelet 02/02/2021: budgets: add missing source terms for SV budgets in LIMA
!  P. Wautelet 03/02/2021: budgets: add new source if LIMA splitting: CORR2
!  P. Wautelet 10/02/2021: budgets: add missing sources for NSV_C2R2BEG+3 budget
!  P. Wautelet 11/02/2021: budgets: add missing term SCAV for NSV_LIMA_SCAVMASS budget
!  P. Wautelet 02/03/2021: budgets: add terms for blowing snow
!  P. Wautelet 04/03/2021: budgets: add terms for drag due to buildings
!  P. Wautelet 17/03/2021: choose source terms for budgets with character strings instead of multiple integer variables
!  C. Barthe   14/03/2022: budgets: add terms for CIBU and RDSF in LIMA
!  M. Taufour  01/07/2022: budgets: add concentration for snow, graupel, hail
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
!
use modd_2d_frc,        only: l2d_adv_frc, l2d_rel_frc
use modd_blowsnow,      only: lblowsnow
use modd_blowsnow_n,    only: lsnowsubl
use modd_budget
use modd_ch_aerosol,    only: lorilam
use modd_conf,          only: l1d, lcartesian, lforcing, lthinshell, nmodel
use modd_dim_n,         only: nimax_ll, njmax_ll, nkmax
use modd_dragbldg_n,    only: ldragbldg
use modd_dust,          only: ldust
use modd_dyn,           only: lcorio, xseglen
use modd_dyn_n,         only: xtstep, locean
use modd_elec_descr,    only: linductive, lrelax2fw_ion
use modd_field,         only: TYPEREAL
use modd_fire_n,        only: lblaze
use modd_nsv,           only: nsv_aerbeg, nsv_aerend, nsv_aerdepbeg, nsv_aerdepend, nsv_c2r2beg, nsv_c2r2end,      &
                              nsv_chembeg, nsv_chemend, nsv_chicbeg, nsv_chicend, nsv_csbeg, nsv_csend,            &
                              nsv_dstbeg, nsv_dstend, nsv_dstdepbeg, nsv_dstdepend, nsv_elecbeg, nsv_elecend,      &
#ifdef MNH_FOREFIRE
                              nsv_ffbeg, nsv_ffend,                                                                &
#endif
                              nsv_lgbeg, nsv_lgend,                                                                &
                              nsv_lima_beg, nsv_lima_end, nsv_lima_ccn_acti, nsv_lima_ccn_free, nsv_lima_hom_haze, &
                              nsv_lima_ifn_free, nsv_lima_ifn_nucl, nsv_lima_imm_nucl,                             &
                              nsv_lima_nc, nsv_lima_nr, nsv_lima_ni, nsv_lima_ns, nsv_lima_ng, nsv_lima_nh,        &
                              nsv_lima_scavmass, nsv_lima_spro,                                                    &
                              nsv_lnoxbeg, nsv_lnoxend, nsv_ppbeg, nsv_ppend,                                      &
                              nsv_sltbeg, nsv_sltend, nsv_sltdepbeg, nsv_sltdepend, nsv_snwbeg, nsv_snwend,        &
                              nsv_user, tsvlist
use modd_parameters,   only: jphext
use modd_param_c2r2,   only: ldepoc_c2r2 => ldepoc, lrain_c2r2 => lrain, lsedc_c2r2 => lsedc, lsupsat_c2r2 => lsupsat
use modd_param_ice,    only: ladj_after, ladj_before, ldeposc_ice => ldeposc, lred, lsedic_ice => lsedic, lwarm_ice => lwarm
use modd_param_n,      only: cactccn, celec
use modd_param_lima,   only: laero_mass_lima => laero_mass, lacti_lima => lacti, ldepoc_lima => ldepoc, &
                             lhhoni_lima => lhhoni, lmeyers_lima => lmeyers, lnucl_lima => lnucl,       &
                             lptsplit,                                                                                       &
                             lscav_lima => lscav, lsedc_lima => lsedc, lsedi_lima => lsedi,             &
                             lspro_lima => lspro,  lcibu, lrdsf,           &
                             nmom_c, nmom_r, nmom_i, nmom_s, nmom_g, nmom_h, nmod_ccn, nmod_ifn, nmod_imm
use modd_ref,          only: lcouples
use modd_salt,         only: lsalt
use modd_turb_n,       only: lsubg_cond
use modd_viscosity,    only: lvisc, lvisc_r, lvisc_sv, lvisc_th, lvisc_uvw

USE MODE_ll

IMPLICIT NONE
!
!*       0.1   declarations of argument
!
!
INTEGER,         INTENT(IN) :: KLUOUT   ! Logical unit number for prints
REAL, INTENT(IN)    :: PTSTEP           ! time step
INTEGER, INTENT(IN) :: KSV              ! number of scalar variables
INTEGER, INTENT(IN) :: KRR              ! number of moist variables
LOGICAL, INTENT(IN) :: ONUMDIFU         ! switch to activate the numerical
                                        ! diffusion for momentum
LOGICAL, INTENT(IN) :: ONUMDIFTH        ! for meteorological scalar variables
LOGICAL, INTENT(IN) :: ONUMDIFSV        ! for tracer scalar variables
LOGICAL, INTENT(IN) :: OHORELAX_UVWTH  ! switch for the
                       ! horizontal relaxation for U,V,W,TH
LOGICAL, INTENT(IN) :: OHORELAX_RV     ! switch for the
                       ! horizontal relaxation for Rv
LOGICAL, INTENT(IN) :: OHORELAX_RC     ! switch for the
                       ! horizontal relaxation for Rc
LOGICAL, INTENT(IN) :: OHORELAX_RR     ! switch for the
                       ! horizontal relaxation for Rr
LOGICAL, INTENT(IN) :: OHORELAX_RI     ! switch for the
                       ! horizontal relaxation for Ri
LOGICAL, INTENT(IN) :: OHORELAX_RS     ! switch for the
                       ! horizontal relaxation for Rs
LOGICAL, INTENT(IN) :: OHORELAX_RG     ! switch for the
                       ! horizontal relaxation for Rg
LOGICAL, INTENT(IN) :: OHORELAX_RH     ! switch for the
                       ! horizontal relaxation for Rh
LOGICAL, INTENT(IN) :: OHORELAX_TKE    ! switch for the
                       ! horizontal relaxation for tke
LOGICAL,DIMENSION(:),INTENT(IN):: OHORELAX_SV     ! switch for the
                       ! horizontal relaxation for scalar variables
LOGICAL, INTENT(IN) :: OVE_RELAX        ! switch to activate the vertical
                                        ! relaxation
logical, intent(in) :: ove_relax_grd    ! switch to activate the vertical
                                        ! relaxation to the lowest verticals
LOGICAL, INTENT(IN) :: OCHTRANS         ! switch to activate convective
                                        !transport for SV
LOGICAL, INTENT(IN) :: ONUDGING         ! switch to activate nudging
LOGICAL, INTENT(IN) :: ODRAGTREE        ! switch to activate vegetation drag
LOGICAL, INTENT(IN) :: ODEPOTREE        ! switch to activate droplet deposition on tree
LOGICAL, INTENT(IN) :: OAERO_EOL        ! switch to activate wind turbine wake
CHARACTER (LEN=*), INTENT(IN) :: HRAD   ! type of the radiation scheme
CHARACTER (LEN=*), INTENT(IN) :: HDCONV ! type of the deep convection scheme
CHARACTER (LEN=*), INTENT(IN) :: HSCONV ! type of the shallow convection scheme
CHARACTER (LEN=*), INTENT(IN) :: HTURB  ! type of the turbulence scheme
CHARACTER (LEN=*), INTENT(IN) :: HTURBDIM! dimensionnality of the turbulence 
                                        ! scheme
CHARACTER (LEN=*), INTENT(IN) :: HCLOUD ! type of microphysical scheme
!
!*       0.2   declarations of local variables
!
real, parameter :: ITOL = 1e-6

INTEGER :: JI, JJ                                         ! loop indices
INTEGER :: IIMAX_ll, IJMAX_ll ! size of the physical global domain
INTEGER :: IIU, IJU                                       ! size along x and y directions
                                                          ! of the extended subdomain
INTEGER :: IBUDIM1                                        ! first dimension of the budget arrays
                                                          ! = NBUIMAX in CART case
                                                          ! = NBUKMAX in MASK case
INTEGER :: IBUDIM2                                        ! second dimension of the budget arrays
                                                          ! = NBUJMAX in CART case
                                                          ! = nbusubwrite in MASK case
INTEGER :: IBUDIM3                                        ! third dimension of the budget arrays
                                                          ! = NBUKMAX in CART case
                                                          ! = NBUMASK in MASK case
INTEGER             :: JSV               ! loop indice for the SVs
INTEGER             :: IINFO_ll ! return status of the interface routine
integer             :: ibudget
logical             :: gtmp
type(tbusourcedata) :: tzsource ! Used to prepare metadate of source terms

call Print_msg( NVERB_DEBUG, 'BUD', 'Ini_budget', 'called' )
!
!*       1.    COMPUTE BUDGET VARIABLES
!              ------------------------
!
NBUSTEP = NINT (XBULEN / PTSTEP)
NBUTSHIFT=0
!
!  common dimension for all CBUTYPE values
!
IF (LBU_KCP) THEN
  NBUKMAX = 1
ELSE
  NBUKMAX = NBUKH - NBUKL +1
END IF
!
if ( cbutype == 'CART' .or. cbutype == 'MASK' ) then
  !Check if xbulen is a multiple of xtstep (within tolerance)
  if ( Abs( Nint( xbulen / xtstep ) * xtstep - xbulen ) > ( ITOL * xtstep ) ) &
    call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget', 'xbulen is not a multiple of xtstep' )

  if ( cbutype == 'CART' ) then
    !Check if xseglen is a multiple of xbulen (within tolerance)
    if ( Abs( Nint( xseglen / xbulen ) * xbulen - xseglen ) > ( ITOL * xseglen ) ) &
      call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget', 'xseglen is not a multiple of xbulen' )

    !Write cartesian budgets every xbulen time period (do not take xbuwri into account)
    xbuwri = xbulen

    nbusubwrite = 1                                      !Number of budget time average periods for each write
    nbutotwrite = nbusubwrite * Nint( xseglen / xbulen ) !Total number of budget time average periods
  else if ( cbutype == 'MASK' ) then
    !Check if xbuwri is a multiple of xtstep (within tolerance)
    if ( Abs( Nint( xbuwri / xtstep ) * xtstep - xbuwri ) > ( ITOL * xtstep ) ) &
      call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget', 'xbuwri is not a multiple of xtstep' )

    !Check if xbuwri is a multiple of xbulen (within tolerance)
    if ( Abs( Nint( xbuwri / xbulen ) * xbulen - xbuwri ) > ( ITOL * xbulen ) ) &
      call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget', 'xbuwri is not a multiple of xbulen' )

    !Check if xseglen is a multiple of xbuwri (within tolerance)
    if ( Abs( Nint( xseglen / xbuwri ) * xbuwri - xseglen ) > ( ITOL * xseglen ) ) &
      call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget', 'xseglen is not a multiple of xbuwri' )

    nbusubwrite = Nint ( xbuwri / xbulen )               !Number of budget time average periods for each write
    nbutotwrite = nbusubwrite * Nint( xseglen / xbuwri ) !Total number of budget time average periods
  end if
end if

IF (CBUTYPE=='CART') THEN              ! cartesian case only
!
  IF ( NBUIL < 1 )        CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUIL too small (<1)' )
  IF ( NBUIL > NIMAX_ll ) CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUIL too large (>NIMAX)' )
  IF ( NBUIH < 1 )        CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUIH too small (<1)' )
  IF ( NBUIH > NIMAX_ll ) CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUIH too large (>NIMAX)' )
  IF ( NBUIH < NBUIL )    CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUIH < NBUIL' )
  IF (LBU_ICP) THEN
    NBUIMAX_ll = 1
  ELSE
    NBUIMAX_ll = NBUIH - NBUIL +1
  END IF

  IF ( NBUJL < 1 )        CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUJL too small (<1)' )
  IF ( NBUJL > NJMAX_ll ) CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUJL too large (>NJMAX)' )
  IF ( NBUJH < 1 )        CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUJH too small (<1)' )
  IF ( NBUJH > NJMAX_ll ) CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUJH too large (>NJMAX)' )
  IF ( NBUJH < NBUJL )    CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUJH < NBUJL' )
  IF (LBU_JCP) THEN
    NBUJMAX_ll = 1
  ELSE
    NBUJMAX_ll = NBUJH - NBUJL +1
  END IF

  IF ( NBUKL < 1 )     CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUKL too small (<1)' )
  IF ( NBUKL > NKMAX ) CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUKL too large (>NKMAX)' )
  IF ( NBUKH < 1 )     CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUKH too small (<1)' )
  IF ( NBUKH > NKMAX ) CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUKH too large (>NKMAX)' )
  IF ( NBUKH < NBUKL ) CALL Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'NBUKH < NBUKL' )

  CALL GET_INTERSECTION_ll(NBUIL+JPHEXT,NBUJL+JPHEXT,NBUIH+JPHEXT,NBUJH+JPHEXT, &
      NBUSIL,NBUSJL,NBUSIH,NBUSJH,"PHYS",IINFO_ll)
  IF ( IINFO_ll /= 1 ) THEN ! 
    IF (LBU_ICP) THEN 
      NBUIMAX = 1
    ELSE
      NBUIMAX = NBUSIH - NBUSIL +1
    END IF
    IF (LBU_JCP) THEN 
      NBUJMAX = 1
    ELSE
      NBUJMAX =  NBUSJH - NBUSJL +1
    END IF
  ELSE ! the intersection is void 
    CBUTYPE='SKIP'  ! no budget on this processor       
    NBUIMAX = 0     ! in order to allocate void arrays
    NBUJMAX = 0
  ENDIF
! three first dimensions of budget arrays in cart and skip cases
   IBUDIM1=NBUIMAX
   IBUDIM2=NBUJMAX
   IBUDIM3=NBUKMAX
! these variables are not be used 
   NBUMASK=-1
!
ELSEIF (CBUTYPE=='MASK') THEN          ! mask case only 
!
  LBU_ENABLE=.TRUE.
                                    ! result on the FM_FILE
  NBUTIME = 1

  CALL GET_DIM_EXT_ll ('B', IIU,IJU)
  ALLOCATE( LBU_MASK( IIU ,IJU, NBUMASK) )
  LBU_MASK(:,:,:)=.FALSE.
  ALLOCATE( NBUSURF( IIU, IJU, NBUMASK, nbusubwrite) )
  NBUSURF(:,:,:,:) = 0
!
! three first dimensions of budget arrays in mask case
!  the order of the dimensions are the order expected in WRITE_DIACHRO routine:
!  x,y,z,time,mask,processus  and in this case x and y are missing
!  first dimension of the arrays : dimension along K
!  second dimension of the arrays : number of the budget time period
!  third dimension of the arrays : number of the budget masks zones
  IBUDIM1=NBUKMAX
  IBUDIM2=nbusubwrite
  IBUDIM3=NBUMASK
! these variables are not used in this case
  NBUIMAX=-1
  NBUJMAX=-1
! the beginning and the end along x and y direction : global extended domain
 ! get dimensions of the physical global domain
   CALL GET_GLOBALDIMS_ll (IIMAX_ll,IJMAX_ll)
   NBUIL=1
   NBUIH=IIMAX_ll + 2 * JPHEXT
   NBUJL=1 
   NBUJH=IJMAX_ll + 2 * JPHEXT
!
ELSE                      ! default case
!
  LBU_ENABLE=.FALSE.
  NBUIMAX = -1
  NBUJMAX = -1
  LBU_RU = .FALSE.
  LBU_RV = .FALSE.
  LBU_RW = .FALSE.
  LBU_RTH= .FALSE.
  LBU_RTKE= .FALSE.
  LBU_RRV= .FALSE.
  LBU_RRC= .FALSE.
  LBU_RRR= .FALSE.
  LBU_RRI= .FALSE.
  LBU_RRS= .FALSE.
  LBU_RRG= .FALSE.
  LBU_RRH= .FALSE.
  LBU_RSV= .FALSE.
!
! three first dimensions of budget arrays in default case
  IBUDIM1=0
  IBUDIM2=0
  IBUDIM3=0
!
END IF  
!
!
!-------------------------------------------------------------------------------
!
!*       2.    ALLOCATE MEMORY FOR BUDGET ARRAYS AND INITIALIZE
!              ------------------------------------------------
!
LBU_BEG =.TRUE. 
!
!-------------------------------------------------------------------------------
!
!*       3.    INITALIZE VARIABLES
!              -------------------
!
!Create intermediate variable to store rhodj for scalar variables
if ( lbu_rth .or. lbu_rtke .or. lbu_rrv .or. lbu_rrc .or. lbu_rrr .or. &
     lbu_rri .or. lbu_rrs  .or. lbu_rrg .or. lbu_rrh .or. lbu_rsv      ) then
  allocate( tburhodj )

  tburhodj%cmnhname  = 'RhodJS'
  tburhodj%cstdname  = ''
  tburhodj%clongname = 'RhodJS'
  tburhodj%cunits    = 'kg'
  tburhodj%ccomment  = 'RhodJ for Scalars variables'
  tburhodj%ngrid     = 1
  tburhodj%ntype     = TYPEREAL
  tburhodj%ndims     = 3

  allocate( tburhodj%xdata(ibudim1, ibudim2, ibudim3) )
  tburhodj%xdata(:, :, :) = 0.
end if


tzsource%ntype    = TYPEREAL
tzsource%ndims    = 3

! Budget of RU
tbudgets(NBUDGET_U)%lenabled = lbu_ru

if ( lbu_ru ) then
  allocate( tbudgets(NBUDGET_U)%trhodj )

  tbudgets(NBUDGET_U)%trhodj%cmnhname  = 'RhodJX'
  tbudgets(NBUDGET_U)%trhodj%cstdname  = ''
  tbudgets(NBUDGET_U)%trhodj%clongname = 'RhodJX'
  tbudgets(NBUDGET_U)%trhodj%cunits    = 'kg'
  tbudgets(NBUDGET_U)%trhodj%ccomment  = 'RhodJ for momentum along X axis'
  tbudgets(NBUDGET_U)%trhodj%ngrid     = 2
  tbudgets(NBUDGET_U)%trhodj%ntype     = TYPEREAL
  tbudgets(NBUDGET_U)%trhodj%ndims     = 3

  allocate( tbudgets(NBUDGET_U)%trhodj%xdata(ibudim1, ibudim2, ibudim3) )
  tbudgets(NBUDGET_U)%trhodj%xdata(:, :, :) = 0.

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_U)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_U)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_U)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_U)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of momentum along X axis'
  tzsource%ngrid    = 2

  tzsource%cunits   = 'm s-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 'm s-2'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'NUD'
  tzsource%clongname  = 'nudging'
  tzsource%lavailable = onudging
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'CURV'
  tzsource%clongname  = 'curvature'
  tzsource%lavailable = .not.l1d .and. .not.lcartesian
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'COR'
  tzsource%clongname  = 'Coriolis'
  tzsource%lavailable = lcorio
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifu
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_uvwth .or. ove_relax .or. ove_relax_grd
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'DRAG'
  tzsource%clongname  = 'drag force due to trees'
  tzsource%lavailable = odragtree
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'DRAGEOL'
  tzsource%clongname  = 'drag force due to wind turbine'
  tzsource%lavailable = OAERO_EOL
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'DRAGB'
  tzsource%clongname  = 'drag force due to buildings'
  tzsource%lavailable = ldragbldg
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'VTURB'
  tzsource%clongname  = 'vertical turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL'
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'HTURB'
  tzsource%clongname  = 'horizontal turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL' .and. HTURBDIM == '3DIM'
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'MAFL'
  tzsource%clongname  = 'mass flux'
  tzsource%lavailable = hsconv == 'EDKF'
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_uvw
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )

  tzsource%cmnhname   = 'PRES'
  tzsource%clongname  = 'pressure'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_U), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_U) )

  call Sourcelist_scan( tbudgets(NBUDGET_U), cbulist_ru )
end if

! Budget of RV
tbudgets(NBUDGET_V)%lenabled = lbu_rv

if ( lbu_rv ) then
  allocate( tbudgets(NBUDGET_V)%trhodj )

  tbudgets(NBUDGET_V)%trhodj%cmnhname  = 'RhodJY'
  tbudgets(NBUDGET_V)%trhodj%cstdname  = ''
  tbudgets(NBUDGET_V)%trhodj%clongname = 'RhodJY'
  tbudgets(NBUDGET_V)%trhodj%cunits    = 'kg'
  tbudgets(NBUDGET_V)%trhodj%ccomment  = 'RhodJ for momentum along Y axis'
  tbudgets(NBUDGET_V)%trhodj%ngrid     = 3
  tbudgets(NBUDGET_V)%trhodj%ntype     = TYPEREAL
  tbudgets(NBUDGET_V)%trhodj%ndims     = 3

  allocate( tbudgets(NBUDGET_V)%trhodj%xdata(ibudim1, ibudim2, ibudim3) )
  tbudgets(NBUDGET_V)%trhodj%xdata(:, :, :) = 0.

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_V)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_V)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_V)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_V)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of momentum along Y axis'
  tzsource%ngrid    = 3

  tzsource%cunits   = 'm s-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 'm s-2'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'NUD'
  tzsource%clongname  = 'nudging'
  tzsource%lavailable = onudging
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'CURV'
  tzsource%clongname  = 'curvature'
  tzsource%lavailable = .not.l1d .and. .not.lcartesian
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'COR'
  tzsource%clongname  = 'Coriolis'
  tzsource%lavailable = lcorio
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifu
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_uvwth .or. ove_relax .or. ove_relax_grd
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'DRAG'
  tzsource%clongname  = 'drag force due to trees'
  tzsource%lavailable = odragtree
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'DRAGEOL'
  tzsource%clongname  = 'drag force due to wind turbine'
  tzsource%lavailable = OAERO_EOL
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'DRAGB'
  tzsource%clongname  = 'drag force due to buildings'
  tzsource%lavailable = ldragbldg
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'VTURB'
  tzsource%clongname  = 'vertical turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL'
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'HTURB'
  tzsource%clongname  = 'horizontal turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL' .and. HTURBDIM == '3DIM'
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'MAFL'
  tzsource%clongname  = 'mass flux'
  tzsource%lavailable = hsconv == 'EDKF'
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_uvw
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )

  tzsource%cmnhname   = 'PRES'
  tzsource%clongname  = 'pressure'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_V), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_V) )

  call Sourcelist_scan( tbudgets(NBUDGET_V), cbulist_rv )
end if

! Budget of RW
tbudgets(NBUDGET_W)%lenabled = lbu_rw

if ( lbu_rw ) then
  allocate( tbudgets(NBUDGET_W)%trhodj )

  tbudgets(NBUDGET_W)%trhodj%cmnhname  = 'RhodJZ'
  tbudgets(NBUDGET_W)%trhodj%cstdname  = ''
  tbudgets(NBUDGET_W)%trhodj%clongname = 'RhodJZ'
  tbudgets(NBUDGET_W)%trhodj%cunits    = 'kg'
  tbudgets(NBUDGET_W)%trhodj%ccomment  = 'RhodJ for momentum along Z axis'
  tbudgets(NBUDGET_W)%trhodj%ngrid     = 4
  tbudgets(NBUDGET_W)%trhodj%ntype     = TYPEREAL
  tbudgets(NBUDGET_W)%trhodj%ndims     = 3

  allocate( tbudgets(NBUDGET_W)%trhodj%xdata(ibudim1, ibudim2, ibudim3) )
  tbudgets(NBUDGET_W)%trhodj%xdata(:, :, :) = 0.

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_W)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_W)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_W)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_W)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of momentum along Z axis'
  tzsource%ngrid    = 4

  tzsource%cunits   = 'm s-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 'm s-2'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'NUD'
  tzsource%clongname  = 'nudging'
  tzsource%lavailable = onudging
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'CURV'
  tzsource%clongname  = 'curvature'
  tzsource%lavailable = .not.l1d .and. .not.lcartesian .and. .not.lthinshell
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'COR'
  tzsource%clongname  = 'Coriolis'
  tzsource%lavailable = lcorio .and. .not.l1d .and. .not.lthinshell
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifu
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_uvwth .or. ove_relax .or. ove_relax_grd
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'VTURB'
  tzsource%clongname  = 'vertical turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL'
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'HTURB'
  tzsource%clongname  = 'horizontal turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL' .and. HTURBDIM == '3DIM'
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_uvw
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'GRAV'
  tzsource%clongname  = 'gravity'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'PRES'
  tzsource%clongname  = 'pressure'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  tzsource%cmnhname   = 'DRAGEOL'
  tzsource%clongname  = 'drag force due to wind turbine'
  tzsource%lavailable = OAERO_EOL
  call Budget_source_add( tbudgets(NBUDGET_W), tzsource )

  call Sourcelist_sort_compact( tbudgets(NBUDGET_W) )

  call Sourcelist_scan( tbudgets(NBUDGET_W), cbulist_rw )
end if

! Budget of RTH
tbudgets(NBUDGET_TH)%lenabled = lbu_rth

if ( lbu_rth ) then
  tbudgets(NBUDGET_TH)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_TH)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_TH)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_TH)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_TH)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of potential temperature'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'K'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 'K s-1'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = '2DADV'
  tzsource%clongname  = 'advective forcing'
  tzsource%lavailable = l2d_adv_frc
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = '2DREL'
  tzsource%clongname  = 'relaxation forcing'
  tzsource%lavailable = l2d_rel_frc
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'NUD'
  tzsource%clongname  = 'nudging'
  tzsource%lavailable = onudging
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'PREF'
  tzsource%clongname  = 'reference pressure'
  tzsource%lavailable = krr > 0 .and. .not.l1d
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_uvwth .or. ove_relax .or. ove_relax_grd
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'RAD'
  tzsource%clongname  = 'radiation'
  tzsource%lavailable = hrad /= 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DCONV'
  tzsource%clongname  = 'KAFR convection'
  tzsource%lavailable = hdconv == 'KAFR' .OR. hsconv == 'KAFR'

  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'BLAZE'
  tzsource%clongname  = 'blaze fire model contribution'
  tzsource%lavailable = lblaze
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'VTURB'
  tzsource%clongname  = 'vertical turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HTURB'
  tzsource%clongname  = 'horizontal turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL' .and. HTURBDIM == '3DIM'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DISSH'
  tzsource%clongname  = 'dissipation'
  tzsource%lavailable = hturb == 'TKEL'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'NETUR'
  tzsource%clongname  = 'negativity correction induced by turbulence'
  tzsource%lavailable = hturb == 'TKEL' .and. (      hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                                                .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' )
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'MAFL'
  tzsource%clongname  = 'mass flux'
  tzsource%lavailable = hsconv == 'EDKF'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'SNSUB'
  tzsource%clongname  = 'blowing snow sublimation'
  tzsource%lavailable = lblowsnow .and. lsnowsubl
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_th
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'OCEAN'
  tzsource%clongname  = 'radiative tendency due to SW penetrating ocean'
  tzsource%lavailable = locean .and. (.not. lcouples)
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'NEADV'
  tzsource%clongname  = 'negativity correction induced by advection'
  tzsource%lavailable =       hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                         .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'NEGA'
  tzsource%clongname  = 'negativity correction'
  tzsource%lavailable =       hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                         .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'SEDI'
  tzsource%clongname  = 'heat transport by hydrometeors sedimentation'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HENU'
  tzsource%clongname  = 'heterogeneous nucleation'
  gtmp = cactccn == 'ABRK' .and. (lorilam .or. ldust .or. lsalt )
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. nmom_c.ge.1 .and. lacti_lima .and. nmod_ccn >= 1   &
                                                     .and. ( .not.lptsplit .or. .not.lsubg_cond )          ) &
                        .or. ( hcloud      == 'C2R2' .and. ( gtmp .or. .not.lsupsat_c2r2 ) )                 &
                        .or. ( hcloud      == 'KHKO' .and. ( gtmp .or. .not.lsupsat_c2r2 ) )
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'REVA'
  tzsource%clongname  = 'rain evaporation'
  tzsource%lavailable =     ( hcloud      == 'LIMA' .and. ( ( .not. lptsplit .and. nmom_r.ge.1 ) .or. lptsplit ) ) &
                       .or. ( hcloud(1:3) == 'ICE'  .and. lwarm_ice )                                              &
                       .or. ( hcloud      == 'C2R2' .and. lrain_c2r2 )                                             &
                       .or. ( hcloud      == 'KHKO' .and. lrain_c2r2 )                                             &
                       .or.   hcloud      == 'KESS'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HIN'
  tzsource%clongname  = 'heterogeneous ice nucleation'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .or. (hcloud == 'LIMA' .and. nmom_i == 1)
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HIND'
  tzsource%clongname  = 'heterogeneous nucleation by deposition'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HINC'
  tzsource%clongname  = 'heterogeneous nucleation by contact'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HON'
  tzsource%clongname  = 'homogeneous nucleation'
  tzsource%lavailable = hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HONH'
  tzsource%clongname  = 'haze homogeneous nucleation'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima .and. lhhoni_lima .and. nmod_ccn >= 1
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HONC'
  tzsource%clongname  = 'droplet homogeneous freezing'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. lnucl_lima ) )
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HONR'
  tzsource%clongname  = 'raindrop homogeneous freezing'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. lnucl_lima .and. nmom_r.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'SFR'
  tzsource%clongname  = 'spontaneous freezing'
  tzsource%lavailable = hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DEPS'
  tzsource%clongname  = 'deposition on snow'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DEPG'
  tzsource%clongname  = 'deposition on graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DEPH'
  tzsource%clongname  = 'deposition on hail'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. nmom_h.ge.1 ) ) &
                        .or. hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'IMLT'
  tzsource%clongname  = 'melting of ice'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'BERFI'
  tzsource%clongname  = 'Bergeron-Findeisen'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'RIM'
  tzsource%clongname  = 'riming of cloud water'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'ACC'
  tzsource%clongname  = 'accretion of rain on aggregates'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. (        lptsplit                                                            &
                                                        .or. ( nmom_s.ge.1 .and. nmom_r.ge.1 ) ) ) &
                        .or.   hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'CFRZ'
  tzsource%clongname  = 'conversion freezing of rain'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'WETG'
  tzsource%clongname  = 'wet growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DRYG'
  tzsource%clongname  = 'dry growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'GMLT'
  tzsource%clongname  = 'graupel melting'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'WETH'
  tzsource%clongname  = 'wet growth of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1 &
                                                .and. nmom_c.ge.1    .and. nmom_s.ge.1 )                   &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4' 
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DRYH'
  tzsource%clongname  = 'dry growth of hail'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'HMLT'
  tzsource%clongname  = 'melting of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1 &
                                                .and. nmom_c.ge.1    .and. nmom_s.ge.1 ) &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'CORR'
  tzsource%clongname  = 'correction'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'CEDS'
  tzsource%clongname  = 'adjustment to saturation'
  tzsource%lavailable = hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'ADJU'
  tzsource%clongname  = 'adjustment to saturation'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. ladj_before .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'DEPI'
  tzsource%clongname  = 'deposition on ice'
  tzsource%lavailable =      ( hcloud(1:3) == 'ICE' .and. ( .not. lred .or. ( lred .and. ladj_after ) .or. celec /= 'NONE') ) &
                        .or. ( hcloud == 'LIMA' .and. lptsplit )
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'COND'
  tzsource%clongname  = 'vapor condensation or cloud water evaporation'
  tzsource%lavailable = hcloud == 'C2R2' .or. hcloud == 'KHKO' .or. hcloud == 'KESS' .or. hcloud == 'REVE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )

  tzsource%cmnhname   = 'NECON'
  tzsource%clongname  = 'negativity correction induced by condensation'
  tzsource%lavailable = (      hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4'   &
                          .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) &
                        .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_TH), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_TH) )

  call Sourcelist_scan( tbudgets(NBUDGET_TH), cbulist_rth )
end if

! Budget of RTKE
tbudgets(NBUDGET_TKE)%lenabled = lbu_rtke

if ( lbu_rtke ) then
  tbudgets(NBUDGET_TKE)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_TKE)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_TKE)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_TKE)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_TKE)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of turbulent kinetic energy'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'm2 s-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 'm2 s-3'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_tke
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'DRAG'
  tzsource%clongname  = 'drag force'
  tzsource%lavailable = odragtree
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'DRAGB'
  tzsource%clongname  = 'drag force due to buildings'
  tzsource%lavailable = ldragbldg
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'DP'
  tzsource%clongname  = 'dynamic production'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'TP'
  tzsource%clongname  = 'thermal production'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'DISS'
  tzsource%clongname  = 'dissipation of TKE'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'TR'
  tzsource%clongname  = 'turbulent transport'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_TKE), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_TKE) )

  call Sourcelist_scan( tbudgets(NBUDGET_TKE), cbulist_rtke )
end if

! Budget of RRV
tbudgets(NBUDGET_RV)%lenabled = lbu_rrv .and. krr >= 1

if ( tbudgets(NBUDGET_RV)%lenabled ) then
  tbudgets(NBUDGET_RV)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_RV)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_RV)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_RV)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_RV)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of water vapor mixing ratio'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'kg kg-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 's-1'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = '2DADV'
  tzsource%clongname  = 'advective forcing'
  tzsource%lavailable = l2d_adv_frc
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = '2DREL'
  tzsource%clongname  = 'relaxation forcing'
  tzsource%lavailable = l2d_rel_frc
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'NUD'
  tzsource%clongname  = 'nudging'
  tzsource%lavailable = onudging
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_rv
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'DCONV'
  tzsource%clongname  = 'KAFR convection'
  tzsource%lavailable = hdconv == 'KAFR' .OR. hsconv == 'KAFR'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'BLAZE'
  tzsource%clongname  = 'blaze fire model contribution'
  tzsource%lavailable = lblaze
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'VTURB'
  tzsource%clongname  = 'vertical turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'HTURB'
  tzsource%clongname  = 'horizontal turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL' .and. HTURBDIM == '3DIM'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'NETUR'
  tzsource%clongname  = 'negativity correction induced by turbulence'
  tzsource%lavailable = hturb == 'TKEL' .and. (      hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                                                .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' )
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'MAFL'
  tzsource%clongname  = 'mass flux'
  tzsource%lavailable = hsconv == 'EDKF'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'SNSUB'
  tzsource%clongname  = 'blowing snow sublimation'
  tzsource%lavailable = lblowsnow .and. lsnowsubl
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_r
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'NEADV'
  tzsource%clongname  = 'negativity correction induced by advection'
  tzsource%lavailable =       hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                         .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'NEGA'
  tzsource%clongname  = 'negativity correction'
  tzsource%lavailable =       hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                         .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'HENU'
  tzsource%clongname  = 'heterogeneous nucleation'
  gtmp = cactccn == 'ABRK' .and. (lorilam .or. ldust .or. lsalt )
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. nmom_c.ge.1 .and. lacti_lima .and. nmod_ccn >= 1   &
                                                     .and. ( .not.lptsplit .or. .not.lsubg_cond )          ) &
                        .or. ( hcloud      == 'C2R2' .and. ( gtmp .or. .not.lsupsat_c2r2 ) )                 &
                        .or. ( hcloud      == 'KHKO' .and. ( gtmp .or. .not.lsupsat_c2r2 ) )
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'REVA'
  tzsource%clongname  = 'rain evaporation'
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. ( ( .not. lptsplit .and. nmom_c.ge.1 .and. nmom_r.ge.1 ) &
                                                             .or.    lptsplit ) )                                 &
                        .or. ( hcloud(1:3) == 'ICE'  .and. lwarm_ice )                                            &
                        .or. ( hcloud      == 'C2R2' .and. lrain_c2r2 )                                           &
                        .or. ( hcloud      == 'KHKO' .and. lrain_c2r2 )                                           &
                        .or.   hcloud      == 'KESS'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'HIN'
  tzsource%clongname  = 'heterogeneous ice nucleation'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .or. ( hcloud == 'LIMA' .and. nmom_i == 1 )
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'HIND'
  tzsource%clongname  = 'heterogeneous nucleation by deposition'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'HONH'
  tzsource%clongname  = 'haze homogeneous nucleation'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima .and. lhhoni_lima .and. nmod_ccn >= 1
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'DEPS'
  tzsource%clongname  = 'deposition on snow'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'DEPG'
  tzsource%clongname  = 'deposition on graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'DEPH'
  tzsource%clongname  = 'deposition on HAIL'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. nmom_h.ge.1 )  &
                        .or. hcloud == 'ICE4' )
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'CEDS'
  tzsource%clongname  = 'adjustment to saturation'
  tzsource%lavailable = hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'ADJU'
  tzsource%clongname  = 'adjustment to saturation'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. ladj_before .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'COND'
  tzsource%clongname  = 'vapor condensation or cloud water evaporation'
  tzsource%lavailable = hcloud == 'C2R2' .or. hcloud == 'KHKO' .or. hcloud == 'KESS' .or. hcloud == 'REVE'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'CORR'
  tzsource%clongname  = 'correction'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'DEPI'
  tzsource%clongname  = 'deposition on ice'
  tzsource%lavailable =      ( hcloud(1:3) == 'ICE' .and. ( .not. lred .or. ( lred .and. ladj_after ) .or. celec /= 'NONE') ) &
                        .or. ( hcloud == 'LIMA' .and. lptsplit )
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'CORR2'
  tzsource%clongname  = 'supplementary correction inside LIMA splitting'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )

  tzsource%cmnhname   = 'NECON'
  tzsource%clongname  = 'negativity correction induced by condensation'
  tzsource%lavailable = (      hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4'   &
                          .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) &
                        .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RV), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_RV) )

  call Sourcelist_scan( tbudgets(NBUDGET_RV), cbulist_rrv )
end if

! Budget of RRC
tbudgets(NBUDGET_RC)%lenabled = lbu_rrc .and. krr >= 2

if ( tbudgets(NBUDGET_RC)%lenabled ) then
  if ( hcloud(1:3) == 'ICE' .and. lred .and. lsedic_ice .and. ldeposc_ice ) &
    call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget', 'lred=T + lsedic=T + ldeposc=T:'// &
                                                        'DEPO and SEDI source terms are mixed and stored in SEDI' )

  tbudgets(NBUDGET_RC)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_RC)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_RC)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_RC)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_RC)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of cloud water mixing ratio'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'kg kg-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 's-1'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_rc
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'DCONV'
  tzsource%clongname  = 'KAFR convection'
  tzsource%lavailable = hdconv == 'KAFR' .OR. hsconv == 'KAFR'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'DEPOTR'
  tzsource%clongname  = 'tree droplet deposition'
  tzsource%lavailable = odragtree .and. odepotree
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'VTURB'
  tzsource%clongname  = 'vertical turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'HTURB'
  tzsource%clongname  = 'horizontal turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL' .and. HTURBDIM == '3DIM'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'NETUR'
  tzsource%clongname  = 'negativity correction induced by turbulence'
  tzsource%lavailable = hturb == 'TKEL' .and. (      hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                                                .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_r
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'NEADV'
  tzsource%clongname  = 'negativity correction induced by advection'
  tzsource%lavailable =       hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                         .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'NEGA'
  tzsource%clongname  = 'negativity correction'
  tzsource%lavailable =       hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                         .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'CORR'
  tzsource%clongname  = 'correction'
!   tzsource%lavailable =       ( hcloud      == 'LIMA' .and. lptsplit .and. nmom_c.ge.1 .and. nmom_r.ge.1 ) &
!                          .or. ( hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE' )
  tzsource%lavailable =  hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'SEDI'
  tzsource%clongname  = 'sedimentation of cloud'
  tzsource%lavailable =    ( hcloud      == 'LIMA' .and. nmom_c.ge.1 .and. lsedc_lima ) &
                      .or. ( hcloud(1:3) == 'ICE'  .and. lsedic_ice )                  &
                      .or. ( hcloud      == 'C2R2' .and. lsedc_c2r2 )                  &
                      .or. ( hcloud      == 'KHKO' .and. lsedc_c2r2 )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'DEPO'
  tzsource%clongname  = 'surface droplet deposition'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. nmom_c.ge.1 .and. ldepoc_lima ) &
                        .or. ( hcloud      == 'C2R2' .and. ldepoc_c2r2 )             &
                        .or. ( hcloud      == 'KHKO' .and. ldepoc_c2r2 )             &
                        .or. ( hcloud(1:3) == 'ICE'  .and. ldeposc_ice .and. celec == 'NONE' )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'R2C1'
  tzsource%clongname  = 'rain to cloud change after sedimentation'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit .and. nmom_c.ge.1 .and. nmom_r.ge.1
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'HENU'
  tzsource%clongname  = 'CCN activation'
  gtmp = cactccn == 'ABRK' .and. (lorilam .or. ldust .or. lsalt )
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. nmom_c.ge.1 .and. lacti_lima .and. nmod_ccn >= 1   &
                                                     .and. ( .not.lptsplit .or. .not.lsubg_cond )          ) &
                        .or. ( hcloud      == 'C2R2' .and. ( gtmp .or. .not.lsupsat_c2r2 ) )                 &
                        .or. ( hcloud      == 'KHKO' .and. ( gtmp .or. .not.lsupsat_c2r2 ) )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'HINC'
  tzsource%clongname  = 'heterogeneous nucleation by contact'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'ADJU'
  tzsource%clongname  = 'adjustment to saturation'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. ladj_before .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'HON'
  tzsource%clongname  = 'homogeneous nucleation'
  tzsource%lavailable = hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'AUTO'
  tzsource%clongname  = 'autoconversion into rain'
  tzsource%lavailable =       ( hcloud      == 'LIMA' .and. ( lptsplit .or. ( nmom_c.ge.1 .and. nmom_r.ge.1 ) ) ) &
                         .or.   hcloud      == 'KESS'                                                           &
                         .or. ( hcloud(1:3) == 'ICE'  .and. lwarm_ice )                                         &
                         .or. ( hcloud      == 'C2R2' .and. lrain_c2r2 )                                        &
                         .or. ( hcloud      == 'KHKO' .and. lrain_c2r2 )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'ACCR'
  tzsource%clongname  = 'accretion of cloud droplets'
  tzsource%lavailable =       ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_c.ge.1 .and. nmom_r.ge.1 ) ) ) &
                         .or.   hcloud      == 'KESS'                                                      &
                         .or. ( hcloud(1:3) == 'ICE'  .and. lwarm_ice )                                    &
                         .or. ( hcloud      == 'C2R2' .and. lrain_c2r2 )                                   &
                         .or. ( hcloud      == 'KHKO' .and. lrain_c2r2 )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'REVA'
  tzsource%clongname  = 'rain evaporation'
  tzsource%lavailable =  hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_c.ge.1 .and. nmom_r.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'HONC'
  tzsource%clongname  = 'droplet homogeneous freezing'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. lnucl_lima ) )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'IMLT'
  tzsource%clongname  = 'melting of ice'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'BERFI'
  tzsource%clongname  = 'Bergeron-Findeisen'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'RIM'
  tzsource%clongname  = 'riming of cloud water'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'CMEL'
  tzsource%clongname  = 'collection by snow and conversion into rain with T>XTT on ice'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'WETG'
  tzsource%clongname  = 'wet growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'DRYG'
  tzsource%clongname  = 'dry growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'CVRC'
  tzsource%clongname  = 'rain to cloud change after other microphysical processes'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'WETH'
  tzsource%clongname  = 'wet growth of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1  &
                                                .and. nmom_c.ge.1    .and. nmom_s.ge.1 )                 &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'DRYH'
  tzsource%clongname  = 'dry growth of hail'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'CEDS'
  tzsource%clongname  = 'adjustment to saturation'
  tzsource%lavailable = hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'DEPI'
  tzsource%clongname  = 'condensation/deposition on ice'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. ( .not. lred .or. ( lred .and. ladj_after ) .or. celec /= 'NONE' )
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'COND'
  tzsource%clongname  = 'vapor condensation or cloud water evaporation'
  tzsource%lavailable = hcloud == 'C2R2' .or. hcloud == 'KHKO' .or. hcloud == 'KESS' .or. hcloud == 'REVE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'CORR2'
  tzsource%clongname  = 'supplementary correction inside LIMA splitting'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )

  tzsource%cmnhname   = 'NECON'
  tzsource%clongname  = 'negativity correction induced by condensation'
  tzsource%lavailable = (      hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4'   &
                          .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) &
                        .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RC), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_RC) )

  call Sourcelist_scan( tbudgets(NBUDGET_RC), cbulist_rrc )
end if

! Budget of RRR
tbudgets(NBUDGET_RR)%lenabled = lbu_rrr .and. krr >= 3

if ( tbudgets(NBUDGET_RR)%lenabled ) then
  tbudgets(NBUDGET_RR)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_RR)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_RR)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_RR)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_RR)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of rain water mixing ratio'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'kg kg-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 's-1'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_rr
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'NETUR'
  tzsource%clongname  = 'negativity correction induced by turbulence'
  tzsource%lavailable = hturb == 'TKEL' .and. ( hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' )
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_r
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'NEADV'
  tzsource%clongname  = 'negativity correction induced by advection'
  tzsource%lavailable =       hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                         .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'NEGA'
  tzsource%clongname  = 'negativity correction'
  tzsource%lavailable =       hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4' &
                         .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'CORR'
  tzsource%clongname  = 'correction'
!   tzsource%lavailable =       ( hcloud      == 'LIMA' .and. lptsplit .and. nmom_c.ge.1 .and. nmom_r.ge.1 ) &
!                          .or. ( hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE' )
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'SEDI'
  tzsource%clongname  = 'sedimentation of rain drops'
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. nmom_c.ge.1 .and. nmom_r.ge.1 ) &
                        .or.   hcloud      == 'KESS'                                     &
                        .or.   hcloud(1:3) == 'ICE'                                      &
                        .or.   hcloud      == 'C2R2'                                     &
                        .or.   hcloud      == 'KHKO'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'R2C1'
  tzsource%clongname  = 'rain to cloud change after sedimentation'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit .and. nmom_c.ge.1 .and. nmom_r.ge.1
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'AUTO'
  tzsource%clongname  = 'autoconversion into rain'
  tzsource%lavailable =       ( hcloud      == 'LIMA' .and. ( lptsplit .or. ( nmom_c.ge.1 .and. nmom_r.ge.1 ) ) ) &
                         .or.   hcloud      == 'KESS'                                                           &
                         .or. ( hcloud(1:3) == 'ICE'  .and. lwarm_ice )                                         &
                         .or. ( hcloud      == 'C2R2' .and. lrain_c2r2 )                                        &
                         .or. ( hcloud      == 'KHKO' .and. lrain_c2r2 )
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'ACCR'
  tzsource%clongname  = 'accretion of cloud droplets'
  tzsource%lavailable =       ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_c.ge.1 .and. nmom_r.ge.1 ) ) ) &
                         .or.   hcloud      == 'KESS'                                                      &
                         .or. ( hcloud(1:3) == 'ICE'  .and. lwarm_ice )                                    &
                         .or. ( hcloud      == 'C2R2' .and. lrain_c2r2 )                                   &
                         .or. ( hcloud      == 'KHKO' .and. lrain_c2r2 )
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'REVA'
  tzsource%clongname  = 'rain evaporation'
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. ( lptsplit .or. ( nmom_c.ge.1 .and. nmom_r.ge.1 ) ) ) &
                        .or.   hcloud      == 'KESS'                                                           &
                        .or. ( hcloud(1:3) == 'ICE'  .and. lwarm_ice )                                         &
                        .or. ( hcloud      == 'C2R2' .and. lrain_c2r2 )                                        &
                        .or. ( hcloud      == 'KHKO' .and. lrain_c2r2 )
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'HONR'
  tzsource%clongname  = 'rain homogeneous freezing'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. lnucl_lima .and. nmom_r.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )


  tzsource%cmnhname   = 'ACC'
  tzsource%clongname  = 'accretion of rain on aggregates'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. ( lptsplit .or. (       nmom_i.ge.1 .and. nmom_c.ge.1      &
                                                                        .and. nmom_s.ge.1 .and. nmom_r.ge.1) ) ) &
                        .or.   hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'CMEL'
  tzsource%clongname  = 'collection of droplets by snow and conversion into rain'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'CFRZ'
  tzsource%clongname  = 'conversion freezing of rain'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'WETG'
  tzsource%clongname  = 'wet growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'DRYG'
  tzsource%clongname  = 'dry growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'GMLT'
  tzsource%clongname  = 'graupel melting'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'CVRC'
  tzsource%clongname  = 'rain to cloud change after other microphysical processes'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'WETH'
  tzsource%clongname  = 'wet growth of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1 &
                                                .and. nmom_c.ge.1    .and. nmom_s.ge.1 )                &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'DRYH'
  tzsource%clongname  = 'dry growth of hail'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'HMLT'
  tzsource%clongname  = 'melting of hail'
  tzsource%lavailable =       ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1 &
                                                .and. nmom_c.ge.1    .and. nmom_s.ge.1 )                &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                         .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'SFR'
  tzsource%clongname  = 'spontaneous freezing'
  tzsource%lavailable = hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

!PW: a documenter
  tzsource%cmnhname   = 'CORR2'
  tzsource%clongname  = 'supplementary correction inside LIMA splitting'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )

  tzsource%cmnhname   = 'NECON'
  tzsource%clongname  = 'negativity correction induced by condensation'
  tzsource%lavailable = (      hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4'   &
                          .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) &
                        .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RR), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_RR) )

  call Sourcelist_scan( tbudgets(NBUDGET_RR), cbulist_rrr )
end if

! Budget of RRI
tbudgets(NBUDGET_RI)%lenabled = lbu_rri .and. krr >= 4

if ( tbudgets(NBUDGET_RI)%lenabled ) then
  tbudgets(NBUDGET_RI)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_RI)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_RI)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_RI)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_RI)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of cloud ice mixing ratio'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'kg kg-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 's-1'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_ri
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'DCONV'
  tzsource%clongname  = 'KAFR convection'
  tzsource%lavailable = hdconv == 'KAFR' .OR. hsconv == 'KAFR'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'VTURB'
  tzsource%clongname  = 'vertical turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HTURB'
  tzsource%clongname  = 'horizontal turbulent diffusion'
  tzsource%lavailable = hturb == 'TKEL' .and. HTURBDIM == '3DIM'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'NETUR'
  tzsource%clongname  = 'negativity correction induced by turbulence'
  tzsource%lavailable = hturb == 'TKEL' .and. ( hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA' )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_r
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'NEADV'
  tzsource%clongname  = 'negativity correction induced by advection'
  tzsource%lavailable =  .true.
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'NEGA'
  tzsource%clongname  = 'negativity correction'
  tzsource%lavailable =  .true.
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'CORR'
  tzsource%clongname  = 'correction'
!   tzsource%lavailable =       ( hcloud      == 'LIMA' .and. lptsplit .and. nmom_i.ge.1 .and. nmom_s.ge.1 ) &
!                          .or. ( hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE' )
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'ADJU'
  tzsource%clongname  = 'adjustment to saturation'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. ladj_before .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'SEDI'
  tzsource%clongname  = 'sedimentation of rain drops'
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. nmom_i.ge.1 .and. lsedi_lima ) &
                        .or.   hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HIN'
  tzsource%clongname  = 'heterogeneous ice nucleation'
  tzsource%lavailable =  hcloud(1:3) == 'ICE' .or. ( hcloud == 'LIMA' .and. nmom_i == 1)
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HIND'
  tzsource%clongname  = 'heterogeneous nucleation by deposition'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HINC'
  tzsource%clongname  = 'heterogeneous nucleation by contact'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HON'
  tzsource%clongname  = 'homogeneous nucleation'
  tzsource%lavailable = hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HONH'
  tzsource%clongname  = 'haze homogeneous nucleation'
  tzsource%lavailable = hcloud == 'LIMA' .and. nmom_i.ge.1 .and. lnucl_lima .and. lhhoni_lima .and. nmod_ccn >= 1
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HONC'
  tzsource%clongname  = 'droplet homogeneous freezing'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. lnucl_lima ) )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'CNVI'
  tzsource%clongname  = 'conversion of snow to cloud ice'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'CNVS'
  tzsource%clongname  = 'conversion of pristine ice to snow'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'AGGS'
  tzsource%clongname  = 'aggregation of snow'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'AUTS'
  tzsource%clongname  = 'autoconversion of ice'
  tzsource%lavailable = hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'IMLT'
  tzsource%clongname  = 'melting of ice'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'BERFI'
  tzsource%clongname  = 'Bergeron-Findeisen'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HMS'
  tzsource%clongname  = 'Hallett-Mossop ice multiplication process due to snow riming'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'CIBU'
  tzsource%clongname  = 'ice multiplication process due to ice collisional breakup'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 .and. lcibu ) )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'CFRZ'
  tzsource%clongname  = 'conversion freezing of rain'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'RDSF'
  tzsource%clongname  = 'ice multiplication process following rain contact freezing'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 .and. lrdsf ) )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'WETG'
  tzsource%clongname  = 'wet growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'DRYG'
  tzsource%clongname  = 'dry growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'HMG'
  tzsource%clongname  = 'Hallett-Mossop ice multiplication process due to graupel riming'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'WETH'
  tzsource%clongname  = 'wet growth of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1 &
                                                .and. nmom_c.ge.1    .and. nmom_s.ge.1 )                &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'DRYH'
  tzsource%clongname  = 'dry growth of hail'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'CEDS'
  tzsource%clongname  = 'adjustment to saturation'
  tzsource%lavailable = hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'DEPI'
  tzsource%clongname  = 'condensation/deposition on ice'
  tzsource%lavailable =      ( hcloud(1:3) == 'ICE' .and. ( .not. lred .or. ( lred .and. ladj_after ) .or. celec /= 'NONE') ) &
                        .or. ( hcloud == 'LIMA' .and. lptsplit )
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'CORR2'
  tzsource%clongname  = 'supplementary correction inside LIMA splitting'
  tzsource%lavailable = hcloud == 'LIMA' .and. lptsplit
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )

  tzsource%cmnhname   = 'NECON'
  tzsource%clongname  = 'negativity correction induced by condensation'
  tzsource%lavailable = celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RI), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_RI) )

  call Sourcelist_scan( tbudgets(NBUDGET_RI), cbulist_rri )
end if

! Budget of RRS
tbudgets(NBUDGET_RS)%lenabled = lbu_rrs .and. krr >= 5

if ( tbudgets(NBUDGET_RS)%lenabled ) then
  tbudgets(NBUDGET_RS)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_RS)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_RS)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_RS)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_RS)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of snow/aggregate mixing ratio'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'kg kg-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 's-1'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_rs
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

!   tzsource%cmnhname   = 'NETUR'
!   tzsource%clongname  = 'negativity correction induced by turbulence'
!   tzsource%lavailable = hturb == 'TKEL' .and. (      hcloud == 'ICE3' .or. hcloud == 'ICE4' &
!                                   .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' )
!   call Budget_source_add( tbudgets(NBUDGET_RS), tzsource nneturrs )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_r
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'NEADV'
  tzsource%clongname  = 'negativity correction induced by advection'
  tzsource%lavailable =  .true.
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'NEGA'
  tzsource%clongname  = 'negativity correction'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'CORR'
  tzsource%clongname  = 'correction'
!   tzsource%lavailable =       ( hcloud      == 'LIMA' .and. lptsplit .and. nmom_i.ge.1 .and. nmom_s.ge.1 ) &
!                          .or. ( hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE' )
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'SEDI'
  tzsource%clongname  = 'sedimentation'
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. nmom_i.ge.1 .and. nmom_s.ge.1 ) &
                        .or.   hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'CNVI'
  tzsource%clongname  = 'conversion of snow to cloud ice'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'DEPS'
  tzsource%clongname  = 'deposition on snow'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'CNVS'
  tzsource%clongname  = 'conversion of pristine ice to snow'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'AGGS'
  tzsource%clongname  = 'aggregation of snow'
  tzsource%lavailable = ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 ) ) ) .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'AUTS'
  tzsource%clongname  = 'autoconversion of ice'
  tzsource%lavailable = hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'RIM'
  tzsource%clongname  = 'riming of cloud water'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'HMS'
  tzsource%clongname  = 'Hallett-Mossop ice multiplication process due to snow riming'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'CIBU'
  tzsource%clongname  = 'ice multiplication process due to ice collisional breakup'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 .and. lcibu ) )
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'ACC'
  tzsource%clongname  = 'accretion of rain on snow'
  tzsource%lavailable =       ( hcloud == 'LIMA' .and. ( lptsplit .or. (       nmom_i.ge.1 .and. nmom_c.ge.1      &
                                                                         .and. nmom_s.ge.1 .and. nmom_r.ge.1) ) ) &
                         .or.   hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'CMEL'
  tzsource%clongname  = 'conversion melting'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'WETG'
  tzsource%clongname  = 'wet growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'DRYG'
  tzsource%clongname  = 'dry growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'WETH'
  tzsource%clongname  = 'wet growth of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1 &
                                                .and. nmom_c.ge.1    .and. nmom_s.ge.1 )                &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'DRYH'
  tzsource%clongname  = 'dry growth of hail'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )

  tzsource%cmnhname   = 'NECON'
  tzsource%clongname  = 'negativity correction induced by condensation'
  tzsource%lavailable = celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RS), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_RS) )

  call Sourcelist_scan( tbudgets(NBUDGET_RS), cbulist_rrs )
end if

! Budget of RRG
tbudgets(NBUDGET_RG)%lenabled = lbu_rrg .and. krr >= 6

if ( tbudgets(NBUDGET_RG)%lenabled ) then
  tbudgets(NBUDGET_RG)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_RG)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_RG)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_RG)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_RG)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of graupel mixing ratio'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'kg kg-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 's-1'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_rg
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

!   tzsource%cmnhname   = 'NETUR'
!   tzsource%clongname  = 'negativity correction induced by turbulence'
!   tzsource%lavailable = hturb == 'TKEL' .and. (      hcloud == 'ICE3' .or. hcloud == 'ICE4' &
!                                   .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' )
!   call Budget_source_add( tbudgets(NBUDGET_RG), tzsource nneturrg )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_r
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'NEADV'
  tzsource%clongname  = 'negativity correction induced by advection'
  tzsource%lavailable =  hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'NEGA'
  tzsource%clongname  = 'negativity correction'
  tzsource%lavailable =  hcloud == 'ICE3' .or. hcloud == 'ICE4' .or. hcloud == 'LIMA'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'CORR'
  tzsource%clongname  = 'correction'
  tzsource%lavailable = hcloud(1:3) == 'ICE' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'SEDI'
  tzsource%clongname  = 'sedimentation'
  tzsource%lavailable =      ( hcloud      == 'LIMA' .and. nmom_i.ge.1 .and. nmom_s.ge.1 ) &
                        .or.   hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'HONR'
  tzsource%clongname  = 'rain homogeneous freezing'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. lnucl_lima .and. nmom_r.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'SFR'
  tzsource%clongname  = 'spontaneous freezing'
  tzsource%lavailable = hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'DEPG'
  tzsource%clongname  = 'deposition on graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'RIM'
  tzsource%clongname  = 'riming of cloud water'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'ACC'
  tzsource%clongname  = 'accretion of rain on graupel'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. ( lptsplit .or. (       nmom_i.ge.1 .and. nmom_c.ge.1      &
                                                                        .and. nmom_s.ge.1 .and. nmom_r.ge.1) ) ) &
                        .or.   hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'CMEL'
  tzsource%clongname  = 'conversion melting of snow'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'CFRZ'
  tzsource%clongname  = 'conversion freezing of rain'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'RDSF'
  tzsource%clongname  = 'ice multiplication process following rain contact freezing'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 .and. lrdsf ) )
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )  

  tzsource%cmnhname   = 'WETG'
  tzsource%clongname  = 'wet growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'GHCV'
  tzsource%clongname  = 'graupel to hail conversion'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'DRYG'
  tzsource%clongname  = 'dry growth of graupel'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'HMG'
  tzsource%clongname  = 'Hallett-Mossop ice multiplication process due to graupel riming'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) )
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'GMLT'
  tzsource%clongname  = 'graupel melting'
  tzsource%lavailable =    ( hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) ) ) &
                        .or. hcloud(1:3) == 'ICE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'WETH'
  tzsource%clongname  = 'wet growth of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1 &
                                                .and. nmom_c.ge.1    .and. nmom_s.ge.1 )                &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'COHG'
  tzsource%clongname  = 'conversion of hail to graupel'
  tzsource%lavailable = hcloud == 'LIMA' .and. (lptsplit .or. (nmom_h.ge.1 .and. nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1) )
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'HGCV'
  tzsource%clongname  = 'hail to graupel conversion'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'DRYH'
  tzsource%clongname  = 'dry growth of hail'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )

  tzsource%cmnhname   = 'NECON'
  tzsource%clongname  = 'negativity correction induced by condensation'
  tzsource%lavailable = (      hcloud == 'KESS' .or. hcloud == 'ICE3' .or. hcloud == 'ICE4'   &
                          .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' ) &
                        .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RG), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_RG) )

  call Sourcelist_scan( tbudgets(NBUDGET_RG), cbulist_rrg )
end if

! Budget of RRH
tbudgets(NBUDGET_RH)%lenabled = lbu_rrh .and. krr >= 7

if ( tbudgets(NBUDGET_RH)%lenabled ) then
  tbudgets(NBUDGET_RH)%trhodj => tburhodj

  !Allocate all basic source terms (used or not)
  !The size should be large enough (bigger than necessary is OK)
  tbudgets(NBUDGET_RH)%nsourcesmax = NSOURCESMAX
  allocate( tbudgets(NBUDGET_RH)%tsources(NSOURCESMAX) )

  allocate( tbudgets(NBUDGET_RH)%xtmpstore(ibudim1, ibudim2, ibudim3) )

  tbudgets(NBUDGET_RH)%tsources(:)%ngroup = 0

  tzsource%ccomment = 'Budget of hail mixing ratio'
  tzsource%ngrid    = 1

  tzsource%cunits   = 'kg kg-1'

  tzsource%cmnhname   = 'INIF'
  tzsource%clongname  = 'initial state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'ENDF'
  tzsource%clongname  = 'final state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource, odonotinit = .true., ooverwrite = .true. )

  tzsource%cmnhname   = 'AVEF'
  tzsource%clongname  = 'averaged state'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource, odonotinit = .true., ooverwrite = .false. )

  tzsource%cunits   = 's-1'

  tzsource%cmnhname   = 'ASSE'
  tzsource%clongname  = 'time filter (Asselin)'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'NEST'
  tzsource%clongname  = 'nesting'
  tzsource%lavailable = nmodel > 1
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'FRC'
  tzsource%clongname  = 'forcing'
  tzsource%lavailable = lforcing
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'DIF'
  tzsource%clongname  = 'numerical diffusion'
  tzsource%lavailable = onumdifth
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'REL'
  tzsource%clongname  = 'relaxation'
  tzsource%lavailable = ohorelax_rh
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

!   tzsource%cmnhname   = 'NETUR'
!   tzsource%clongname  = 'negativity correction induced by turbulence'
!   tzsource%lavailable = hturb == 'TKEL' .and. (      hcloud == 'ICE3' .or. hcloud == 'ICE4' &
!                                   .or. hcloud == 'KHKO' .or. hcloud == 'C2R2' .or. hcloud == 'LIMA' )
!   call Budget_source_add( tbudgets(NBUDGET_RH), tzsource nneturrh )

  tzsource%cmnhname   = 'VISC'
  tzsource%clongname  = 'viscosity'
  tzsource%lavailable = lvisc .and. lvisc_r
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'ADV'
  tzsource%clongname  = 'total advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'NEADV'
  tzsource%clongname  = 'negativity correction induced by advection'
  tzsource%lavailable = .true.
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'NEGA'
  tzsource%clongname  = 'negativity correction'
  tzsource%lavailable =  .true.
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'SEDI'
  tzsource%clongname  = 'sedimentation'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. nmom_i.ge.1 .and. nmom_h.ge.1 ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'DEPH'
  tzsource%clongname  = 'deposition on hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. nmom_h.ge.1 ) 
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )


  tzsource%cmnhname   = 'GHCV'
  tzsource%clongname  = 'graupel to hail conversion'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'WETG'
  tzsource%clongname  = 'wet growth of graupel'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. nmom_h.ge.1                                                            &
                                                .and. ( lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) ) ) &
                        .or. ( hcloud == 'ICE4' .and. ( .not. lred .or. celec /= 'NONE' ) )
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'WETH'
  tzsource%clongname  = 'wet growth of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not.lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1    &
                                                                    .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'COHG'
  tzsource%clongname  = 'conversion from hail to graupel'
  tzsource%lavailable = hcloud == 'LIMA' .and. ( lptsplit .or. (nmom_h.ge.1 .and. nmom_i.ge.1 &
                                                             .and. nmom_c.ge.1 .and. nmom_s.ge.1) )
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'HGCV'
  tzsource%clongname  = 'hail to graupel conversion'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'DRYH'
  tzsource%clongname  = 'dry growth of hail'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'HMLT'
  tzsource%clongname  = 'melting of hail'
  tzsource%lavailable =      ( hcloud == 'LIMA' .and. .not. lptsplit .and. nmom_h.ge.1 .and. nmom_i.ge.1   &
                                                                     .and. nmom_c.ge.1 .and. nmom_s.ge.1 ) &
                        .or. ( hcloud == 'LIMA' .and. lptsplit ) &
                        .or.   hcloud == 'ICE4'
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'CORR'
  tzsource%clongname  = 'correction'
  tzsource%lavailable = hcloud == 'ICE4' .and. lred .and. celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )

  tzsource%cmnhname   = 'NECON'
  tzsource%clongname  = 'negativity correction induced by condensation'
  tzsource%lavailable = celec == 'NONE'
  call Budget_source_add( tbudgets(NBUDGET_RH), tzsource )


  call Sourcelist_sort_compact( tbudgets(NBUDGET_RH) )

  call Sourcelist_scan( tbudgets(NBUDGET_RH), cbulist_rrh )
end if

! Budgets of RSV (scalar variables)

if ( ksv > 999 ) call Print_msg( NVERB_FATAL, 'BUD', 'Ini_budget', 'number of scalar variables > 999' )

SV_BUDGETS: do jsv = 1, ksv
  ibudget = NBUDGET_SV1 - 1 + jsv

  tbudgets(ibudget)%lenabled = lbu_rsv

  if ( lbu_rsv ) then
    tbudgets(ibudget)%trhodj => tburhodj

    !Allocate all basic source terms (used or not)
    !The size should be large enough (bigger than necessary is OK)
    tbudgets(ibudget)%nsourcesmax = NSOURCESMAX
    allocate( tbudgets(ibudget)%tsources(NSOURCESMAX) )

    allocate( tbudgets(ibudget)%xtmpstore(ibudim1, ibudim2, ibudim3) )

    tbudgets(ibudget)%tsources(:)%ngroup = 0

    tzsource%ccomment = 'Budget of scalar variable ' // tsvlist(jsv)%cmnhname
    tzsource%ngrid    = 1

    tzsource%cunits   = '1'

    tzsource%cmnhname   = 'INIF'
    tzsource%clongname  = 'initial state'
    tzsource%lavailable = .true.
    call Budget_source_add( tbudgets(ibudget), tzsource, odonotinit = .true., ooverwrite = .true. )

    tzsource%cmnhname   = 'ENDF'
    tzsource%clongname  = 'final state'
    tzsource%lavailable = .true.
    call Budget_source_add( tbudgets(ibudget), tzsource, odonotinit = .true., ooverwrite = .true. )

    tzsource%cmnhname   = 'AVEF'
    tzsource%clongname  = 'averaged state'
    tzsource%lavailable = .true.
    call Budget_source_add( tbudgets(ibudget), tzsource, odonotinit = .true., ooverwrite = .false. )

    tzsource%cunits   = 's-1'

    tzsource%cmnhname   = 'ASSE'
    tzsource%clongname  = 'time filter (Asselin)'
    tzsource%lavailable = .true.
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'NEST'
    tzsource%clongname  = 'nesting'
    tzsource%lavailable = nmodel > 1
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'FRC'
    tzsource%clongname  = 'forcing'
    tzsource%lavailable = lforcing
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'DIF'
    tzsource%clongname  = 'numerical diffusion'
    tzsource%lavailable = onumdifsv
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'REL'
    tzsource%clongname  = 'relaxation'
    tzsource%lavailable = ohorelax_sv( jsv ) .or. ( celec /= 'NONE' .and. lrelax2fw_ion                 &
                                                    .and. (jsv == nsv_elecbeg .or. jsv == nsv_elecend ) )
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'DCONV'
    tzsource%clongname  = 'KAFR convection'
    tzsource%lavailable = ( hdconv == 'KAFR' .or. hsconv == 'KAFR' ) .and. ochtrans
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'VTURB'
    tzsource%clongname  = 'vertical turbulent diffusion'
    tzsource%lavailable = hturb == 'TKEL'
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'HTURB'
    tzsource%clongname  = 'horizontal turbulent diffusion'
    tzsource%lavailable = hturb == 'TKEL' .and. HTURBDIM == '3DIM'
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'MAFL'
    tzsource%clongname  = 'mass flux'
    tzsource%lavailable = hsconv == 'EDKF'
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'VISC'
    tzsource%clongname  = 'viscosity'
    tzsource%lavailable = lvisc .and. lvisc_sv
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'ADV'
    tzsource%clongname  = 'total advection'
    tzsource%lavailable = .true.
    call Budget_source_add( tbudgets(ibudget), tzsource )

    tzsource%cmnhname   = 'NEGA2'
    tzsource%clongname  = 'negativity correction'
    tzsource%lavailable = .true.
    call Budget_source_add( tbudgets(ibudget), tzsource )

    ! Add specific source terms to different scalar variables
    SV_VAR: if ( jsv <= nsv_user ) then
      ! nsv_user case
      ! Nothing to do

    else if ( jsv >= nsv_c2r2beg .and. jsv <= nsv_c2r2end ) then SV_VAR
      ! C2R2 or KHKO Case

      ! Source terms in common for all C2R2/KHKO budgets
      tzsource%cmnhname   = 'NETUR'
      tzsource%clongname  = 'negativity correction induced by turbulence'
      tzsource%lavailable = hturb == 'TKEL'
      call Budget_source_add( tbudgets(ibudget), tzsource )

      tzsource%cmnhname   = 'NEADV'
      tzsource%clongname  = 'negativity correction induced by advection'
      tzsource%lavailable = .true.
      call Budget_source_add( tbudgets(ibudget), tzsource )

      tzsource%cmnhname   = 'NEGA'
      tzsource%clongname  = 'negativity correction'
      tzsource%lavailable = .true.
      call Budget_source_add( tbudgets(ibudget), tzsource )

      tzsource%cmnhname   = 'NECON'
      tzsource%clongname  = 'negativity correction induced by condensation'
      tzsource%lavailable = .true.
      call Budget_source_add( tbudgets(ibudget), tzsource )

      ! Source terms specific to each budget
      SV_C2R2: select case( jsv - nsv_c2r2beg + 1 )
        case ( 1 ) SV_C2R2
          ! Concentration of activated nuclei
          tzsource%cmnhname   = 'HENU'
          tzsource%clongname  = 'CCN activation'
          gtmp = cactccn == 'ABRK' .and. (lorilam .or. ldust .or. lsalt )
          tzsource%lavailable =  gtmp .or. ( .not.gtmp .and. .not.lsupsat_c2r2 )
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'CEVA'
          tzsource%clongname  = 'evaporation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )


        case ( 2 ) SV_C2R2
          ! Concentration of cloud droplets
          tzsource%cmnhname   = 'DEPOTR'
          tzsource%clongname  = 'tree droplet deposition'
          tzsource%lavailable = odragtree .and. odepotree
          call Budget_source_add( tbudgets(ibudget), tzsource)

          tzsource%cmnhname   = 'HENU'
          tzsource%clongname  = 'CCN activation'
          gtmp = cactccn == 'ABRK' .and. (lorilam .or. ldust .or. lsalt )
          tzsource%lavailable =  gtmp .or. ( .not.gtmp .and. .not.lsupsat_c2r2 )
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SELF'
          tzsource%clongname  = 'self-collection of cloud droplets'
          tzsource%lavailable = lrain_c2r2
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'ACCR'
          tzsource%clongname  = 'accretion of cloud droplets'
          tzsource%lavailable = lrain_c2r2
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SEDI'
          tzsource%clongname  = 'sedimentation'
          tzsource%lavailable = lsedc_c2r2
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DEPO'
          tzsource%clongname  = 'surface droplet deposition'
          tzsource%lavailable = ldepoc_c2r2
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'CEVA'
          tzsource%clongname  = 'evaporation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )


        case ( 3 ) SV_C2R2
          ! Concentration of raindrops
          tzsource%cmnhname   = 'AUTO'
          tzsource%clongname  = 'autoconversion into rain'
          tzsource%lavailable = lrain_c2r2
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SCBU'
          tzsource%clongname  = 'self collection - coalescence/break-up'
          tzsource%lavailable = hcloud /= 'KHKO'
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'REVA'
          tzsource%clongname  = 'rain evaporation'
          tzsource%lavailable = lrain_c2r2
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'BRKU'
          tzsource%clongname  = 'spontaneous break-up'
          tzsource%lavailable = lrain_c2r2
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SEDI'
          tzsource%clongname  = 'sedimentation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )


        case ( 4 ) SV_C2R2
          ! Supersaturation
          tzsource%cmnhname   = 'CEVA'
          tzsource%clongname  = 'evaporation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

      end select SV_C2R2


    else if ( jsv >= nsv_lima_beg .and. jsv <= nsv_lima_end ) then SV_VAR
      ! LIMA case

      ! Source terms in common for all LIMA budgets (except supersaturation)
      if ( jsv /= nsv_lima_spro ) then
        tzsource%cmnhname   = 'NETUR'
        tzsource%clongname  = 'negativity correction induced by turbulence'
        tzsource%lavailable = hturb == 'TKEL'
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'NEADV'
        tzsource%clongname  = 'negativity correction induced by advection'
        tzsource%lavailable = .true.
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'NEGA'
        tzsource%clongname  = 'negativity correction'
        tzsource%lavailable = .true.
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'NECON'
        tzsource%clongname  = 'negativity correction induced by condensation'
        tzsource%lavailable = .true.
        call Budget_source_add( tbudgets(ibudget), tzsource )
      end if


      ! Source terms specific to each budget
      SV_LIMA: if ( jsv == nsv_lima_nc ) then
        ! Cloud droplets concentration
        tzsource%cmnhname   = 'DEPOTR'
        tzsource%clongname  = 'tree droplet deposition'
        tzsource%lavailable = odragtree .and. odepotree
        call Budget_source_add( tbudgets(ibudget), tzsource )

!         tzsource%cmnhname   = 'CORR'
!         tzsource%clongname  = 'correction'
!         tzsource%lavailable = lptsplit .and. nmom_c.ge.1  .and. nmom_r.ge.1
!         call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'SEDI'
        tzsource%clongname  = 'sedimentation'
        tzsource%lavailable = nmom_c.ge.1  .and. lsedc_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'DEPO'
        tzsource%clongname  = 'surface droplet deposition'
        tzsource%lavailable = nmom_c.ge.1  .and. ldepoc_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'R2C1'
        tzsource%clongname  = 'rain to cloud change after sedimentation'
        tzsource%lavailable = lptsplit .and. nmom_c.ge.1  .and. nmom_r.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HENU'
        tzsource%clongname  = 'CCN activation'
        tzsource%lavailable = nmom_c.ge.1  .and. lacti_lima .and. nmod_ccn >= 1 .and. ( .not.lptsplit .or. .not.lsubg_cond )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HINC'
        tzsource%clongname  = 'heterogeneous nucleation by contact'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'SELF'
        tzsource%clongname  = 'self-collection of cloud droplets'
        tzsource%lavailable = lptsplit .or. (nmom_c.ge.1 .and. nmom_r.ge.1)
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'AUTO'
        tzsource%clongname  = 'autoconversion into rain'
        tzsource%lavailable = lptsplit .or. ( nmom_c.ge.1  .and. nmom_r.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'ACCR'
        tzsource%clongname  = 'accretion of cloud droplets'
        tzsource%lavailable = lptsplit .or. ( nmom_c.ge.1  .and. nmom_r.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'REVA'
        tzsource%clongname  = 'rain evaporation'
        tzsource%lavailable = lptsplit .or. ( nmom_c.ge.1  .and. nmom_r.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HONC'
        tzsource%clongname  = 'droplet homogeneous freezing'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. lnucl_lima )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'IMLT'
        tzsource%clongname  = 'melting of ice'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'RIM'
        tzsource%clongname  = 'riming of cloud water'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETG'
        tzsource%clongname  = 'wet growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'DRYG'
        tzsource%clongname  = 'dry growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CVRC'
        tzsource%clongname  = 'rain to cloud change after other microphysical processes'
        tzsource%lavailable = lptsplit
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETH'
        tzsource%clongname  = 'wet growth of hail'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CORR2'
        tzsource%clongname  = 'supplementary correction inside LIMA splitting'
        tzsource%lavailable = lptsplit
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_c.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv == nsv_lima_nr ) then SV_LIMA
        ! Rain drops concentration
!         tzsource%cmnhname   = 'CORR'
!         tzsource%clongname  = 'correction'
!         tzsource%lavailable = lptsplit .and. nmom_c.ge.1  .and. nmom_r.ge.1
!         call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'SEDI'
        tzsource%clongname  = 'sedimentation'
        tzsource%lavailable = nmom_c.ge.1  .and. nmom_r.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'R2C1'
        tzsource%clongname  = 'rain to cloud change after sedimentation'
        tzsource%lavailable = lptsplit .and. nmom_c.ge.1  .and. nmom_r.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'AUTO'
        tzsource%clongname  = 'autoconversion into rain'
        tzsource%lavailable = lptsplit .or. (nmom_c.ge.1  .and. nmom_r.ge.1)
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'SCBU'
        tzsource%clongname  = 'self collection - coalescence/break-up'
        tzsource%lavailable = lptsplit .or. (nmom_c.ge.1  .and. nmom_r.ge.1)
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'REVA'
        tzsource%clongname  = 'rain evaporation'
        tzsource%lavailable = lptsplit .or. (nmom_c.ge.1  .and. nmom_r.ge.1)
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'BRKU'
        tzsource%clongname  = 'spontaneous break-up'
        tzsource%lavailable = lptsplit .or. (nmom_c.ge.1  .and. nmom_r.ge.1)
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HONR'
        tzsource%clongname  = 'rain homogeneous freezing'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_r.ge.1 .and. lnucl_lima )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'ACC'
        tzsource%clongname  = 'accretion of rain  on aggregates'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 .and. nmom_r.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CFRZ'
        tzsource%clongname  = 'conversion freezing of rain'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETG'
        tzsource%clongname  = 'wet growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'DRYG'
        tzsource%clongname  = 'dry growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'GMLT'
        tzsource%clongname  = 'graupel melting'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CVRC'
        tzsource%clongname  = 'rain to cloud change after other microphysical processes'
        tzsource%lavailable = lptsplit
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETH'
        tzsource%clongname  = 'wet growth of hail'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HMLT'
        tzsource%clongname  = 'melting of hail'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CORR2'
        tzsource%clongname  = 'supplementary correction inside LIMA splitting'
        tzsource%lavailable = lptsplit
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv >= nsv_lima_ccn_free .and. jsv <= nsv_lima_ccn_free + nmod_ccn - 1 ) then SV_LIMA
        ! Free CCN concentration
        tzsource%cmnhname   = 'HENU'
        tzsource%clongname  = 'CCN activation'
        tzsource%lavailable = nmom_c.ge.1  .and. lacti_lima .and. nmod_ccn >= 1 .and. ( .not.lptsplit .or. .not.lsubg_cond )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HONH'
        tzsource%clongname  = 'haze homogeneous nucleation'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima .and. lhhoni_lima .and. nmod_ccn >= 1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_c.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'SCAV'
        tzsource%clongname  = 'scavenging'
        tzsource%lavailable = lscav_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv >= nsv_lima_ccn_acti .and. jsv <= nsv_lima_ccn_acti + nmod_ccn - 1 ) then SV_LIMA
        ! Activated CCN concentration
        tzsource%cmnhname   = 'HENU'
        tzsource%clongname  = 'CCN activation'
        tzsource%lavailable = nmom_c.ge.1  .and. lacti_lima .and. nmod_ccn >= 1 .and. ( .not.lptsplit .or. .not.lsubg_cond )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HINC'
        tzsource%clongname  = 'heterogeneous nucleation by contact'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima .and. .not. lmeyers_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_c.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv == nsv_lima_scavmass ) then SV_LIMA
        ! Scavenged mass variable
        tzsource%cmnhname   = 'SCAV'
        tzsource%clongname  = 'scavenging'
        tzsource%lavailable = lscav_lima .and. laero_mass_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = lscav_lima .and. laero_mass_lima .and. .not.lspro_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv == nsv_lima_ni ) then SV_LIMA
        ! Pristine ice crystals concentration
!         tzsource%cmnhname   = 'CORR'
!         tzsource%clongname  = 'correction'
!         tzsource%lavailable = lptsplit .and. nmom_i.ge.1 .and. nmom_s.ge.1
!         call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'SEDI'
        tzsource%clongname  = 'sedimentation'
        tzsource%lavailable = nmom_i.ge.1 .and. lsedi_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HIND'
        tzsource%clongname  = 'heterogeneous nucleation by deposition'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HINC'
        tzsource%clongname  = 'heterogeneous nucleation by contact'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HONH'
        tzsource%clongname  = 'haze homogeneous nucleation'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima .and. lhhoni_lima .and. nmod_ccn >= 1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HONC'
        tzsource%clongname  = 'droplet homogeneous freezing'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. lnucl_lima )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CNVI'
        tzsource%clongname  = 'conversion of snow to cloud ice'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CNVS'
        tzsource%clongname  = 'conversion of pristine ice to snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'AGGS'
        tzsource%clongname  = 'aggregation of snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'IMLT'
        tzsource%clongname  = 'melting of ice'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HMS'
        tzsource%clongname  = 'Hallett-Mossop ice multiplication process due to snow riming'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CIBU'
        tzsource%clongname  = 'ice multiplication process due to ice collisional breakup'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 .and. lcibu )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CFRZ'
        tzsource%clongname  = 'conversion freezing of rain'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'RDSF'
        tzsource%clongname  = 'ice multiplication process following rain contact freezing'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 .and. lrdsf )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETG'
        tzsource%clongname  = 'wet growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'DRYG'
        tzsource%clongname  = 'dry growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HMG'
        tzsource%clongname  = 'Hallett-Mossop ice multiplication process due to graupel riming'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETH'
        tzsource%clongname  = 'wet growth of hail'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_i.ge.1 .and. .not.lptsplit .and. .not.lspro_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CORR2'
        tzsource%clongname  = 'supplementary correction inside LIMA splitting'
        tzsource%lavailable = lptsplit
        call Budget_source_add( tbudgets(ibudget), tzsource )

      else if ( jsv == nsv_lima_ns ) then SV_LIMA

        tzsource%cmnhname   = 'SEDI'
        tzsource%clongname  = 'sedimentation'
        tzsource%lavailable = lptsplit .or. ( nmom_s.ge.2 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CNVI'
        tzsource%clongname  = 'conversion of snow to cloud ice'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CNVS'
        tzsource%clongname  = 'conversion of pristine ice to snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'BRKU'
        tzsource%clongname  = 'break up of snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'RIM'
        tzsource%clongname  = 'heavy riming of cloud droplet on snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'ACC'
        tzsource%clongname  = 'accretion of rain on snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CMEL'
        tzsource%clongname  = 'conversion melting of snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETG'
        tzsource%clongname  = 'wet growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'DRYG'
        tzsource%clongname  = 'dry growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETH'
        tzsource%clongname  = 'wet growth of hail'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'SSC'
        tzsource%clongname  = 'snow self collection'
        tzsource%lavailable = lptsplit  .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )
        
        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_i.ge.1 .and. .not.lptsplit .and. .not.lspro_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

      else if ( jsv == nsv_lima_ng ) then SV_LIMA

        tzsource%cmnhname   = 'SEDI'
        tzsource%clongname  = 'sedimentation'
        tzsource%lavailable = nmom_i.ge.1 .or. ( nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'RIM'
        tzsource%clongname  = 'heavy riming of cloud droplet on snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'ACC'
        tzsource%clongname  = 'accretion of rain on snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CMEL'
        tzsource%clongname  = 'conversion melting of snow'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CFRZ'
        tzsource%clongname  = 'conversion freezing of raindrop'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETG'
        tzsource%clongname  = 'wet growth of graupel'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'GMLT'
        tzsource%clongname  = 'graupel melting'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1  .and. nmom_s.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETH'
        tzsource%clongname  = 'wet growth of hail'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'COHG'
        tzsource%clongname  = 'conversion hail graupel'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )
        
        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_i.ge.1 .and. .not.lptsplit .and. .not.lspro_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )   

      else if ( jsv == nsv_lima_nh .and. nmom_h.ge.1) then SV_LIMA

        tzsource%cmnhname   = 'SEDI'
        tzsource%clongname  = 'sedimentation'
        tzsource%lavailable = nmom_i.ge.1 .or. ( nmom_i.ge.1 .and. nmom_s.ge.1 .and. nmom_h.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'WETG'
        tzsource%clongname  = 'wet growth of graupel'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'COHG'
        tzsource%clongname  = 'conversion hail graupel'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HMLT'
        tzsource%clongname  = 'hail melting'
        tzsource%lavailable = lptsplit .or. nmom_h.ge.1
        call Budget_source_add( tbudgets(ibudget), tzsource )
        
      else if ( jsv >= nsv_lima_ifn_free .and. jsv <= nsv_lima_ifn_free + nmod_ifn - 1 ) then SV_LIMA
        ! Free IFN concentration
        tzsource%cmnhname   = 'HIND'
        tzsource%clongname  = 'heterogeneous nucleation by deposition'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima .and. .not. lmeyers_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_i.ge.1 .and. .not.lptsplit .and. .not.lspro_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'SCAV'
        tzsource%clongname  = 'scavenging'
        tzsource%lavailable = lscav_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv >= nsv_lima_ifn_nucl .and. jsv <= nsv_lima_ifn_nucl + nmod_ifn - 1 ) then SV_LIMA
        ! Nucleated IFN concentration
        tzsource%cmnhname   = 'HIND'
        tzsource%clongname  = 'heterogeneous nucleation by deposition'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima                                                     &
                              .and. ( ( lmeyers_lima .and. jsv == nsv_lima_ifn_nucl ) .or. .not. lmeyers_lima )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'HINC'
        tzsource%clongname  = 'heterogeneous nucleation by contact'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima .and. lmeyers_lima .and. jsv == nsv_lima_ifn_nucl
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'IMLT'
        tzsource%clongname  = 'melting of ice'
        tzsource%lavailable = lptsplit .or. ( nmom_i.ge.1 .and. nmom_c.ge.1 )
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_i.ge.1 .and. .not.lptsplit .and. .not.lspro_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv >= nsv_lima_imm_nucl .and. jsv <= nsv_lima_imm_nucl + nmod_imm - 1 ) then SV_LIMA
        ! Nucleated IMM concentration
        tzsource%cmnhname   = 'HINC'
        tzsource%clongname  = 'heterogeneous nucleation by contact'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima .and. .not. lmeyers_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )

        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = nmom_i.ge.1 .and. .not.lptsplit .and. .not.lspro_lima
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv == nsv_lima_hom_haze ) then SV_LIMA
        ! Homogeneous freezing of CCN
        tzsource%cmnhname   = 'HONH'
        tzsource%clongname  = 'haze homogeneous nucleation'
        tzsource%lavailable = nmom_i.ge.1 .and. lnucl_lima .and.                                             &
                              ( ( lhhoni_lima .and. nmod_ccn >= 1 ) .or. ( .not.lptsplit .and. nmom_c.ge.1 ) )
        call Budget_source_add( tbudgets(ibudget), tzsource )


      else if ( jsv == nsv_lima_spro ) then SV_LIMA
        ! Supersaturation
        tzsource%cmnhname   = 'CEDS'
        tzsource%clongname  = 'adjustment to saturation'
        tzsource%lavailable = .true.
        call Budget_source_add( tbudgets(ibudget), tzsource )

      end if SV_LIMA


    else if ( jsv >= nsv_elecbeg .and. jsv <= nsv_elecend ) then SV_VAR
      ! Electricity case
      tzsource%cmnhname   = 'NEGA'
      tzsource%clongname  = 'negativity correction'
      tzsource%lavailable = .true.
      call Budget_source_add( tbudgets(ibudget), tzsource )

      SV_ELEC: select case( jsv - nsv_elecbeg + 1 )
        case ( 1 ) SV_ELEC
          ! volumetric charge of water vapor
          tzsource%cmnhname   = 'DRIFT'
          tzsource%clongname  = 'ion drift motion'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'CORAY'
          tzsource%clongname  = 'cosmic ray source'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DEPS'
          tzsource%clongname  = 'deposition on snow'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DEPG'
          tzsource%clongname  = 'deposition on graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'REVA'
          tzsource%clongname  = 'rain evaporation'
          tzsource%lavailable = lwarm_ice
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DEPI'
          tzsource%clongname  = 'condensation/deposition on ice'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'NEUT'
          tzsource%clongname  = 'neutralization'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )


        case ( 2 ) SV_ELEC
          ! volumetric charge of cloud droplets
          tzsource%cmnhname   = 'HON'
          tzsource%clongname  = 'homogeneous nucleation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'AUTO'
          tzsource%clongname  = 'autoconversion into rain'
          tzsource%lavailable = lwarm_ice
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'ACCR'
          tzsource%clongname  = 'accretion of cloud droplets'
          tzsource%lavailable = lwarm_ice
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'RIM'
          tzsource%clongname  = 'riming of cloud water'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETG'
          tzsource%clongname  = 'wet growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DRYG'
          tzsource%clongname  = 'dry growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'INCG'
          tzsource%clongname  = 'inductive charge transfer between cloud droplets and graupel'
          tzsource%lavailable = linductive
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETH'
          tzsource%clongname  = 'wet growth of hail'
          tzsource%lavailable = hcloud == 'ICE4'
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'IMLT'
          tzsource%clongname  = 'melting of ice'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'BERFI'
          tzsource%clongname  = 'Bergeron-Findeisen'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SEDI'
          tzsource%clongname  = 'sedimentation'
          tzsource%lavailable = lsedic_ice
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DEPI'
          tzsource%clongname  = 'condensation/deposition on ice'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'NEUT'
          tzsource%clongname  = 'neutralization'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )


        case ( 3 ) SV_ELEC
          ! volumetric charge of rain drops
          tzsource%cmnhname   = 'SFR'
          tzsource%clongname  = 'spontaneous freezing'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'AUTO'
          tzsource%clongname  = 'autoconversion into rain'
          tzsource%lavailable = lwarm_ice
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'ACCR'
          tzsource%clongname  = 'accretion of cloud droplets'
          tzsource%lavailable = lwarm_ice
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'REVA'
          tzsource%clongname  = 'rain evaporation'
          tzsource%lavailable = lwarm_ice
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'ACC'
          tzsource%clongname  = 'accretion of rain  on aggregates'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'CFRZ'
          tzsource%clongname  = 'conversion freezing of rain'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETG'
          tzsource%clongname  = 'wet growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DRYG'
          tzsource%clongname  = 'dry growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'GMLT'
          tzsource%clongname  = 'graupel melting'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETH'
          tzsource%clongname  = 'wet growth of hail'
          tzsource%lavailable = hcloud == 'ICE4'
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'HMLT'
          tzsource%clongname  = 'melting of hail'
          tzsource%lavailable = hcloud == 'ICE4'
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SEDI'
          tzsource%clongname  = 'sedimentation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'NEUT'
          tzsource%clongname  = 'neutralization'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

        case ( 4 ) SV_ELEC
          ! volumetric charge of ice crystals
          tzsource%cmnhname   = 'HON'
          tzsource%clongname  = 'homogeneous nucleation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'AGGS'
          tzsource%clongname  = 'aggregation of snow'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'AUTS'
          tzsource%clongname  = 'autoconversion of ice'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'CFRZ'
          tzsource%clongname  = 'conversion freezing of rain'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETG'
          tzsource%clongname  = 'wet growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DRYG'
          tzsource%clongname  = 'dry growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETH'
          tzsource%clongname  = 'wet growth of hail'
          tzsource%lavailable = hcloud == 'ICE4'
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'IMLT'
          tzsource%clongname  = 'melting of ice'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'BERFI'
          tzsource%clongname  = 'Bergeron-Findeisen'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'NIIS'
          tzsource%clongname  = 'non-inductive charge separation due to ice-snow collisions'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SEDI'
          tzsource%clongname  = 'sedimentation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DEPI'
          tzsource%clongname  = 'condensation/deposition on ice'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'NEUT'
          tzsource%clongname  = 'neutralization'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )


        case ( 5 ) SV_ELEC
          ! volumetric charge of snow
          tzsource%cmnhname   = 'DEPS'
          tzsource%clongname  = 'deposition on snow'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'AGGS'
          tzsource%clongname  = 'aggregation of snow'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'AUTS'
          tzsource%clongname  = 'autoconversion of ice'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'RIM'
          tzsource%clongname  = 'riming of cloud water'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'ACC'
          tzsource%clongname  = 'accretion of rain on snow'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'CMEL'
          tzsource%clongname  = 'conversion melting'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETG'
          tzsource%clongname  = 'wet growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DRYG'
          tzsource%clongname  = 'dry growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'NIIS'
          tzsource%clongname  = 'non-inductive charge separation due to ice-snow collisions'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETH'
          tzsource%clongname  = 'wet growth of hail'
          tzsource%lavailable = hcloud == 'ICE4'
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SEDI'
          tzsource%clongname  = 'sedimentation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'NEUT'
          tzsource%clongname  = 'neutralization'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )


        case ( 6 ) SV_ELEC
          ! volumetric charge of graupel
          tzsource%cmnhname   = 'SFR'
          tzsource%clongname  = 'spontaneous freezing'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DEPG'
          tzsource%clongname  = 'deposition on graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'RIM'
          tzsource%clongname  = 'riming of cloud water'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'ACC'
          tzsource%clongname  = 'accretion of rain  on graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'CMEL'
          tzsource%clongname  = 'conversion melting'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'CFRZ'
          tzsource%clongname  = 'conversion freezing of rain'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETG'
          tzsource%clongname  = 'wet growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'DRYG'
          tzsource%clongname  = 'dry growth of graupel'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'INCG'
          tzsource%clongname  = 'inductive charge transfer between cloud droplets and graupel'
          tzsource%lavailable = linductive
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'GMLT'
          tzsource%clongname  = 'graupel melting'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'WETH'
          tzsource%clongname  = 'wet growth of hail'
          tzsource%lavailable = hcloud == 'ICE4'
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'SEDI'
          tzsource%clongname  = 'sedimentation'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )

          tzsource%cmnhname   = 'NEUT'
          tzsource%clongname  = 'neutralization'
          tzsource%lavailable = .true.
          call Budget_source_add( tbudgets(ibudget), tzsource )


        case ( 7: ) SV_ELEC
          if ( ( hcloud == 'ICE4' .and. ( jsv - nsv_elecbeg + 1 ) == 7 ) ) then
            ! volumetric charge of hail
            tzsource%cmnhname   = 'WETG'
            tzsource%clongname  = 'wet growth of graupel'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'WETH'
            tzsource%clongname  = 'wet growth of hail'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'HMLT'
            tzsource%clongname  = 'melting of hail'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'SEDI'
            tzsource%clongname  = 'sedimentation'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'NEUT'
            tzsource%clongname  = 'neutralization'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

          else if (      ( hcloud == 'ICE3' .and. ( jsv - nsv_elecbeg + 1 ) == 7 ) &
                    .or. ( hcloud == 'ICE4' .and. ( jsv - nsv_elecbeg + 1 ) == 8 ) ) then
            ! Negative ions (NSV_ELECEND case)
            tzsource%cmnhname   = 'DRIFT'
            tzsource%clongname  = 'ion drift motion'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'CORAY'
            tzsource%clongname  = 'cosmic ray source'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'DEPS'
            tzsource%clongname  = 'deposition on snow'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'DEPG'
            tzsource%clongname  = 'deposition on graupel'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'REVA'
            tzsource%clongname  = 'rain evaporation'
            tzsource%lavailable = lwarm_ice
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'DEPI'
            tzsource%clongname  = 'condensation/deposition on ice'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

            tzsource%cmnhname   = 'NEUT'
            tzsource%clongname  = 'neutralization'
            tzsource%lavailable = .true.
            call Budget_source_add( tbudgets(ibudget), tzsource )

          else
            call Print_msg( NVERB_FATAL, 'BUD', 'Ini_budget', 'unknown electricity budget' )
          end if

      end select SV_ELEC


    else if ( jsv >= nsv_lgbeg .and. jsv <= nsv_lgend ) then SV_VAR
      !Lagrangian variables


    else if ( jsv >= nsv_ppbeg .and. jsv <= nsv_ppend ) then SV_VAR
      !Passive pollutants


#ifdef MNH_FOREFIRE
    else if ( jsv >= nsv_ffbeg .and. jsv <= nsv_ffend ) then SV_VAR
      !Forefire

#endif
    else if ( jsv >= nsv_csbeg .and. jsv <= nsv_csend ) then SV_VAR
      !Conditional sampling


    else if ( jsv >= nsv_chembeg .and. jsv <= nsv_chemend ) then SV_VAR
      !Chemical case
      tzsource%cmnhname   = 'CHEM'
      tzsource%clongname  = 'chemistry activity'
      tzsource%lavailable = .true.
      call Budget_source_add( tbudgets(ibudget), tzsource )

      tzsource%cmnhname   = 'NEGA'
      tzsource%clongname  = 'negativity correction'
      tzsource%lavailable = .true.
      call Budget_source_add( tbudgets(ibudget), tzsource )


    else if ( jsv >= nsv_chicbeg .and. jsv <= nsv_chicend ) then SV_VAR
      !Ice phase chemistry


    else if ( jsv >= nsv_aerbeg .and. jsv <= nsv_aerend ) then SV_VAR
      !Chemical aerosol case
      tzsource%cmnhname   = 'NEGA'
      tzsource%clongname  = 'negativity correction'
      tzsource%lavailable = lorilam
      call Budget_source_add( tbudgets(ibudget), tzsource )

    else if ( jsv >= nsv_aerdepbeg .and. jsv <= nsv_aerdepend ) then SV_VAR
      !Aerosol wet deposition

    else if ( jsv >= nsv_dstbeg .and. jsv <= nsv_dstend ) then SV_VAR
      !Dust

    else if ( jsv >= nsv_dstdepbeg .and. jsv <= nsv_dstdepend ) then SV_VAR
      !Dust wet deposition

    else if ( jsv >= nsv_sltbeg .and. jsv <= nsv_sltend ) then SV_VAR
      !Salt

    else if ( jsv >= nsv_sltdepbeg .and. jsv <= nsv_sltdepend ) then SV_VAR
      !Salt wet deposition

    else if ( jsv >= nsv_snwbeg .and. jsv <= nsv_snwend ) then SV_VAR
      !Snow
      tzsource%cmnhname   = 'SNSUB'
      tzsource%clongname  = 'blowing snow sublimation'
      tzsource%lavailable = lblowsnow .and. lsnowsubl
      call Budget_source_add( tbudgets(ibudget), tzsource )

      tzsource%cmnhname   = 'SNSED'
      tzsource%clongname  = 'blowing snow sedimentation'
      tzsource%lavailable = lblowsnow
      call Budget_source_add( tbudgets(ibudget), tzsource )


    else if ( jsv >= nsv_lnoxbeg .and. jsv <= nsv_lnoxend ) then SV_VAR
      !LiNOX passive tracer

    else SV_VAR
      call Print_msg( NVERB_FATAL, 'BUD', 'Ini_budget', 'unknown scalar variable' )
    end if SV_VAR


    call Sourcelist_sort_compact( tbudgets(ibudget) )

    call Sourcelist_scan( tbudgets(ibudget), cbulist_rsv )
  end if
end do SV_BUDGETS

call Ini_budget_groups( tbudgets, ibudim1, ibudim2, ibudim3 )

if ( tbudgets(NBUDGET_U)  %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_U),   cbulist_ru   )
if ( tbudgets(NBUDGET_V)  %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_V),   cbulist_rv   )
if ( tbudgets(NBUDGET_W)  %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_W),   cbulist_rw   )
if ( tbudgets(NBUDGET_TH) %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_TH),  cbulist_rth  )
if ( tbudgets(NBUDGET_TKE)%lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_TKE), cbulist_rtke )
if ( tbudgets(NBUDGET_RV) %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_RV),  cbulist_rrv  )
if ( tbudgets(NBUDGET_RC) %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_RC),  cbulist_rrc  )
if ( tbudgets(NBUDGET_RR) %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_RR),  cbulist_rrr  )
if ( tbudgets(NBUDGET_RI) %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_RI),  cbulist_rri  )
if ( tbudgets(NBUDGET_RS) %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_RS),  cbulist_rrs  )
if ( tbudgets(NBUDGET_RG) %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_RG),  cbulist_rrg  )
if ( tbudgets(NBUDGET_RH) %lenabled ) call Sourcelist_nml_compact( tbudgets(NBUDGET_RH),  cbulist_rrh  )
if ( lbu_rsv )                        call Sourcelist_sv_nml_compact( cbulist_rsv  )
end subroutine Ini_budget


subroutine Budget_source_add( tpbudget, tpsource, odonotinit, ooverwrite )
  use modd_budget, only: tbudgetdata, tbusourcedata

  type(tbudgetdata),   intent(inout) :: tpbudget
  type(tbusourcedata), intent(in)    :: tpsource ! Metadata basis
  logical, optional,   intent(in)    :: odonotinit
  logical, optional,   intent(in)    :: ooverwrite

  character(len=4) :: ynum
  integer          :: isourcenumber

  call Print_msg( NVERB_DEBUG, 'BUD', 'Budget_source_add', 'called for ' // Trim( tpbudget%cname ) &
                  // ': ' // Trim( tpsource%cmnhname ) )

  isourcenumber = tpbudget%nsources + 1
  if ( isourcenumber > tpbudget%nsourcesmax ) then
    Write( ynum, '( i4 )' ) tpbudget%nsourcesmax
    cmnhmsg(1) = 'Insufficient max number of source terms (' // Trim(ynum) // ') for budget ' // Trim( tpbudget%cname )
    cmnhmsg(2) = 'Please increaze value of parameter NSOURCESMAX'
    call Print_msg( NVERB_FATAL, 'BUD', 'Budget_source_add' )
  else
    tpbudget%nsources = tpbudget%nsources + 1
  end if

  ! Copy metadata from provided tpsource
  ! Modifications to source term metadata done with the other dummy arguments
  tpbudget%tsources(isourcenumber) = tpsource

  if ( present( odonotinit ) ) tpbudget%tsources(isourcenumber)%ldonotinit = odonotinit

  if ( present( ooverwrite ) ) tpbudget%tsources(isourcenumber)%loverwrite = ooverwrite
end subroutine Budget_source_add


subroutine Ini_budget_groups( tpbudgets, kbudim1, kbudim2, kbudim3 )
  use modd_budget,     only: tbudgetdata
  use modd_field,      only: TYPEINT, TYPEREAL
  use modd_parameters, only: NMNHNAMELGTMAX, NSTDNAMELGTMAX, NLONGNAMELGTMAX, NUNITLGTMAX, NCOMMENTLGTMAX

  use mode_tools,  only: Quicksort

  type(tbudgetdata), dimension(:), intent(inout) :: tpbudgets
  integer,                         intent(in)    :: kbudim1
  integer,                         intent(in)    :: kbudim2
  integer,                         intent(in)    :: kbudim3

  character(len=NMNHNAMELGTMAX)      :: ymnhname
  character(len=NSTDNAMELGTMAX)      :: ystdname
  character(len=NLONGNAMELGTMAX)     :: ylongname
  character(len=NUNITLGTMAX)         :: yunits
  character(len=NCOMMENTLGTMAX)      :: ycomment
  integer                            :: ji, jj, jk
  integer                            :: isources ! Number of source terms in a budget
  integer                            :: inbgroups ! Number of budget groups
  integer                            :: ival
  integer                            :: icount
  integer                            :: ivalmax, ivalmin
  integer                            :: igrid
  integer                            :: itype
  integer                            :: idims
  integer, dimension(:), allocatable :: igroups ! Temporary array to store sorted group numbers
  integer, dimension(:), allocatable :: ipos    ! Temporary array to store initial position of group numbers
  real                               :: zval
  real                               :: zvalmax, zvalmin

  call Print_msg( NVERB_DEBUG, 'BUD', 'Ini_budget_groups', 'called' )

  BUDGETS: do ji = 1, size( tpbudgets )
    ENABLED: if ( tpbudgets(ji)%lenabled ) then
      isources = size( tpbudgets(ji)%tsources )
      do jj = 1, isources
        ! Check if ngroup is an allowed value
        if ( tpbudgets(ji)%tsources(jj)%ngroup < 0 ) then
          call Print_msg( NVERB_ERROR, 'BUD', 'Ini_budget', 'negative group value is not allowed' )
          tpbudgets(ji)%tsources(jj)%ngroup = 0
        end if

        if ( tpbudgets(ji)%tsources(jj)%ngroup > 0 ) tpbudgets(ji)%tsources(jj)%lenabled = .true.
      end do

      !Count the number of groups of source terms
      !ngroup=1 is for individual entries, >1 values are groups
      allocate( igroups(isources ) )
      allocate( ipos   (isources ) )
      igroups(:) = tpbudgets(ji)%tsources(:)%ngroup
      ipos(:) = [ ( jj, jj = 1, isources ) ]

      !Sort the group list number
      call Quicksort( igroups, 1, isources, ipos )

      !Count the number of different groups
      !and renumber the entries (from 1 to inbgroups)
      inbgroups = 0
      ival = igroups(1)
      if ( igroups(1) /= 0 ) then
        inbgroups = 1
        igroups(1) = inbgroups
      end if
      do jj = 2, isources
        if ( igroups(jj) == 1 ) then
          inbgroups = inbgroups + 1
          igroups(jj) = inbgroups
        else if ( igroups(jj) > 0 ) then
          if ( igroups(jj) /= ival ) then
            ival = igroups(jj)
            inbgroups = inbgroups + 1
          end if
          igroups(jj) = inbgroups
        end if
      end do

      !Write the igroups values to the budget structure
      do jj = 1, isources
        tpbudgets(ji)%tsources(ipos(jj))%ngroup = igroups(jj)
      end do

      !Allocate the group structure + populate it
      tpbudgets(ji)%ngroups = inbgroups
      allocate( tpbudgets(ji)%tgroups(inbgroups) )

      do jj = 1, inbgroups
        !Search the list of sources for each group
        !not the most efficient algorithm but do the job
        icount = 0
        do jk = 1, isources
          if ( tpbudgets(ji)%tsources(jk)%ngroup == jj ) then
            icount = icount + 1
            ipos(icount) = jk !ipos is reused as a temporary work array
          end if
        end do
        tpbudgets(ji)%tgroups(jj)%nsources = icount

        allocate( tpbudgets(ji)%tgroups(jj)%nsourcelist(icount) )
        tpbudgets(ji)%tgroups(jj)%nsourcelist(:) = ipos(1 : icount)

        ! Set the name of the field
        ymnhname = tpbudgets(ji)%tsources(ipos(1))%cmnhname
        do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
          ymnhname = trim( ymnhname ) // '_' // trim( tpbudgets(ji)%tsources(ipos(jk))%cmnhname )
        end do
        tpbudgets(ji)%tgroups(jj)%cmnhname = ymnhname

        ! Set the standard name (CF convention)
        if ( tpbudgets(ji)%tgroups(jj)%nsources == 1 ) then
          ystdname = tpbudgets(ji)%tsources(ipos(1))%cstdname
        else
          ! The CF standard name is probably wrong if combining several source terms => set to ''
          ystdname = ''
        end if
        tpbudgets(ji)%tgroups(jj)%cstdname = ystdname

        ! Set the long name (CF convention)
        ylongname = tpbudgets(ji)%tsources(ipos(1))%clongname
        do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
          ylongname = trim( ylongname ) // ' + ' // tpbudgets(ji)%tsources(ipos(jk))%clongname
        end do
        tpbudgets(ji)%tgroups(jj)%clongname = ylongname

        ! Set the units
        yunits = tpbudgets(ji)%tsources(ipos(1))%cunits
        do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
          if ( trim( yunits ) /= trim( tpbudgets(ji)%tsources(ipos(jk))%cunits ) ) then
            call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget',                               &
                            'incompatible units for the different source terms of the group ' &
                            //trim( tpbudgets(ji)%tgroups(jj)%cmnhname ) )
            yunits = 'unknown'
          end if
        end do
        tpbudgets(ji)%tgroups(jj)%cunits = yunits

        ! Set the comment
        ! It is composed of the source comment followed by the clongnames of the different sources
        ycomment = trim( tpbudgets(ji)%tsources(ipos(1))%ccomment ) // ': '// trim( tpbudgets(ji)%tsources(ipos(1))%clongname )
        do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
          ycomment = trim( ycomment ) // ', ' // trim( tpbudgets(ji)%tsources(ipos(jk))%clongname )
        end do
        ycomment = trim( ycomment ) // ' source term'
        if ( tpbudgets(ji)%tgroups(jj)%nsources > 1 ) ycomment = trim( ycomment ) // 's'
        tpbudgets(ji)%tgroups(jj)%ccomment = ycomment

        ! Set the Arakawa grid
        igrid = tpbudgets(ji)%tsources(ipos(1))%ngrid
        do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
          if ( igrid /= tpbudgets(ji)%tsources(ipos(jk))%ngrid ) then
            call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget',                                             &
                            'different Arakawa grid positions for the different source terms of the group ' &
                            //trim( tpbudgets(ji)%tgroups(jj)%cmnhname ) )
          end if
        end do
        tpbudgets(ji)%tgroups(jj)%ngrid = igrid

        ! Set the data type
        itype = tpbudgets(ji)%tsources(ipos(1))%ntype
        do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
          if ( itype /= tpbudgets(ji)%tsources(ipos(jk))%ntype ) then
            call Print_msg( NVERB_FATAL, 'BUD', 'Ini_budget',                                             &
                            'incompatible data types for the different source terms of the group ' &
                            //trim( tpbudgets(ji)%tgroups(jj)%cmnhname ) )
          end if
        end do
        tpbudgets(ji)%tgroups(jj)%ntype = itype

        ! Set the number of dimensions
        idims = tpbudgets(ji)%tsources(ipos(1))%ndims
        do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
          if ( idims /= tpbudgets(ji)%tsources(ipos(jk))%ndims ) then
            call Print_msg( NVERB_FATAL, 'BUD', 'Ini_budget',                                             &
                            'incompatible number of dimensions for the different source terms of the group ' &
                            //trim( tpbudgets(ji)%tgroups(jj)%cmnhname ) )
          end if
        end do
        tpbudgets(ji)%tgroups(jj)%ndims = idims

        ! Set the fill values
        if ( tpbudgets(ji)%tgroups(jj)%ntype == TYPEINT ) then
          ival = tpbudgets(ji)%tsources(ipos(1))%nfillvalue
          do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
            if ( ival /= tpbudgets(ji)%tsources(ipos(jk))%nfillvalue ) then
              call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget',                                             &
                              'different (integer) fill values for the different source terms of the group ' &
                              //trim( tpbudgets(ji)%tgroups(jj)%cmnhname ) )
            end if
          end do
          tpbudgets(ji)%tgroups(jj)%nfillvalue = ival
        end if

        if ( tpbudgets(ji)%tgroups(jj)%ntype == TYPEREAL ) then
          zval = tpbudgets(ji)%tsources(ipos(1))%xfillvalue
          do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
            if ( zval /= tpbudgets(ji)%tsources(ipos(jk))%xfillvalue ) then
              call Print_msg( NVERB_WARNING, 'BUD', 'Ini_budget',                                             &
                              'different (real) fill values for the different source terms of the group ' &
                              //trim( tpbudgets(ji)%tgroups(jj)%cmnhname ) )
            end if
          end do
          tpbudgets(ji)%tgroups(jj)%xfillvalue = zval
        end if

        ! Set the valid min/max values
        ! Take the min or max of all the sources
        ! Maybe, it would be better to take the sum? (if same sign, if not already the maximum allowed value for this type)
        if ( tpbudgets(ji)%tgroups(jj)%ntype == TYPEINT ) then
          ivalmin = tpbudgets(ji)%tsources(ipos(1))%nvalidmin
          ivalmax = tpbudgets(ji)%tsources(ipos(1))%nvalidmax
          do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
            ivalmin = min( ivalmin, tpbudgets(ji)%tsources(ipos(jk))%nvalidmin )
            ivalmax = max( ivalmax, tpbudgets(ji)%tsources(ipos(jk))%nvalidmax )
          end do
          tpbudgets(ji)%tgroups(jj)%nvalidmin = ivalmin
          tpbudgets(ji)%tgroups(jj)%nvalidmax = ivalmax
        end if

        if ( tpbudgets(ji)%tgroups(jj)%ntype == TYPEREAL ) then
          zvalmin = tpbudgets(ji)%tsources(ipos(1))%xvalidmin
          zvalmax = tpbudgets(ji)%tsources(ipos(1))%xvalidmax
          do jk = 2, tpbudgets(ji)%tgroups(jj)%nsources
            zvalmin = min( zvalmin, tpbudgets(ji)%tsources(ipos(jk))%xvalidmin )
            zvalmax = max( zvalmax, tpbudgets(ji)%tsources(ipos(jk))%xvalidmax )
          end do
          tpbudgets(ji)%tgroups(jj)%xvalidmin = zvalmin
          tpbudgets(ji)%tgroups(jj)%xvalidmax = zvalmax
        end if

        allocate( tpbudgets(ji)%tgroups(jj)%xdata(kbudim1, kbudim2, kbudim3 ) )
        tpbudgets(ji)%tgroups(jj)%xdata(:, :, :) = 0.
      end do

      deallocate( igroups )
      deallocate( ipos )

      !Check that a group does not contain more than 1 source term with ldonotinit=.true.
      do jj = 1, inbgroups
        if ( tpbudgets(ji)%tgroups(jj)%nsources > 1 ) then
          do jk = 1, tpbudgets(ji)%tgroups(jj)%nsources
            if ( tpbudgets(ji)%tsources(tpbudgets(ji)%tgroups(jj)%nsourcelist(jk) )%ldonotinit ) &
              call Print_msg( NVERB_FATAL, 'BUD', 'Ini_budget', &
                              'a group with more than 1 source term may not contain sources with ldonotinit=true' )
            if ( tpbudgets(ji)%tsources(tpbudgets(ji)%tgroups(jj)%nsourcelist(jk) )%loverwrite ) &
              call Print_msg( NVERB_FATAL, 'BUD', 'Ini_budget', &
                              'a group with more than 1 source term may not contain sources with loverwrite=true' )
          end do
        end if
      end do

    end if ENABLED
  end do BUDGETS

end subroutine Ini_budget_groups


subroutine Sourcelist_sort_compact( tpbudget )
  !Sort the list of sources to put the non-available source terms at the end of the list
  !and compact the list
  use modd_budget, only: tbudgetdata, tbusourcedata

  type(tbudgetdata), intent(inout) :: tpbudget

  integer                                        :: ji
  integer                                        :: isrc_avail, isrc_notavail
  type(tbusourcedata), dimension(:), allocatable :: tzsources_avail
  type(tbusourcedata), dimension(:), allocatable :: tzsources_notavail

  isrc_avail    = 0
  isrc_notavail = 0

  Allocate( tzsources_avail   (tpbudget%nsources) )
  Allocate( tzsources_notavail(tpbudget%nsources) )

  !Separate source terms available or not during the execution
  !(based on the criteria provided to Budget_source_add and stored in lavailable field)
  do ji = 1, tpbudget%nsources
    if ( tpbudget%tsources(ji)%lavailable ) then
      isrc_avail = isrc_avail + 1
      tzsources_avail(isrc_avail) = tpbudget%tsources(ji)
    else
      isrc_notavail = isrc_notavail + 1
      tzsources_notavail(isrc_notavail) = tpbudget%tsources(ji)
    end if
  end do

  !Reallocate/compact the source list
  if ( Allocated( tpbudget%tsources ) ) Deallocate( tpbudget%tsources )
  Allocate( tpbudget%tsources( tpbudget%nsources ) )

  tpbudget%nsourcesmax = tpbudget%nsources
  !Limit the number of sources to the available list
  tpbudget%nsources    = isrc_avail

  !Fill the source list beginning with the available sources and finishing with the non-available ones
  do ji = 1, isrc_avail
    tpbudget%tsources(ji) = tzsources_avail(ji)
  end do

  do ji = 1, isrc_notavail
    tpbudget%tsources(isrc_avail + ji) = tzsources_notavail(ji)
  end do

end subroutine Sourcelist_sort_compact


subroutine Sourcelist_scan( tpbudget, hbulist )
  use modd_budget, only: tbudgetdata

  type(tbudgetdata),              intent(inout) :: tpbudget
  character(len=*), dimension(:), intent(in)    :: hbulist

  character(len=:), allocatable :: yline
  character(len=:), allocatable :: ysrc
  character(len=:), dimension(:), allocatable :: ymsg
  integer                       :: idx
  integer                       :: igroup
  integer                       :: igroup_idx
  integer                       :: ipos
  integer                       :: istart
  integer                       :: ji

  istart = 1

  ! Case 'LIST_AVAIL': list all the available source terms
  if ( Size( hbulist ) > 0 ) then
    if ( Trim( hbulist(1) ) == 'LIST_AVAIL' ) then
      Allocate( character(len=65) :: ymsg(tpbudget%nsources + 1) )
      ymsg(1) = '---------------------------------------------------------------------'
      ymsg(2) = 'Available source terms for budget ' // Trim( tpbudget%cname )
      Write( ymsg(3), '( A32, " ", A32 )' ) 'Name', 'Long name'
      idx = 3
      do ji = 1, tpbudget%nsources
        if ( All( tpbudget%tsources(ji)%cmnhname /= [ 'INIF' , 'ENDF', 'AVEF' ] ) ) then
          idx = idx + 1
          Write( ymsg(idx), '( A32, " ", A32 )' ) tpbudget%tsources(ji)%cmnhname, tpbudget%tsources(ji)%clongname
        end if
      end do
      ymsg(tpbudget%nsources + 1 ) = '---------------------------------------------------------------------'
      call Print_msg_multi( NVERB_WARNING, 'BUD', 'Sourcelist_scan', ymsg )
      !To not read the 1st line again
      istart = 2
    end if
  end if

  ! Case 'LIST_ALL': list all the source terms
  if ( Size( hbulist ) > 0 ) then
    if ( Trim( hbulist(1) ) == 'LIST_ALL' ) then
      Allocate( character(len=65) :: ymsg(tpbudget%nsourcesmax + 1) )
      ymsg(1) = '---------------------------------------------------------------------'
      ymsg(2) = 'Source terms for budget ' // Trim( tpbudget%cname )
      Write( ymsg(3), '( A32, " ", A32 )' ) 'Name', 'Long name'
      idx = 3
      do ji = 1, tpbudget%nsourcesmax
        if ( All( tpbudget%tsources(ji)%cmnhname /= [ 'INIF' , 'ENDF', 'AVEF' ] ) ) then
          idx = idx + 1
          Write( ymsg(idx), '( A32, " ", A32 )' ) tpbudget%tsources(ji)%cmnhname, tpbudget%tsources(ji)%clongname
        end if
      end do
      ymsg(tpbudget%nsourcesmax + 1 ) = '---------------------------------------------------------------------'
      call Print_msg_multi( NVERB_WARNING, 'BUD', 'Sourcelist_scan', ymsg )
      !To not read the 1st line again
      istart = 2
    end if
  end if

  ! Case 'ALL': enable all available source terms
  if ( Size( hbulist ) > 0 ) then
    if ( Trim( hbulist(1) ) == 'ALL' ) then
      do ji = 1, tpbudget%nsources
        tpbudget%tsources(ji)%ngroup = 1
      end do
      return
    end if
  end if

  !Always enable INIF, ENDF and AVEF terms
  ipos = Source_find( tpbudget, 'INIF' )
  if ( ipos < 1 ) call Print_msg( NVERB_FATAL, 'BUD', 'Sourcelist_scan', 'source term ' // Trim( tpbudget%cname ) &
                            // ': INIF not found' )
  tpbudget%tsources(ipos)%ngroup = 1

  ipos = Source_find( tpbudget, 'ENDF' )
  if ( ipos < 1 ) call Print_msg( NVERB_FATAL, 'BUD', 'Sourcelist_scan', 'source term ' // Trim( tpbudget%cname ) &
                            // ': ENDF not found' )
  tpbudget%tsources(ipos)%ngroup = 1

  ipos = Source_find( tpbudget, 'AVEF' )
  if ( ipos < 1 ) call Print_msg( NVERB_FATAL, 'BUD', 'Sourcelist_scan', 'source term ' // Trim( tpbudget%cname ) &
                            // ': AVEF not found' )
  tpbudget%tsources(ipos)%ngroup = 1

  !igroup_idx start at 2 because 1 is reserved for individually stored source terms
  igroup_idx = 2

  do ji = istart, Size( hbulist )
    if ( Len_trim( hbulist(ji) ) > 0 ) then
      ! Scan the line and separate the different sources (separated by + signs)
      yline = Trim(hbulist(ji))

      idx = Index( yline, '+' )
      if ( idx < 1 ) then
        igroup = 1
      else
        igroup = igroup_idx
        igroup_idx = igroup_idx + 1
      end if

      do
        idx = Index( yline, '+' )
        if ( idx < 1 ) then
          ysrc = yline
        else
          ysrc = yline(1 : idx - 1)
          yline = yline(idx + 1 :)
        end if

        !Check if the source is known
        if ( Len_trim( ysrc ) > 0 ) then
          ipos = Source_find( tpbudget, ysrc )

          if ( ipos > 0 ) then
            call Print_msg( NVERB_DEBUG, 'BUD', 'Sourcelist_scan', 'source term ' // Trim( tpbudget%cname ) &
                            // ': ' // ysrc // ' found' )

            if ( .not.  tpbudget%tsources(ipos)%lavailable ) then
              call Print_msg( NVERB_WARNING, 'BUD', 'Sourcelist_scan', 'source term ' // Trim( tpbudget%cname ) &
                              // ': ' // ysrc // ' not available' )
              tpbudget%tsources(ipos)%ngroup = 0
            else
              tpbudget%tsources(ipos)%ngroup = igroup
            end if
          else
            call Print_msg( NVERB_ERROR, 'BUD', 'Sourcelist_scan', 'source term ' // Trim( tpbudget%cname ) &
                            // ': ' // ysrc // ' not found' )
          end if
        end if

        if ( idx < 1 ) exit
      end do
    end if
  end do
end subroutine Sourcelist_scan


subroutine Sourcelist_nml_compact( tpbudget, hbulist )
  !This subroutine reduce the size of the hbulist to the minimum
  !The list is generated from the group list
  use modd_budget, only: NBULISTMAXLEN, tbudgetdata

  type(tbudgetdata),                                       intent(in)    :: tpbudget
  character(len=NBULISTMAXLEN), dimension(:), allocatable, intent(inout) :: hbulist

  integer :: idx
  integer :: isource
  integer :: jg
  integer :: js

  if ( Allocated( hbulist ) ) Deallocate( hbulist )

  if ( tpbudget%ngroups < 3 ) then
    call Print_msg( NVERB_ERROR, 'BUD', 'Sourcelist_nml_compact', 'ngroups is too small' )
    return
  end if

  Allocate( character(len=NBULISTMAXLEN) :: hbulist(tpbudget%ngroups - 3) )
  hbulist(:) = ''

  idx = 0
  do jg = 1, tpbudget%ngroups
    if ( tpbudget%tgroups(jg)%nsources < 1 ) then
      call Print_msg( NVERB_ERROR, 'BUD', 'Sourcelist_nml_compact', 'no source for group' )
      cycle
    end if

    !Do not put 'INIF', 'ENDF', 'AVEF' in hbulist because their presence is automatic if the corresponding budget is enabled
    isource = tpbudget%tgroups(jg)%nsourcelist(1)
    if ( Any( tpbudget%tsources(isource)%cmnhname ==  [ 'INIF', 'ENDF', 'AVEF' ] ) ) cycle

    idx = idx + 1
#if 0
    !Do not do this way because the group cmnhname may be truncated (NMNHNAMELGTMAX is smaller than NBULISTMAXLEN)
    !and the name separator is different ('_')
    hbulist(idx) = Trim( tpbudget%tgroups(jg)%cmnhname )
#else
    do js = 1, tpbudget%tgroups(jg)%nsources
      isource = tpbudget%tgroups(jg)%nsourcelist(js)
      hbulist(idx) = Trim( hbulist(idx) ) // Trim( tpbudget%tsources(isource)%cmnhname )
      if ( js < tpbudget%tgroups(jg)%nsources ) hbulist(idx) = Trim( hbulist(idx) ) // '+'
    end do
#endif
  end do
end subroutine Sourcelist_nml_compact


subroutine Sourcelist_sv_nml_compact( hbulist )
  !This subroutine reduce the size of the hbulist
  !For SV variables the reduction is simpler than for other variables
  !because it is too complex to do this cleanly (the enabled source terms are different for each scalar variable)
  use modd_budget, only: NBULISTMAXLEN, tbudgetdata

  character(len=*), dimension(:), allocatable, intent(inout) :: hbulist

  character(len=NBULISTMAXLEN), dimension(:), allocatable :: ybulist_new
  integer :: ilines
  integer :: ji

  ilines = 0
  do ji = 1, Size( hbulist )
    if ( Len_trim(hbulist(ji)) > 0 ) ilines = ilines + 1
  end do

  Allocate( ybulist_new(ilines) )

  ilines = 0
  do ji = 1, Size( hbulist )
    if ( Len_trim(hbulist(ji)) > 0 ) then
      ilines = ilines + 1
      ybulist_new(ilines) = Trim( hbulist(ji) )
    end if
  end do

  call Move_alloc( from = ybulist_new, to = hbulist )
end subroutine Sourcelist_sv_nml_compact


pure function Source_find( tpbudget, hsource ) result( ipos )
  use modd_budget,     only: tbudgetdata

  type(tbudgetdata), intent(in) :: tpbudget
  character(len=*),  intent(in) :: hsource
  integer :: ipos

  integer :: ji
  logical :: gfound

  ipos = -1
  gfound = .false.
  do ji = 1, tpbudget%nsourcesmax
    if ( Trim( hsource ) == Trim ( tpbudget%tsources(ji)%cmnhname ) ) then
      gfound = .true.
      ipos = ji
      exit
    end if
  end do

end function Source_find

end module mode_ini_budget
