!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ##################
      MODULE MODD_BUDGET
!     ##################
!
!!****  *MODD_BUDGET* - declaration of budget variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the budget
!     variables
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_BUDGET)
!!          
!!    AUTHOR
!!    ------
!!	P. Hereil   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        23/02/95 
!!      J.-P. Lafore    10/02/98    adding of rhodj declaration for budget  
!!      V. Ducrocq      4/06/99     //
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 19/07/2019: parameters to identify budget number
!  P. Wautelet 15/11/2019: remove unused CBURECORD variable
!  P. Wautelet 17/01/2020: add new budget data types
!  P. Wautelet 27/01/2020: use the tfieldmetadata_base abstract datatype
!  P. Wautelet 28/01/2020: add trhodj in tbudgetdata datatype
!  P. Wautelet 09/03/2020: add tburhodj variable
!  P. Wautelet 17/04/2020: set default values for budgets switch values
!  P. Wautelet 23/04/2020: add nid in tbudgetdata datatype
!  P. Wautelet 17/08/2020: add xtmplesstore in tbudgetdata datatype
!  P. Wautelet 08/10/2020: add clessource in tbudgetdata datatype
!  P. Wautelet 08/12/2020: add nbusubwrite and nbutotwrite
!  P. Wautelet 11/01/2021: remove nbuwrnb (replaced by nbusubwrite)
!  P. Wautelet 14/01/2021: change xbusurf type to integer (+ rename it to nbusurf)
!  P. Wautelet 03/03/2021: add tbudiachrometadata type (useful to pass more information to Write_diachro)
!  P. Wautelet 17/03/2021: choose source terms for budgets with character strings instead of multiple integer variables
!  P. Wautelet 30/03/2021: budgets: cartesian subdomain limits are defined in the physical domain
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------

use modd_field,      only: tfieldmetadata_base
use modd_parameters, only: NBUNAMELGTMAX, NCOMMENTLGTMAX

implicit none

public

integer, parameter :: NBULISTMAXLEN   = 128
integer, parameter :: NBULISTMAXLINES = 50

integer, parameter :: NBUDGET_RHO = 0  ! Reference number for budget of RhoJ
integer, parameter :: NBUDGET_U   = 1  ! Reference number for budget of RhoJu  and/or LES budgets with u
integer, parameter :: NBUDGET_V   = 2  ! Reference number for budget of RhoJv  and/or LES budgets with u
integer, parameter :: NBUDGET_W   = 3  ! Reference number for budget of RhoJw  and/or LES budgets with u
integer, parameter :: NBUDGET_TH  = 4  ! Reference number for budget of RhoJTh and/or LES budgets with th
integer, parameter :: NBUDGET_TKE = 5  ! Reference number for budget of RhoJTke and/or LES budgets with Tke
integer, parameter :: NBUDGET_RV  = 6  ! Reference number for budget of RhoJrv and/or LES budgets with rv
integer, parameter :: NBUDGET_RC  = 7  ! Reference number for budget of RhoJrc and/or LES budgets with rc
integer, parameter :: NBUDGET_RR  = 8  ! Reference number for budget of RhoJrr and/or LES budgets with rr
integer, parameter :: NBUDGET_RI  = 9  ! Reference number for budget of RhoJri and/or LES budgets with ri
integer, parameter :: NBUDGET_RS  = 10 ! Reference number for budget of RhoJrs and/or LES budgets with rs
integer, parameter :: NBUDGET_RG  = 11 ! Reference number for budget of RhoJrg and/or LES budgets with rg
integer, parameter :: NBUDGET_RH  = 12 ! Reference number for budget of RhoJrh and/or LES budgets with rh
integer, parameter :: NBUDGET_SV1 = 13 ! Reference number for 1st budget of RhoJsv and/or LES budgets with sv

integer, parameter :: NMAXLEVELS       = 7
integer, parameter :: NLVL_ROOT        = 0
integer, parameter :: NLVL_CATEGORY    = 1
integer, parameter :: NLVL_SUBCATEGORY = 2
integer, parameter :: NLVL_GROUP       = 3
integer, parameter :: NLVL_SHAPE       = 4
integer, parameter :: NLVL_TIMEAVG     = 5
integer, parameter :: NLVL_NORM        = 6
integer, parameter :: NLVL_MASK        = 7

#ifdef MNH_IOCDF4
character(len=*), dimension(NMAXLEVELS), parameter :: CNCGROUPNAMES = [ &
                                         'category   ', &  !Name of the different type of groups/levels in the netCDF file
                                         'subcategory', &
                                         'group      ', &
                                         'shape      ', &
                                         'timeavg    ', &
                                         'norm       ', &
                                         'mask       '  ]
#endif

integer :: nbudgets ! Number of budget categories


type, extends( tfieldmetadata_base ) :: tbusourcedata
  integer :: ngroup = 0 ! Number of the source term group in which storing the source term
                        !  (0: no store, 1: individual store, >1: number of the group)
  logical :: lavailable = .false. ! If true, the source is available in the run (conditions to access it are met),
                                  ! but it doesn't mean it is used (see lenabled field)
  logical :: lenabled   = .false.
  logical :: ldonotinit = .false. ! if true, does not need a call to Budget_store_init
                                  ! It may be true only if the source term is in a group not containing other sources
  logical :: loverwrite = .false. ! if true, source term values will overwrite the previous ones
                                  ! It may be true only if the source term is in a group not containing other sources
end type tbusourcedata

type, extends( tfieldmetadata_base ) :: tbugroupdata
  integer :: nsources = 0 ! Number of source terms composing this group
  integer, dimension(:),     allocatable :: nsourcelist ! List of the source terms composing this group
  real,    dimension(:,:,:), allocatable :: xdata ! Array to store the budget data
end type tbugroupdata

type, extends( tfieldmetadata_base ) :: tburhodata
  real, dimension(:,:,:), allocatable :: xdata ! Array to store the budget data
end type tburhodata

type :: tbudiachrometadata
  character(len=NBUNAMELGTMAX),  dimension(NMAXLEVELS) :: clevels  = '' !Name of the different groups/levels in the netCDF file
  character(len=NCOMMENTLGTMAX), dimension(NMAXLEVELS) :: ccomments ='' !Comments for the different groups/levels in the netCDF file
  character(len=1)              :: cdirection   = ''                    !Used for 2pt correlation and spectrum
  logical :: lmobile    = .false.                                       !Is the domain moving? (ie for aircrafts and balloons)
  logical :: licompress = .false.
  logical :: ljcompress = .false.
  logical :: lkcompress = .false.
  logical :: ltcompress = .false. ! true if values are time averaged (can be on multiple time periods)
  logical :: lnorm      = .false. ! true if values are normalized
  logical, dimension(NMAXLEVELS) :: lleveluse = .false.
  integer :: nil = -1 !Cartesian box boundaries in physical domain coordinates
  integer :: nih = -1
  integer :: njl = -1
  integer :: njh = -1
  integer :: nkl = -1
  integer :: nkh = -1
  integer :: nsv = -1 !Reference number of the corresponding scalar variable
end type tbudiachrometadata

type tbudgetdata
  character(len=NBUNAMELGTMAX)  :: cname    = ''
  character(len=NCOMMENTLGTMAX) :: ccomment = ''
  character(len=100)            :: clessource = '' ! Last source stored
  integer :: nid         = -1 !Identifier number (based on parameters NBUDGET_*)
  integer :: ngroups     = 0 !Number of groups of source terms to store
  integer :: nsources    = 0 !Number of available source terms
  integer :: nsourcesmax = 0 !Maximum number of source terms
  integer :: ntmpstoresource = 0 !Reference of the source term using the xtmpstore array
  logical :: lenabled = .false. ! True if corresponding budget flag is set to true
  real, dimension(:,:,:), allocatable :: xtmpstore ! Array to store temporary data
                                                   !  (to allow to store the difference between 2 places)
  real, dimension(:,:,:), allocatable :: xtmplesstore ! Array to store temporary data for LES budgets
                                                      !  (to allow to store the difference between 2 places)
  type(tbusourcedata), dimension(:), allocatable :: tsources ! Full list of source terms (used or not)
  type(tbugroupdata),  dimension(:), allocatable :: tgroups  ! Full list of groups of source terms (to be written)
  type(tburhodata),    pointer                   :: trhodj => null() ! Budget array for rhodj
end type tbudgetdata

TYPE TBUDGETCONF_t
 LOGICAL :: LBU_ENABLE=.FALSE.
 LOGICAL :: LBUDGET_U=.FALSE.  ! flag to compute budget of RhoJu  and/or LES budgets with u
 LOGICAL :: LBUDGET_V=.FALSE.  ! flag to compute budget of RhoJv  and/or LES budgets with u
 LOGICAL :: LBUDGET_W=.FALSE.  ! flag to compute budget of RhoJw  and/or LES budgets with u
 LOGICAL :: LBUDGET_TH=.FALSE. ! flag to compute budget of RhoJTh and/or LES budgets with th
 LOGICAL :: LBUDGET_TKE=.FALSE.! flag to compute budget of RhoJTke and/or LES budgets with Tke
 LOGICAL :: LBUDGET_RV=.FALSE. ! flag to compute budget of RhoJrv and/or LES budgets with rv
 LOGICAL :: LBUDGET_RC=.FALSE. ! flag to compute budget of RhoJrc and/or LES budgets with rc
 LOGICAL :: LBUDGET_RR=.FALSE. ! flag to compute budget of RhoJrr and/or LES budgets with rr
 LOGICAL :: LBUDGET_RI=.FALSE. ! flag to compute budget of RhoJri and/or LES budgets with ri
 LOGICAL :: LBUDGET_RS=.FALSE. ! flag to compute budget of RhoJrs and/or LES budgets with rs
 LOGICAL :: LBUDGET_RG=.FALSE. ! flag to compute budget of RhoJrg and/or LES budgets with rg
 LOGICAL :: LBUDGET_RH=.FALSE. ! flag to compute budget of RhoJrh and/or LES budgets with rh
 LOGICAL :: LBUDGET_SV=.FALSE. ! flag to compute budget of RhoJsv and/or LES budgets with sv
END TYPE TBUDGETCONF_t
!
TYPE(TBUDGETCONF_t), TARGET :: TBUCONF
!
type(tbudgetdata), dimension(:), allocatable, save :: tbudgets
type(tburhodata),                pointer,     save :: tburhodj => null() ! Budget array for rhodj used inside some tbudgets


!                       General variables
LOGICAL, POINTER :: LBU_ENABLE=>NULL()
!
CHARACTER (LEN=4), SAVE :: CBUTYPE         ! type of desired budget 'CART'
                                           ! (cartesian box) or 'MASK' (budget
                                           ! zone defined by a mask) or 'NONE'
                                           ! (no budget)
INTEGER, SAVE :: NBUMOD                    ! model in which budget is 
                                           ! calculated
!
LOGICAL, SAVE :: LBU_BEG                   ! switch for budget beginning
!
REAL, SAVE    :: XBULEN                    ! length in seconds of the budget 
                                           ! temporal average
!
INTEGER, SAVE :: NBUSTEP                   ! number of model timesteps required 
                                           ! for the budget time average
REAL, SAVE    :: XBUWRI                    ! period in seconds between
                                           ! budget writing for budget masks
INTEGER, SAVE :: NBUTSHIFT                 ! temporal shift for budgets writing
integer, save :: nbusubwrite = 0           ! Number of budget time average periods for each write
integer, save :: nbutotwrite = 0           ! Total number of budget time average periods
!
INTEGER, SAVE :: NBUKL, NBUKH              ! lowest and highest K indice values 
                                           ! of the budget box in the physical domain
LOGICAL, SAVE :: LBU_KCP                   ! switch for compression in K
                                           ! direction
!
!                Variables used by the cartesian box case ('CART') only
!
INTEGER, SAVE :: NBUIL, NBUIH              ! lowest and highest I indice values 
                                           ! of the cartesian box in the physical domain
INTEGER, SAVE :: NBUJL, NBUJH              ! lowest and highest J indice values 
                                           ! of the cartesian box in the physical domain
LOGICAL, SAVE :: LBU_ICP                   ! switch for compression in I
                                           ! direction
LOGICAL, SAVE :: LBU_JCP                   ! switch for comppression in J
                                           ! direction
!
!                Variables used by the  mask case ('MASK') only
!
INTEGER, SAVE :: NBUMASK                   ! number of MASK zones for which 
                                           ! budgets are performed 
LOGICAL, SAVE, DIMENSION(:,:,:),         & ! define the zone where the MASK 
           ALLOCATABLE :: LBU_MASK         ! is True 
!                                          
INTEGER, SAVE, DIMENSION(:,:,:,:),       & ! surface for each mask at each
           ALLOCATABLE :: NBUSURF          ! budget step
!             
INTEGER, SAVE :: NBUTIME                   ! number of budget time periods
!
!                       Variables for budget storage 
!
!                       General variables
INTEGER, SAVE :: NBUSIL, NBUSIH      ! lowest and highest I indices of the intersection
                                     ! of the cartesian box with the sub-domain
INTEGER, SAVE :: NBUSJL, NBUSJH      ! lowest and highest J indices of the intersection
                                     ! of the global cartesian box
INTEGER, SAVE :: NBUIMAX_ll                ! second dimension of the budget
INTEGER, SAVE :: NBUJMAX_ll                ! second dimension of the budget
                                           ! array in the global domain (in CART case)
!
INTEGER, SAVE :: NBUIMAX                   ! first dimension of the budget
                                           ! tabular
INTEGER, SAVE :: NBUJMAX                   ! second dimension of the budget
                                           ! tabular
INTEGER, SAVE :: NBUKMAX                   ! dimension along K of the budget
                                           ! tabular
!
!      Allowed processes for the budget of the x scalar variables
!        (transport part only)
!
! For each budget, the switches values for budgets
! activation may be set by the user in a namelist. Their default value is 0.
! In the following declaration, the corresponding process names  are given 
! beside as comments.
!     
!      Allowed processes for the budget of RU (wind component along x)
!
! Current namelist: NAM_BU_RU
!
LOGICAL, SAVE :: LBU_RU = .FALSE. ! True when the budget of RU is performed
!                         
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RU
!
!      Allowed processes for the budget of RV (wind component along y)
!                                                  
! Current namelist: NAM_BU_RV
!
LOGICAL, SAVE :: LBU_RV = .FALSE. ! True when the budget of RV is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RV
!
!      Allowed processes for the budget of RW (wind vertical component)
!                                                  
! Current namelist: NAM_BU_RW
!
LOGICAL, SAVE :: LBU_RW = .FALSE. ! True when the budget of RW is performed
!                                                  
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RW
!
!      Allowed processes for the budget of RTH (potential temperature)
!                                                  
! Current namelist: NAM_BU_RTH
!
LOGICAL, SAVE :: LBU_RTH = .FALSE. ! True when the budget of RTH is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RTH
!
!      Allowed processes for the budget of RTKE (kinetic energy)
!                                                  
! Current namelist: NAM_BU_RTKE
!
LOGICAL, SAVE :: LBU_RTKE = .FALSE. ! True when the budget of RTKE is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RTKE
!
!      Allowed processes for the budget of moist variable RRV (water vapor)
!                                                  
! Current namelist: NAM_BU_RRV
!
LOGICAL, SAVE :: LBU_RRV = .FALSE. ! true when the budget of RRV is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RRV
!
!      Allowed processes for the budget of moist variable RRC (cloud water)
!                                                  
! Current namelist: NAM_BU_RRC
!
LOGICAL, SAVE :: LBU_RRC = .FALSE. ! True when the budget of RRC is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RRC
!
!      Allowed processes for the budget of moist variable RRR (rain water)
!
! Current namelist: NAM_BU_RRR
!
LOGICAL, SAVE :: LBU_RRR = .FALSE. ! True when the budget of RRR is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RRR
!
!      Allowed processes for the budget of moist variable RRI (ice)
!
! Current namelist: NAM_BU_RRI
!
LOGICAL, SAVE :: LBU_RRI = .FALSE. ! True when the budget of RRI is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RRI
!
!      Allowed processes for the budget of moist variable RRS (snow)
!
! Current namelist: NAM_BU_RRS
!
LOGICAL, SAVE :: LBU_RRS = .FALSE. ! True when the budget of RRS is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RRS
!
!      Allowed processes for the budget of moist variable RRG (graupel)
!
! Current namelist: NAM_BU_RRG
!
LOGICAL, SAVE :: LBU_RRG = .FALSE. ! True when the budget of RRG is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RRG
!
!      Allowed processes for the budget of moist variable RRH (hail)
!
! Current namelist: NAM_BU_RRH
!
LOGICAL, SAVE :: LBU_RRH = .FALSE. ! True when the budget of RRH is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RRH
!
! Current namelist: NAM_BU_RSV
!
LOGICAL, SAVE :: LBU_RSV = .FALSE. ! True when the budget of RSVx is performed
!
CHARACTER(LEN=NBULISTMAXLEN), DIMENSION(:), ALLOCATABLE :: CBULIST_RSV
!
!
REAL :: XTIME_BU          ! budget time in this time-step
REAL :: XTIME_BU_PROCESS  ! budget time per process for this time-step
!
LOGICAL, POINTER :: LBUDGET_U=>NULL() ! flag to compute budget of RhoJu  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_V=>NULL()  ! flag to compute budget of RhoJv  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_W=>NULL()  ! flag to compute budget of RhoJw  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_TH=>NULL() ! flag to compute budget of RhoJTh and/or LES budgets with th
LOGICAL, POINTER :: LBUDGET_TKE=>NULL() ! flag to compute budget of RhoJTke and/or LES budgets with Tke
LOGICAL, POINTER :: LBUDGET_RV=>NULL() ! flag to compute budget of RhoJrv and/or LES budgets with rv
LOGICAL, POINTER :: LBUDGET_RC=>NULL() ! flag to compute budget of RhoJrc and/or LES budgets with rc
LOGICAL, POINTER :: LBUDGET_RR=>NULL() ! flag to compute budget of RhoJrr and/or LES budgets with rr
LOGICAL, POINTER :: LBUDGET_RI=>NULL() ! flag to compute budget of RhoJri and/or LES budgets with ri
LOGICAL, POINTER :: LBUDGET_RS=>NULL() ! flag to compute budget of RhoJrs and/or LES budgets with rs
LOGICAL, POINTER :: LBUDGET_RG=>NULL() ! flag to compute budget of RhoJrg and/or LES budgets with rg
LOGICAL, POINTER :: LBUDGET_RH=>NULL() ! flag to compute budget of RhoJrh and/or LES budgets with rh
LOGICAL, POINTER :: LBUDGET_SV=>NULL() ! flag to compute budget of RhoJsv and/or LES budgets with sv
!
CONTAINS
SUBROUTINE TBUCONF_ASSOCIATE()
  IMPLICIT NONE
  LBU_ENABLE=>TBUCONF%LBU_ENABLE

  LBUDGET_U=>TBUCONF%LBUDGET_U
  LBUDGET_V=>TBUCONF%LBUDGET_V
  LBUDGET_W=>TBUCONF%LBUDGET_W
  LBUDGET_TH=>TBUCONF%LBUDGET_TH
  LBUDGET_TKE=>TBUCONF%LBUDGET_TKE
  LBUDGET_RV=>TBUCONF%LBUDGET_RV
  LBUDGET_RC=>TBUCONF%LBUDGET_RC
  LBUDGET_RR=>TBUCONF%LBUDGET_RR
  LBUDGET_RI=>TBUCONF%LBUDGET_RI
  LBUDGET_RS=>TBUCONF%LBUDGET_RS
  LBUDGET_RG=>TBUCONF%LBUDGET_RG
  LBUDGET_RH=>TBUCONF%LBUDGET_RH
  LBUDGET_SV=>TBUCONF%LBUDGET_SV
END SUBROUTINE TBUCONF_ASSOCIATE

END MODULE MODD_BUDGET
