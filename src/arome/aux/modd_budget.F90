!     ######spl
      MODULE MODD_BUDGET
!     ##################
!
!!****  *MODD_BUDGET* - declaration of budget variables
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to specify  the budget
!     variables.
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_PARAMETERS: JPBUMAX, JPBUPROCMAX
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (module MODD_BUDGET)
!!
!!    AUTHOR
!!    ------
!!      P. Hereil   *Meteo France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        23/02/95
!!      J.-P. Lafore    10/02/98    adding of rhodj declaration for budget
!!      V. Ducrocq      4/06/99     //
!!      J.-P. Pinty     25/09/00    additional budget terms for C2R2 scheme
!!      D. Gazen        22/01/01    add NCHEMSV
!!      V. Masson       06/11/02    new flags for budget calls and time counters
!!      V. Masson       27/11/02    add 2way nesting effect
!!      P. Jabouille    07/07/04    add budget terms for microphysics
!!      C. Barthe       19/11/09    add budget terms for electricity
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
USE MODD_PARAMETERS, ONLY :JPBUMAX, JPBUPROMAX
USE DDH_MIX, ONLY : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
IMPLICIT NONE

SAVE
!
INTEGER, PARAMETER:: NBUDGET_RHO = 0  ! Reference number for budget of RhoJ
INTEGER, PARAMETER:: NBUDGET_U   = 1  ! Reference number for budget of RhoJu  and/or LES budgets with u
INTEGER, PARAMETER:: NBUDGET_V   = 2  ! Reference number for budget of RhoJv  and/or LES budgets with u
INTEGER, PARAMETER:: NBUDGET_W   = 3  ! Reference number for budget of RhoJw  and/or LES budgets with u
INTEGER, PARAMETER:: NBUDGET_TH  = 4  ! Reference number for budget of RhoJTh and/or LES budgets with th
INTEGER, PARAMETER:: NBUDGET_TKE = 5  ! Reference number for budget of RhoJTke and/or LES budgets with Tke
INTEGER, PARAMETER:: NBUDGET_RV  = 6  ! Reference number for budget of RhoJrv and/or LES budgets with rv
INTEGER, PARAMETER:: NBUDGET_RC  = 7  ! Reference number for budget of RhoJrc and/or LES budgets with rc
INTEGER, PARAMETER:: NBUDGET_RR  = 8  ! Reference number for budget of RhoJrr and/or LES budgets with rr
INTEGER, PARAMETER:: NBUDGET_RI  = 9  ! Reference number for budget of RhoJri and/or LES budgets with ri
INTEGER, PARAMETER:: NBUDGET_RS  = 10 ! Reference number for budget of RhoJrs and/or LES budgets with rs
INTEGER, PARAMETER:: NBUDGET_RG  = 11 ! Reference number for budget of RhoJrg and/or LES budgets with rg
INTEGER, PARAMETER:: NBUDGET_RH  = 12 ! Reference number for budget of RhoJrh and/or LES budgets with rh
INTEGER, PARAMETER:: NBUDGET_SV1 = 13 ! Reference number for 1st budget of RhoJsv and/or LES budgets with sv
!
TYPE TBUDGETDATA
  INTEGER :: NBUDGET
  TYPE(TYP_DDH), POINTER :: YDDDH=>NULL()
  TYPE(TLDDH), POINTER :: YDLDDH=>NULL()
  TYPE(TMDDH), POINTER :: YDMDDH=>NULL()
ENDTYPE TBUDGETDATA
!
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
!                       General variables
LOGICAL, POINTER :: LBU_ENABLE=>NULL()
!
INTEGER, DIMENSION(JPBUMAX,JPBUPROMAX)  &  ! number of processes to be
                             :: NBUINC=0           ! avoided for every budget
                                                 ! between one active
                                                 ! source to the next one
INTEGER, DIMENSION(JPBUMAX)             &  ! counter for all the processes
                             :: NBUCTR_ACTV=0      ! activated or not
!
CHARACTER (LEN=4) :: CBUTYPE         ! type of desired budget 'CART'
                                           ! (cartesian box) or 'MASK' (budget
                                           ! zone defined by a mask) or 'NONE'
                                           ! (no budget)
INTEGER :: NBUMOD=0                    ! model in which budget is
                                           ! calculated
INTEGER, DIMENSION(:),             & ! number of processes for each
                 ALLOCATABLE :: NBUPROCNBR ! budget
!
INTEGER, DIMENSION(:),             & ! process counter linked to each
                 ALLOCATABLE :: NBUPROCCTR ! budget
!
CHARACTER(LEN=2), DIMENSION(:,:),  & ! resulting string character of the
        ALLOCATABLE :: CBUACTION           ! transcription of the budget actions
                                           ! (integer) read in  namelists or
                                           ! set by default
CHARACTER (LEN=16), DIMENSION(:,:),& ! names of records on the FM file
                 ALLOCATABLE :: CBURECORD  ! for the budgets
!
CHARACTER (LEN=99), DIMENSION(:,:),& ! name of a process for a budget. It
                 ALLOCATABLE :: CBUCOMMENT ! will appear in the comment part of
                                           ! the previous record
!
LOGICAL :: LBU_BEG=.FALSE.           ! switch for budget beginning
!
REAL    :: XBULEN=0.                    ! length in seconds of the budget
                                           ! temporal average
!
INTEGER :: NBUSTEP=0                   ! number of model timesteps required
                                           ! for the budget time average
REAL    :: XBUWRI=0.                       ! period in seconds of
                                           ! budget writing on FM-files
INTEGER :: NBUWRNB=0                   ! number of budget periods when storage
                                           ! arrays are written on FM-files
INTEGER :: NBUTSHIFT=0                 ! temporal shift for budgets writing
!
INTEGER :: NBUKH=0              ! lowest and highest K indice values
INTEGER :: NBUKL=0              ! lowest and highest K indice values
                                           ! of the budget box
LOGICAL :: LBU_KCP=.FALSE.           ! switch for compression in K
                                           ! direction
!
!                Variables used by the cartesian box case ('CART') only
!
INTEGER :: NBUIL=0              ! lowest and highest I indice values
INTEGER :: NBUIH=0              ! lowest and highest I indice values
                                           ! of the cartesian box
INTEGER :: NBUJL=0              ! lowest and highest J indice values
INTEGER :: NBUJH=0              ! lowest and highest J indice values
                                           ! of the cartesian box
LOGICAL :: LBU_ICP=.FALSE.           ! switch for compression in I
                                           ! direction
LOGICAL :: LBU_JCP=.FALSE.           ! switch for comppression in J
                                           ! direction
!
!                Variables used by the  mask case ('MASK') only
!
INTEGER :: NBUMASK=0                   ! number of MASK zones for which
                                           ! budgets are performed
LOGICAL, DIMENSION(:,:,:),         & ! define the zone where the MASK
           ALLOCATABLE :: LBU_MASK         ! is True
!
REAL, DIMENSION(:,:,:,:),          & ! surface for each mask at each
           ALLOCATABLE :: XBUSURF          ! budget step
!
INTEGER :: NBUTIME=0                   ! number of budget time periods
!
!                       Variables for budget storage
!
!                       General variables
INTEGER :: NBUSIL=0      ! lowest and highest I indices of the intersection
INTEGER :: NBUSIH=0      ! lowest and highest I indices of the intersection
                                     ! of the cartesian box with the sub-domain
INTEGER :: NBUSJL=0      ! lowest and highest J indices of the intersection
INTEGER :: NBUSJH=0      ! lowest and highest J indices of the intersection
                                     ! of the global cartesian box
INTEGER :: NBUIMAX_ll=0                ! second dimension of the budget
INTEGER :: NBUJMAX_ll=0                ! second dimension of the budget
                                           ! array in the global domain (in CART case)
!
INTEGER :: NBUIMAX=0                   ! first dimension of the budget
                                           ! tabular
INTEGER :: NBUJMAX=0                   ! second dimension of the budget
                                           ! tabular
INTEGER :: NBUKMAX=0                   ! dimension along K of the budget
                                           ! tabular
REAL, DIMENSION(:,:,:,:),          & ! budget arrays for RU, RV and
        ALLOCATABLE :: XBURU, XBURV, XBURW ! RW (wind components) respectively
REAL, DIMENSION(:,:,:,:),          & ! budget arrays for RTH (potential
        ALLOCATABLE :: XBURTH, XBURTKE     ! temperature) and RTKE (kinetic
                                           ! energy)
REAL, DIMENSION(:,:,:,:),          & ! budget arrays for RRV (water vapor)
        ALLOCATABLE :: XBURRV, XBURRC      ! and RRC (cloud water)
REAL, DIMENSION(:,:,:,:),          & ! budget arrays for RRR (rain water)
        ALLOCATABLE :: XBURRR, XBURRI      ! and RRI (ice)
REAL, DIMENSION(:,:,:,:),          & ! budget arrays for RRS (snow)
        ALLOCATABLE :: XBURRS, XBURRG      ! and RRG (graupel)
REAL, DIMENSION(:,:,:,:),          & ! budget array for RRH (hail)
        ALLOCATABLE :: XBURRH              !
REAL, DIMENSION(:,:,:,:,:), &
                 ALLOCATABLE :: XBURSV       ! Budget of the SVx
REAL, DIMENSION(:,:,:),            & ! budget arrays for RHODJ at
               ALLOCATABLE :: XBURHODJ , & !   scalar localization
                              XBURHODJU, & !        U localization
                              XBURHODJV, & !        V localization
                              XBURHODJW    !    and W localization
!
!      Allowed processes for the budget of the x scalar variables
!        (transport part only)
!
! For each budget, the switches values (from 0 to JPBUPROMAX) for budgets
! activation may be set by the user in a namelist. Their default value is 0.
! In the following declaration, the corresponding process names  are given
! beside as comments.
!
!      Allowed processes for the budget of RU (wind component along x)
!
! Courant namelist: NAM_BURU
!
LOGICAL :: LBU_RU=.FALSE.     ! True when the budget of RU is performed
!
INTEGER :: NASSEU=0     ! time filter
INTEGER :: NNESTU=0     ! Efffect of 2way nesting on U
INTEGER :: NADVXU=0     ! advection along X
INTEGER :: NADVYU=0     ! advection along Y
INTEGER :: NADVZU=0     ! advection along Z
INTEGER :: NFRCU=0      ! forcing
INTEGER :: NNUDU=0      ! nudging
INTEGER :: NCURVU=0     ! curvature
INTEGER :: NCORU=0      ! Coriolis terms
INTEGER :: NDIFU=0      ! numerical diffusion
INTEGER :: NRELU=0      ! relaxation
INTEGER :: NHTURBU=0    ! horizontal TURBulence
INTEGER :: NVTURBU=0    ! vertical turbulence
INTEGER :: NDRAGU=0     ! vegetation drag
INTEGER :: NMAFLU=0     ! mass flux
INTEGER :: NPRESU=0     ! pressure term
!
!      Allowed processes for the budget of RV (wind component along y)
!
! Courant namelist: NAM_BURV
!
LOGICAL :: LBU_RV=.FALSE.     ! True when the budget of RV is performed
!
INTEGER :: NASSEV=0     ! time filter
INTEGER :: NNESTV=0     ! Efffect of 2way nesting on V
INTEGER :: NADVXV=0     ! advection along X
INTEGER :: NADVYV=0     ! advection along Y
INTEGER :: NADVZV=0     ! advection along Z
INTEGER :: NFRCV=0      ! forcing
INTEGER :: NNUDV=0      ! nudging
INTEGER :: NCURVV=0     ! curvature
INTEGER :: NCORV=0      ! Coriolis terms
INTEGER :: NDIFV=0      ! numerical diffusion
INTEGER :: NRELV=0      ! relaxation
INTEGER :: NHTURBV=0    ! horizontal turbulence
INTEGER :: NVTURBV=0    ! vertical turbulence
INTEGER :: NDRAGV=0     ! vegetation drag
INTEGER :: NMAFLV=0     ! mass flux
INTEGER :: NPRESV=0     ! pressure term
!
!      Allowed processes for the budget of RW (wind vertical component)
!
! Courant namelist: NAM_BURW
!
LOGICAL :: LBU_RW=.FALSE.     ! True when the budget of RW is performed
!
INTEGER :: NASSEW=0     ! time filter
INTEGER :: NNESTW=0     ! Efffect of 2way nesting on W
INTEGER :: NADVXW=0     ! advection along X
INTEGER :: NADVYW=0     ! advection along Y
INTEGER :: NADVZW=0     ! advection along Z
INTEGER :: NFRCW=0      ! forcing
INTEGER :: NNUDW=0      ! nudging
INTEGER :: NCURVW=0     ! curvature
INTEGER :: NCORW=0      ! Coriolis terms
INTEGER :: NGRAVW=0     ! gravity term
INTEGER :: NDIFW=0      ! numerical diffusion
INTEGER :: NRELW=0      ! relaxation
INTEGER :: NHTURBW=0    ! horizontal turbulence
INTEGER :: NVTURBW=0    ! vertical turbulence
INTEGER :: NPRESW=0     ! pressure term
!
!      Allowed processes for the budget of RTH (potential temperature)
!
! Courant namelist: NAM_BURTH
!
LOGICAL :: LBU_RTH=.FALSE.    ! True when the budget of RTH is performed
!
INTEGER :: NASSETH=0    ! time filter
INTEGER :: NNESTTH=0    ! Efffect of 2way nesting on Th
INTEGER :: NADVTH=0     ! Total advection for PPM
INTEGER :: NADVXTH=0    ! advection along X (all except PPM)
INTEGER :: NADVYTH=0    ! advection along Y (all except PPM)
INTEGER :: NADVZTH=0    ! advection along Z (all except PPM)
INTEGER :: NFRCTH=0     ! forcing
INTEGER :: N2DADVTH=0  ! 2d advecting forcing
INTEGER :: N2DRELTH=0   ! 2d relaxation forcing
INTEGER :: NNUDTH=0     ! nudging
INTEGER :: NPREFTH=0    ! theta source term due to the reference pressure
                            ! (Dyn. Sources) only present if KRR>0
INTEGER :: NDIFTH=0     ! numerical diffusion
INTEGER :: NRELTH=0     ! relaxation
INTEGER :: NRADTH=0     ! RADiation
INTEGER :: NDCONVTH=0   ! KAFR CONVection
INTEGER :: NMAFLTH=0    ! Mass flux
INTEGER :: NHTURBTH=0   ! horizontal turbulence
INTEGER :: NVTURBTH=0   ! vertical turbulence
INTEGER :: NDISSHTH=0   ! dissipative heating
INTEGER :: NNEGATH=0    ! negative correction induced by hydrometeors
INTEGER :: NREVATH=0    ! rain evaporation
INTEGER :: NCONDTH=0    ! evaporation/condensation
INTEGER :: NHENUTH=0    ! HEterogenous NUcleation ICE3
INTEGER :: NHONTH=0     ! HOmogeneous Nucleation  ICE3
INTEGER :: NSFRTH=0     ! Spontaneous FReezing    ICE3
INTEGER :: NDEPSTH=0    ! DEPosition on Snow      ICE3
INTEGER :: NDEPGTH=0    ! DEPosition on Graupel   ICE3
INTEGER :: NRIMTH=0     ! RIMing of cloudwater    ICE3
INTEGER :: NACCTH=0     ! ACCretion of rainwater  ICE3
INTEGER :: NCFRZTH=0    ! Conversion FReeZing     ICE3
INTEGER :: NWETGTH=0    ! WET Growth of graupel   ICE3
INTEGER :: NDRYGTH=0    ! DRY Growth of graupel   ICE3
INTEGER :: NGMLTTH=0    ! Graupel MeLTing         ICE3
INTEGER :: NIMLTTH=0    ! Ice MeLTing             ICE3
INTEGER :: NBERFITH=0   ! BERgeron-FIndeisen gth. ICE3
INTEGER :: NCDEPITH=0   ! Cond./DEPosition on ice ICE3
INTEGER :: NWETHTH=0    ! wet growth of hail      ICE4
INTEGER :: NHMLTTH=0    ! melting of hail         ICE4
INTEGER :: NHINDTH=0    ! Heterogeneous Nucleation by Deposition C3R5
INTEGER :: NHINCTH=0    ! Heterogeneous Nucleation by Contact    C3R5
INTEGER :: NHONHTH=0    ! Haze Homogeneous Nucleation            C3R5
INTEGER :: NHONCTH=0    ! droplet homogeneous nucleation         C3R5
INTEGER :: NHONRTH=0    ! drop homogeneous nucleation            C3R5
INTEGER :: NCEDSTH=0    ! adjustment                             C3R5
!
!      Allowed processes for the budget of RTKE (kinetic energy)
!
! Courant namelist: NAM_BURTKE
!
LOGICAL :: LBU_RTKE=.FALSE.   ! True when the budget of RTKE is performed
!
INTEGER :: NASSETKE=0   ! time filter
INTEGER :: NADVTKE=0    ! Total advection for PPM
INTEGER :: NADVXTKE=0   ! advection along X (all except PPM)
INTEGER :: NADVYTKE=0   ! advection along Y (all except PPM)
INTEGER :: NADVZTKE=0   ! advection along Z (all except PPM)
INTEGER :: NFRCTKE=0    ! forcing
INTEGER :: NDIFTKE=0    ! numerical diffusion
INTEGER :: NRELTKE=0    ! relaxation
INTEGER :: NDPTKE=0     ! dynamic production of TKE
INTEGER :: NTPTKE=0     ! thermal production of TKE
INTEGER :: NDRAGTKE=0   ! vegetation drag
INTEGER :: NDISSTKE=0   ! dissipation of TKE
INTEGER :: NTRTKE=0     ! turbulent transport of TKE
!
!
!      Allowed processes for the budget of moist variable RRV (water vapor)
!
! Courant namelist: NAM_BURRV
!
LOGICAL :: LBU_RRV=.FALSE.   ! true when the budget of RRV is performed
!
INTEGER :: NASSERV=0   ! time filter
INTEGER :: NNESTRV=0   ! Effect of 2way nesting on Rv
INTEGER :: NADVRV=0    ! Total advection for PPM
INTEGER :: NADVXRV=0   ! advection along X (all except PPM)
INTEGER :: NADVYRV=0   ! advection along Y (all except PPM)
INTEGER :: NADVZRV=0   ! advection along Z (all except PPM)
INTEGER :: NFRCRV=0    ! forcing
INTEGER :: N2DADVRV=0  ! 2d advecting forcing
INTEGER :: N2DRELRV=0  ! 2d relaxation forcing
INTEGER :: NNUDRV=0    ! nudging
INTEGER :: NDIFRV=0    ! numerical diffusion
INTEGER :: NRELRV=0    ! relaxation
INTEGER :: NDCONVRV=0  ! KAFR CONVection
INTEGER :: NMAFLRV=0   ! Mass flux
INTEGER :: NHTURBRV=0  ! horizontal turbulence
INTEGER :: NVTURBRV=0  ! vertical turbulence
INTEGER :: NNEGARV=0   ! negative correction
INTEGER :: NREVARV=0   ! rain evaporation
INTEGER :: NCONDRV=0   ! evaporation/condensation
INTEGER :: NHENURV=0   ! HEterogenous NUcleation ICE3
INTEGER :: NDEPSRV=0   ! DEPosition on Snow      ICE3
INTEGER :: NDEPGRV=0   ! DEPosition on Graupel   ICE3
INTEGER :: NCDEPIRV=0  ! Cond./DEPosition on ice ICE3
INTEGER :: NHINDRV=0   ! Heterogeneous Nucleation by Deposition C3R5
INTEGER :: NHONHRV=0   ! Haze Homogeneous Nucleation            C3R5
INTEGER :: NCEDSRV=0   ! adjustement                            C3R5
!
!      Allowed processes for the budget of moist variable RRC (cloud water)
!
! Courant namelist: NAM_BURRC
!
LOGICAL :: LBU_RRC=.FALSE.    ! True when the budget of RRC is performed
!
INTEGER :: NASSERC=0    ! time filter
INTEGER :: NNESTRC=0    ! Efffect of 2way nesting on Rc
INTEGER :: NADVRC=0     ! Total advection for PPM
INTEGER :: NADVXRC=0    ! advection along X (all except PPM)
INTEGER :: NADVYRC=0    ! advection along Y (all except PPM)
INTEGER :: NADVZRC=0    ! advection along Z (all except PPM)
INTEGER :: NFRCRC=0     ! forcing
INTEGER :: NDIFRC=0     ! numerical diffusion
INTEGER :: NRELRC=0     ! relaxation
INTEGER :: NDCONVRC=0   ! Deep CONVection
INTEGER :: NHTURBRC=0   ! horizontal turbulence
INTEGER :: NVTURBRC=0   ! vertical turbulence
INTEGER :: NNEGARC=0    ! negative correction
INTEGER :: NACCRRC=0    ! accretion
INTEGER :: NAUTORC=0    ! autoconversion
INTEGER :: NCONDRC=0    ! evaporation/condensation
INTEGER :: NHONRC=0     ! HOmogeneous Nucleation  ICE3
INTEGER :: NRIMRC=0     ! RIMing of cloudwater    ICE3
INTEGER :: NWETGRC=0    ! WET Growth of graupel   ICE3
INTEGER :: NDRYGRC=0    ! DRY Growth of graupel   ICE3
INTEGER :: NIMLTRC=0    ! Ice MeLTing             ICE3
INTEGER :: NBERFIRC=0   ! BERgeron-FIndeisen gth. ICE3
INTEGER :: NCDEPIRC=0   ! Cond./DEPosition on ice ICE3
INTEGER :: NHENURC=0    ! CCN Activation C2R2
INTEGER :: NSEDIRC=0    ! sedimentation  C2R2
INTEGER :: NWETHRC=0    ! wet growth of hail
INTEGER :: NHINCRC=0    ! Heterogeneous Nucleation by Contact C3R5
INTEGER :: NHONCRC=0    ! droplet homogeneous nucleation      C3R5
INTEGER :: NCEDSRC=0    ! adjustment                          C3R5
INTEGER :: NREVARC=0    ! evaporation of rain drops
INTEGER :: NDEPORC=0    ! ground deposition     
INTEGER :: NDEPOTRRC=0  ! deposition on tree
!
!      Allowed processes for the budget of moist variable RRR (rain water)
!
! Courant namelist: NAM_BURRR
!
LOGICAL :: LBU_RRR=.FALSE.    ! True when the budget of RRR is performed
!
INTEGER :: NASSERR=0    ! time filter
INTEGER :: NNESTRR=0    ! Efffect of 2way nesting on Rr
INTEGER :: NADVRR=0     ! Total advection for PPM
INTEGER :: NADVXRR=0    ! advection along X (all except PPM)
INTEGER :: NADVYRR=0    ! advection along Y (all except PPM)
INTEGER :: NADVZRR=0    ! advection along Z (all except PPM)
INTEGER :: NFRCRR=0     ! forcing
INTEGER :: NDIFRR=0     ! numerical diffusion
INTEGER :: NRELRR=0     ! relaxation
INTEGER :: NNEGARR=0    ! negative correction
INTEGER :: NACCRRR=0    ! accretion
INTEGER :: NAUTORR=0    ! autoconversion
INTEGER :: NREVARR=0    ! rain evaporation
INTEGER :: NSEDIRR=0    ! sedimentation
INTEGER :: NSFRRR=0     ! Spontaneous FReezing    ICE3
INTEGER :: NACCRR=0     ! ACCretion of rainwater  ICE3
INTEGER :: NCFRZRR=0    ! Conversion FReeZing     ICE3
INTEGER :: NWETGRR=0    ! WET Growth of graupel   ICE3
INTEGER :: NDRYGRR=0    ! DRY Growth of graupel   ICE3
INTEGER :: NGMLTRR=0    ! Graupel MeLTing         ICE3
INTEGER :: NWETHRR=0    ! wet growth of hail      ICE4
INTEGER :: NHMLTRR=0    ! melting of hail         ICE4
INTEGER :: NHONRRR=0    ! drop homogeneous nucleation C3R5
!
!      Allowed processes for the budget of moist variable RRI (ice)
!
! Courant namelist: NAM_BURRI
!
LOGICAL :: LBU_RRI=.FALSE.    ! True when the budget of RRI is performed
!
INTEGER :: NASSERI=0    ! time filter
INTEGER :: NNESTRI=0    ! Efffect of 2way nesting on Ri
INTEGER :: NADVRI=0     ! Total advection for PPM
INTEGER :: NADVXRI=0    ! advection along X (all except PPM)
INTEGER :: NADVYRI=0    ! advection along Y (all except PPM)
INTEGER :: NADVZRI=0    ! advection along Z (all except PPM)
INTEGER :: NFRCRI=0     ! forcing
INTEGER :: NDIFRI=0     ! numerical diffusion
INTEGER :: NRELRI=0     ! relaxation
INTEGER :: NDCONVRI=0   ! Deep CONVection
INTEGER :: NHTURBRI=0   ! horizontal turbulence
INTEGER :: NVTURBRI=0   ! vertical turbulence
INTEGER :: NNEGARI=0    ! negative correction
INTEGER :: NSEDIRI=0    ! SEDImentation           ICE3
INTEGER :: NHENURI=0    ! HEterogenous NUcleation ICE3
INTEGER :: NHONRI=0     ! HOmogeneous Nucleation  ICE3
INTEGER :: NAGGSRI=0    ! AGGregation of snow     ICE3
INTEGER :: NAUTSRI=0    ! AUToconversion of ice   ICE3
INTEGER :: NCFRZRI=0    ! Conversion FReeZing     ICE3
INTEGER :: NWETGRI=0    ! WET Growth of graupel   ICE3
INTEGER :: NDRYGRI=0    ! DRY Growth of graupel   ICE3
INTEGER :: NIMLTRI=0    ! Ice MeLTing             ICE3
INTEGER :: NBERFIRI=0   ! BERgeron-FIndeisen gth. ICE3
INTEGER :: NCDEPIRI=0   ! Cond./DEPosition on ice ICE3
INTEGER :: NWETHRI=0    ! wet growth of hail      ICE4
INTEGER :: NHINDRI=0 ! heterogeneous nucleation by deposition C3R5
INTEGER :: NHINCRI=0 ! heterogeneous nucleation by contact    C3R5
INTEGER :: NHONHRI=0 ! haze homogeneous nucleation source     C3R5
INTEGER :: NHONCRI=0 ! droplet homogeneous nucleation         C3R5
INTEGER :: NCNVIRI=0 ! Conversion of snow to r_i              C3R5
INTEGER :: NCNVSRI=0 ! Conversion of pristine ice to r_s      C3R5
INTEGER :: NHMSRI=0  ! Hallett-Mossop ice multiplication process due to snow riming C3R5
INTEGER :: NHMGRI=0  ! Hallett-Mossop ice multiplication process due to graupel riming C3R5
INTEGER :: NCEDSRI=0 ! adjustement                            C3R5
!
!      Allowed processes for the budget of moist variable RRS (snow)
!
! Courant namelist: NAM_BURRS
!
LOGICAL :: LBU_RRS=.FALSE.    ! True when the budget of RRS is performed
!
INTEGER :: NASSERS=0    ! time filter
INTEGER :: NNESTRS=0    ! Efffect of 2way nesting on Rs
INTEGER :: NADVRS=0     ! Total advection for PPM
INTEGER :: NADVXRS=0    ! advection along X (all except PPM)
INTEGER :: NADVYRS=0    ! advection along Y (all except PPM)
INTEGER :: NADVZRS=0    ! advection along Z (all except PPM)
INTEGER :: NFRCRS=0     ! forcing
INTEGER :: NDIFRS=0     ! numerical diffusion
INTEGER :: NRELRS=0     ! relaxation
INTEGER :: NNEGARS=0    ! negative correction
INTEGER :: NSEDIRS=0    ! SEDImentation           ICE3
INTEGER :: NDEPSRS=0    ! DEPosition on Snow      ICE3
INTEGER :: NAGGSRS=0    ! AGGregation of snow     ICE3
INTEGER :: NAUTSRS=0    ! AUToconversion of ice   ICE3
INTEGER :: NRIMRS=0     ! RIMing of cloudwater    ICE3
INTEGER :: NACCRS=0     ! ACCretion of rainwater  ICE3
INTEGER :: NCMELRS=0    ! Conversion MeLTing      ICE3
INTEGER :: NWETGRS=0    ! WET Growth of graupel   ICE3
INTEGER :: NDRYGRS=0    ! DRY Growth of graupel   ICE3
INTEGER :: NWETHRS=0    ! wet growth of hail      ICE4
INTEGER :: NCNVIRS=0   ! Conversion of snow to r_i         C3R5
INTEGER :: NCNVSRS=0   ! Conversion of pristine ice to r_s C3R5
INTEGER :: NHMSRS=0    ! Hallett-Mossop ice multiplication process due to snow riming C3R5
!
!      Allowed processes for the budget of moist variable RRG (graupel)
!
! Courant namelist: NAM_BURRG
!
LOGICAL :: LBU_RRG=.FALSE.    ! True when the budget of RRG is performed
!
INTEGER :: NASSERG=0    ! time filter
INTEGER :: NNESTRG=0    ! Efffect of 2way nesting on Rg
INTEGER :: NADVRG=0    ! Total advection for PPM
INTEGER :: NADVXRG=0    ! advection along X (all except PPM)
INTEGER :: NADVYRG=0    ! advection along Y (all except PPM)
INTEGER :: NADVZRG=0    ! advection along Z (all except PPM)
INTEGER :: NFRCRG=0     ! forcing
INTEGER :: NDIFRG=0     ! numerical diffusion
INTEGER :: NRELRG=0     ! relaxation
INTEGER :: NNEGARG=0    ! negative correction
INTEGER :: NSEDIRG=0    ! SEDImentation           ICE3
INTEGER :: NSFRRG=0     ! Spontaneous FReezing    ICE3
INTEGER :: NDEPGRG=0    ! DEPosition on Snow      ICE3
INTEGER :: NRIMRG=0     ! RIMing of cloudwater    ICE3
INTEGER :: NACCRG=0     ! ACCretion of rainwater  ICE3
INTEGER :: NCMELRG=0    ! Conversion MeLTing      ICE3
INTEGER :: NCFRZRG=0    ! Conversion FReeZing     ICE3
INTEGER :: NWETGRG=0    ! WET Growth of graupel   ICE3
INTEGER :: NDRYGRG=0    ! DRY Growth of graupel   ICE3
INTEGER :: NGMLTRG=0    ! Graupel MeLTing         ICE3
INTEGER :: NWETHRG=0    ! wet growth of hail      ICE4
INTEGER :: NHONRRG=0    ! drop homogeneous nucleation C3R5
INTEGER :: NHMGRG=0     ! Hallett-Mossop ice multiplication process due to graupel riming
INTEGER :: NCOHGRG=0    ! conversion of hail to graupel
!
!      Allowed processes for the budget of moist variable RRH (hail)
!
! Courant namelist: NAM_BURRH
!
LOGICAL :: LBU_RRH=.FALSE.    ! True when the budget of RRH is performed
!
INTEGER :: NASSERH=0    ! time filter
INTEGER :: NNESTRH=0    ! Efffect of 2way nesting on Rh
INTEGER :: NADVRH=0     ! Total advection for PPM
INTEGER :: NADVXRH=0    ! advection along X (all except PPM)
INTEGER :: NADVYRH=0    ! advection along Y (all except PPM)
INTEGER :: NADVZRH=0    ! advection along Z (all except PPM)
INTEGER :: NFRCRH=0     ! forcing
INTEGER :: NDIFRH=0     ! numerical diffusion
INTEGER :: NRELRH=0     ! relaxation
INTEGER :: NNEGARH=0    ! negative correction
INTEGER :: NSEDIRH=0    ! sedimentation
INTEGER :: NWETGRH=0    ! wet growth of graupel
INTEGER :: NWETHRH=0    ! wet growth of hail
INTEGER :: NHMLTRH=0    ! melting
INTEGER :: NCOHGRH=0    ! conversion of hail to graupel
!
! Courant namelist: NAM_BURSV
!
LOGICAL :: LBU_RSV=.FALSE.    ! True when the budget of RSVx is performed
!
INTEGER :: NASSESV=0    ! Asselin-Robert time filter
INTEGER :: NNESTSV=0    ! Efffect of 2way nesting on Sv
INTEGER :: NADVSV=0     ! Total advection for PPM
INTEGER :: NADVXSV=0    ! advection along X (all except PPM)
INTEGER :: NADVYSV=0    ! advection along Y (all except PPM)
INTEGER :: NADVZSV=0    ! advection along Z (all except PPM)
INTEGER :: NFRCSV=0     ! forcing
INTEGER :: NDIFSV=0     ! numerical diffusion
INTEGER :: NRELSV=0     ! relaxation
INTEGER :: NDCONVSV=0   !  Deep CONVection
INTEGER :: NMAFLSV=0    ! mass flux
INTEGER :: NHTURBSV=0   ! horizontal turbulence
INTEGER :: NVTURBSV=0   ! vertical turbulence
INTEGER :: NCHEMSV=0    ! chemistry activity
!
INTEGER :: NNEGASV=0
!
! Allowed processes for the budget of electric charge carried by water vapor
INTEGER :: NDEPSQV=0
INTEGER :: NDEPGQV=0
INTEGER :: NREVAQV=0
INTEGER :: NDEPIQV=0
INTEGER :: NNEUTQV=0
!
! Allowed processes for the budget of electric charge carried by cloud droplets
INTEGER :: NAUTOQC=0
INTEGER :: NACCRQC=0
INTEGER :: NRIMQC=0
INTEGER :: NWETGQC=0
INTEGER :: NDRYGQC=0
INTEGER :: NIMLTQC=0
INTEGER :: NBERFIQC=0
INTEGER :: NDEPIQC=0
INTEGER :: NINDQC=0  ! inductive process
INTEGER :: NSEDIQC=0
INTEGER :: NNEUTQC=0
!
! Allowed processes for the budget of electric charge carried by rain drops
INTEGER :: NAUTOQR=0
INTEGER :: NACCRQR=0
INTEGER :: NREVAQR=0
INTEGER :: NACCQR=0
INTEGER :: NCFRZQR=0
INTEGER :: NWETGQR=0
INTEGER :: NDRYGQR=0
INTEGER :: NGMLTQR=0
INTEGER :: NSEDIQR=0
INTEGER :: NNEUTQR=0
!
! Allowed processes for the budget of electric charge carried by ice crystals
INTEGER :: NAGGSQI=0
INTEGER :: NAUTSQI=0
INTEGER :: NCFRZQI=0
INTEGER :: NWETGQI=0
INTEGER :: NDRYGQI=0
INTEGER :: NIMLTQI=0
INTEGER :: NBERFIQI=0
INTEGER :: NDEPIQI=0
INTEGER :: NNIISQI=0 ! non-inductive I-S
INTEGER :: NSEDIQI=0
INTEGER :: NNEUTQI=0
!
! Allowed processes for the budget of electric charge carried by snow
INTEGER :: NDEPSQS=0
INTEGER :: NAGGSQS=0
INTEGER :: NAUTSQS=0
INTEGER :: NRIMQS=0
INTEGER :: NACCQS=0
INTEGER :: NCMELQS=0
INTEGER :: NWETGQS=0
INTEGER :: NDRYGQS=0
INTEGER :: NNIISQS=0  ! non-inductive I-S
INTEGER :: NSEDIQS=0
INTEGER :: NNEUTQS=0
!
! Allowed processes for the budget of electric charge carried by graupel
INTEGER :: NDEPGQG=0
INTEGER :: NRIMQG=0
INTEGER :: NACCQG=0
INTEGER :: NCMELQG=0
INTEGER :: NCFRZQG=0
INTEGER :: NWETGQG=0
INTEGER :: NDRYGQG=0
INTEGER :: NGMLTQG=0
INTEGER :: NINDQG=0  ! inductive process
INTEGER :: NSEDIQG=0
INTEGER :: NNEUTQG=0
!
! must add processes for electric charge carried by hail
!
!
REAL :: XTIME_BU=0.          ! budget time in this time-step
REAL :: XTIME_BU_PROCESS=0.  ! budget time per process for this time-step
!
LOGICAL, POINTER :: LBUDGET_U=>NULL()    ! flag to compute budget of RhoJu  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_V=>NULL()    ! flag to compute budget of RhoJv  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_W=>NULL()    ! flag to compute budget of RhoJw  and/or LES budgets with u
LOGICAL, POINTER :: LBUDGET_TH=>NULL()   ! flag to compute budget of RhoJTh and/or LES budgets with th
LOGICAL, POINTER :: LBUDGET_TKE=>NULL()  ! flag to compute budget of RhoJTke and/or LES budgets with Tke
LOGICAL, POINTER :: LBUDGET_RV=>NULL()   ! flag to compute budget of RhoJrv and/or LES budgets with rv
LOGICAL, POINTER :: LBUDGET_RC=>NULL()   ! flag to compute budget of RhoJrc and/or LES budgets with rc
LOGICAL, POINTER :: LBUDGET_RR=>NULL()   ! flag to compute budget of RhoJrr and/or LES budgets with rr
LOGICAL, POINTER :: LBUDGET_RI=>NULL()   ! flag to compute budget of RhoJri and/or LES budgets with ri
LOGICAL, POINTER :: LBUDGET_RS=>NULL()   ! flag to compute budget of RhoJrs and/or LES budgets with rs
LOGICAL, POINTER :: LBUDGET_RG=>NULL()   ! flag to compute budget of RhoJrg and/or LES budgets with rg
LOGICAL, POINTER :: LBUDGET_RH=>NULL()   ! flag to compute budget of RhoJrh and/or LES budgets with rh
LOGICAL, POINTER :: LBUDGET_SV=>NULL()   ! flag to compute budget of RhoJsv and/or LES budgets with sv

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
!
END MODULE MODD_BUDGET
