!     ######spl
      SUBROUTINE ARO_SUBUDGET(KLON,KLEV,PTSTEP)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##########################################################################

!**** *ARO_SUBUDGET*   - Initialize common meso_NH MODD_ used in BUDGET for DDH in AROME

!     Purpose.
!     --------
!            Set implicit values for MODD_BUDGET in the way needed for DDH in AROME 

!**   Interface.
!     ----------
!        *CALL* *ARO_SUBUDGET

!        Explicit arguments :
!        --------------------
!       None

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        To use budgets in DDH for AROME, budget must have type CART. 
!        First dimension is NPROMA and second dimension is 1. Budgets 
!        are reset after each tipe step. Processes not used in AROME are
!        assigned 3.
!
!     Externals.
!     ----------

!     Reference.
!     ----------
!        Documentation AROME 

!     Author.
!     -------
!        T. Kovacic 
!
!     Modifications.
!     --------------
!        Original : 05-04-27
!     ------------------------------------------------------------------
!USE MODD_PARAMETERS, ONLY : LWARM
USE MODD_BUDGET
USE MODD_DYN, ONLY : LCORIO

IMPLICIT NONE
!     ------------------------------------------------------------------
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER, INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER, INTENT(IN)   :: KLEV     !Number of vertical levels 
REAL,    INTENT(IN)   :: PTSTEP   ! time step

!*********************************************************************************


REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_SUBUDGET',0,ZHOOK_HANDLE)

CALL TBUCONF_ASSOCIATE()

LBU_ENABLE = .TRUE.

LCORIO = .TRUE.

CBUTYPE = 'CART'

NBUMOD= 1
!
XBULEN = PTSTEP
!
NBUSTEP = 1
                                           
XBUWRI = -1  ! not used

NBUWRNB = -1 ! not used, in INI_BUDGET
                                           
NBUTSHIFT = -1 ! not used
!
NBUKL = 2
NBUKH = KLEV + 1
LBU_KCP = .FALSE.

!
!  Variables used by the cartesian box case ('CART') only
!
NBUIL = 1
NBUIH = KLON
NBUJL = 1 
NBUJH = 1
LBU_ICP = .FALSE.
LBU_JCP = .FALSE.
!
!                Variables used by the  mask case ('MASK') only
!
NBUMASK = -1 ! not used, in INI_BUDGET
!             
NBUTIME = -1 ! not used
!

!
!*******************************************************************************
!****    W A R N I N G ! ! !
! Budgets are initialized with subroutine AROINI_BUDGET which calls MNH 
! subroutine INI_BUDGET. In this subroutine budgets are defined so that 4th 
! process must be Asselin filter and it can be activated only with 0 or 1. 
! Our choice is 0. With this choise we don't have to make a budget for Asselin 
! filter in ARO_STARTBU as it is obligatory in INITIAL_GUESS, it's equivalent in MNH. 
!*********************************************************************************


!  Budget array XBURU      : RU budget
! MNH namelist: NAM_BURU

  LBU_RU=.TRUE.
  NASSEU  =0  ! time filter
  NNESTU  =3  ! Efffect of 2way nesting on U
  NADVXU  =3  ! advection along X
  NADVYU  =3  ! advection along Y
  NADVZU  =3  ! advection along Z 
  NFRCU   =3  ! forcing
  NCURVU  =3  ! curvature
  NCORU   =3  ! Coriolis terms 
  NDIFU   =3  ! numerical diffusion
  NRELU   =3  ! relaxation
  NHTURBU =4  !horizontal TURBulence
  NVTURBU =1  !vertical turbulence 
  NPRESU  =3  ! pressure term

!  Budget array XBURV, 

  LBU_RV=.TRUE. !  when the budget of RV is  not performed
!
  NASSEV   =0  ! time filter
  NNESTV   =3  ! Efffect of 2way nesting on V
  NADVXV   =3  ! advection along X 
  NADVYV   =3  ! advection along Y 
  NADVZV   =3  ! advection along Z 
  NFRCV    =3  ! forcing
  NCURVV   =3  ! curvature
  NCORV    =3  ! Coriolis terms 
  NDIFV    =3  ! numerical diffusion
  NRELV    =3  ! relaxation
  NHTURBV  =4  ! horizontal turbulence
  NVTURBV  =1  ! vertical turbulence 
  NPRESV   =3  ! pressure term

!  Budget array XBURW 
!
!      Allowed processes for the budget of RW (wind vertical component)
!                                                  
! MNH namelist: NAM_BURW
!
  LBU_RW=.TRUE. ! when the budget of RW is not performed 
!                                                  
  NASSEW   = 0  ! time filter
  NNESTW   = 3 ! Efffect of 2way nesting on W
  NADVXW   = 3 ! advection along X
  NADVYW   = 3 ! advection along Y
  NADVZW   = 3 ! advection along Z 
  NFRCW    = 3 ! forcing
  NCURVW   = 3 ! curvature
  NCORW    = 3 ! Coriolis terms 
  NGRAVW   = 3 ! gravity term
  NDIFW    = 3 ! numerical diffusion
  NRELW    = 3 ! relaxation
  NHTURBW  = 4 ! horizontal turbulence 
  NVTURBW  = 1 ! vertical turbulence 
  NPRESW   = 3 ! pressure term
!
!      Allowed processes for the budget of RTH (potential temperature)
!                                                  
!  Budget array XBURTH 
! MNH namelist: NAM_BURTH
!
  LBU_RTH=.TRUE. ! when the budget of RTH is not performed
!
  NASSETH  = 0  ! time filter
  NNESTTH  = 3 ! Efffect of 2way nesting on Th
  NADVXTH  = 3 ! advection along X
  NADVYTH  = 3 ! advection along Y
  NADVZTH  = 3 ! advection along Z 
  NFRCTH   = 3 ! forcing
  NPREFTH  = 3 ! theta source term due to the reference pressure
                            ! (Dyn. Sources) only present if KRR>0
  NDIFTH   = 3 ! numerical diffusion
  NRELTH   = 3 ! relaxation
  NRADTH   = 1 ! RADiation
  NDCONVTH = 1 ! Deep CONVection
  NHTURBTH = 4 ! horizontal turbulence
  NVTURBTH = 1 ! vertical turbulence
  NDISSHTH = 1 ! dissipative heating
  NNEGATH  = 1 ! negative correction induced by hydrometeors
  NREVATH  = 1 ! rain evaporation
  NCONDTH  = 3 ! evaporation/condensation
  NHENUTH  = 1 ! HEterogenous NUcleation ICE3
  NHONTH   = 1 ! HOmogeneous Nucleation  ICE3
  NSFRTH   = 1 ! Spontaneous FReezing    ICE3
  NDEPSTH  = 1 ! DEPosition on Snow      ICE3
  NDEPGTH  = 1 ! DEPosition on Graupel   ICE3
  NRIMTH   = 1 ! RIMing of cloudwater    ICE3
  NACCTH   = 1 ! ACCretion of rainwater  ICE3
  NCFRZTH  = 1 ! Conversion FReeZing     ICE3
  NWETGTH  = 1 ! WET Growth of graupel   ICE3
  NDRYGTH  = 1 ! DRY Growth of graupel   ICE3
  NGMLTTH  = 1 ! Graupel MeLTing         ICE3
  NIMLTTH  = 1 ! Ice MeLTing             ICE3
  NBERFITH = 1 ! BERgeron-FIndeisen gth. ICE3
  NCDEPITH = 1 ! Cond./DEPosition on ice ICE3
!
!      Allowed processes for the budget of RTKE (kinetic energy)
!                                                  
!  Budget array XBURTKE     
! MNH namelist: NAM_BURTKE
!
  LBU_RTKE=.TRUE. ! when the budget of RTKE is not performed
!
  NASSETKE = 0  ! time filter
  NADVXTKE = 3 ! advection along X
  NADVYTKE = 3 ! advection along Y
  NADVZTKE = 3 ! advection along Z 
  NFRCTKE  = 3 ! forcing
  NDIFTKE  = 3 ! numerical diffusion
  NRELTKE  = 3 ! relaxation
  NDPTKE   = 1 ! dynamic production of TKE
  NTPTKE   = 1 ! thermal production of TKE
  NDISSTKE = 1 ! dissipation of TKE
  NTRTKE   = 1 ! turbulent transport of TKE
!
!
!      Allowed processes for the budget of moist variable RRV (water vapor)
!                                                  
!  Budget array XBURRV 
! MNH namelist: NAM_BURRV
!
  LBU_RRV=.TRUE. ! when the budget of RRV is not performed
!
  NASSERV  = 0  ! time filter
  NNESTRV  = 3 ! Efffect of 2way nesting on Rv
  NADVXRV  = 3 ! advection along X
  NADVYRV  = 3 ! advection along Y
  NADVZRV  = 3 ! advection along Z 
  NFRCRV   = 3 ! forcing
  NDIFRV   = 3 ! numerical diffusion
  NRELRV   = 3 ! relaxation
  NDCONVRV = 1 ! Deep CONVection
  NHTURBRV = 4 ! horizontal turbulence 
  NVTURBRV = 1 ! vertical turbulence
  NNEGARV  = 1 ! negative correction                            
  NREVARV  = 1 ! rain evaporation
  NCONDRV  = 3 ! evaporation/condensation
  NHENURV  = 1 ! HEterogenous NUcleation ICE3
  NDEPSRV  = 1 ! DEPosition on Snow      ICE3
  NDEPGRV  = 1 ! DEPosition on Graupel   ICE3
  NCDEPIRV = 1 ! Cond./DEPosition on ice ICE3
!
!  Budget array XBURRC      
!      Allowed processes for the budget of moist variable RRC (cloud water)
!                                                  
! MNH namelist: NAM_BURRC
!
  LBU_RRC=.TRUE. ! when the budget of RRC is not performed
!
  NASSERC  = 0  ! time filter
  NNESTRC  = 3 ! Efffect of 2way nesting on Rc
  NADVXRC  = 3 ! advection along X
  NADVYRC  = 3 ! advection along Y
  NADVZRC  = 3 ! advection along Z 
  NFRCRC   = 3 ! forcing
  NDIFRC   = 3 ! numerical diffusion
  NRELRC   = 3 ! relaxation
  NDCONVRC = 1 ! Deep CONVection
  NHTURBRC = 4 ! horizontal turbulence 
  NVTURBRC = 1 ! vertical turbulence
  NNEGARC  = 1 ! negative correction                            
  NACCRRC  = 1 ! accretion
  NAUTORC  = 1 ! autoconversion
  NCONDRC  = 3 ! evaporation/condensation
  NHONRC   = 1 ! HOmogeneous Nucleation  ICE3
  NRIMRC   = 1 ! RIMing of cloudwater    ICE3
  NWETGRC  = 1 ! WET Growth of graupel   ICE3
  NDRYGRC  = 1 ! DRY Growth of graupel   ICE3
  NIMLTRC  = 1 ! Ice MeLTing             ICE3
  NBERFIRC = 1 ! BERgeron-FIndeisen gth. ICE3
  NCDEPIRC = 1 ! Cond./DEPosition on ice ICE3
  NHENURC  = 3 ! CCN Activation C2R2
  NSEDIRC  = 3 ! sedimentation  C2R2
!
!      Allowed processes for the budget of moist variable RRR (rain water)
!
!  Budget array XBURRR 
! MNH namelist: NAM_BURRR
!
  LBU_RRR=.TRUE. ! when the budget of RRR is not performed
!
  NASSERR  = 0  ! time filter
  NNESTRR  = 3 ! Efffect of 2way nesting on Rr
  NADVXRR  = 3 ! advection along X
  NADVYRR  = 3 ! advection along Y
  NADVZRR  = 3 ! advection along Z 
  NFRCRR   = 3 ! forcing
  NDIFRR   = 3 ! numerical diffusion
  NRELRR   = 3 ! relaxation
  NNEGARR  = 1 ! negative correction                            
  NACCRRR  = 1 ! accretion
  NAUTORR  = 1 ! autoconversion
  NREVARR  = 1 ! rain evaporation
  NSEDIRR  = 1 ! sedimentation
  NSFRRR   = 1 ! Spontaneous FReezing    ICE3
  NACCRR   = 1 ! ACCretion of rainwater  ICE3
  NCFRZRR  = 1 ! Conversion FReeZing     ICE3
  NWETGRR  = 1 ! WET Growth of graupel   ICE3
  NDRYGRR  = 1 ! DRY Growth of graupel   ICE3
  NGMLTRR  = 1 ! Graupel MeLTing         ICE3
!
!      Allowed processes for the budget of moist variable RRI (ice)
!
!  Budget array XBURRI      
! MNH namelist: NAM_BURRI
!
  LBU_RRI=.TRUE. ! when the budget of RRI is not performed
!
  NASSERI  = 0  ! time filter
  NNESTRI  = 3 ! Efffect of 2way nesting on Ri
  NADVXRI  = 3 ! advection along X
  NADVYRI  = 3 ! advection along Y
  NADVZRI  = 3 ! advection along Z 
  NFRCRI   = 3 ! forcing
  NDIFRI   = 3 ! numerical diffusion
  NRELRI   = 3 ! relaxation
  NDCONVRI = 1 ! Deep CONVection
  NHTURBRI = 4 ! horizontal turbulence
  NVTURBRI = 1 ! vertical turbulence
  NNEGARI  = 1 ! negative correction                            
  NSEDIRI  = 1 ! SEDImentation           ICE3
  NHENURI  = 1 ! HEterogenous NUcleation ICE3
  NHONRI   = 1 ! HOmogeneous Nucleation  ICE3
  NAGGSRI  = 1 ! AGGregation of snow     ICE3
  NAUTSRI  = 1 ! AUToconversion of ice   ICE3
  NCFRZRI  = 1 ! Conversion FReeZing     ICE3
  NWETGRI  = 1 ! WET Growth of graupel   ICE3
  NDRYGRI  = 1 ! DRY Growth of graupel   ICE3
  NIMLTRI  = 1 ! Ice MeLTing             ICE3
  NBERFIRI = 1 ! BERgeron-FIndeisen gth. ICE3
  NCDEPIRI = 1 ! Cond./DEPosition on ice ICE3
!
!      Allowed processes for the budget of moist variable RRS (snow)
!
!  Budget array XBURRS 
! MNH namelist: NAM_BURRS
!
  LBU_RRS=.TRUE. ! when the budget of RRS is not performed
!
  NASSERS  = 0  ! time filter
  NNESTRS  = 3 ! Efffect of 2way nesting on Rs
  NADVXRS  = 3 ! advection along X
  NADVYRS  = 3 ! advection along Y
  NADVZRS  = 3 ! advection along Z 
  NFRCRS   = 3 ! forcing
  NDIFRS   = 3 ! numerical diffusion
  NRELRS   = 3 ! relaxation
  NNEGARS  = 1 ! negative correction                            
  NSEDIRS  = 1 ! SEDImentation           ICE3
  NDEPSRS  = 1 ! DEPosition on Snow      ICE3
  NAGGSRS  = 1 ! AGGregation of snow     ICE3
  NAUTSRS  = 1 ! AUToconversion of ice   ICE3
  NRIMRS   = 1 ! RIMing of cloudwater    ICE3
  NACCRS   = 1 ! ACCretion of rainwater  ICE3
  NCMELRS  = 1 ! Conversion MeLTing      ICE3
  NWETGRS  = 1 ! WET Growth of graupel   ICE3
  NDRYGRS  = 1 ! DRY Growth of graupel   ICE3
!
!      Allowed processes for the budget of moist variable RRG (graupel)
!
!  Budget array XBURRG     
! MNH namelist: NAM_BURRG
!
  LBU_RRG=.TRUE. ! when the budget of RRG is not performed
!
  NASSERG  = 0  ! time filter
  NNESTRG  = 3 ! Efffect of 2way nesting on Rg
  NADVXRG  = 3 ! advection along X
  NADVYRG  = 3 ! advection along Y
  NADVZRG  = 3 ! advection along Z 
  NFRCRG   = 3 ! forcing
  NDIFRG   = 3 ! numerical diffusion
  NRELRG   = 3 ! relaxation
  NNEGARG  = 1 ! negative correction                            
  NSEDIRG  = 1 ! SEDImentation           ICE3
  NSFRRG   = 1 ! Spontaneous FReezing    ICE3
  NDEPGRG  = 1 ! DEPosition on Snow      ICE3
  NRIMRG   = 1 ! RIMing of cloudwater    ICE3
  NACCRG   = 1 ! ACCretion of rainwater  ICE3
  NCMELRG  = 1 ! Conversion MeLTing      ICE3
  NCFRZRG  = 1 ! Conversion FReeZing     ICE3
  NWETGRG  = 1 ! WET Growth of graupel   ICE3
  NDRYGRG  = 1 ! DRY Growth of graupel   ICE3
  NGMLTRG  = 1 ! Graupel MeLTing         ICE3
!
!      Allowed processes for the budget of moist variable RRH (hail)
!
!  Budget array XBURRH              ! 
! MNH namelist: NAM_BURRH
!
  LBU_RRH=.FALSE. ! when the budget of RRH is not performed
!
  NASSERH  = 0  ! time filter
  NNESTRH  = 3 ! Efffect of 2way nesting on Rh
  NADVXRH  = 3 ! advection along X
  NADVYRH  = 3 ! advection along Y
  NADVZRH  = 3 ! advection along Z 
  NFRCRH   = 3 ! forcing
  NDIFRH   = 3 ! numerical diffusion
  NRELRH   = 3 ! relaxation
  NNEGARH  = 3 ! negative correction                            
!
! MNH namelist: NAM_BURSV
!
  LBU_RSV=.FALSE. ! the budget of RSVx is not performed
!
  NASSESV  = 0 ! Asselin-Robert time filter
  NNESTSV  = 3 ! Efffect of 2way nesting on Sv
  NADVXSV  = 3 ! advection along X
  NADVYSV  = 3 ! advection along Y
  NADVZSV  = 3 ! advection along Z
  NFRCSV   = 3 ! forcing
  NDIFSV   = 3 ! numerical diffusion
  NRELSV   = 3 ! relaxation
  NDCONVSV = 1 !  Deep CONVection
  NHTURBSV = 4 ! horizontal turbulence
  NVTURBSV = 1 ! vertical turbulence
  NCHEMSV  = 1 ! chemistry activity
!
!
!
LBUDGET_U  =.TRUE.  ! flag to compute budget of RhoJu  and/or LES budgets with u
LBUDGET_V  =.TRUE.  ! flag to compute budget of RhoJv  and/or LES budgets with u
LBUDGET_W  =.TRUE.  ! flag to compute budget of RhoJw  and/or LES budgets with u
LBUDGET_TH =.TRUE.  ! flag to compute budget of RhoJTh and/or LES budgets with th
LBUDGET_TKE=.TRUE.  ! flag to compute budget of RhoJTke and/or LES budgets with Tke
LBUDGET_RV =.TRUE.  ! flag to compute budget of RhoJrv and/or LES budgets with rv
LBUDGET_RC =.TRUE.  ! flag to compute budget of RhoJrc and/or LES budgets with rc
LBUDGET_RR =.TRUE.  ! flag to compute budget of RhoJrr and/or LES budgets with rr
LBUDGET_RI =.TRUE.  ! flag to compute budget of RhoJri and/or LES budgets with ri
LBUDGET_RS =.TRUE.  ! flag to compute budget of RhoJrs and/or LES budgets with rs
LBUDGET_RG =.TRUE.  ! flag to compute budget of RhoJrg and/or LES budgets with rg
LBUDGET_RH =.TRUE. ! flag to compute budget of RhoJrh and/or LES budgets with rh
LBUDGET_SV =.FALSE. ! flag to compute budget of RhoJsv and/or LES budgets with sv
!
IF (LHOOK) CALL DR_HOOK('ARO_SUBUDGET',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE ARO_SUBUDGET
               
