MODULE YEMLBC_MODEL

! Purpose :
! -------
!    Forcing a LAM model by another model: part 0B
!    - forcing by lateral boundary conditions
!    - pressure tendency coupling
!    - spectral nudging

! Interface :
! ---------
!    Empty.

! External :
! --------
!    None.

! Method :
! ------
!    See Documentation.

! Reference :
! ---------

! Author :
! ------
!    K. YESSAD (CNRM/GMAP) after YEMBICU, YEMDYN, YEMGT3B, SUEBICU, SUEDYN, SUESC2.
! Original : December 2010

! Modifications :
! -------------
! Daan Degrauwe: Feb 2012 Boyd biperiodization         
! T. Wilhelmsson and K. Yessad (Oct 2013) Geometry and setup refactoring.
! B. Bochenek (Oct 2013): Weights for LBC interpolation
! K. Yessad (July 2014): Move some variables.
! F. Taillefer (Aug 2015)  Add control of no coupling at all in canari by namelist
! M. Hortal (Dec 2014): Upper boundary relaxation
! P. Marguinaud (Oct 2016) : Port to single precision
! J. Vivoda (Mar 2017): Fixing bug in options LQCPL and LCCPL
! WARNING! The bug is in swapping LBC buffers P*CPL entering ESC2R. Fix does
! not remove it, but adjusts interpolation weights EWB accordingly. Clean
! solution is to correct swapping of LBC buffers, otherwise the code will
! remain difficult to understand.
! H. Dhouioui (Sep 2017) renamed from elbc0b_mod.F90
! O. Vignes (Feb 2020): Upper boundary relaxation fixes
! M. Hamrud (Oct 2021) incorporated YEMLBC_INIT. Forced by move of YOMDYNA 
! into model object
!-----------------------------------------------------------------------------

USE PARKIND1 , ONLY : JPIM, JPRB, JPRD
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMHOOK  , ONLY : LHOOK, DR_HOOK, JPHOOK
USE YOMDFI   , ONLY : NSTDFI, NSTDFIA, RTDFI, RTDFIA
USE YOMCT0   , ONLY : LRPLANE, LALLOPR, LCANARI,NCONF,LELAM
USE YOMGMV   , ONLY : TGMV
USE YOMINI   , ONLY : LDFI
USE YOMLUN   , ONLY : NULOUT, NULNAM
USE YOMMP0   , ONLY : MYSETB, MYSETV, LOUTPUT

IMPLICIT NONE
SAVE

!=============================================================================

!      1.    TYPE DEFINITION
!            ---------------

! Moved form YEMLBC_FIELDS
! Structure for GMVS coupled fields in LTENC option.
TYPE TGMVSTENC
INTEGER(KIND=JPIM) :: MSP        ! surface pressure variable
INTEGER(KIND=JPIM) :: MSPL       ! zonal component of grad(surface pressure variable)
INTEGER(KIND=JPIM) :: MSPM       ! meridian component of grad(surface pressure variable)
INTEGER(KIND=JPIM) :: NDIM       ! number of coupled fields (includes derivatives)
INTEGER(KIND=JPIM) :: NDIMT      ! number of temporally interpolated fields (includes derivatives)
END TYPE TGMVSTENC

TYPE TGMVCPL
! coupled GMV
INTEGER(KIND=JPIM) :: MU         ! U-wind
INTEGER(KIND=JPIM) :: MV         ! V-wind
INTEGER(KIND=JPIM) :: MT         ! temperature
INTEGER(KIND=JPIM) :: MSPD       ! pressure departure variable
INTEGER(KIND=JPIM) :: MSVD       ! vertical divergence variable
INTEGER(KIND=JPIM) :: MNHX       ! NHX term
! derivatives required for linear terms calculation (in ESEIMPLS)
INTEGER(KIND=JPIM) :: MDIV       ! horizontal divergence
INTEGER(KIND=JPIM) :: MTL        ! zonal component of grad(temperature)
INTEGER(KIND=JPIM) :: MTM        ! meridian component of grad(temperature)
INTEGER(KIND=JPIM) :: MSPDL      ! zonal component of grad(pressure departure variable)
INTEGER(KIND=JPIM) :: MSPDM      ! meridian component of grad(pressure departure variable)
INTEGER(KIND=JPIM) :: NDIM       ! number of coupled fields (does not include derivatives)
INTEGER(KIND=JPIM) :: NDIMT      ! number of temporally interpolated fields (includes derivatives)
END TYPE TGMVCPL

TYPE TGMVSCPL
INTEGER(KIND=JPIM) :: MSP        ! surface pressure variable
INTEGER(KIND=JPIM) :: MSPL       ! zonal component of grad(surface pressure variable)
INTEGER(KIND=JPIM) :: MSPM       ! meridian component of grad(surface pressure variable)
INTEGER(KIND=JPIM) :: NDIM       ! number of coupled fields (does not include derivatives)
INTEGER(KIND=JPIM) :: NDIMT      ! number of temporally interpolated fields (includes derivatives)
END TYPE TGMVSCPL

! ** IMPORTED FROM YEMLBC_INIT
INTEGER(KIND=JPIM), PARAMETER :: JPLSGT=20
INTEGER(KIND=JPIM), PARAMETER :: JPALFNM=31

TYPE :: TELBC_MODEL

!      1.1   Coupling of surface pressure tendency

! LTENC   : TRUE if tendency coupling of surface pressure is switched on
! LALLTC  : used together with LTENC when LTENC=.T.
!           - no meaning for quadratic tendency coupling, where just t1 coupling
!             is applied at every NEFRCL time step
!           - for lin. tendency coupling:
!             TRUE if tendency coupling of surf. pres. at every step
!             FALSE if tend. coupl., except at every NEFRCL time steps
!             when just t1 coupling

LOGICAL :: LTENC
LOGICAL :: LALLTC
LOGICAL :: LRFIRST   ! Force reading of first coupling file (usually, it is the
                     ! same as the initial conditions file)

!      1.2   Lateral forcing

!   NBICOU  : controls coupling of wind components (GMV)
!   NBICOT  : controls coupling of temperature (GMV)
!   NBICPD  : controls coupling of pressure departure variable (GMV)
!   NBICVD  : controls coupling of vertical divergence variable (GMV)
!   NBICNHX : controls coupling of "NHX" term (GMV)
!   NBICOP  : controls coupling of surface pressure (GMVS)
!   Possible value for the NBIC[X] variables:
!   * 0: no coupling
!   * 1: default coupling
!   * 2: specific coupling function

!   NECRIPL : controls timelevel of coupling
!   * 0: coupling at t
!   * 1: coupling at t+dt
!   * 2: coupling at t and t+dt
!   LQCPL   : if T (resp. F), quadratic (resp. linear) temporal interpolation
!   LCCPL   : if T (resp. F), cubic (resp. linear) temporal interpolation
!   NECOTL  : Controls the coupling in the tangent linear model:
!             0   - no coupling
!             < 0 - coupling with the Davies relaxation
!             -1  - the linear time interpolation is switched on
!             > 0 - coupling other than Davies relaxation (to be implemented)
!   NECOAD  : cf. NECOTL but for the adjoint model.
!   LE0COTA : TRUE if the relaxation in the I+E zone is towards
!             nullified boundary conditions (assumed exact);
!             FALSE if the relaxation is towards a predefined forcing.
!             LE0COTA has the same meaning in the TL and AD models
!             (respectively, relaxation to 0 perturbation or sensitivity)
!   LEREADINI: TRUE if initial historical file has to be read
!              (used for E501, E801 - not a namelist parameter)
!   N1LSG   : if =1, the gradient with respect to the large scale
!             coupling data has to be written into a file.
!   NFRLSG  : frequency of writing large scale gradients (time or steps)
!   NLSGTS  : array containing large scale gradients timesteps
!             * if NLSGTS(0)=0 action if MOD(JSTEP,NFRLSG)=0
!             * if NLSGTS(0)>0 NLSGTS(0) significant numbers in NLSGTS
!               are then considered and action for JSTEP=NLSGTS(.)*NFRLSG
!             * if NLSGTS(0)<0 action for JSTEP=(NLSGTS(.)/DELTAT)*NFRLSG
!   LRDLSG  : switch for using boundary data perturbation in conf 1
!             for sensitivity forecast run.
!   JPALFNM : Dimension for reading alpha function parameters

INTEGER(KIND=JPIM) :: NBICOU
INTEGER(KIND=JPIM) :: NBICOT
INTEGER(KIND=JPIM) :: NBICPD
INTEGER(KIND=JPIM) :: NBICVD
INTEGER(KIND=JPIM) :: NBICNHX
INTEGER(KIND=JPIM) :: NBICOP
INTEGER(KIND=JPIM) :: NECRIPL
LOGICAL :: LQCPL 
LOGICAL :: LCCPL
INTEGER(KIND=JPIM) :: NECOTL
INTEGER(KIND=JPIM) :: NECOAD
LOGICAL :: LE0COTA
LOGICAL :: LEREADINI
INTEGER(KIND=JPIM) :: N1LSG
INTEGER(KIND=JPIM) :: NFRLSG
INTEGER(KIND=JPIM) :: NLSGTS(0:JPLSGT)
LOGICAL :: LRDLSG

!      1.3   Spectral nudging

!   LESPCPL  : control of spectral nudging

LOGICAL :: LESPCPL
LOGICAL :: LSPTENC

!      1.4   Upper nesting boundary conditions

!   LUNBC    : controls upper nesting boundary conditions

LOGICAL :: LUNBC


!      2.1   LECOBI.

! LECOBI  : T if there is coupling and biperiodicisation

LOGICAL :: LECOBI

! ** END OF IMPORTED FROM YEMLBC_INIT

! Moved from YEMLBC_FIELDS
!      1.1   Number of coupled fields, structures for coupled fields.

! YYTGMVSTENC  : contains pointers and number of coupled fields for GMVS in LTENC option
! YYTGMVCPL    : contains pointers and number of coupled fields for GMV
! YYTGMVSCPL   : contains pointers and number of coupled fields for GMVS
! NDIMCPL      : number of GFL fields with true LCOUPLING attribute.
! NGALEF       : total number of coupled fields.

TYPE(TGMVSTENC) :: YYTGMVSTENC
TYPE(TGMVCPL) :: YYTGMVCPL
TYPE(TGMVSCPL) :: YYTGMVSCPL
INTEGER(KIND=JPIM) :: NDIMCPL
INTEGER(KIND=JPIM) :: NGALEF
! End ofMoved from YEMLBC_FIELDS

!      2.    DECLARATIONS
!            ------------

!      2.0   Control frequency of LBC.

! NEFRCL   : frequency of updating the lateral boundary coupling fields.
!            The LBC fields will be updated every NEFRCL time steps.
! NETLS1   : Time step of the first set of lateral boundary  fields.
! TEFRCL   : time interval between two updatings of the lateral boundary fields

INTEGER(KIND=JPIM) :: NEFRCL
INTEGER(KIND=JPIM) :: NETLS1
REAL(KIND=JPRB) :: TEFRCL

!      2.2   Namelist variables for relaxation coefficients (resp for GMV, GMVS, GFL).

REAL(KIND=JPRB) :: EPA_GMV(JPALFNM)
REAL(KIND=JPRB) :: EPA_GMVS(JPALFNM)
REAL(KIND=JPRB) :: EPA_GFL(JPALFNM)

!      2.3   Relaxation coefficients.

! EALFA_GMV    : relaxation coefficients alpha for GMV.
! EALFA_GMVS   : relaxation coefficients alpha for GMVS.
! EALFA_GFL    : relaxation coefficients alpha for GFL.
! EALFA_TENC   : relaxation coefficients alpha for LTENC (GMVS only).
! EALFAGT3GMV  : ALFA (relax. coef.) of coupling points for GMV
! EALFAGT3GMVS : ALFA (relax. coef.) of coupling points for GMVS
! EALFAGT3GFL  : ALFA (relax. coef.) of coupling points for GFL
! EALFAU_GMV   : relaxation coefficients alpha for GMV (upper boundary).
! EALFAU_GFL   : relaxation coefficients alpha for GFL (upper boundary).

REAL(KIND=JPRB),ALLOCATABLE:: EALFA_GMV(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EALFA_GMVS(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EALFA_GFL(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EALFA_TENC(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EALFAGT3GMV(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EALFAGT3GMVS(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EALFAGT3GFL(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EALFAU_GMV(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EALFAU_GFL(:,:)

!      2.4   Other variables for grid-point coupling.

! GMGT3        : GM array of coupling points
! GMGT4        : GM array of coupling points (upper boundary).
! EWB          : weights for couplings
! EWBDFIFW     : weights for forward DFI
! EWBDFIBW     : weights for backward DFI
! RTENC        : multiplier of EALFA in the tendency coupling scheme
!                for stability reasons (RTENC<=1. close to 1)

REAL(KIND=JPRB),ALLOCATABLE:: GMGT3(:,:)
REAL(KIND=JPRB),ALLOCATABLE:: GMGT4(:)
REAL(KIND=JPRB),ALLOCATABLE:: EWB(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EWBDFIFW(:,:,:,:)
REAL(KIND=JPRB),ALLOCATABLE:: EWBDFIBW(:,:,:,:)
REAL(KIND=JPRB) :: RTENC

!      2.5   Other variables for spectral nudging.

! LSPNUSPDL  : .TRUE. if spectral nudging on Ps is relevant on this MPI task
! RNUDTFRAC  : Time fraction for spectral nudging
! NEFRSPCPL  : frequency of spectral nudging
! NEK0,NEK1  : lower and upper limits for total wavenumber for spectral nudging
! NEN1,NEN2  : lower and upper model levels for spectral nudging
! SPNUDVOR   : spectral nudging coeficient for vorticity
! SPNUDDIV   : spectral nudging coeficient for divergence
! SPNUDT     : spectral nudging coeficient for temperature
! SPNUDQ     : spectral nudging coeficient for specific humidity
! SPNUDSP    : spectral nudging coeficient for surface pressure
! LNUDSPGFL  : An array to control if any spectral GFL, for nudging

LOGICAL, ALLOCATABLE :: LNUDSPGFL(:)
LOGICAL :: LSPNUSPDL
REAL(KIND=JPRB) :: RNUDTFRAC
INTEGER(KIND=JPIM) :: NEFRSPCPL
INTEGER(KIND=JPIM) :: NEK0
INTEGER(KIND=JPIM) :: NEK1
INTEGER(KIND=JPIM) :: NEN1
INTEGER(KIND=JPIM) :: NEN2
REAL(KIND=JPRB) :: SPNUDVOR
REAL(KIND=JPRB) :: SPNUDDIV
REAL(KIND=JPRB) :: SPNUDT
REAL(KIND=JPRB) :: SPNUDQ
REAL(KIND=JPRB) :: SPNUDSP
REAL(KIND=JPRB) :: RNUTENC

END TYPE TELBC_MODEL


!=============================================================================

CONTAINS

SUBROUTINE SUELBC_FIELDS_DIM(YDML_LBC,YDGEOMETRY,YDDYNA,YGFL,KNFD2D,KNS3D)
!--------------------------------------------------------------------------
! Sets-up part 0C of forcing a LAM model by another model
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------

USE GEOMETRY_MOD , ONLY : GEOMETRY
USE YOM_YGFL     , ONLY : TYPE_GFLD
USE YOMDYNA      , ONLY : TDYNA
TYPE(TELBC_MODEL)     , INTENT(INOUT) :: YDML_LBC 
TYPE(GEOMETRY)    , INTENT(IN)    :: YDGEOMETRY
TYPE(TDYNA)       , INTENT(IN)    :: YDDYNA
TYPE(TYPE_GFLD)   , INTENT(INOUT) :: YGFL
INTEGER(KIND=JPIM), INTENT(IN)    :: KNFD2D
INTEGER(KIND=JPIM), INTENT(IN)    :: KNS3D
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER(KIND=JPIM) :: IWS_TENC, ICPL, JGFL, IW, IW_TENC

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YEMLBC_FIELDS:SUELBC_FIELDS_DIM',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------

!      Part A: calculations and allocations.

! * Calculation of YYGMVSTENC, YYTGMVCPL, YYGMVSCPL, NDIMCPL, NGALEF:


! YYGMVSTENC:
IF (YDML_LBC%LTENC) THEN
  YDML_LBC%YYTGMVSTENC%MSP=1
  YDML_LBC%YYTGMVSTENC%MSPL=2
  YDML_LBC%YYTGMVSTENC%MSPM=3
  YDML_LBC%YYTGMVSTENC%NDIM=3
  YDML_LBC%YYTGMVSTENC%NDIMT=3
ELSE
  YDML_LBC%YYTGMVSTENC%MSP=1
  YDML_LBC%YYTGMVSTENC%MSPL=1
  YDML_LBC%YYTGMVSTENC%MSPM=1
  YDML_LBC%YYTGMVSTENC%NDIM=1
  YDML_LBC%YYTGMVSTENC%NDIMT=1
ENDIF

! YYTGMVCPL:
IF (NCONF == 701) THEN
  ! derivatives are useless because ESEIMPLS is not called.
  IF (YDDYNA%LNHDYN.AND.YDDYNA%LNHX) THEN
    YDML_LBC%YYTGMVCPL%MU   = 1
    YDML_LBC%YYTGMVCPL%MV   = 2
    YDML_LBC%YYTGMVCPL%MT   = 3
    YDML_LBC%YYTGMVCPL%MSPD = 4
    YDML_LBC%YYTGMVCPL%MSVD = 5
    YDML_LBC%YYTGMVCPL%MNHX = 6
    YDML_LBC%YYTGMVCPL%NDIM = 6
    YDML_LBC%YYTGMVCPL%NDIMT= 6
  ELSEIF (YDDYNA%LNHDYN.AND.(.NOT.YDDYNA%LNHX)) THEN
    YDML_LBC%YYTGMVCPL%MU   = 1
    YDML_LBC%YYTGMVCPL%MV   = 2
    YDML_LBC%YYTGMVCPL%MT   = 3
    YDML_LBC%YYTGMVCPL%MSPD = 4
    YDML_LBC%YYTGMVCPL%MSVD = 5
    YDML_LBC%YYTGMVCPL%MNHX = 5
    YDML_LBC%YYTGMVCPL%NDIM = 5
    YDML_LBC%YYTGMVCPL%NDIMT= 5
  ELSE
    YDML_LBC%YYTGMVCPL%MU   = 1
    YDML_LBC%YYTGMVCPL%MV   = 2
    YDML_LBC%YYTGMVCPL%MT   = 3
    YDML_LBC%YYTGMVCPL%MSPD = 3
    YDML_LBC%YYTGMVCPL%MSVD = 3
    YDML_LBC%YYTGMVCPL%MNHX = 3
    YDML_LBC%YYTGMVCPL%NDIM = 3
    YDML_LBC%YYTGMVCPL%NDIMT= 3
  ENDIF
ELSE
! derivatives are useful because ESEIMPLS is called.
  IF (YDDYNA%LNHDYN.AND.YDDYNA%LNHX) THEN
    YDML_LBC%YYTGMVCPL%MU   = 1
    YDML_LBC%YYTGMVCPL%MV   = 2
    YDML_LBC%YYTGMVCPL%MT   = 3
    YDML_LBC%YYTGMVCPL%MSPD = 4
    YDML_LBC%YYTGMVCPL%MSVD = 5
    YDML_LBC%YYTGMVCPL%MNHX = 6
    YDML_LBC%YYTGMVCPL%NDIM = 6
    YDML_LBC%YYTGMVCPL%MDIV = 7
    YDML_LBC%YYTGMVCPL%MTL  = 8
    YDML_LBC%YYTGMVCPL%MTM  = 9
    YDML_LBC%YYTGMVCPL%MSPDL=10
    YDML_LBC%YYTGMVCPL%MSPDM=11
    YDML_LBC%YYTGMVCPL%NDIMT=11
  ELSEIF (YDDYNA%LNHDYN.AND.(.NOT.YDDYNA%LNHX)) THEN
    YDML_LBC%YYTGMVCPL%MU   = 1
    YDML_LBC%YYTGMVCPL%MV   = 2
    YDML_LBC%YYTGMVCPL%MT   = 3
    YDML_LBC%YYTGMVCPL%MSPD = 4
    YDML_LBC%YYTGMVCPL%MSVD = 5
    YDML_LBC%YYTGMVCPL%MNHX = 5
    YDML_LBC%YYTGMVCPL%NDIM = 5
    YDML_LBC%YYTGMVCPL%MDIV = 6
    YDML_LBC%YYTGMVCPL%MTL  = 7
    YDML_LBC%YYTGMVCPL%MTM  = 8
    YDML_LBC%YYTGMVCPL%MSPDL= 9
    YDML_LBC%YYTGMVCPL%MSPDM=10
    YDML_LBC%YYTGMVCPL%NDIMT=10
  ELSE
    YDML_LBC%YYTGMVCPL%MU   = 1
    YDML_LBC%YYTGMVCPL%MV   = 2
    YDML_LBC%YYTGMVCPL%MT   = 3
    YDML_LBC%YYTGMVCPL%MSPD = 3
    YDML_LBC%YYTGMVCPL%MSVD = 3
    YDML_LBC%YYTGMVCPL%MNHX = 3
    YDML_LBC%YYTGMVCPL%NDIM = 3
    YDML_LBC%YYTGMVCPL%MDIV = 4
    YDML_LBC%YYTGMVCPL%MTL  = 5
    YDML_LBC%YYTGMVCPL%MTM  = 6
    YDML_LBC%YYTGMVCPL%MSPDL= 6
    YDML_LBC%YYTGMVCPL%MSPDM= 6
    YDML_LBC%YYTGMVCPL%NDIMT= 6
  ENDIF
ENDIF

! YYGMVSCPL:
IF (NCONF == 701) THEN
  ! derivatives are useless because ESEIMPLS is not called.
  YDML_LBC%YYTGMVSCPL%MSP=1
  YDML_LBC%YYTGMVSCPL%NDIM=1
  YDML_LBC%YYTGMVSCPL%NDIMT=1
ELSE
  ! derivatives are useful because ESEIMPLS is called.
  YDML_LBC%YYTGMVSCPL%MSP=1
  YDML_LBC%YYTGMVSCPL%MSPL=2
  YDML_LBC%YYTGMVSCPL%MSPM=3
  YDML_LBC%YYTGMVSCPL%NDIM=1
  YDML_LBC%YYTGMVSCPL%NDIMT=3
ENDIF

! YDML_LBC%NDIMCPL:
ICPL=0
DO JGFL=1,YGFL%NUMFLDS
  IF (YGFL%YCOMP(JGFL)%NCOUPLING /= 0) THEN
    ICPL=ICPL+1
  ENDIF
ENDDO
YDML_LBC%NDIMCPL=ICPL

! YDML_LBC%NGALEF:
YDML_LBC%NGALEF=YDML_LBC%YYTGMVCPL%NDIM+YDML_LBC%YYTGMVSCPL%NDIM+YDML_LBC%NDIMCPL


!--------------------------------------------------------------------------

!      Part B: printings.

WRITE(NULOUT,*) ' --- PRINTINGS IN SUELBC_FIELDS_DIM --- '

! * B1: Number of coupled fields.
WRITE(UNIT=NULOUT,FMT='('' Grid-point coupling: '')')
WRITE(UNIT=NULOUT,FMT='(''  nb of GMV fields with temporal interpolation: YDML_LBC%YYTGMVCPL%NDIMT = '',I4)')&
  & YDML_LBC%YYTGMVCPL%NDIMT
WRITE(UNIT=NULOUT,FMT='(''  nb of coupled GMV fields: YDML_LBC%YYTGMVCPL%NDIM = '',I4)') YDML_LBC%YYTGMVCPL%NDIM
WRITE(UNIT=NULOUT,FMT='(''  nb of GMVS fields with temporal interpolation: YDML_LBC%YYTGMVSCPL%NDIMT = '',I4)')&
 & YDML_LBC%YYTGMVSCPL%NDIMT
WRITE(UNIT=NULOUT,FMT='(''  nb of coupled GMVS fields: YDML_LBC%YYTGMVSCPL%NDIM = '',I4)') YDML_LBC%YYTGMVSCPL%NDIM
WRITE(UNIT=NULOUT,FMT='(''  nb of coupled GFL fields: YDML_LBC%NDIMCPL = '',I4)') YDML_LBC%NDIMCPL
WRITE(UNIT=NULOUT,FMT='(''  total nb of coupled fields: YDML_LBC%NGALEF = '',I4)') YDML_LBC%NGALEF


!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YEMLBC_FIELDS:SUELBC_FIELDS_DIM',1,ZHOOK_HANDLE)
END SUBROUTINE SUELBC_FIELDS_DIM 
!=====================================================================================
SUBROUTINE SUELBC_INIT(YDDYNA,YDML_LBC)
USE YOMDYNA      , ONLY : TDYNA

!--------------------------------------------------------------------------
! Sets-up part 0A of forcing a LAM model by another model
!--------------------------------------------------------------------------
TYPE(TDYNA)       , INTENT(IN)    :: YDDYNA
TYPE(TELBC_MODEL),TARGET,INTENT(INOUT):: YDML_LBC

!--------------------------------------------------------------------------

INTEGER(KIND=JPIM) :: J
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!----POINTERS FOR ASSOCIATION --------------------------
LOGICAL,POINTER :: LTENC,LALLTC,LRFIRST
INTEGER(KIND=JPIM),POINTER :: NBICOU,NBICOT,NBICPD,NBICVD,NBICNHX,NBICOP,NECRIPL
LOGICAL,POINTER :: LQCPL,LCCPL
INTEGER(KIND=JPIM),POINTER :: NECOTL,NECOAD
LOGICAL,POINTER :: LE0COTA,LEREADINI
INTEGER(KIND=JPIM),POINTER :: N1LSG,NFRLSG,NLSGTS(:)
LOGICAL,POINTER :: LRDLSG
LOGICAL,POINTER :: LESPCPL,LSPTENC,LUNBC
!--------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "posnam.intfb.h"

#include "nemelbc0a.nam.h"


!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YEMLBC_INIT:SUELBC_INIT',0,ZHOOK_HANDLE)
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------

LTENC=>YDML_LBC%LTENC
LALLTC=>YDML_LBC%LALLTC
LRFIRST=>YDML_LBC%LRFIRST
NBICOU=>YDML_LBC%NBICOU
NBICOT=>YDML_LBC%NBICOT
NBICPD=>YDML_LBC%NBICPD
NBICVD=>YDML_LBC%NBICVD
NBICNHX=>YDML_LBC%NBICNHX
NBICOP=>YDML_LBC%NBICOP
NECRIPL=>YDML_LBC%NECRIPL
LQCPL=>YDML_LBC%LQCPL
LCCPL=>YDML_LBC%LCCPL
NECOTL=>YDML_LBC%NECOTL
NECOAD=>YDML_LBC%NECOAD
LE0COTA=>YDML_LBC%LE0COTA
LEREADINI=>YDML_LBC%LEREADINI
N1LSG=>YDML_LBC%N1LSG
NFRLSG=>YDML_LBC%NFRLSG
NLSGTS=>YDML_LBC%NLSGTS
LRDLSG=>YDML_LBC%LRDLSG
LESPCPL=>YDML_LBC%LESPCPL
LSPTENC=>YDML_LBC%LSPTENC
LUNBC=>YDML_LBC%LUNBC
!      Part A: default values.

! * Tendency coupling
LTENC=.FALSE.
LALLTC=.FALSE.
LRFIRST=.TRUE.

! * Lateral forcing
IF (LELAM) THEN
  NBICOU=1
  NBICOT=1
  IF (YDDYNA%LNHDYN) THEN
    NBICPD=1
    NBICVD=1
    IF(YDDYNA%LNHX) THEN
      NBICNHX=1
    ELSE
      NBICNHX=0
    ENDIF
  ELSE
    NBICPD=0
    NBICVD=0
    NBICNHX=0
  ENDIF
  NBICOP=1
  NECRIPL=1
  LQCPL=.FALSE.
  LCCPL=.FALSE.
  NECOAD=0
  NECOTL=0
  LE0COTA=.FALSE.
  LEREADINI=.TRUE.
  N1LSG=0
  NFRLSG=1
  DO J=0,JPLSGT
    NLSGTS(J)=0
  ENDDO
  LRDLSG=.FALSE.
ELSE
  NBICOU=0
  NBICOT=0
  NBICOP=0
  NBICPD=0
  NBICVD=0
  NBICNHX=0
  NECRIPL=1
  LQCPL=.FALSE.
  LCCPL=.FALSE.
  NECOAD=0
  NECOTL=0
  LE0COTA=.FALSE.
  LEREADINI=.FALSE.
  N1LSG=0
  NFRLSG=1
  DO J=0,JPLSGT
    NLSGTS(J)=0
  ENDDO
  LRDLSG=.FALSE.
ENDIF

! * Spectral nudging
LESPCPL=.FALSE.
LSPTENC=.FALSE.

! * Upper nesting boundary conditions
LUNBC=.FALSE.


!--------------------------------------------------------------------------

!      Part B: namelist reading.

IF (LELAM) THEN
  CALL POSNAM(NULNAM,'NEMELBC0A')
  READ(NULNAM,NEMELBC0A)
ENDIF

!--------------------------------------------------------------------------

!      Part C: checkings.

! * checkings on tendency coupling
IF (.NOT.LTENC.AND.LALLTC) THEN
  CALL ABOR1('SUELBC_INIT: ABOR1 CALLED: .NOT.LTENC.AND.LALLTC')
ENDIF
IF (LTENC.AND.LRDLSG) THEN
  CALL ABOR1('SUELBC_INIT: ABOR1 CALLED: LTENC.AND.LRDLSG')
ENDIF

! * checkings on lateral boundary forcing
IF (NECRIPL /= 0.AND.NECRIPL /= 1.AND.NECRIPL /= 2) THEN
  CALL ABOR1('SUELBC_INIT: IMPROPER VALUE FOR NECRIPL')
ENDIF
IF (LQCPL.AND.NCONF == 701) THEN
  CALL ABOR1('SUELBC_INIT: LQCPL=.T. NOT PREPARED FOR NCONF=701')
ENDIF
IF (NECRIPL /= 0.AND.( NBICOU /= NBICOT.OR.NBICOU /= NBICOP ) ) THEN
  WRITE(NULOUT,*) 'T1 COUPLING DOESNT LET COUPLING FUNCTIONS DIFFER'
  CALL ABOR1 ('SUELBC_INIT: NECRIPL /= 0 and improper values of some NBIC..')
ENDIF

!--------------------------------------------------------------------------

!      Part D: printings.

IF (LELAM) THEN
  WRITE(NULOUT,*) ' --- PRINTINGS IN SUELBC_INIT --- '
  WRITE(NULOUT,*) ' LTENC = ',LTENC,' LALLTC = ',LALLTC
  WRITE(NULOUT,*) ' LSPTENC = ',LSPTENC
  WRITE(NULOUT,*) ' NBICOU = ',NBICOU,' NBICOT = ',NBICOT,' NBICOP = ',NBICOP
  WRITE(NULOUT,*) ' NBICPD = ',NBICPD,' NBICVD = ',NBICVD,&
   & ' NBICNHX = ',NBICNHX,' LQCPL = ',LQCPL,' LCCPL = ',LCCPL
  WRITE(UNIT=NULOUT,FMT='('' NECRIPL = '',I2)') NECRIPL
  WRITE(NULOUT,*) ' NECOAD = ',NECOAD,' NECOTL = ',NECOTL,' LE0COTA = ',LE0COTA
  WRITE(NULOUT,FMT='('' NFRLSG = '',I2,'' N1LSG = '',I2)') NFRLSG,N1LSG
  WRITE(NULOUT,*) ' NLSGTS =  ',NLSGTS(0),(NLSGTS(J),J=1,ABS(NLSGTS(0)))
  WRITE(NULOUT,FMT='('' LRDLSG = '',L2)') LRDLSG
  WRITE(NULOUT,*) 'LESPCPL = ',LESPCPL
  WRITE(NULOUT,*) 'LUNBC = ',LUNBC
  WRITE(NULOUT,*) ' '
ENDIF

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YEMLBC_INIT:SUELBC_INIT',1,ZHOOK_HANDLE)
END SUBROUTINE SUELBC_INIT

SUBROUTINE SUELBC_MODEL(YDML_LBC,YDGEOMETRY,YDDYNA,YDGMV,YDML_GCONF)

!--------------------------------------------------------------------------
! Sets-up part 0B of forcing a LAM model by another model
!--------------------------------------------------------------------------

!--------------------------------------------------------------------------

USE GEOMETRY_MOD           , ONLY : GEOMETRY
USE MODEL_GENERAL_CONF_MOD , ONLY : MODEL_GENERAL_CONF_TYPE
USE YOMDYNA                , ONLY : TDYNA

TYPE(TELBC_MODEL) ,TARGET,INTENT(INOUT)    :: YDML_LBC 
TYPE(GEOMETRY),INTENT(IN)                  :: YDGEOMETRY
TYPE(TDYNA)   , INTENT(IN)                 :: YDDYNA
TYPE(TGMV)    ,INTENT(INOUT)               :: YDGMV
TYPE(MODEL_GENERAL_CONF_TYPE),INTENT(INOUT):: YDML_GCONF

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
REAL(KIND=JPRB),ALLOCATABLE :: ZB(:,:),ZBU(:,:),ZEALP(:),ZREPA(:),ZEALFA(:,:,:)
REAL(KIND=JPRB),ALLOCATABLE :: ZGP(:,:,:),ZSPM(:,:)
INTEGER(KIND=JPIM),ALLOCATABLE :: INEAL(:),INNAL(:),INMAL(:)
INTEGER(KIND=JPIM),ALLOCATABLE :: IBICO(:)

INTEGER(KIND=JPIM) :: IENDLON, ISTLON, IBZONGL, IGPTOT
INTEGER(KIND=JPIM) :: IA, IGLG, IIA, IJA, IROF, ISUP, ICPL
INTEGER(KIND=JPIM) :: JFLD, JK, JA, JIA, JJA, JLON, JGL
INTEGER(KIND=JPIM) :: INSTEP, ICANCPL, JLEV
INTEGER(KIND=JPIM) :: IBL, JKGLO, JROF, IBLOCKCPL

REAL(KIND=JPRB) :: ZA, ZE, ZO, ZRZONG, ZRZONL, ZRZONU, ZXA, ZYA, ZDIV, ZREM, ZYAOXA
REAL(KIND=JPRB) :: ZT, ZT1, ZT2, ZT3, ZT4
REAL(KIND=JPRB) :: ZEPS=1.E-4_JPRB
REAL(KIND=JPRB) :: ZEPS2=1.E-10_JPRB
INTEGER(KIND=JPIM) :: INEFRCLDFI, INEFRCLDFIA, ISWP

!----POINTERS FOR ASSOCIATION --------------------------
INTEGER(KIND=JPIM),POINTER :: NEFRSPCPL,NEK0,NEK1,NEN1,NEN2
REAL(KIND=JPRB),POINTER :: EPA_GMV(:),EPA_GMVS(:),EPA_GFL(:)
REAL(KIND=JPRB),POINTER :: TEFRCL,SPNUDVOR,SPNUDDIV,SPNUDT,SPNUDQ,SPNUDSP,RNUTENC,RTENC


!--------------------------------------------------------------------------

#include "abor1.intfb.h"
#include "ereespe.intfb.h"
#include "esperee.intfb.h"
#include "posnam.intfb.h"
#include "fctez.func.h"

#include "nemelbc0b.nam.h"

!--------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YEMLBC_MODEL:SUELBC_MODEL',0,ZHOOK_HANDLE)


EPA_GMV=>YDML_LBC%EPA_GMV
EPA_GMVS=>YDML_LBC%EPA_GMVS
EPA_GFL=>YDML_LBC%EPA_GFL
TEFRCL=>YDML_LBC%TEFRCL
NEFRSPCPL=>YDML_LBC%NEFRSPCPL
NEK0=>YDML_LBC%NEK0
NEK1=>YDML_LBC%NEK1
NEN1=>YDML_LBC%NEN1
NEN2=>YDML_LBC%NEN2
SPNUDVOR=>YDML_LBC%SPNUDVOR
SPNUDDIV=>YDML_LBC%SPNUDDIV
SPNUDT=>YDML_LBC%SPNUDT
SPNUDQ=>YDML_LBC%SPNUDQ
SPNUDSP=>YDML_LBC%SPNUDSP
RNUTENC=>YDML_LBC%RNUTENC
RTENC=>YDML_LBC%RTENC

!--------------------------------------------------------------------------
ASSOCIATE(YDGSGEOM_NB=>YDGEOMETRY%YRGSGEOM_NB, &
 & NBZONG=>YDGEOMETRY%YREDIM%NBZONG, NBZONL=>YDGEOMETRY%YREDIM%NBZONL, NBZONU=>YDGEOMETRY%YREDIM%NBZONU,&
 & NBIPINCIX=>YDGEOMETRY%YREDIM%NBIPINCIX, NBIPINCIY=>YDGEOMETRY%YREDIM%NBIPINCIY, &
 & NDGLG=>YDGEOMETRY%YRDIM%NDGLG, NDLUXG=>YDGEOMETRY%YRDIM%NDLUXG, NDGUNG=>YDGEOMETRY%YRDIM%NDGUNG, &
 & NDGUXG=>YDGEOMETRY%YRDIM%NDGUXG, NDLUNG=>YDGEOMETRY%YRDIM%NDLUNG, NDGENL=>YDGEOMETRY%YRDIM%NDGENL, &
 & NPROMA=>YDGEOMETRY%YRDIM%NPROMA, NGPBLKS=>YDGEOMETRY%YRDIM%NGPBLKS, NSPEC2=>YDGEOMETRY%YRDIM%NSPEC2, &
 & NDLON=>YDGEOMETRY%YRDIM%NDLON, NDGSAL=>YDGEOMETRY%YRDIM%NDGSAL, &
 & NSTA=>YDGEOMETRY%YRMP%NSTA, NPTRFLOFF=>YDGEOMETRY%YRMP%NPTRFLOFF, MYLATS=>YDGEOMETRY%YRMP%MYLATS, &
 & NONL=>YDGEOMETRY%YRMP%NONL, NBSETSP=>YDGEOMETRY%YRMP%NBSETSP, &
 & NGPTOT=>YDGEOMETRY%YRGEM%NGPTOT, NSTAGP=>YDGEOMETRY%YRGEM%NSTAGP, &
 & NFD2D=>YDML_GCONF%YRDIMF%NFD2D, NS3D=>YDML_GCONF%YRDIMF%NS3D, NFLEVG=>YDGEOMETRY%YRDIMV%NFLEVG, &
 & NFLEVL=>YDGEOMETRY%YRDIMV%NFLEVL, YT0=>YDGMV%YT0, YT1=>YDGMV%YT1 , &
 & NSTOP=>YDML_GCONF%YRRIP%NSTOP, TSTEP=>YDML_GCONF%YRRIP%TSTEP)
!--------------------------------------------------------------------------

!      Part A: default values.

! * Relaxation coefficients
YDML_LBC%EPA_GMV(:)=2.16_JPRB
YDML_LBC%EPA_GMVS(:)=2.16_JPRB
YDML_LBC%EPA_GFL(:)=5.52_JPRB

! * Control frequency of LBC
TEFRCL=TSTEP

! * Spectral nudging
YDML_LBC%NEFRSPCPL=1
YDML_LBC%NEK0=0
YDML_LBC%NEK1=0
YDML_LBC%NEN1=0
YDML_LBC%NEN2=0
YDML_LBC%SPNUDVOR=0._JPRB
YDML_LBC%SPNUDDIV=0._JPRB
YDML_LBC%SPNUDT=0._JPRB
YDML_LBC%SPNUDQ=0._JPRB
YDML_LBC%SPNUDSP=0._JPRB
YDML_LBC%RNUTENC=0.0_JPRB

! * Tendency coupling
YDML_LBC%RTENC=1.0_JPRB

!--------------------------------------------------------------------------

!      Part B: namelist reading.

CALL POSNAM(NULNAM,'NEMELBC0B')
READ(NULNAM,NEMELBC0B)

!--------------------------------------------------------------------------

!      Part C: checkings, modify values.

! * C1: Control frequency of LBC (compute NEFRCL)

! TEFRCL, NEFRCL:
IF (TSTEP > ZEPS) THEN 
  ZDIV = TEFRCL/TSTEP
  ZREM = ZDIV - REAL(NINT(ZDIV),JPRB)
  ! 1 second/coupling-interval error is tolerated
  IF(ABS(ZREM) >= 1._JPRB/TSTEP) THEN
    WRITE(UNIT=NULOUT,FMT='('' TEFRCL MUST BE A MULTIPLE '',''OF TSTEP'')')
    CALL ABOR1('SUELBC_MODEL: TEFRCL MUST BE A MULTIPLE OF TSTEP')
  ELSE
    YDML_LBC%NEFRCL=NINT(TEFRCL/TSTEP)
  ENDIF
ELSE
  YDML_LBC%NEFRCL=1
ENDIF

! * C2: Weights for LBC interpolation (compute EWB, EWBDFIFW, EWBDFIBW)

! Fill EWB:
IF ((NSTOP/=0).AND.(YDML_LBC%NEFRCL/=0)) THEN

  ALLOCATE(YDML_LBC%EWB(0:NSTOP,1:4,0:9))
  YDML_LBC%EWB(:,:,:)=0.0_JPRB

  IF (YDML_LBC%LQCPL) THEN
!DEC$ IVDEP
    DO INSTEP=0,NSTOP
      IF (INSTEP<YDML_LBC%NEFRCL) THEN
        ! KTIMLEV=0
        YDML_LBC%EWB(INSTEP,1,0)=REAL(YDML_LBC%NEFRCL-MOD(INSTEP,YDML_LBC%NEFRCL),JPRB)/REAL(YDML_LBC%NEFRCL,JPRB)
        YDML_LBC%EWB(INSTEP,2,0)=REAL(MOD(INSTEP,YDML_LBC%NEFRCL),JPRB)/REAL(YDML_LBC%NEFRCL,JPRB)
        YDML_LBC%EWB(INSTEP,3,:)=0.0_JPRB

        ! KTIMLEV=1
        YDML_LBC%EWB(INSTEP,1,1)=YDML_LBC%EWB(INSTEP,1,0)-1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)
        YDML_LBC%EWB(INSTEP,2,1)=YDML_LBC%EWB(INSTEP,2,0)+1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)

        ! KTIMLEV=9
        YDML_LBC%EWB(INSTEP,1,9)=YDML_LBC%EWB(INSTEP,1,0)+1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)
        YDML_LBC%EWB(INSTEP,2,9)=YDML_LBC%EWB(INSTEP,2,0)-1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)
      ELSE
        ISWP=MOD(INSTEP,YDML_LBC%NEFRCL)+1
        ZT  =REAL(              ISWP,JPRB)
        ZT1 =REAL(-1*YDML_LBC%NEFRCL,JPRB)
        ZT2 =REAL( 0*YDML_LBC%NEFRCL,JPRB)
        ZT3 =REAL(+1*YDML_LBC%NEFRCL,JPRB)

        YDML_LBC%EWB(INSTEP,1,1)=(ZT -ZT2)*(ZT -ZT3)/ &
         &                      ((ZT1-ZT2)*(ZT1-ZT3))
        YDML_LBC%EWB(INSTEP,2,1)=(ZT -ZT1)*(ZT -ZT3)/ &
         &                      ((ZT2-ZT1)*(ZT2-ZT3))
        YDML_LBC%EWB(INSTEP,3,1)=(ZT -ZT1)*(ZT -ZT2)/ &
         &                      ((ZT3-ZT1)*(ZT3-ZT2))
      ENDIF
    ENDDO
  ELSEIF (YDML_LBC%LCCPL) THEN
!DEC$ IVDEP
    DO INSTEP=0,NSTOP
      IF (INSTEP<YDML_LBC%NEFRCL) THEN
        ! KTIMLEV=0
        YDML_LBC%EWB(INSTEP,1,0)=REAL(YDML_LBC%NEFRCL-MOD(INSTEP,YDML_LBC%NEFRCL),JPRB)/REAL(YDML_LBC%NEFRCL,JPRB)
        YDML_LBC%EWB(INSTEP,2,0)=REAL(MOD(INSTEP,YDML_LBC%NEFRCL),JPRB)/REAL(YDML_LBC%NEFRCL,JPRB)
        YDML_LBC%EWB(INSTEP,3,:)=0.0_JPRB
        YDML_LBC%EWB(INSTEP,4,:)=0.0_JPRB

        ! KTIMLEV=1
        YDML_LBC%EWB(INSTEP,1,1)=YDML_LBC%EWB(INSTEP,1,0)-1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)
        YDML_LBC%EWB(INSTEP,2,1)=YDML_LBC%EWB(INSTEP,2,0)+1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)


        ! KTIMLEV=9
        YDML_LBC%EWB(INSTEP,1,9)=YDML_LBC%EWB(INSTEP,1,0)+1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)
        YDML_LBC%EWB(INSTEP,2,9)=YDML_LBC%EWB(INSTEP,2,0)-1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)

      ELSEIF (INSTEP<2*YDML_LBC%NEFRCL) THEN
        ISWP=MOD(INSTEP,YDML_LBC%NEFRCL)+1
        ZT  =REAL(              ISWP,JPRB)
        ZT1 =REAL(-1*YDML_LBC%NEFRCL,JPRB)
        ZT2 =REAL( 0*YDML_LBC%NEFRCL,JPRB)
        ZT3 =REAL(+1*YDML_LBC%NEFRCL,JPRB)

        YDML_LBC%EWB(INSTEP,1,1)=(ZT -ZT2)*(ZT -ZT3)/ &
         &                      ((ZT1-ZT2)*(ZT1-ZT3))
        YDML_LBC%EWB(INSTEP,2,1)=(ZT -ZT1)*(ZT -ZT3)/ &
         &                      ((ZT2-ZT1)*(ZT2-ZT3))
        YDML_LBC%EWB(INSTEP,3,1)=(ZT -ZT1)*(ZT -ZT2)/ &
         &                      ((ZT3-ZT1)*(ZT3-ZT2))
      ELSE
        ISWP=MOD(INSTEP,YDML_LBC%NEFRCL)+1
        ZT  =REAL(ISWP,JPRB)
        IF (INSTEP >= NSTOP-YDML_LBC%NEFRCL) THEN
          ZT1=REAL(-1*YDML_LBC%NEFRCL,JPRB)
          ZT2=REAL( 0*YDML_LBC%NEFRCL,JPRB)
          ZT3=REAL(+1*YDML_LBC%NEFRCL,JPRB)

          YDML_LBC%EWB(INSTEP,1,1)=0.0_JPRB
          YDML_LBC%EWB(INSTEP,2,1)=(ZT -ZT2)*(ZT -ZT3)/ &
           &                      ((ZT1-ZT2)*(ZT1-ZT3))
          YDML_LBC%EWB(INSTEP,3,1)=(ZT -ZT1)*(ZT -ZT3)/ &
           &                      ((ZT2-ZT1)*(ZT2-ZT3))
          YDML_LBC%EWB(INSTEP,4,1)=(ZT -ZT1)*(ZT -ZT2)/ &
           &                      ((ZT3-ZT1)*(ZT3-ZT2))
        ELSE
          ZT1=REAL(-1*YDML_LBC%NEFRCL,JPRB)
          ZT2=REAL( 0*YDML_LBC%NEFRCL,JPRB)
          ZT3=REAL(+1*YDML_LBC%NEFRCL,JPRB)
          ZT4=REAL(+2*YDML_LBC%NEFRCL,JPRB)

          YDML_LBC%EWB(INSTEP,1,1)=(ZT -ZT2)*(ZT -ZT3)*(ZT -ZT4)/ &
           &                      ((ZT1-ZT2)*(ZT1-ZT3)*(ZT1-ZT4))
          YDML_LBC%EWB(INSTEP,2,1)=(ZT -ZT1)*(ZT -ZT3)*(ZT -ZT4)/ &
           &                      ((ZT2-ZT1)*(ZT2-ZT3)*(ZT2-ZT4))
          YDML_LBC%EWB(INSTEP,3,1)=(ZT -ZT1)*(ZT -ZT2)*(ZT -ZT4)/ &
           &                      ((ZT3-ZT1)*(ZT3-ZT2)*(ZT3-ZT4))
          YDML_LBC%EWB(INSTEP,4,1)=(ZT -ZT1)*(ZT -ZT2)*(ZT -ZT3)/ &
           &                      ((ZT4-ZT1)*(ZT4-ZT2)*(ZT4-ZT3))
        ENDIF
      ENDIF
    ENDDO
  ELSE
!DEC$ IVDEP
    DO INSTEP=0,NSTOP

      ! KTIMLEV=0
      YDML_LBC%EWB(INSTEP,1,0)=REAL(YDML_LBC%NEFRCL-MOD(INSTEP,YDML_LBC%NEFRCL),JPRB)/REAL(YDML_LBC%NEFRCL,JPRB)
      YDML_LBC%EWB(INSTEP,2,0)=REAL(MOD(INSTEP,YDML_LBC%NEFRCL),JPRB)/REAL(YDML_LBC%NEFRCL,JPRB)

      ! KTIMLEV=1
      YDML_LBC%EWB(INSTEP,1,1)=YDML_LBC%EWB(INSTEP,1,0)-1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)
      YDML_LBC%EWB(INSTEP,2,1)=YDML_LBC%EWB(INSTEP,2,0)+1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)

      ! KTIMLEV=9
      YDML_LBC%EWB(INSTEP,1,9)=YDML_LBC%EWB(INSTEP,1,0)+1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)
      YDML_LBC%EWB(INSTEP,2,9)=YDML_LBC%EWB(INSTEP,2,0)-1.0_JPRB/REAL(YDML_LBC%NEFRCL,JPRB)

    ENDDO

  ENDIF
ENDIF

IF (LDFI) THEN

  ! Fill EWBDFIFW (DFI forward weights):
  ALLOCATE(YDML_LBC%EWBDFIFW(0:2*NSTDFI,1:2,0:9,0:1))
  INEFRCLDFI=YDML_LBC%TEFRCL/RTDFI

  DO INSTEP=0,2*NSTDFI
    ! KTIMLEV=0
    YDML_LBC%EWBDFIFW(INSTEP,1,0,:)=REAL(INEFRCLDFI-MOD(INSTEP,INEFRCLDFI)+NSTDFI,JPRB)/&
     & REAL(INEFRCLDFI,JPRB)
    YDML_LBC%EWBDFIFW(INSTEP,2,0,:)=REAL(MOD(INSTEP,INEFRCLDFI)-NSTDFI,JPRB)/&
     & REAL(INEFRCLDFI,JPRB)
    ! KTIMLEV=1, LBIAS=F
    YDML_LBC%EWBDFIFW(INSTEP,1,1,0)=YDML_LBC%EWBDFIFW(INSTEP,1,0,0)-1.0_JPRB/REAL(INEFRCLDFI,JPRB)
    YDML_LBC%EWBDFIFW(INSTEP,2,1,0)=YDML_LBC%EWBDFIFW(INSTEP,2,0,0)+1.0_JPRB/REAL(INEFRCLDFI,JPRB)
    ! KTIMLEV=1, LBIAS=T
    YDML_LBC%EWBDFIFW(INSTEP,1,1,1)=REAL(INEFRCLDFI+MOD(INSTEP,INEFRCLDFI)-NSTDFI,JPRB)/&
     & REAL(INEFRCLDFI,JPRB)+1.0_JPRB/REAL(INEFRCLDFI,JPRB)
    YDML_LBC%EWBDFIFW(INSTEP,2,1,1)=REAL(-MOD(INSTEP,INEFRCLDFI)+NSTDFI,JPRB)/&
     & REAL(INEFRCLDFI,JPRB)-1.0_JPRB/REAL(INEFRCLDFI,JPRB)
    ! KTIMLEV=9, LBIAS=F
    YDML_LBC%EWBDFIFW(INSTEP,1,9,0)=YDML_LBC%EWBDFIFW(INSTEP,1,0,0)+1.0_JPRB/REAL(INEFRCLDFI,JPRB)
    YDML_LBC%EWBDFIFW(INSTEP,2,9,0)=YDML_LBC%EWBDFIFW(INSTEP,2,0,0)-1.0_JPRB/REAL(INEFRCLDFI,JPRB)
    ! KTIMLEV=9, LBIAS=T
    YDML_LBC%EWBDFIFW(INSTEP,1,9,1)=REAL(INEFRCLDFI+MOD(INSTEP,INEFRCLDFI)-NSTDFI,JPRB)/&
     & REAL(INEFRCLDFI,JPRB)-1.0_JPRB/REAL(INEFRCLDFI,JPRB)
    YDML_LBC%EWBDFIFW(INSTEP,2,9,1)=REAL(-MOD(INSTEP,INEFRCLDFI)+NSTDFI,JPRB)/&
     & REAL(INEFRCLDFI,JPRB)+1.0_JPRB/REAL(INEFRCLDFI,JPRB)
  ENDDO

  ! Fill EWBDFIBW (DFI backward weights):
  ALLOCATE(YDML_LBC%EWBDFIBW(0:2*NSTDFIA,1:2,0:9,0:1))
  INEFRCLDFIA=YDML_LBC%TEFRCL/RTDFIA

  DO INSTEP=0,2*NSTDFIA
    YDML_LBC%EWBDFIBW(INSTEP,1,0,:)=REAL(INEFRCLDFIA+MOD(INSTEP,INEFRCLDFIA),JPRB)/REAL(INEFRCLDFIA,JPRB)
    YDML_LBC%EWBDFIBW(INSTEP,2,0,:)=-REAL(MOD(INSTEP,INEFRCLDFIA),JPRB)/REAL(INEFRCLDFIA,JPRB)

    ! KTIMLEV=1, LBIAS=F
    YDML_LBC%EWBDFIBW(INSTEP,1,1,0)=YDML_LBC%EWBDFIBW(INSTEP,1,0,0)+1.0_JPRB/REAL(INEFRCLDFIA,JPRB)
    YDML_LBC%EWBDFIBW(INSTEP,2,1,0)=YDML_LBC%EWBDFIBW(INSTEP,2,0,0)-1.0_JPRB/REAL(INEFRCLDFIA,JPRB)
    ! KTIMLEV=1, LBIAS=T
    YDML_LBC%EWBDFIBW(INSTEP,1,1,1)=REAL(INEFRCLDFIA-MOD(INSTEP,INEFRCLDFIA),JPRB)/REAL(INEFRCLDFIA,JPRB)&
     & -1.0_JPRB/REAL(INEFRCLDFIA,JPRB)
    YDML_LBC%EWBDFIBW(INSTEP,2,1,1)=REAL(MOD(INSTEP,INEFRCLDFIA),JPRB)/REAL(INEFRCLDFIA,JPRB)&
     & +1.0_JPRB/REAL(INEFRCLDFIA,JPRB)
    ! KTIMLEV=9, LBIAS=F
    YDML_LBC%EWBDFIBW(INSTEP,1,9,0)=YDML_LBC%EWBDFIBW(INSTEP,1,0,0)-1.0_JPRB/REAL(INEFRCLDFIA,JPRB)
    YDML_LBC%EWBDFIBW(INSTEP,2,9,0)=YDML_LBC%EWBDFIBW(INSTEP,2,0,0)+1.0_JPRB/REAL(INEFRCLDFIA,JPRB)
    ! KTIMLEV=9, LBIAS=T
    YDML_LBC%EWBDFIBW(INSTEP,1,9,1)=REAL(INEFRCLDFIA-MOD(INSTEP,INEFRCLDFIA),JPRB)/REAL(INEFRCLDFIA,JPRB)&
     & +1.0_JPRB/REAL(INEFRCLDFIA,JPRB)
    YDML_LBC%EWBDFIBW(INSTEP,2,9,1)=REAL(MOD(INSTEP,INEFRCLDFIA),JPRB)/REAL(INEFRCLDFIA,JPRB)&
     & -1.0_JPRB/REAL(INEFRCLDFIA,JPRB)
  ENDDO

ENDIF

! * C3: Calculation of IBICO and LECOBI.

ALLOCATE(IBICO(YDML_LBC%NGALEF))
! GMV:
IBICO(YDML_LBC%YYTGMVCPL%MU)=YDML_LBC%NBICOU
IBICO(YDML_LBC%YYTGMVCPL%MV)=YDML_LBC%NBICOU
IBICO(YDML_LBC%YYTGMVCPL%MT)=YDML_LBC%NBICOT
IF (YDDYNA%LNHDYN) IBICO(YDML_LBC%YYTGMVCPL%MSPD)=YDML_LBC%NBICPD
IF (YDDYNA%LNHDYN) IBICO(YDML_LBC%YYTGMVCPL%MSVD)=YDML_LBC%NBICVD
IF (YDDYNA%LNHDYN.AND.YDDYNA%LNHX) IBICO(YDML_LBC%YYTGMVCPL%MNHX)=YDML_LBC%NBICNHX
! GMVS:
IBICO(YDML_LBC%YYTGMVCPL%NDIM+YDML_LBC%YYTGMVSCPL%MSP)=YDML_LBC%NBICOP
! GFL:
DO JFLD=1,YDML_LBC%NDIMCPL
  IBICO(YDML_LBC%YYTGMVCPL%NDIM+YDML_LBC%YYTGMVSCPL%NDIM+JFLD)=1
ENDDO

IF((NBZONL /= 0).OR.(NBZONG /= 0).OR.(NDLUXG /= NDLON).OR.(NDGUXG /= NDGLG)) THEN
  IF ( MAXVAL(IBICO(1:YDML_LBC%NGALEF))==0 .AND. MINVAL(IBICO(1:YDML_LBC%NGALEF))==0 ) THEN
    ! no field coupled; LECOBI set to F.
    YDML_LBC%LECOBI=.FALSE.
  ELSE
    ! at least one field coupled; non-empty coupling zone.
    YDML_LBC%LECOBI=.TRUE.
  ENDIF
ELSE
  ! empty coupling zone.
  YDML_LBC%LECOBI=.FALSE.
ENDIF
! For canari
ICANCPL=YDML_LBC%NBICOU+YDML_LBC%NBICOT+YDML_LBC%NBICOP
IF (LCANARI .AND. ICANCPL==0) YDML_LBC%LECOBI=.FALSE.

! * C4: Relaxation coefficients EALFA_GMV, EALFA_GMVS, EALFA_GFL, EALFA_TENC (former SUEBICU).

ALLOCATE(YDML_LBC%EALFA_GMV(NGPTOT+1,YDML_LBC%YYTGMVCPL%NDIM))
IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'EALFA_GMV ',&
 & SIZE(YDML_LBC%EALFA_GMV ),SHAPE(YDML_LBC%EALFA_GMV )

ALLOCATE(YDML_LBC%EALFA_GMVS(NGPTOT+1,YDML_LBC%YYTGMVSCPL%NDIM))
IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'EALFA_GMVS',&
 & SIZE(YDML_LBC%EALFA_GMVS),SHAPE(YDML_LBC%EALFA_GMVS)

ALLOCATE(YDML_LBC%EALFA_GFL(NGPTOT+1,YDML_LBC%NDIMCPL))
IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'EALFA_GFL ',&
& SIZE(YDML_LBC%EALFA_GFL ),SHAPE(YDML_LBC%EALFA_GFL )

IF(YDML_LBC%LUNBC) THEN
  ALLOCATE(YDML_LBC%EALFAU_GMV(NFLEVG,YDML_LBC%YYTGMVCPL%NDIM))
  IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'EALFAU_GMV ',&
   & SIZE(YDML_LBC%EALFAU_GMV ),SHAPE(YDML_LBC%EALFAU_GMV )

  ALLOCATE(YDML_LBC%EALFAU_GFL(NFLEVG,YDML_LBC%NDIMCPL))
  IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'EALFAU_GFL ',&
  & SIZE(YDML_LBC%EALFAU_GFL ),SHAPE(YDML_LBC%EALFAU_GFL )
ENDIF

IF (YDML_LBC%LTENC) THEN
  ALLOCATE(YDML_LBC%EALFA_TENC(NGPTOT+1,YDML_LBC%YYTGMVSTENC%NDIM))
  IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'EALFA_TENC',&
 & SIZE(YDML_LBC%EALFA_TENC),SHAPE(YDML_LBC%EALFA_TENC)
ENDIF

IF (YDML_LBC%LECOBI) THEN

  ! * C4.1: allocations.

  ALLOCATE(ZREPA(YDML_LBC%NGALEF))
  ALLOCATE(ZEALP(YDML_LBC%NGALEF))
  ALLOCATE(ZEALFA(NDLON,YDML_LBC%NGALEF,NDGLG))
  ALLOCATE(INEAL(YDML_LBC%NGALEF))
  ALLOCATE(INNAL(YDML_LBC%NGALEF))
  ALLOCATE(INMAL(YDML_LBC%NGALEF))
  ALLOCATE(ZB(NBZONL+1,NBZONG+1))
  IF (YDML_LBC%LUNBC) THEN
    IF (NBZONU > 0 .AND. NBZONU < NFLEVG) THEN
      ALLOCATE(ZBU(NBZONU+1,YDML_LBC%NGALEF))
    ELSE
      CALL ABOR1('SUELBC_MODEL: NBZONU IS OUT OF BOUNDS')
    ENDIF
  ENDIF

  ISUP=100
  IBZONGL=MAX(NBZONL-NBIPINCIX,NBZONG-NBIPINCIY)

  ! * C4.2: fill ZREPA.

  ICPL=0
  IF (SIZE(YDML_LBC%EPA_GMV) < YDML_LBC%YYTGMVCPL%NDIM) THEN
    CALL ABOR1('SUELBC_MODEL:YDML_LBC%EPA_GMV TOO SMALL !')
  ELSE
    ZREPA(ICPL+1:ICPL+YDML_LBC%YYTGMVCPL%NDIM)=YDML_LBC%EPA_GMV(1:YDML_LBC%YYTGMVCPL%NDIM)
  ENDIF
  ICPL=YDML_LBC%YYTGMVCPL%NDIM
  IF (SIZE(YDML_LBC%EPA_GMVS) < YDML_LBC%YYTGMVSCPL%NDIM) THEN
    CALL ABOR1('SUELBC_MODEL:YDML_LBC%EPA_GMVS TOO SMALL !')
  ELSE
    ZREPA(ICPL+1:ICPL+YDML_LBC%YYTGMVSCPL%NDIM)=YDML_LBC%EPA_GMVS(1:YDML_LBC%YYTGMVSCPL%NDIM)
  ENDIF
  ICPL=YDML_LBC%YYTGMVCPL%NDIM+YDML_LBC%YYTGMVSCPL%NDIM
  IF (SIZE(YDML_LBC%EPA_GFL) < YDML_LBC%NDIMCPL) THEN
    CALL ABOR1('SUELBC_MODEL:YDML_LBC%EPA_GFL TOO SMALL !')
  ELSE
    ZREPA(ICPL+1:ICPL+YDML_LBC%NDIMCPL)=YDML_LBC%EPA_GFL(1:YDML_LBC%NDIMCPL)
  ENDIF

  DO JFLD =1,YDML_LBC%NGALEF
    IF((ZREPA(JFLD) > -2.0_JPRB).AND.(ZREPA(JFLD) < 2.0_JPRB))THEN
      WRITE(UNIT=NULOUT,FMT='("ERROR ZREPA CANT BE",F5.2)')ZREPA(JFLD)
      CALL ABOR1('SUELBC_MODEL: ZREPA MUST NOT HAVE A VALUE STRICTLY BETWEEN -2.0 AND 2.0')
    ENDIF
  ENDDO

  ! * C4.3: fill INEAL,INMAL,INNAL (identical for all coupled fields).

  INEAL(1:YDML_LBC%NGALEF)=2
  INMAL(1:YDML_LBC%NGALEF)=1

  IF (IBZONGL>=11 .AND. IBZONGL<=26) THEN
    INNAL(1:YDML_LBC%NGALEF)=2
    WRITE(NULOUT,*) 'INNAL FORCED TO 2 FOR CONVERGENCE'
  ELSEIF (IBZONGL>=27) THEN
    INNAL(1:YDML_LBC%NGALEF)=1
    WRITE(NULOUT,*) 'INNAL FORCED TO 1 FOR CONVERGENCE'
  ELSE
    INNAL(1:YDML_LBC%NGALEF)=3
    WRITE(NULOUT,*) 'INNAL SET TO 3'
  ENDIF

  ! * C4.4: compute auxilary variables ZRZONL, ZRZONG, ZRZONU, ZEALP (identical for all
  ! coupled fields).

  IF (NBZONL > NBIPINCIX) THEN
    ZRZONL=1.0_JPRB/REAL(NBZONL-NBIPINCIX,JPRB)
  ENDIF
  IF (NBZONG > NBIPINCIY) THEN
    ZRZONG=1.0_JPRB/REAL(NBZONG-NBIPINCIY,JPRB)
  ENDIF

  IF (YDML_LBC%LUNBC) THEN
    ZRZONU = 1.0_JPRB/(YDGEOMETRY%YRVETA%VETAF(NBZONU+1)-YDGEOMETRY%YRVETA%VETAF(1))
  ENDIF

  DO JFLD =1,YDML_LBC%NGALEF
    ZEALP(JFLD)=&
     & REAL((INMAL(JFLD)+INNAL(JFLD))**(INMAL(JFLD)+INNAL(JFLD)),JPRB)/&
     & REAL((INNAL(JFLD)**INNAL(JFLD))*(INMAL(JFLD)**INMAL(JFLD))*INEAL(JFLD))
  ENDDO

  ! * C4.5: compute ZEALFA.

  DO JFLD =1,YDML_LBC%NGALEF

    IF (IBICO(JFLD) == 0) THEN

      ! --- no coupling applied to this field; we simply set ZEALFA=0 everywhere.
      ZEALFA(1:NDLON,JFLD,1:NDGLG)=0.0_JPRB

    ELSE

      ! --- coupling applied to this field.

      ! * ZEALFA: initialize the center domain to 0. and the outer domain to 1.

      ZEALFA(NDLUNG+NBZONL:NDLUXG-NBZONL,JFLD,NDGUNG+NBZONG:NDGUXG-NBZONG)=0.0_JPRB
      ZEALFA(1:NDLON,JFLD,1:NDGUNG-1)=1.0_JPRB
      ZEALFA(1:NDLON,JFLD,NDGUXG+1:NDGLG)=1.0_JPRB
      ZEALFA(1:NDLUNG-1,JFLD,NDGUNG:NDGUXG)=1.0_JPRB
      ZEALFA(NDLUXG+1:NDLON,JFLD,NDGUNG:NDGUXG)=1.0_JPRB

      ! * compute ZEALFA in the relaxation area.

      IF ((NBZONL > 0).OR.(NBZONG > 0))THEN

        ! Compute ZB:
        IF (NBZONL > 0) THEN
          DO JA=2,NBZONL
            ! relaxation function is 1 in 1:NBIPINCIX
            IF (JA<=NBZONL-NBIPINCIX) THEN
              ZA=REAL(JA-1,JPRB)*ZRZONL
              IF(ZREPA(JFLD) <= -2.0_JPRB) THEN
                ZB(JA,NBZONG+1)=FEZBM(ZA,-ZREPA(JFLD))
              ELSE
                ZB(JA,NBZONG+1)=FEZBP(ZA,ZREPA(JFLD))
              ENDIF
            ELSE
              ZB(JA,NBZONG+1)=1._JPRB
            ENDIF
          ENDDO
        ENDIF

        IF (NBZONG > 0) THEN
          DO JA=2,NBZONG
            ! relaxation function is 1 in 1:NBIPINCIY
            IF (JA<=NBZONG-NBIPINCIY) THEN
              ZA=REAL(JA-1,JPRB)*ZRZONG
              IF(ZREPA(JFLD) <= -2.0_JPRB) THEN
                ZB(NBZONL+1,JA)=FEZBM(ZA,-ZREPA(JFLD))
              ELSE
                ZB(NBZONL+1,JA)=FEZBP(ZA,ZREPA(JFLD))
              ENDIF
            ELSE
              ZB(NBZONL+1,JA)=1._JPRB
            ENDIF
          ENDDO
        ENDIF

        IF ((NBZONL > 0).AND.(NBZONG > 0)) THEN
          DO JIA=2,NBZONL
            IF (JIA<=NBZONL-NBIPINCIX) THEN
              ZXA=REAL(JIA-1,JPRB)*ZRZONL
              DO JJA=2,NBZONG
                IF (JJA<=NBZONG-NBIPINCIY) THEN
                  ZYA=REAL(JJA-1,JPRB)*ZRZONG
                  ZYAOXA=ZYA/ZXA
                  ZA=MAX(ZXA,ZYA)
                  ZO=ZA
                  DO JK=1,ISUP
                    ZE=FEZE(ZA,ZEALP(JFLD),INNAL(JFLD),INMAL(JFLD))
                    IF (JPRB == JPRD) THEN
                      ZA=(ZXA**ZE+ZYA**ZE)**(1.0_JPRB/ZE)
                    ELSE
                      ZA=ZXA*(1._JPRB+ZYAOXA**ZE)**(1.0_JPRB/ZE)
                    ENDIF
                    ZT=ABS(ZA-ZO)/ZO
                    IF (ZT <= ZEPS) EXIT
                    IF (JK == ISUP) THEN
                      WRITE(NULOUT,*) 'NO CONVERGENCE FOR EALFA'
                      CALL ABOR1('SUELBC_MODEL: NO CONVERGENCE FOR EALFA')
                    ENDIF
                    ZA=0.5_JPRB*(ZA+ZO)
                    ZO=ZA
                  ENDDO
                  IF(ZREPA(JFLD) <= -2.0_JPRB) THEN
                    ZB(JIA,JJA)=FEZBM(ZA,-ZREPA(JFLD))
                  ELSE
                    ZB(JIA,JJA)=FEZBP(ZA,ZREPA(JFLD))
                  ENDIF
                ELSE
                  ZB(JIA,JJA)=1._JPRB
                ENDIF
              ENDDO
            ELSE
              DO JJA=2,NBZONG
                ZB(JIA,JJA)=1._JPRB
              ENDDO
            ENDIF
          ENDDO
        ENDIF

        ! Initialize ZEALFA on the relaxation area

        IF (NBZONG > 0)THEN
          DO JLON=NDLUNG+NBZONL,NDLUXG-NBZONL
            ZEALFA(JLON,JFLD,NDGUNG)=1.0_JPRB
            DO JGL=NDGUNG+1,NDGUNG+NBZONG-1
              IA=NDGUNG+NBZONG+1-JGL
              ZEALFA(JLON,JFLD,JGL)=ZB(NBZONL+1,IA)
            ENDDO
            ZEALFA(JLON,JFLD,NDGUXG)=1.0_JPRB
            DO JGL=NDGUXG-NBZONG+1,NDGUXG-1
              IA=JGL-NDGUXG+NBZONG+1
              ZEALFA(JLON,JFLD,JGL)=ZB(NBZONL+1,IA)
            ENDDO
          ENDDO
        ENDIF

        IF(NBZONL > 0)THEN
          DO JGL=NDGUNG+NBZONG,NDGUXG-NBZONG
            ZEALFA(NDLUNG,JFLD,JGL)=1.0_JPRB
            DO JLON=NDLUNG+1,NDLUNG+NBZONL-1
              IA=NDLUNG+NBZONL+1-JLON
              ZEALFA(JLON,JFLD,JGL)=ZB(IA,NBZONG+1)
            ENDDO
            ZEALFA(NDLUXG,JFLD,JGL)=1.0_JPRB
            DO JLON=NDLUXG-NBZONL+1,NDLUXG-1
              IA=JLON-NDLUXG+NBZONL+1
              ZEALFA(JLON,JFLD,JGL)=ZB(IA,NBZONG+1)
            ENDDO
          ENDDO
        ENDIF

        IF((NBZONL > 0).AND.(NBZONG > 0)) THEN
          DO JLON=NDLUNG+1,NDLUNG+NBZONL-1
            IIA=NDLUNG+NBZONL+1-JLON
            DO JGL=NDGUNG+1,NDGUNG+NBZONG-1
              IJA=NDGUNG+NBZONG+1-JGL
              ZEALFA(JLON,JFLD,JGL)=ZB(IIA,IJA)
            ENDDO
            DO JGL=NDGUXG-NBZONG+1,NDGUXG-1
              IJA=JGL-NDGUXG+NBZONG+1
              ZEALFA(JLON,JFLD,JGL)=ZB(IIA,IJA)
            ENDDO
          ENDDO
          DO JLON=NDLUXG-NBZONL+1,NDLUXG-1
            IIA=JLON-NDLUXG+NBZONL+1
            DO JGL=NDGUNG+1,NDGUNG+NBZONG-1
              IJA=NDGUNG+NBZONG+1-JGL
              ZEALFA(JLON,JFLD,JGL)=ZB(IIA,IJA)
            ENDDO
            DO JGL=NDGUXG-NBZONG+1,NDGUXG-1
              IJA=JGL-NDGUXG+NBZONG+1
              ZEALFA(JLON,JFLD,JGL)=ZB(IIA,IJA)
            ENDDO
          ENDDO
          ZEALFA(NDLUNG:NDLUNG+NBZONL-1,JFLD,NDGUNG)=1.0_JPRB
          ZEALFA(NDLUNG:NDLUNG+NBZONL-1,JFLD,NDGUXG)=1.0_JPRB
          ZEALFA(NDLUXG-NBZONL+1:NDLUXG,JFLD,NDGUNG)=1.0_JPRB
          ZEALFA(NDLUXG-NBZONL+1:NDLUXG,JFLD,NDGUXG)=1.0_JPRB
          ZEALFA(NDLUNG,JFLD,NDGUNG:NDGUNG+NBZONG-1)=1.0_JPRB
          ZEALFA(NDLUXG,JFLD,NDGUNG:NDGUNG+NBZONG-1)=1.0_JPRB
          ZEALFA(NDLUNG,JFLD,NDGUXG-NBZONG+1:NDGUXG)=1.0_JPRB
          ZEALFA(NDLUXG,JFLD,NDGUXG-NBZONG+1:NDGUXG)=1.0_JPRB
        ENDIF

        ! Feed the extra-longitudes and latitudes
        ZEALFA(NDLUNG+NBZONL:NDLUXG-NBZONL,JFLD,NDGUNG+NBZONG:NDGUXG-NBZONG)=0.0_JPRB

      ENDIF ! ((NBZONL > 0).OR.(NBZONG > 0))

      IF (YDML_LBC%LUNBC) THEN
        ! Compute ZBU:
        DO JA=2,NBZONU
          ZA=(YDGEOMETRY%YRVETA%VETAF(JA)-YDGEOMETRY%YRVETA%VETAF(1))*ZRZONU
          IF(ZREPA(JFLD) <= -2.0_JPRB) THEN
            ZBU(JA,JFLD)=FEZBM(ZA,-ZREPA(JFLD))
          ELSE
            ZBU(JA,JFLD)=FEZBP(ZA,ZREPA(JFLD))
          ENDIF
        ENDDO
      ENDIF

    ENDIF ! IBICO(JFLD)

  ENDDO ! JFLD

  ! * C4.6: compute EALFA_GMV, EALFA_GMVS, EALFA_GFL, EALFA_TENC from ZEALFA.

  ! EALFA_GMV:
  ICPL=0
  DO JFLD=1,YDML_LBC%YYTGMVCPL%NDIM
    IROF=1
    DO JGL=1,NDGENL
      IGLG=MYLATS(JGL)
      ISTLON=NSTA(NPTRFLOFF+JGL,MYSETB)
      IENDLON=NSTA(NPTRFLOFF+JGL,MYSETB)+NONL(NPTRFLOFF+JGL,MYSETB)-1
!DEC$ IVDEP
      DO JLON=ISTLON,IENDLON
        YDML_LBC%EALFA_GMV(IROF,JFLD)=ZEALFA(JLON,ICPL+JFLD,IGLG)
        IROF=IROF+1
      ENDDO
    ENDDO
  ENDDO

  ! EALFA_GMVS:
  ICPL=YDML_LBC%YYTGMVCPL%NDIM
  DO JFLD=1,YDML_LBC%YYTGMVSCPL%NDIM
    IROF=1
    DO JGL=1,NDGENL
      IGLG=MYLATS(JGL)
      ISTLON=NSTA(NPTRFLOFF+JGL,MYSETB)
      IENDLON=NSTA(NPTRFLOFF+JGL,MYSETB)+NONL(NPTRFLOFF+JGL,MYSETB)-1
!DEC$ IVDEP
      DO JLON=ISTLON,IENDLON
        YDML_LBC%EALFA_GMVS(IROF,JFLD)=ZEALFA(JLON,ICPL+JFLD,IGLG)
        IROF=IROF+1
      ENDDO
    ENDDO
  ENDDO

  ! EALFA_GFL:
  ICPL=YDML_LBC%YYTGMVCPL%NDIM+YDML_LBC%YYTGMVSCPL%NDIM
  DO JFLD=1,YDML_LBC%NDIMCPL
    IROF=1
    DO JGL=1,NDGENL
      IGLG=MYLATS(JGL)
      ISTLON=NSTA(NPTRFLOFF+JGL,MYSETB)
      IENDLON=NSTA(NPTRFLOFF+JGL,MYSETB)+NONL(NPTRFLOFF+JGL,MYSETB)-1
!DEC$ IVDEP
      DO JLON=ISTLON,IENDLON
        YDML_LBC%EALFA_GFL(IROF,JFLD)=ZEALFA(JLON,ICPL+JFLD,IGLG)
        IROF=IROF+1
      ENDDO
    ENDDO
  ENDDO

  IF(YDML_LBC%LUNBC) THEN
  ! * C4.6.2: compute EALFAU_GMV, EALFAU_GFL from ZBU.

  ! EALFAU_GMV:
    DO JFLD=1,YDML_LBC%YYTGMVCPL%NDIM
      YDML_LBC%EALFAU_GMV(1:NFLEVG,JFLD)=0._JPRB
      IF (NBZONU > 0) YDML_LBC%EALFAU_GMV(1,JFLD)=1._JPRB

      DO JLEV=2,NBZONU
        YDML_LBC%EALFAU_GMV(JLEV,JFLD)=ZBU(NBZONU+2-JLEV,JFLD)
      ENDDO
    ENDDO

  ! EALFAU_GFL:
    ICPL=YDML_LBC%YYTGMVCPL%NDIM+YDML_LBC%YYTGMVSCPL%NDIM
    DO JFLD=1,YDML_LBC%NDIMCPL
      YDML_LBC%EALFAU_GFL(1:NFLEVG,JFLD)=0._JPRB
      IF (NBZONU > 0) YDML_LBC%EALFAU_GFL(1,JFLD)=1._JPRB

      DO JLEV=2,NBZONU
        YDML_LBC%EALFAU_GFL(JLEV,JFLD)=ZBU(NBZONU+2-JLEV,JFLD+ICPL)
      ENDDO
    ENDDO
  ENDIF

 ! EALFA_TENC:
  IF (YDML_LBC%LTENC .AND. LRPLANE) THEN

    ! EALFA_TENC for surface pressure variable:
    IROF=1
    DO JGL=1,NDGENL
      IGLG=MYLATS(JGL)
      ISTLON=NSTA(NPTRFLOFF+JGL,MYSETB)
      IENDLON=NSTA(NPTRFLOFF+JGL,MYSETB)+NONL(NPTRFLOFF+JGL,MYSETB)-1
      DO JLON=ISTLON,IENDLON
        YDML_LBC%EALFA_TENC(IROF,YDML_LBC%YYTGMVSTENC%MSP)=YDML_LBC%EALFA_GMVS(IROF,YDML_LBC%YYTGMVSCPL%MSP)
        IROF=IROF+1
      ENDDO
    ENDDO

    ! EALFA_TENC for horizontal derivatives of surface pressure variable:
    ALLOCATE(ZSPM(1,NSPEC2))
    ALLOCATE(ZGP(NGPTOT,1,3))
    ZGP(1:NGPTOT,1,1)=YDML_LBC%EALFA_TENC(1:NGPTOT,YDML_LBC%YYTGMVSTENC%MSP)
    CALL EREESPE(YDGEOMETRY,1,1,ZSPM,ZGP(1,1,1))
    CALL ESPEREE(YDGEOMETRY,1,1,ZSPM,ZGP(1,1,1),PREELL=ZGP(1,1,2),PREELM=ZGP(1,1,3))
    YDML_LBC%EALFA_TENC(1:NGPTOT,YDML_LBC%YYTGMVSTENC%MSPL)=ZGP(1:NGPTOT,1,2)
    YDML_LBC%EALFA_TENC(1:NGPTOT,YDML_LBC%YYTGMVSTENC%MSPM)=ZGP(1:NGPTOT,1,3)
    DEALLOCATE(ZGP)
    DEALLOCATE(ZSPM)
    DO JLON=1,NGPTOT
      IF (YDML_LBC%EALFA_TENC(JLON,YDML_LBC%YYTGMVSTENC%MSP) == 1.0_JPRB&
       & .OR.YDML_LBC%EALFA_TENC(JLON,YDML_LBC%YYTGMVSTENC%MSP) == 0.0_JPRB) THEN
        YDML_LBC%EALFA_TENC(JLON,YDML_LBC%YYTGMVSTENC%MSPL) = 0.0_JPRB
        YDML_LBC%EALFA_TENC(JLON,YDML_LBC%YYTGMVSTENC%MSPM) = 0.0_JPRB
      ELSEIF (YDML_LBC%EALFA_TENC(JLON,YDML_LBC%YYTGMVSTENC%MSP) > 1.0_JPRB.OR.&
       & YDML_LBC%EALFA_TENC(JLON,YDML_LBC%YYTGMVSTENC%MSP) < 0.0_JPRB) THEN
        CALL ABOR1('SUELBC_MODEL: EALFA_TENC IS OUT OF [0,1]')
      ENDIF
    ENDDO
  ENDIF

  ! * C4.7: deallocations.

  IF (ALLOCATED(ZEALP)) DEALLOCATE(ZEALP)
  IF (ALLOCATED(ZREPA)) DEALLOCATE(ZREPA)
  IF (ALLOCATED(ZEALFA)) DEALLOCATE(ZEALFA)
  IF (ALLOCATED(INEAL)) DEALLOCATE(INEAL)
  IF (ALLOCATED(INNAL)) DEALLOCATE(INNAL)
  IF (ALLOCATED(INMAL)) DEALLOCATE(INMAL)
  IF (ALLOCATED(ZB)) DEALLOCATE(ZB)
  IF (ALLOCATED(ZBU)) DEALLOCATE(ZBU)

ELSE

  YDML_LBC%EALFA_GMV(:,:)=0.0_JPRB
  YDML_LBC%EALFA_GMVS(:,:)=0.0_JPRB
  YDML_LBC%EALFA_GFL(:,:)=0.0_JPRB
  IF (YDML_LBC%LTENC) YDML_LBC%EALFA_TENC(:,:)=0.0_JPRB

ENDIF ! LECOBI

IF (ALLOCATED(IBICO)) DEALLOCATE(IBICO)

! * C5: Other variables for grid-point coupling:
!       Allocation and computation of GMGT3, EALFAGT3GMV, EALFAGT3GMVS, EALFAGT3GFL.

ALLOCATE(YDML_LBC%GMGT3(NPROMA,YDGEOMETRY%YRELBC_GEO%NCPLBLKS))
IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'GMGT3 ',&
 & SIZE(YDML_LBC%GMGT3 ),SHAPE(YDML_LBC%GMGT3 )

ALLOCATE(YDML_LBC%EALFAGT3GMV(NPROMA,YDML_LBC%YYTGMVCPL%NDIM,YDGEOMETRY%YRELBC_GEO%NCPLBLKS))
IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A12,' ALLOCATED ',8I8)") 'EALFAGT3GMV ',&
 & SIZE(YDML_LBC%EALFAGT3GMV),SHAPE(YDML_LBC%EALFAGT3GMV)

ALLOCATE(YDML_LBC%EALFAGT3GMVS(NPROMA,YDML_LBC%YYTGMVSCPL%NDIM,YDGEOMETRY%YRELBC_GEO%NCPLBLKS))
IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A12,' ALLOCATED ',8I8)") 'EALFAGT3GMVS',&
 & SIZE(YDML_LBC%EALFAGT3GMVS),SHAPE(YDML_LBC%EALFAGT3GMVS)

ALLOCATE(YDML_LBC%EALFAGT3GFL(NPROMA,YDML_LBC%NDIMCPL,YDGEOMETRY%YRELBC_GEO%NCPLBLKS))
IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A12,' ALLOCATED ',8I8)") 'EALFAGT3GFL ',&
 & SIZE(YDML_LBC%EALFAGT3GFL),SHAPE(YDML_LBC%EALFAGT3GFL)

IF(YDML_LBC%LUNBC) THEN
  ALLOCATE(YDML_LBC%GMGT4(NGPTOT))
  IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'GMGT4 ',&
   & SIZE(YDML_LBC%GMGT4 ),SHAPE(YDML_LBC%GMGT4 )
ENDIF

DO JKGLO=1,NGPTOT,NPROMA
  IBL=(JKGLO-1)/NPROMA+1
  IBLOCKCPL=YDGEOMETRY%YRELBC_GEO%MPTRCPLBLK(IBL)
  IF (IBLOCKCPL > 0) THEN ! This block contains coupling points
    DO JROF=1,YDGEOMETRY%YRELBC_GEO%NIND_LEN(IBLOCKCPL)
      IGPTOT=JKGLO+YDGEOMETRY%YRELBC_GEO%NIND_LIST(JROF,IBLOCKCPL)-1
      YDML_LBC%EALFAGT3GMV(JROF,:,IBLOCKCPL)=YDML_LBC%EALFA_GMV(IGPTOT,:)
      YDML_LBC%EALFAGT3GMVS(JROF,:,IBLOCKCPL)=YDML_LBC%EALFA_GMVS(IGPTOT,:)
      YDML_LBC%EALFAGT3GFL(JROF,:,IBLOCKCPL)=YDML_LBC%EALFA_GFL(IGPTOT,:)
      YDML_LBC%GMGT3(JROF,IBLOCKCPL)=YDGSGEOM_NB%GM(IGPTOT)
    ENDDO
  ENDIF
  IF (YDML_LBC%LUNBC) THEN
    DO JROF=1,MIN(NPROMA,NGPTOT-JKGLO+1)
      IGPTOT=JKGLO+JROF-1
      YDML_LBC%GMGT4(IGPTOT)=YDGSGEOM_NB%GM(IGPTOT)
    ENDDO
  ENDIF
ENDDO

! * C6: Other variables for spectral nudging

IF (YDML_LBC%LESPCPL) THEN
  YDML_LBC%RNUDTFRAC=SIGN(1._JPRB,YDML_GCONF%YRRIP%TSTEP)/REAL(YDML_LBC%NEFRSPCPL,JPRB)
  YDML_LBC%LSPNUSPDL=(SPNUDSP > ZEPS2).AND.(MYSETV==NBSETSP)
ELSE
  YDML_LBC%RNUDTFRAC=0._JPRB
  YDML_LBC%LSPNUSPDL=.FALSE.
ENDIF

IF (YDML_LBC%LESPCPL) THEN
  ALLOCATE(YDML_LBC%LNUDSPGFL(MAX(1,YDML_GCONF%YGFL%NUMSPFLDS)))
  IF (LALLOPR) WRITE(NULOUT,"(1X,'ARRAY ',A10,' ALLOCATED ',8I8)") 'LNUDSPGFL',&
 & SIZE(YDML_LBC%LNUDSPGFL),SHAPE(YDML_LBC%LNUDSPGFL)
  YDML_LBC%LNUDSPGFL(:)=.FALSE.
  IF ((SPNUDQ>ZEPS2).AND.YDML_GCONF%YGFL%YQ%MPSP>0.AND.YDML_GCONF%YGFL%YQ_NL%LSP) THEN
    YDML_LBC%LNUDSPGFL(YDML_GCONF%YGFL%YQ%MPSP)=.TRUE.
  ENDIF
ENDIF

!--------------------------------------------------------------------------

!      Part D: printings.

WRITE(NULOUT,*) ' --- PRINTINGS IN SUELBC_MODEL --- '

! * D1: Control frequency of LBC.
WRITE(UNIT=NULOUT,FMT='('' Frequency of LBC: '')')
WRITE(UNIT=NULOUT,FMT='(''  TEFRCL = '',F10.1,'' NEFRCL = '',I15)') YDML_LBC%TEFRCL,YDML_LBC%NEFRCL

! * D2: LECOBI.
WRITE(NULOUT,'(''  LECOBI = '',L2)') YDML_LBC%LECOBI

! * D3: Relaxation coefficients EALFA_GMV, EALFA_GMVS, EALFA_GFL, EALFA_TENC.
IF (YDML_LBC%LECOBI) THEN
  IF (LOUTPUT) THEN
    WRITE(UNIT=NULOUT,FMT='('' Relaxation coefficients: '')')
    DO JFLD=1,YDML_LBC%YYTGMVCPL%NDIM
      WRITE(UNIT=NULOUT,FMT=' (''  JFLD = '',I3)') JFLD
      WRITE(UNIT=NULOUT,FMT=' (''  EALFA_GMV FOR JFLD'')')
      DO JGL=1,NDGENL,100
        IGLG=MYLATS(JGL)
        WRITE(UNIT=NULOUT,FMT='(2X,14(1X,E8.2))') (YDML_LBC%EALFA_GMV(JLON+NSTAGP(JGL),JFLD),JLON=1,NDLON,100)
      ENDDO
    ENDDO
    DO JFLD=1,YDML_LBC%YYTGMVSCPL%NDIM
      WRITE(UNIT=NULOUT,FMT=' (''  JFLD = '',I3)') JFLD
      WRITE(UNIT=NULOUT,FMT=' (''  EALFA_GMVS FOR JFLD'')')
      DO JGL=1,NDGENL,100
        IGLG=MYLATS(JGL)
        WRITE(UNIT=NULOUT,FMT='(2X,14(1X,E8.2))') (YDML_LBC%EALFA_GMVS(JLON+NSTAGP(JGL),JFLD),JLON=1,NDLON,100)
      ENDDO
    ENDDO
    DO JFLD=1,YDML_LBC%NDIMCPL
      WRITE(UNIT=NULOUT,FMT=' (''  JFLD = '',I3)') JFLD
      WRITE(UNIT=NULOUT,FMT=' (''  EALFA_GFL FOR JFLD'')')
      DO JGL=1,NDGENL,100
        IGLG=MYLATS(JGL)
        WRITE(UNIT=NULOUT,FMT='(2X,14(1X,E8.2))') (YDML_LBC%EALFA_GFL(JLON+NSTAGP(JGL),JFLD),JLON=1,NDLON,100)
      ENDDO
    ENDDO
    IF(YDML_LBC%LUNBC) THEN
      WRITE(UNIT=NULOUT,FMT='('' Upper Relaxation coefficients: '')')
      DO JFLD=1,YDML_LBC%YYTGMVCPL%NDIM
        WRITE(UNIT=NULOUT,FMT=' (''  JFLD = '',I3)') JFLD
        WRITE(UNIT=NULOUT,FMT=' (''  EALFAU_GMV FOR JFLD'')')
        WRITE(UNIT=NULOUT,FMT='(2X,14(1X,E8.2))')&
          &YDML_LBC%EALFAU_GMV(1:NBZONU+1,JFLD)
      ENDDO
    ENDIF
    IF (YDML_LBC%LTENC.AND.LRPLANE) THEN
      DO JFLD=1,YDML_LBC%YYTGMVSTENC%NDIM
        WRITE(UNIT=NULOUT,FMT=' (''  JFLD = '',I3)') JFLD
        WRITE(UNIT=NULOUT,FMT=' (''  EALFA_TENC FOR JFLD'')')
        DO JGL=1,NDGENL,100
          IGLG=MYLATS(JGL)
          WRITE(UNIT=NULOUT,FMT='(2X,14(1X,E8.2))') (YDML_LBC%EALFA_TENC(JLON+NSTAGP(JGL),JFLD),JLON=1,NDLON,100)
        ENDDO
      ENDDO
    ENDIF
  ENDIF
ENDIF

! * D4: Other variables for spectral nudging.
IF (YDML_LBC%LESPCPL) THEN
  WRITE(UNIT=NULOUT,FMT='('' Other variables for spectral nudging: '')')
  WRITE(UNIT=NULOUT,FMT='(''  RNUDTFRAC = '',E20.14)') YDML_LBC%RNUDTFRAC
  WRITE(UNIT=NULOUT,FMT='(''  LSPNUSPDL = '',L2)') YDML_LBC%LSPNUSPDL
  WRITE(UNIT=NULOUT,FMT='(''  LNUDSPGFL = '',20(1X,L2))') YDML_LBC%LNUDSPGFL(:)
  WRITE(UNIT=NULOUT,FMT='(''  NEFRSPCPL = '',I5)') YDML_LBC%NEFRSPCPL
  WRITE(UNIT=NULOUT,FMT='(''   NEK0     = '',I5)') YDML_LBC%NEK0
  WRITE(UNIT=NULOUT,FMT='(''   NEK1     = '',I5)') YDML_LBC%NEK1
  WRITE(UNIT=NULOUT,FMT='(''   NEN1     = '',I5)') YDML_LBC%NEN1
  WRITE(UNIT=NULOUT,FMT='(''   NEN2     = '',I5)') YDML_LBC%NEN2
  WRITE(UNIT=NULOUT,FMT='(''  RNUTENC   = '',E20.14)') YDML_LBC%RNUTENC
  WRITE(UNIT=NULOUT,FMT='(''   RTENC    = '',E20.14)') YDML_LBC%RTENC
  WRITE(NULOUT, FMT='(&
  & '' SPNUDVOR = '',E11.4,'' SPNUDDIV = '',E11.4,'' SPNUDT = '',E11.4,&
  & '' SPNUDQ = '',E11.4,'' SPNUDSP = '',E11.4&
  & )') YDML_LBC%SPNUDVOR,YDML_LBC%SPNUDDIV,YDML_LBC%SPNUDT,YDML_LBC%SPNUDQ,YDML_LBC%SPNUDSP
ENDIF

WRITE(NULOUT,*) ' '

!--------------------------------------------------------------------------
END ASSOCIATE
IF (LHOOK) CALL DR_HOOK('YEMLBC_MODEL:SUELBC_MODEL',1,ZHOOK_HANDLE)
END SUBROUTINE SUELBC_MODEL

SUBROUTINE DEALLOCATE_ELBC0B(YDML_LBC)

!--------------------------------------------------------------------------
! deallocates 'ELBC0B' arrays
!--------------------------------------------------------------------------

IMPLICIT NONE

TYPE(TELBC_MODEL) ,INTENT(INOUT) :: YDML_LBC 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('YEMLBC_MODEL:DEALLOCATE_ELBC0B',0,ZHOOK_HANDLE)
!     ------------------------------------------------------------------

IF (ALLOCATED(YDML_LBC%EALFA_GMV)) DEALLOCATE(YDML_LBC%EALFA_GMV)
IF (ALLOCATED(YDML_LBC%EALFA_GMVS)) DEALLOCATE(YDML_LBC%EALFA_GMVS)
IF (ALLOCATED(YDML_LBC%EALFA_GFL)) DEALLOCATE(YDML_LBC%EALFA_GFL)
IF (ALLOCATED(YDML_LBC%EALFA_TENC)) DEALLOCATE(YDML_LBC%EALFA_TENC)
IF (ALLOCATED(YDML_LBC%EALFAGT3GMV)) DEALLOCATE(YDML_LBC%EALFAGT3GMV)
IF (ALLOCATED(YDML_LBC%EALFAGT3GMVS)) DEALLOCATE(YDML_LBC%EALFAGT3GMVS)
IF (ALLOCATED(YDML_LBC%EALFAGT3GFL)) DEALLOCATE(YDML_LBC%EALFAGT3GFL)
IF (ALLOCATED(YDML_LBC%EALFAU_GMV)) DEALLOCATE(YDML_LBC%EALFAU_GMV)
IF (ALLOCATED(YDML_LBC%EALFAU_GFL)) DEALLOCATE(YDML_LBC%EALFAU_GFL)

IF (ALLOCATED(YDML_LBC%GMGT3   )) DEALLOCATE(YDML_LBC%GMGT3)
IF (ALLOCATED(YDML_LBC%GMGT4   )) DEALLOCATE(YDML_LBC%GMGT4)
IF (ALLOCATED(YDML_LBC%EWB     )) DEALLOCATE(YDML_LBC%EWB)
IF (ALLOCATED(YDML_LBC%EWBDFIBW)) DEALLOCATE(YDML_LBC%EWBDFIBW)
IF (ALLOCATED(YDML_LBC%EWBDFIFW)) DEALLOCATE(YDML_LBC%EWBDFIFW)
IF (ALLOCATED(YDML_LBC%LNUDSPGFL)) DEALLOCATE(YDML_LBC%LNUDSPGFL)

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('YEMLBC_MODEL:DEALLOCATE_ELBC0B',1,ZHOOK_HANDLE)
END SUBROUTINE DEALLOCATE_ELBC0B

!=============================================================================

END MODULE YEMLBC_MODEL
