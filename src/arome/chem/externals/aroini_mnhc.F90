!!    ##############################
      SUBROUTINE AROINI_MNHC(OUSECHEM, OORILAM, ODUST, ODEPOS,&
                             OINITCHEM, OINITDUST, OINITORILAM,&
                             KDAY, KMONTH, KYEAR,&
                             KLUOUT, KPROC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ##############################
!!
!!*** *AROINI_MNHC*
!!
!!    PURPOSE
!!    -------
!        initialize chemical variables
!        initialize chemical core system
!!
!!
!!    AUTHOR
!!    ------
!!    P. Tulet  *CNRM / GMEI* and contributors of MesoNH-C (K. Shure, C. Mari, V. Crassier)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 02/04/05
!!
!!    EXTERNAL
!!    --------
USE MODI_CH_INIT_SCHEME
USE MODI_CH_INIT_CCS
USE MODI_CH_INIT_JVALUES
USE MODI_CH_AER_INIT_SOA
USE MODI_CH_AER_MOD_INIT
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CONF,        ONLY: CCONF
USE MODD_CH_M9
USE MODD_CH_MNHC_n 
USE MODD_CH_SOLVER_n
USE MODD_CH_AERO_n
USE MODD_CH_AEROSOL 
USE MODD_DUST
USE MODE_DUSTOPT

!!
!!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
!
LOGICAL, INTENT(IN)  :: OUSECHEM
LOGICAL, INTENT(IN)  :: ODUST
LOGICAL, INTENT(IN)  :: OORILAM
LOGICAL, INTENT(IN)  :: OINITCHEM
LOGICAL, INTENT(IN)  :: OINITDUST
LOGICAL, INTENT(IN)  :: OINITORILAM
LOGICAL, INTENT(IN)  :: ODEPOS
INTEGER, INTENT(IN)  :: KLUOUT   ! output listing channel
INTEGER, INTENT(IN)  :: KYEAR    ! year of initialization
INTEGER, INTENT(IN)  :: KMONTH   ! month of initialization
INTEGER, INTENT(IN)  :: KDAY     ! day of initialization
INTEGER, INTENT(IN)  :: KPROC    ! proc number

!
!*      0.2    pointer to model 1
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AROINI_MNHC',0,ZHOOK_HANDLE)
CALL CH_MNHC_GOTO_MODEL(1, 1)
CALL CH_SOLVER_GOTO_MODEL(1, 1)
CALL CH_AERO_GOTO_MODEL(1, 1)
!-------------------------------------------------------------------------------
!
!
!*       0.    PROLOGUE (default chemical namelist )
!              --------
! 0.1 MODD_CH_MNHC_n

LUSECHEM            = OUSECHEM
LCH_INIT_FIELD      = OINITCHEM
LCH_SURFACE_FLUX    = .TRUE.
LCH_CONV_SCAV       = .FALSE.
CCH_EXPLICIT_SCAV   = 'NONE'
CCHEM_INPUT_FILE    = 'chimie.nam'
CCH_TDISCRETIZATION = 'SPLIT' ! not use in AROME
NCH_SUBSTEPS        = 5
LCH_TUV_ONLINE      = .FALSE.
CCH_TUV_LOOKUP      = 'PHOTO.TUV39'
CCH_TUV_CLOUDS      = 'NONE'
XCH_TUV_ALBNEW      = 0.02
XCH_TUV_DOBNEW      = 320.
XCH_TUV_TUPDATE     = 1200. 
CCH_VEC_METHOD      = 'HOR'  
NCH_VEC_LENGTH      = 1000   ! not use in AROME
XCH_TS1D_TSTEP      = 600.   ! not use in AROME
CCH_TS1D_COMMENT    = 'no comment' ! not use in AROME
CCH_TS1D_FILENAME   = 'IO1D' ! not use in AROME
CCH_SCHEME          = 'ReLACS'

! 0.2 MODD_CH_SOLVER_n

CSOLVER = 'EXQSSA' ! default values of each solver can be found/modified in MODD_CH_SOLVER_n
XSLOW    = 100.
XFAST    = 0.1
XDTMIN   = 10.
XDTMAX   = 10.
XDTFIRST = 10.
XRTOL    = 0.001
XATOL    = 1E2 

! 0.3 MODD_CH_AEROSOL

LORILAM   = OORILAM
CRGUNIT   = 'NUMB'    ! type of log-normal geometric mean radius given
!                                     ! in nameliste (mass on number)

LVARSIGI  = .FALSE.   ! switch to active pronostic dispersion for I mode
LVARSIGJ  = .FALSE.   ! switch to active pronostic dispersion for J mode
LHETEROSO4 = .FALSE.  ! switch to active sulfates heteronegeous production
LSEDIMAERO = .FALSE.  ! switch to active aerosol sedimentation
LAERINIT   = OINITORILAM  ! switch to active thermodynamics
                      ! mineral equilibrium for each mode
XN0IMIN    = 10.       ! minimum particule number value for I mode
XN0JMIN    = 1.        ! minimum particule number value for J mode
XINIRADIUSI= 0.05      ! mean radius initialization for I mode (um)
XINIRADIUSJ= 0.2       ! mean radius initialization for J mode (um)
XINISIGI   = 1.8       ! dispersion initialization for I mode 
XINISIGJ   = 2.0       ! dispersion initialization for J mode
XSIGIMIN   = 1.20      ! minimum dispersion value for I mode
XSIGJMIN   = 1.20      ! minimum dispersion value for J mode
XSIGIMAX   = 3.60      ! maximum dispersion value for I mode
XSIGJMAX   = 3.60      ! maximum dispersion value for J mode
XCOEFRADIMAX  = 10.    ! maximum increasement for Rg mode I
XCOEFRADIMIN  = 0.1    ! maximum decreasement for Rg mode I
XCOEFRADJMAX  = 10.    ! maximum increasement for Rg mode J
XCOEFRADJMIN  = 0.1    ! maximum decreasement for Rg mode J
CMINERAL      = "EQSAM" ! mineral equilibrium scheme
CORGANIC      = "NONE" ! organic equilibrium scheme
CNUCLEATION   = "NONE" ! sulfates nucleation scheme

! 0.4 MODD_DUST

LDUST      = ODUST
LDSTINIT   = OINITDUST
LVARSIG    = .FALSE.   ! switch to active pronostic dispersion for all modes
LSEDIMDUST = .FALSE.   ! switch to active aerosol sedimentation
XSIGMIN   = 1.20      ! minimum dispersion value for dust mode
XSIGMAX   = 3.60      ! maximum dispersion value for dust mode
XCOEFRADMAX  = 10.    ! maximum increasement for Rg mode dust
XCOEFRADMIN  = 0.1    ! maximum decreasement for Rg mode dust

LDEPOS_DST(:) = ODEPOS

IF (LUSECHEM) THEN
!*       1.1    INITIALIZE CHEMICAL CORE SYSTEM
!              -------------------------------

!
  CALL CH_INIT_SCHEME(KLUOUT)
  CALL CH_INIT_CCS(NEQ, KLUOUT, 5)
  CCH_SCHEME = "NONE"

!
!
!
!-------------------------------------------------------------------------------
!
!*       1.2    INITIALIZE JVALUES
!              ------------------
 CALL CH_INIT_JVALUES(KDAY, KMONTH, KYEAR, KLUOUT)
!
END IF

!-------------------------------------------------------------------------------
!
!*       2.    INITIALIZE SOA
!              ---------------
IF ((LUSECHEM).AND.(LORILAM)) THEN
 CALL CH_AER_INIT_SOA(KLUOUT, 5)
 CALL CH_AER_MOD_INIT
END IF


!*       3.     Read in look up tables of dust optical properties
!                -------------------------------------------------
IF (LDUST) THEN
  CALL DUST_OPT_LKT_SET1()
  CALL DUST_OPT_LKT_SET2()
  CALL DUST_OPT_LKT_SET3()
  CALL DUST_OPT_LKT_SET4()
  CALL DUST_OPT_LKT_SET5()
  CALL DUST_OPT_LKT_SET6()
  CALL DUST_OPT_LKT_SET7()
  CALL DUST_OPT_LKT_SET8()
  CALL DUST_OPT_LKT_SET9()
  CALL DUST_OPT_LKT_SET10()
ENDIF

!-------------------------------------------------------------------------------
!

IF (LHOOK) CALL DR_HOOK('AROINI_MNHC',1,ZHOOK_HANDLE)
END SUBROUTINE AROINI_MNHC
