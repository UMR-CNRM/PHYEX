!     ######spl
SUBROUTINE AROINI_BUDGET(LDBU_ENABLE, HLUOUT,KULOUT, PTSTEP,KSV,KRR, &
                         HRAD,HDCONV,HSCONV,HTURB,HTURBDIM,HCLOUD,   &
                         HMET_ADV_SCHEME, HSV_ADV_SCHEME)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
!!!!old definition SUBROUTINE AROINI_BUDGET 
!**** *AROINI_BUDGET*   - Initialize common meso_NH MODD_ used in BUDGET for AROME

!     Purpose.
!     --------
!            Set implicit default values for MODD_BUDGET for the use in AROME 

!**   Interface.
!     ----------
!        *CALL* *AROINI_BUDGET

!        Explicit arguments :
!        --------------------
!       None

!        Implicit arguments :
!        --------------------
!       None

!     Method.
!     -------
!        See documentation
!        To use budgets in DDH for AROME, budget must have type CART. First dimension is NPROMA and
!        second dimension is 1. Budgets are reset after each tipe step. Processes not used in AROME are
!        marked with 3

!     Externals.
!     ----------

!     Reference.
!     ----------
!
!     Author.
!     -------
!        Y. Seity 

!     Modifications.
!     --------------
!        Original :    03-12-12
!        T. Kovacic    05-04-27   Initialization for DDH
!     ------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------ 
!
USE MODD_BUDGET  !!!!, ONLY : NBUMOD, LBUDGET_U, LBUDGET_V, LBUDGET_W, LBUDGET_SV&
!                !!!!     &, LBUDGET_RV, LBUDGET_RC, LBUDGET_RI, LBUDGET_TH
USE MODI_INI_BUDGET
!
IMPLICIT NONE
!
!*       0.1   declarations of argument
!
!
LOGICAL :: LDBU_ENABLE
CHARACTER (LEN=*), INTENT(IN) :: HLUOUT ! name of output listing
INTEGER, INTENT(IN)   :: KULOUT
REAL, INTENT(IN) :: PTSTEP              ! time step
INTEGER, INTENT(IN) :: KSV              ! number of scalar variables
INTEGER, INTENT(IN) :: KRR              ! number of moist variables
CHARACTER (LEN=*), INTENT(IN) :: HRAD   ! type of the radiation scheme
CHARACTER (LEN=*), INTENT(IN) :: HDCONV ! type of the deep convection scheme
CHARACTER (LEN=*), INTENT(IN) :: HSCONV ! type of the shallow convection scheme
CHARACTER (LEN=*), INTENT(IN) :: HTURB  ! type of the turbulence scheme


CHARACTER (LEN=*), INTENT(IN) :: HTURBDIM! dimensionnality of the turbulence 
                                        ! scheme
CHARACTER (LEN=*), INTENT(IN) :: HCLOUD ! type of microphysical scheme
CHARACTER (LEN=*), INTENT(IN) :: HMET_ADV_SCHEME ! type of advection scheme for meteorological scalar variables
CHARACTER (LEN=*), INTENT(IN) :: HSV_ADV_SCHEME ! type of advection scheme for tracer scalar variables
!
!*       0.2   Declarations of local variables :
!
!
LOGICAL  :: GNUMDIFU        ! switch to activate momentum numerical diffusion
LOGICAL  :: GNUMDIFTH       ! for theta
LOGICAL  :: GNUMDIFSV       ! for scalars
LOGICAL  :: GHORELAX_UVWTH  ! switch for the horizontal relaxation for U,V,W,TH
LOGICAL  :: GHORELAX_RV     ! switch for the horizontal relaxation for Rv
LOGICAL  :: GHORELAX_RC     ! switch for the horizontal relaxation for Rc
LOGICAL  :: GHORELAX_RR     ! switch for the horizontal relaxation for Rr
LOGICAL  :: GHORELAX_RI     ! switch for the horizontal relaxation for Ri
LOGICAL  :: GHORELAX_RS     ! switch for the horizontal relaxation for Rs
LOGICAL  :: GHORELAX_RG     ! switch for the horizontal relaxation for Rg
LOGICAL  :: GHORELAX_RH     ! switch for the horizontal relaxation for Rh
LOGICAL  :: GHORELAX_TKE    ! switch for the horizontal relaxation for tke
LOGICAL,DIMENSION(:), ALLOCATABLE :: GHORELAX_SV     ! switch for the
                              ! horizontal relaxation for scalar variables
LOGICAL  :: GVE_RELAX        ! switch to activate the vertical 
                                        ! relaxation
LOGICAL  :: GCHTRANS         ! switch to activate convective 
                                        !transport for SV
LOGICAL  :: GDRAGTREE        ! switch to activate vegetation drag
LOGICAL  :: LLNUDGING
REAL  :: ZTSTEP 
!
!-------------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('AROINI_BUDGET',0,ZHOOK_HANDLE)
LBU_ENABLE = LDBU_ENABLE
IF (.NOT.LBU_ENABLE) THEN
!
!*       1.    RUN WITHOUT BUDGETS, early versions of AROME
!              --------------------------------------------
!
  NBUMOD=5       !pour qu'il ne fasse pas de budget dans rain_ice
  !pour qu'il ne fasse pas de budget dans Ice_adjust
  LBUDGET_U =.FALSE.
  LBUDGET_V =.FALSE.
  LBUDGET_W =.FALSE.
  LBUDGET_TH=.FALSE.
  LBUDGET_TKE=.FALSE.
  LBUDGET_RV=.FALSE.   
  LBUDGET_RC=.FALSE.
  
  LBUDGET_RR =.FALSE.
  LBUDGET_RI =.FALSE.
  LBUDGET_RS =.FALSE.
  LBUDGET_RG =.FALSE.
  LBUDGET_RH =.FALSE.
  LBUDGET_SV=.FALSE.

ELSE
!
!-------------------------------------------------------------------------------
!
!*       2.    SETUP 
!              -----
!
!
!*       2.1  SETUP FOR DDH
!
!
!-- For DDH in AROME, variables from MODD_BUDGERT are not read from namelists.
!   Proper values are given in subroutine DDH_BUDTET_SETUP
!
!  CALL DDH_BUDTET_SETUP
!
!*       3.  INITIALIZATION FOR AROME
!
!
!
!*       3.1.  Setting switches for not needed processes to FALS
!
!
  GNUMDIFU        = .FALSE.  
  GNUMDIFTH       = .FALSE.
  GNUMDIFSV       = .FALSE.
  GHORELAX_UVWTH  = .FALSE.  
  GHORELAX_RV     = .FALSE.  
  GHORELAX_RC     = .FALSE.    
  GHORELAX_RR     = .FALSE.  
  GHORELAX_RI     = .FALSE.  
  GHORELAX_RS     = .FALSE.  
  GHORELAX_RG     = .FALSE.  
  GHORELAX_RH     = .FALSE.  
  GHORELAX_TKE    = .FALSE.  
  GVE_RELAX       = .FALSE.  
  GCHTRANS        = .FALSE.  
  ALLOCATE(GHORELAX_SV(1))
  GHORELAX_SV(1)  = .FALSE.  
  GDRAGTREE       = .FALSE.

!
!*       3.1.  Initialization
!
!
LLNUDGING=.FALSE.
ZTSTEP = PTSTEP
  CALL INI_BUDGET(KULOUT, HLUOUT, ZTSTEP,KSV,KRR,                      &
      &GNUMDIFU,GNUMDIFTH,GNUMDIFSV,                                   &
      &GHORELAX_UVWTH,GHORELAX_RV,GHORELAX_RC,GHORELAX_RR,             &
      &GHORELAX_RI,GHORELAX_RS, GHORELAX_RG, GHORELAX_RH,GHORELAX_TKE, &
      &GHORELAX_SV,GVE_RELAX,GCHTRANS,LLNUDGING,GDRAGTREE,             &
      &HRAD,HSCONV,HDCONV,HTURB,HTURBDIM,HCLOUD,                       &
      &HMET_ADV_SCHEME,HSV_ADV_SCHEME)

!  LBU_ENABLE=.FALSE.


ENDIF

IF (LHOOK) CALL DR_HOOK('AROINI_BUDGET',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE AROINI_BUDGET
