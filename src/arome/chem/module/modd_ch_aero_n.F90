!     ######spl
       MODULE MODD_CH_AERO_n
       USE PARKIND1, ONLY : JPRB
       USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!     ######################
!!
!!     PURPOSE
!!     -------
!!
!!     declaration of variables and types for the aerosol system
!!
!!     METHOD
!!     ------
!!
!!
!!     REFERENCE
!!     ---------
!!     none
!!
!!
!!     AUTHOR
!!     ------
!!     P. Tulet (LA)
!!
!!
!!     MODIFICATIONS
!!     -------------
!!     (30-01-01) P.Tulet (LA) * modifications for secondary biogenics aerosols
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
USE MODD_PARAMETERS, ONLY: JPMODELMAX
IMPLICIT NONE

TYPE CH_AERO_t
!
!* normalisation parameters
!
  REAL, DIMENSION(:,:), POINTER :: XN0=>NULL()     ! Number concentration
  REAL, DIMENSION(:,:), POINTER :: XRG0=>NULL()    ! Geometric mean size
  REAL, DIMENSION(:,:), POINTER :: XSIG0=>NULL()   ! Dispersion ln(sigma)
  REAL, DIMENSION(:,:,:,:), POINTER   :: XN3D=>NULL()      ! Number concentration
  REAL, DIMENSION(:,:,:,:), POINTER   :: XRG3D=>NULL()     ! Geometric mean size
  REAL, DIMENSION(:,:,:,:), POINTER   :: XSIG3D=>NULL()    ! dispersion (sigma)
  REAL, DIMENSION(:,:,:,:), POINTER   :: XM3D=>NULL()      ! moments
  REAL, DIMENSION(:,:,:,:), POINTER   :: XSEDA=>NULL()     ! sedimentation
  REAL, DIMENSION(:,:,:),   POINTER   :: XVDEPAERO=>NULL() ! aerosol dry deposition
  REAL, DIMENSION(:,:,:,:,:), POINTER :: XCTOTA3D=>NULL()  ! Total concentration of species
!
  REAL, DIMENSION(:,:,:), POINTER :: XFTEST=>NULL()
  REAL, DIMENSION(:,:,:), POINTER :: XCTOTA=>NULL() ! Total concentration of species
                                                    ! (HNO3, ! H2SO4, NH3) present in
                                                    ! each of the aerosol mode (ug/m3)
  REAL, DIMENSION(:,:,:), POINTER :: XCCTOT=>NULL() ! Composition of 3rd Moment (%)
  REAL, DIMENSION(:,:),   POINTER :: XCTOTG=>NULL() ! Total concentration of volatile
                                                    ! species (HNO3, NH3) (ug/m3) in
                                                    ! the air
  REAL, DIMENSION(:,:,:,:), POINTER :: XFRAC=>NULL() ! Gas fraction into organic species
  REAL, DIMENSION(:,:,:,:), POINTER :: XMI=>NULL()   ! Molar mass of aerosol species (g/mol)
  REAL, DIMENSION(:,:,:,:), POINTER :: XSOLORG=>NULL() ! Solubility fraction of SOA (%)
  REAL, DIMENSION(:,:), POINTER :: XRHOP0=>NULL()    ! Condensed phase density (kg/m3)
  REAL, DIMENSION(:,:,:,:), POINTER :: XRHOP3D=>NULL()    ! Condensed phase density (kg/m3)
  REAL, DIMENSION(:),   POINTER :: XLAMBDA=>NULL()  ! Mean free path of background
                                              ! gas molecules
  REAL, DIMENSION(:),   POINTER :: XMU=>NULL()      ! Gas viscosity (kg/(ms))
!
!--------------------------------------------------------------------------
!
!* Growth parameters
!
  REAL, DIMENSION(:,:), POINTER :: XOM=>NULL()

!
!----------------------------------------------------------------------------
!
!* Nucleation/cond. growth parameters
!
  REAL, DIMENSION(:), POINTER :: XSO4RAT=>NULL()
                                ! Rate of gas phase production of
                                ! sulfuric acid (molec./cm3/s)
!
!----------------------------------------------------------------------------
!
  LOGICAL       :: GSEDFIX = .TRUE. ! flag used in CH_AER_SEDIM_n routine
!
END TYPE CH_AERO_t

TYPE(CH_AERO_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: CH_AERO_MODEL

REAL, DIMENSION(:,:), POINTER :: XN0=>NULL()
REAL, DIMENSION(:,:), POINTER :: XRG0=>NULL()
REAL, DIMENSION(:,:), POINTER :: XSIG0=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER   :: XN3D=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER   :: XRG3D=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER   :: XSIG3D=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER   :: XM3D=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER   :: XSEDA=>NULL()
REAL, DIMENSION(:,:,:),   POINTER   :: XVDEPAERO=>NULL()
REAL, DIMENSION(:,:,:,:,:), POINTER :: XCTOTA3D=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XFTEST=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCTOTA=>NULL()
REAL, DIMENSION(:,:,:), POINTER :: XCCTOT=>NULL()
REAL, DIMENSION(:,:),   POINTER :: XCTOTG=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XFRAC=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XMI=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XSOLORG=>NULL()
REAL, DIMENSION(:,:), POINTER :: XRHOP0=>NULL()
REAL, DIMENSION(:,:,:,:), POINTER :: XRHOP3D=>NULL()
REAL, DIMENSION(:),   POINTER :: XLAMBDA=>NULL()
REAL, DIMENSION(:),   POINTER :: XMU=>NULL()
REAL, DIMENSION(:,:), POINTER :: XOM=>NULL()
REAL, DIMENSION(:), POINTER :: XSO4RAT=>NULL()
LOGICAL,            POINTER :: GSEDFIX=>NULL()

CONTAINS

SUBROUTINE CH_AERO_GOTO_MODEL(KFROM, KTO)
INTEGER, INTENT(IN) :: KFROM, KTO
!
! Save current state for allocated arrays
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MODD_CH_AERO_N:CH_AERO_GOTO_MODEL',0,ZHOOK_HANDLE)
CH_AERO_MODEL(KFROM)%XN0=>XN0
CH_AERO_MODEL(KFROM)%XRG0=>XRG0
CH_AERO_MODEL(KFROM)%XSIG0=>XSIG0
CH_AERO_MODEL(KFROM)%XN3D=>XN3D
CH_AERO_MODEL(KFROM)%XRG3D=>XRG3D
CH_AERO_MODEL(KFROM)%XSIG3D=>XSIG3D
CH_AERO_MODEL(KFROM)%XM3D=>XM3D
CH_AERO_MODEL(KFROM)%XSEDA=>XSEDA
CH_AERO_MODEL(KFROM)%XVDEPAERO=>XVDEPAERO
CH_AERO_MODEL(KFROM)%XCTOTA3D=>XCTOTA3D
CH_AERO_MODEL(KFROM)%XFTEST=>XFTEST
CH_AERO_MODEL(KFROM)%XCTOTA=>XCTOTA
CH_AERO_MODEL(KFROM)%XCCTOT=>XCCTOT
CH_AERO_MODEL(KFROM)%XCTOTG=>XCTOTG
CH_AERO_MODEL(KFROM)%XFRAC=>XFRAC
CH_AERO_MODEL(KFROM)%XMI=>XMI
CH_AERO_MODEL(KFROM)%XSOLORG=>XSOLORG
CH_AERO_MODEL(KFROM)%XRHOP0=>XRHOP0
CH_AERO_MODEL(KFROM)%XRHOP3D=>XRHOP3D
CH_AERO_MODEL(KFROM)%XLAMBDA=>XLAMBDA
CH_AERO_MODEL(KFROM)%XMU=>XMU
CH_AERO_MODEL(KFROM)%XOM=>XOM
CH_AERO_MODEL(KFROM)%XSO4RAT=>XSO4RAT
!
! Current model is set to model KTO
XN0=>CH_AERO_MODEL(KTO)%XN0
XRG0=>CH_AERO_MODEL(KTO)%XRG0
XSIG0=>CH_AERO_MODEL(KTO)%XSIG0
XN3D=>CH_AERO_MODEL(KTO)%XN3D
XRG3D=>CH_AERO_MODEL(KTO)%XRG3D
XSIG3D=>CH_AERO_MODEL(KTO)%XSIG3D
XM3D=>CH_AERO_MODEL(KTO)%XM3D
XSEDA=>CH_AERO_MODEL(KTO)%XSEDA
XVDEPAERO=>CH_AERO_MODEL(KTO)%XVDEPAERO
XCTOTA3D=>CH_AERO_MODEL(KTO)%XCTOTA3D
XFTEST=>CH_AERO_MODEL(KTO)%XFTEST
XCTOTA=>CH_AERO_MODEL(KTO)%XCTOTA
XCCTOT=>CH_AERO_MODEL(KTO)%XCCTOT
XCTOTG=>CH_AERO_MODEL(KTO)%XCTOTG
XFRAC=>CH_AERO_MODEL(KTO)%XFRAC
XMI=>CH_AERO_MODEL(KTO)%XMI
XSOLORG=>CH_AERO_MODEL(KTO)%XSOLORG
XRHOP0=>CH_AERO_MODEL(KTO)%XRHOP0
XRHOP3D=>CH_AERO_MODEL(KTO)%XRHOP3D
XLAMBDA=>CH_AERO_MODEL(KTO)%XLAMBDA
XMU=>CH_AERO_MODEL(KTO)%XMU
XOM=>CH_AERO_MODEL(KTO)%XOM
XSO4RAT=>CH_AERO_MODEL(KTO)%XSO4RAT
GSEDFIX=>CH_AERO_MODEL(KTO)%GSEDFIX
IF (LHOOK) CALL DR_HOOK('MODD_CH_AERO_N:CH_AERO_GOTO_MODEL',1,ZHOOK_HANDLE)
END SUBROUTINE CH_AERO_GOTO_MODEL

END MODULE MODD_CH_AERO_n
