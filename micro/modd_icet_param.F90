
MODULE MODD_ICET_PARAM

!!    AUTHOR
!!    ------
!!      B.J.K.Engdahl   2022 add ICE-T
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        2022
!!      Rolf H. Myhre   2024 Separated out to adapt to CY49T1
!-------------------------------------------------------------------------------

  USE MODD_PRECISION,  ONLY: MNHREAL64
  USE MODD_PARAMETERS, ONLY: JPMODELMAX

  IMPLICIT NONE

  ! Ice initiates with this mass (kg), corresponding diameter calc.
  ! Min diameters and mass of cloud, rain, snow, and graupel (m, kg).
  REAL, PARAMETER :: XD0C = 1.E-6
  REAL, PARAMETER :: XD0R = 50.E-6
  REAL, PARAMETER :: XD0G = 250.E-6
  REAL, PARAMETER :: XD0S = 200.E-6

  ! Minimum microphys values
  ! R1 value, 1.E-12, cannot be set lower because of numerical
  ! problems with Paul Field's moments and should not be set larger
  ! because of truncation problems in snow/ice growth.
  REAL, PARAMETER :: XICET_EPS = 1.E-15
  REAL, PARAMETER :: XICET_R1  = 1.E-12

  ! Lookup table dimensions
  INTEGER, PARAMETER :: NBINS = 100
  INTEGER, PARAMETER :: NBC = NBINS
  INTEGER, PARAMETER :: NBR = NBINS
  INTEGER, PARAMETER :: NBS = NBINS
  INTEGER, PARAMETER :: NTB_C = 37
  INTEGER, PARAMETER :: NTB_R = 37
  INTEGER, PARAMETER :: NTB_G = 28
  INTEGER, PARAMETER :: NTB_R1 = 37
  INTEGER, PARAMETER :: NTB_IN = 55

  !************ BJKE for ice nucleation 28.11.16 *************************
  ! Constants in Cooper curve relation for cloud ice number.
  REAL, PARAMETER :: XTNO = 5.0
  REAL, PARAMETER :: XATO = 0.304

  !******************** BJKE Bigg (1953) freezing ++ *********************
  ! Generalized gamma distributions for rain, graupel and cloud ice.
  ! N(D) = N_0 * D**mu * exp(-lamda*D);  mu=0 is exponential.
  REAL, PARAMETER :: XMU_R = 0.0
  REAL, PARAMETER :: XMU_G = 0.0
  REAL, PARAMETER :: XMU_I = 0.0

  ! Sum of two gamma distrib for snow (Field et al. 2005).
  ! N(D) = M2**4/M3**3 * [Kap0*exp(-M2*Lam0*D/M3)
  !      + Kap1*(M2/M3)**mu_s * D**mu_s * exp(-M2*Lam1*D/M3)]
  ! M2 and M3 are the (XBS)th and (XBS+1)th moments respectively
  ! calculated as function of ice water content and temperature.
  REAL, PARAMETER :: XMU_S = 0.6357

  ! Densities of rain, snow, graupel, and cloud ice.
  REAL, PARAMETER :: XRHO_G = 500.0
  REAL, PARAMETER :: XRHO_I = 890.0

  REAL, PARAMETER :: XBM_G = 3.0
  REAL, PARAMETER :: XBM_I = 3.0

  ! Rho_not used in fallspeed relations (rho_not/rho)**.5 adjustment.
  REAL, PARAMETER :: XRHO_NOT = 101325.0/(287.05*298.0)

  ! Fallspeed power laws relations:  v = (av*D**bv)*exp(-fv*D)
  ! Rain from Ferrier (1994), ice, snow, and graupel from
  ! Thompson et al (2008). Coefficient fv is zero for graupel/ice.
  REAL, PARAMETER :: XBV_R = 1.0
  REAL, PARAMETER :: XFV_S = 100.0
  REAL, PARAMETER :: XBV_I = 1.0
  REAL, PARAMETER :: XBV_C = 2.0

  ! Water vapor and air gas constants at constant pressure
  REAL, PARAMETER :: XR_THOM = 287.04

  ! Y-intercept parameter for graupel is not constant and depends on
  ! mixing ratio.  Also, when mu_g is non-zero, these become equiv
  ! y-intercept for an exponential distrib and proper values are
  ! computed based on same mixing ratio and total number concentration.
  REAL, PARAMETER :: XGONV_MIN = 1.E4
  REAL, PARAMETER :: XGONV_MAX = 3.E6

  REAL, PARAMETER, DIMENSION(2) :: XTHVREFZ=300. ! Thetav(z) for reference
                                                 ! state without orography

  ! Lookup tables for cloud water content (kg/m**3).
  REAL, DIMENSION(NTB_C), PARAMETER :: &
    XR_C = (/1.E-6,2.E-6,3.E-6,4.E-6,5.E-6,6.E-6,7.E-6,8.E-6,9.E-6, &
             1.E-5,2.E-5,3.E-5,4.E-5,5.E-5,6.E-5,7.E-5,8.E-5,9.E-5, &
             1.E-4,2.E-4,3.E-4,4.E-4,5.E-4,6.E-4,7.E-4,8.E-4,9.E-4, &
             1.E-3,2.E-3,3.E-3,4.E-3,5.E-3,6.E-3,7.E-3,8.E-3,9.E-3, &
             1.E-2/)

  ! Lookup tables for rain content (kg/m**3).
  REAL, DIMENSION(NTB_R), PARAMETER :: &
    XR_R = (/1.E-6,2.E-6,3.E-6,4.E-6,5.E-6,6.E-6,7.E-6,8.E-6,9.E-6, &
             1.E-5,2.E-5,3.E-5,4.E-5,5.E-5,6.E-5,7.E-5,8.E-5,9.E-5, &
             1.E-4,2.E-4,3.E-4,4.E-4,5.E-4,6.E-4,7.E-4,8.E-4,9.E-4, &
             1.E-3,2.E-3,3.E-3,4.E-3,5.E-3,6.E-3,7.E-3,8.E-3,9.E-3, &
             1.E-2/)

  ! Lookup tables for graupel content (kg/m**3).
  REAL, DIMENSION(NTB_G), PARAMETER :: &
    XR_G = (/1.E-5,2.E-5,3.E-5,4.E-5,5.E-5,6.E-5,7.E-5,8.E-5,9.E-5, &
             1.E-4,2.E-4,3.E-4,4.E-4,5.E-4,6.E-4,7.E-4,8.E-4,9.E-4, &
             1.E-3,2.E-3,3.E-3,4.E-3,5.E-3,6.E-3,7.E-3,8.E-3,9.E-3, &
             1.E-2/)

  ! Lookup tables for IN concentration (/m**3) from 0.001 to 1000/Liter.
  REAL, DIMENSION(NTB_IN), PARAMETER :: &
    XNT_IN = (/1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0,9.0, &
              1.E1,2.E1,3.E1,4.E1,5.E1,6.E1,7.E1,8.E1,9.E1, &
              1.E2,2.E2,3.E2,4.E2,5.E2,6.E2,7.E2,8.E2,9.E2, &
              1.E3,2.E3,3.E3,4.E3,5.E3,6.E3,7.E3,8.E3,9.E3, &
              1.E4,2.E4,3.E4,4.E4,5.E4,6.E4,7.E4,8.E4,9.E4, &
              1.E5,2.E5,3.E5,4.E5,5.E5,6.E5,7.E5,8.E5,9.E5, &
              1.E6/)

  ! Lookup tables for rain y-intercept parameter (/m**4).
  REAL, DIMENSION(NTB_R1), PARAMETER :: &
    XN0R_EXP = (/1.E6,2.E6,3.E6,4.E6,5.E6,6.E6,7.E6,8.E6,9.E6, &
                 1.E7,2.E7,3.E7,4.E7,5.E7,6.E7,7.E7,8.E7,9.E7, &
                 1.E8,2.E8,3.E8,4.E8,5.E8,6.E8,7.E8,8.E8,9.E8, &
                 1.E9,2.E9,3.E9,4.E9,5.E9,6.E9,7.E9,8.E9,9.E9, &
                 1.E10/)

  ! For snow moments conversions (from Field et al. 2005)
  REAL, DIMENSION(10), PARAMETER :: &
    XSA = (/ 5.065339, -0.062659, -3.032362, 0.029469, -0.000285, &
             0.31255,   0.000204,  0.003199, 0.0,      -0.015952/)
  REAL, DIMENSION(10), PARAMETER :: &
    XSB = (/ 0.476221, -0.015896,  0.165977, 0.007468, -0.000141, &
             0.060366,  0.000079,  0.000594, 0.0, -0.003577/)

  TYPE ICET_PARAM_t

    REAL :: XM0G

    REAL :: XODTS, XOBMR
    REAL, DIMENSION(5,15) :: XCC_EX, XCC_GM
    REAL, DIMENSION(15)   :: XOCG1, XOCG2
    REAL, DIMENSION(12)   :: XCG_E, XCG_G
    REAL :: XOGE1, XOGG1, XOGG2, XOBMG, XT1_QG_QC

    ! Mass power law relations:  mass = am*D**bm
    ! Snow from Field et al. (2005), others assume spherical form.
    REAL :: XAM_G, XAM_R, XAM_I
    INTEGER :: NIR2, NIR3
    INTEGER :: NIIN2
    INTEGER :: NIC2, NIC1
    REAL(KIND=MNHREAL64), DIMENSION(NBINS+1) :: XDX
    REAL(KIND=MNHREAL64), DIMENSION(NBC)     :: XT_NC
    REAL(KIND=MNHREAL64), DIMENSION(NBS)     :: XITDS

    ! Variables holding a bunch of exponents and gamma values (cloud water,
    ! cloud ice, rain, snow, then graupel).
    REAL, DIMENSION(13) :: XCR_EX, XCR_GM
    REAL, DIMENSION(18) :: XCS_EX, XCS_GM
    REAL, DIMENSION(7)  :: XCI_EX, XCI_GM
    REAL :: XOAMS, XT1_QS_QC
    REAL :: XORE1, XORG1, XORG2

    ! Lookup tables
    REAL(KIND=MNHREAL64), ALLOCATABLE, DIMENSION(:,:,:,:) :: XTPI_QRFZ, XTPG_QRFZ
    REAL(KIND=MNHREAL64), ALLOCATABLE, DIMENSION(:,:,:,:) :: XTPI_QCFZ
    REAL(KIND=MNHREAL64), ALLOCATABLE, DIMENSION(:,:)     :: XT_EFRW
    REAL(KIND=MNHREAL64), ALLOCATABLE, DIMENSION(:,:)     :: XT_EFSW
    REAL(KIND=MNHREAL64), ALLOCATABLE, DIMENSION(:)       :: XITDC, XITDTC
    REAL(KIND=MNHREAL64), ALLOCATABLE, DIMENSION(:)       :: XITDR, XITDTR

  CONTAINS

    PROCEDURE, PUBLIC :: ALLOC
    PROCEDURE, PUBLIC :: DEALLOC

  END TYPE ICET_PARAM_t

  TYPE(ICET_PARAM_t), DIMENSION(JPMODELMAX), TARGET, SAVE :: ICE_T_PARAMETERS_MODEL
  TYPE(ICET_PARAM_t), POINTER, SAVE :: ICE_T_PARAMETERS => NULL()

CONTAINS

  SUBROUTINE ICET_PARAM_GOTO_MODEL(KTO)

    !! Set ICE_T_PARAMETERS pointer to element KTO in MODEL list

    INTEGER, INTENT(IN) :: KTO

    ICE_T_PARAMETERS => ICE_T_PARAMETERS_MODEL(KTO)

  END SUBROUTINE ICET_PARAM_GOTO_MODEL

  SUBROUTINE ALLOC(SELF)

    !! Allocate arrays in ICET_PARAM_t type

    CLASS(ICET_PARAM_t), INTENT(INOUT) :: SELF

    ALLOCATE(SELF%XTPI_QRFZ(NTB_R, NTB_R1, 45, NTB_IN))
    ALLOCATE(SELF%XTPG_QRFZ(NTB_R, NTB_R1, 45, NTB_IN))

    ALLOCATE(SELF%XTPI_QCFZ(NTB_C, NBC, 45, NTB_IN))

    ALLOCATE(SELF%XT_EFRW(NBR, NBC))
    ALLOCATE(SELF%XT_EFSW(NBS, NBC))

    ALLOCATE(SELF%XITDC(NBC))
    ALLOCATE(SELF%XITDTC(NBC))

    ALLOCATE(SELF%XITDR(NBR))
    ALLOCATE(SELF%XITDTR(NBR))

  END SUBROUTINE ALLOC

  SUBROUTINE DEALLOC(SELF)

    !! Dellocate arrays in ICET_PARAM_t type
    !! Not used anywhere at the moment,
    !! but everything will be automatically deallocated when
    !! object goes out of scope.

    CLASS(ICET_PARAM_t), INTENT(INOUT) :: SELF

    DEALLOCATE(SELF%XTPI_QRFZ)
    DEALLOCATE(SELF%XTPG_QRFZ)

    DEALLOCATE(SELF%XTPI_QCFZ)

    DEALLOCATE(SELF%XT_EFRW)
    DEALLOCATE(SELF%XT_EFSW)

    DEALLOCATE(SELF%XITDC)
    DEALLOCATE(SELF%XITDTC)

    DEALLOCATE(SELF%XITDR)
    DEALLOCATE(SELF%XITDTR)

  END SUBROUTINE DEALLOC

END MODULE MODD_ICET_PARAM
