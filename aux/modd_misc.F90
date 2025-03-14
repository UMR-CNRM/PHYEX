MODULE MODD_MISC
!
USE MODD_BUDGET, ONLY: TBUDGETCONF_t
IMPLICIT NONE
! 
!> @file  
!!    MODD_MISC - definition of a structure containing all the control parameters for arome
!!
!!    This is a structure specifically built for arome.
!!    All the constants needed to call the parametrisations in arome which are not
!!    contained in the different structures entering the parametrisation are put here.
!!    These keys are mainly in modd_conf in mesonh.
TYPE MISC_t
  !These values are not (yet) tuneable in arome.
  LOGICAL                        :: LMFCONV=.TRUE.       !< Use convective mass flux in the condensation scheme
  LOGICAL                        :: OCOMPUTE_SRC=.TRUE.  !< Compute s'r'
  INTEGER                        :: KHALO=1              !< Size of the halo for parallel distribution (used in turb)
  CHARACTER(LEN=6)               :: CPROGRAM='AROME'     !< Name of the model
  LOGICAL                        :: ONOMIXLG=.FALSE.     !< Turbulence for lagrangian variables
  LOGICAL                        :: OOCEAN=.FALSE.       !< Ocean version of the turbulence scheme
  LOGICAL                        :: ODEEPOC=.FALSE.      !< Ocean version of the turbulence scheme
  LOGICAL                        :: OCOUPLES=.FALSE.     !< Ocean-atmo LES interactive coupling
  LOGICAL                        :: OBLOWSNOW=.FALSE.    !< Blowsnow
  REAL                           :: XRSNOW=1.            !< Blowing factor
  CHARACTER(LEN=4), DIMENSION(2) :: HLBCX='CYCL'         !< Boundary condition
  CHARACTER(LEN=4), DIMENSION(2) :: HLBCY='CYCL'         !< Boundary condition
  LOGICAL                        :: OIBM=.FALSE.         !< Run with IBM
  LOGICAL                        :: OFLYER=.FALSE.       !< MesoNH flyer diagnostic
  LOGICAL                        :: ODIAG_IN_RUN=.FALSE. !< LES diagnostics
  LOGICAL                        :: O2D=.FALSE.          !< 2D version of the turbulence
  CHARACTER(LEN=4)               :: CELEC='NONE'         !< Name of the electricity scheme
  LOGICAL                        :: OELEC=.FALSE.        !< Lightning prognostic scheme
  LOGICAL                        :: OSEDIM_BEARD=.FALSE. !< Switch for effect of electrical forces on sedim.
  !These values are computed from the model setup
  LOGICAL                        :: OFLAT                !< Flat configuration

  !Budget configuration
  TYPE(TBUDGETCONF_t)            :: TBUCONF              !< Budget configuration
END TYPE MISC_t
END MODULE MODD_MISC
