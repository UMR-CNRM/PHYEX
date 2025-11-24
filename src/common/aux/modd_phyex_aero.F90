MODULE MODD_PHYEX_AERO
  
IMPLICIT NONE

TYPE PHYEX_AERO_t

LOGICAL           :: LORILAM     ! switch to active aerosols fluxes
REAL              :: XINISIGI    ! dispersion initialization for I mode 
REAL              :: XINISIGJ    ! dispersion initialization for J mode
REAL              :: XINIRADIUSJ ! mean radius initialization for J mode (um)
CHARACTER(LEN=4)  :: CRGUNIT     ! type of log-normal geometric mean radius

REAL, DIMENSION(:), POINTER ::  XRHOI ! volumar mass of species in [kg/m3]
INTEGER :: NCARB     ! number of chemically inert species (like black carbon)
INTEGER :: NSOA      ! number of condensable species that may form secondary aerosols

INTEGER :: NSP       ! number of chemical species for ARES or isorropia NSP=4 these are

INTEGER :: JP_AER_OC
INTEGER :: JP_AER_BC
INTEGER :: JP_AER_DST
INTEGER :: JP_AER_H2O
INTEGER :: JP_AER_SO4

! modd_salt
LOGICAL      :: LSALT   ! switch to active pronostic sea salts
INTEGER      :: NMODE_SLT  ! number of sea salt modes (max 3; default = 3)
!Initial dry number median radius (um) from Schultz et al., 2004
REAL, DIMENSION(:), POINTER :: XINIRADIUS_SLT
!Initial, standard deviation from Vignati et al., 2001
REAL, DIMENSION(:), POINTER :: XINISIG_SLT
CHARACTER(LEN=4)  :: CRGUNITS  ! type of log-normal geometric mean radius


! modd_dust
LOGICAL      :: LDUST  ! switch to active pronostic dusts
!2 modes will be mode 2 & 3, whereas 3 modes will modes 1, 2 and 3
INTEGER, DIMENSION(3) :: JPDUSTORDER
! NEW PARAMETERIZATION FROM AMMA, default
!Initial dry number median radius (um) 
REAL, DIMENSION(:), POINTER :: XINIRADIUS
!Initial, standard deviation
REAL, DIMENSION(:), POINTER :: XINISIG
CHARACTER(LEN=4)  :: CRGUNITD  ! type of log-normal geometric mean radius
!                              ! given in namelist (mass on number)


! modd_csts_salt and modd_csts_dust
REAL  :: XDENSITY_SALT  ![kg/m3] density of dust
REAL  :: XDENSITY_DUST  ![kg/m3] density of dust

END TYPE PHYEX_AERO_t

END MODULE MODD_PHYEX_AERO
