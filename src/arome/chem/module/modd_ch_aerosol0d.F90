!     ######spl
       MODULE MODD_CH_AEROSOL0D
!!     ######################
!!
!!     PURPOSE
!!     -------
!!
!!     declaration of variables and types for the aerosol system in the
!!     surface meteo-france scheme
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
!!     P. Tulet (CNRM)
!!
!!
!!     MODIFICATIONS
!!     -------------
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
IMPLICIT NONE
!
! aerosol mode parameters
LOGICAL      :: LCH_AERO_FLUX     = .FALSE. ! switch to active pronostic aerosols
!
LOGICAL      :: LCO2PM            = .FALSE. ! switch to active primary emission derived from CO 

LOGICAL      :: LVARSIGI  = .FALSE.   ! switch to active pronostic dispersion for I mode
LOGICAL      :: LVARSIGJ  = .FALSE.   ! switch to active pronostic dispersion for J mode

!
INTEGER, PARAMETER         :: JPMODE=2      ! number of modes
INTEGER, PARAMETER         :: JPIN=JPMODE*3 ! number of differential equations
INTEGER, SAVE, DIMENSION(JPMODE) :: NM0,NM3,NM6   ! index of the moments in arrays
!
!* indices of Aerosol chemical parameters
!
INTEGER, PARAMETER :: NSP=4        ! number of chemical species
                                   ! for ARES or isorropia NSP=4 these are
INTEGER, PARAMETER :: JP_AER_SO4 = 1
INTEGER, PARAMETER :: JP_AER_NO3 = 2
INTEGER, PARAMETER :: JP_AER_NH3 = 3
INTEGER, PARAMETER :: JP_AER_H2O = 4
!
INTEGER, PARAMETER :: JP_AER_SO4g = JP_AER_SO4
INTEGER, PARAMETER :: JP_AER_NO3g = JP_AER_NO3
INTEGER, PARAMETER :: JP_AER_NH3g = JP_AER_NH3
!
INTEGER, PARAMETER :: NCARB=2     ! number of chemically inert species
                                  ! (like black carbon)
INTEGER, PARAMETER :: JP_AER_OC = 5
INTEGER, PARAMETER :: JP_AER_BC = 6

INTEGER, PARAMETER :: NSOA=10      ! number of condensable species that may form
                                   ! secondary aerosols
INTEGER, PARAMETER :: JP_AER_SOA1 = 7
INTEGER, PARAMETER :: JP_AER_SOA2 = 8
INTEGER, PARAMETER :: JP_AER_SOA3 = 9
INTEGER, PARAMETER :: JP_AER_SOA4 = 10
INTEGER, PARAMETER :: JP_AER_SOA5 = 11
INTEGER, PARAMETER :: JP_AER_SOA6 = 12
INTEGER, PARAMETER :: JP_AER_SOA7 = 13
INTEGER, PARAMETER :: JP_AER_SOA8 = 14
INTEGER, PARAMETER :: JP_AER_SOA9 = 15
INTEGER, PARAMETER :: JP_AER_SOA10 = 16
                                  ! Dusts
INTEGER, PARAMETER :: JP_AER_DUST = 17
! switch to use Pun et al. (if equal
! to 1)
INTEGER, PARAMETER :: JP_AER_SOA = 1 ! 
! to compile ch_monitorn (to be change if use RACM or ReLAS)
INTEGER, PARAMETER :: JP_AER_SOAA = 1 ! 
INTEGER, PARAMETER :: JP_AER_SOAB = 1 ! 
INTEGER, PARAMETER :: JP_SOAAi = 1 ! 
INTEGER, PARAMETER :: JP_SOAAj = 1 ! 
INTEGER, PARAMETER :: JP_SOABi = 1 ! 
INTEGER, PARAMETER :: JP_SOABj = 1 ! 

INTEGER, PARAMETER :: JPNN=NSP+NSOA+NCARB+1

INTEGER :: JP_CH_SO4i, JP_CH_SO4j, &
           JP_CH_NO3i, JP_CH_NO3j, JP_CH_NH3i, JP_CH_NH3j, &
           JP_CH_BCi,  JP_CH_BCj,  JP_CH_OCi,  JP_CH_OCj,  &
           JP_CH_H2Oi, JP_CH_H2Oj, JP_CH_DUST, JP_CH_M0i,  &
           JP_CH_M0j,  JP_CH_M6i,  JP_CH_M6j
!
REAL, SAVE, DIMENSION(JPNN)  :: XMI    ! molecular mass of species i (g/mol)

REAL, SAVE, DIMENSION(JPNN)  :: XRHOI  ! volumar mass of species i (kg/m3)

REAL, SAVE, DIMENSION(JPNN)  :: XFAC   ! conversion factor um3/m3 -> ug/cm3
                                       ! for each chemical species
!
REAL         :: XEMISRADIUSI  = 0.01   ! mean radius of primary aerosol
                                       ! emission for I mode
REAL         :: XEMISRADIUSJ  = 0.5    ! mean radius of primary aerosol
                                       ! emission for J mode
REAL         :: XEMISRADIUSK  = 4.     ! mean radius of primary aerosol
                                       ! emission for K mode
REAL         :: XEMISSIGI     = 1.60   ! dispersion of primary aerosol
                                       ! emission for I mode
REAL         :: XEMISSIGJ     = 1.60   ! dispersion of primary aerosol
                                       ! emission for J mode
REAL         :: XEMISSIGK     = 1.60   ! dispersion of primary aerosol
                                       ! emission for K mode
CHARACTER*4  :: CRGUNIT   = 'MASS'    ! type of log-normal geometric mean radius given
!                                     ! in nameliste (mass on number)



!----------------------------------------------------------------------------
!
!*  constants
!
REAL, PARAMETER :: XPBOLTZ=1.380658e-23    ! Boltzmann constant (J/K)
REAL, PARAMETER :: XAVOGADRO=6.0221367E+23 ! Avogadro constant
REAL, PARAMETER :: XMD    = 28.9644E-3     ! Air mass molarity

!
END MODULE MODD_CH_AEROSOL0D
