!     ######spl
       MODULE MODD_CH_AEROSOL
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
!!     Vincent Crassier (LA)
!!
!!
!!     MODIFICATIONS
!!     -------------
!!     (30-01-01) P.Tulet (LA) * modifications for secondary biogenics aerosols
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
USE MODD_CST, ONLY :    &
       XPI              & !Definition of pi
      ,XBOLTZ           & ! Boltzman constant 
      ,XAVOGADRO        & ![molec/mol] avogadros number
      ,XG               & ! Gravity constant
      ,XP00             & ! Reference pressure
      ,XMD              & ![kg/mol] molar weight of air
      ,XRD              & ! Gaz constant for dry air
      ,XCPD               !  Cpd (dry air)
!
IMPLICIT NONE
!
! aerosol mode parameters
!
LOGICAL      :: LORILAM     = .FALSE.       ! switch to active aerosols fluxes
LOGICAL      :: LINITPM     = .TRUE.        ! switch to initialize BC
LOGICAL      :: LAERINIT    = .FALSE.       ! switch to initialize aerosols

                                            ! and OC from CO concentration (real_case)
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

INTEGER            :: NSOA = 10    ! number of condensable species that may form
                                   ! secondary aerosols
INTEGER            :: NM6_AER = 2  ! number of condensable species that may form
                                   ! secondary aerosols
INTEGER            :: JP_AER_SOA1 = 7 
INTEGER            :: JP_AER_SOA2 = 8
INTEGER            :: JP_AER_SOA3 = 9
INTEGER            :: JP_AER_SOA4 = 10
INTEGER            :: JP_AER_SOA5 = 11
INTEGER            :: JP_AER_SOA6 = 12
INTEGER            :: JP_AER_SOA7 = 13
INTEGER            :: JP_AER_SOA8 = 14
INTEGER            :: JP_AER_SOA9 = 15 
INTEGER            :: JP_AER_SOA10 = 16

CHARACTER(LEN=32),DIMENSION(:), ALLOCATABLE :: CAERONAMES

INTEGER            :: JP_CH_SO4I = 1  
INTEGER            :: JP_CH_SO4J = 2  
INTEGER            :: JP_CH_NO3I = 3  
INTEGER            :: JP_CH_NO3J = 4  
INTEGER            :: JP_CH_NH3I = 5  
INTEGER            :: JP_CH_NH3J = 6  
INTEGER            :: JP_CH_H2OI = 7  
INTEGER            :: JP_CH_H2OJ = 8  
INTEGER            :: JP_CH_OCI  = 9  
INTEGER            :: JP_CH_OCJ  = 10  
INTEGER            :: JP_CH_BCI  = 11  
INTEGER            :: JP_CH_BCJ  = 12 
INTEGER            :: JP_CH_SOA1I  = 13  
INTEGER            :: JP_CH_SOA1J  = 14
INTEGER            :: JP_CH_SOA2I  = 15
INTEGER            :: JP_CH_SOA2J  = 16
INTEGER            :: JP_CH_SOA3I  = 17
INTEGER            :: JP_CH_SOA3J  = 18  
INTEGER            :: JP_CH_SOA4I  = 19  
INTEGER            :: JP_CH_SOA4J  = 20 
INTEGER            :: JP_CH_SOA5I  = 21  
INTEGER            :: JP_CH_SOA5J  = 22  
INTEGER            :: JP_CH_SOA6I  = 23  
INTEGER            :: JP_CH_SOA6J  = 24  
INTEGER            :: JP_CH_SOA7I  = 25  
INTEGER            :: JP_CH_SOA7J  = 26  
INTEGER            :: JP_CH_SOA8I  = 27  
INTEGER            :: JP_CH_SOA8J  = 28  
INTEGER            :: JP_CH_SOA9I  = 29  
INTEGER            :: JP_CH_SOA9J  = 30  
INTEGER            :: JP_CH_SOA10I = 31  
INTEGER            :: JP_CH_SOA10J = 32  
INTEGER            :: JP_CH_M0I = 33  
INTEGER            :: JP_CH_M0J = 34  
INTEGER            :: JP_CH_M6I = 35  
INTEGER            :: JP_CH_M6J = 36  

!INTEGER, PARAMETER :: JPNN=NSP+NSOA+NCARB
                                   
! Index for gas species which interact with aerosols
INTEGER :: JP_CH_HNO3,  JP_CH_H2SO4, JP_CH_NH3, JP_CH_O3, JP_CH_CO,  &
           JP_CH_URG1, JP_CH_URG2, JP_CH_RPG2, JP_CH_RP18, JP_CH_UR26,&
           JP_CH_RPG3, JP_CH_URG4, JP_CH_UR8, JP_CH_UR17, JP_CH_UR7, JP_CH_URG6, &
           JP_CH_ARAC, JP_CH_URG7, JP_CH_RPG7, JP_CH_RPR7, JP_CH_URG8, JP_CH_UR19, &
           JP_CH_URG9, JP_CH_URG10, JP_CH_PAN8, JP_CH_UR22, JP_CH_RPR4, JP_CH_AP7, &
           JP_CH_RPR3, JP_CH_UR21, JP_CH_UR28, JP_CH_UR29,  JP_CH_UR30, &
           JP_CH_RPR9, JP_CH_RP12, JP_CH_UR3, JP_CH_UR23, JP_CH_UR31, JP_CH_AP1, &
           JP_CH_AP6, JP_CH_ADAC, JP_CH_UR2, JP_CH_UR14, JP_CH_UR27, JP_CH_RP14, &
           JP_CH_RP19, JP_CH_UR11, JP_CH_UR15, JP_CH_AP10, JP_CH_UR20, JP_CH_UR34, &
           JP_CH_AP11, JP_CH_AP12, JP_CH_UR5, JP_CH_UR6, JP_CH_AP8, JP_CH_RP17,&
           JP_CH_RP13

INTEGER :: JP_CH_H2O2,  JP_CH_SO2, JP_CH_SO42M
!
REAL, SAVE, DIMENSION(:), ALLOCATABLE ::  XRHOI  ! volumar mass of species i (kg/m3)

REAL, SAVE, DIMENSION(:), ALLOCATABLE ::  XFAC   ! conversion factor um3/m3 -> ug/cm3
                                       ! for each chemical species


! Molar mass of each aerosols parents (in kg/mol)
REAL, PARAMETER :: XHNO3=63.01287
REAL, PARAMETER :: XNH3 =17.03061
REAL, PARAMETER :: XURG1=88.
REAL, PARAMETER :: XURG2=1.76981E+02
REAL, PARAMETER :: XRPG2=1.68000E+02
REAL, PARAMETER :: XRP18=1.84000E+02
REAL, PARAMETER :: XRPG3=1.53772E+02
REAL, PARAMETER :: XURG4=1.95867E+02
REAL, PARAMETER :: XUR17=1.72000E+02
REAL, PARAMETER :: XRPR3=1.86000E+02
REAL, PARAMETER :: XAP7 =2.33000E+02
REAL, PARAMETER :: XURG6=1.89153E+02
REAL, PARAMETER :: XUR22=2.12000E+02
REAL, PARAMETER :: XURG7=1.56781E+02
REAL, PARAMETER :: XADAC=1.56781E+02
REAL, PARAMETER :: XRPR4=1.67000E+02
REAL, PARAMETER :: XRPR7=1.50000E+02
REAL, PARAMETER :: XRPG7=1.96059E+02
REAL, PARAMETER :: XURG8=1.73777E+02
REAL, PARAMETER :: XURG9=2.61676E+02
REAL, PARAMETER :: XUR26=1.68000E+02
REAL, PARAMETER :: XURG10=2.14834E+02
REAL, PARAMETER :: XUR7=1.68000E+02
REAL, PARAMETER :: XUR8=1.84000E+02
REAL, PARAMETER :: XPAN8=2.63000E+02
REAL, PARAMETER :: XARAC=1.66000E+02
REAL, PARAMETER :: XUR19=1.70000E+02
REAL, PARAMETER :: XUR21=88.
REAL, PARAMETER :: XUR28=90.
REAL, PARAMETER :: XUR29=186.0
REAL, PARAMETER :: XUR30=200.0 
REAL, PARAMETER :: XRP13=168.
REAL, PARAMETER :: XRP17=170.0
REAL, PARAMETER :: XRPR9=154.0
REAL, PARAMETER :: XRP12=152.0
REAL, PARAMETER :: XUR3=202.0
REAL, PARAMETER :: XUR23=144.0
REAL, PARAMETER :: XUR31=220.0
REAL, PARAMETER :: XAP1=183.0
REAL, PARAMETER :: XAP6=197.0
REAL, PARAMETER :: XRP14=188.0
REAL, PARAMETER :: XRP19=204.0
REAL, PARAMETER :: XUR2=152.0
REAL, PARAMETER :: XUR14=181.0
REAL, PARAMETER :: XUR27=164.0
REAL, PARAMETER :: XUR11=172.0
REAL, PARAMETER :: XUR15=201.0
REAL, PARAMETER :: XAP10=217.0
REAL, PARAMETER :: XUR20=256.0
REAL, PARAMETER :: XUR34=240.0
REAL, PARAMETER :: XAP11=287.0
REAL, PARAMETER :: XAP12=303.0
REAL, PARAMETER :: XUR5=170.0
REAL, PARAMETER :: XUR6=170.0
REAL, PARAMETER :: XAP8=215.0


!----------------------------------------------------------------------------
!
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XSURF
REAL, SAVE, DIMENSION(:,:), ALLOCATABLE :: XDP

! Declaration for the Bessagnet tabulation
REAL, SAVE, DIMENSION(:), ALLOCATABLE :: rhi
REAL, SAVE, DIMENSION(:), ALLOCATABLE :: tempi
REAL, SAVE, DIMENSION(:), ALLOCATABLE :: zsu, znh, zni
REAL, SAVE, DIMENSION(:,:,:,:,:,:), ALLOCATABLE :: zf

! Declaration of  the neuronal coefficients

!     .. weights 
REAL, SAVE, DIMENSION(100,100) :: W1IJA,W1JKA,W2IJA,W2JKA
REAL, SAVE, DIMENSION(100,100) :: W1IJB,W1JKB,W2IJB,W2JKB
REAL, SAVE, DIMENSION(100,100) :: W1IJC,W1JKC,W2IJC,W2JKC
REAL, SAVE, DIMENSION(2,100) :: X1MINA,X1MAXA,X1MODA,X2MINA,X2MAXA,X2MODA
REAL, SAVE, DIMENSION(2,100) :: X1MINB,X1MAXB,X1MODB,X2MINB,X2MAXB,X2MODB
REAL, SAVE, DIMENSION(2,100) :: X1MINC,X1MAXC,X1MODC,X2MINC,X2MAXC,X2MODC

!     .. counters and indices
INTEGER, SAVE :: I1IA,J1JA,K1KA,I2IA,J2JA,K2KA
INTEGER, SAVE :: I1IB,J1JB,K1KB,I2IB,J2JB,K2KB
INTEGER, SAVE :: I1IC,J1JC,K1KC,I2IC,J2JC,K2KC

!----------------------------------------------------------------------------
! aerosol lognormal parameterization
CHARACTER(LEN=4)  :: CRGUNIT   = 'NUMB'    ! type of log-normal geometric mean radius
!                                     !given in namelist (mass on number)

LOGICAL      :: LVARSIGI  = .FALSE.   ! switch to active pronostic dispersion for I mode
LOGICAL      :: LVARSIGJ  = .FALSE.   ! switch to active pronostic dispersion for J mode
LOGICAL      :: LVARSIGK  = .FALSE.   ! switch to active pronostic dispersion for J mode
LOGICAL      :: LHETEROSO4 = .FALSE.   ! switch to active sulfates heteronegeous
                                       ! production
LOGICAL      :: LRGFIX     = .FALSE.    ! switch to active aerosol sedimentation
LOGICAL      :: LSEDIMAERO = .FALSE.   ! switch to active aerosol sedimentation
REAL         :: XN0IMIN    = 10.       ! minimum particule number value for I mode
REAL         :: XN0JMIN    = 1.        ! minimum particule number value for J mode
REAL         :: XINIRADIUSI= 0.05      ! mean radius initialization for I mode (um)
REAL         :: XINIRADIUSJ= 0.2       ! mean radius initialization for J mode (um)
REAL         :: XINISIGI   = 1.80      ! dispersion initialization for I mode 
REAL         :: XINISIGJ   = 2.00      ! dispersion initialization for J mode
REAL         :: XSIGIMIN   = 1.20      ! minimum dispersion value for I mode
REAL         :: XSIGJMIN   = 1.20      ! minimum dispersion value for J mode
REAL         :: XSIGIMAX   = 3.60      ! maximum dispersion value for I mode
REAL         :: XSIGJMAX   = 3.60      ! maximum dispersion value for J mode
REAL         :: XCOEFRADIMAX  = 10.    ! maximum increasement for Rg mode I
REAL         :: XCOEFRADIMIN  = 0.1    ! maximum decreasement for Rg mode I
REAL         :: XCOEFRADJMAX  = 10.    ! maximum increasement for Rg mode J
REAL         :: XCOEFRADJMIN  = 0.1    ! maximum decreasement for Rg mode J
CHARACTER(LEN=5)  :: CMINERAL      = "NONE"   ! mineral equilibrium scheme
CHARACTER(LEN=5)  :: CORGANIC      = "NONE"   ! organic equilibrium scheme
CHARACTER(LEN=80) :: CNUCLEATION   = "NONE"   ! sulfates nucleation scheme
!
!
END MODULE MODD_CH_AEROSOL
