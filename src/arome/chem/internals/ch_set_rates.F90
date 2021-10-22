!     ######spl
      SUBROUTINE CH_SET_RATES(PTIME,PCONC,TPM,KMI,KOUT,KVERB,KVECNPT,KEQ)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    #######################################################
!! This code has been created automatically by preprocessor m10,
!! version: 9.7, copyright 1995-1999 by Meteo France/Universite Paul Sabatier.
!! Please report all bugs to K. Suhre (Lab. d'Aerologie UPS/CNRS).
!!*** *MODD_CH_SET_RATES*
!!
!!    PURPOSE
!!    -------
!       set or calculate reaction rates
!!
!!**  METHOD
!!    ------
!!      simple
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    Karsten Suhre (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 26/07/96
!!    Modified 05/05/98: Vectorization (Vincent Crassier & KS)
!!
!!----------------------------------------------------------------------
!!
!!    EXTERNAL
!!    --------
!!    none
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
USE MODD_CH_M9_SCHEME
USE MODD_CH_M9,        ONLY : METEOTRANSTYPE
USE MODI_CH_ALLOCATE_TACCS
USE MODI_CH_DEALLOCATE_TACCS
!     USER DEFINED FUNCTIONS
USE MODI_TROE
USE MODI_TROE_EQUIL
!!
!!    EXPLICIT ARGUMENTS
!!    ------------------
IMPLICIT NONE
REAL,    INTENT(IN)                      :: PTIME
INTEGER, INTENT(IN)                      :: KVECNPT
INTEGER, INTENT(IN)                      :: KEQ
REAL,    INTENT(IN),  DIMENSION(KVECNPT,KEQ)        :: PCONC
TYPE(METEOTRANSTYPE), DIMENSION(KVECNPT), INTENT(IN):: TPM
INTEGER, INTENT(IN)                      :: KMI
INTEGER, INTENT(IN)                      :: KOUT,KVERB
!!
! /BEGIN_SET_RATES/
 !
 ! transfer of meteo-variables into variables used by the 
 ! chemical core system (and some unit conversion):
 !
 !   molecular weight of air:   m_mol^air = 28.8 g/mol
 !   molecular weight of H2O:   m_mol^H2O = 18.0 g/mol
 !                        ==>   m_mol^air / m_mol^H2O = 1.6
 !   density conversion factor: 1 g/cm3 = 1E+3 kg/m3
 !   n_molec (moelc./cm3):      M = 1E-3*RHO(kg/m3) * Navo / m_mol
 !   n_water:                   H2O = M * 1.6 * Rv
 !   pressure = RHO * T * R * 1E-5 (in atm)
 !   assuming 20.95 vol% O2
 !
  TYPE(CCSTYPE),  POINTER     :: TPK
 !!
 !!----------------------------------------------------------------------
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_SET_RATES',0,ZHOOK_HANDLE)
IF (.NOT. ASSOCIATED(TACCS(KMI)%NVERB)) THEN
   CALL CH_ALLOCATE_TACCS(KMI,KVECNPT)
END IF
 
 IF (SIZE(TACCS(KMI)%NVERB) .NE. KVECNPT) THEN
  CALL CH_DEALLOCATE_TACCS(KMI)
  CALL CH_ALLOCATE_TACCS(KMI,KVECNPT)
 END IF
 
TPK=>TACCS(KMI)

 TPK%MODELLEVEL = TPM%XMETEOVAR(1)
 TPK%M          = 1E-3*TPM%XMETEOVAR(2) * 6.0221367E+23 / 28.9644
 TPK%T          = TPM%XMETEOVAR(3)
 TPK%H2O        = TPK%M*1.6077*TPM%XMETEOVAR(4)
 TPK%CLOUDWATER = TPM%XMETEOVAR(5)
 TPK%LAT        = TPM%XMETEOVAR(6)
 TPK%LON        = TPM%XMETEOVAR(7)
 TPK%YEAR       = INT(TPM%XMETEOVAR(8))
 TPK%MONTH      = INT(TPM%XMETEOVAR(9))
 TPK%DAY        = INT(TPM%XMETEOVAR(10))
 ! derived variables
 TPK%PRESSURE   = TPM%XMETEOVAR(2) * TPK%T * 288.290947 * 1E-5
 TPK%O2         = 0.2095 * TPK%M
 TPK%N2         = TPK%M - TPK%O2
 TPK%H2         = 1.23e13
 !
 ! the following prints will be erased
 IF (KVERB >= 15) THEN
   WRITE(KOUT,*) "CH_SET_RATES: the following variables have been updated"
   WRITE(KOUT,*) "MODELLEVEL:   ", TPK%MODELLEVEL(1)
   WRITE(KOUT,*) "M:            ", TPK%M(1)          , "molec/cm3"
   WRITE(KOUT,*) "T:            ", TPK%T(1)          , "K"
   WRITE(KOUT,*) "H2O           ", TPK%H2O(1)        , "molec/cm3"
   WRITE(KOUT,*) "CLOUDWATER:   ", TPK%CLOUDWATER(1) , "kg/kg"
   WRITE(KOUT,*) "LATITUDE:     ", TPK%LAT(1)        , "degree"
   WRITE(KOUT,*) "LONGITUDE:    ", TPK%LON(1)        , "degree"
   WRITE(KOUT,*) "YEAR:         ", TPK%YEAR(1) 
   WRITE(KOUT,*) "MONTH:        ", TPK%MONTH(1) 
   WRITE(KOUT,*) "DAY:          ", TPK%DAY(1) 
   WRITE(KOUT,*) "PRESSURE:     ", TPK%PRESSURE(1)   , "atm"
   WRITE(KOUT,*) "O2:           ", TPK%O2(1)         , "molec/cm3"
   WRITE(KOUT,*) "N2:           ", TPK%N2(1)         , "molec/cm3"
   WRITE(KOUT,*) "H2:           ", TPK%H2(1)         , "molec/cm3"
 END IF
 !
! /END_SET_RATES/
 TPK%K018=TPK%M*6.00E-34*(TPK%T/300)**(-2.3)
 TPK%K019=8.00E-12*exp(-(2060.0/TPK%T))
 TPK%K020=1.80E-11*exp(-(-110.0/TPK%T))
 TPK%K021=3.20E-11*exp(-(-70.0/TPK%T))
 TPK%K022=2.20E-10
 TPK%K023=1.60E-12*exp(-(940.0/TPK%T))
 TPK%K024=1.10E-14*exp(-(500.0/TPK%T))
 TPK%K025=4.80E-11*exp(-(-250.0/TPK%T))
 TPK%K026=2.90E-12*exp(-(160.0/TPK%T))
 TPK%K027=2.3E-13*EXP(600./TPK%T)+1.7E-33*TPK%M*EXP(1000./TPK%T)
 TPK%K028=3.22E-34*EXP(2800./TPK%T)+2.38E-54*TPK%M*EXP(3200./TPK%T)
 TPK%K029=TROE(1.,9.00E-32,1.5,3.00E-11,0.0,TPK%M,TPK%T,KVECNPT)
 TPK%K030=6.50E-12*exp(-(-120.0/TPK%T))
 TPK%K031=TROE(1.,9.00E-32,2.0,2.20E-11,0.0,TPK%M,TPK%T,KVECNPT)
 TPK%K032=TROE(1.,7.00E-31,2.6,1.50E-11,0.5,TPK%M,TPK%T,KVECNPT)
 TPK%K033=TROE(1.,2.60E-30,3.2,2.40E-11,1.3,TPK%M,TPK%T,KVECNPT)
 TPK%K034=2.20E-11
 TPK%K035=3.70E-12*exp(-(-250.0/TPK%T))
 TPK%K036=TROE(1.,1.80E-31,3.2,4.70E-12,1.4,TPK%M,TPK%T,KVECNPT)
 TPK%K037=TROE_EQUIL(1.80E-31,3.2,4.70E-12,1.4,4.76E+26,10900.,TPK%M,TPK%T,KVEC&
&NPT)
 TPK%K038=3.50E-12
 TPK%K039=1.80E-11*exp(-(390.0/TPK%T))
 TPK%K040=(7.2E-15*EXP(785/TPK%T))+(1.9E-33*EXP(725/TPK%T)*TPK%M)/(1+(1.9E-33*E&
&XP(725/TPK%T)*TPK%M)/(4.1E-16*EXP(1440/TPK%T)))
 TPK%K041=1.30E-12*exp(-(-380.0/TPK%T))
 TPK%K042=2.00E-12*exp(-(1400.0/TPK%T))
 TPK%K043=1.20E-13*exp(-(2450.0/TPK%T))
 TPK%K044=3.30E-39*exp(-(-530.0/TPK%T))
 TPK%K045=1.50E-11*exp(-(-170.0/TPK%T))
 TPK%K046=4.50E-14*exp(-(1260.0/TPK%T))
 TPK%K047=TROE(1.,2.20E-30,3.9,1.50E-12,0.7,TPK%M,TPK%T,KVECNPT)
 TPK%K048=TROE_EQUIL(2.20E-30,3.9,1.50E-12,0.7,3.70E+26,11000.0,TPK%M,TPK%T,KVE&
&CNPT)
 TPK%K049=8.50E-13*exp(-(2450.0/TPK%T))
 TPK%K050=5.50E-12*exp(-(2000.0/TPK%T))
 TPK%K051=TROE(1.,3.00E-31,3.3,1.50E-12,0.0,TPK%M,TPK%T,KVECNPT)
 TPK%K052=1.5E-13*(1.+2.439E-20*TPK%M)
 TPK%K053=6.00E-11
 TPK%K054=0.00E-01*exp(-(-13.0/TPK%T))
 TPK%K055=TPK%T*TPK%T*7.44E-18*exp(-(1361./TPK%T))
 TPK%K056=1.51E-17*TPK%T*TPK%T*exp(-(492./TPK%T))
 TPK%K057=3.76E-12*exp(-(260.0/TPK%T))+1.70E-12*exp(-(155.0/TPK%T))+1.21E-12*ex&
&p(-(125.0/TPK%T))
 TPK%K058=1.78E-12*exp(-(-438.0/TPK%T))+6.07E-13*exp(-(-500.0/TPK%T))+0.00E-01*&
&exp(-(-448.0/TPK%T))
 TPK%K059=2.54E-11*exp(-(-410.0/TPK%T))+0.00E-01*exp(-(-444.0/TPK%T))+0.00E-01
 TPK%K060=3.31E-12*exp(-(-355.0/TPK%T))+3.45E-13
 TPK%K061=1.00E-11
 TPK%K062=5.55E-12*exp(-(-331.0/TPK%T))
 TPK%K063=TPK%T*TPK%T*5.68E-18*exp(-(-92.0/TPK%T))
 TPK%K064=1.32E-11+1.88E-12*exp(-(-175.0/TPK%T))
 TPK%K065=2.93E-12*exp(-(-190.0/TPK%T))
 TPK%K066=3.36E-12*exp(-(-190.0/TPK%T))
 TPK%K067=3.80E-14+1.59E-14*exp(-(-500.0/TPK%T))
 TPK%K068=5.31E-12*exp(-(260.0/TPK%T))
 TPK%K069=3.40E-13*exp(-(1900.0/TPK%T))
 TPK%K070=1.40E-12*exp(-(1900.0/TPK%T))
 TPK%K071=1.62E-12*exp(-(1900.0/TPK%T))+0.00E-01*exp(-(150.0/TPK%T))+1.94E-14*e&
&xp(-(1000.0/TPK%T))
 TPK%K072=4.92E-16
 TPK%K073=4.35E-18*TPK%T*TPK%T*exp(-(2282.0/TPK%T))+1.91E-14*exp(-(450.0/TPK%T)&
&)+1.08E-15*exp(-(-450.0/TPK%T))+0.00E-01
 TPK%K074=4.00E-12*exp(-(446.0/TPK%T))+0.00E-01*exp(-(-490.0/TPK%T))+0.00E-01
 TPK%K075=3.76E-16*exp(-(500.0/TPK%T))
 TPK%K076=8.17E-15*exp(-(2580.0/TPK%T))+4.32E-16*exp(-(1800.0/TPK%T))+2.87E-17*&
&exp(-(845.0/TPK%T))+0.00E-01*exp(-(2283.0/TPK%T))
 TPK%K077=7.86E-15*exp(-(1913.0/TPK%T))+0.00E-01*exp(-(732.0/TPK%T))+0.00E-01
 TPK%K078=0.00E-01*exp(-(2112.0/TPK%T))+1.38E-19
 TPK%K079=7.20E-17*exp(-(1700.0/TPK%T))
 TPK%K080=2.00E-11
 TPK%K081=1.00E-11
 TPK%K082=3.60E-11
 TPK%K083=1.66E-17*exp(-(-1044.0/TPK%T))
 TPK%K084=2.80E-11
 TPK%K085=TROE(5.86E-01,9.70E-29,5.6,9.30E-12,1.5,TPK%M,TPK%T,KVECNPT)
 TPK%K086=TROE_EQUIL(9.70E-29,5.6,9.30E-12,1.5,1.16E+28,13954.,TPK%M,TPK%T,KVEC&
&NPT)
 TPK%K087=4.20E-12*exp(-(-180.0/TPK%T))
 TPK%K088=4.36E-12
 TPK%K089=6.93E-12
 TPK%K090=4.00E-12
 TPK%K091=4.00E-12
 TPK%K092=1.22E-11
 TPK%K093=4.00E-12
 TPK%K094=3.80E-13*exp(-(-800.0/TPK%T))
 TPK%K095=6.16E-14*exp(-(-700.0/TPK%T))+1.52E-13*exp(-(-1300.0/TPK%T))
 TPK%K096=1.81E-13*exp(-(-1300.0/TPK%T))
 TPK%K097=1.28E-13*exp(-(-1300.0/TPK%T))+0.00E-01
 TPK%K098=3.75E-13*exp(-(-980.0/TPK%T))
 TPK%K099=5.94E-13*exp(-(-550.0/TPK%T))+1.99E-16*exp(-(-2640.0/TPK%T))+5.56E-14&
&*exp(-(-1300.0/TPK%T))
 TPK%K100=1.66E-13*exp(-(-1300.0/TPK%T))
 TPK%K101=9.10E-14*exp(-(-416.0/TPK%T))
 TPK%K102=1.03E-14*exp(-(-158.0/TPK%T))+6.24E-14*exp(-(-431.0/TPK%T))+1.53E-14*&
&exp(-(-467.0/TPK%T))+4.34E-15*exp(-(-633.0/TPK%T))
 TPK%K103=1.57E-13*exp(-(-708.0/TPK%T))
 TPK%K104=1.36E-13*exp(-(-708.0/TPK%T))
 TPK%K105=3.56E-14*exp(-(-708.0/TPK%T))
 TPK%K106=1.77E-11*exp(-(440.0/TPK%T))+1.48E-16*exp(-(-2510.0/TPK%T))+3.10E-13*&
&exp(-(-508.0/TPK%T))
 TPK%K107=1.12E-13*exp(-(-708.0/TPK%T))
 TPK%K108=4.44E-14*exp(-(-211.0/TPK%T))+2.23E-13*exp(-(-460.0/TPK%T))+4.10E-14*&
&exp(-(-522.0/TPK%T))+1.17E-14*exp(-(-683.0/TPK%T))
 TPK%K109=4.36E-13*exp(-(-765.0/TPK%T))
 TPK%K110=7.60E-13*exp(-(-765.0/TPK%T))
 TPK%K111=3.63E-13*exp(-(-765.0/TPK%T))
 TPK%K112=7.73E-13*exp(-(-530.0/TPK%T))+1.70E-13*exp(-(-565.0/TPK%T))
 TPK%K113=4.85E-13*exp(-(-765.0/TPK%T))
 TPK%K114=4.19E-15*exp(-(-1000.0/TPK%T))
 TPK%K115=2.48E-14*exp(-(-1000.0/TPK%T))
 TPK%K116=1.20E-12
 TPK%K117=1.20E-12
 TPK%K118=1.20E-12
 TPK%K119=1.20E-12
 TPK%K120=1.20E-12
 TPK%K121=3.48E-12
 TPK%K122=1.20E-12
 TPK%K123=1.66E-13*exp(-(-1300.0/TPK%T))
 TPK%K124=5.99E-15*exp(-(-1510.0/TPK%T))
 TPK%K125=1.69E-14*exp(-(-1560.0/TPK%T))
 TPK%K126=7.13E-17*exp(-(-2950.0/TPK%T))
 TPK%K127=4.00E-12
 TPK%K128=1.20E-12
TPK%NOUT(:)  = KOUT
TPK%NVERB(:) = KVERB
IF (LHOOK) CALL DR_HOOK('CH_SET_RATES',1,ZHOOK_HANDLE)
RETURN
END SUBROUTINE CH_SET_RATES
