!MNH_LIC Copyright 2002-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #############################
       MODULE MODI_INI_RAIN_ICE_ELEC
!      #############################
!
INTERFACE
      SUBROUTINE INI_RAIN_ICE_ELEC (KLUOUT, PTSTEP, PDZMIN, KSPLITR, HCLOUD, &
                                    KINTVL, PFDINFTY                         )
!
INTEGER,           INTENT(IN) :: KLUOUT    ! Logical unit number for prints
INTEGER,           INTENT(OUT):: KSPLITR   ! Number of small time step
                                           ! integration for  rain
                                           ! sedimendation
REAL,              INTENT(IN) :: PTSTEP    ! Effective Time step 
REAL,              INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
CHARACTER (LEN=4), INTENT(IN) :: HCLOUD    ! Indicator of the cloud scheme
INTEGER,           INTENT(INOUT) :: KINTVL   ! Number of interval to integrate the kernels
REAL,              INTENT(INOUT) :: PFDINFTY ! Factor used to define the "infinite" diameter
!
END SUBROUTINE INI_RAIN_ICE_ELEC
END INTERFACE
END MODULE MODI_INI_RAIN_ICE_ELEC
!
!     ########################################################################
      SUBROUTINE INI_RAIN_ICE_ELEC (KLUOUT, PTSTEP, PDZMIN, KSPLITR, HCLOUD, &
                                    KINTVL, PFDINFTY                         )
!     ########################################################################
!
!!****  *INI_RAIN_ICE_ELEC * - initialize the constants necessary for the warm and
!!                             cold microphysical schemes, 
!!                             and for the electrical scheme
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used to 
!!    resolve the mixed phase microphysical scheme. The collection kernels of 
!!    the precipitating particles are recomputed if necessary if some parameters
!!    defining the ice categories have been modified. The number of small
!!    time steps leading to stable scheme for the rain, ice, snow and ggraupeln
!!    sedimentation is also computed (time-splitting technique).     
!!
!!**  METHOD
!!    ------
!!      The constants are initialized to their numerical values and the number
!!    of small time step is computed by dividing the 2* Deltat time interval of
!!    the Leap-frog scheme so that the stability criterion for the rain
!!    sedimentation is fulfilled for a Raindrop maximal fall velocity equal
!!    VTRMAX. The parameters defining the collection kernels are read and are
!!    checked against the new ones. If any change occurs, these kernels are
!!    recomputed and their numerical values are written in the output listiing. 
!!
!!    EXTERNAL
!!    --------
!!      GAMMA    :  gamma function
!!     
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_CST
!!        XPI                  !
!!        XP00                 ! Reference pressure
!!        XRD                  ! Gaz constant for dry air
!!        XRHOLW               ! Liquid water density
!!      Module MODD_REF
!!        XTHVREFZ             ! Reference virtual pot.temp. without orography
!!      Module MODD_PARAMETERS
!!        JPVEXT               !
!!      Module MODD_RAIN_ICE_DESCR
!!      Module MODD_RAIN_ICE_PARAM
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine INI_RAIN_ICE )
!!
!!    AUTHOR
!!    ------
!!
!!    * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!    Original: 2002
!!    Modifications:
!!      C. Barthe   20/11/09   update to version 4.8.1
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_LUNIT
USE MODD_PARAMETERS
USE MODD_PARAM_ICE
USE MODD_RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM
USE MODD_REF
USE MODD_ELEC_PARAM, ONLY : XGAMINC_RIM3, XFCI
USE MODD_ELEC_DESCR, ONLY : XFS
!
USE MODI_MOMG
USE MODI_GAMMA
USE MODI_GAMMA_INC
USE MODI_RRCOLSS
USE MODI_RZCOLX
USE MODI_RSCOLRG
USE MODI_READ_XKER_RACCS
USE MODI_READ_XKER_SDRYG
USE MODI_READ_XKER_RDRYG
USE MODI_READ_XKER_SWETH
USE MODI_READ_XKER_GWETH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,           INTENT(IN) :: KLUOUT   ! Logical unit number for prints
INTEGER,           INTENT(OUT):: KSPLITR   ! Number of small time step
                                           ! integration for  rain
                                           ! sedimendation
REAL,              INTENT(IN) :: PTSTEP    ! Effective Time step 
REAL,              INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
CHARACTER (LEN=4), INTENT(IN) :: HCLOUD    ! Indicator of the cloud scheme
INTEGER,           INTENT(INOUT) :: KINTVL   ! Number of interval to integrate the kernels
REAL,              INTENT(INOUT) :: PFDINFTY ! Factor used to define the "infinite" diameter
!
!
!*       0.2   Declarations of local variables :
!
REAL :: ZT                    ! Work variable
REAL :: ZVTRMAX               ! Raindrop maximal fall velocity
REAL :: ZRHO00                ! Surface reference air density
REAL :: ZRATE                 ! Geometrical growth of Lbda in the tabulated
                              ! functions and kernels
REAL :: ZBOUND                ! XDCSLIM*Lbda_s: upper bound for the partial
                              ! integration of the riming rate of the aggregates
REAL :: ZEGS, ZEGR, ZEHS, ZEHG! Bulk collection efficiencies
REAL :: ZALPHA, &  !      Parameters to compute
        ZNU,    &  !    the value of the p_moment
        ZP         ! of the generalized Gamma function
REAL :: ZESR                  ! Mean efficiency of rain-aggregate collection
REAL :: ZFDINFTY              ! Factor used to define the "infinite" diameter
REAL :: ZCONC_MAX     ! Maximal concentration for snow
REAL :: ZGAMC, ZGAMC2 ! parameters involving various moments of the generalized gamma law
REAL :: ZFACT_NUCL    ! Amplification factor for the minimal ice concentration
REAL :: ZXR           ! Value of x_r in N_r = C_r lambda_r ** x_r
!
INTEGER :: IKB     ! Coordinates of the first physical points along z
INTEGER :: J1, J2  ! Internal loop indexes
INTEGER :: IND     ! Number of interval to integrate the kernels
!
LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
                    ! listing
!  
INTEGER  :: KND
INTEGER  :: KACCLBDAS, KACCLBDAR, KDRYLBDAG, KDRYLBDAS, KDRYLBDAR
INTEGER  :: KWETLBDAS, KWETLBDAG, KWETLBDAH
REAL     :: PALPHAR, PALPHAS, PALPHAG, PALPHAH
REAL     :: PNUR, PNUS, PNUG, PNUH
REAL     :: PBR, PBS, PBG, PBH
REAL     :: PCR, PCS, PCG, PCH
REAL     :: PDR, PDS, PDG, PDH
REAL     :: PESR, PEGS, PEGR, PEHS, PEHG
REAL     :: PACCLBDAS_MAX, PACCLBDAR_MAX, PACCLBDAS_MIN, PACCLBDAR_MIN
REAL     :: PDRYLBDAG_MAX, PDRYLBDAS_MAX, PDRYLBDAG_MIN, PDRYLBDAS_MIN
REAL     :: PDRYLBDAR_MAX, PDRYLBDAR_MIN
REAL     :: PWETLBDAS_MAX, PWETLBDAG_MAX, PWETLBDAS_MIN, PWETLBDAG_MIN
REAL     :: PWETLBDAH_MAX, PWETLBDAH_MIN
!
!-------------------------------------------------------------------------------
!
!*       0.     FUNCTION STATEMENTS
!   	        -------------------
!
!*       0.1    p_moment of the Generalized GAMMA function
!
!
!
!        1.     COMPUTE KSPLTR FOR EACH MODEL
!               -----------------------------
!
!*       1.1    Set the hailstones maximum fall velocity
!
IF (CSEDIM == 'SPLI') THEN
  IF (HCLOUD == 'ICE4') THEN
    ZVTRMAX = 40.
  ELSE IF (HCLOUD == 'ICE3') THEN
    ZVTRMAX = 10.
  END IF 
END IF
!
!*       1.2    Compute the number of small time step integration
!
KSPLITR = 1
IF (CSEDIM == 'SPLI') THEN
  SPLIT : DO
    ZT = PTSTEP / REAL(KSPLITR)
    IF (ZT * ZVTRMAX / PDZMIN .LT. 1.) EXIT SPLIT
    KSPLITR = KSPLITR + 1
  END DO SPLIT
END IF
!
IF (ALLOCATED(XRTMIN)) THEN       ! In case of nesting microphysics constants of
                                  ! MODD_RAIN_ICE_PARAM are computed only once,
                                  ! but if INI_RAIN_ICE has been called already
                                  ! one must change the XRTMIN size.
  DEALLOCATE(XRTMIN)
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.     CHARACTERISTICS OF THE SPECIES
!   	        ------------------------------
!
!*       2.1    Cloud droplet and Raindrop characteristics
!
XAC = (XPI / 6.0) * XRHOLW
XBC = 3.0
XCC = XRHOLW * XG / (18.0 * 1.7E-5) ! Stokes flow (Pruppacher p 322 for T=273K)
XDC = 2.0
!
XAR = (XPI / 6.0) * XRHOLW
XBR = 3.0
XCR = 842.
XDR = 0.8
!
XCCR = 8.E6    ! N0_r =  XCXR * lambda_r ** ZXR
ZXR  = -1.     !
!
XF0R = 1.00
XF1R = 0.26
!
XC1R = 1. / 2.
!
!
!*       2.2    Ice crystal characteristics
!
SELECT CASE (CPRISTINE_ICE)
  CASE('PLAT')
    XAI = 0.82      ! Plates
    XBI = 2.5       ! Plates 
    XC_I = 800.     ! Plates
    XDI = 1.0       ! Plates
    XC1I = 1./XPI   ! Plates
  CASE('COLU')
    XAI = 2.14E-3   ! Columns
    XBI = 1.7       ! Columns
    XC_I = 2.1E5    ! Columns
    XDI = 1.585     ! Columns
    XC1I = 0.8      ! Columns 
  CASE('BURO')
    XAI = 44.0      ! Bullet rosettes
    XBI = 3.0       ! Bullet rosettes
    XC_I = 4.3E5    ! Bullet rosettes
    XDI = 1.663     ! Bullet rosettes
    XC1I = 0.5      ! Bullet rosettes
END SELECT
!
!  Note that XCCI=N_i (a locally predicted value) and XCXI=0.0, implicitly
!
XF0I = 1.00
XF2I = 0.14
!
!
!*       2.3    Snowflakes/aggregates characteristics
!
XAS = 0.02
XBS = 1.9
XCS = 5.1
XDS = 0.27
!
XCCS = 5.0
XCXS = 1.0
!
XF0S = 0.86
XF1S = 0.28
!
XC1S = 1. / XPI
!
!
!*       2.4    Graupel/Frozen drop characteristics
!
XAG = 19.6  ! Lump graupel case
XBG = 2.8   ! Lump graupel case
XCG = 124.  ! Lump graupel case
XDG = 0.66  ! Lump graupel case
!
XCCG = 5.E5
XCXG = -0.5
! XCCG = 4.E4 ! Test of Ziegler (1988)
! XCXG = -1.0 ! Test of Ziegler (1988)
!
XF0G = 0.86
XF1G = 0.28
!
XC1G = 1. / 2.
!
!
!*       2.5    Hailstone characteristics
!
XAH = 470.
XBH = 3.0
XCH = 207.
XDH = 0.64
!
XCCH = 4.E4 ! Test of Ziegler (1988)
XCXH = -1.0 ! Test of Ziegler (1988)
!
XF0H = 0.86
XF1H = 0.28
!
XC1H = 1./2.
!
!-------------------------------------------------------------------------------
!
!*       3.     DIMENSIONAL DISTRIBUTIONS OF THE SPECIES
!	            ----------------------------------------
!
!        3.1    Cloud droplet distribution
!
! Over land
XALPHAC = 1.0  ! Gamma law of the Cloud droplet (here volume-like distribution)
XNUC    = 3.0  ! Gamma law with little dispersion
!
!
! Over sea
XALPHAC2 = 3.0  ! Gamma law of the Cloud droplet (here volume-like distribution)
XNUC2    = 1.0  ! Gamma law with little dispersion
!
!*       3.2    Raindrops distribution
!
XALPHAR = 1.0  ! Exponential law
XNUR    = 1.0  ! Exponential law
!
!*       3.3    Ice crystal distribution
!
XALPHAI = 3.0  ! Gamma law for the ice crystal volume
XNUI    = 3.0  ! Gamma law with little dispersion
!
XALPHAS = 1.0  ! Exponential law
XNUS    = 1.0  ! Exponential law
!
XALPHAG = 1.0  ! Exponential law
XNUG    = 1.0  ! Exponential law
!
XALPHAH = 1.0  ! Gamma law
XNUH    = 8.0  ! Gamma law with little dispersion
!
!*       3.4    Constants for shape parameter
!
ZGAMC  = MOMG(XALPHAC,XNUC,3.)
ZGAMC2 = MOMG(XALPHAC2,XNUC2,3.)
XLBC(1) = XAR * ZGAMC
XLBC(2) = XAR * ZGAMC2
XLBEXC = 1.0 / XBC
!
XLBEXR = 1.0 / (-1.0 - XBR)
XLBR   = (XAR * XCCR * MOMG(XALPHAR,XNUR,XBR))**(-XLBEXR)
!
XLBEXI = 1.0 / (-XBI)
XLBI   = (XAI * MOMG(XALPHAI,XNUI,XBI))**(-XLBEXI)
!
XLBEXS = 1.0 / (XCXS - XBS)
XLBS   = (XAS * XCCS * MOMG(XALPHAS,XNUS,XBS))**(-XLBEXS)
!
XLBEXG = 1.0 / (XCXG - XBG)
XLBG   = (XAG * XCCG * MOMG(XALPHAG,XNUG,XBG))**(-XLBEXG)
!
XLBEXH = 1.0/(XCXH-XBH)
XLBH   = (XAH * XCCH * MOMG(XALPHAH,XNUH,XBH))**(-XLBEXH)
!
!*       3.5    Minimal values allowed for the mixing ratios
!
XLBDAR_MAX = 100000.0
XLBDAS_MAX = 100000.0
XLBDAG_MAX = 100000.0
!
ZCONC_MAX  = 1.E6 ! Maximal concentration for falling particules set to 1 per cc
XLBDAS_MAX = (ZCONC_MAX / XCCS)**(1./XCXS)
!
IF (HCLOUD == 'ICE4') THEN
  ALLOCATE( XRTMIN(7) )
ELSE IF (HCLOUD == 'ICE3') THEN
  ALLOCATE( XRTMIN(6) )
END IF
!
XRTMIN(1) = 1.0E-20
XRTMIN(2) = 1.0E-20
XRTMIN(3) = 1.0E-20
XRTMIN(4) = 1.0E-20
XRTMIN(5) = 1.0E-15
XRTMIN(6) = 1.0E-15
IF (HCLOUD == 'ICE4') XRTMIN(7) = 1.0E-15
!
XCONC_SEA = 1.E8 ! 100/cm3
XCONC_LAND = 3.E8 ! 300/cm3
XCONC_URBAN = 5.E8 ! 500/cm3
!
!-------------------------------------------------------------------------------
!
!*       4.     CONSTANTS FOR THE SEDIMENTATION
!   	        -------------------------------
!
!*       4.1    Exponent of the fall-speed air density correction
!
XCEXVT = 0.4
!
IKB = 1 + JPVEXT
ZRHO00 = XP00 / (XRD * XTHVREFZ(IKB))
!
!*       4.2    Constants for sedimentation
!
XFSEDC(1)  = GAMMA(XNUC+(XDC+3.)/XALPHAC) / GAMMA(XNUC+3./XALPHAC) * &
            (ZRHO00)**XCEXVT
XFSEDC(2)  = GAMMA(XNUC2+(XDC+3.)/XALPHAC2) / GAMMA(XNUC2+3./XALPHAC2)* &
            (ZRHO00)**XCEXVT
!
XEXSEDR = (XBR + XDR + 1.0) / (XBR + 1.0)
XFSEDR  = XCR * XAR * XCCR * MOMG(XALPHAR,XNUR,XBR+XDR) *  &
         (XAR * XCCR * MOMG(XALPHAR,XNUR,XBR))**(-XEXSEDR) * (ZRHO00)**XCEXVT
!
XEXRSEDI = (XBI + XDI) / XBI
XEXCSEDI = 1.0 - XEXRSEDI
XFSEDI   = (4. * XPI * 900.)**(-XEXCSEDI) *                   &
           XC_I * XAI * MOMG(XALPHAI,XNUI,XBI+XDI) *          &
           ((XAI * MOMG(XALPHAI,XNUI,XBI)))**(-XEXRSEDI) *    &
           (ZRHO00)**XCEXVT
XFCI     = (4. * XPI * 900.)**(-1) 
!
!  Computations made for Columns
!
XEXRSEDI = 1.9324
XEXCSEDI =-0.9324
XFSEDI   = 3.89745E11 * MOMG(XALPHAI,XNUI,3.285) *                          &
                        MOMG(XALPHAI,XNUI,1.7)**(-XEXRSEDI)*(ZRHO00)**XCEXVT
XEXCSEDI =-0.9324 * 3.0
WRITE (KLUOUT,FMT=*)' PRISTINE ICE SEDIMENTATION for columns XFSEDI =',XFSEDI
!
XEXSEDS = (XBS + XDS - XCXS) / (XBS - XCXS)
XFSEDS  = XCS * XAS * XCCS * MOMG(XALPHAS,XNUS,XBS+XDS) *   &
         (XAS * XCCS * MOMG(XALPHAS,XNUS,XBS))**(-XEXSEDS) * (ZRHO00)**XCEXVT
!
XEXSEDG = (XBG + XDG - XCXG) / (XBG - XCXG)
XFSEDG  = XCG * XAG * XCCG * MOMG(XALPHAG,XNUG,XBG+XDG) *   &
         (XAG * XCCG * MOMG(XALPHAG,XNUG,XBG))**(-XEXSEDG) * (ZRHO00)**XCEXVT
!
XEXSEDH = (XBH + XDH - XCXH) / (XBH - XCXH)
XFSEDH  = XCH * XAH * XCCH * MOMG(XALPHAH,XNUH,XBH+XDH) *   &
         (XAH * XCCH * MOMG(XALPHAH,XNUH,XBH))**(-XEXSEDH) * (ZRHO00)**XCEXVT
!
!
!-------------------------------------------------------------------------------
!
!*       5.     CONSTANTS FOR THE SLOW COLD PROCESSES
!      	        -------------------------------------
!
!*       5.1    Constants for ice nucleation
!
SELECT CASE (CPRISTINE_ICE)
  CASE('PLAT')
    ZFACT_NUCL =  1.0    ! Plates
  CASE('COLU')
    ZFACT_NUCL = 25.0    ! Columns
  CASE('BURO')
    ZFACT_NUCL = 17.0    ! Bullet rosettes
END SELECT
!
XNU10 = 50. * ZFACT_NUCL
XALPHA1 = 4.5
XBETA1 = 0.6
!
XNU20 = 1000. * ZFACT_NUCL
XALPHA2 = 12.96
XBETA2 = 0.639
!
XMNU0 = 6.88E-13
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      Heterogeneous nucleation")')
  WRITE(UNIT=KLUOUT,FMT='(" NU10=",E13.6," ALPHA1=",E13.6," BETA1=",E13.6)') &
                                                      XNU10,XALPHA1,XBETA1
  WRITE(UNIT=KLUOUT,FMT='(" NU20=",E13.6," ALPHA2=",E13.6," BETA2=",E13.6)') &
                                                      XNU20,XALPHA2,XBETA2
  WRITE(UNIT=KLUOUT,FMT='(" mass of embryo XMNU0=",E13.6)') XMNU0
END IF
!
XALPHA3 = -3.075
XBETA3 = 81.00356
XHON = (XPI / 6.) * ((2.0 * 3.0 * 4.0 * 5.0 * 6.0) / &
                     (2.0 * 3.0)) * (1.1E5)**(-3.0)
                                       ! Pi/6 * (G_c(6)/G_c(3)) * (1/Lbda_c**3)
                                       ! avec Lbda_c=1.1E5 m^-1
                                       !     the formula is equivalent to
                                       !        rho_dref * r_c     G(6)
                                       ! Pi/6 * -------------- * ---------
                                       !         rho_lw * N_c    G(3)*G(3)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      Homogeneous nucleation")')
  WRITE(UNIT=KLUOUT,FMT='(" ALPHA3=",E13.6," BETA3=",E13.6)') XALPHA3,XBETA3
  WRITE(UNIT=KLUOUT,FMT='(" constant XHON=",E13.6)') XHON
END IF
!
!
!*       5.2    Constants for vapor deposition on ice
!
XSCFAC = (0.63**(1./3.)) * SQRT((ZRHO00)**XCEXVT) ! One assumes Sc=0.63
!
X0DEPI = (4.0 * XPI) * XC1I * XF0I * MOMG(XALPHAI,XNUI,1.)
X2DEPI = (4.0 * XPI) * XC1I * XF2I * XC_I * MOMG(XALPHAI,XNUI,XDI+2.0)
!
X0DEPS = (4.0 * XPI) * XCCS * XC1S * XF0S * MOMG(XALPHAS,XNUS,1.)
X1DEPS = (4.0 * XPI) * XCCS * XC1S * XF1S * SQRT(XCS) * MOMG(XALPHAS,XNUS,0.5*XDS+1.5)
XEX0DEPS = XCXS - 1.0
XEX1DEPS = XCXS - 0.5 * (XDS + 3.0)
!
X0DEPG = (4.0 * XPI) * XCCG * XC1G * XF0G * MOMG(XALPHAG,XNUG,1.)
X1DEPG = (4.0 * XPI) * XCCG * XC1G * XF1G * SQRT(XCG) * MOMG(XALPHAG,XNUG,0.5*XDG+1.5)
XEX0DEPG = XCXG - 1.0
XEX1DEPG = XCXG - 0.5 * (XDG + 3.0)
!
X0DEPH = (4.0 * XPI) * XCCH * XC1H * XF0H * MOMG(XALPHAH,XNUH,1.)
X1DEPH = (4.0 * XPI) * XCCH * XC1H * XF1H * SQRT(XCH) * MOMG(XALPHAH,XNUH,0.5*XDH+1.5)
XEX0DEPH = XCXH - 1.0
XEX1DEPH = XCXH - 0.5 * (XDH + 3.0)
!
!
!*       5.3    Constants for pristine ice autoconversion
!
XTIMAUTI = 1.E-3  !  Time constant at T=T_t
XTEXAUTI = 0.015  !  Temperature factor of the I+I collection efficiency
XCRIAUTI = 0.2E-4 !  Critical ice content for the autoconversion to occur
                  !  Revised value by Chaboureau et al. (2001)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      pristine ice autoconversion")')
  WRITE(UNIT=KLUOUT,FMT='(" Time constant   XTIMAUTI=",E13.6)') XTIMAUTI
  WRITE(UNIT=KLUOUT,FMT='(" Temp. factor    XTEXAUTI=",E13.6)') XTEXAUTI
  WRITE(UNIT=KLUOUT,FMT='(" Crit. ice cont. XCRIAUTI=",E13.6)') XCRIAUTI
END IF
!
!
!*       5.4    Constants for snow aggregation
!
XCOLIS   = 0.25 ! Collection efficiency of I+S
XCOLEXIS = 0.05 ! Temperature factor of the I+S collection efficiency
XFIAGGS  = (XPI / 4.0) * XCOLIS * XCCS * XCS * (ZRHO00**XCEXVT) * &
            MOMG(XALPHAS,XNUS,XDS+2.0)
XEXIAGGS = XCXS - XDS - 2.0
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      snow aggregation")')
  WRITE(UNIT=KLUOUT,FMT='(" Coll. efficiency XCOLIS=",E13.6)') XCOLIS
  WRITE(UNIT=KLUOUT,FMT='(" Temp. factor     XCOLEXIS=",E13.6)') XCOLEXIS
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       6.     CONSTANTS FOR THE SLOW WARM PROCESSES
!   	        -------------------------------------
!
!*       6.1    Constants for the cloud droplets autoconversion
!
XTIMAUTC = 1.E-3
XCRIAUTC = 0.5E-3 
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      cloud droplets autoconversion")')
  WRITE(UNIT=KLUOUT,FMT='(" Time constant   XTIMAUTC=",E13.6)') XTIMAUTC
  WRITE(UNIT=KLUOUT,FMT='(" Crit. ice cont. XCRIAUTC=",E13.6)') XCRIAUTC
END IF
!
!*       6.2    Constants for the accretion of cloud droplets by raindrops
!
XFCACCR  = (XPI / 4.0) * XCCR * XCR * (ZRHO00**XCEXVT) * MOMG(XALPHAR,XNUR,XDR+2.0)
XEXCACCR = -XDR - 3.0
!
!*       6.3    Constants for the evaporation of the raindrops
!
X0EVAR = (4.0 * XPI) * XCCR * XC1R * XF0R * MOMG(XALPHAR,XNUR,1.)
X1EVAR = (4.0 * XPI) * XCCR * XC1R * XF1R * SQRT(XCR)*MOMG(XALPHAR,XNUR,0.5*XDR+1.5)
XEX0EVAR = -2.0
XEX1EVAR = -1.0 - 0.5 * (XDR + 3.0)
!
!
!-------------------------------------------------------------------------------
!
!*       7.     CONSTANTS FOR THE FAST COLD PROCESSES FOR THE AGGREGATES
!	        --------------------------------------------------------
!
!
!*       7.1    Constants for the riming of the aggregates
!
XDCSLIM  = 0.007 ! D_cs^lim = 7 mm as suggested by Farley et al. (1989)
XCOLCS   = 1.0
XEXCRIMSS= XCXS - XDS - 2.0
XCRIMSS  = (XPI / 4.0) * XCOLCS * XCCS * XCS * (ZRHO00**XCEXVT) * MOMG(XALPHAS,XNUS,XDS+2.0)
XEXCRIMSG= XEXCRIMSS
XCRIMSG  = XCRIMSS
XSRIMCG  = XCCS * XAS * MOMG(XALPHAS,XNUS,XBS)
XEXSRIMCG= XCXS - XBS
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      riming of the aggregates")')
  WRITE(UNIT=KLUOUT,FMT='(" D_cs^lim (Farley et al.) XDCSLIM=",E13.6)') XDCSLIM
  WRITE(UNIT=KLUOUT,FMT='(" Coll. efficiency          XCOLCS=",E13.6)') XCOLCS
END IF
!
NGAMINC = 80
XGAMINC_BOUND_MIN = 1.0E-1 ! Minimal value of (Lbda * D_cs^lim)**alpha
XGAMINC_BOUND_MAX = 1.0E7  ! Maximal value of (Lbda * D_cs^lim)**alpha
ZRATE = EXP(LOG(XGAMINC_BOUND_MAX/XGAMINC_BOUND_MIN)/REAL(NGAMINC-1))
!
IF( .NOT.ALLOCATED(XGAMINC_RIM1) ) ALLOCATE( XGAMINC_RIM1(NGAMINC) )
IF( .NOT.ALLOCATED(XGAMINC_RIM2) ) ALLOCATE( XGAMINC_RIM2(NGAMINC) )
IF( .NOT.ALLOCATED(XGAMINC_RIM3) ) ALLOCATE( XGAMINC_RIM3(NGAMINC) )
!
DO J1 = 1, NGAMINC
  ZBOUND = XGAMINC_BOUND_MIN * ZRATE**(J1-1)
  XGAMINC_RIM1(J1) = GAMMA_INC(XNUS+(2.0+XDS)/XALPHAS,ZBOUND)
  XGAMINC_RIM2(J1) = GAMMA_INC(XNUS+XBS/XALPHAS      ,ZBOUND)
  XFS = 1.3  ! cf values initiated in ini_param_elec
  XGAMINC_RIM3(J1) = GAMMA_INC(XNUS+XFS/XALPHAS      ,ZBOUND)
END DO
!
XRIMINTP1 = XALPHAS / LOG(ZRATE)
XRIMINTP2 = 1.0 + XRIMINTP1 * LOG(XDCSLIM/(XGAMINC_BOUND_MIN)**(1.0/XALPHAS))
!
!*       7.2    Constants for the accretion of raindrops onto aggregates
!
XFRACCSS = ((XPI**2) / 24.0) * XCCS * XCCR * XRHOLW * (ZRHO00**XCEXVT)
!
XLBRACCS1   =      MOMG(XALPHAS,XNUS,2.) * MOMG(XALPHAR,XNUR,3.)
XLBRACCS2   = 2. * MOMG(XALPHAS,XNUS,1.) * MOMG(XALPHAR,XNUR,4.)
XLBRACCS3   =                              MOMG(XALPHAR,XNUR,5.)
!
XFSACCRG = (XPI / 4.0) * XAS * XCCS * XCCR * (ZRHO00**XCEXVT)
!
XLBSACCR1   =      MOMG(XALPHAR,XNUR,2.) * MOMG(XALPHAS,XNUS,XBS)
XLBSACCR2   = 2. * MOMG(XALPHAR,XNUR,1.) * MOMG(XALPHAS,XNUS,XBS+1.)
XLBSACCR3   =                              MOMG(XALPHAS,XNUS,XBS+2.)
!
!*       7.2.1  Defining the ranges for the computation of the kernels
!
! Notice: One magnitude of lambda discretized over 10 points for rain
! Notice: One magnitude of lambda discretized over 10 points for snow
!
NACCLBDAS = 40
XACCLBDAS_MIN = 5.0E1 ! Minimal value of Lbda_s to tabulate XKER_RACCS
XACCLBDAS_MAX = 5.0E5 ! Maximal value of Lbda_s to tabulate XKER_RACCS
ZRATE = LOG(XACCLBDAS_MAX/XACCLBDAS_MIN)/REAL(NACCLBDAS-1)
XACCINTP1S = 1.0 / ZRATE
XACCINTP2S = 1.0 - LOG( XACCLBDAS_MIN ) / ZRATE
!
NACCLBDAR = 40
XACCLBDAR_MIN = 1.0E3 ! Minimal value of Lbda_r to tabulate XKER_RACCS
XACCLBDAR_MAX = 1.0E7 ! Maximal value of Lbda_r to tabulate XKER_RACCS
ZRATE = LOG(XACCLBDAR_MAX/XACCLBDAR_MIN)/REAL(NACCLBDAR-1)
XACCINTP1R = 1.0 / ZRATE
XACCINTP2R = 1.0 - LOG( XACCLBDAR_MIN ) / ZRATE
!
!*       7.2.2  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZESR     = 1.0   ! factor used to integrate the dimensional distributions when
ZFDINFTY = 20.0  ! computing the kernels XKER_RACCSS, XKER_RACCS and XKER_SACCRG
!
IF( .NOT.ALLOCATED(XKER_RACCSS) ) ALLOCATE( XKER_RACCSS(NACCLBDAS,NACCLBDAR) )
IF( .NOT.ALLOCATED(XKER_RACCS ) ) ALLOCATE( XKER_RACCS (NACCLBDAS,NACCLBDAR) )
IF( .NOT.ALLOCATED(XKER_SACCRG) ) ALLOCATE( XKER_SACCRG(NACCLBDAR,NACCLBDAS) )
!
CALL READ_XKER_RACCS (KACCLBDAS,KACCLBDAR,KND,                                &
                      PALPHAS,PNUS,PALPHAR,PNUR,PESR,PBS,PBR,PCS,PDS,PCR,PDR, &
                      PACCLBDAS_MAX,PACCLBDAR_MAX,PACCLBDAS_MIN,PACCLBDAR_MIN,&
                      PFDINFTY                                                )
IF( (KACCLBDAS/=NACCLBDAS) .OR. (KACCLBDAR/=NACCLBDAR) .OR. (KND/=IND) .OR. &
    (PALPHAS/=XALPHAS) .OR. (PNUS/=XNUS)                               .OR. &
    (PALPHAR/=XALPHAR) .OR. (PNUR/=XNUR)                               .OR. &
    (PESR/=ZESR) .OR. (PBS/=XBS) .OR. (PBR/=XBR)                       .OR. &
    (PCS/=XCS) .OR. (PDS/=XDS) .OR. (PCR/=XCR) .OR. (PDR/=XDR)         .OR. &
    (PACCLBDAS_MAX/=XACCLBDAS_MAX) .OR. (PACCLBDAR_MAX/=XACCLBDAR_MAX) .OR. &
    (PACCLBDAS_MIN/=XACCLBDAS_MIN) .OR. (PACCLBDAR_MIN/=XACCLBDAR_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RRCOLSS ( IND, XALPHAS, XNUS, XALPHAR, XNUR,                          &
                 ZESR, XBR, XCS, XDS, XCR, XDR,                              &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_RACCSS, XAG, XBS, XAS                        )
  CALL RZCOLX  ( IND, XALPHAS, XNUS, XALPHAR, XNUR,                          &
                 ZESR, XBR, XCS, XDS, XCR, XDR,                              &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_RACCS                                        )
  CALL RSCOLRG ( IND, XALPHAS, XNUS, XALPHAR, XNUR,                          &
                 ZESR, XBS, XCS, XDS, XCR, XDR,                              &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_SACCRG,  XAG, XBS, XAS                       )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF RACSS KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF RACS  KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF SACRG KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KACCLBDAS=",I3)') NACCLBDAS
  WRITE(UNIT=KLUOUT,FMT='("KACCLBDAR=",I3)') NACCLBDAR
  WRITE(UNIT=KLUOUT,FMT='("PALPHAS=",E13.6)') XALPHAS
  WRITE(UNIT=KLUOUT,FMT='("PNUS=",E13.6)') XNUS
  WRITE(UNIT=KLUOUT,FMT='("PALPHAR=",E13.6)') XALPHAR
  WRITE(UNIT=KLUOUT,FMT='("PNUR=",E13.6)') XNUR
  WRITE(UNIT=KLUOUT,FMT='("PESR=",E13.6)') ZESR
  WRITE(UNIT=KLUOUT,FMT='("PBS=",E13.6)') XBS
  WRITE(UNIT=KLUOUT,FMT='("PBR=",E13.6)') XBR
  WRITE(UNIT=KLUOUT,FMT='("PCS=",E13.6)') XCS
  WRITE(UNIT=KLUOUT,FMT='("PDS=",E13.6)') XDS
  WRITE(UNIT=KLUOUT,FMT='("PCR=",E13.6)') XCR
  WRITE(UNIT=KLUOUT,FMT='("PDR=",E13.6)') XDR
  WRITE(UNIT=KLUOUT,FMT='("PACCLBDAS_MAX=",E13.6)') &
                                                    XACCLBDAS_MAX
  WRITE(UNIT=KLUOUT,FMT='("PACCLBDAR_MAX=",E13.6)') &
                                                    XACCLBDAR_MAX
  WRITE(UNIT=KLUOUT,FMT='("PACCLBDAS_MIN=",E13.6)') &
                                                    XACCLBDAS_MIN
  WRITE(UNIT=KLUOUT,FMT='("PACCLBDAR_MIN=",E13.6)') &
                                                    XACCLBDAR_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_RACCSS) ) THEN")')
  DO J1 = 1 , NACCLBDAS
    DO J2 = 1 , NACCLBDAR
    WRITE(UNIT=KLUOUT,FMT='("  PKER_RACCSS(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RACCSS(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_RACCS ) ) THEN")')
  DO J1 = 1 , NACCLBDAS
    DO J2 = 1 , NACCLBDAR
    WRITE(UNIT=KLUOUT,FMT='("  PKER_RACCS (",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RACCS (J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_SACCRG) ) THEN")')
  DO J1 = 1 , NACCLBDAR
    DO J2 = 1 , NACCLBDAS
    WRITE(UNIT=KLUOUT,FMT='("  PKER_SACCRG(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_SACCRG(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_RACCS (KACCLBDAS,KACCLBDAR,KND,                               &
                       PALPHAS,PNUS,PALPHAR,PNUR,PESR,PBS,PBR,PCS,PDS,PCR,PDR, &
                       PACCLBDAS_MAX,PACCLBDAR_MAX,PACCLBDAS_MIN,PACCLBDAR_MIN,&
                       PFDINFTY,XKER_RACCSS,XKER_RACCS,XKER_SACCRG             )
  WRITE(UNIT=KLUOUT,FMT='(" Read XKER_RACCSS")')
  WRITE(UNIT=KLUOUT,FMT='(" Read XKER_RACCS ")')
  WRITE(UNIT=KLUOUT,FMT='(" Read XKER_SACCRG")')
END IF
!
!*       7.3    Constant for the conversion-melting rate
!
XFSCVMG = 2.0
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      conversion-melting of the aggregates")')
  WRITE(UNIT=KLUOUT,FMT='(" Conv. factor XFSCVMG=",E13.6)') XFSCVMG
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       8.     CONSTANTS FOR THE FAST COLD PROCESSES FOR THE GRAUPELN
!               ------------------------------------------------------
!
!
!*       8.1    Constants for the rain contact freezing
!
XCOLIR    = 1.0
!
XEXRCFRI  = -XDR - 5.0 + ZXR
XRCFRI    = ((XPI**2) / 24.0) * XCCR * XRHOLW * XCOLIR * XCR * &
             (ZRHO00**XCEXVT) * MOMG(XALPHAR,XNUR,XDR+5.0)
XEXICFRR  = -XDR - 2.0 + ZXR
XICFRR    = (XPI / 4.0) * XCOLIR * XCR * (ZRHO00**XCEXVT) * &
             XCCR * MOMG(XALPHAR,XNUR,XDR+2.0)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      rain contact freezing")')
  WRITE(UNIT=KLUOUT,FMT='(" Coll. efficiency          XCOLIR=",E13.6)') XCOLIR
END IF
!
!
!*       8.2    Constants for the dry growth of the graupeln
!
!*       8.2.1  Constants for the cloud droplet collection by the graupeln
!
XFCDRYG = (XPI / 4.0) * XCCG * XCG * (ZRHO00**XCEXVT) * MOMG(XALPHAG,XNUG,XDG+2.0)
!
!*       8.2.2  Constants for the cloud ice collection by the graupeln
!
XCOLIG   = 0.25 ! Collection efficiency of I+G
XCOLEXIG = 0.05 ! Temperature factor of the I+G collection efficiency
XCOLIG   = 0.01 ! Collection efficiency of I+G
XCOLEXIG = 0.1  ! Temperature factor of the I+G collection efficiency
WRITE (KLUOUT, FMT=*) ' NEW Constants for the cloud ice collection by the graupeln'
WRITE (KLUOUT, FMT=*) ' XCOLIG, XCOLEXIG  = ',XCOLIG,XCOLEXIG
!
XFIDRYG = (XPI / 4.0) * XCOLIG * XCCG * XCG * (ZRHO00**XCEXVT) * &
           MOMG(XALPHAG,XNUG,XDG+2.0)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      cloud ice collection by the graupeln")')
  WRITE(UNIT=KLUOUT,FMT='(" Coll. efficiency XCOLIG=",E13.6)') XCOLIG
  WRITE(UNIT=KLUOUT,FMT='(" Temp. factor     XCOLEXIG=",E13.6)') XCOLEXIG
END IF
!
!*       8.2.3  Constants for the aggregate collection by the graupeln
!
XCOLSG   = 0.25 ! Collection efficiency of S+G
XCOLEXSG = 0.05 ! Temperature factor of the S+G collection efficiency
XCOLSG   = 0.01 ! Collection efficiency of S+G
XCOLEXSG = 0.1  ! Temperature factor of the S+G collection efficiency
WRITE (KLUOUT, FMT=*) ' NEW Constants for the aggregate collection by the graupeln'
WRITE (KLUOUT, FMT=*) ' XCOLSG, XCOLEXSG  = ',XCOLSG,XCOLEXSG
!
XFSDRYG = (XPI / 4.0) * XCOLSG * XCCG * XCCS * XAS * (ZRHO00**XCEXVT)
!
XLBSDRYG1   =      MOMG(XALPHAG,XNUG,2.) * MOMG(XALPHAS,XNUS,XBS)
XLBSDRYG2   = 2. * MOMG(XALPHAG,XNUG,1.) * MOMG(XALPHAS,XNUS,XBS+1.)
XLBSDRYG3   =                              MOMG(XALPHAS,XNUS,XBS+2.)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      aggregate collection by the graupeln")')
  WRITE(UNIT=KLUOUT,FMT='(" Coll. efficiency XCOLSG=",E13.6)') XCOLSG
  WRITE(UNIT=KLUOUT,FMT='(" Temp. factor     XCOLEXSG=",E13.6)') XCOLEXSG
END IF
!
!*       8.2.4  Constants for the raindrop collection by the graupeln
!
XFRDRYG = ((XPI**2) / 24.0) * XCCG * XCCR * XRHOLW * (ZRHO00**XCEXVT)
!
XLBRDRYG1   =      MOMG(XALPHAG,XNUG,2.) * MOMG(XALPHAR,XNUR,3.)
XLBRDRYG2   = 2. * MOMG(XALPHAG,XNUG,1.) * MOMG(XALPHAR,XNUR,4.)
XLBRDRYG3   =                              MOMG(XALPHAR,XNUR,5.)
!
! Notice: One magnitude of lambda discretized over 10 points
!
NDRYLBDAR = 40
XDRYLBDAR_MIN = 1.0E3 ! Minimal value of Lbda_r to tabulate XKER_RDRYG
XDRYLBDAR_MAX = 1.0E7 ! Maximal value of Lbda_r to tabulate XKER_RDRYG
ZRATE = LOG(XDRYLBDAR_MAX/XDRYLBDAR_MIN) / REAL(NDRYLBDAR-1)
XDRYINTP1R = 1.0 / ZRATE
XDRYINTP2R = 1.0 - LOG( XDRYLBDAR_MIN ) / ZRATE
!
NDRYLBDAS = 80
XDRYLBDAS_MIN = 2.5E1 ! Minimal value of Lbda_s to tabulate XKER_SDRYG
XDRYLBDAS_MAX = 2.5E9 ! Maximal value of Lbda_s to tabulate XKER_SDRYG
ZRATE = LOG(XDRYLBDAS_MAX/XDRYLBDAS_MIN) / REAL(NDRYLBDAS-1)
XDRYINTP1S = 1.0 / ZRATE
XDRYINTP2S = 1.0 - LOG( XDRYLBDAS_MIN ) / ZRATE
!
NDRYLBDAG = 40
XDRYLBDAG_MIN = 1.0E3 ! Min value of Lbda_g to tabulate XKER_SDRYG,XKER_RDRYG
XDRYLBDAG_MAX = 1.0E7 ! Max value of Lbda_g to tabulate XKER_SDRYG,XKER_RDRYG
ZRATE = LOG(XDRYLBDAG_MAX/XDRYLBDAG_MIN) / REAL(NDRYLBDAG-1)
XDRYINTP1G = 1.0 / ZRATE
XDRYINTP2G = 1.0 - LOG( XDRYLBDAG_MIN ) / ZRATE
!
!*       8.2.5  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZEGS     = 1.0   ! factor used to integrate the dimensional distributions when
ZFDINFTY = 20.0  ! computing the kernels XKER_SDRYG
!
IF( .NOT.ALLOCATED(XKER_SDRYG) ) ALLOCATE( XKER_SDRYG(NDRYLBDAG,NDRYLBDAS) )
!
CALL READ_XKER_SDRYG (KDRYLBDAG,KDRYLBDAS,KND,                              &
                   PALPHAG,PNUG,PALPHAS,PNUS,PEGS,PBS,PCG,PDG,PCS,PDS,      &
                   PDRYLBDAG_MAX,PDRYLBDAS_MAX,PDRYLBDAG_MIN,PDRYLBDAS_MIN, &
                   PFDINFTY                                                 )
IF( (KDRYLBDAG/=NDRYLBDAG) .OR. (KDRYLBDAS/=NDRYLBDAS) .OR. (KND/=IND) .OR. &
    (PALPHAG/=XALPHAG) .OR. (PNUG/=XNUG)                               .OR. &
    (PALPHAS/=XALPHAS) .OR. (PNUS/=XNUS)                               .OR. &
    (PEGS/=ZEGS) .OR. (PBS/=XBS)                                       .OR. &
    (PCG/=XCG) .OR. (PDG/=XDG) .OR. (PCS/=XCS) .OR. (PDS/=XDS)         .OR. &
    (PDRYLBDAG_MAX/=XDRYLBDAG_MAX) .OR. (PDRYLBDAS_MAX/=XDRYLBDAS_MAX) .OR. &
    (PDRYLBDAG_MIN/=XDRYLBDAG_MIN) .OR. (PDRYLBDAS_MIN/=XDRYLBDAS_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAG, XNUG, XALPHAS, XNUS,                          &
                ZEGS, XBS, XCG, XDG, XCS, XDS,                              &
                XDRYLBDAG_MAX, XDRYLBDAS_MAX, XDRYLBDAG_MIN, XDRYLBDAS_MIN, &
                ZFDINFTY, XKER_SDRYG                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF SDRYG KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KDRYLBDAG=",I3)') NDRYLBDAG
  WRITE(UNIT=KLUOUT,FMT='("KDRYLBDAS=",I3)') NDRYLBDAS
  WRITE(UNIT=KLUOUT,FMT='("PALPHAG=",E13.6)') XALPHAG
  WRITE(UNIT=KLUOUT,FMT='("PNUG=",E13.6)') XNUG
  WRITE(UNIT=KLUOUT,FMT='("PALPHAS=",E13.6)') XALPHAS
  WRITE(UNIT=KLUOUT,FMT='("PNUS=",E13.6)') XNUS
  WRITE(UNIT=KLUOUT,FMT='("PEGS=",E13.6)') ZEGS
  WRITE(UNIT=KLUOUT,FMT='("PBS=",E13.6)') XBS
  WRITE(UNIT=KLUOUT,FMT='("PCG=",E13.6)') XCG
  WRITE(UNIT=KLUOUT,FMT='("PDG=",E13.6)') XDG
  WRITE(UNIT=KLUOUT,FMT='("PCS=",E13.6)') XCS
  WRITE(UNIT=KLUOUT,FMT='("PDS=",E13.6)') XDS
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAG_MAX=",E13.6)') &
                                                    XDRYLBDAG_MAX
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAS_MAX=",E13.6)') &
                                                    XDRYLBDAS_MAX
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAG_MIN=",E13.6)') &
                                                    XDRYLBDAG_MIN
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAS_MIN=",E13.6)') &
                                                    XDRYLBDAS_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_SDRYG) ) THEN")')
  DO J1 = 1 , NDRYLBDAG
    DO J2 = 1 , NDRYLBDAS
    WRITE(UNIT=KLUOUT,FMT='("PKER_SDRYG(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_SDRYG(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_SDRYG (KDRYLBDAG,KDRYLBDAS,KND,                              &
                     PALPHAG,PNUG,PALPHAS,PNUS,PEGS,PBS,PCG,PDG,PCS,PDS,      &
                     PDRYLBDAG_MAX,PDRYLBDAS_MAX,PDRYLBDAG_MIN,PDRYLBDAS_MIN, &
                     PFDINFTY,XKER_SDRYG                                      )
  WRITE(UNIT=KLUOUT,FMT='(" Read XKER_SDRYG")')
END IF
!
!
IND      = 50    ! Number of interval used to integrate the dimensional
ZEGR     = 1.0   ! distributions when computing the kernel XKER_RDRYG
ZFDINFTY = 20.0
!
IF( .NOT.ALLOCATED(XKER_RDRYG) ) ALLOCATE( XKER_RDRYG(NDRYLBDAG,NDRYLBDAR) )
!
CALL READ_XKER_RDRYG (KDRYLBDAG,KDRYLBDAR,KND,                              &
                   PALPHAG,PNUG,PALPHAR,PNUR,PEGR,PBR,PCG,PDG,PCR,PDR,      &
                   PDRYLBDAG_MAX,PDRYLBDAR_MAX,PDRYLBDAG_MIN,PDRYLBDAR_MIN, &
                   PFDINFTY                                                 )
IF( (KDRYLBDAG/=NDRYLBDAG) .OR. (KDRYLBDAR/=NDRYLBDAR) .OR. (KND/=IND) .OR. &
    (PALPHAG/=XALPHAG) .OR. (PNUG/=XNUG)                               .OR. &
    (PALPHAR/=XALPHAR) .OR. (PNUR/=XNUR)                               .OR. &
    (PEGR/=ZEGR) .OR. (PBR/=XBR)                                       .OR. &
    (PCG/=XCG) .OR. (PDG/=XDG) .OR. (PCR/=XCR) .OR. (PDR/=XDR)         .OR. &
    (PDRYLBDAG_MAX/=XDRYLBDAG_MAX) .OR. (PDRYLBDAR_MAX/=XDRYLBDAR_MAX) .OR. &
    (PDRYLBDAG_MIN/=XDRYLBDAG_MIN) .OR. (PDRYLBDAR_MIN/=XDRYLBDAR_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAG, XNUG, XALPHAR, XNUR,                          &
                ZEGR, XBR, XCG, XDG, XCR, XDR,                              &
                XDRYLBDAG_MAX, XDRYLBDAR_MAX, XDRYLBDAG_MIN, XDRYLBDAR_MIN, &
                ZFDINFTY, XKER_RDRYG                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF RDRYG KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KDRYLBDAG=",I3)') NDRYLBDAG
  WRITE(UNIT=KLUOUT,FMT='("KDRYLBDAR=",I3)') NDRYLBDAR
  WRITE(UNIT=KLUOUT,FMT='("PALPHAG=",E13.6)') XALPHAG
  WRITE(UNIT=KLUOUT,FMT='("PNUG=",E13.6)') XNUG
  WRITE(UNIT=KLUOUT,FMT='("PALPHAR=",E13.6)') XALPHAR
  WRITE(UNIT=KLUOUT,FMT='("PNUR=",E13.6)') XNUR
  WRITE(UNIT=KLUOUT,FMT='("PEGR=",E13.6)') ZEGR
  WRITE(UNIT=KLUOUT,FMT='("PBR=",E13.6)') XBR
  WRITE(UNIT=KLUOUT,FMT='("PCG=",E13.6)') XCG
  WRITE(UNIT=KLUOUT,FMT='("PDG=",E13.6)') XDG
  WRITE(UNIT=KLUOUT,FMT='("PCR=",E13.6)') XCR
  WRITE(UNIT=KLUOUT,FMT='("PDR=",E13.6)') XDR
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAG_MAX=",E13.6)') &
                                                    XDRYLBDAG_MAX
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAR_MAX=",E13.6)') &
                                                    XDRYLBDAR_MAX
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAG_MIN=",E13.6)') &
                                                    XDRYLBDAG_MIN
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAR_MIN=",E13.6)') &
                                                    XDRYLBDAR_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_RDRYG) ) THEN")')
  DO J1 = 1 , NDRYLBDAG
    DO J2 = 1 , NDRYLBDAR
    WRITE(UNIT=KLUOUT,FMT='("PKER_RDRYG(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RDRYG(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_RDRYG (KDRYLBDAG,KDRYLBDAR,KND,                              &
                     PALPHAG,PNUG,PALPHAR,PNUR,PEGR,PBR,PCG,PDG,PCR,PDR,      &
                     PDRYLBDAG_MAX,PDRYLBDAR_MAX,PDRYLBDAG_MIN,PDRYLBDAR_MIN, &
                     PFDINFTY,XKER_RDRYG                                      )
  WRITE(UNIT=KLUOUT,FMT='(" Read XKER_RDRYG")')
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       9.     CONSTANTS FOR THE FAST COLD PROCESSES FOR THE HAILSTONES
!               --------------------------------------------------------
!
!*       9.2    Constants for the wet growth of the hailstones
!
!
!*       9.2.1  Constant for the cloud droplet and cloud ice collection
!               by the hailstones
!
XFWETH = (XPI / 4.0) * XCCH * XCH * (ZRHO00**XCEXVT) * MOMG(XALPHAH,XNUH,XDH+2.0)
!
!*       9.2.2  Constants for the aggregate collection by the hailstones
!
XFSWETH = (XPI/4.0) * XCCH * XCCS * XAS * (ZRHO00**XCEXVT)
!
XLBSWETH1   =      MOMG(XALPHAH,XNUH,2.) * MOMG(XALPHAS,XNUS,XBS)
XLBSWETH2   = 2. * MOMG(XALPHAH,XNUH,1.) * MOMG(XALPHAS,XNUS,XBS+1.)
XLBSWETH3   =                              MOMG(XALPHAS,XNUS,XBS+2.)
!
!*       9.2.3  Constants for the graupel collection by the hailstones
!
XFGWETH = (XPI / 4.0) * XCCH * XCCG * XAG * (ZRHO00**XCEXVT)
!
XLBGWETH1   =      MOMG(XALPHAH,XNUH,2.) * MOMG(XALPHAG,XNUG,XBG)
XLBGWETH2   = 2. * MOMG(XALPHAH,XNUH,1.) * MOMG(XALPHAG,XNUG,XBG+1.)
XLBGWETH3   =                              MOMG(XALPHAG,XNUG,XBG+2.)
!
! Notice: One magnitude of lambda discretized over 10 points
!
NWETLBDAS = 80
XWETLBDAS_MIN = 2.5E1 ! Minimal value of Lbda_s to tabulate XKER_SWETH
XWETLBDAS_MAX = 2.5E9 ! Maximal value of Lbda_s to tabulate XKER_SWETH
ZRATE = LOG(XWETLBDAS_MAX/XWETLBDAS_MIN) / REAL(NWETLBDAS-1)
XWETINTP1S = 1.0 / ZRATE
XWETINTP2S = 1.0 - LOG( XWETLBDAS_MIN ) / ZRATE
NWETLBDAG = 40
XWETLBDAG_MIN = 1.0E3 ! Min value of Lbda_g to tabulate XKER_GWETH
XWETLBDAG_MAX = 1.0E7 ! Max value of Lbda_g to tabulate XKER_GWETH
ZRATE = LOG(XWETLBDAG_MAX/XWETLBDAG_MIN) / REAL(NWETLBDAG-1)
XWETINTP1G = 1.0 / ZRATE
XWETINTP2G = 1.0 - LOG( XWETLBDAG_MIN ) / ZRATE
NWETLBDAH = 40
XWETLBDAH_MIN = 1.0E3 ! Min value of Lbda_h to tabulate XKER_SWETH,XKER_GWETH
XWETLBDAH_MAX = 1.0E7 ! Max value of Lbda_h to tabulate XKER_SWETH,XKER_GWETH
ZRATE = LOG(XWETLBDAH_MAX/XWETLBDAH_MIN) / REAL(NWETLBDAH-1)
XWETINTP1H = 1.0 / ZRATE
XWETINTP2H = 1.0 - LOG( XWETLBDAH_MIN ) / ZRATE
!
!*       9.2.4  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZEHS     = 1.0   ! factor used to integrate the dimensional distributions when
ZFDINFTY = 20.0  ! computing the kernels XKER_SWETH
!
IF( .NOT.ALLOCATED(XKER_SWETH) ) ALLOCATE( XKER_SWETH(NWETLBDAH,NWETLBDAS) )
!
CALL READ_XKER_SWETH (KWETLBDAH,KWETLBDAS,KND,                              &
                   PALPHAH,PNUH,PALPHAS,PNUS,PEHS,PBS,PCH,PDH,PCS,PDS,      &
                   PWETLBDAH_MAX,PWETLBDAS_MAX,PWETLBDAH_MIN,PWETLBDAS_MIN, &
                   PFDINFTY                                                 )
IF( (KWETLBDAH/=NWETLBDAH) .OR. (KWETLBDAS/=NWETLBDAS) .OR. (KND/=IND) .OR. &
    (PALPHAH/=XALPHAH) .OR. (PNUH/=XNUH)                               .OR. &
    (PALPHAS/=XALPHAS) .OR. (PNUS/=XNUS)                               .OR. &
    (PEHS/=ZEHS) .OR. (PBS/=XBS)                                       .OR. &
    (PCH/=XCH) .OR. (PDH/=XDH) .OR. (PCS/=XCS) .OR. (PDS/=XDS)         .OR. &
    (PWETLBDAH_MAX/=XWETLBDAH_MAX) .OR. (PWETLBDAS_MAX/=XWETLBDAS_MAX) .OR. &
    (PWETLBDAH_MIN/=XWETLBDAH_MIN) .OR. (PWETLBDAS_MIN/=XWETLBDAS_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAH, XNUH, XALPHAS, XNUS,                          &
                ZEHS, XBS, XCH, XDH, XCS, XDS,                              &
                XWETLBDAH_MAX, XWETLBDAS_MAX, XWETLBDAH_MIN, XWETLBDAS_MIN, &
                ZFDINFTY, XKER_SWETH                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF SWETH KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAH=",I3)') NWETLBDAH
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAS=",I3)') NWETLBDAS
  WRITE(UNIT=KLUOUT,FMT='("PALPHAH=",E13.6)') XALPHAH
  WRITE(UNIT=KLUOUT,FMT='("PNUH=",E13.6)') XNUH
  WRITE(UNIT=KLUOUT,FMT='("PALPHAS=",E13.6)') XALPHAS
  WRITE(UNIT=KLUOUT,FMT='("PNUS=",E13.6)') XNUS
  WRITE(UNIT=KLUOUT,FMT='("PEHS=",E13.6)') ZEHS
  WRITE(UNIT=KLUOUT,FMT='("PBS=",E13.6)') XBS
  WRITE(UNIT=KLUOUT,FMT='("PCH=",E13.6)') XCH
  WRITE(UNIT=KLUOUT,FMT='("PDH=",E13.6)') XDH
  WRITE(UNIT=KLUOUT,FMT='("PCS=",E13.6)') XCS
  WRITE(UNIT=KLUOUT,FMT='("PDS=",E13.6)') XDS
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MAX=",E13.6)') &
                                                    XWETLBDAH_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAS_MAX=",E13.6)') &
                                                    XWETLBDAS_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MIN=",E13.6)') &
                                                    XWETLBDAH_MIN
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAS_MIN=",E13.6)') &
                                                    XWETLBDAS_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_SWETH) ) THEN")')
  DO J1 = 1 , NWETLBDAH
    DO J2 = 1 , NWETLBDAS
    WRITE(UNIT=KLUOUT,FMT='("PKER_SWETH(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_SWETH(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_SWETH (KWETLBDAH,KWETLBDAS,KND,                              &
                     PALPHAH,PNUH,PALPHAS,PNUS,PEHS,PBS,PCH,PDH,PCS,PDS,      &
                     PWETLBDAH_MAX,PWETLBDAS_MAX,PWETLBDAH_MIN,PWETLBDAS_MIN, &
                     PFDINFTY,XKER_SWETH                                      )
  WRITE(UNIT=KLUOUT,FMT='(" Read XKER_SWETH")')
END IF
!
!
IND      = 50    ! Number of interval used to integrate the dimensional
ZEHG     = 1.0   ! distributions when computing the kernel XKER_GWETH
ZFDINFTY = 20.0
!
IF( .NOT.ALLOCATED(XKER_GWETH) ) ALLOCATE( XKER_GWETH(NWETLBDAH,NWETLBDAG) )
!
CALL READ_XKER_GWETH (KWETLBDAH,KWETLBDAG,KND,                              &
                   PALPHAH,PNUH,PALPHAG,PNUG,PEHG,PBG,PCH,PDH,PCG,PDG,      &
                   PWETLBDAH_MAX,PWETLBDAG_MAX,PWETLBDAH_MIN,PWETLBDAG_MIN, &
                   PFDINFTY                                                 )
IF( (KWETLBDAH/=NWETLBDAH) .OR. (KWETLBDAG/=NWETLBDAG) .OR. (KND/=IND) .OR. &
    (PALPHAH/=XALPHAH) .OR. (PNUH/=XNUH)                               .OR. &
    (PALPHAG/=XALPHAG) .OR. (PNUG/=XNUG)                               .OR. &
    (PEHG/=ZEHG) .OR. (PBG/=XBG)                                       .OR. &
    (PCH/=XCH) .OR. (PDH/=XDH) .OR. (PCG/=XCG) .OR. (PDG/=XDG)         .OR. &
    (PWETLBDAH_MAX/=XWETLBDAH_MAX) .OR. (PWETLBDAG_MAX/=XWETLBDAG_MAX) .OR. &
    (PWETLBDAH_MIN/=XWETLBDAH_MIN) .OR. (PWETLBDAG_MIN/=XWETLBDAG_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAH, XNUH, XALPHAG, XNUG,                          &
                ZEHG, XBG, XCH, XDH, XCG, XDG,                              &
                XWETLBDAH_MAX, XWETLBDAG_MAX, XWETLBDAH_MIN, XWETLBDAG_MIN, &
                ZFDINFTY, XKER_GWETH                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF GWETH KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAH=",I3)') NWETLBDAH
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAG=",I3)') NWETLBDAG
  WRITE(UNIT=KLUOUT,FMT='("PALPHAH=",E13.6)') XALPHAH
  WRITE(UNIT=KLUOUT,FMT='("PNUH=",E13.6)') XNUH
  WRITE(UNIT=KLUOUT,FMT='("PALPHAG=",E13.6)') XALPHAG
  WRITE(UNIT=KLUOUT,FMT='("PNUG=",E13.6)') XNUG
  WRITE(UNIT=KLUOUT,FMT='("PEHG=",E13.6)') ZEHG
  WRITE(UNIT=KLUOUT,FMT='("PBG=",E13.6)') XBG
  WRITE(UNIT=KLUOUT,FMT='("PCH=",E13.6)') XCH
  WRITE(UNIT=KLUOUT,FMT='("PDH=",E13.6)') XDH
  WRITE(UNIT=KLUOUT,FMT='("PCG=",E13.6)') XCG
  WRITE(UNIT=KLUOUT,FMT='("PDG=",E13.6)') XDG
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MAX=",E13.6)') &
                                                    XWETLBDAH_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAG_MAX=",E13.6)') &
                                                    XWETLBDAG_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MIN=",E13.6)') &
                                                    XWETLBDAH_MIN
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAG_MIN=",E13.6)') &
                                                    XWETLBDAG_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_GWETH) ) THEN")')
  DO J1 = 1 , NWETLBDAH
    DO J2 = 1 , NWETLBDAG
    WRITE(UNIT=KLUOUT,FMT='("PKER_GWETH(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_GWETH(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_GWETH (KWETLBDAH,KWETLBDAG,KND,                              &
                     PALPHAH,PNUH,PALPHAG,PNUG,PEHG,PBG,PCH,PDH,PCG,PDG,      &
                     PWETLBDAH_MAX,PWETLBDAG_MAX,PWETLBDAH_MIN,PWETLBDAG_MIN, &
                     PFDINFTY,XKER_GWETH                                      )
  WRITE(UNIT=KLUOUT,FMT='(" Read XKER_GWETH")')
END IF
!
!
!-------------------------------------------------------------------------------
!
!*      10.     SOME PRINTS FOR CONTROL
!               -----------------------
!
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='(" Summary of the ice particule characteristics")')
  WRITE(UNIT=KLUOUT,FMT='("      PRISTINE ICE")')
  WRITE(UNIT=KLUOUT,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAI,XBI
  WRITE(UNIT=KLUOUT,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XC_I,XDI
  WRITE(UNIT=KLUOUT,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAI,XNUI
  WRITE(UNIT=KLUOUT,FMT='("              SNOW")')
  WRITE(UNIT=KLUOUT,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAS,XBS
  WRITE(UNIT=KLUOUT,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XCS,XDS
  WRITE(UNIT=KLUOUT,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
                                                      XCCS,XCXS
  WRITE(UNIT=KLUOUT,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAS,XNUS
  WRITE(UNIT=KLUOUT,FMT='("            GRAUPEL")')
  WRITE(UNIT=KLUOUT,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAG,XBG
  WRITE(UNIT=KLUOUT,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XCG,XDG
  WRITE(UNIT=KLUOUT,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
                                                      XCCG,XCXG
  WRITE(UNIT=KLUOUT,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAG,XNUG
  WRITE(UNIT=KLUOUT,FMT='("               HAIL")')
  WRITE(UNIT=KLUOUT,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAH,XBH
  WRITE(UNIT=KLUOUT,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XCH,XDH
  WRITE(UNIT=KLUOUT,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
                                                      XCCH,XCXH
  WRITE(UNIT=KLUOUT,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAH,XNUH
END IF
!
KINTVL = IND
PFDINFTY = ZFDINFTY
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_RAIN_ICE_ELEC
