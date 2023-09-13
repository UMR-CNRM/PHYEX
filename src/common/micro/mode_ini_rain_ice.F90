!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
MODULE MODE_INI_RAIN_ICE
IMPLICIT NONE
CONTAINS
      SUBROUTINE INI_RAIN_ICE ( KLUOUT, PTSTEP, PDZMIN, KSPLITR, HCLOUD )
!     ###########################################################
!
!!****  *INI_RAIN_ICE * - initialize the constants necessary for the warm and
!!                        cold microphysical schemes.
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
!!    recomputed and their numerical values are written in the output listing.
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
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/95
!!      J.-P. Pinty 05/04/96 Add automatic control and regeneration of the
!!                           collection kernels
!!      J.-P. Pinty 10/05/96 Correction of ZRATE and computations of RIM
!!      J.-P. Pinty 24/11/97 Sedimentation of ice made for Columns and bug for XAG
!!      J.-P. Lafore 23/11/98 Back to Lin et al. 83 formulation for RIAUTS
!!                            with a Critical ice content set to .5 g/Kg
!!      N. Asencio  13/08/98 parallel code: PDZMIN is computed outside in ini_modeln
!!      J.-P. Lafore 12/8/98 In case of nesting microphysics constants of
!!                           MODD_RAIN_ICE_PARAM are computed only once.
!!                           Only KSPLTR is computed for each model.
!!      J. Stein    20/04/99 remove 2 unused local variables
!!      G Molinie   21/05/99 Bug in XEXRCFRI and XRCFRI
!!      J.-P. Pinty 24/06/00 Bug in RCRIMS
!!      J.-P. Pinty 24/12/00 Update hail case
!!      J.-P. Chaboureau & J.-P. Pinty
!!                  24/03/01 Update XCRIAUTI for cirrus cases
!!      J.-P. Pinty 24/11/01 Update ICE3/ICE4 options
!!      S. Riette 2016-11: new ICE3/ICE4 options
!!      P. Wautelet 22/01/2019 bug correction: incorrect write
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!!      S. Riette 2022-03: use of RAIN_ICE_PARAM structure for some variables
!!                         to reproduce results on belenos. The reason why
!!                         those variables must have a specifi treatment was
!!                         not understood
!  J. Wurtz       03/2022: New snow characteristics with LSNOW_T
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_PARAM_ICE_n, ONLY: LSNOW_T, CSEDIM, LRED, CPRISTINE_ICE, &
                        & LCRIAUTI, XACRIAUTI_NAM, XBCRIAUTI_NAM, XCRIAUTC_NAM, XCRIAUTI_NAM, XT0CRIAUTI_NAM, &
                        & XFRACM90, XFRMIN_NAM, XRDEPGRED_NAM, XRDEPSRED_NAM
USE MODD_RAIN_ICE_DESCR_n
USE MODD_RAIN_ICE_PARAM_n
!
USE MODI_GAMMA
USE MODI_GAMMA_INC
USE MODE_RRCOLSS, ONLY: RRCOLSS
USE MODE_RZCOLX, ONLY: RZCOLX
USE MODE_RSCOLRG, ONLY: RSCOLRG
USE MODE_READ_XKER_RACCS, ONLY: READ_XKER_RACCS
USE MODE_READ_XKER_SDRYG, ONLY: READ_XKER_SDRYG
USE MODE_READ_XKER_RDRYG, ONLY: READ_XKER_RDRYG
USE MODE_READ_XKER_SWETH, ONLY: READ_XKER_SWETH
USE MODE_READ_XKER_GWETH, ONLY: READ_XKER_GWETH
USE MODE_READ_XKER_RWETH, ONLY: READ_XKER_RWETH
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                 INTENT(IN) :: KLUOUT   ! Logical unit number for prints
INTEGER,                 INTENT(OUT):: KSPLITR   ! Number of small time step
                                                 ! integration for  rain
                                                 ! sedimendation
!
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step
!
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
CHARACTER (LEN=4), INTENT(IN)       :: HCLOUD    ! Indicator of the cloud scheme
!
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: J1,J2              ! Internal loop indexes
REAL :: ZT                    ! Work variable
REAL :: ZVTRMAX               ! Raindrop maximal fall velocity
REAL :: ZRHO00                ! Surface reference air density
REAL :: ZE, ZRV               ! Work array for ZRHO00 computation
REAL :: ZRATE                 ! Geometrical growth of Lbda in the tabulated
                              ! functions and kernels
REAL :: ZBOUND                ! XDCSLIM*Lbda_s: upper bound for the partial
                              ! integration of the riming rate of the aggregates
REAL :: ZEGS, ZEGR, ZEHS, &   ! Bulk collection efficiencies
      & ZEHG, ZEHR
!
INTEGER :: IND                ! Number of interval to integrate the kernels
REAL :: ZESR                  ! Mean efficiency of rain-aggregate collection
REAL :: ZFDINFTY              ! Factor used to define the "infinite" diameter
!
REAL :: ZCRI0, ZTCRI0         ! Second point to determine 10**(aT+b) law of ri->rs autoconversion
!
!
!
LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
                    ! listing
REAL     :: ZCONC_MAX ! Maximal concentration for snow
REAL     :: ZGAMC,ZGAMC2 ! parameters
                    ! involving various moments of the generalized gamma law
REAL     :: ZFACT_NUCL! Amplification factor for the minimal ice concentration
REAL     :: ZXR     ! Value of x_r in N_r = C_r lambda_r ** x_r
!
INTEGER  :: KND
INTEGER  :: KACCLBDAS,KACCLBDAR,KDRYLBDAG,KDRYLBDAS,KDRYLBDAR
INTEGER  :: KWETLBDAS,KWETLBDAG,KWETLBDAR,KWETLBDAH
REAL     :: PALPHAR,PALPHAS,PALPHAG,PALPHAH
REAL     :: PNUR,PNUS,PNUG,PNUH
REAL     :: PBR,PBS,PBG
REAL     :: PCR,PCS,PCG,PCH
REAL     :: PDR,PDS,PFVELOS,PDG,PDH
REAL     :: PESR,PEGS,PEGR,PEHS,PEHG,PEHR
REAL     :: PFDINFTY
REAL     :: PACCLBDAS_MAX,PACCLBDAR_MAX,PACCLBDAS_MIN,PACCLBDAR_MIN
REAL     :: PDRYLBDAG_MAX,PDRYLBDAS_MAX,PDRYLBDAG_MIN,PDRYLBDAS_MIN
REAL     :: PDRYLBDAR_MAX,PDRYLBDAR_MIN
REAL     :: PWETLBDAS_MAX,PWETLBDAG_MAX,PWETLBDAS_MIN,PWETLBDAG_MIN
REAL     :: PWETLBDAR_MAX,PWETLBDAH_MAX,PWETLBDAR_MIN,PWETLBDAH_MIN
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('INI_RAIN_ICE',0,ZHOOK_HANDLE)
!
!
!
!*       0.     FUNCTION STATEMENTS
!               -------------------
!
!
!*       0.1    p_moment of the Generalized GAMMA function
!
!
!
!        1.     COMPUTE KSPLTR FOR EACH MODEL
!               ---------------------------------------------------------
!
!*       1.1    Set the hailstones maximum fall velocity
!
IF (CSEDIM == 'SPLI' .AND. .NOT. LRED) THEN
 IF (HCLOUD == 'ICE4' .OR. HCLOUD=='LIMA') THEN
  ZVTRMAX = 40.
 ELSE IF (HCLOUD == 'ICE3') THEN
  ZVTRMAX = 10.
 END IF
END IF
!
!*       1.2    Compute the number of small time step integration
!
KSPLITR = 1
IF (CSEDIM == 'SPLI' .AND. .NOT. LRED) THEN
 SPLIT : DO
  ZT = PTSTEP / REAL(KSPLITR)
  IF ( ZT * ZVTRMAX / PDZMIN .LT. 1.) EXIT SPLIT
  KSPLITR = KSPLITR + 1
 END DO SPLIT
END IF
!
IF (HCLOUD == 'ICE4' .OR. HCLOUD=='LIMA') THEN
  CALL RAIN_ICE_DESCR_ALLOCATE(7)
ELSE IF (HCLOUD == 'ICE3') THEN
  CALL RAIN_ICE_DESCR_ALLOCATE(6)
END IF
!
XRTMIN(1) = 1.0E-20
XRTMIN(2) = 1.0E-20
XRTMIN(3) = 1.0E-20
XRTMIN(4) = 1.0E-20
XRTMIN(5) = 1.0E-15
XRTMIN(6) = 1.0E-15
IF (HCLOUD == 'ICE4' .OR. HCLOUD=='LIMA') XRTMIN(7) = 1.0E-15
!
!-------------------------------------------------------------------------------
!
!*       2.     CHARACTERISTICS OF THE SPECIES
!               ------------------------------
!
!
!*       2.1    Cloud droplet and Raindrop characteristics
!
XAC = (XPI/6.0)*XRHOLW
XBC = 3.0
XCC = XRHOLW*XG/(18.0*1.7E-5) ! Stokes flow (Pruppacher p 322 for T=273K)
XDC = 2.0
!
!
XAR = (XPI/6.0)*XRHOLW
XBR = 3.0
XCR = 842.
XDR = 0.8
!
!XCCR = 1.E7    ! N0_r =  XCXR * lambda_r ** ZXR
XCCR = 8.E6    ! N0_r =  XCXR * lambda_r ** ZXR
ZXR  = -1.     !
!
XF0R = 1.00
XF1R = 0.26
!
XC1R = 1./2.
!
!
!*       2.2    Ice crystal characteristics
!
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
!
XAS = 0.02
XBS = 1.9
IF (LSNOW_T) THEN
  !Cas Gamma generalisee
  XCS = 11.52
  XDS = 0.39
  XFVELOS =0.097
  !Cas MP
  !XCS = 13.2
  !XDS = 0.423       
  !XFVELOS = 25.14
ELSE
  XCS = 5.1
  XDS = 0.27
  XFVELOS = 0.
END IF
!
IF (.NOT. LSNOW_T) THEN
  XCCS = 5.0
  XCXS = 1.0
END IF
!
XF0S = 0.86
XF1S = 0.28
!
XC1S = 1./XPI
!
!
!*       2.4    Graupel/Frozen drop characteristics
!
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
XC1G = 1./2.
!
!
!*       2.5    Hailstone characteristics
!
!
XAH = 470.
XBH = 3.0
XCH = 207.
XDH = 0.64
!
!XCCH = 5.E-4
!XCXH = 2.0
!!!!!!!!!!!!
   XCCH = 4.E4 ! Test of Ziegler (1988)
   XCXH = -1.0 ! Test of Ziegler (1988)
!!!    XCCH = 5.E5 ! Graupel_like
!!!    XCXH = -0.5 ! Graupel_like
!!!!!!!!!!!!
!
XF0H = 0.86
XF1H = 0.28
!
XC1H = 1./2.
!
!-------------------------------------------------------------------------------
!
!*       3.     DIMENSIONAL DISTRIBUTIONS OF THE SPECIES
!               ----------------------------------------
!
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
IF (LSNOW_T) THEN
!Cas GAMMAGEN
   XALPHAS = .214   ! Generalized gamma law
   XNUS    = 43.7   ! Generalized gamma law
   XTRANS_MP_GAMMAS = SQRT( ( GAMMA(XNUS + 2./XALPHAS)*GAMMA(XNUS + 4./XALPHAS) ) / &
                            ( 8.* GAMMA(XNUS + 1./XALPHAS)*GAMMA(XNUS + 3./XALPHAS) ) )
ELSE
   XALPHAS = 1.0  ! Exponential law
   XNUS    = 1.0  ! Exponential law
   XTRANS_MP_GAMMAS = 1.
END IF
!
XALPHAG = 1.0  ! Exponential law
XNUG    = 1.0  ! Exponential law
!
XALPHAH = 1.0  ! Gamma law
XNUH    = 8.0  ! Gamma law with little dispersion
!
!*       3.4    Constants for shape parameter
!
ZGAMC = MOMG(XALPHAC,XNUC,3.)
ZGAMC2 = MOMG(XALPHAC2,XNUC2,3.)
XLBC(1)   = XAR*ZGAMC
XLBC(2)   = XAR*ZGAMC2
XLBEXC = 1.0/XBC
!
XLBEXR = 1.0/(-1.0-XBR)
XLBR   = ( XAR*XCCR*MOMG(XALPHAR,XNUR,XBR) )**(-XLBEXR)
!
XLBEXI = 1.0/(-XBI)
XLBI   = ( XAI*MOMG(XALPHAI,XNUI,XBI) )**(-XLBEXI)
!
#ifdef REPRO48
#else
XNS   = 1.0/(XAS*MOMG(XALPHAS,XNUS,XBS))
#endif
XLBEXS = 1.0/(XCXS-XBS)
XLBS   = ( XAS*XCCS*MOMG(XALPHAS,XNUS,XBS) )**(-XLBEXS)
!
XLBEXG = 1.0/(XCXG-XBG)
XLBG   = ( XAG*XCCG*MOMG(XALPHAG,XNUG,XBG) )**(-XLBEXG)
!
XLBEXH = 1.0/(XCXH-XBH)
XLBH   = ( XAH*XCCH*MOMG(XALPHAH,XNUH,XBH) )**(-XLBEXH)
!
!*       3.5    Minimal values allowed for the mixing ratios
!
XLBDAR_MAX = 100000.0
XLBDAS_MAX = 100000.0
XLBDAG_MAX = 100000.0
!
ZCONC_MAX  = 1.E6 ! Maximal concentration for falling particules set to 1 per cc
IF(XCCS>0. .AND. XCXS>0. )XLBDAS_MAX = ( ZCONC_MAX/XCCS )**(1./XCXS)
#ifdef REPRO48
#else
IF (LSNOW_T) XLBDAS_MAX = 1.E6
XLBDAS_MIN = 1.E-10
#endif
!
XCONC_SEA=1E8 ! 100/cm3
XCONC_LAND=3E8 ! 300/cm3
XCONC_URBAN=5E8 ! 500/cm3
!
!-------------------------------------------------------------------------------
!
!*       4.     CONSTANTS FOR THE SEDIMENTATION
!               -------------------------------
!
!
!*       4.1    Exponent of the fall-speed air density correction
!
XCEXVT = 0.4
!
!ZRHO00 = XP00/(XRD*XTHVREFZ(1+JPVEXT))
!According to Foote and Du Toit (1969) and List (1958), ZRHO00 must be computed for Hu=50%, P=101325Pa and T=293.15K
ZE = (50./100.) * EXP(XALPW-XBETAW/293.15-XGAMW*LOG(293.15))
ZRV = (XRD/XRV) * ZE / (101325.-ZE)
ZRHO00 = 101325.*(1.+ZRV)/(XRD+ZRV*XRV)/293.15
!
!*       4.2    Constants for sedimentation
!
XFSEDC(1)  = GAMMA(XNUC+(XDC+3.)/XALPHAC)/GAMMA(XNUC+3./XALPHAC)*     &
            (ZRHO00)**XCEXVT
XFSEDC(2)  = GAMMA(XNUC2+(XDC+3.)/XALPHAC2)/GAMMA(XNUC2+3./XALPHAC2)*     &
            (ZRHO00)**XCEXVT
!
XEXSEDR = (XBR+XDR+1.0)/(XBR+1.0)
XFSEDR  = XCR*XAR*XCCR*MOMG(XALPHAR,XNUR,XBR+XDR)*                         &
            (XAR*XCCR*MOMG(XALPHAR,XNUR,XBR))**(-XEXSEDR)*(ZRHO00)**XCEXVT
!
XEXRSEDI = (XBI+XDI)/XBI
XEXCSEDI = 1.0-XEXRSEDI
XFSEDI   = (4.*XPI*900.)**(-XEXCSEDI) *                         &
           XC_I*XAI*MOMG(XALPHAI,XNUI,XBI+XDI) *                &
           ((XAI*MOMG(XALPHAI,XNUI,XBI)))**(-XEXRSEDI) *        &
           (ZRHO00)**XCEXVT
!When we do not use computations for columns, I think we must uncomment line just below
!XEXCSEDI = XEXCSEDI * 3. to be checked
!
!  Computations made for Columns
!
XEXRSEDI = 1.9324
XEXCSEDI =-0.9324
XFSEDI   = 3.89745E11*MOMG(XALPHAI,XNUI,3.285)*                          &
                      MOMG(XALPHAI,XNUI,1.7)**(-XEXRSEDI)*(ZRHO00)**XCEXVT
XEXCSEDI =-0.9324*3.0
WRITE (KLUOUT,FMT=*)' PRISTINE ICE SEDIMENTATION for columns XFSEDI =',XFSEDI
!
!
#ifdef REPRO48
XEXSEDS = (XBS+XDS-XCXS)/(XBS-XCXS)
XFSEDS  = XCS*XAS*XCCS*MOMG(XALPHAS,XNUS,XBS+XDS)*                         &
         (XAS*XCCS*MOMG(XALPHAS,XNUS,XBS))**(-XEXSEDS)*(ZRHO00)**XCEXVT
#else
IF (LRED) THEN
   XEXSEDS = -XDS-XBS
   XFSEDS  = XCS*MOMG(XALPHAS,XNUS,XBS+XDS)/(MOMG(XALPHAS,XNUS,XBS))    &
            *(ZRHO00)**XCEXVT
ELSE
   XEXSEDS = (XBS+XDS-XCXS)/(XBS-XCXS)
   XFSEDS  = XCS*XAS*XCCS*MOMG(XALPHAS,XNUS,XBS+XDS)*                         &
            (XAS*XCCS*MOMG(XALPHAS,XNUS,XBS))**(-XEXSEDS)*(ZRHO00)**XCEXVT
END IF
#endif
!
XEXSEDG = (XBG+XDG-XCXG)/(XBG-XCXG)
XFSEDG  = XCG*XAG*XCCG*MOMG(XALPHAG,XNUG,XBG+XDG)*                         &
            (XAG*XCCG*MOMG(XALPHAG,XNUG,XBG))**(-XEXSEDG)*(ZRHO00)**XCEXVT
!
XEXSEDH = (XBH+XDH-XCXH)/(XBH-XCXH)
XFSEDH  = XCH*XAH*XCCH*MOMG(XALPHAH,XNUH,XBH+XDH)*                         &
            (XAH*XCCH*MOMG(XALPHAH,XNUH,XBH))**(-XEXSEDH)*(ZRHO00)**XCEXVT
!
!
!-------------------------------------------------------------------------------
!
!*       5.     CONSTANTS FOR THE SLOW COLD PROCESSES
!               -------------------------------------
!
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
XNU10 = 50.*ZFACT_NUCL
XALPHA1 = 4.5
XBETA1 = 0.6
!
XNU20 = 1000.*ZFACT_NUCL
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
XHON = (XPI/6.)*((2.0*3.0*4.0*5.0*6.0)/(2.0*3.0))*(1.1E5)**(-3.0) !
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
XSCFAC = (0.63**(1./3.))*SQRT((ZRHO00)**XCEXVT) ! One assumes Sc=0.63
!
X0DEPI = (4.0*XPI)*XC1I*XF0I*MOMG(XALPHAI,XNUI,1.)
X2DEPI = (4.0*XPI)*XC1I*XF2I*XC_I*MOMG(XALPHAI,XNUI,XDI+2.0)
!
#ifdef REPRO48
X0DEPS = (4.0*XPI)*XCCS*XC1S*XF0S*MOMG(XALPHAS,XNUS,1.)
X1DEPS = (4.0*XPI)*XCCS*XC1S*XF1S*SQRT(XCS)*MOMG(XALPHAS,XNUS,0.5*XDS+1.5)
XEX0DEPS = XCXS-1.0
XEX1DEPS = XCXS-0.5*(XDS+3.0)
#else
X0DEPS = XNS*(4.0*XPI)*XC1S*XF0S*MOMG(XALPHAS,XNUS,1.)
X1DEPS = XNS*(4.0*XPI)*XC1S*XF1S*SQRT(XCS)*MOMG(XALPHAS,XNUS,0.5*XDS+1.5)
XEX0DEPS = -1.0
XEX1DEPS = -0.5*(XDS+3.0)
#endif
XRDEPSRED = XRDEPSRED_NAM
!
X0DEPG = (4.0*XPI)*XCCG*XC1G*XF0G*MOMG(XALPHAG,XNUG,1.)
X1DEPG = (4.0*XPI)*XCCG*XC1G*XF1G*SQRT(XCG)*MOMG(XALPHAG,XNUG,0.5*XDG+1.5)
XEX0DEPG = XCXG-1.0
XEX1DEPG = XCXG-0.5*(XDG+3.0)
XRDEPGRED = XRDEPGRED_NAM
!
X0DEPH = (4.0*XPI)*XCCH*XC1H*XF0H*MOMG(XALPHAH,XNUH,1.)
X1DEPH = (4.0*XPI)*XCCH*XC1H*XF1H*SQRT(XCH)*MOMG(XALPHAH,XNUH,0.5*XDH+1.5)
XEX0DEPH = XCXH-1.0
XEX1DEPH = XCXH-0.5*(XDH+3.0)

GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      factors sublimation snow/groupel")')
  WRITE(UNIT=KLUOUT,FMT='(" mod sublim snow =",E13.6)') XRDEPSRED
  WRITE(UNIT=KLUOUT,FMT='(" mod sublim graupel =",E13.6)') XRDEPGRED
END IF

!
!*       5.3    Constants for pristine ice autoconversion
!
XTIMAUTI = 1.E-3  !  Time constant at T=T_t
XTEXAUTI = 0.015  !  Temperature factor of the I+I collection efficiency
XCRIAUTI = XCRIAUTI_NAM
IF(LCRIAUTI) THEN
  XT0CRIAUTI = XT0CRIAUTI_NAM
  !second point to determine 10**(aT+b) law
  ZTCRI0=-40.0
  ZCRI0=1.25E-6
  XBCRIAUTI=-( LOG10(XCRIAUTI) - LOG10(ZCRI0)*XT0CRIAUTI/ZTCRI0 )&
                   *ZTCRI0/(XT0CRIAUTI-ZTCRI0)
  XACRIAUTI=(LOG10(ZCRI0)-XBCRIAUTI)/ZTCRI0
ELSE
  XACRIAUTI=XACRIAUTI_NAM
  XBCRIAUTI=XBCRIAUTI_NAM
  XT0CRIAUTI=(LOG10(XCRIAUTI)-XBCRIAUTI)/0.06
ENDIF
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      pristine ice autoconversion")')
  WRITE(UNIT=KLUOUT,FMT='(" Time constant   XTIMAUTI=",E13.6)') XTIMAUTI
  WRITE(UNIT=KLUOUT,FMT='(" Temp. factor    XTEXAUTI=",E13.6)') XTEXAUTI
  WRITE(UNIT=KLUOUT,FMT='(" Crit. ice cont. XCRIAUTI=",E13.6)') XCRIAUTI
  WRITE(UNIT=KLUOUT,FMT='(" A Coef. for cirrus law XACRIAUTI=",E13.6)')XACRIAUTI
  WRITE(UNIT=KLUOUT,FMT='(" B Coef. for cirrus law XBCRIAUTI=",E13.6)')XBCRIAUTI
  WRITE(UNIT=KLUOUT,FMT='(" Temp degC at which cirrus law starts to be used=",E13.6)') XT0CRIAUTI
END IF
!
!
!*       5.4    Constants for snow aggregation
!
XCOLIS   = 0.25 ! Collection efficiency of I+S
XCOLEXIS = 0.05 ! Temperature factor of the I+S collection efficiency
#ifdef REPRO48
XFIAGGS  = (XPI/4.0)*XCOLIS*XCCS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
XEXIAGGS = XCXS-XDS-2.0
#else
XFIAGGS  = XNS*(XPI/4.0)*XCOLIS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
XEXIAGGS = -XDS - 2.0 ! GAMMGEN LH_EXTENDED
#endif
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
!               -------------------------------------
!
!
!*       6.1    Constants for the cloud droplets autoconversion
!
XTIMAUTC = 1.E-3
XCRIAUTC = XCRIAUTC_NAM
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
XFCACCR  = (XPI/4.0)*XCCR*XCR*(ZRHO00**XCEXVT)*MOMG(XALPHAR,XNUR,XDR+2.0)
XEXCACCR = -XDR-3.0
!
!*       6.3    Constants for the evaporation of the raindrops
!
X0EVAR = (4.0*XPI)*XCCR*XC1R*XF0R*MOMG(XALPHAR,XNUR,1.)
X1EVAR = (4.0*XPI)*XCCR*XC1R*XF1R*SQRT(XCR)*MOMG(XALPHAR,XNUR,0.5*XDR+1.5)
XEX0EVAR = -2.0
XEX1EVAR = -1.0-0.5*(XDR+3.0)
!
!
!-------------------------------------------------------------------------------
!
!*       7.     CONSTANTS FOR THE FAST COLD PROCESSES FOR THE AGGREGATES
!               --------------------------------------------------------
!
!
!*       7.1    Constants for the riming of the aggregates
!
XDCSLIM  = 0.007 ! D_cs^lim = 7 mm as suggested by Farley et al. (1989)
XCOLCS   = 1.0
#ifdef REPRO48
XEXCRIMSS= XCXS-XDS-2.0
XCRIMSS  = (XPI/4.0)*XCOLCS*XCCS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
#else
XEXCRIMSS= -XDS-2.0
XCRIMSS  = XNS * (XPI/4.0)*XCOLCS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
#endif
XEXCRIMSG= XEXCRIMSS
XCRIMSG  = XCRIMSS
#ifdef REPRO48
XSRIMCG  = XCCS*XAS*MOMG(XALPHAS,XNUS,XBS)
XEXSRIMCG= XCXS-XBS
XSRIMCG2 = XCCS*XAG*MOMG(XALPHAS,XNUS,XBG)
XSRIMCG3 = XFRACM90
XEXSRIMCG2=XCXS-XBG
#else
XSRIMCG  = XNS*XAS*MOMG(XALPHAS,XNUS,XBS)
XEXSRIMCG = -XBS
XSRIMCG2 = XNS*XAG*MOMG(XALPHAS,XNUS,XBG)
XSRIMCG3 = XFRACM90
XEXSRIMCG2=XBS-XBG
#endif
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=KLUOUT,FMT='("      riming of the aggregates")')
  WRITE(UNIT=KLUOUT,FMT='(" D_cs^lim (Farley et al.) XDCSLIM=",E13.6)') XDCSLIM
  WRITE(UNIT=KLUOUT,FMT='(" Coll. efficiency          XCOLCS=",E13.6)') XCOLCS
END IF
!
RAIN_ICE_PARAMN%NGAMINC = 80
RAIN_ICE_PARAMN%XGAMINC_BOUND_MIN = 1.0E-1 ! Minimal value of (Lbda * D_cs^lim)**alpha
RAIN_ICE_PARAMN%XGAMINC_BOUND_MAX = 1.0E7  ! Maximal value of (Lbda * D_cs^lim)**alpha
ZRATE = EXP(LOG(RAIN_ICE_PARAMN%XGAMINC_BOUND_MAX/RAIN_ICE_PARAMN%XGAMINC_BOUND_MIN)/REAL(RAIN_ICE_PARAMN%NGAMINC-1))
!
IF( .NOT.ASSOCIATED(XGAMINC_RIM1) ) CALL RAIN_ICE_PARAM_ALLOCATE('XGAMINC_RIM1', RAIN_ICE_PARAMN%NGAMINC)
IF( .NOT.ASSOCIATED(XGAMINC_RIM2) ) CALL RAIN_ICE_PARAM_ALLOCATE('XGAMINC_RIM2', RAIN_ICE_PARAMN%NGAMINC)
IF( .NOT.ASSOCIATED(XGAMINC_RIM4) ) CALL RAIN_ICE_PARAM_ALLOCATE('XGAMINC_RIM4', RAIN_ICE_PARAMN%NGAMINC)
!
DO J1=1,RAIN_ICE_PARAMN%NGAMINC
  ZBOUND = RAIN_ICE_PARAMN%XGAMINC_BOUND_MIN*ZRATE**(J1-1)
  XGAMINC_RIM1(J1) = GAMMA_INC(XNUS+(2.0+XDS)/XALPHAS,ZBOUND)
  XGAMINC_RIM2(J1) = GAMMA_INC(XNUS+XBS/XALPHAS      ,ZBOUND)
  XGAMINC_RIM4(J1) = GAMMA_INC(XNUS+XBG/XALPHAS      ,ZBOUND)
END DO
!
RAIN_ICE_PARAMN%XRIMINTP1 = XALPHAS / LOG(ZRATE)
RAIN_ICE_PARAMN%XRIMINTP2 = 1.0 + RAIN_ICE_PARAMN%XRIMINTP1*LOG( XDCSLIM/(RAIN_ICE_PARAMN%XGAMINC_BOUND_MIN)**(1.0/XALPHAS) )
!
!*       7.2    Constants for the accretion of raindrops onto aggregates
!
#ifdef REPRO48
XFRACCSS = ((XPI**2)/24.0)*XCCS*XCCR*XRHOLW*(ZRHO00**XCEXVT)
#else
XFRACCSS = XNS*((XPI**2)/24.0)*XCCR*XRHOLW*(ZRHO00**XCEXVT)
#endif
!
XLBRACCS1   =    MOMG(XALPHAS,XNUS,2.)*MOMG(XALPHAR,XNUR,3.)
XLBRACCS2   = 2.*MOMG(XALPHAS,XNUS,1.)*MOMG(XALPHAR,XNUR,4.)
XLBRACCS3   =                          MOMG(XALPHAR,XNUR,5.)
!
#ifdef REPRO48
XFSACCRG = (XPI/4.0)*XAS*XCCS*XCCR*(ZRHO00**XCEXVT)
#else
XFSACCRG = XNS*(XPI/4.0)*XAS*XCCR*(ZRHO00**XCEXVT)
#endif
!
XLBSACCR1   =    MOMG(XALPHAR,XNUR,2.)*MOMG(XALPHAS,XNUS,XBS)
XLBSACCR2   = 2.*MOMG(XALPHAR,XNUR,1.)*MOMG(XALPHAS,XNUS,XBS+1.)
XLBSACCR3   =                          MOMG(XALPHAS,XNUS,XBS+2.)
!
!*       7.2.1  Defining the ranges for the computation of the kernels
!
! Notice: One magnitude of lambda discretized over 10 points for rain
! Notice: One magnitude of lambda discretized over 10 points for snow
!
RAIN_ICE_PARAMN%NACCLBDAS = 40
RAIN_ICE_PARAMN%XACCLBDAS_MIN = 5.0E1 ! Minimal value of Lbda_s to tabulate XKER_RACCS
RAIN_ICE_PARAMN%XACCLBDAS_MAX = 5.0E5 ! Maximal value of Lbda_s to tabulate XKER_RACCS
ZRATE = LOG(RAIN_ICE_PARAMN%XACCLBDAS_MAX/RAIN_ICE_PARAMN%XACCLBDAS_MIN)/REAL(RAIN_ICE_PARAMN%NACCLBDAS-1)
RAIN_ICE_PARAMN%XACCINTP1S = 1.0 / ZRATE
RAIN_ICE_PARAMN%XACCINTP2S = 1.0 - LOG( RAIN_ICE_PARAMN%XACCLBDAS_MIN ) / ZRATE
RAIN_ICE_PARAMN%NACCLBDAR = 40
RAIN_ICE_PARAMN%XACCLBDAR_MIN = 1.0E3 ! Minimal value of Lbda_r to tabulate XKER_RACCS
RAIN_ICE_PARAMN%XACCLBDAR_MAX = 1.0E7 ! Maximal value of Lbda_r to tabulate XKER_RACCS
ZRATE = LOG(RAIN_ICE_PARAMN%XACCLBDAR_MAX/RAIN_ICE_PARAMN%XACCLBDAR_MIN)/REAL(RAIN_ICE_PARAMN%NACCLBDAR-1)
RAIN_ICE_PARAMN%XACCINTP1R = 1.0 / ZRATE
RAIN_ICE_PARAMN%XACCINTP2R = 1.0 - LOG( RAIN_ICE_PARAMN%XACCLBDAR_MIN ) / ZRATE
!
!*       7.2.2  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZESR     = 1.0   ! factor used to integrate the dimensional distributions when
ZFDINFTY = 20.0  ! computing the kernels XKER_RACCSS, XKER_RACCS and XKER_SACCRG
!
IF( .NOT.ASSOCIATED(XKER_RACCSS) ) CALL RAIN_ICE_PARAM_ALLOCATE('XKER_RACCSS', RAIN_ICE_PARAMN%NACCLBDAS,RAIN_ICE_PARAMN%NACCLBDAR)
IF( .NOT.ASSOCIATED(XKER_RACCS ) ) CALL RAIN_ICE_PARAM_ALLOCATE('XKER_RACCS', RAIN_ICE_PARAMN%NACCLBDAS,RAIN_ICE_PARAMN%NACCLBDAR)
IF( .NOT.ASSOCIATED(XKER_SACCRG) ) CALL RAIN_ICE_PARAM_ALLOCATE('XKER_SACCRG', RAIN_ICE_PARAMN%NACCLBDAR,RAIN_ICE_PARAMN%NACCLBDAS)
!
CALL READ_XKER_RACCS (KACCLBDAS,KACCLBDAR,KND,                                        &
                      PALPHAS,PNUS,PALPHAR,PNUR,PESR,PBS,PBR,PCS,PDS,PFVELOS,PCR,PDR, &
                      PACCLBDAS_MAX,PACCLBDAR_MAX,PACCLBDAS_MIN,PACCLBDAR_MIN,        &
                      PFDINFTY                                                        )
IF( (KACCLBDAS/=RAIN_ICE_PARAMN%NACCLBDAS) .OR. (KACCLBDAR/=RAIN_ICE_PARAMN%NACCLBDAR) .OR. (KND/=IND) .OR. &
    (PALPHAS/=XALPHAS) .OR. (PNUS/=XNUS)                               .OR. &
    (PALPHAR/=XALPHAR) .OR. (PNUR/=XNUR)                               .OR. &
    (PESR/=ZESR) .OR. (PBS/=XBS) .OR. (PBR/=XBR)                       .OR. &
    (PCS/=XCS) .OR. (PDS/=XDS) .OR. (PFVELOS/=XFVELOS) .OR. (PCR/=XCR) .OR. (PDR/=XDR) .OR. &
    (PACCLBDAS_MAX/=RAIN_ICE_PARAMN%XACCLBDAS_MAX) .OR. (PACCLBDAR_MAX/=RAIN_ICE_PARAMN%XACCLBDAR_MAX) .OR. &
    (PACCLBDAS_MIN/=RAIN_ICE_PARAMN%XACCLBDAS_MIN) .OR. (PACCLBDAR_MIN/=RAIN_ICE_PARAMN%XACCLBDAR_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RRCOLSS ( IND, XALPHAS, XNUS, XALPHAR, XNUR,                          &
                 ZESR, XBR, XCS, XDS, XFVELOS, XCR, XDR,                     &
                 RAIN_ICE_PARAMN%XACCLBDAS_MAX, RAIN_ICE_PARAMN%XACCLBDAR_MAX, &
                 RAIN_ICE_PARAMN%XACCLBDAS_MIN, RAIN_ICE_PARAMN%XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_RACCSS, XAG, XBS, XAS                        )
  CALL RZCOLX  ( IND, XALPHAS, XNUS, XALPHAR, XNUR,                          &
                 ZESR, XBR, XCS, XDS, XFVELOS, XCR, XDR, 0.,                 &
                 RAIN_ICE_PARAMN%XACCLBDAS_MAX, RAIN_ICE_PARAMN%XACCLBDAR_MAX, &
                 RAIN_ICE_PARAMN%XACCLBDAS_MIN, RAIN_ICE_PARAMN%XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_RACCS                                        )
  CALL RSCOLRG ( IND, XALPHAS, XNUS, XALPHAR, XNUR,                          &
                 ZESR, XBS, XCS, XDS, XFVELOS, XCR, XDR,                     & 
                 RAIN_ICE_PARAMN%XACCLBDAS_MAX, RAIN_ICE_PARAMN%XACCLBDAR_MAX, &
                 RAIN_ICE_PARAMN%XACCLBDAS_MIN, RAIN_ICE_PARAMN%XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_SACCRG,  XAG, XBS, XAS                       )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF RACSS KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF RACS  KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF SACRG KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KACCLBDAS=",I3)') RAIN_ICE_PARAMN%NACCLBDAS
  WRITE(UNIT=KLUOUT,FMT='("KACCLBDAR=",I3)') RAIN_ICE_PARAMN%NACCLBDAR
  WRITE(UNIT=KLUOUT,FMT='("PALPHAS=",E13.6)') XALPHAS
  WRITE(UNIT=KLUOUT,FMT='("PNUS=",E13.6)') XNUS
  WRITE(UNIT=KLUOUT,FMT='("PALPHAR=",E13.6)') XALPHAR
  WRITE(UNIT=KLUOUT,FMT='("PNUR=",E13.6)') XNUR
  WRITE(UNIT=KLUOUT,FMT='("PESR=",E13.6)') ZESR
  WRITE(UNIT=KLUOUT,FMT='("PBS=",E13.6)') XBS
  WRITE(UNIT=KLUOUT,FMT='("PBR=",E13.6)') XBR
  WRITE(UNIT=KLUOUT,FMT='("PCS=",E13.6)') XCS
  WRITE(UNIT=KLUOUT,FMT='("PDS=",E13.6)') XDS
  WRITE(UNIT=KLUOUT,FMT='("PFVELOS=",E13.6)') XFVELOS
  WRITE(UNIT=KLUOUT,FMT='("PCR=",E13.6)') XCR
  WRITE(UNIT=KLUOUT,FMT='("PDR=",E13.6)') XDR
  WRITE(UNIT=KLUOUT,FMT='("PACCLBDAS_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XACCLBDAS_MAX
  WRITE(UNIT=KLUOUT,FMT='("PACCLBDAR_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XACCLBDAR_MAX
  WRITE(UNIT=KLUOUT,FMT='("PACCLBDAS_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XACCLBDAS_MIN
  WRITE(UNIT=KLUOUT,FMT='("PACCLBDAR_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XACCLBDAR_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_RACCSS) ) THEN")')
  DO J1 = 1 , RAIN_ICE_PARAMN%NACCLBDAS
    DO J2 = 1 , RAIN_ICE_PARAMN%NACCLBDAR
    WRITE(UNIT=KLUOUT,FMT='("  PKER_RACCSS(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RACCSS(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_RACCS ) ) THEN")')
  DO J1 = 1 , RAIN_ICE_PARAMN%NACCLBDAS
    DO J2 = 1 , RAIN_ICE_PARAMN%NACCLBDAR
    WRITE(UNIT=KLUOUT,FMT='("  PKER_RACCS (",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RACCS (J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_SACCRG) ) THEN")')
  DO J1 = 1 , RAIN_ICE_PARAMN%NACCLBDAR
    DO J2 = 1 , RAIN_ICE_PARAMN%NACCLBDAS
    WRITE(UNIT=KLUOUT,FMT='("  PKER_SACCRG(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_SACCRG(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_RACCS (KACCLBDAS,KACCLBDAR,KND,                                       &
                       PALPHAS,PNUS,PALPHAR,PNUR,PESR,PBS,PBR,PCS,PDS,PFVELOS,PCR,PDR, &
                       PACCLBDAS_MAX,PACCLBDAR_MAX,PACCLBDAS_MIN,PACCLBDAR_MIN,        &
                       PFDINFTY,XKER_RACCSS,XKER_RACCS,XKER_SACCRG                     )
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
XEXRCFRI  = -XDR-5.0+ZXR
XRCFRI    = ((XPI**2)/24.0)*XCCR*XRHOLW*XCOLIR*XCR*(ZRHO00**XCEXVT)     &
                                                     *MOMG(XALPHAR,XNUR,XDR+5.0)
XEXICFRR  = -XDR-2.0+ZXR
XICFRR    = (XPI/4.0)*XCOLIR*XCR*(ZRHO00**XCEXVT)          &
                                   *XCCR*MOMG(XALPHAR,XNUR,XDR+2.0)
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
XFCDRYG = (XPI/4.0)*XCCG*XCG*(ZRHO00**XCEXVT)*MOMG(XALPHAG,XNUG,XDG+2.0)
!
!*       8.2.2  Constants for the cloud ice collection by the graupeln
!
XCOLIG    = 0.25 ! Collection efficiency of I+G
XCOLEXIG  = 0.05 ! Temperature factor of the I+G collection efficiency
XCOLIG   = 0.01 ! Collection efficiency of I+G
XCOLEXIG = 0.1  ! Temperature factor of the I+G collection efficiency
WRITE (KLUOUT, FMT=*) ' NEW Constants for the cloud ice collection by the graupeln'
WRITE (KLUOUT, FMT=*) ' XCOLIG, XCOLEXIG  = ',XCOLIG,XCOLEXIG
XFIDRYG = (XPI/4.0)*XCOLIG*XCCG*XCG*(ZRHO00**XCEXVT)*MOMG(XALPHAG,XNUG,XDG+2.0)
XEXFIDRYG=(XCXG-XDG-2.)/(XCXG-XBG)
XFIDRYG2=XFIDRYG/XCOLIG*(XAG*XCCG*MOMG(XALPHAG,XNUG,XBG))**(-XEXFIDRYG)
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
XCOLSG    = 0.25 ! Collection efficiency of S+G
XCOLEXSG  = 0.05 ! Temperature factor of the S+G collection efficiency
XCOLSG   = 0.01 ! Collection efficiency of S+G
XCOLEXSG = 0.1  ! Temperature factor of the S+G collection efficiency
WRITE (KLUOUT, FMT=*) ' NEW Constants for the aggregate collection by the graupeln'
WRITE (KLUOUT, FMT=*) ' XCOLSG, XCOLEXSG  = ',XCOLSG,XCOLEXSG
#ifdef REPRO48
XFSDRYG = (XPI/4.0)*XCOLSG*XCCG*XCCS*XAS*(ZRHO00**XCEXVT)
#else
XFSDRYG = XNS*(XPI/4.0)*XCOLSG*XCCG*XAS*(ZRHO00**XCEXVT)
#endif
!
XLBSDRYG1   =    MOMG(XALPHAG,XNUG,2.)*MOMG(XALPHAS,XNUS,XBS)
XLBSDRYG2   = 2.*MOMG(XALPHAG,XNUG,1.)*MOMG(XALPHAS,XNUS,XBS+1.)
XLBSDRYG3   =                          MOMG(XALPHAS,XNUS,XBS+2.)
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
XFRDRYG = ((XPI**2)/24.0)*XCCG*XCCR*XRHOLW*(ZRHO00**XCEXVT)
!
XLBRDRYG1   =    MOMG(XALPHAG,XNUG,2.)*MOMG(XALPHAR,XNUR,3.)
XLBRDRYG2   = 2.*MOMG(XALPHAG,XNUG,1.)*MOMG(XALPHAR,XNUR,4.)
XLBRDRYG3   =                          MOMG(XALPHAR,XNUR,5.)
!
! Notice: One magnitude of lambda discretized over 10 points
!
RAIN_ICE_PARAMN%NDRYLBDAR = 40
RAIN_ICE_PARAMN%XDRYLBDAR_MIN = 1.0E3 ! Minimal value of Lbda_r to tabulate XKER_RDRYG
RAIN_ICE_PARAMN%XDRYLBDAR_MAX = 1.0E7 ! Maximal value of Lbda_r to tabulate XKER_RDRYG
ZRATE = LOG(RAIN_ICE_PARAMN%XDRYLBDAR_MAX/RAIN_ICE_PARAMN%XDRYLBDAR_MIN)/REAL(RAIN_ICE_PARAMN%NDRYLBDAR-1)
RAIN_ICE_PARAMN%XDRYINTP1R = 1.0 / ZRATE
RAIN_ICE_PARAMN%XDRYINTP2R = 1.0 - LOG( RAIN_ICE_PARAMN%XDRYLBDAR_MIN ) / ZRATE
RAIN_ICE_PARAMN%NDRYLBDAS = 80
RAIN_ICE_PARAMN%XDRYLBDAS_MIN = 2.5E1 ! Minimal value of Lbda_s to tabulate XKER_SDRYG
RAIN_ICE_PARAMN%XDRYLBDAS_MAX = 2.5E9 ! Maximal value of Lbda_s to tabulate XKER_SDRYG
ZRATE = LOG(RAIN_ICE_PARAMN%XDRYLBDAS_MAX/RAIN_ICE_PARAMN%XDRYLBDAS_MIN)/REAL(RAIN_ICE_PARAMN%NDRYLBDAS-1)
RAIN_ICE_PARAMN%XDRYINTP1S = 1.0 / ZRATE
RAIN_ICE_PARAMN%XDRYINTP2S = 1.0 - LOG( RAIN_ICE_PARAMN%XDRYLBDAS_MIN ) / ZRATE
RAIN_ICE_PARAMN%NDRYLBDAG = 40
RAIN_ICE_PARAMN%XDRYLBDAG_MIN = 1.0E3 ! Min value of Lbda_g to tabulate XKER_SDRYG,XKER_RDRYG
RAIN_ICE_PARAMN%XDRYLBDAG_MAX = 1.0E7 ! Max value of Lbda_g to tabulate XKER_SDRYG,XKER_RDRYG
ZRATE = LOG(RAIN_ICE_PARAMN%XDRYLBDAG_MAX/RAIN_ICE_PARAMN%XDRYLBDAG_MIN)/REAL(RAIN_ICE_PARAMN%NDRYLBDAG-1)
RAIN_ICE_PARAMN%XDRYINTP1G = 1.0 / ZRATE
RAIN_ICE_PARAMN%XDRYINTP2G = 1.0 - LOG( RAIN_ICE_PARAMN%XDRYLBDAG_MIN ) / ZRATE
!
!*       8.2.5  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZEGS     = 1.0   ! factor used to integrate the dimensional distributions when
ZFDINFTY = 20.0  ! computing the kernels XKER_SDRYG
!
IF( .NOT.ASSOCIATED(XKER_SDRYG) ) CALL RAIN_ICE_PARAM_ALLOCATE('XKER_SDRYG', RAIN_ICE_PARAMN%NDRYLBDAG,RAIN_ICE_PARAMN%NDRYLBDAS)
!
CALL READ_XKER_SDRYG (KDRYLBDAG,KDRYLBDAS,KND,                              &
                   PALPHAG,PNUG,PALPHAS,PNUS,PEGS,PBS,PCG,PDG,PCS,PDS,PFVELOS, &
                   PDRYLBDAG_MAX,PDRYLBDAS_MAX,PDRYLBDAG_MIN,PDRYLBDAS_MIN, &
                   PFDINFTY                                                 )
IF( (KDRYLBDAG/=RAIN_ICE_PARAMN%NDRYLBDAG) .OR. (KDRYLBDAS/=RAIN_ICE_PARAMN%NDRYLBDAS) .OR. (KND/=IND) .OR. &
    (PALPHAG/=XALPHAG) .OR. (PNUG/=XNUG)                               .OR. &
    (PALPHAS/=XALPHAS) .OR. (PNUS/=XNUS)                               .OR. &
    (PEGS/=ZEGS) .OR. (PBS/=XBS)                                       .OR. &
    (PCG/=XCG) .OR. (PDG/=XDG) .OR. (PCS/=XCS) .OR. (PDS/=XDS) .OR. (PFVELOS/=XFVELOS) .OR. &
    (PDRYLBDAG_MAX/=RAIN_ICE_PARAMN%XDRYLBDAG_MAX) .OR. (PDRYLBDAS_MAX/=RAIN_ICE_PARAMN%XDRYLBDAS_MAX) .OR. &
    (PDRYLBDAG_MIN/=RAIN_ICE_PARAMN%XDRYLBDAG_MIN) .OR. (PDRYLBDAS_MIN/=RAIN_ICE_PARAMN%XDRYLBDAS_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAG, XNUG, XALPHAS, XNUS,                          &
                ZEGS, XBS, XCG, XDG, 0., XCS, XDS, XFVELOS,                 &
                RAIN_ICE_PARAMN%XDRYLBDAG_MAX, RAIN_ICE_PARAMN%XDRYLBDAS_MAX, &
                RAIN_ICE_PARAMN%XDRYLBDAG_MIN, RAIN_ICE_PARAMN%XDRYLBDAS_MIN, &
                ZFDINFTY, XKER_SDRYG                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF SDRYG KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KDRYLBDAG=",I3)') RAIN_ICE_PARAMN%NDRYLBDAG
  WRITE(UNIT=KLUOUT,FMT='("KDRYLBDAS=",I3)') RAIN_ICE_PARAMN%NDRYLBDAS
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
  WRITE(UNIT=KLUOUT,FMT='("PFVELOS=",E13.6)') XFVELOS
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAG_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XDRYLBDAG_MAX
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAS_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XDRYLBDAS_MAX
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAG_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XDRYLBDAG_MIN
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAS_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XDRYLBDAS_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_SDRYG) ) THEN")')
  DO J1 = 1 , RAIN_ICE_PARAMN%NDRYLBDAG
    DO J2 = 1 , RAIN_ICE_PARAMN%NDRYLBDAS
    WRITE(UNIT=KLUOUT,FMT='("PKER_SDRYG(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_SDRYG(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_SDRYG (KDRYLBDAG,KDRYLBDAS,KND,                              &
                     PALPHAG,PNUG,PALPHAS,PNUS,PEGS,PBS,PCG,PDG,PCS,PDS,PFVELOS, &
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
IF( .NOT.ASSOCIATED(XKER_RDRYG) ) CALL RAIN_ICE_PARAM_ALLOCATE('XKER_RDRYG', RAIN_ICE_PARAMN%NDRYLBDAG,RAIN_ICE_PARAMN%NDRYLBDAR)
!
CALL READ_XKER_RDRYG (KDRYLBDAG,KDRYLBDAR,KND,                              &
                   PALPHAG,PNUG,PALPHAR,PNUR,PEGR,PBR,PCG,PDG,PCR,PDR,      &
                   PDRYLBDAG_MAX,PDRYLBDAR_MAX,PDRYLBDAG_MIN,PDRYLBDAR_MIN, &
                   PFDINFTY                                                 )
IF( (KDRYLBDAG/=RAIN_ICE_PARAMN%NDRYLBDAG) .OR. (KDRYLBDAR/=RAIN_ICE_PARAMN%NDRYLBDAR) .OR. (KND/=IND) .OR. &
    (PALPHAG/=XALPHAG) .OR. (PNUG/=XNUG)                               .OR. &
    (PALPHAR/=XALPHAR) .OR. (PNUR/=XNUR)                               .OR. &
    (PEGR/=ZEGR) .OR. (PBR/=XBR)                                       .OR. &
    (PCG/=XCG) .OR. (PDG/=XDG) .OR. (PCR/=XCR) .OR. (PDR/=XDR)         .OR. &
    (PDRYLBDAG_MAX/=RAIN_ICE_PARAMN%XDRYLBDAG_MAX) .OR. (PDRYLBDAR_MAX/=RAIN_ICE_PARAMN%XDRYLBDAR_MAX) .OR. &
    (PDRYLBDAG_MIN/=RAIN_ICE_PARAMN%XDRYLBDAG_MIN) .OR. (PDRYLBDAR_MIN/=RAIN_ICE_PARAMN%XDRYLBDAR_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAG, XNUG, XALPHAR, XNUR,                          &
                ZEGR, XBR, XCG, XDG, 0., XCR, XDR, 0.,                      &
                RAIN_ICE_PARAMN%XDRYLBDAG_MAX, RAIN_ICE_PARAMN%XDRYLBDAR_MAX, &
                RAIN_ICE_PARAMN%XDRYLBDAG_MIN, RAIN_ICE_PARAMN%XDRYLBDAR_MIN, &
                ZFDINFTY, XKER_RDRYG                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF RDRYG KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KDRYLBDAG=",I3)') RAIN_ICE_PARAMN%NDRYLBDAG
  WRITE(UNIT=KLUOUT,FMT='("KDRYLBDAR=",I3)') RAIN_ICE_PARAMN%NDRYLBDAR
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
                                                    RAIN_ICE_PARAMN%XDRYLBDAG_MAX
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAR_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XDRYLBDAR_MAX
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAG_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XDRYLBDAG_MIN
  WRITE(UNIT=KLUOUT,FMT='("PDRYLBDAR_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XDRYLBDAR_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_RDRYG) ) THEN")')
  DO J1 = 1 , RAIN_ICE_PARAMN%NDRYLBDAG
    DO J2 = 1 , RAIN_ICE_PARAMN%NDRYLBDAR
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
         
!          8.2.6 Constants for possible modifying some processes related to 
!                graupeln in XFRMIN(1:8),  IN - concentration in XFRMIN(9) and Kogan 
!                autoconversion in XFRMIN(10:11). May be used for e.g. ensemble spread
  XFRMIN=XFRMIN_NAM
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
XCOLIH   = 0.01 ! Collection efficiency of I+H
XCOLEXIH = 0.1  ! Temperature factor of the I+H collection efficiency
XFWETH = (XPI/4.0)*XCCH*XCH*(ZRHO00**XCEXVT)*MOMG(XALPHAH,XNUH,XDH+2.0)
!
!*       9.2.2  Constants for the aggregate collection by the hailstones
!
XCOLSH   = 0.01 ! Collection efficiency of S+H
XCOLEXSH = 0.1  ! Temperature factor of the S+H collection efficiency
#ifdef REPRO48
XFSWETH = (XPI/4.0)*XCCH*XCCS*XAS*(ZRHO00**XCEXVT)
#else
XFSWETH = XNS*(XPI/4.0)*XCCH*XAS*(ZRHO00**XCEXVT) ! Wurtz
#endif
!
XLBSWETH1   =    MOMG(XALPHAH,XNUH,2.)*MOMG(XALPHAS,XNUS,XBS)
XLBSWETH2   = 2.*MOMG(XALPHAH,XNUH,1.)*MOMG(XALPHAS,XNUS,XBS+1.)
XLBSWETH3   =                          MOMG(XALPHAS,XNUS,XBS+2.)
!
!*       9.2.3  Constants for the graupel collection by the hailstones
!
XCOLGH   = 0.01 ! Collection efficiency of G+H
XCOLEXGH = 0.1  ! Temperature factor of the G+H collection efficiency
XFGWETH = (XPI/4.0)*XCCH*XCCG*XAG*(ZRHO00**XCEXVT)
!
XLBGWETH1   =    MOMG(XALPHAH,XNUH,2.)*MOMG(XALPHAG,XNUG,XBG)
XLBGWETH2   = 2.*MOMG(XALPHAH,XNUH,1.)*MOMG(XALPHAG,XNUG,XBG+1.)
XLBGWETH3   =                          MOMG(XALPHAG,XNUG,XBG+2.)
!
!*       9.2.3 bis Constants for the rain collection by the hailstones
!
XFRWETH = (XPI/4.0)*XCCH*XCCR*XAR*(ZRHO00**XCEXVT)
!
XLBRWETH1   =    MOMG(XALPHAH,XNUH,2.)*MOMG(XALPHAR,XNUR,XBR)
XLBRWETH2   = 2.*MOMG(XALPHAH,XNUH,1.)*MOMG(XALPHAR,XNUR,XBR+1.)
XLBRWETH3   =                          MOMG(XALPHAR,XNUR,XBR+2.)
!
! Notice: One magnitude of lambda discretized over 10 points
!
RAIN_ICE_PARAMN%NWETLBDAS = 80
RAIN_ICE_PARAMN%XWETLBDAS_MIN = 2.5E1 ! Minimal value of Lbda_s to tabulate XKER_SWETH
RAIN_ICE_PARAMN%XWETLBDAS_MAX = 2.5E9 ! Maximal value of Lbda_s to tabulate XKER_SWETH
ZRATE = LOG(RAIN_ICE_PARAMN%XWETLBDAS_MAX/RAIN_ICE_PARAMN%XWETLBDAS_MIN)/REAL(RAIN_ICE_PARAMN%NWETLBDAS-1)
RAIN_ICE_PARAMN%XWETINTP1S = 1.0 / ZRATE
RAIN_ICE_PARAMN%XWETINTP2S = 1.0 - LOG( RAIN_ICE_PARAMN%XWETLBDAS_MIN ) / ZRATE
RAIN_ICE_PARAMN%NWETLBDAG = 40
RAIN_ICE_PARAMN%XWETLBDAG_MIN = 1.0E3 ! Min value of Lbda_g to tabulate XKER_GWETH
RAIN_ICE_PARAMN%XWETLBDAG_MAX = 1.0E7 ! Max value of Lbda_g to tabulate XKER_GWETH
ZRATE = LOG(RAIN_ICE_PARAMN%XWETLBDAG_MAX/RAIN_ICE_PARAMN%XWETLBDAG_MIN)/REAL(RAIN_ICE_PARAMN%NWETLBDAG-1)
RAIN_ICE_PARAMN%XWETINTP1G = 1.0 / ZRATE
RAIN_ICE_PARAMN%XWETINTP2G = 1.0 - LOG( RAIN_ICE_PARAMN%XWETLBDAG_MIN ) / ZRATE
RAIN_ICE_PARAMN%NWETLBDAR = 40
RAIN_ICE_PARAMN%XWETLBDAR_MIN = 1.0E3 ! Minimal value of Lbda_r to tabulate XKER_RWETH
RAIN_ICE_PARAMN%XWETLBDAR_MAX = 1.0E7 ! Maximal value of Lbda_r to tabulate XKER_RWETH
ZRATE = LOG(RAIN_ICE_PARAMN%XWETLBDAR_MAX/RAIN_ICE_PARAMN%XWETLBDAR_MIN)/REAL(RAIN_ICE_PARAMN%NWETLBDAR-1)
RAIN_ICE_PARAMN%XWETINTP1R = 1.0 / ZRATE
RAIN_ICE_PARAMN%XWETINTP2R = 1.0 - LOG( RAIN_ICE_PARAMN%XWETLBDAR_MIN ) / ZRATE
RAIN_ICE_PARAMN%NWETLBDAH = 40
RAIN_ICE_PARAMN%XWETLBDAH_MIN = 1.0E3 ! Min value of Lbda_h to tabulate XKER_SWETH,XKER_GWETH,XKER_RWETH
RAIN_ICE_PARAMN%XWETLBDAH_MAX = 1.0E7 ! Max value of Lbda_h to tabulate XKER_SWETH,XKER_GWETH,XKER_RWETH
ZRATE = LOG(RAIN_ICE_PARAMN%XWETLBDAH_MAX/RAIN_ICE_PARAMN%XWETLBDAH_MIN)/REAL(RAIN_ICE_PARAMN%NWETLBDAH-1)
RAIN_ICE_PARAMN%XWETINTP1H = 1.0 / ZRATE
RAIN_ICE_PARAMN%XWETINTP2H = 1.0 - LOG( RAIN_ICE_PARAMN%XWETLBDAH_MIN ) / ZRATE
!
!*       9.2.4  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZEHS     = 1.0   ! factor used to integrate the dimensional distributions when
ZFDINFTY = 20.0  ! computing the kernels XKER_SWETH
!
IF( .NOT.ASSOCIATED(XKER_SWETH) ) CALL RAIN_ICE_PARAM_ALLOCATE('XKER_SWETH', RAIN_ICE_PARAMN%NWETLBDAH,RAIN_ICE_PARAMN%NWETLBDAS)
!
CALL READ_XKER_SWETH (KWETLBDAH,KWETLBDAS,KND,                              &
                   PALPHAH,PNUH,PALPHAS,PNUS,PEHS,PBS,PCH,PDH,PCS,PDS,PFVELOS, &
                   PWETLBDAH_MAX,PWETLBDAS_MAX,PWETLBDAH_MIN,PWETLBDAS_MIN, &
                   PFDINFTY                                                 )
IF( (KWETLBDAH/=RAIN_ICE_PARAMN%NWETLBDAH) .OR. (KWETLBDAS/=RAIN_ICE_PARAMN%NWETLBDAS) .OR. (KND/=IND) .OR. &
    (PALPHAH/=XALPHAH) .OR. (PNUH/=XNUH)                               .OR. &
    (PALPHAS/=XALPHAS) .OR. (PNUS/=XNUS)                               .OR. &
    (PEHS/=ZEHS) .OR. (PBS/=XBS)                                       .OR. &
    (PCH/=XCH) .OR. (PDH/=XDH) .OR. (PCS/=XCS) .OR. (PDS/=XDS) .OR. (PFVELOS/=XFVELOS) .OR. &
    (PWETLBDAH_MAX/=RAIN_ICE_PARAMN%XWETLBDAH_MAX) .OR. (PWETLBDAS_MAX/=RAIN_ICE_PARAMN%XWETLBDAS_MAX) .OR. &
    (PWETLBDAH_MIN/=RAIN_ICE_PARAMN%XWETLBDAH_MIN) .OR. (PWETLBDAS_MIN/=RAIN_ICE_PARAMN%XWETLBDAS_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAH, XNUH, XALPHAS, XNUS,                          &
                ZEHS, XBS, XCH, XDH, 0., XCS, XDS, XFVELOS,                 &
                RAIN_ICE_PARAMN%XWETLBDAH_MAX, RAIN_ICE_PARAMN%XWETLBDAS_MAX, &
                RAIN_ICE_PARAMN%XWETLBDAH_MIN, RAIN_ICE_PARAMN%XWETLBDAS_MIN, &
                ZFDINFTY, XKER_SWETH                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF SWETH KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAH=",I3)') RAIN_ICE_PARAMN%NWETLBDAH
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAS=",I3)') RAIN_ICE_PARAMN%NWETLBDAS
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
  WRITE(UNIT=KLUOUT,FMT='("PFVELOS=",E13.6)') XFVELOS
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAH_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAS_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAS_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAH_MIN
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAS_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAS_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_SWETH) ) THEN")')
  DO J1 = 1 , RAIN_ICE_PARAMN%NWETLBDAH
    DO J2 = 1 , RAIN_ICE_PARAMN%NWETLBDAS
    WRITE(UNIT=KLUOUT,FMT='("PKER_SWETH(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_SWETH(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_SWETH (KWETLBDAH,KWETLBDAS,KND,                              &
                     PALPHAH,PNUH,PALPHAS,PNUS,PEHS,PBS,PCH,PDH,PCS,PDS,PFVELOS, &
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
IF( .NOT.ASSOCIATED(XKER_GWETH) ) CALL RAIN_ICE_PARAM_ALLOCATE('XKER_GWETH', RAIN_ICE_PARAMN%NWETLBDAH,RAIN_ICE_PARAMN%NWETLBDAG)
!
CALL READ_XKER_GWETH (KWETLBDAH,KWETLBDAG,KND,                              &
                   PALPHAH,PNUH,PALPHAG,PNUG,PEHG,PBG,PCH,PDH,PCG,PDG,      &
                   PWETLBDAH_MAX,PWETLBDAG_MAX,PWETLBDAH_MIN,PWETLBDAG_MIN, &
                   PFDINFTY                                                 )
IF( (KWETLBDAH/=RAIN_ICE_PARAMN%NWETLBDAH) .OR. (KWETLBDAG/=RAIN_ICE_PARAMN%NWETLBDAG) .OR. (KND/=IND) .OR. &
    (PALPHAH/=XALPHAH) .OR. (PNUH/=XNUH)                               .OR. &
    (PALPHAG/=XALPHAG) .OR. (PNUG/=XNUG)                               .OR. &
    (PEHG/=ZEHG) .OR. (PBG/=XBG)                                       .OR. &
    (PCH/=XCH) .OR. (PDH/=XDH) .OR. (PCG/=XCG) .OR. (PDG/=XDG)         .OR. &
    (PWETLBDAH_MAX/=RAIN_ICE_PARAMN%XWETLBDAH_MAX) .OR. (PWETLBDAG_MAX/=RAIN_ICE_PARAMN%XWETLBDAG_MAX) .OR. &
    (PWETLBDAH_MIN/=RAIN_ICE_PARAMN%XWETLBDAH_MIN) .OR. (PWETLBDAG_MIN/=RAIN_ICE_PARAMN%XWETLBDAG_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAH, XNUH, XALPHAG, XNUG,                          &
                ZEHG, XBG, XCH, XDH, 0., XCG, XDG, 0.,                      &
                RAIN_ICE_PARAMN%XWETLBDAH_MAX, RAIN_ICE_PARAMN%XWETLBDAG_MAX, &
                RAIN_ICE_PARAMN%XWETLBDAH_MIN, RAIN_ICE_PARAMN%XWETLBDAG_MIN, &
                ZFDINFTY, XKER_GWETH                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF GWETH KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAH=",I3)') RAIN_ICE_PARAMN%NWETLBDAH
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAG=",I3)') RAIN_ICE_PARAMN%NWETLBDAG
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
                                                    RAIN_ICE_PARAMN%XWETLBDAH_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAG_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAG_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAH_MIN
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAG_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAG_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_GWETH) ) THEN")')
  DO J1 = 1 , RAIN_ICE_PARAMN%NWETLBDAH
    DO J2 = 1 , RAIN_ICE_PARAMN%NWETLBDAG
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
IND      = 50    ! Number of interval used to integrate the dimensional
ZEHR     = 1.0   ! distributions when computing the kernel XKER_RWETH
ZFDINFTY = 20.0
!
IF( .NOT.ASSOCIATED(XKER_RWETH) ) CALL RAIN_ICE_PARAM_ALLOCATE('XKER_RWETH', RAIN_ICE_PARAMN%NWETLBDAH,RAIN_ICE_PARAMN%NWETLBDAR)
!
CALL READ_XKER_RWETH (KWETLBDAH,KWETLBDAR,KND,                              &
                   PALPHAH,PNUH,PALPHAR,PNUR,PEHR,PBR,PCH,PDH,PCR,PDR,      &
                   PWETLBDAH_MAX,PWETLBDAR_MAX,PWETLBDAH_MIN,PWETLBDAR_MIN, &
                   PFDINFTY                                                 )
IF( (KWETLBDAH/=RAIN_ICE_PARAMN%NWETLBDAH) .OR. (KWETLBDAR/=RAIN_ICE_PARAMN%NWETLBDAR) .OR. (KND/=IND) .OR. &
    (PALPHAH/=XALPHAH) .OR. (PNUH/=XNUH)                               .OR. &
    (PALPHAR/=XALPHAR) .OR. (PNUR/=XNUR)                               .OR. &
    (PEHR/=ZEHR) .OR. (PBR/=XBR)                                       .OR. &
    (PCH/=XCH) .OR. (PDH/=XDH) .OR. (PCR/=XCR) .OR. (PDR/=XDR)         .OR. &
    (PWETLBDAH_MAX/=RAIN_ICE_PARAMN%XWETLBDAH_MAX) .OR. (PWETLBDAR_MAX/=RAIN_ICE_PARAMN%XWETLBDAR_MAX) .OR. &
    (PWETLBDAH_MIN/=RAIN_ICE_PARAMN%XWETLBDAH_MIN) .OR. (PWETLBDAR_MIN/=RAIN_ICE_PARAMN%XWETLBDAR_MIN) .OR. &
    (PFDINFTY/=ZFDINFTY)                                               ) THEN
  CALL RZCOLX ( IND, XALPHAH, XNUH, XALPHAR, XNUR,                          &
                ZEHR, XBR, XCH, XDH, 0., XCR, XDR, 0.,                      &
                RAIN_ICE_PARAMN%XWETLBDAH_MAX, RAIN_ICE_PARAMN%XWETLBDAR_MAX, &
                RAIN_ICE_PARAMN%XWETLBDAH_MIN, RAIN_ICE_PARAMN%XWETLBDAR_MIN, &
                ZFDINFTY, XKER_RWETH                                        )
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("**** UPDATE NEW SET OF RWETH KERNELS ****")')
  WRITE(UNIT=KLUOUT,FMT='("*****************************************")')
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("KND=",I3)') IND
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAH=",I3)') RAIN_ICE_PARAMN%NWETLBDAH
  WRITE(UNIT=KLUOUT,FMT='("KWETLBDAR=",I3)') RAIN_ICE_PARAMN%NWETLBDAR
  WRITE(UNIT=KLUOUT,FMT='("PALPHAH=",E13.6)') XALPHAH
  WRITE(UNIT=KLUOUT,FMT='("PNUH=",E13.6)') XNUH
  WRITE(UNIT=KLUOUT,FMT='("PALPHAR=",E13.6)') XALPHAR
  WRITE(UNIT=KLUOUT,FMT='("PNUR=",E13.6)') XNUR
  WRITE(UNIT=KLUOUT,FMT='("PEHR=",E13.6)') ZEHR
  WRITE(UNIT=KLUOUT,FMT='("PBR=",E13.6)') XBR
  WRITE(UNIT=KLUOUT,FMT='("PCH=",E13.6)') XCH
  WRITE(UNIT=KLUOUT,FMT='("PDH=",E13.6)') XDH
  WRITE(UNIT=KLUOUT,FMT='("PCR=",E13.6)') XCR
  WRITE(UNIT=KLUOUT,FMT='("PDR=",E13.6)') XDR
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAH_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAR_MAX=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAR_MAX
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAH_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAH_MIN
  WRITE(UNIT=KLUOUT,FMT='("PWETLBDAR_MIN=",E13.6)') &
                                                    RAIN_ICE_PARAMN%XWETLBDAR_MIN
  WRITE(UNIT=KLUOUT,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=KLUOUT,FMT='("!")')
  WRITE(UNIT=KLUOUT,FMT='("IF( PRESENT(PKER_RWETH) ) THEN")')
  DO J1 = 1 , RAIN_ICE_PARAMN%NWETLBDAH
    DO J2 = 1 , RAIN_ICE_PARAMN%NWETLBDAR
    WRITE(UNIT=KLUOUT,FMT='("PKER_RWETH(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RWETH(J1,J2)
    END DO
  END DO
  WRITE(UNIT=KLUOUT,FMT='("END IF")')
  ELSE
  CALL READ_XKER_RWETH (KWETLBDAH,KWETLBDAR,KND,                              &
                     PALPHAH,PNUH,PALPHAR,PNUR,PEHR,PBR,PCH,PDH,PCR,PDR,      &
                     PWETLBDAH_MAX,PWETLBDAR_MAX,PWETLBDAH_MIN,PWETLBDAR_MIN, &
                     PFDINFTY,XKER_RWETH                                      )
  WRITE(UNIT=KLUOUT,FMT='(" Read XKER_RWETH")')
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
IF (LHOOK) CALL DR_HOOK('INI_RAIN_ICE',1,ZHOOK_HANDLE)
CONTAINS
!
!------------------------------------------------------------------------------
!
  FUNCTION MOMG(PALPHA,PNU,PP) RESULT (PMOMG)
!
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law
!
  USE MODI_GAMMA
!
  IMPLICIT NONE
!
  REAL, INTENT(IN)     :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL, INTENT(IN)     :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL, INTENT(IN)     :: PP     ! order of the moment
  REAL     :: PMOMG  ! result: moment of order ZP
!
!------------------------------------------------------------------------------
!
!
  PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)
!
  END FUNCTION MOMG
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INI_RAIN_ICE
END MODULE MODE_INI_RAIN_ICE
