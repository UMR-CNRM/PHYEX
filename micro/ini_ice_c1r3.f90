!MNH_LIC Copyright 2000-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ########################
       MODULE MODI_INI_ICE_C1R3 
!      ########################
!
INTERFACE
      SUBROUTINE INI_ICE_C1R3 ( PTSTEP, PDZMIN, KSPLITG )
!
INTEGER,                 INTENT(OUT):: KSPLITG   ! Number of small time step
                                                 ! integration for  rain
                                                 ! sedimendation
!
REAL,                    INTENT(IN) :: PTSTEP    ! Time step
!
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
!
END SUBROUTINE INI_ICE_C1R3
!
END INTERFACE
!
END MODULE MODI_INI_ICE_C1R3
!     ###################################################
      SUBROUTINE INI_ICE_C1R3 ( PTSTEP, PDZMIN, KSPLITG )
!     ###################################################
!
!!****  *INI_ICE_C1R3 * - initialize the constants necessary for the warm and
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
!!      Module MODD_ICE_C1R3_DESCR
!!      Module MODD_ICE_C1R3_PARAM
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine INI_ICE_C1R3 )
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/2000
!!      J.-P. Pinty 28/05/2001 Correction for RHONI
!!      J.-P. Pinty 31/05/2001 Correction for ICNVS factors
!!      J.-P. Pinty 29/06/2001 Bug in RCHONI and RVHNCI
!!      J.-P. Pinty 29/06/2001 Add RHHONI process (freezing haze part.)
!!      J.-P. Pinty 23/09/2001 Review the HM process constants
!!      J.-P. Pinty 23/10/2001 Add XRHORSMIN
!!      J.-P. Pinty 05/04/2002 Add computation of the effective radius
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  J. Wurtz       03/2022: new snow characteristics
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_ICE_C1R3_DESCR
USE MODD_ICE_C1R3_PARAM
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAMETERS
USE MODD_PARAM_C1R3
USE MODD_PARAM_C2R2,      ONLY : XALPHAC,XNUC,XALPHAR,XNUR
USE MODD_RAIN_C2R2_DESCR, ONLY : XAR,XBR,XCR,XDR,XF0R,XF1R,XAC,XBC,XCC,XDC, &
                                 XLBC,XLBEXC,XLBR,XLBEXR
USE MODD_REF
!
use mode_msg
!
USE MODD_RAIN_ICE_DESCR, ONLY : XFVELOS
!
USE MODI_GAMMA
USE MODI_GAMMA_INC
USE MODE_READ_XKER_RACCS, ONLY: READ_XKER_RACCS
USE MODE_READ_XKER_RDRYG, ONLY: READ_XKER_RDRYG
USE MODE_READ_XKER_SDRYG, ONLY: READ_XKER_SDRYG
USE MODE_RRCOLSS, ONLY: RRCOLSS
USE MODE_RSCOLRG, ONLY: RSCOLRG
USE MODE_RZCOLX,  ONLY: RZCOLX
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                 INTENT(OUT):: KSPLITG   ! Number of small time step
                                                 ! integration for  rain
                                                 ! sedimendation
!
REAL,                    INTENT(IN) :: PTSTEP    ! Time step
!
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!                                                                           diat
!
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB                ! Coordinates of the first  physical 
                              ! points along z
INTEGER :: J1,J2              ! Internal loop indexes
!
REAL, DIMENSION(8)  :: ZGAMI  ! parameters involving various moments
REAL, DIMENSION(2)  :: ZGAMS  ! of the generalized gamma law
!
REAL :: ZT                    ! Work variable
REAL :: ZVTRMAX               ! Raindrop maximal fall velocity
REAL :: ZRHO00                ! Surface reference air density
REAL :: ZRATE                 ! Geometrical growth of Lbda in the tabulated
                              ! functions and kernels
REAL :: ZBOUND                ! XDCSLIM*Lbda_s: upper bound for the partial
                              ! integration of the riming rate of the aggregates
REAL :: ZEGS, ZEGR            ! Bulk collection efficiencies
!
INTEGER :: IND                ! Number of interval to integrate the kernels
REAL :: ZESR                  ! Mean efficiency of rain-aggregate collection
REAL :: ZFDINFTY              ! Factor used to define the "infinite" diameter
!
!
INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
                    ! listing
REAL     :: ZCONC_MAX ! Maximal concentration for snow
REAL     :: ZFACT_NUCL! Amplification factor for the minimal ice concentration
!  
INTEGER  :: KND
INTEGER  :: KACCLBDAS,KACCLBDAR,KDRYLBDAG,KDRYLBDAS,KDRYLBDAR
REAL     :: PALPHAR,PALPHAS,PALPHAG
REAL     :: PNUR,PNUS,PNUG
REAL     :: PBR,PBS
REAL     :: PCR,PCS,PCG
REAL     :: PDR,PDS,PFVELOS,PDG
REAL     :: PESR,PEGS,PEGR
REAL     :: PFDINFTY
REAL     :: PACCLBDAS_MAX,PACCLBDAR_MAX,PACCLBDAS_MIN,PACCLBDAR_MIN
REAL     :: PDRYLBDAG_MAX,PDRYLBDAS_MAX,PDRYLBDAG_MIN,PDRYLBDAS_MIN
REAL     :: PDRYLBDAR_MAX,PDRYLBDAR_MIN
!
REAL     :: ZFAC_ZRNIC ! Zrnic factor used to decrease Long Kernels
!
!-------------------------------------------------------------------------------
!
!
!*       0.     FUNCTION STATEMENTS
!   	        -------------------
!
!
!*       0.1    G(p) for p_moment of the Generalized GAMMA function
!
!
! recall that MOMG(ZALPHA,ZNU,ZP)=GAMMA(ZNU+ZP/ZALPHA)/GAMMA(ZNU)
!
!
!        1.     INTIALIZE OUTPUT LISTING AND COMPUTE KSPLITG FOR EACH MODEL
!               -----------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
!*       1.1    Set the graupel maximum fall velocity
!
ZVTRMAX = 30.                          
IF( CHEVRIMED_ICE_C1R3 == 'HAIL' ) THEN
  ZVTRMAX = 60. ! Hail case 
END IF
!
!*       1.2    Compute the number of small time step integration
!
KSPLITG = 1
SPLIT : DO
  ZT = 2.* PTSTEP / REAL(KSPLITG)
  IF ( ZT * ZVTRMAX / PDZMIN .LT. 1.) EXIT SPLIT
  KSPLITG = KSPLITG + 1
END DO SPLIT
!
IF (ALLOCATED(XRTMIN)) RETURN     ! In case of nesting microphysics constants of
!                                 ! MODD_ICE_C1R3_PARAM are computed only once.
!
!-------------------------------------------------------------------------------
!
!*       2.     CHARACTERISTICS OF THE SPECIES
!   	        ------------------------------
!
!
!*       2.1    Raindrops characteristics
!
!
!*       2.2    Ice crystal characteristics
!
SELECT CASE (CPRISTINE_ICE_C1R3)
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
XF2I = 0.103
XF0IS = 0.86
XF1IS = 0.28
!
!
!*       2.3    Snowflakes/aggregates characteristics
!
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
XC1S = 1./XPI
!
!
!*       2.4    Heavily rimed crystals characteristics
!
!
SELECT CASE (CHEVRIMED_ICE_C1R3)
  CASE('GRAU')
    XAG = 19.6  ! Lump graupel case
    XBG = 2.8   ! Lump graupel case
    XCG = 124.  ! Lump graupel case
    XDG = 0.66  ! Lump graupel case
  CASE('HAIL')
    XAG = 470.  ! Hail case
    XBG = 3.0   ! Hail case
    XCG = 207.  ! Hail case
    XDG = 0.64  ! Hail case
END SELECT
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
!-------------------------------------------------------------------------------
!
!*       3.     DIMENSIONAL DISTRIBUTIONS OF THE SPECIES
!	            ----------------------------------------
!
!
!*       3.2    Ice crystal distribution
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
!*       3.3    Constants for shape parameter
!
XLBEXI = 1.0/XBI
XLBI   = XAI*MOMG(XALPHAI,XNUI,XBI)
!
XLBEXS = 1.0/(XCXS-XBS)
XLBS   = ( XAS*XCCS*MOMG(XALPHAS,XNUS,XBS) )**(-XLBEXS)
!
XLBEXG = 1.0/(XCXG-XBG)
XLBG   = ( XAG*XCCG*MOMG(XALPHAG,XNUG,XBG))**(-XLBEXG)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      Shape Parameters")')
  WRITE(UNIT=ILUOUT0,FMT='(" XLBEXI =",E13.6," XLBI =",E13.6)') XLBEXI,XLBI
  WRITE(UNIT=ILUOUT0,FMT='(" XLBEXS =",E13.6," XLBS =",E13.6)') XLBEXS,XLBS
  WRITE(UNIT=ILUOUT0,FMT='(" XLBEXG =",E13.6," XLBG =",E13.6)') XLBEXG,XLBG
END IF
!
!*       3.4    Minimal values allowed for the mixing ratios
!
XLBDAS_MAX = 100000.0
XLBDAG_MAX = 100000.0
!
ZCONC_MAX  = 1.E6 ! Maximal concentration for falling particules set to 1 per cc
XLBDAS_MAX = ( ZCONC_MAX/XCCS )**(1./XCXS) 
!
ALLOCATE( XRTMIN(6) )
XRTMIN(1) = 1.0E-20
XRTMIN(2) = 1.0E-20
XRTMIN(3) = 1.0E-20
XRTMIN(4) = 1.0E-20
XRTMIN(5) = 1.0E-15
XRTMIN(6) = 1.0E-15
ALLOCATE( XCTMIN(6) )
XCTMIN(1) = 1.0
XCTMIN(2) = 1.0
XCTMIN(3) = 1.0E-3
XCTMIN(4) = 1.0E-3
XCTMIN(5) = 1.0E-3
XCTMIN(6) = 1.0E-3
!
!-------------------------------------------------------------------------------
!
!*       4.     CONSTANTS FOR THE SEDIMENTATION
!   	        -------------------------------
!
!
!*       4.1    Exponent of the fall-speed air density correction
!
XCEXVT = 0.4
!
IKB = 1 + JPVEXT
ZRHO00 = XP00/(XRD*XTHVREFZ(IKB))
!
!*       4.2    Constants for sedimentation
!
!! XEXRSEDI = (XBI+XDI)/XBI
!! XEXCSEDI = 1.0-XEXRSEDI
!! XFSEDI   = (4.*XPI*900.)**(-XEXCSEDI) *                         &
!!            XC_I*XAI*MOMG(XALPHAI,XNUI,XBI+XDI) *                &
!!            ((XAI*MOMG(XALPHAI,XNUI,XBI)))**(-XEXRSEDI) *        &
!!            (ZRHO00)**XCEXVT
!! !
!! !  Computations made for Columns
!! !
!! XEXRSEDI = 1.9324
!! XEXCSEDI =-0.9324
!! XFSEDI   = 3.89745E11*MOMG(XALPHAI,XNUI,3.285)*                          &
!!                       MOMG(XALPHAI,XNUI,1.7)**(-XEXRSEDI)*(ZRHO00)**XCEXVT
!! XEXCSEDI =-0.9324*3.0
!! WRITE (ILUOUT0,FMT=*)' PRISTINE ICE SEDIMENTATION for columns XFSEDI=',XFSEDI
!
!
XEXSEDS = (XBS+XDS-XCXS)/(XBS-XCXS)
XFSEDS  = XCS*XAS*XCCS*MOMG(XALPHAS,XNUS,XBS+XDS)*                         &
            (XAS*XCCS*MOMG(XALPHAS,XNUS,XBS))**(-XEXSEDS)*(ZRHO00)**XCEXVT
!
XEXSEDG = (XBG+XDG-XCXG)/(XBG-XCXG)
XFSEDG  = XCG*XAG*XCCG*MOMG(XALPHAG,XNUG,XBG+XDG)*                         &
            (XAG*XCCG*MOMG(XALPHAG,XNUG,XBG))**(-XEXSEDG)*(ZRHO00)**XCEXVT
!
!
!-------------------------------------------------------------------------------
!
!*       5.     CONSTANTS FOR THE SLOW COLD PROCESSES
!      	        -------------------------------------
!
!
!*       5.1    Constants for ice nucleation
!
SELECT CASE (CPRISTINE_ICE_C1R3)
  CASE('PLAT')
    ZFACT_NUCL =  1.0    ! Plates
  CASE('COLU')
    ZFACT_NUCL = 25.0    ! Columns
  CASE('BURO')
    ZFACT_NUCL = 17.0    ! Bullet rosettes
END SELECT
!
!*       5.1.1  Constants for nucleation from ice nuclei
!
XNUC_DEP  = XFACTNUC_DEP*1000.*ZFACT_NUCL
XEXSI_DEP = 12.96
XEX_DEP   = -0.639
!
XCONCI_MAX = 100.E3  !  Assume a maximum concentration of 100 per liter
XNUC_CON   = XFACTNUC_CON*1000.*ZFACT_NUCL
XEXTT_CON  = -0.262
XEX_CON    = -2.8
!
XMNU0 = 6.88E-13
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      Heterogeneous nucleation")')
  WRITE(UNIT=ILUOUT0,FMT='(" XNUC_DEP=",E13.6," XEXSI=",E13.6," XEX=",E13.6)') &
                                                      XNUC_DEP,XEXSI_DEP,XEX_DEP
  WRITE(UNIT=ILUOUT0,FMT='(" XNUC_CON=",E13.6," XEXTT=",E13.6," XEX=",E13.6)') &
                                                      XNUC_CON,XEXTT_CON,XEX_CON
  WRITE(UNIT=ILUOUT0,FMT='(" mass of embryo XMNU0=",E13.6)') XMNU0
END IF
!
!*       5.1.2  Constants for homogeneous nucleation from haze particules
!
XRHOI_HONH = 925.0
XCEXP_DIFVAP_HONH = 1.94
XCOEF_DIFVAP_HONH = (2.0*XPI)*0.211E-4*XP00/XTT**XCEXP_DIFVAP_HONH
XCRITSAT1_HONH = 2.583
XCRITSAT2_HONH = 207.83
XTMIN_HONH = 180.0
XTMAX_HONH = 240.0
XDLNJODT1_HONH = 4.37
XDLNJODT2_HONH = 0.03
XC1_HONH = 100.0
XC2_HONH = 22.6
XC3_HONH = 0.1
XRCOEF_HONH = (XPI/6.0)*XRHOI_HONH
!
!*       5.1.3  Constants for homogeneous nucleation from cloud droplets
!
XTEXP1_HONC = -606.3952*LOG(10.0)
XTEXP2_HONC =  -52.6611*LOG(10.0)
XTEXP3_HONC =   -1.7439*LOG(10.0)
XTEXP4_HONC =   -0.0265*LOG(10.0)
XTEXP5_HONC = -1.536E-4*LOG(10.0)
IF (XALPHAC == 3.0) THEN
  XC_HONC   = XPI/6.0
  XR_HONC   = XPI/6.0
ELSE
  WRITE(UNIT=ILUOUT0,FMT='("      Homogeneous nucleation")')
  WRITE(UNIT=ILUOUT0,FMT='(" XALPHAC=",E13.6," IS NOT 3.0")') XALPHAC
  WRITE(UNIT=ILUOUT0,FMT='(" No algorithm yet developed in this case !")')
  call Print_msg(NVERB_FATAL,'GEN','INI_ICE_C1R3','')
END IF
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      Homogeneous nucleation")')
  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP1_HONC=",E13.6)') XTEXP1_HONC
  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP2_HONC=",E13.6)') XTEXP2_HONC
  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP3_HONC=",E13.6)') XTEXP3_HONC
  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP4_HONC=",E13.6)') XTEXP4_HONC
  WRITE(UNIT=ILUOUT0,FMT='(" XTEXP5_HONC=",E13.6)') XTEXP5_HONC
  WRITE(UNIT=ILUOUT0,FMT='("XC_HONC=",E13.6," XR_HONC=",E13.6)') XC_HONC,XR_HONC
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
! Harrington parameterization for ice to snow conversion
!
XDICNVS_LIM = 125.E-6  ! size in microns
XLBDAICNVS_LIM = (50.0**(1.0/(XALPHAI)))/XDICNVS_LIM  ! ZLBDAI Limitation
XC0DEPIS = ((4.0*XPI)/(XAI*XBI))*XC1I*XF0IS*                        &
           (XALPHAI/GAMMA(XNUI))*XDICNVS_LIM**(1.0-XBI)
XC1DEPIS = ((4.0*XPI)/(XAI*XBI))*XC1I*XF1IS*SQRT(XC_I)*             &
           (XALPHAI/GAMMA(XNUI))*XDICNVS_LIM**(1.0-XBI+(XDI+1.0)/2.0)
XR0DEPIS = XC0DEPIS *(XAI*XDICNVS_LIM**XBI)
XR1DEPIS = XC1DEPIS *(XAI*XDICNVS_LIM**XBI)
!
! Harrington parameterization for snow to ice conversion
!
XLBDASCNVI_MAX = 6000. ! lbdas max after Field (1999)
XCSCNVI_MAX    = 1000. ! estimated ice conc. due to S->I conversion
XRHORSMIN      = (XLBDASCNVI_MAX/XLBS)**(1.0/XLBEXS)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='(" snow is converted into pristine ice with ")')
  WRITE(UNIT=ILUOUT0,FMT='(" XRHORSMIN=",E13.6)') XRHORSMIN
END IF
!
XDSCNVI_LIM    = 125.E-6  ! size in microns
XLBDASCNVI_LIM = (50.0**(1.0/(XALPHAS)))/XDSCNVI_LIM  ! ZLBDAS Limitation
XC0DEPSI = ((4.0*XPI)/(XAS*XBS))*XC1S*XF0IS*                        &
           (XALPHAS/GAMMA(XNUS))*XDSCNVI_LIM**(1.0-XBS)
XC1DEPSI = ((4.0*XPI)/(XAS*XBS))*XC1S*XF1IS*SQRT(XCS)*             &
           (XALPHAS/GAMMA(XNUS))*XDSCNVI_LIM**(1.0-XBS+(XDS+1.0)/2.0)
XR0DEPSI = XC0DEPSI *(XAS*XDSCNVI_LIM**XBS)
XR1DEPSI = XC1DEPSI *(XAS*XDSCNVI_LIM**XBS)
!
! Vapor deposition on the snow and the graupels
!
X0DEPS = (4.0*XPI)*XCCS*XC1S*XF0S*MOMG(XALPHAS,XNUS,1.)
X1DEPS = (4.0*XPI)*XCCS*XC1S*XF1S*SQRT(XCS)*MOMG(XALPHAS,XNUS,0.5*XDS+1.5)
XEX0DEPS = XCXS-1.0
XEX1DEPS = XCXS-0.5*(XDS+3.0)
!
X0DEPG = (4.0*XPI)*XCCG*XC1G*XF0G*MOMG(XALPHAG,XNUG,1.)
X1DEPG = (4.0*XPI)*XCCG*XC1G*XF1G*SQRT(XCG)*MOMG(XALPHAG,XNUG,0.5*XDG+1.5)
XEX0DEPG = XCXG-1.0
XEX1DEPG = XCXG-0.5*(XDG+3.0)
!
!-------------------------------------------------------------------------------
!
!*       6.     CONSTANTS FOR THE COALESCENCE PROCESSES
!               --------------------------------------
!
!
!*       6.0    Precalculation of the gamma function momentum
!
!
ZGAMI(1) = GAMMA(XNUI)
ZGAMI(2) = MOMG(XALPHAI,XNUI,3.)
ZGAMI(3) = MOMG(XALPHAI,XNUI,6.)
ZGAMI(4) = ZGAMI(3)-ZGAMI(2)**2  ! useful for Sig_I
ZGAMI(5) = MOMG(XALPHAI,XNUI,9.)
ZGAMI(6) = MOMG(XALPHAI,XNUI,3.+XBI)
ZGAMI(7) = MOMG(XALPHAI,XNUI,XBI)
ZGAMI(8) = MOMG(XALPHAI,XNUI,3.)/MOMG(XALPHAI,XNUI,2.)
!
ZGAMS(1) = GAMMA(XNUS)
ZGAMS(2) = MOMG(XALPHAS,XNUS,3.)
!
!*       6.1    Csts for the coalescence processes
!
ZFAC_ZRNIC = 0.1
XKER_ZRNIC_A1 = 2.59E15*ZFAC_ZRNIC**2! From Long  a1=9.44E9 cm-3 
                                     ! so XKERA1= 9.44E9*1E6*(PI/6)**2
XKER_ZRNIC_A2 = 3.03E3*ZFAC_ZRNIC    ! From Long  a2=5.78E3      
                                     ! so XKERA2= 5.78E3*    (PI/6)
!
!*       6.2    Csts for the pristine ice selfcollection process
!
XSELFI = XKER_ZRNIC_A1*ZGAMI(3)
XCOLEXII = 0.025   !  Temperature factor of the I+I collection efficiency
!
!*       6.3    Constants for pristine ice autoconversion
!
XTEXAUTI = 0.025   !  Temperature factor of the I+I collection efficiency
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      pristine ice autoconversion")')
  WRITE(UNIT=ILUOUT0,FMT='(" Temp. factor    XTEXAUTI=",E13.6)') XTEXAUTI
END IF
!
XAUTO3 = 6.25E18*(ZGAMI(2))**(1./3.)*SQRT(ZGAMI(4))
XAUTO4 = 0.5E6*(ZGAMI(4))**(1./6.)
XLAUTS = 2.7E-2
XLAUTS_THRESHOLD  = 0.4
XITAUTS= 0.27 ! (Notice that T2 of BR74 is uncorrect and that 0.27=1./3.7
XITAUTS_THRESHOLD = 7.5
!
!*       6.4    Constants for snow aggregation
!
XCOLEXIS = 0.05    ! Temperature factor of the I+S collection efficiency
XAGGS_CLARGE1 = XKER_ZRNIC_A2*ZGAMI(2)
XAGGS_CLARGE2 = XKER_ZRNIC_A2*ZGAMS(2)
XAGGS_RLARGE1 = XKER_ZRNIC_A2*ZGAMI(6)*XAI
XAGGS_RLARGE2 = XKER_ZRNIC_A2*ZGAMI(7)*ZGAMS(2)*XAI
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      snow aggregation")')
  WRITE(UNIT=ILUOUT0,FMT='(" Temp. factor     XCOLEXIS=",E13.6)') XCOLEXIS
END IF
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
XEXCRIMSS= XCXS-XDS-2.0
XCRIMSS  = (XPI/4.0)*XCOLCS*XCCS*XCS*(ZRHO00**XCEXVT)*MOMG(XALPHAS,XNUS,XDS+2.0)
XEXCRIMSG= XEXCRIMSS
XCRIMSG  = XCRIMSS
XSRIMCG  = XCCS*XAS*MOMG(XALPHAS,XNUS,XBS)
XEXSRIMCG= XCXS-XBS
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      riming of the aggregates")')
  WRITE(UNIT=ILUOUT0,FMT='(" D_cs^lim (Farley et al.) XDCSLIM=",E13.6)') XDCSLIM
  WRITE(UNIT=ILUOUT0,FMT='(" Coll. efficiency          XCOLCS=",E13.6)') XCOLCS
END IF
!
NGAMINC = 80
XGAMINC_BOUND_MIN = 1.0E-1 ! Minimal value of (Lbda * D_cs^lim)**alpha
XGAMINC_BOUND_MAX = 1.0E7  ! Maximal value of (Lbda * D_cs^lim)**alpha
ZRATE = EXP(LOG(XGAMINC_BOUND_MAX/XGAMINC_BOUND_MIN)/REAL(NGAMINC-1))
!
ALLOCATE( XGAMINC_RIM1(NGAMINC) )
ALLOCATE( XGAMINC_RIM2(NGAMINC) )
!
DO J1=1,NGAMINC
  ZBOUND = XGAMINC_BOUND_MIN*ZRATE**(J1-1)
  XGAMINC_RIM1(J1) = GAMMA_INC(XNUS+(2.0+XDS)/XALPHAS,ZBOUND)
  XGAMINC_RIM2(J1) = GAMMA_INC(XNUS+XBS/XALPHAS      ,ZBOUND)
END DO
!
XRIMINTP1 = XALPHAS / LOG(ZRATE)
XRIMINTP2 = 1.0 + XRIMINTP1*LOG( XDCSLIM/(XGAMINC_BOUND_MIN)**(1.0/XALPHAS) )
!
!*       7.1.1  Defining the constants for the Hallett-Mossop 
!               secondary ice nucleation process
!
XHMTMIN = XTT - 8.0
XHMTMAX = XTT - 3.0
XHM1 = 9.3E-3            ! Obsolete parameterization
XHM2 = 1.5E-3/LOG(10.0)  ! from Ferrier (1995)
XHM_YIELD = 5.E-3 ! A splinter is produced after the riming of 200 droplets
XHM_COLLCS= 1.0   ! Collision efficiency snow/droplet (with Dc>25 microns)
XHM_FACTS = XHM_YIELD*(XHM_COLLCS/XCOLCS)
!
! Notice: One magnitude of lambda discretized over 10 points for the droplets 
!
XGAMINC_HMC_BOUND_MIN = 1.0E-3 ! Min value of (Lbda * (12,25) microns)**alpha
XGAMINC_HMC_BOUND_MAX = 1.0E5  ! Max value of (Lbda * (12,25) microns)**alpha
ZRATE = EXP(LOG(XGAMINC_HMC_BOUND_MAX/XGAMINC_HMC_BOUND_MIN)/REAL(NGAMINC-1))
!
ALLOCATE( XGAMINC_HMC(NGAMINC) )
!
DO J1=1,NGAMINC
  ZBOUND = XGAMINC_HMC_BOUND_MIN*ZRATE**(J1-1)
  XGAMINC_HMC(J1) = GAMMA_INC(XNUC,ZBOUND)
END DO
!
XHMSINTP1 = XALPHAC / LOG(ZRATE)
XHMSINTP2 = 1.0 + XHMSINTP1*LOG( 12.E-6/(XGAMINC_HMC_BOUND_MIN)**(1.0/XALPHAC) )
XHMLINTP1 = XALPHAC / LOG(ZRATE)
XHMLINTP2 = 1.0 + XHMLINTP1*LOG( 25.E-6/(XGAMINC_HMC_BOUND_MIN)**(1.0/XALPHAC) )
!
!*       7.2    Constants for the accretion of raindrops onto aggregates
!
XFRACCSS = ((XPI**2)/24.0)*XCCS*XRHOLW*(ZRHO00**XCEXVT)
!
XLBRACCS1   =    MOMG(XALPHAS,XNUS,2.)*MOMG(XALPHAR,XNUR,3.)
XLBRACCS3   =                          MOMG(XALPHAR,XNUR,5.)
!
XFSACCRG = (XPI/4.0)*XAS*XCCS*(ZRHO00**XCEXVT)
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
NACCLBDAS = 40
XACCLBDAS_MIN = 5.0E1 ! Minimal value of Lbda_s to tabulate XKER_RACCS
XACCLBDAS_MAX = 5.0E5 ! Maximal value of Lbda_s to tabulate XKER_RACCS
ZRATE = LOG(XACCLBDAS_MAX/XACCLBDAS_MIN)/REAL(NACCLBDAS-1)
XACCINTP1S = 1.0 / ZRATE
XACCINTP2S = 1.0 - LOG( XACCLBDAS_MIN ) / ZRATE
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
ALLOCATE( XKER_RACCSS(NACCLBDAS,NACCLBDAR) )
ALLOCATE( XKER_RACCS (NACCLBDAS,NACCLBDAR) )
ALLOCATE( XKER_SACCRG(NACCLBDAR,NACCLBDAS) )
!
CALL READ_XKER_RACCS (KACCLBDAS,KACCLBDAR,KND,                                &
                      PALPHAS,PNUS,PALPHAR,PNUR,PESR,PBS,PBR,PCS,PDS,PFVELOS,PCR,PDR, &
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
                 ZESR, XBR, XCS, XDS, XFVELOS, XCR, XDR,                     &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_RACCSS, XAG, XBS, XAS                        )
  CALL RZCOLX  ( IND, XALPHAS, XNUS, XALPHAR, XNUR,                          &
                 ZESR, XBR, XCS, XDS, XFVELOS, XCR, XDR, 0.,                             &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_RACCS                                        )
  CALL RSCOLRG ( IND, XALPHAS, XNUS, XALPHAR, XNUR,                          &
                 ZESR, XBS, XCS, XDS, XFVELOS, XCR, XDR,                              &
                 XACCLBDAS_MAX, XACCLBDAR_MAX, XACCLBDAS_MIN, XACCLBDAR_MIN, &
                 ZFDINFTY, XKER_SACCRG,XAG, XBS, XAS                         )
  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF RACSS KERNELS ****")')
  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF RACS  KERNELS ****")')
  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF SACRG KERNELS ****")')
  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
  WRITE(UNIT=ILUOUT0,FMT='("!")')
  WRITE(UNIT=ILUOUT0,FMT='("KND=",I3)') IND
  WRITE(UNIT=ILUOUT0,FMT='("KACCLBDAS=",I3)') NACCLBDAS
  WRITE(UNIT=ILUOUT0,FMT='("KACCLBDAR=",I3)') NACCLBDAR
  WRITE(UNIT=ILUOUT0,FMT='("PALPHAS=",E13.6)') XALPHAS
  WRITE(UNIT=ILUOUT0,FMT='("PNUS=",E13.6)') XNUS
  WRITE(UNIT=ILUOUT0,FMT='("PALPHAR=",E13.6)') XALPHAR
  WRITE(UNIT=ILUOUT0,FMT='("PNUR=",E13.6)') XNUR
  WRITE(UNIT=ILUOUT0,FMT='("PESR=",E13.6)') ZESR
  WRITE(UNIT=ILUOUT0,FMT='("PBS=",E13.6)') XBS
  WRITE(UNIT=ILUOUT0,FMT='("PBR=",E13.6)') XBR
  WRITE(UNIT=ILUOUT0,FMT='("PCS=",E13.6)') XCS
  WRITE(UNIT=ILUOUT0,FMT='("PDS=",E13.6)') XDS
  WRITE(UNIT=ILUOUT0,FMT='("PCR=",E13.6)') XCR
  WRITE(UNIT=ILUOUT0,FMT='("PDR=",E13.6)') XDR
  WRITE(UNIT=ILUOUT0,FMT='("PACCLBDAS_MAX=",E13.6)') &
                                                    XACCLBDAS_MAX
  WRITE(UNIT=ILUOUT0,FMT='("PACCLBDAR_MAX=",E13.6)') &
                                                    XACCLBDAR_MAX
  WRITE(UNIT=ILUOUT0,FMT='("PACCLBDAS_MIN=",E13.6)') &
                                                    XACCLBDAS_MIN
  WRITE(UNIT=ILUOUT0,FMT='("PACCLBDAR_MIN=",E13.6)') &
                                                    XACCLBDAR_MIN
  WRITE(UNIT=ILUOUT0,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=ILUOUT0,FMT='("!")')
  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_RACCSS) ) THEN")')
  DO J1 = 1 , NACCLBDAS
    DO J2 = 1 , NACCLBDAR
    WRITE(UNIT=ILUOUT0,FMT='("  PKER_RACCSS(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RACCSS(J1,J2)
    END DO
  END DO
  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  WRITE(UNIT=ILUOUT0,FMT='("!")')
  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_RACCS ) ) THEN")')
  DO J1 = 1 , NACCLBDAS
    DO J2 = 1 , NACCLBDAR
    WRITE(UNIT=ILUOUT0,FMT='("  PKER_RACCS (",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RACCS (J1,J2)
    END DO
  END DO
  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  WRITE(UNIT=ILUOUT0,FMT='("!")')
  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_SACCRG) ) THEN")')
  DO J1 = 1 , NACCLBDAR
    DO J2 = 1 , NACCLBDAS
    WRITE(UNIT=ILUOUT0,FMT='("  PKER_SACCRG(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_SACCRG(J1,J2)
    END DO
  END DO
  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  ELSE
  CALL READ_XKER_RACCS (KACCLBDAS,KACCLBDAR,KND,                               &
                       PALPHAS,PNUS,PALPHAR,PNUR,PESR,PBS,PBR,PCS,PDS,PFVELOS,PCR,PDR, &
                       PACCLBDAS_MAX,PACCLBDAR_MAX,PACCLBDAS_MIN,PACCLBDAR_MIN,&
                       PFDINFTY,XKER_RACCSS,XKER_RACCS,XKER_SACCRG             )
  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_RACCSS")')
  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_RACCS ")')
  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_SACCRG")')
END IF
!
!*       7.3    Constant for the conversion-melting rate
!
XFSCVMG = 2.0
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      conversion-melting of the aggregates")')
  WRITE(UNIT=ILUOUT0,FMT='(" Conv. factor XFSCVMG=",E13.6)') XFSCVMG
END IF
!
!
!-------------------------------------------------------------------------------
!
!*       8.     CONSTANTS FOR THE FAST COLD PROCESSES FOR THE GRAUPELN
!	            ------------------------------------------------------
!
!
!*       8.1    Constants for the rain contact freezing
!
XCOLIR    = 1.0
!
! values of these coeficients differ from the single-momemt rain_ice case
!
XEXRCFRI  = -XDR-5.0
XRCFRI    = ((XPI**2)/24.0)*XRHOLW*XCOLIR*XCR*(ZRHO00**XCEXVT)     &
                                                     *MOMG(XALPHAR,XNUR,XDR+5.0)
XEXICFRR  = -XDR-2.0
XICFRR    = (XPI/4.0)*XCOLIR*XCR*(ZRHO00**XCEXVT)                  &
                                                     *MOMG(XALPHAR,XNUR,XDR+2.0)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      rain contact freezing")')
  WRITE(UNIT=ILUOUT0,FMT='(" Coll. efficiency          XCOLIR=",E13.6)') XCOLIR
END IF
!
!
!*       8.2    Constants for the dry growth of the graupeln
!
!*       8.2.1  Constants for the cloud droplet collection by the graupeln
!               and for the Hallett-Mossop process
!
XCOLCG  = 0.6  !  Estimated from Cober and List (1993)
XFCDRYG = (XPI/4.0)*XCOLCG*XCCG*XCG*(ZRHO00**XCEXVT)*MOMG(XALPHAG,XNUG,XDG+2.0)
!
XHM_COLLCG= 0.9   ! Collision efficiency graupel/droplet (with Dc>25 microns)
XHM_FACTG = XHM_YIELD*(XHM_COLLCG/XCOLCG)
!
!*       8.2.2  Constants for the cloud ice collection by the graupeln
!
XCOLIG    = 0.25 ! Collection efficiency of I+G
XCOLEXIG  = 0.05 ! Temperature factor of the I+G collection efficiency
XCOLIG   = 0.01 ! Collection efficiency of I+G
XCOLEXIG = 0.1  ! Temperature factor of the I+G collection efficiency
WRITE (ILUOUT0, FMT=*) ' NEW Constants for the cloud ice collection by the graupeln'
WRITE (ILUOUT0, FMT=*) ' XCOLIG, XCOLEXIG  = ',XCOLIG,XCOLEXIG
XFIDRYG = (XPI/4.0)*XCOLIG*XCCG*XCG*(ZRHO00**XCEXVT)*MOMG(XALPHAG,XNUG,XDG+2.0)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      cloud ice collection by the graupeln")')
  WRITE(UNIT=ILUOUT0,FMT='(" Coll. efficiency XCOLIG=",E13.6)') XCOLIG
  WRITE(UNIT=ILUOUT0,FMT='(" Temp. factor     XCOLEXIG=",E13.6)') XCOLEXIG
END IF
!
!*       8.2.3  Constants for the aggregate collection by the graupeln
!
XCOLSG    = 0.25 ! Collection efficiency of S+G
XCOLEXSG  = 0.05 ! Temperature factor of the S+G collection efficiency
XCOLSG   = 0.01 ! Collection efficiency of S+G
XCOLEXSG = 0.1  ! Temperature factor of the S+G collection efficiency
WRITE (ILUOUT0, FMT=*) ' NEW Constants for the aggregate collection by the graupeln'
WRITE (ILUOUT0, FMT=*) ' XCOLSG, XCOLEXSG  = ',XCOLSG,XCOLEXSG
XFSDRYG = (XPI/4.0)*XCOLSG*XCCG*XCCS*XAS*(ZRHO00**XCEXVT)
!
XLBSDRYG1   =    MOMG(XALPHAG,XNUG,2.)*MOMG(XALPHAS,XNUS,XBS)
XLBSDRYG2   = 2.*MOMG(XALPHAG,XNUG,1.)*MOMG(XALPHAS,XNUS,XBS+1.)
XLBSDRYG3   =                          MOMG(XALPHAS,XNUS,XBS+2.)
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='("      aggregate collection by the graupeln")')
  WRITE(UNIT=ILUOUT0,FMT='(" Coll. efficiency XCOLSG=",E13.6)') XCOLSG
  WRITE(UNIT=ILUOUT0,FMT='(" Temp. factor     XCOLEXSG=",E13.6)') XCOLEXSG
END IF
!
!*       8.2.4  Constants for the raindrop collection by the graupeln
!
XFRDRYG = ((XPI**2)/24.0)*XCCG*XRHOLW*(ZRHO00**XCEXVT)
!
XLBRDRYG1   =    MOMG(XALPHAG,XNUG,2.)*MOMG(XALPHAR,XNUR,3.)
XLBRDRYG2   = 2.*MOMG(XALPHAG,XNUG,1.)*MOMG(XALPHAR,XNUR,4.)
XLBRDRYG3   =                          MOMG(XALPHAR,XNUR,5.)
!
! Notice: One magnitude of lambda discretized over 10 points
!
NDRYLBDAR = 40
XDRYLBDAR_MIN = 1.0E3 ! Minimal value of Lbda_r to tabulate XKER_RDRYG
XDRYLBDAR_MAX = 1.0E7 ! Maximal value of Lbda_r to tabulate XKER_RDRYG
ZRATE = LOG(XDRYLBDAR_MAX/XDRYLBDAR_MIN)/REAL(NDRYLBDAR-1)
XDRYINTP1R = 1.0 / ZRATE
XDRYINTP2R = 1.0 - LOG( XDRYLBDAR_MIN ) / ZRATE
NDRYLBDAS = 80
XDRYLBDAS_MIN = 2.5E1 ! Minimal value of Lbda_s to tabulate XKER_SDRYG
XDRYLBDAS_MAX = 2.5E9 ! Maximal value of Lbda_s to tabulate XKER_SDRYG
ZRATE = LOG(XDRYLBDAS_MAX/XDRYLBDAS_MIN)/REAL(NDRYLBDAS-1)
XDRYINTP1S = 1.0 / ZRATE
XDRYINTP2S = 1.0 - LOG( XDRYLBDAS_MIN ) / ZRATE
NDRYLBDAG = 40
XDRYLBDAG_MIN = 1.0E3 ! Min value of Lbda_g to tabulate XKER_SDRYG,XKER_RDRYG
XDRYLBDAG_MAX = 1.0E7 ! Max value of Lbda_g to tabulate XKER_SDRYG,XKER_RDRYG
ZRATE = LOG(XDRYLBDAG_MAX/XDRYLBDAG_MIN)/REAL(NDRYLBDAG-1)
XDRYINTP1G = 1.0 / ZRATE
XDRYINTP2G = 1.0 - LOG( XDRYLBDAG_MIN ) / ZRATE
!
!*       8.2.5  Computations of the tabulated normalized kernels
!
IND      = 50    ! Interval number, collection efficiency and infinite diameter
ZEGS     = 1.0   ! factor used to integrate the dimensional distributions when
ZFDINFTY = 20.0  ! computing the kernels XKER_SDRYG
!
ALLOCATE( XKER_SDRYG(NDRYLBDAG,NDRYLBDAS) )
!
CALL READ_XKER_SDRYG (KDRYLBDAG,KDRYLBDAS,KND,                              &
                   PALPHAG,PNUG,PALPHAS,PNUS,PEGS,PBS,PCG,PDG,PCS,PDS,PFVELOS,      &
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
                ZEGS, XBS, XCG, XDG, 0., XCS, XDS, XFVELOS,                             &
                XDRYLBDAG_MAX, XDRYLBDAS_MAX, XDRYLBDAG_MIN, XDRYLBDAS_MIN, &
                ZFDINFTY, XKER_SDRYG                                        )
  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF SDRYG KERNELS ****")')
  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
  WRITE(UNIT=ILUOUT0,FMT='("!")')
  WRITE(UNIT=ILUOUT0,FMT='("KND=",I3)') IND
  WRITE(UNIT=ILUOUT0,FMT='("KDRYLBDAG=",I3)') NDRYLBDAG
  WRITE(UNIT=ILUOUT0,FMT='("KDRYLBDAS=",I3)') NDRYLBDAS
  WRITE(UNIT=ILUOUT0,FMT='("PALPHAG=",E13.6)') XALPHAG
  WRITE(UNIT=ILUOUT0,FMT='("PNUG=",E13.6)') XNUG
  WRITE(UNIT=ILUOUT0,FMT='("PALPHAS=",E13.6)') XALPHAS
  WRITE(UNIT=ILUOUT0,FMT='("PNUS=",E13.6)') XNUS
  WRITE(UNIT=ILUOUT0,FMT='("PEGS=",E13.6)') ZEGS
  WRITE(UNIT=ILUOUT0,FMT='("PBS=",E13.6)') XBS
  WRITE(UNIT=ILUOUT0,FMT='("PCG=",E13.6)') XCG
  WRITE(UNIT=ILUOUT0,FMT='("PDG=",E13.6)') XDG
  WRITE(UNIT=ILUOUT0,FMT='("PCS=",E13.6)') XCS
  WRITE(UNIT=ILUOUT0,FMT='("PDS=",E13.6)') XDS
  WRITE(UNIT=ILUOUT0,FMT='("PDRYLBDAG_MAX=",E13.6)') &
                                                    XDRYLBDAG_MAX
  WRITE(UNIT=ILUOUT0,FMT='("PDRYLBDAS_MAX=",E13.6)') &
                                                    XDRYLBDAS_MAX
  WRITE(UNIT=ILUOUT0,FMT='("PDRYLBDAG_MIN=",E13.6)') &
                                                    XDRYLBDAG_MIN
  WRITE(UNIT=ILUOUT0,FMT='("PDRYLBDAS_MIN=",E13.6)') &
                                                    XDRYLBDAS_MIN
  WRITE(UNIT=ILUOUT0,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=ILUOUT0,FMT='("!")')
  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_SDRYG) ) THEN")')
  DO J1 = 1 , NDRYLBDAG
    DO J2 = 1 , NDRYLBDAS
    WRITE(UNIT=ILUOUT0,FMT='("PKER_SDRYG(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_SDRYG(J1,J2)
    END DO
  END DO
  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  ELSE
  CALL READ_XKER_SDRYG (KDRYLBDAG,KDRYLBDAS,KND,                              &
                     PALPHAG,PNUG,PALPHAS,PNUS,PEGS,PBS,PCG,PDG,PCS,PDS,PFVELOS,      &
                     PDRYLBDAG_MAX,PDRYLBDAS_MAX,PDRYLBDAG_MIN,PDRYLBDAS_MIN, &
                     PFDINFTY,XKER_SDRYG                                      )
  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_SDRYG")')
END IF
!
!
IND      = 50    ! Number of interval used to integrate the dimensional
ZEGR     = 1.0   ! distributions when computing the kernel XKER_RDRYG
ZFDINFTY = 20.0
!
ALLOCATE( XKER_RDRYG(NDRYLBDAG,NDRYLBDAR) )
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
                ZEGR, XBR, XCG, XDG, 0., XCR, XDR, 0.,                      &
                XDRYLBDAG_MAX, XDRYLBDAR_MAX, XDRYLBDAG_MIN, XDRYLBDAR_MIN, &
                ZFDINFTY, XKER_RDRYG                                        )
  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
  WRITE(UNIT=ILUOUT0,FMT='("**** UPDATE NEW SET OF RDRYG KERNELS ****")')
  WRITE(UNIT=ILUOUT0,FMT='("*****************************************")')
  WRITE(UNIT=ILUOUT0,FMT='("!")')
  WRITE(UNIT=ILUOUT0,FMT='("KND=",I3)') IND
  WRITE(UNIT=ILUOUT0,FMT='("KDRYLBDAG=",I3)') NDRYLBDAG
  WRITE(UNIT=ILUOUT0,FMT='("KDRYLBDAR=",I3)') NDRYLBDAR
  WRITE(UNIT=ILUOUT0,FMT='("PALPHAG=",E13.6)') XALPHAG
  WRITE(UNIT=ILUOUT0,FMT='("PNUG=",E13.6)') XNUG
  WRITE(UNIT=ILUOUT0,FMT='("PALPHAR=",E13.6)') XALPHAR
  WRITE(UNIT=ILUOUT0,FMT='("PNUR=",E13.6)') XNUR
  WRITE(UNIT=ILUOUT0,FMT='("PEGR=",E13.6)') ZEGR
  WRITE(UNIT=ILUOUT0,FMT='("PBR=",E13.6)') XBR
  WRITE(UNIT=ILUOUT0,FMT='("PCG=",E13.6)') XCG
  WRITE(UNIT=ILUOUT0,FMT='("PDG=",E13.6)') XDG
  WRITE(UNIT=ILUOUT0,FMT='("PCR=",E13.6)') XCR
  WRITE(UNIT=ILUOUT0,FMT='("PDR=",E13.6)') XDR
  WRITE(UNIT=ILUOUT0,FMT='("PDRYLBDAG_MAX=",E13.6)') &
                                                    XDRYLBDAG_MAX
  WRITE(UNIT=ILUOUT0,FMT='("PDRYLBDAR_MAX=",E13.6)') &
                                                    XDRYLBDAR_MAX
  WRITE(UNIT=ILUOUT0,FMT='("PDRYLBDAG_MIN=",E13.6)') &
                                                    XDRYLBDAG_MIN
  WRITE(UNIT=ILUOUT0,FMT='("PDRYLBDAR_MIN=",E13.6)') &
                                                    XDRYLBDAR_MIN
  WRITE(UNIT=ILUOUT0,FMT='("PFDINFTY=",E13.6)') ZFDINFTY
  WRITE(UNIT=ILUOUT0,FMT='("!")')
  WRITE(UNIT=ILUOUT0,FMT='("IF( PRESENT(PKER_RDRYG) ) THEN")')
  DO J1 = 1 , NDRYLBDAG
    DO J2 = 1 , NDRYLBDAR
    WRITE(UNIT=ILUOUT0,FMT='("PKER_RDRYG(",I3,",",I3,") = ",E13.6)') &
                        J1,J2,XKER_RDRYG(J1,J2)
    END DO
  END DO
  WRITE(UNIT=ILUOUT0,FMT='("END IF")')
  ELSE
  CALL READ_XKER_RDRYG (KDRYLBDAG,KDRYLBDAR,KND,                              &
                     PALPHAG,PNUG,PALPHAR,PNUR,PEGR,PBR,PCG,PDG,PCR,PDR,      &
                     PDRYLBDAG_MAX,PDRYLBDAR_MAX,PDRYLBDAG_MIN,PDRYLBDAR_MIN, &
                     PFDINFTY,XKER_RDRYG                                      )
  WRITE(UNIT=ILUOUT0,FMT='(" Read XKER_RDRYG")')
END IF
!
!-------------------------------------------------------------------------------
!
!*       9.     SET-UP RADIATIVE PARAMETERS
!               ---------------------------
!
! R_eff_i = XFREFFI * (rho*r_i/N_i)**(1/3)
!
!
XFREFFI = 0.5 * ZGAMI(8) * (1.0/XLBI)**XLBEXI
!
!-------------------------------------------------------------------------------
!
!*       10.     SOME PRINTS FOR CONTROL
!                -----------------------
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='(" Summary of the ice particule characteristics")')
  WRITE(UNIT=ILUOUT0,FMT='("      PRISTINE ICE")')
  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAI,XBI
  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XC_I,XDI
  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAI,XNUI
  WRITE(UNIT=ILUOUT0,FMT='("              SNOW")')
  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAS,XBS
  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XCS,XDS
  WRITE(UNIT=ILUOUT0,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
                                                      XCCS,XCXS
  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAS,XNUS
  WRITE(UNIT=ILUOUT0,FMT='("            GRAUPEL")')
  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAG,XBG
  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XCG,XDG
  WRITE(UNIT=ILUOUT0,FMT='("           concentration:CC=",E13.6," x=",E13.6)') &
                                                      XCCG,XCXG
  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAG,XNUG
END IF
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!------------------------------------------------------------------------------
!
  FUNCTION MOMG (PALPHA,PNU,PP) RESULT (PMOMG)
!
! auxiliary routine used to compute the Pth moment order of the generalized
! gamma law
!
  USE MODI_GAMMA
!
  IMPLICIT NONE
!
  REAL     :: PALPHA ! first shape parameter of the dimensionnal distribution
  REAL     :: PNU    ! second shape parameter of the dimensionnal distribution
  REAL     :: PP     ! order of the moment
  REAL     :: PMOMG  ! result: moment of order ZP
!
!------------------------------------------------------------------------------
!
!
  PMOMG = GAMMA(PNU+PP/PALPHA)/GAMMA(PNU)
!
  END FUNCTION MOMG
!
!------------------------------------------------------------------------------
!
!
END SUBROUTINE INI_ICE_C1R3
