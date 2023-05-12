!MNH_LIC Copyright 1996-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #########################
       MODULE MODI_INI_RAIN_C2R2 
!      #########################
!
IMPLICIT NONE
INTERFACE
      SUBROUTINE INI_RAIN_C2R2 ( PTSTEP, PDZMIN, KSPLITR, HCLOUD )
IMPLICIT NONE
!
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
END SUBROUTINE INI_RAIN_C2R2
!
END INTERFACE
!
END MODULE MODI_INI_RAIN_C2R2
!     ####################################################
      SUBROUTINE INI_RAIN_C2R2 ( PTSTEP, PDZMIN, KSPLITR, HCLOUD )
!     ####################################################
!
!!****  *INI_RAIN_C2R2 * - initialize the constants for the two-moment scheme
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used in the 
!!    warm microphysical scheme C2R2. The routine allows for the choice of
!!    several activation schemes CPB, TFH and TWO. The CPB scheme can be 
!!    initialized either with a CCN shape function or directly from the 
!!    specification of aerosol properties.
!!    The cloud droplets and rain drops are assumed to follow a generalized
!!    gamma law.
!!
!!**  METHOD
!!    ------
!!      The constants are initialized to their numerical values and the number
!!    of small time step in the sedimentation scheme is computed by dividing 
!!    the 2* Deltat time interval of the leap-frog scheme so that the stability 
!!    criterion for the rain sedimentation is fulfilled for a raindrop maximal 
!!    fall velocity equal VTRMAX. 
!!
!!    EXTERNAL
!!    --------
!!      GAMMA    :  gamma function
!!      HYPGEO   :  hypergeometric function
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
!!      Module MODD_RAIN_C2R2_DESCR
!!      Module MODD_RAIN_C2R2_KHKO_PARAM
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation ( routine INI_RAIN_C2R2 )
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    31/12/96
!!      J.-P. Pinty 07/07/00  In revised form
!!      J.-P. Pinty 05/04/02  Add computation of the effective radius
!!      J.-P. Pinty 29/11/02  Add cloud doplet fall speed parameters
!!      O.Geoffroy  03/2006   Add KHKO scheme
!!      G.Delautier 09/2014   fusion MODD_RAIN_C2R2_PARAM et MODD_RAIN_KHKO_PARAM
!!      M.Mazoyer   10/2016 Constants for Droplet sedimentation adapted to fog for KHKO
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_LUNIT, ONLY: TLUOUT0
USE MODD_PARAM_C2R2
USE MODD_PARAMETERS
USE MODD_RAIN_C2R2_DESCR
USE MODD_RAIN_C2R2_KHKO_PARAM
USE MODD_REF
!
USE MODI_GAMMA
USE MODI_HYPGEO
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
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
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB                ! Coordinates of the first and last physical 
                              ! points along z
INTEGER :: J1                 ! Internal loop indexes
!
REAL, DIMENSION(6)  :: ZGAMC, ZGAMR ! parameters involving various moments of
    ! the generalized gamma law
!
REAL :: ZT                          ! Work variable
REAL :: ZTT                   ! Temperature in Celsius
REAL :: ZLV                   ! Latent heat of vaporization
REAL :: ZSS                   ! Supersaturation
REAL :: ZPSI1, ZG             ! Psi1 and G functions
REAL :: ZAHENR                ! r_star (FH92)
REAL :: ZVTRMAX               ! Raindrop maximal fall velocity
REAL :: ZRHO00                ! Surface reference air density
REAL :: ZSURF_TEN             ! Water drop surface tension
REAL :: ZSMIN, ZSMAX          ! Minimal and maximal supersaturation used to
      ! discretize the HYP functions
!
!
INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
                    ! listing
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
!        1.     INTIALIZE OUTPUT LISTING AND COMPUTE KSPLITR FOR EACH MODEL
!               -----------------------------------------------------------
!
ILUOUT0 = TLUOUT0%NLU
!
!*       1.1    Set the raindrop maximum fall velocity
!
ZVTRMAX = 30.                          
!
!*       1.2    Compute the number of small time step integration
!
KSPLITR = 1
SPLIT : DO
  ZT = PTSTEP / REAL(KSPLITR)
  IF ( ZT * ZVTRMAX / PDZMIN < 1.0) EXIT SPLIT
  KSPLITR = KSPLITR + 1
END DO SPLIT
!
IF (ALLOCATED(XRTMIN)) RETURN     ! In case of nesting microphysics constants of
!                                 ! MODD_RAIN_C2R2_KHKO_PARAM are computed only once.
!
!-------------------------------------------------------------------------------
!
!*       2.     CHARACTERISTICS OF THE SPECIES
!	        ------------------------------
!
!
!*       2.1    Cloud droplet characteristics
!
XAC = (XPI/6.0)*XRHOLW
XBC = 3.0
IF (HCLOUD=='KHKO') THEN
 XCC = XRHOLW*XG/(18.0*1.816E-5) ! Stokes flow (Pruppacher p 322 for T=293K)
ELSE
 XCC = XRHOLW*XG/(18.0*1.7E-5) ! Stokes flow (Pruppacher p 322 for T=273K)
ENDIF
XDC = 2.0
!
XF0C = 1.00
XF2C = 0.08
!
XC1C = 1./2.
!
!*       2.2    Raindrops characteristics
!
XAR = (XPI/6.0)*XRHOLW
XBR = 3.0
XCR = 842.
XDR = 0.8
!
XF0R = 0.780
XF1R = 0.265
!
!
!
!-------------------------------------------------------------------------------
!
!*       3.     DIMENSIONAL DISTRIBUTIONS OF THE SPECIES
!               ----------------------------------------
!
!*       3.1    Cloud droplet distribution
!
!XALPHAC = 3.0  ! Gamma law of the Cloud droplet (here volume-like distribution)
!XNUC    = 3.0  ! Gamma law with little dispersion
!
!*       3.2    Raindrop distribution
!
!XALPHAR = 3.0  ! Gamma law of the raindrops (here volume-like distribution)
!XNUR    = 3.0  ! Gamma law for the raindrops 
!XNUR    = 0.1
!
!*       3.3    Precalculation of the gamma function momentum
!
!
  ZGAMC(1) = GAMMA(XNUC)
  ZGAMC(2) = MOMG(XALPHAC,XNUC,3.)
  ZGAMC(3) = MOMG(XALPHAC,XNUC,6.)
  ZGAMC(4) = ZGAMC(3)-ZGAMC(2)**2  ! useful for Sig_c
  ZGAMC(5) = MOMG(XALPHAC,XNUC,9.)
  ZGAMC(6) = MOMG(XALPHAC,XNUC,3.)**(2./3.)/MOMG(XALPHAC,XNUC,2.)
!
  ZGAMR(1) = GAMMA(XNUR)
  ZGAMR(2) = MOMG(XALPHAR,XNUR,3.)
  ZGAMR(3) = MOMG(XALPHAR,XNUR,6.)
  ZGAMR(4) = MOMG(XALPHAR,XNUR,6.)
  ZGAMR(5) = MOMG(XALPHAR,XNUR,9.)
  ZGAMR(6) = MOMG(XALPHAR,XNUR,3.)**(2./3.)/MOMG(XALPHAR,XNUR,2.)
!
!
!*       3.4    Set bounds
!
ALLOCATE( XRTMIN(3) )
ALLOCATE( XCTMIN(3) )
IF (HCLOUD == 'C2R2') THEN
XRTMIN(1) = 1.0E-20
XRTMIN(2) = 1.0E-20
XRTMIN(3) = 1.0E-17
ELSE
XRTMIN(1) = 1.0E-20
XRTMIN(2) = 1.E-7
XRTMIN(3) = 1.E-8
ENDIF
!
XCTMIN(1) = 1.0
XCTMIN(2) = 1.0 
XCTMIN(3) = 1.0E-3
!
!*       3.4    Csts for the shape parameter
!
XLBC   = XAR*ZGAMC(2)
XLBEXC = 1.0/XBC
XLBR   = XAR*ZGAMR(2)
XLBEXR = 1.0/XBR
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
ZRHO00 = XP00/(XRD*XTHVREFZ(IKB))
!
!*       4.2    Constants for sedimentation
!
XFSEDRR  = XCR*GAMMA(XNUR+(XDR+3.)/XALPHAR)/GAMMA(XNUR+3./XALPHAR)*     &
    (ZRHO00)**XCEXVT
XFSEDCR  = XCR*GAMMA(XNUR+XDR/XALPHAR)/GAMMA(XNUR)*     &
     (ZRHO00)**XCEXVT
XFSEDRC  = XCC*GAMMA(XNUC+(XDC+3.)/XALPHAC)/GAMMA(XNUC+3./XALPHAC)*     &
    (ZRHO00)**XCEXVT
XFSEDCC  = XCC*GAMMA(XNUC+XDC/XALPHAC)/GAMMA(XNUC)*     &
    (ZRHO00)**XCEXVT
!
!
!-------------------------------------------------------------------------------
!
!*       5.     CONSTANTS FOR THE NUCLEATION PROCESS
!      	        -------------------------------------
!
!               
!		Compute CCN spectra parameters from CCN characteristics
!
IF (HPARAM_CCN == 'CPB' .AND. HINI_CCN == 'AER') THEN
  SELECT CASE (HTYPE_CCN)
    CASE('M')  ! NaCl maritime case
      XKHEN    = 3.251*(XLOGSIG_CCN/0.4835)**(-1.297)
      XMUHEN   = 2.589*(XLOGSIG_CCN/0.4835)**(-1.511)
      XBETAHEN = 621.689*(XR_MEAN_CCN/0.133E-6)**(3.002)      &
                        *EXP(1.081*((XLOGSIG_CCN/0.4835)-1.)) &
                        *XFSOLUB_CCN                          &
                        *(XACTEMP_CCN/290.16)**(2.995)
    CASE('C')  ! (NH4)2SO4 continental case
      XKHEN    = 1.403*(XLOGSIG_CCN/1.16)**(-1.172)
      XMUHEN   = 0.834*(XLOGSIG_CCN/1.16)**(-1.350)
      XBETAHEN = 25.499*(XR_MEAN_CCN/0.0218E-6)**(3.057)      &
                        *EXP(4.092*((XLOGSIG_CCN/1.16)-1.))   &
                        *XFSOLUB_CCN**(1.011)                 &
                        *(XACTEMP_CCN/290.16)**(3.076)
  END SELECT
  XCHEN = XCONC_CCN*(XBETAHEN**(0.5*XKHEN)*GAMMA(XMUHEN))       &
                   /(GAMMA(0.5*XKHEN+1.)*GAMMA(XMUHEN-0.5*XKHEN))
END IF
!
XWMIN = 0.01 ! Minimal positive vertical velocity required 
             ! for the activation process in Twomey and CPB scheme
XTMIN = -0.000278 ! Minimal cooling required 1K/h

!
XDIVA = 226.E-7 ! Diffusivity of water vapor in the air
XTHCO = 24.3E-3 ! Air thermal conductivity
!
!                        ( 8 Mw (Sigma)sw )3  Pi*Rho_l
!            XCSTDCRIT = ( -------------- ) * --------
!                        ( 3    Ra Rhow   )     6
!
ZSURF_TEN = 76.1E-3  ! Surface tension of a water drop at T=0 C
XCSTDCRIT = (XPI/6.)*XRHOLW*( (8.0*ZSURF_TEN )/( 3.0*XRV*XRHOLW ) )**3
!
!		Tabulation of the hypergeometric functions
!
!               F(mu,k/2, k/2+1 ,-Beta s**2) and
!               F(mu,k/2,(k+3)/2,-Beta s**2) as a function of s
!
NHYP = 200
ALLOCATE (XHYPF12(NHYP))
ALLOCATE (XHYPF32(NHYP))
!
ZSMIN = 1.0E-5  ! soit Smin=0.001 %
ZSMAX = 1.0E-1  ! soit Smax=   10 %
XHYPINTP1 = REAL(NHYP-1)/LOG(ZSMAX/ZSMIN)
XHYPINTP2 = REAL(NHYP)-XHYPINTP1*LOG(ZSMAX)
IF (HPARAM_CCN == 'CPB') THEN ! CPB98's case 
  TAB_HYP : DO J1 = 1,NHYP    ! tabulation using a logarithmic scale for the
                      ! supersaturations (0.00001<S<0.1 in "no unit")
              ZSS =ZSMAX*(ZSMIN/ZSMAX)**(REAL(NHYP-J1)/REAL(NHYP-1))
              XHYPF12(J1) = HYPGEO(XMUHEN,XKHEN/2.0,(XKHEN+2.0)/2.0,XBETAHEN, &
                                   100.*ZSS)
              XHYPF32(J1) = HYPGEO(XMUHEN,XKHEN/2.0,(XKHEN+3.0)/2.0,XBETAHEN*100**2, &
                                   ZSS)
            END DO TAB_HYP
  IF (HINI_CCN == 'CCN') THEN
    XCONC_CCN = XCHEN*(GAMMA(0.5*XKHEN+1.)*GAMMA(XMUHEN-0.5*XKHEN)) &
                     /(XBETAHEN**(0.5*XKHEN)*GAMMA(XMUHEN)) 
  END IF
ELSE                          ! other cases (but not used)
  XHYPF12(:) = 1.0
  XHYPF32(:) = 1.0
! XCONC_CCN = -1.0 ! Negative value to recall that CCN spectra, other than those
                   ! defined by CPB98, are unbounded
END IF
!
!		Compute the tabulation of function of T :
!
!                              (Psi1)**(3/2)
!               XAHENG = -----------------------
!                                G**(3/2)
!
!    		XAHENY = a2 C**p k**q            as given by Feingold
!
NAHEN = 81 ! Tabulation for each Kelvin degree in the range XTT-40 to XTT+40
XAHENINTP1 = 1.0
XAHENINTP2 = 0.5*REAL(NAHEN-1) - XTT
IF (HPARAM_CCN == 'TFH') THEN
  ALLOCATE (XAHENY(NAHEN))
  ALLOCATE (XAHENF(NAHEN))
!
!            Compute constants for the calculation of Smax.
!            XCSTHEN = 1/(rho_l 4 pi C (100)^k
!
  XCSTHEN = 1.0 / ( XRHOLW*4.0*XPI*XCHEN*(100.0)**XKHEN )
  DO J1 = 1,NAHEN
    ZTT = XTT + REAL(J1-(NAHEN-1)/2)                     ! T
    ZLV = XLVTT+(XCPV-XCL)*(ZTT-XTT)                     ! Lv
    ZPSI1 = (XG/(XRD*ZTT))*(XMV*ZLV/(XMD*XCPD*ZTT)-1.)   ! Psi1
    ZG    = 1.E-4*(6.224E-7 + 0.281E-7 * ZTT + 2.320E-10 * ZTT**2) *  & ! G
           XCHEN**(-0.127   + 2.668E-3 * ZTT + 7.583E-7  * ZTT**2) *  &
           XKHEN**(-0.214   + 9.416E-3 * ZTT - 1.173E-4  * ZTT**2)
    ZAHENR     = 1.E-2*(2.124E-3 + 3.373E-5 * ZTT + 9.632E-8 * ZTT**2) *  & ! r_star
                XCHEN**(-0.321   - 3.333E-4 * ZTT - 9.972E-6 * ZTT**2) *  &
                XKHEN**(-0.464   + 9.253E-3 * ZTT - 2.066E-5 * ZTT**2)
    XAHENF(J1) = XCSTHEN*(ZPSI1/(ZG*ZAHENR))
!
    XAHENY(J1) =  (7.128E-5 + 1.094E-6 * ZTT + 4.314E-9 * ZTT**2) *  & ! y_bar
           XCHEN**( 0.230   - 1.200E-4 * ZTT + 1.607E-5 * ZTT**2) *  &
           XKHEN**( 1.132   - 9.083E-3 * ZTT - 1.482E-5 * ZTT**2)
  END DO
!
! Additional coefficients for the dependence on W
!
  XWCOEF_F1 =-0.149
  XWCOEF_F2 = 1.514E-3
  XWCOEF_F3 = 4.375E-6
  XWCOEF_Y1 = 0.132
  XWCOEF_Y2 =-2.191E-3
  XWCOEF_Y3 = 3.934E-5
ELSE
  ALLOCATE (XAHENG(NAHEN))
  ALLOCATE (XPSI1(NAHEN))
  ALLOCATE (XPSI3(NAHEN))
!
!            Compute constants for the calculation of Smax.
!            XCSTHEN = 1/(rho_l 2 pi k C B(k/2,3/2))
!
  XCSTHEN = 1.0 / ( XRHOLW*2.0*XPI*XKHEN*XCHEN*(100.0)**XKHEN * &
                    GAMMA(XKHEN/2.0)*GAMMA(3.0/2.0)/GAMMA((XKHEN+3.0)/2.0) )
  DO J1 = 1,NAHEN
    ZTT = XTT + REAL(J1-(NAHEN-1)/2)                     ! T
    ZLV = XLVTT+(XCPV-XCL)*(ZTT-XTT)                     ! Lv
    XPSI1(J1) = (XG/(XRD*ZTT))*(XMV*ZLV/(XMD*XCPD*ZTT)-1.)   ! Psi1
    XPSI3(J1) = -1*XMV*ZLV/(XMD*XRD*(ZTT**2))   ! Psi3
    ZG    = 1./( XRHOLW*( (XRV*ZTT)/                                        & !G
                          (XDIVA*EXP(XALPW-(XBETAW/ZTT)-(XGAMW*ALOG(ZTT)))) &
                        + (ZLV/ZTT)**2/(XTHCO*XRV) ) )
    XAHENG(J1) = XCSTHEN/(ZG)**(3./2.)
  END DO
END IF
!
!
!-------------------------------------------------------------------------------
!
! Parameters used to initialise the droplet and drop concentration
! from the respective mixing ratios (used in RESTART_RAIN_C2R2)
!
! Droplet case
!
IF( HPARAM_CCN=='CPB' ) THEN
  XCONCC_INI = 0.8 * XCONC_CCN       ! 80% of the maximum CCN conc. is assumed
ELSE
  XCONCC_INI = XCHEN * (0.1)**XKHEN  ! 0.1% supersaturation is assumed
END IF
!
! Raindrop case
!
XCONCR_PARAM_INI = (1.E7)**3/(XPI*XRHOLW) ! MP law with N_O=1.E7 m-1 is assumed
!
!-------------------------------------------------------------------------------
!
!*       6.     CONSTANTS FOR THE COALESCENCE PROCESSES
!               --------------------------------------
!
!
!*       6.1    Csts for the coalescence processes
!
XKERA1 = 2.59E15   ! From Long  a1=9.44E9 cm-3 so XKERA1= 9.44E9*1E6*(PI/6)**2
XKERA2 = 3.03E3    ! From Long  a2=5.78E3      so XKERA2= 5.78E3*    (PI/6)
!
! Cst for the cloud droplet selfcollection process
!
XSELFC = XKERA1*ZGAMC(3)
!
! Cst for the autoconversion process
!
XAUTO1 = 6.25E18*(ZGAMC(2))**(1./3.)*SQRT(ZGAMC(4))
XAUTO2 = 0.5E6*(ZGAMC(4))**(1./6.)
XLAUTR = 2.7E-2
XLAUTR_THRESHOLD  = 0.4
XITAUTR= 0.27 ! (Notice that T2 of BR74 is uncorrect and that 0.27=1./3.7
XITAUTR_THRESHOLD = 7.5
XCAUTR = 3.5E9
!
! Cst for the accretion process
!
XACCR1 = ZGAMR(2)**(1./3.)
XACCR2 = 5.0E-6
XACCR3 = 12.6E-4
XACCR4 = XAUTO2
XACCR5 = 3.5
XACCR6 = 1.2*XCAUTR
XACCR_CLARGE1 = XKERA2*ZGAMC(2)
XACCR_CLARGE2 = XKERA2*ZGAMR(2)
XACCR_RLARGE1 = XKERA2*ZGAMC(3)*XRHOLW*(XPI/6.0)
XACCR_RLARGE2 = XKERA2*ZGAMC(2)*ZGAMR(2)*XRHOLW*(XPI/6.0)
XACCR_CSMALL1 = XKERA1*ZGAMC(3)
XACCR_CSMALL2 = XKERA1*ZGAMR(3)
XACCR_RSMALL1 = XKERA1*ZGAMC(5)*XRHOLW*(XPI/6.0)
XACCR_RSMALL2 = XKERA1*ZGAMC(2)*ZGAMR(3)*XRHOLW*(XPI/6.0)
!
! Cst for the raindrop self-collection/breakup process
!
XSCBU2 = XKERA2*ZGAMR(2)
XSCBU3 = XKERA1*ZGAMR(3)
XSCBU_EFF1 = 0.6E-3
XSCBU_EFF2 = 2.0E-3
XSCBUEXP1 = -2500.0
!
!
!-------------------------------------------------------------------------------
!
!*       7.     CONSTANTS FOR THE "SONTANEOUS" BREAK-UP
!               ---------------------------------------
!
XSPONBUD1 = 3.0E-3
XSPONBUD2 = 4.0E-3
XSPONBUD3 = 5.0E-3
XSPONCOEF2 = ((XSPONBUD3/XSPONBUD2)**3 - 1.0)/(XSPONBUD3-XSPONBUD1)**2
!
!
!------------------------------------------------------------------------------
!
!*        8.    CONSTANTS FOR EVAPORATION PROCESS
!               ---------------------------------
!
X0CNDC = (4.0*XPI)*XC1C*XF0C*MOMG(XALPHAC,XNUC,1.)
X2CNDC = (4.0*XPI)*XC1C*XF2C*XCC*MOMG(XALPHAC,XNUC,XDC+2.0)
!
XEX0EVAR = -1.0
XEX1EVAR = -1.0 - (XDR+1.0)*0.5
XEX2EVAR = -0.5*XCEXVT
!
X0EVAR = (2.0*XPI)*XF0R*GAMMA(XNUR+1./XALPHAR)/GAMMA(XNUR)
X1EVAR = (2.0*XPI)*XF1R*((ZRHO00)**(XCEXVT)*(XCR/0.15E-4))**0.5*    &
	   GAMMA(XNUR+(XDR+3.0)/(2.0*XALPHAR))/GAMMA(XNUR)
!
XEX0EVAR = 2.0
XEX1EVAR = 2.0 - (XDR+1.0)*0.5
XEX2EVAR = -0.5*XCEXVT
!
X0EVAR = (12.0)*XF0R*GAMMA(XNUR+1./XALPHAR)/GAMMA(XNUR+3./XALPHAR)
X1EVAR = (12.0)*XF1R*((ZRHO00)**(XCEXVT)*(XCR/0.15E-4))**0.5*    &
	   GAMMA(XNUR+(XDR+3.0)/(2.0*XALPHAR))/GAMMA(XNUR+3./XALPHAR)
!
!-------------------------------------------------------------------------------
!
!*       9.     SET-UP RADIATIVE PARAMETERS
!               ---------------------------
!
! R_eff_c = XFREFFC * (rho*r_c/N_c)**(1/3)
!
!
XFREFFC = 0.5 * ZGAMC(6) * (1.0/XAC)**(1.0/3.0)
XFREFFR = 0.5 * ZGAMR(6) * (1.0/XAR)**(1.0/3.0)
!
! Coefficients used to compute reff when both cloud and rain are present
!
XCREC = 1.0/ (ZGAMC(6) * XAC**(2.0/3.0))
XCRER = 1.0/ (ZGAMR(6) * XAR**(2.0/3.0))
!
!-------------------------------------------------------------------------------
!
!*       10.     SOME PRINTS FOR CONTROL
!                -----------------------
!
!
GFLAG = .TRUE.
IF (GFLAG) THEN
  WRITE(UNIT=ILUOUT0,FMT='(" Summary of the cloud particule characteristics")')
  WRITE(UNIT=ILUOUT0,FMT='("             CLOUD")')
  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAR,XBR
  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XCC,XDC
  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAC,XNUC
  WRITE(UNIT=ILUOUT0,FMT='("               RAIN")')
  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
                                                      XAR,XBR
  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
                                                      XCR,XDR
  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
                                                      XALPHAR,XNUR
  WRITE(UNIT=ILUOUT0,FMT='(" Description of the nucleation spectrum")')
  WRITE(UNIT=ILUOUT0,FMT='("         C=",E13.6,"  k=",E13.6)') XCHEN, XKHEN
  WRITE(UNIT=ILUOUT0,FMT='("      Beta=",E13.6," MU=",E13.6)') XBETAHEN, XMUHEN
  WRITE(UNIT=ILUOUT0,FMT='("      CCN max=",E13.6)') XCONC_CCN
END IF
!
!-------------------------------------------------------------------------------
!
!*       11.     Constants only for KHKO scheme
!               ---------------------------
!
!*       11.1    Cst for the coalescence processes
!
XR0 = 25.0E-6
!
!*       11.2    Cst for evaporation processes
!
XCEVAP = 0.86
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
!------------------------------------------------------------------------------
!
!
END SUBROUTINE INI_RAIN_C2R2
