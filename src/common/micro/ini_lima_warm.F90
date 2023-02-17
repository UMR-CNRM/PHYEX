!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #########################
       MODULE MODI_INI_LIMA_WARM
!      #########################
!
INTERFACE
      SUBROUTINE INI_LIMA_WARM (PTSTEP, PDZMIN)
!
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step 
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
END SUBROUTINE INI_LIMA_WARM
!
END INTERFACE
!
END MODULE MODI_INI_LIMA_WARM
!     #########################################
      SUBROUTINE INI_LIMA_WARM (PTSTEP, PDZMIN)
!     #########################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to initialize the constants used in the 
!!    microphysical scheme LIMA for the warm phase species and processes. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST
USE MODD_REF
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_WARM
USE MODD_PARAMETERS
!USE MODD_LUNIT, ONLY : TLUOUT0
!
USE MODE_LIMA_FUNCTIONS, ONLY: MOMG
USE MODI_HYPGEO
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step 
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB                ! Coordinates of the first and last physical 
                              ! points along z
INTEGER :: J1                 ! Internal loop indexes
INTEGER :: JMOD               ! Internal loop to index the CCN modes
!
REAL, DIMENSION(6)  :: ZGAMC, ZGAMR ! parameters involving various moments of
                              ! the generalized gamma law
!
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
!INTEGER  :: ILUOUT0 ! Logical unit number for output-listing
!INTEGER  :: IRESP   ! Return code of FM-routines
!LOGICAL  :: GFLAG   ! Logical flag for printing the constatnts on the output
                    ! listing
!  
!-------------------------------------------------------------------------------
!
!
!*       1.     CHARACTERISTICS OF THE SPECIES
!   	        ------------------------------
!
!
!*       1.1    Cloud droplet characteristics
!
XAC = (XPI/6.0)*XRHOLW
XBC = 3.0
XCC = XRHOLW*XG/(18.0*1.816E-5) ! Stokes flow (Pruppacher eq. 10-138 for T=293K)
!XCC = XRHOLW*XG/(18.0*1.7E-5) ! Stokes flow (Pruppacher eq. 10-138 for T=273K)
XDC = 2.0
!
XF0C = 1.00
XF2C = 0.108
!
XC1C = 1./2.
!
!*       1.2    Raindrops characteristics
!
XAR = (XPI/6.0)*XRHOLW
XBR = 3.0
XCR = 842.
XDR = 0.8
!
XF0R = 0.780
!Correction BVIE Pruppacher 1997 eq. 13-61
!XF1R = 0.265
XF1R = 0.308
!
!
!------------------------------------------------------------------------------
!
!
!*       2.     DIMENSIONAL DISTRIBUTIONS OF THE SPECIES
!   	        ----------------------------------------
!
!
!*       2.1    Cloud droplet distribution
!
!XALPHAC = 3.0  ! Gamma law of the Cloud droplet (here volume-like distribution)
!XNUC    = 3.0  ! Gamma law with little dispersion
!
!*       2.2    Raindrop distribution
!
!XALPHAR = 3.0  ! Gamma law of the raindrops (here volume-like distribution)
!XNUR    = 3.0  ! Gamma law for the raindrops 
!XNUR    = 0.1
!
!*       2.3    Precalculation of the gamma function momentum
!
!
ZGAMC(1) = GAMMA_X0D(XNUC)
ZGAMC(2) = MOMG(XALPHAC,XNUC,3.)
ZGAMC(3) = MOMG(XALPHAC,XNUC,6.)
ZGAMC(4) = ZGAMC(3)-ZGAMC(2)**2  ! useful for Sig_c
ZGAMC(5) = MOMG(XALPHAC,XNUC,9.)
ZGAMC(6) = MOMG(XALPHAC,XNUC,3.)**(2./3.)/MOMG(XALPHAC,XNUC,2.)
!
ZGAMR(1) = GAMMA_X0D(XNUR)
ZGAMR(2) = MOMG(XALPHAR,XNUR,3.)
ZGAMR(3) = MOMG(XALPHAR,XNUR,6.)
ZGAMR(4) = MOMG(XALPHAR,XNUR,6.)
ZGAMR(5) = MOMG(XALPHAR,XNUR,9.)
ZGAMR(6) = MOMG(XALPHAR,XNUR,3.)**(2./3.)/MOMG(XALPHAR,XNUR,2.)
!
!*       2.4    Csts for the shape parameter
!
XLBC   = XAR*ZGAMC(2)
XLBEXC = 1.0/XBC
!
XNR = 1.0/(XAR*MOMG(XALPHAR,XNUR,XBR))
XCCR   = 8.E6
XCXR   = -1.
IF (NMOM_R.EQ.1) THEN
   XLBEXR = 1.0/(XCXR-XBR)
   XLBR   = ( XAR*XCCR*MOMG(XALPHAR,XNUR,XBR) )**(-XLBEXR)
ELSE
   XLBR   = XAR*ZGAMR(2)
   XLBEXR = 1.0/XBR
END IF
!
!
!------------------------------------------------------------------------------
!
!
!*       3.     CONSTANTS FOR THE SEDIMENTATION
!   	        -------------------------------
!
!
!*       4.1    Exponent of the fall-speed air density correction
!
IKB = 1 + JPVEXT
! Correction
!ZRHO00 = XP00/(XRD*XTHVREFZ(IKB))
ZRHO00 = 1.2041 ! at P=1013.25hPa and T=20°C
!
!*       4.2    Constants for sedimentation
!
XFSEDRR  = XCR*GAMMA_X0D(XNUR+(XDR+3.)/XALPHAR)/GAMMA_X0D(XNUR+3./XALPHAR)*     &
            (ZRHO00)**XCEXVT
XFSEDCR  = XCR*GAMMA_X0D(XNUR+XDR/XALPHAR)/GAMMA_X0D(XNUR)*     &
            (ZRHO00)**XCEXVT
XFSEDRC  = XCC*GAMMA_X0D(XNUC+(XDC+3.)/XALPHAC)/GAMMA_X0D(XNUC+3./XALPHAC)*     &
            (ZRHO00)**XCEXVT
XFSEDCC  = XCC*GAMMA_X0D(XNUC+XDC/XALPHAC)/GAMMA_X0D(XNUC)*     &
            (ZRHO00)**XCEXVT

!
XLB(2)    = XLBC
XLBEX(2)  = XLBEXC
XD(2)     = XDC
XFSEDR(2) = XFSEDRC
XFSEDC(2) = XFSEDCC
!
XLB(3)    = XLBR
XLBEX(3)  = XLBEXR
XD(3)     = XDR
XFSEDR(3) = XFSEDRR
XFSEDC(3) = XFSEDCR
!
!------------------------------------------------------------------------------
!
!
!*       4.     CONSTANTS FOR THE NUCLEATION PROCESS
!   	        ------------------------------------
!
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
!
!
!	4.1   Tabulation of the hypergeometric functions in 'no units'
!             --------------------------------------------------------
!
!               In LIMA's nucleation parameterization,
!               supersaturation is not in % : Smax=0.01 for a 1% supersaturation.
!               This is accounted for in the modified Beta and C values.
!
!               Here, we tabulate the
!               F(mu,k/2, k/2+1 ,-Beta S**2) -> XHYPF12
!               F(mu,k/2,(k+3)/2,-Beta S**2) -> XHYPF32 functions
!               using a logarithmic scale for S
!
NHYP = 500 ! Number of points for the tabulation
ALLOCATE (XHYPF12( NHYP, NMOD_CCN ))
ALLOCATE (XHYPF32( NHYP, NMOD_CCN ))
!
ZSMIN = 1.0E-5  ! Minimum supersaturation set at 0.001 % 
ZSMAX = 5.0E-2  ! Maximum supersaturation set at 5 %
XHYPINTP1 = REAL(NHYP-1)/LOG(ZSMAX/ZSMIN)
XHYPINTP2 = REAL(NHYP)-XHYPINTP1*LOG(ZSMAX)
!
DO JMOD = 1,NMOD_CCN 
   DO J1 = 1,NHYP
      ZSS =ZSMAX*(ZSMIN/ZSMAX)**(REAL(NHYP-J1)/REAL(NHYP-1))
      XHYPF12(J1,JMOD) = HYPGEO(XMUHEN_MULTI(JMOD),0.5*XKHEN_MULTI(JMOD),&
                                0.5*XKHEN_MULTI(JMOD)+1.0,XBETAHEN_MULTI(JMOD),ZSS)
      XHYPF32(J1,JMOD) = HYPGEO(XMUHEN_MULTI(JMOD),0.5*XKHEN_MULTI(JMOD),&
                                0.5*XKHEN_MULTI(JMOD)+1.5,XBETAHEN_MULTI(JMOD),ZSS)
   END DO
ENDDO
!
NAHEN = 81 ! Tabulation for each Kelvin degree in the range XTT-40 to XTT+40
XAHENINTP1 = 1.0
XAHENINTP2 = 0.5*REAL(NAHEN-1) - XTT
!
!		Compute the tabulation of function of T :
!
!                                   1
!               XAHENG = -----------------------
!                          XCSTHEN * G**(3/2)
!
!            Compute constants for the calculation of Smax.
!            XCSTHEN = 1/(rho_l 2 pi)   
!            PSI1
!            PSI3
!            T
!            Lv
!            G
!
ALLOCATE (XAHENG(NAHEN))
ALLOCATE (XAHENG2(NAHEN))
ALLOCATE (XAHENG3(NAHEN))
ALLOCATE (XPSI1(NAHEN))
ALLOCATE (XPSI3(NAHEN))
XCSTHEN = 1.0 / ( XRHOLW*2.0*XPI )
DO J1 = 1,NAHEN
   ZTT = XTT + REAL(J1-(NAHEN-1)/2)                                          ! T
   ZLV = XLVTT+(XCPV-XCL)*(ZTT-XTT)                                          ! Lv
   XPSI1(J1) = (XG/(XRD*ZTT))*(XMV*ZLV/(XMD*XCPD*ZTT)-1.)                    ! Psi1
   XPSI3(J1) = -1*XMV*ZLV/(XMD*XRD*(ZTT**2))                                 ! Psi3
   ZG    = 1./( XRHOLW*( (XRV*ZTT)/                                        & ! G
                         (XDIVA*EXP(XALPW-(XBETAW/ZTT)-(XGAMW*ALOG(ZTT)))) &
                       + (ZLV/ZTT)**2/(XTHCO*XRV) ) )              
   XAHENG(J1) = XCSTHEN/(ZG)**(3./2.)
   XAHENG2(J1) = 1/(ZG)**(1./2.) * GAMMA_X0D(XNUC+1./XALPHAC)/GAMMA_X0D(XNUC)
   XAHENG3(J1) = (ZG) * GAMMA_X0D(XNUC+1./XALPHAC)/GAMMA_X0D(XNUC)
END DO
!-------------------------------------------------------------------------------
!
! Parameters used to initialise the droplet and drop concentration
! from the respective mixing ratios (used in RESTART_RAIN_C2R2)
!
! Droplet case
!
!!ALLOCATE(XCONCC_INI(SIZE(PNFS,1),SIZE(PNFS,2),SIZE(PNFS,3),SIZE(PNFS,4))) !NMOD_CCN))
!!  XCONCC_INI(:,:,:,:) = 0.8 * PNFS(:,:,:,:) ! 80% of the maximum CCN conc. is assumed
!
! Raindrop case
!
XCONCR_PARAM_INI = (1.E7)**3/(XPI*XRHOLW) ! MP law with N_O=1.E7 m-1 is assumed
!
!
!------------------------------------------------------------------------------
!
!
!*       5.     CONSTANTS FOR THE COALESCENCE PROCESSES
!   	        ---------------------------------------
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
XR0 = 25.0E-6
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
! ICE3 accretion of cloud droplets by rain drops
!
XFCACCR  = (XPI/4.0)*XCCR*XCR*(ZRHO00**XCEXVT)*MOMG(XALPHAR,XNUR,XDR+2.0)
XEXCACCR = -XDR-3.0
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
!------------------------------------------------------------------------------
!
!
!*       6.     CONSTANTS FOR THE "SONTANEOUS" BREAK-UP
!   	        ---------------------------------------
!
!
XSPONBUD1 = 3.0E-3
XSPONBUD2 = 4.0E-3
XSPONBUD3 = 5.0E-3
XSPONCOEF2 = ((XSPONBUD3/XSPONBUD2)**3 - 1.0)/(XSPONBUD3-XSPONBUD1)**2
!
!
!------------------------------------------------------------------------------
!
!
!*        7.    CONSTANTS FOR EVAPORATION PROCESS
!   	        ---------------------------------------
!
!
X0CNDC = (4.0*XPI)*XC1C*XF0C*MOMG(XALPHAC,XNUC,1.)
X2CNDC = (4.0*XPI)*XC1C*XF2C*XCC*MOMG(XALPHAC,XNUC,XDC+2.0)
!
! Valeurs utiles pour le calcul de l'évaporation en fonction de N_r
!
!XEX0EVAR = -1.0
!XEX1EVAR = -1.0 - (XDR+1.0)*0.5
!XEX2EVAR = -0.5*XCEXVT
!
!X0EVAR = (2.0*XPI)*XF0R*GAMMA_X0D(XNUR+1./XALPHAR)/GAMMA_X0D(XNUR)
!X1EVAR = (2.0*XPI)*XF1R*((ZRHO00)**(XCEXVT)*(XCR/0.15E-4))**0.5*    &
!           GAMMA_X0D(XNUR+(XDR+3.0)/(2.0*XALPHAR))/GAMMA_X0D(XNUR)
!
!
! Valeurs utiles pour le calcul de l'évaporation en fonction de r_r
!
XEX0EVAR = 2.0
XEX1EVAR = 2.0 - (XDR+1.0)*0.5
XEX2EVAR = -0.5*XCEXVT
!
X0EVAR = (12.0)*XF0R*GAMMA_X0D(XNUR+1./XALPHAR)/GAMMA_X0D(XNUR+3./XALPHAR)
X1EVAR = (12.0)*XF1R*((ZRHO00)**(XCEXVT)*(XCR/0.15E-4))**0.5*    &
           GAMMA_X0D(XNUR+(XDR+3.0)/(2.0*XALPHAR))/GAMMA_X0D(XNUR+3./XALPHAR)
!
XCEVAP = 0.86
!
!------------------------------------------------------------------------------
!
!
!*       8.     SET-UP RADIATIVE PARAMETERS
!   	        ---------------------------
!
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
!
!------------------------------------------------------------------------------
!
!
!*       9.     SOME PRINTS FOR CONTROL
!   	        -----------------------
!
!
!!$GFLAG = .TRUE.
!!$IF (GFLAG) THEN
!!$  ILUOUT0 = TLUOUT0%NLU
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Summary of the cloud particule characteristics")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("             CLOUD")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
!!$                                                      XAR,XBR
!!$  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
!!$                                                      XCC,XDC
!!$  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
!!$                                                      XALPHAC,XNUC
!!$  WRITE(UNIT=ILUOUT0,FMT='("               RAIN")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("                   masse: A=",E13.6," B=",E13.6)') &
!!$                                                      XAR,XBR
!!$  WRITE(UNIT=ILUOUT0,FMT='("                 vitesse: C=",E13.6," D=",E13.6)') &
!!$                                                      XCR,XDR
!!$  WRITE(UNIT=ILUOUT0,FMT='("            distribution:AL=",E13.6,"NU=",E13.6)') &
!!$                                                      XALPHAR,XNUR
!!$  WRITE(UNIT=ILUOUT0,FMT='(" Description of the nucleation spectrum")')
!!$  WRITE(UNIT=ILUOUT0,FMT='("         C=",E13.6,"  k=",E13.6)') XCHEN, XKHEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      Beta=",E13.6," MU=",E13.6)') XBETAHEN, XMUHEN
!!$  WRITE(UNIT=ILUOUT0,FMT='("      CCN max=",E13.6)') XCONC_CCN
!!$END IF
!
!------------------------------------------------------------------------------
!
END SUBROUTINE INI_LIMA_WARM
