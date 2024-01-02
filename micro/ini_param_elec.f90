!MNH_LIC Copyright 1994-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     	##########################
        MODULE MODI_INI_PARAM_ELEC
!      	##########################
!
IMPLICIT NONE
INTERFACE
!
        SUBROUTINE INI_PARAM_ELEC (TPINIFILE, HGETSVT, HCLOUD, HELEC, &
                                   PRHO00, KRR, IIU, IJU, IKU )
!
USE MODD_IO, ONLY : TFILEDATA
IMPLICIT NONE
!
TYPE(TFILEDATA),                INTENT(IN) :: TPINIFILE ! Initial file
CHARACTER (LEN=*), DIMENSION(:),INTENT(IN) :: HGETSVT
CHARACTER (LEN=4),              INTENT(IN) :: HCLOUD   ! microphysics scheme
CHARACTER (LEN=4),              INTENT(IN) :: HELEC    ! electrical scheme
INTEGER, INTENT(IN) :: KRR      ! Number of moist variables
REAL,    INTENT(IN) :: PRHO00   ! Pressure at ground level
INTEGER, INTENT(IN) :: IIU      ! Upper dimension in x direction (local)
INTEGER, INTENT(IN) :: IJU      ! Upper dimension in y direction (local)
INTEGER, INTENT(IN) :: IKU      ! Upper dimension in z direction
!
END SUBROUTINE INI_PARAM_ELEC
END INTERFACE
END MODULE MODI_INI_PARAM_ELEC
!
!       ###############################################################
        SUBROUTINE INI_PARAM_ELEC (TPINIFILE, HGETSVT, HCLOUD, HELEC, &
                                   PRHO00, KRR, IIU, IJU, IKU)
!       ###############################################################
!
!!****  *INI_PARAM_ELEC* -  initialize the constants necessary 
!!                          for the electrical scheme.
!!
!!    PURPOSE
!!    -------
!!        The purpose of this routine is to initialize the constants used to 
!!	resolve the electrical scheme.
!!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      None
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!	Helsdon and Farley, 1987: A numerical study of a Montana thunderstorm: 
!!	  2. Model results versus observations involving electrical aspects.
!!	  J. Geophys. Res., 92, 5661-5675.
!!
!!	Takahashi, 1978: Riming electrification as a charge generation 
!!	  mechanism in thunderstorms. J. Atmos. Sci., 35, 1536-1548.
!!
!!	Gardiner et al., 1985: Measurements of initial potential gradient and
!!	  particles charges in a Montana supercell thunderstorm.
!!	  J. Geophys. Res., 90, 6079-6086.
!!
!!	Saunders et al., 1991: The effect of liquid water on thunderstorm
!!	  charging. J. Geophys. Res., 96, 11007-11017.
!!
!!    AUTHOR
!!    ------
!!      Gilles Molinie      * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!        C. Barthe  01/02/2004 coefficients f/b
!!        C. Barthe  21/05/2004 Add limitations for the NI processes
!!        C. Barthe  31/05/2004 Add constants for the inductive process
!!        C. Barthe  10/11/2009 Update to Masdev 4.8.1
!!        M. Chong      26/01/10  Small ions parameters 
!!                               +Fair weather field from Helsdon-Farley
!!                                (JGR, 1987, 5661-5675)
!!        J.-P. Pinty jan 2015  tabulate the equations for Saunders
!!        J. Escobar 8/01/2016 bug , missing YDIR='XY' in READ 
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  J. Wurtz       03/2022: new snow characteristics
!!        C. Barthe  04/02/2022 Add XGAMINC_RIM3 (no more initialized in ini_rain_ice_elec)
!!        C. Barthe  07/06/2022 Add parameters for charge sedimentation in LIMA
!!        C. Barthe  30/11/2022 Remove the section about charge neutralization ; 
!!                              already done in ini_flash_geom_elec
!!        C. Barthe  28/03/2023 Add parameters for sedimentation of cloud droplets charge
!!        C. Barthe  13/07/2023 Modify parameters that contain C_x and x_x for Ns, Ng and Nh
!
!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
USE MODD_CST
USE MODD_ELEC_n
USE MODD_ELEC_DESCR
USE MODD_ELEC_PARAM
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_NSV,            ONLY: NSV_ELECEND
USE MODD_PARAMETERS
USE MODD_PARAM_ICE_n
USE MODD_PARAM_LIMA,      ONLY : XALPHAC_L=>XALPHAC, XNUC_L=>XNUC, XALPHAR_L=>XALPHAR, XNUR_L=>XNUR, &
                                 XALPHAI_L=>XALPHAI, XNUI_L=>XNUI, XALPHAS_L=>XALPHAS, XNUS_L=>XNUS, &
                                 XALPHAG_L=>XALPHAG, XNUG_L=>XNUG,                                   &
                                 XCEXVT_L=>XCEXVT
USE MODD_PARAM_LIMA_COLD, ONLY : XAI_L=>XAI, XBI_L=>XBI, XC_I_L=>XC_I, XDI_L=>XDI, &
                                 XAS_L=>XAS, XBS_L=>XBS, XCS_L=>XCS, XDS_L=>XDS, XCCS_L=>XCCS, XCXS_L=>XCXS
USE MODD_PARAM_LIMA_MIXED,ONLY : XAG_L=>XAG, XBG_L=>XBG, XCG_L=>XCG, XDG_L=>XDG, XCCG_L=>XCCG, XCXG_L=>XCXG,     &
                                 XAH_L=>XAH, XBH_L=>XBH, XCH_L=>XCH, XDH_L=>XDH, XCCH_L=>XCCH, XCXH_L=>XCXH,     &
                                 XALPHAH_L=>XALPHAH, XNUH_L=>XNUH,                                               &
                                 XGAMINC_BOUND_MIN_L=>XGAMINC_BOUND_MIN, XGAMINC_BOUND_MAX_L=>XGAMINC_BOUND_MAX, &
                                 NGAMINC_L=>NGAMINC, NACCLBDAR_L=>NACCLBDAR, NACCLBDAS_L=>NACCLBDAS,             &
                                 XACCLBDAR_MIN_L=>XACCLBDAR_MIN, XACCLBDAR_MAX_L=>XACCLBDAR_MAX,                 &
                                 XACCLBDAS_MIN_L=>XACCLBDAS_MIN, XACCLBDAS_MAX_L=>XACCLBDAS_MAX,                 &
                                 NDRYLBDAR_L=>NDRYLBDAR, NDRYLBDAS_L=>NDRYLBDAS, NDRYLBDAG_L=>NDRYLBDAG,         &
                                 XDRYLBDAR_MIN_L=>XDRYLBDAR_MIN, XDRYLBDAR_MAX_L=>XDRYLBDAR_MAX,                 &
                                 XDRYLBDAS_MIN_L=>XDRYLBDAS_MIN, XDRYLBDAS_MAX_L=>XDRYLBDAS_MAX,                 &
                                 XDRYLBDAG_MIN_L=>XDRYLBDAG_MIN, XDRYLBDAG_MAX_L=>XDRYLBDAG_MAX
USE MODD_PARAM_LIMA_WARM, ONLY : XAR_L=>XAR, XBR_L=>XBR, XCR_L=>XCR, XDR_L=>XDR, &
                                 XCC_L=>XCC, XDC_L=>XDC, XCCR_L=>XCCR
USE MODD_RAIN_ICE_DESCR_n,ONLY : XALPHAC_I=>XALPHAC, XNUC_I=>XNUC, XALPHAR_I=>XALPHAR, XNUR_I=>XNUR,         &
                                 XALPHAI_I=>XALPHAI, XNUI_I=>XNUI, XALPHAS_I=>XALPHAS, XNUS_I=>XNUS,         &
                                 XALPHAG_I=>XALPHAG, XNUG_I=>XNUG, XALPHAH_I=>XALPHAH, XNUH_I=>XNUH,         &
                                 XCC_I=>XCC, XDC_I=>XDC,                                                     &
                                 XAR_I=>XAR, XBR_I=>XBR, XCR_I=>XCR, XDR_I=>XDR, XCCR_I=>XCCR,               &
                                 XAI_I=>XAI, XBI_I=>XBI, XC_I_I=>XC_I, XDI_I=>XDI,                           &
                                 XAS_I=>XAS, XBS_I=>XBS, XCS_I=>XCS, XDS_I=>XDS, XCCS_I=>XCCS, XCXS_I=>XCXS, &
                                 XAG_I=>XAG, XBG_I=>XBG, XCG_I=>XCG, XDG_I=>XDG, XCCG_I=>XCCG, XCXG_I=>XCXG, &
                                 XAH_I=>XAH, XBH_I=>XBH, XCH_I=>XCH, XDH_I=>XDH, XCCH_I=>XCCH, XCXH_I=>XCXH, &
                                 XCEXVT_I=>XCEXVT
USE MODD_RAIN_ICE_PARAM_n,ONLY : XGAMINC_BOUND_MIN_I=>XGAMINC_BOUND_MIN, XGAMINC_BOUND_MAX_I=>XGAMINC_BOUND_MAX, &
                                 NGAMINC_I=>NGAMINC, NACCLBDAR_I=>NACCLBDAR, NACCLBDAS_I=>NACCLBDAS,             &
                                 XACCLBDAR_MIN_I=>XACCLBDAR_MIN, XACCLBDAR_MAX_I=>XACCLBDAR_MAX,                 &
                                 XACCLBDAS_MIN_I=>XACCLBDAS_MIN, XACCLBDAS_MAX_I=>XACCLBDAS_MAX,                 &
                                 NDRYLBDAR_I=>NDRYLBDAR, NDRYLBDAS_I=>NDRYLBDAS, NDRYLBDAG_I=>NDRYLBDAG,         &
                                 XDRYLBDAR_MIN_I=>XDRYLBDAR_MIN, XDRYLBDAR_MAX_I=>XDRYLBDAR_MAX,                 &
                                 XDRYLBDAS_MIN_I=>XDRYLBDAS_MIN, XDRYLBDAS_MAX_I=>XDRYLBDAS_MAX,                 &
                                 XDRYLBDAG_MIN_I=>XDRYLBDAG_MIN, XDRYLBDAG_MAX_I=>XDRYLBDAG_MAX
USE MODD_VAR_ll
!
USE MODE_IO_FIELD_READ,  only: IO_Field_read
!
USE MODI_GAMMA_INC
USE MODI_MOMG
USE MODE_RRCOLSS, ONLY: RRCOLSS
USE MODE_RSCOLRG, ONLY: RSCOLRG
USE MODE_RZCOLX,  ONLY: RZCOLX
USE MODI_VQZCOLX
!
IMPLICIT NONE
!
!*	0.1	 Declaration of dummy arguments
!
TYPE(TFILEDATA),                INTENT(IN) :: TPINIFILE ! Initial file
CHARACTER (LEN=*), DIMENSION(:),INTENT(IN) :: HGETSVT
CHARACTER (LEN=4),              INTENT(IN) :: HCLOUD   ! microphysics scheme
CHARACTER (LEN=4),              INTENT(IN) :: HELEC    ! electrical scheme
INTEGER, INTENT(IN) :: KRR      ! Number of moist variables
REAL,    INTENT(IN) :: PRHO00   ! Pressure at ground level
INTEGER, INTENT(IN) :: IIU      ! Upper dimension in x direction (local)
INTEGER, INTENT(IN) :: IJU      ! Upper dimension in y direction (local)
INTEGER, INTENT(IN) :: IKU      ! Upper dimension in z direction
!
!*	0.2	 Declaration of local variables
!
INTEGER :: IND                ! Number of intervals to integrate kernels
INTEGER :: J1, JLWC, JTEMP
INTEGER :: IGAMINC, IACCLBDAR, IACCLBDAS, IDRYLBDAR, IDRYLBDAS, IDRYLBDAG
!
REAL    :: ZRATE              ! Geometrical growth of Lbda in the tabulated
                              ! functions and kernels
REAL    :: ZBOUND             ! XDCSLIM*Lbda_s: upper bound for the partial
                              ! integration of the riming rate of the aggregates
!
REAL    :: ZESR               ! Mean efficiency of rain-aggregate collection
REAL    :: ZEGS               !  
REAL    :: ZEGR
REAL    :: ZFDINFTY           ! Factor used to define the "infinite" diameter
! variables used to cope with the module variables common to icex and lima
REAL    :: ZCEXVT,                                                     &
           ZCC, ZDC, ZALPHAC, ZNUC,                                    &
           ZAR, ZBR, ZCR, ZDR, ZCCR, ZALPHAR, ZNUR,                    &
           ZAI, ZBI, ZCI, ZDI, ZALPHAI, ZNUI,                          &
           ZAS, ZBS, ZCS, ZDS, ZCCS, ZCXS, ZALPHAS, ZNUS,              &
           ZAG, ZBG, ZCG, ZDG, ZCCG, ZCXG, ZALPHAG, ZNUG,              &
           ZAH, ZBH, ZCH, ZDH, ZCCH, ZCXH, ZALPHAH, ZNUH,              &
           ZGAMINC_BOUND_MIN, ZGAMINC_BOUND_MAX,                       &
           ZACCLBDAR_MIN, ZACCLBDAR_MAX, ZACCLBDAS_MIN, ZACCLBDAS_MAX, &
           ZDRYLBDAR_MIN, ZDRYLBDAR_MAX, ZDRYLBDAS_MIN, ZDRYLBDAS_MAX, &
           ZDRYLBDAG_MIN, ZDRYLBDAG_MAX
!
REAL, DIMENSION(:,:), ALLOCATABLE :: ZMANSELL1, ZMANSELL2 ! Used to initialize
                                                          ! XMANSELL array
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZT, ZLWCC, ZEW
!
!-------------------------------------------------------------------------------
!
!*	1.	PRELIMINARIES
!		-------------
!
!*      1.1     Constants for electricity
!
XEPSILON = 8.85E-12     ! Dielectric permittivity of the air
XECHARGE = 1.6E-19      ! Elementary charge (C) 
!
!
!*      1.2     Address module variables common to ICEx and LIMA
!
IF (HCLOUD(1:3) == 'ICE') THEN
  ZCEXVT = XCEXVT_I
  !
  ZCC = XCC_I
  ZDC = XDC_I
  ZALPHAC = XALPHAC_I
  ZNUC    = XNUC_I
  !
  ZAR     = XAR_I
  ZBR     = XBR_I
  ZCR     = XCR_I
  ZDR     = XDR_I
  ZCCR    = XCCR_I
  ZALPHAR = XALPHAR_I
  ZNUR    = XNUR_I
  !
  ZAI     = XAI_I
  ZBI     = XBI_I
  ZCI     = XC_I_I
  ZDI     = XDI_I
  ZALPHAI = XALPHAI_I
  ZNUI    = XNUI_I
  !
  ZAS     = XAS_I
  ZBS     = XBS_I
  ZCS     = XCS_I
  ZDS     = XDS_I
  ZCCS    = XCCS_I
  ZCXS    = XCXS_I
  ZALPHAS = XALPHAS_I
  ZNUS    = XNUS_I
  !
  ZAG     = XAG_I
  ZBG     = XBG_I
  ZCG     = XCG_I
  ZDG     = XDG_I
  ZCCG    = XCCG_I
  ZCXG    = XCXG_I
  ZALPHAG = XALPHAG_I
  ZNUG    = XNUG_I
  !
  ZAH     = XAH_I
  ZBH     = XBH_I
  ZCH     = XCH_I
  ZDH     = XDH_I
  ZCCH    = XCCH_I
  ZCXH    = XCXH_I
  ZALPHAH = XALPHAH_I
  ZNUH    = XNUH_I
  !
  IGAMINC = NGAMINC_I
  ZGAMINC_BOUND_MIN = XGAMINC_BOUND_MIN_I
  ZGAMINC_BOUND_MAX = XGAMINC_BOUND_MAX_I
  !
  IACCLBDAR = NACCLBDAR_I
  IACCLBDAS = NACCLBDAS_I
  ZACCLBDAR_MIN = XACCLBDAR_MIN_I
  ZACCLBDAR_MAX = XACCLBDAR_MAX_I
  ZACCLBDAS_MIN = XACCLBDAS_MIN_I
  ZACCLBDAS_MAX = XACCLBDAS_MAX_I
  !
  IDRYLBDAR = NDRYLBDAR_I
  IDRYLBDAS = NDRYLBDAS_I
  IDRYLBDAG = NDRYLBDAG_I
  ZDRYLBDAR_MIN = XDRYLBDAR_MIN_I
  ZDRYLBDAR_MAX = XDRYLBDAR_MAX_I
  ZDRYLBDAS_MIN = XDRYLBDAS_MIN_I
  ZDRYLBDAS_MAX = XDRYLBDAS_MAX_I
  ZDRYLBDAG_MIN = XDRYLBDAG_MIN_I
  ZDRYLBDAG_MAX = XDRYLBDAG_MAX_I
  !
ELSE IF (HCLOUD == 'LIMA') THEN
  ZCEXVT = XCEXVT_L
  !
  ZCC = XCC_L
  ZDC = XDC_L
  ZALPHAC = XALPHAC_L
  ZNUC    = XNUC_L
  !
  ZAR     = XAR_L
  ZBR     = XBR_L
  ZCR     = XCR_L
  ZDR     = XDR_L
  ZCCR    = XCCR_L
  ZALPHAR = XALPHAR_L
  ZNUR    = XNUR_L
  !
  ZAI     = XAI_L
  ZBI     = XBI_L
  ZCI     = XC_I_L
  ZDI     = XDI_L
  ZALPHAI = XALPHAI_L
  ZNUI    = XNUI_L
  !
  ZAS     = XAS_L
  ZBS     = XBS_L
  ZCS     = XCS_L
  ZDS     = XDS_L
  ZCCS    = XCCS_L
  ZCXS    = XCXS_L
  ZALPHAS = XALPHAS_L
  ZNUS    = XNUS_L
  !
  ZAG     = XAG_L
  ZBG     = XBG_L
  ZCG     = XCG_L
  ZDG     = XDG_L
  ZCCG    = XCCG_L
  ZCXG    = XCXG_L
  ZALPHAG = XALPHAG_L
  ZNUG    = XNUG_L
  !
  ZAH     = XAH_L
  ZBH     = XBH_L
  ZCH     = XCH_L
  ZDH     = XDH_L
  ZCCH    = XCCH_L
  ZCXH    = XCXH_L
  ZALPHAH = XALPHAH_L
  ZNUH    = XNUH_L
  !
  IGAMINC = NGAMINC_L
  ZGAMINC_BOUND_MIN = XGAMINC_BOUND_MIN_L
  ZGAMINC_BOUND_MAX = XGAMINC_BOUND_MAX_L
  !
  IACCLBDAR = NACCLBDAR_L
  IACCLBDAS = NACCLBDAS_L
  ZACCLBDAR_MIN = XACCLBDAR_MIN_L
  ZACCLBDAR_MAX = XACCLBDAR_MAX_L
  ZACCLBDAS_MIN = XACCLBDAS_MIN_L
  ZACCLBDAS_MAX = XACCLBDAS_MAX_L
  !
  IDRYLBDAR = NDRYLBDAR_L
  IDRYLBDAS = NDRYLBDAS_L
  IDRYLBDAG = NDRYLBDAG_L
  ZDRYLBDAR_MIN = XDRYLBDAR_MIN_L
  ZDRYLBDAR_MAX = XDRYLBDAR_MAX_L
  ZDRYLBDAS_MIN = XDRYLBDAS_MIN_L
  ZDRYLBDAS_MAX = XDRYLBDAS_MAX_L
  ZDRYLBDAG_MIN = XDRYLBDAG_MIN_L
  ZDRYLBDAG_MAX = XDRYLBDAG_MAX_L
END IF  
!
!-------------------------------------------------------------------------------
!
!*	2.	SHAPE PARAMETERS
!		----------------
!
XCXR = -1.0 ! Raindrop characteristic : XCXR (not initialized in ini_rain_ice.f90)
!
! Individual charge  q(d) = e_x * d ** f_x with f_x = XFx 
!
XFC = 0.5   ! cloud 
XFR = 1.3   ! rain 
XFI = 0.5   ! pristine ice 
XFS = 1.3   ! snow
XFG = 2.0   ! graupel
XFH = 2.0   ! hail
!
! Min/max values of e_x
!
XEGMIN = 1.E-12
XEGMAX = 1.E-3 
!
XESMIN = 1.E-14
XESMAX = 1.E-4 
!
XEIMIN = 1.E-12
XEIMAX = 1.E-3 
!
XECMIN = 1.E-12
XECMAX = 1.E-3
! 
XERMIN = 1.E-14
XERMAX = 1.E-4
!
XEHMIN = 1.E-14
XEHMAX = 1.E-4 
!
!  E=E_0 * exp(k_e*z)
XE_0 = -100.
XKEF = -2.E-4            !  229.E-6
!
!  E=E_0 (b1 exp(-a1 z) + b2 exp(-a2 z) + b3 exp(-a3 z) : Helsdon-Farley, 1987
XE0_HF = -80.
XA1_HF = 4.5E-3
XB1_HF = 0.5
XA2_HF = 3.8E-4
XB2_HF = 0.65
XA3_HF = 1.0E-4
XB3_HF = 0.1
!
XIONCOMB = 1.6E-12
XF_POS = 1.4E-4
XF_NEG = 1.9E-4
XEXPMOB = 1.4E-4
!
XFCORONA = 2.E-20
XECORONA = 5000.
!
XJCURR_FW = -2.7E-12
!
!
!-------------------------------------------------------------------------------
!
!*	3.	COEFFICIENTS FOR CHARGE TRANSFERS
!		---------------------------------
!
! proportionality coefficient between mass transfer and charge transfer rates
! the mixing ratio is proportional to the volume of the particle
! the electric charge is proportional to the surface of the particle
!
XCOEF_RQ_V = 1
XCOEF_RQ_C = XFC / 3.0  ! XBC=3
XCOEF_RQ_R = XFR / ZBR
XCOEF_RQ_I = XFI / ZBI
XCOEF_RQ_S = XFS / ZBS
XCOEF_RQ_G = XFG / ZBG
XCOEF_RQ_H = XFH / ZBH
!
!
!-------------------------------------------------------------------------------
!
!*	3.	HOMOGENEOUS NUCLEATION
!		----------------------
!
XALPHACQ = 3.         !>
XNUCQ    = 1.         ! >--- generic values
XLBDACQ  = 1.1E5      !>
!
XQHON = (1. / XRHOLW)
XQHON = XQHON * MOMG(XALPHACQ,XNUCQ,XFC+3.)
XQHON = XQHON / MOMG(XALPHACQ,XNUCQ,3.)
XQHON = XQHON / (XLBDACQ**XFC)
!
!
!-------------------------------------------------------------------------------
!
!*	4.	SEDIMENTATION
!		-------------
!
IF (ALLOCATED(XQTMIN)) DEALLOCATE(XQTMIN)       
IF (ALLOCATED(XRTMIN_ELEC)) DEALLOCATE(XRTMIN_ELEC)
!
IF (KRR == 6) THEN
  ALLOCATE( XQTMIN(6) )
  ALLOCATE( XRTMIN_ELEC(6) )
ELSE IF (KRR == 7) THEN
  ALLOCATE( XQTMIN(7) )
  ALLOCATE( XRTMIN_ELEC(7) )
END IF
!
XQTMIN(1) = 1.0E-17  !
XQTMIN(2) = 1.0E-17  !
XQTMIN(3) = 1.0E-17  ! ten particles per cubic meter that carried a charge
XQTMIN(4) = 1.0E-17  ! of one electron
XQTMIN(5) = 1.0E-17  !
XQTMIN(6) = 1.0E-17  !
IF (KRR == 7) XQTMIN(7) = 1.0E-17  !
!
XRTMIN_ELEC(1) = 1.0E-6  
XRTMIN_ELEC(2) = 1.0E-6 
XRTMIN_ELEC(3) = 1.0E-6
XRTMIN_ELEC(4) = 1.0E-6
XRTMIN_ELEC(5) = 1.0E-6
XRTMIN_ELEC(6) = 1.0E-6
IF (KRR == 7) XRTMIN_ELEC(7) = 1.0E-6
!
XLBDAR_MAXE = 2.E3 ! Less than 10000 particles in cube meter of cloud. 
XLBDAS_MAXE = 2.E3 ! Less than 10000 particles in cube meter of cloud. 
XLBDAG_MAXE = 2.E3 !
XLBDAH_MAXE = 2.E3 !
!
IF (HCLOUD(1:3) == 'ICE') THEN
  !
  ! Cloud droplets
  !
  ZCEXVT = 0.4
  XEXQSEDC = XFC + ZDC
  XFQSEDC = ZCC * MOMG(ZALPHAC,ZNUC,ZDC+XFC) * (PRHO00)**ZCEXVT
  !
  ! Rain
  !
  XEXQSEDR = (XCXR - XFR - ZDR) / (XCXR - ZBR)
  XFQSEDR  = ZCR * (ZCCR**(1 - XEXQSEDR)) * MOMG(ZALPHAR,ZNUR,ZDR+XFR) * &
           ((ZAR * MOMG(ZALPHAR,ZNUR,ZBR))**(-XEXQSEDR)) * (PRHO00)**ZCEXVT
  !
  ! Ice
!!++cb++ 23/02/23 pour la microphysique, calcul fait pour des colonnes
! => on fait pareil ici pour garder la coherence
!XEXQSEDI = (ZDI + XFI) / ZBI
!XFQSEDI  = ZCI * MOMG(ZALPHAI,ZNUI,ZDI+XFI) * (PRHO00**ZCEXVT) * &
!           (ZAI * MOMG(ZALPHAI,ZNUI,ZBI))**(-XEXQSEDI)  
  XEXQSEDI = (1.585 + XFI) / 1.7
  XFQSEDI  = 2.1E5 * MOMG(ZALPHAI,ZNUI,1.585+XFI) * (PRHO00**ZCEXVT) * &
             (2.14E-3 * MOMG(ZALPHAI,ZNUI,1.7))**(-XEXQSEDI)
!--cb--
  XFCI     = (4. * XPI * 900.)**(-1)
  !
  ! Snow
  !
  XEXQSEDS = (ZCXS - XFS - ZDS) / (ZCXS - ZBS)
  XFQSEDS  = ZCS * (ZCCS**(1 - XEXQSEDS)) * MOMG(ZALPHAS,ZNUS,ZDS+XFS) * &
           ((ZAS * MOMG(ZALPHAS,ZNUS,ZBS))**(-XEXQSEDS)) * (PRHO00)**ZCEXVT
  !
  ! Graupeln 
  !
  XEXQSEDG = (ZCXG - XFG - ZDG) / (ZCXG - ZBG)
  XFQSEDG  = ZCG * (ZCCG**(1 - XEXQSEDG)) * MOMG(ZALPHAG,ZNUG,ZDG+XFG) * &
           ((ZAG * MOMG(ZALPHAG,ZNUG,ZBG))**(-XEXQSEDG)) * (PRHO00)**ZCEXVT
  !
  ! Hail
  !
  XEXQSEDH = (ZCXH - XFH - ZDH) / (ZCXH - ZBH)
  XFQSEDH  = ZCH * (ZCCH**(1 - XEXQSEDH)) * MOMG(ZALPHAH,ZNUH,ZDH+XFH) * &
           ((ZAH * MOMG(ZALPHAH,ZNUH,ZBH))**(-XEXQSEDH)) * (PRHO00)**ZCEXVT
!
ELSE IF (HCLOUD == 'LIMA') THEN
  ALLOCATE(XFQSED(KRR))
  XFQSED(:) = 0.
  XFQSED(2) = ZCC * MOMG(ZALPHAC,ZNUC,ZDC+XFC)
  XFQSED(3) = ZCR * MOMG(ZALPHAR,ZNUR,ZDR+XFR)
  XFQSED(4) = ZCI * MOMG(ZALPHAI,ZNUI,ZDI+XFI)
  XFQSED(5) = ZCS * MOMG(ZALPHAS,ZNUS,ZDS+XFS)
  XFQSED(6) = ZCG * MOMG(ZALPHAG,ZNUG,ZDG+XFG)
  IF (KRR == 7) XFQSED(7) = ZCH * MOMG(ZALPHAH,ZNUH,ZDH+XFH)
  !
  ALLOCATE(XDQ(KRR))
  XDQ(:) = 0.
  XDQ(2) = ZDC + XFC
  XDQ(3) = ZDR + XFR
  XDQ(4) = ZDI + XFI
  XDQ(5) = ZDS + XFS
  XDQ(6) = ZDG + XFG
  IF (KRR == 7) XDQ(7) = ZDH + XFH
END IF
!
!-------------------------------------------------------------------------------
!
!*	5.	EVAPORATION OF RAINDROPS
!		------------------------
!
!XQREVAV1 = (2. / XPI) * MOMG(ZALPHAR,ZNUR,XFR) / MOMG(ZALPHAR,ZNUR,2.)
!XQREVAV2 = (XPI / ZAR) * (MOMG(ZALPHAR,ZNUR,2.) / MOMG(ZALPHAR,ZNUR,ZBR)) * &
!           (XCXR - 2.) / (XCXR - ZBR)
!
!
!-------------------------------------------------------------------------------
!
!*	6.	RIMING OF CLOUD DROPLETS ON SNOW
!		--------------------------------
!
IF (HELEC == 'ELE4') THEN
  XEXQSRIMCG = -XFS
  XQSRIMCG   = MOMG(ZALPHAS,ZNUS,XFS)
ELSE
  XEXQSRIMCG = ZCXS - XFS
  XQSRIMCG   = ZCCS * MOMG(ZALPHAS,ZNUS,XFS)
END IF
!
! The array containing the tabulated function M(fs,D_cs^lim)/M(fs) 
! is no more implemented in ini_rain_ice.f90
ZRATE = EXP(LOG(ZGAMINC_BOUND_MAX/ZGAMINC_BOUND_MIN)/REAL(IGAMINC-1))
IF( .NOT.ALLOCATED(XGAMINC_RIM3) ) ALLOCATE( XGAMINC_RIM3(IGAMINC) )
DO J1 = 1, IGAMINC
  ZBOUND = ZGAMINC_BOUND_MIN * ZRATE**(J1-1)
  XGAMINC_RIM3(J1) = GAMMA_INC(ZNUS+XFS/ZALPHAS,ZBOUND)
END DO
!
!
!-------------------------------------------------------------------------------
!
!*	7.	CONTACT FREEZING BETWEEN RAINDROPS AND PRISTINE ICE
!		---------------------------------------------------
!
IF (HELEC == 'ELE4') THEN
  XEXQRCFRIG = - ZDR - XFR - 2.0
  XQRCFRIG   = (XPI / 4.0) * ZCR * MOMG(ZALPHAR,ZNUR,ZDR+XFR+2.) * &
                PRHO00**ZCEXVT
ELSE
  XEXQRCFRIG = XCXR - ZDR - XFR - 2.0
  XQRCFRIG   = (XPI / 4.0) * ZCR * ZCCR * MOMG(ZALPHAR,ZNUR,ZDR+XFR+2.) * &
                PRHO00**ZCEXVT
END IF
!
!
!-------------------------------------------------------------------------------
!
!*	8.	INITIALIZATIONS FOR THE NON INDUCTIVE PROCESSES
!		-----------------------------------------------
!
! arrays allocation for NI charging rate
!
ALLOCATE( XNI_SDRYG(IIU, IJU, IKU) )
ALLOCATE( XNI_IDRYG(IIU, IJU, IKU) )
ALLOCATE( XNI_IAGGS(IIU, IJU, IKU) )
ALLOCATE( XIND_RATE(IIU, IJU, IKU) )
ALLOCATE( XEW(IIU, IJU, IKU) )
XEW(:,:,:) = 0.
!
SELECT CASE(HGETSVT(NSV_ELECEND))
  CASE ('READ')
    CALL IO_Field_read(TPINIFILE,'NI_IAGGS',XNI_IAGGS)
    CALL IO_Field_read(TPINIFILE,'NI_IDRYG',XNI_IDRYG)
    CALL IO_Field_read(TPINIFILE,'NI_SDRYG',XNI_SDRYG)
    CALL IO_Field_read(TPINIFILE,'INDUC_CG',XIND_RATE)
  CASE ('INIT')
    XNI_IAGGS(:,:,:) = 0.
    XNI_IDRYG(:,:,:) = 0.
    XNI_SDRYG(:,:,:) = 0.
    XIND_RATE(:,:,:) = 0.
END SELECT

!
!*      8.1     Gardiner et al. (1985) parameterization
!
IF (CNI_CHARGING == 'GARDI') THEN
  XLWCC = 0.1   ! g.m^-3
END IF
!
!
!*      8.2     Saunders et al. (1991) and 
!*              Saunders and Peck (1998) parameterizations
!
IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR.  &
    CNI_CHARGING == 'SAP98' .OR.                               &
    CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR.  &
    CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
!
! ice particle = the smallest particle (I-S and I-G collisions)
  XIMP = 3.76     ! for positive charge
  XINP = 2.5
  XIKP = 4.92E13
  XIKP_TAK = 6.1E12     ! for Takahashi
  XIMN = 2.54     ! for negative charge   
  XINN = 2.8
  XIKN = 5.25E8
  XIKN_TAK = 4.3E7      ! for Takahashi
! 
! snow = the smallest particle (S-G collisions)
  XSMP = 0.44     ! for positive charge
  XSNP = 2.5
  XSKP = 52.8
  XSKP_TAK = 6.5     ! for Takahashi
  XSMN = 0.5      ! for negative charge
  XSNN = 2.8
  XSKN = 24.
  XSKN_TAK = 2.0        ! for Takahashi
!
  XFQIAGGSP = XIKP * ZCS**(1. + XINP) *                 &
                MOMG(ZALPHAS, ZNUS, 2.+ZDS*(1.+XINP)) * &
                MOMG(ZALPHAI, ZNUI, XIMP)
  XFQIAGGSN = XIKN * ZCS**(1. + XINN) *                 &
                MOMG(ZALPHAS, ZNUS, 2.+ZDS*(1.+XINN)) * &
                MOMG(ZALPHAI, ZNUI, XIMN)
!
  XFQIDRYGBSP = XIKP * ZCG**(1. + XINP) *               &
                MOMG(ZALPHAG, ZNUG, 2.+ZDG*(1.+XINP)) * &
                MOMG(ZALPHAI, ZNUI, XIMP)
  XFQIDRYGBSN = XIKN * ZCG**(1. + XINN) *               &
                MOMG(ZALPHAG, ZNUG, 2.+ZDG*(1.+XINN)) * &
                MOMG(ZALPHAI, ZNUI, XIMN)
!
  XFQIAGGSP_TAK = XFQIAGGSP * XIKP_TAK / XIKP
  XFQIAGGSN_TAK = XFQIAGGSN * XIKN_TAK / XIKN
  XFQIDRYGBSP_TAK = XFQIDRYGBSP * XIKP_TAK / XIKP
  XFQIDRYGBSN_TAK = XFQIDRYGBSN * XIKN_TAK / XIKN
!
  XAIGAMMABI      = ZAI * MOMG(ZALPHAI, ZNUI, ZBI)
!
  XLBQSDRYGB1SP = MOMG(ZALPHAG,ZNUG,2.) * MOMG(ZALPHAS, ZNUS, XSMP)
  XLBQSDRYGB1SN = MOMG(ZALPHAG,ZNUG,2.) * MOMG(ZALPHAS, ZNUS, XSMN)
  XLBQSDRYGB2SP = 2. * MOMG(ZALPHAG,ZNUG,1.) * MOMG(ZALPHAS, ZNUS, 1.+XSMP)
  XLBQSDRYGB2SN = 2. * MOMG(ZALPHAG,ZNUG,1.) * MOMG(ZALPHAS, ZNUS, 1.+XSMN)
  XLBQSDRYGB3SP =                              MOMG(ZALPHAS, ZNUS, 2.+XSMP)
  XLBQSDRYGB3SN =                              MOMG(ZALPHAS, ZNUS, 2.+XSMN)
ENDIF
!
IF (CNI_CHARGING == 'SAP98' .OR. CNI_CHARGING == 'TERAR' .OR. &
    CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2') THEN
  XVSCOEF = ZCS * MOMG(ZALPHAS, ZNUS, ZBS+ZDS) / MOMG(ZALPHAS, ZNUS, ZBS)
  XVGCOEF = ZCG * MOMG(ZALPHAG, ZNUG, ZBG+ZDG) / MOMG(ZALPHAG, ZNUG, ZBG)
END IF
!
!
!*      8.3    Takahashi (1978) parameterization
!
IF (CNI_CHARGING == 'TAKAH') THEN
!
! last column and line are duplicated for interpolation
  NIND_TEMP = 31
  NIND_LWC  = 28
  IF( .NOT.ALLOCATED(XMANSELL)) ALLOCATE( XMANSELL(NIND_LWC+1,NIND_TEMP+1))
  ALLOCATE( ZMANSELL1(29,16) )
  ALLOCATE( ZMANSELL2(29,16) )
  ZMANSELL1 = RESHAPE( &
(/  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  1.65,  3.3 ,  &
    4.95, 13.2 , 19.8 , 26.4 , 31.02, 33.0 , 42.9 , 46.2 , 49.5 ,  &
   52.8 , 49.5 , 39.6 , 36.3 , 32.34, 32.01, 32.01, 31.68, 31.35, 31.35, 31.35,&
    0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.33,  2.64,  4.95,  6.6 ,  &
    6.93, 15.18, 23.1 , 29.7 , 33.0 , 36.3 , 49.5 , 52.8 , 52.8 ,  &
   52.8 , 49.5 , 36.3 , 33.0 , 32.34, 32.01, 32.01, 31.35, 31.35, 31.35, 31.35,&
    0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.33,  3.3 ,  4.95,  6.6 ,  7.59,  &
    9.24, 17.16, 26.4 , 32.34, 36.3 , 39.6 , 52.8 , 56.1 , 56.1 ,  &
   56.1 , 42.9 , 33.0 , 32.01, 32.01, 31.68, 31.68, 31.02, 31.02, 30.36, 30.36,&
    0.0 ,  0.0 ,  0.0 ,  0.33,  2.97,  5.61,  7.26,  8.25,  9.9 ,  &
   10.56, 19.8 , 28.71, 36.3 , 39.6 , 46.2 , 56.1 , 59.4 , 59.4 ,  &
   59.4 , 42.9 , 31.35, 31.35, 31.68, 31.35, 30.69, 30.36, 31.02, 30.03, 30.03,&
    0.0 ,  0.0 ,  0.33,  2.64,  5.61,  7.59,  9.24,  9.9 , 11.22,  &
   11.88, 22.44, 31.35, 39.6 , 47.2 , 52.8 , 59.4 , 62.7 , 62.7 ,  &
   52.8 , 33.0 , 30.69, 30.69, 30.36, 30.69, 30.36, 30.03, 30.36, 30.03, 30.03,&
    0.0 ,  0.0 ,  0.33,  4.29,  7.26,  8.91, 10.23, 11.55, 12.21,  &
   13.2 , 24.09, 33.0 , 39.6 , 49.5 , 56.1 , 56.1 , 56.1 , 66.0 ,  &
   49.5 , 31.35, 30.03, 29.7 , 29.37, 29.37, 29.37, 29.04, 30.03, 29.04, 29.04,&
    0.0 ,  0.0 ,  2.31,  6.6 ,  8.91,  9.9 , 11.22, 13.2 , 13.2 ,  &
   14.85, 25.08, 36.3 , 42.9 , 49.5 , 52.8 , 46.2 , 42.9 , 52.8 ,  &
   42.9 , 28.05, 28.05, 29.71, 28.38, 28.71, 28.71, 28.71, 29.37, 28.71, 28.71,&
    0.0 ,  0.0 ,  4.29,  7.59,  9.9 , 10.23, 11.22, 14.85, 14.85,  &
   16.5 , 26.4 , 36.3 , 42.9 , 42.9 , 46.2 , 42.9 , 39.6 , 46.2 ,  &
   33.0 , 24.75, 26.4 , 27.39, 27.39, 27.72, 28.05, 28.38, 29.04, 28.05, 28.05,&
    0.0 ,  0.0 ,  6.27,  9.24, 10.56, 12.21, 11.88, 15.18, 15.84,  &
   16.5 , 27.39, 33.0 , 39.6 , 39.6 , 39.6 , 36.3 , 33.0 , 39.6 ,  &
   19.8 , 16.5 , 23.1 , 25.74, 26.73, 27.06, 27.39, 27.39, 28.38, 27.39, 27.39,&
    0.0 ,  0.66,  6.93,  9.57, 11.22, 12.87, 12.54, 15.51, 16.5 ,  &
   18.15, 27.39, 31.02, 36.3 , 33.0 , 29.7 , 23.1 , 19.8 , 26.4 ,  &
    9.9 ,  0.99, 20.46, 23.43, 25.41, 26.4 , 27.06, 27.06, 27.39, 27.06, 27.06,&
    0.0 ,  2.64,  7.59, 10.23, 11.88, 13.53, 13.2 , 15.84, 16.5 ,  &
   19.14, 27.06, 28.71, 33.0 , 26.4 , 21.45, 16.5 ,  9.9 ,  3.3 ,  &
    0.0 , -6.6 , 18.15, 21.45, 23.76, 25.08, 26.07, 26.73, 26.73, 26.73, 26.73,&
    0.0 ,  3.3 ,  8.25, 10.89, 11.88, 13.53, 13.86, 15.84, 16.83,  &
   18.81, 26.73, 27.39, 28.71, 21.45, 16.5 ,  3.3 ,  0.99, -0.99,  &
   -4.95,-16.5 , 16.5 , 19.8 , 22.77, 24.09, 25.08, 25.41, 25.74, 26.4 , 26.4 ,&
    0.0 ,  4.62,  8.91, 11.22, 12.21, 13.53, 14.85, 16.17, 16.83,  &
   18.81, 25.74, 26.4 , 26.4 , 16.5 ,  2.97,  0.0 , -1.98, -9.9 ,  &
  -11.55,-19.8 ,  3.3 , 18.15, 20.13, 23.1 , 23.76, 24.42, 25.41, 25.41, 25.41,&
    0.0 ,  5.8 ,  9.24, 11.22, 12.54, 13.86, 15.18, 15.84, 16.83,  &
   18.48, 24.42, 23.43, 22.44,  3.3 ,  0.0 , -1.98, -4.29, -9.9 ,  &
  -13.2 ,-26.4 ,  0.0 , 16.83, 19.14, 21.45, 23.1 , 23.76, 24.42, 24.75, 24.75,&
    0.0 ,  5.94,  9.57, 11.22, 12.54, 13.86, 14.85, 15.84, 16.5 ,  &
   17.82, 23.1 , 20.79, 19.47,  0.33, -1.32, -3.3 , -6.6 , -9.9 ,  &
  -14.85,-33.0 , -0.33, 14.85, 18.48, 19.8 , 21.78, 23.1 , 23.76, 24.42, 24.42,&
    0.0 ,  5.94,  9.57, 11.22, 12.54, 13.53, 14.85, 15.51, 16.5 ,  &
   16.83, 21.45, 17.82,  9.9 , -0.33, -2.31, -3.3 , -8.25, -9.99,  &
  -14.85,-33.0 , -1.32,  4.95, 17.49, 19.47, 20.79, 22.44, 23.43, 24.09, 24.09/)&
    ,(/29, 16/))
  ZMANSELL2 = RESHAPE( &
(/  0.0 ,  6.27,  9.57, 10.89, 12.21, 13.53, 14.52, 15.18, 15.84,  &
   16.17, 19.47,  9.9 ,  0.0 , -0.99, -2.97, -4.95, -9.99, -9.99,  &
  -14.85,-29.7 , -2.97,  3.3 , 16.83, 19.14, 19.47, 21.12, 22.44, 23.43, 23.43,&
    0.0 ,  5.94,  9.57, 10.89, 12.21, 13.2 , 14.19, 14.85, 15.18,  &
   15.51, 17.16,  3.3 , -0.33, -1.65, -3.3 , -6.6 , -9.99, -9.99,  &
  -13.2 ,-28.1 , -3.3 ,  2.64, 16.5 , 17.82, 19.14, 20.13, 21.45, 22.77, 22.77,&
    0.0 ,  5.61,  8.91, 10.89, 11.88, 13.2 , 13.86, 14.52, 14.85,  &
   14.85,  9.9 ,  1.65, -0.99, -1.98, -3.3 , -6.6 , -9.99, -9.99,  &
  -13.2 ,-26.4 , -3.3 ,  1.65, 13.2 , 17.49, 18.81, 19.8 , 20.79, 22.11, 22.11,&
    0.0 ,  5.28,  8.58, 10.56, 11.55, 12.54, 13.2 , 13.86, 14.19,  &
   13.86,  6.6 ,  0.0 , -1.32, -2.31, -3.3 , -6.6 , -9.99, -9.99,  &
  -13.2 ,-24.8 , -3.3 ,  1.32,  6.6 , 17.16, 18.48, 19.47, 20.46, 21.45, 21.45,&
    0.0 ,  4.95,  8.25,  9.9 , 11.22, 11.88, 12.54, 12.87, 13.2 ,  &
   13.2 ,  3.3 , -0.66, -1.98, -2.64, -3.3 , -6.6 , -9.9 , -9.9 ,  &
  -13.2 ,-23.1 , -3.3 ,  0.66,  4.95, 16.5 , 17.82, 18.81, 19.8 , 20.46, 20.46,&
    0.0 ,  4.29,  7.59,  9.57, 10.56, 11.55, 11.88, 11.88, 12.21,  &
   11.88,  2.64, -0.66, -2.31, -2.64, -3.3 , -6.6 , -9.9 , -9.9 ,  &
  -13.2 ,-21.5 , -3.3 ,  0.0 ,  4.29, 14.85, 17.16, 18.48, 19.47, 20.13, 20.13,&
    0.0 ,  3.96,  6.93,  8.91,  9.9 , 10.89, 10.89, 10.89, 10.89,  &
   10.56,  0.99, -0.66, -2.31, -2.64, -3.3 , -6.6 , -9.9 , -9.99,  &
  -13.2 ,-19.8 , -6.6 ,  0.0 ,  3.3 , 13.2 , 16.83, 18.15, 19.47, 19.8 , 19.8 ,&
    0.0 ,  2.97,  6.27,  8.25,  9.24,  9.9 , 10.23, 10.23,  9.9 ,  &
    9.57,  0.0 , -0.99, -2.31, -2.64, -3.3 , -6.6 , -9.9 , -9.99,  &
  -11.55,-18.2 , -9.9 ,  0.0 ,  3.3 , 11.55, 16.5 , 17.82, 18.81, 19.8 , 19.8 ,&
    0.0 ,  0.66,  5.61,  7.59,  8.91,  9.24,  9.24,  9.24,  8.91,  &
    8.58,  0.0 , -0.99, -1.98, -2.64, -3.3 , -6.6 , -9.9 , -9.9 ,  &
  -11.55,-17.5 , -6.6 ,  0.0 ,  2.97,  8.25, 16.5 , 17.16, 18.48, 19.47, 19.47,&
    0.0 ,  0.0 ,  4.29,  6.6 ,  7.59,  8.25,  8.25,  7.59,  7.92,  &
    7.59, -0.33, -1.32, -1.98, -2.64, -3.3 , -4.95, -9.9 , -9.9 ,  &
   -9.9 ,-16.5 , -6.6 ,  0.0 ,  2.97,  6.6 , 14.85, 16.83, 18.15, 19.47, 19.47,&
    0.0 ,  0.0 ,  2.64,  5.28,  6.93,  7.26,  7.59,  7.26,  6.93,  &
    6.6 , -0.66, -1.32, -1.98, -2.64, -3.3 , -4.29, -8.25, -9.9 ,  &
   -9.9 ,-16.5 , -6.6 ,  0.0 ,  2.64,  6.6 , 14.85, 16.5 , 17.82, 19.14, 19.14,&
    0.0 ,  0.0 ,  0.33,  3.63,  5.61,  6.6 ,  6.6 ,  6.6 ,  4.95,  &
    3.63, -0.66, -1.32, -1.98, -2.64, -3.3 , -3.3 , -6.6 , -8.25,  &
   -8.25,-16.5 , -6.6 ,  0.0 ,  2.64,  6.6 , 13.2 , 16.5 , 17.49, 18.81, 18.81,&
    0.0 ,  0.0 ,  0.0 ,  0.99,  3.3 ,  4.29,  4.29,  4.62,  3.3 ,  &
    2.97, -0.66, -1.32, -1.98, -2.64, -2.97, -2.97, -3.3 , -4.95,  &
   -4.95,-15.8 , -4.95,  0.0 ,  2.31,  6.6 , 11.55, 16.5 , 17.49, 18.15, 18.15,&
    0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.99,  3.3 ,  2.64,  2.31,  2.31,  &
    2.31, -0.66, -1.32, -1.98, -2.31, -2.97, -2.64, -2.97, -3.63,  &
   -4.95, -9.99, -4.95,  0.0 ,  2.31,  6.6 ,  9.9 , 14.85, 17.16, 18.15, 18.15,&
    0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.66,  0.99,  1.65,  1.98,  &
    1.65, -0.66, -1.32, -1.65, -2.31, -2.64, -2.64, -2.97, -2.97,  &
   -3.3 , -9.99, -3.3 ,  0.0 ,  2.31,  6.6 ,  9.9 , 14.85, 17.16, 18.15, 18.15,&
    0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.0 ,  0.66,  0.99,  1.65,  1.98,  &
    1.65, -0.66, -1.32, -1.65, -2.31, -2.64, -2.64, -2.97, -2.97,  &
   -3.3 , -9.99, -3.3 ,  0.0 ,  2.31,  6.6 ,  9.9 , 14.85, 17.16, 18.15, 18.15/)&        
    ,(/29, 16/)) 
  XMANSELL(:, 1:16) = ZMANSELL1(:,:)
  XMANSELL(:,17:32) = ZMANSELL2(:,:)
  DEALLOCATE(ZMANSELL1)
  DEALLOCATE(ZMANSELL2)
!
  XMANSELL(:,:) = XMANSELL(:,:) * 1.E-15 ! in C
END IF
!
!
!*      8.4    Saunders et al. (1991) parameterization
!              Idem for Brooks et al. (1997), but with EW = ZRAR/3.
!
!
IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR.  &
    CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2')  THEN
!
  NIND_TEMP = 31
  NIND_LWC  = 28
!
  IF( .NOT.ALLOCATED(XSAUNDER)) ALLOCATE(XSAUNDER(NIND_LWC+1,NIND_TEMP+1))
  ALLOCATE(ZT(NIND_TEMP+1))    ! Kelvin
  ALLOCATE(ZLWCC(NIND_TEMP+1))
  DO JTEMP = 1, NIND_TEMP+1
    ZT(JTEMP) = 1.0 - REAL(JTEMP) + XTT
  END DO
  ZLWCC(:) = MIN( MAX( -0.49 + 6.64E-2*(XTT-ZT(:)),0.22 ),1.1 )   ! (g m^-3)
  ALLOCATE(ZEW(NIND_LWC+1))
!
!                       LWC index (0.01 g.m^-3 --> 10 g.m^-3)
!                                  0.01 to 0.09 every 0.01 (9 values)
!                                  0.10 to 0.90 every 0.10 (9 values)
!                                  1.00 to 10.0 every 1.00 (10 values)
  DO JLWC = 1, 9
    ZEW(JLWC) = 0.01 * REAL(JLWC)
  END DO
  DO JLWC = 10, 18
    ZEW(JLWC) = 0.1 + 0.1 * REAL(JLWC-10)
  END DO
  DO JLWC = 19, NIND_LWC+1
    ZEW(JLWC) = 1.0 + REAL(JLWC-19)
  END DO
!
!
  XSAUNDER(:,:) = 0.0
  DO JTEMP = 1, NIND_TEMP+1
    DO JLWC = 1, NIND_LWC+1
!
! region S4 : positive
      IF (ZT(JTEMP) <= (XTT-7.35) .AND. ZT(JTEMP) > (XTT-23.9458) .AND. &
          ZEW(JLWC) > ZLWCC(JTEMP)) THEN
        XSAUNDER(JLWC,JTEMP) = MAX( 0.,                                        &
                                    20.22*ZEW(JLWC)+1.36*(ZT(JTEMP)-XTT)+10.05 )
      ENDIF
!
! region S1 : positive --> linear interpolation
      IF (ZT(JTEMP) > (XTT-7.35) .AND. ZT(JTEMP) < XTT .AND. &
          ZEW(JLWC) > ZLWCC(JTEMP)) THEN
        XSAUNDER(JLWC,JTEMP) = MAX( 0.,-(2.75*ZEW(JLWC)+0.007)*(ZT(JTEMP)-XTT) )
      ENDIF
!
! region S8 : positive
      IF (ZT(JTEMP) <= (XTT-23.9458) .AND. ZT(JTEMP) > (XTT-40.0) .AND. &
          ZEW(JLWC) > ZLWCC(JTEMP)) THEN
        XSAUNDER(JLWC,JTEMP)  = MAX( 0.,20.22*ZEW(JLWC)-22.26 )
      ENDIF
!
! region S7 : negative
      IF (ZT(JTEMP) <= (XTT-7.35) .AND. ZT(JTEMP) > (XTT-40.0) .AND. &
          ZEW(JLWC) >= 0.104149   .AND. ZEW(JLWC) < ZLWCC(JTEMP)) THEN
        XSAUNDER(JLWC,JTEMP) = MIN( 0.,3.02-31.76*ZEW(JLWC)+26.53*ZEW(JLWC)**2 )
      ENDIF
    END DO
  END DO
END IF
!
! SAUN1 doesn't take into account marginal positive and negative regions at
! low LWC
!
IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'BSMP1') THEN
  DO JTEMP = 1, NIND_TEMP+1
    DO JLWC = 1, NIND_LWC+1
!
! region S1 : negative --> linear interpolation
      IF (ZT(JTEMP) > (XTT-7.35)   .AND. ZT(JTEMP) < XTT .AND. &
          ZEW(JLWC) < ZLWCC(JTEMP) .AND. ZEW(JLWC) >= 0.104149) THEN
        XSAUNDER(JLWC,JTEMP) = MIN( 0.,                                        &
                      (-0.41+4.32*ZEW(JLWC)-3.61*ZEW(JLWC)**2)*(ZT(JTEMP)-XTT) )
      ENDIF
    END DO
  END DO
!
  XSAUNDER(:,:) = XSAUNDER(:,:) * 1.E-15 ! in C
!
END IF
!
! SAUN2 takes into account marginal positive and negative regions at low LWC
!
IF (CNI_CHARGING == 'SAUN2' .OR. CNI_CHARGING == 'BSMP2') THEN
  DO JTEMP = 1, NIND_TEMP+1
    DO JLWC = 1, NIND_LWC+1
!
! region S2 : negative
      IF (ZT(JTEMP) <= (XTT-7.35) .AND. ZT(JTEMP) > (XTT-16.0) .AND. &
          ZEW(JLWC) >= 0.026      .AND. ZEW(JLWC) < 0.14) THEN
        XSAUNDER(JLWC,JTEMP) = MIN( 0.,-314.4*ZEW(JLWC) + 7.92 )
      ENDIF
!
! region S3 : negative
      IF (ZT(JTEMP) <= (XTT-7.35) .AND. ZT(JTEMP) > (XTT-16.0) .AND. &
          ZEW(JLWC) >= 0.14       .AND. ZEW(JLWC) < 0.22) THEN
        XSAUNDER(JLWC,JTEMP) = MIN( 0.,419.4 * ZEW(JLWC) - 92.64 )
      ENDIF
!
! region S5 : positive
      IF (ZT(JTEMP) < (XTT-20.0) .AND. ZT(JTEMP) > (XTT-40.0) .AND. &
          ZEW(JLWC) >= 0.063034  .AND. ZEW(JLWC) < 0.12) THEN
        XSAUNDER(JLWC,JTEMP) = MAX( 0.,2041.76*ZEW(JLWC) - 128.7 )
      ENDIF
!
! region S6 : positive
      IF (ZT(JTEMP) < (XTT-20.0) .AND. ZT(JTEMP) > (XTT-40.0) .AND. &
          ZEW(JLWC)  >= 0.12      .AND. ZEW(JLWC)  < 0.1596) THEN
        XSAUNDER(JLWC,JTEMP) = MAX( 0.,-2900.22*ZEW(JLWC) + 462.91 )
      ENDIF
!
! region S1 : negative --> linear interpolation of S3
      IF (ZT(JTEMP) > (XTT-7.35) .AND. ZT(JTEMP) < XTT .AND. &
          ZEW(JLWC)  >= 0.14      .AND. ZEW(JLWC)  < ZLWCC(JTEMP)) THEN
        XSAUNDER(JLWC,JTEMP) = MIN( 0.,(-57.06*ZEW(JLWC)+12.6)*(ZT(JTEMP)-XTT) )
      ENDIF
!
! region S1 : negative --> linear interpolation of S2
      IF (ZT(JTEMP) > (XTT-7.35) .AND. ZT(JTEMP) < XTT .AND. &
          ZEW(JLWC)  >= 0.026     .AND. ZEW(JLWC)  < 0.14) THEN
        XSAUNDER(JLWC,JTEMP) = MIN( 0.,(42.8*ZEW(JLWC)-1.08)*(ZT(JTEMP)-XTT) )
      ENDIF
    END DO
  END DO
!
  XSAUNDER(:,:) = XSAUNDER(:,:) * 1.E-15 ! in C
!
END IF
!
!*      8.5    Takahashi with EW or ZRAR (Tsenova and Mitzeva, 2009, 2011)
!                       here ZRAR = 9 * EW
!                       Temperature index (0C --> -30C)
!                       LWC index (0.01 g.m^-3 --> 10 g.m^-3)
!                                  0.01 to 0.09 every 0.01 (9 values)
!                                  0.10 to 0.90 every 0.10 (9 values)
!                                  1.00 to 10.0 every 1.00 (10 values)
!
IF (CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
!
  NIND_TEMP = 31
  NIND_LWC  = 28
!
  IF( .NOT.ALLOCATED(XTAKA_TM)) ALLOCATE(XTAKA_TM(NIND_LWC+1,NIND_TEMP+1))
  ALLOCATE(ZT(NIND_TEMP+1))    ! Kelvin
  ALLOCATE(ZEW(NIND_LWC+1))
  DO JTEMP = 1, NIND_TEMP+1
    ZT(JTEMP) = 1.0 - REAL(JTEMP) + XTT
  END DO

  DO JLWC = 1, 9
    ZEW(JLWC) = 0.01 * REAL(JLWC)
  END DO
  DO JLWC = 10, 18
    ZEW(JLWC) = 0.1 + 0.1 * REAL(JLWC-10)
  END DO
  DO JLWC = 19, NIND_LWC+1
    ZEW(JLWC) = 1.0 + REAL(JLWC-19)
  END DO
!
  XTAKA_TM(:,:) = 0.0
  DO JTEMP = 1, NIND_TEMP+1
    DO JLWC = 1, NIND_LWC+1
!
! Eq. 1: >0
      IF ( ZT(JTEMP) > (XTT - 10.) .AND. ZEW(JLWC) <= 1.6) THEN
        XTAKA_TM(JLWC, JTEMP) = 146.981 * ZEW(JLWC) - 116.37 * ZEW(JLWC)**2  &
                               + 29.76 * ZEW(JLWC)**3                        &
                               - 0.03 * (ZT(JTEMP) - XTT)**3 * ZEW(JLWC)     &
                               - 2.58 * (ZT(JTEMP) - XTT)                    &
                               - 0.21 * (ZT(JTEMP) - XTT)**3 * ZEW(JLWC)**3  &
                               + 0.36 * (ZT(JTEMP) - XTT)**3 * ZEW(JLWC)**2  &
                               + 0.15 * (ZT(JTEMP) - XTT)**2                 &
                               + 2.92 * (ZT(JTEMP) - XTT)  * ZEW(JLWC)**3    &
                               - 4.22 * (ZT(JTEMP) - XTT)  * ZEW(JLWC) - 8.506
      END IF
!
!  Eq. 2: >0
      IF ( ZT(JTEMP) > (XTT - 10.) .AND. &
             ZEW(JLWC) > 1.6 .AND. ZEW(JLWC) <= 8.) THEN
        XTAKA_TM(JLWC, JTEMP) = 4.179 * (ZT(JTEMP) - XTT)                    &
                               - 0.005 * (ZT(JTEMP) - XTT)**2 * ZEW(JLWC)**2 &
                               + 0.916 * ZEW(JLWC)**2                        &
                               - 1.333 * (ZT(JTEMP) - XTT)    * ZEW(JLWC)    &
                               - 7.465 * ZEW(JLWC)                           &
                               + 0.109 * (ZT(JTEMP) - XTT)    * ZEW(JLWC)**2 &
                               + 0.001 * (ZT(JTEMP) - XTT)**2 * ZEW(JLWC)**3 &
                               - 0.035 * ZEW(JLWC)**3  + 50.84454
      END IF
!
!  Eq. 8: > 0
      IF ( ZEW(JLWC) <= 0.4 .AND. &
              ZT(JTEMP) <= (XTT - 10.) .AND. ZT(JTEMP) >= (XTT - 40.)) THEN
        XTAKA_TM(JLWC, JTEMP) = - 3.3515 * (ZT(JTEMP) - XTT)                   &
                                + 95.957 * (ZT(JTEMP) - XTT)    * ZEW(JLWC)**2 &
                                + 511.83 * ZEW(JLWC)                           &
                                + 17.448 * (ZT(JTEMP) - XTT)**2 * ZEW(JLWC)**3 &
                                - 0.0007 * (ZT(JTEMP) - XTT)**3                &
                                + 20.570 * (ZT(JTEMP) - XTT)    * ZEW(JLWC)    &
                                + 0.1656 * (ZT(JTEMP) - XTT)**2 * ZEW(JLWC)    &
                                + 0.4954 * (ZT(JTEMP) - XTT)**3 * ZEW(JLWC)**3 &
                                - 0.0975 * (ZT(JTEMP) - XTT)**3 * ZEW(JLWC)**2 &
                                + 67.457 * (ZT(JTEMP) - XTT)    * ZEW(JLWC)**3 &
                                - 0.1066 * (ZT(JTEMP) - XTT)**2 - 24.5715
      END IF
!
! Eq. 9: < 0
      IF ( ZT(JTEMP) <= (XTT - 10.) .AND. ZT(JTEMP) >= (XTT - 40.) .AND. &
             ZEW(JLWC) > 0.4 .AND. ZEW(JLWC) <= 3.2) THEN
        XTAKA_TM(JLWC, JTEMP) = - 1.5676 * (ZT(JTEMP) - XTT) * ZEW(JLWC)      &
                                + 0.2484 * (ZT(JTEMP) - XTT)   * ZEW(JLWC)**3 &
                                + 0.0112 * (ZT(JTEMP) - XTT)**3               &
                                + 19.199 * (ZT(JTEMP) - XTT)                  &
                                + 0.8051 * (ZT(JTEMP) - XTT)**2               &
                                - 83.4 * ZEW(JLWC)                            &
                                + 15.4 * ZEW(JLWC)**2                         &
                                + 5.97 * ZEW(JLWC)**3 + 167.9278
      END IF
!
!  Eq. 10: > 0
      IF ( ZT(JTEMP) <= (XTT - 10.) .AND. ZT(JTEMP) >= (XTT - 40.) .AND. &
           ZEW(JLWC) > 3.2 .AND. ZEW(JLWC) <= 8. ) THEN
        XTAKA_TM(JLWC, JTEMP) = 4.2127 * (ZT(JTEMP) - XTT)                 &
                              - 0.8311 * (ZT(JTEMP) - XTT) * ZEW(JLWC)     &
                              + 0.0670 * (ZT(JTEMP) - XTT) * ZEW(JLWC) **2 &
                              + 0.0042 * (ZT(JTEMP) - XTT)**2 * ZEW(JLWC)  &
                              + 40.9642
      END IF
    END DO
  END DO
!
  XTAKA_TM(:,:) = XTAKA_TM(:,:) * 1.E-15 ! in C
!
END IF
!
!
!-------------------------------------------------------------------------------
!
!*	9.	NON_INDUCTIVE PROCESS: AGGREGATION OF ICE ON SNOW
!		-------------------------------------------------
!
!*      9.1     Helsdon and Farley (1987) parameterization
!
XFQIAGGSBH = 2.E-14  ! (C.) Constant for ice-snow charging process 
!
!
!*      9.2     Gardiner et al. (1985) parameterization
!
IF (HELEC == 'ELE4') THEN
  XFQIAGGSBG = (XPI / 4.0) * ZCS**4. * PRHO00**(4. * ZCEXVT) * &
                MOMG(ZALPHAS,ZNUS,2.+4.*ZDS) * 7.3 *          &
                MOMG(ZALPHAI,ZNUI,4.)
ELSE
  XFQIAGGSBG = (XPI / 4.0) * ZCCS * ZCS**4. * PRHO00**(4. * ZCEXVT) * &
                MOMG(ZALPHAS,ZNUS,2.+4.*ZDS) * 7.3 *          &
                MOMG(ZALPHAI,ZNUI,4.)
END IF
!
!
!*      9.3     Saunders et al.(1991) parameterization
!
IF (HELEC == 'ELE4') THEN
  XFQIAGGSBS = XPI / 4.0
ELSE
  XFQIAGGSBS = (XPI / 4.0) * ZCCS
END IF
!
!
!*      9.4     Takahashi (1978) parameterization
!
IF (CNI_CHARGING == 'TAKAH') THEN
  IF (HELEC == 'ELE4') THEN
    XFQIAGGSBT1 = (XPI / 4.0) * ZCS
  ELSE
    XFQIAGGSBT1 = (XPI / 4.0) * ZCCS * ZCS
  END IF
  XFQIAGGSBT2 = 10. * MOMG(ZALPHAS,ZNUS,2.+ZDS)
  XFQIAGGSBT3 = 5. * ZCS * MOMG(ZALPHAI,ZNUI,2.) *               &
                MOMG(ZALPHAS,ZNUS,2.+2.*ZDS) / ((1.E-4)**2 * 8. * & 
                (ZAI * MOMG(ZALPHAI,ZNUI,ZBI))**(2. / ZBI))  
END IF
!
!
!-------------------------------------------------------------------------------
!
!*	10.	ACCRETION OF RAINDROPS ON SNOW
!		------------------------------
!
IF( .NOT.ALLOCATED(XKER_Q_RACCSS)) ALLOCATE( XKER_Q_RACCSS(IACCLBDAS,IACCLBDAR) )
IF( .NOT.ALLOCATED(XKER_Q_RACCS)) ALLOCATE( XKER_Q_RACCS (IACCLBDAS,IACCLBDAR) )
IF( .NOT.ALLOCATED(XKER_Q_SACCRG)) ALLOCATE( XKER_Q_SACCRG(IACCLBDAR,IACCLBDAS) )
!
IF (HELEC == 'ELE4') THEN
  XFQRACCS = (XPI / 4.0) * PRHO00**ZCEXVT
ELSE
  XFQRACCS = (XPI / 4.0) * ZCCS * ZCCR * (PRHO00**ZCEXVT)
END IF
!
XLBQRACCS1  =      MOMG(ZALPHAR,ZNUR,2.+XFR)
XLBQRACCS2  = 2. * MOMG(ZALPHAR,ZNUR,1.+XFR) * MOMG(ZALPHAS,ZNUS,1.)
XLBQRACCS3  =      MOMG(ZALPHAR,ZNUR,XFR)    * MOMG(ZALPHAS,ZNUS,2.)
!
XLBQSACCRG1 =      MOMG(ZALPHAS,ZNUS,2.+XFS)
XLBQSACCRG2 = 2. * MOMG(ZALPHAS,ZNUS,1.+XFS) * MOMG(ZALPHAR,ZNUR,1.)
XLBQSACCRG3 =      MOMG(ZALPHAS,ZNUS,XFS)    * MOMG(ZALPHAR,ZNUR,2.)
!
! These values are pasted from ini_rain_ice (7.2.2)
IND      = 50   ! Interval number, collection efficiency and infinite diameter
ZESR     = 1.0  ! factor used to integrate the dimensional distributions when
ZFDINFTY = 20.0 ! computing the kernels XKER_RACCSS, XKER_RACCS and XKER_SACCRG
!
CALL RRCOLSS (IND, ZALPHAS, ZNUS, ZALPHAR, ZNUR,                          &
              ZESR, XFR, ZCS, ZDS, 0., ZCR, ZDR,                          &
              ZACCLBDAS_MAX, ZACCLBDAR_MAX, ZACCLBDAS_MIN, ZACCLBDAR_MIN, &
              ZFDINFTY, XKER_Q_RACCSS, ZAG, ZBS, ZAS                      )
!
CALL RZCOLX  (IND, ZALPHAS, ZNUS, ZALPHAR, ZNUR,                          &
              ZESR, XFR, ZCS, ZDS, 0., ZCR, ZDR, 0.,                      &
              ZACCLBDAS_MAX, ZACCLBDAR_MAX, ZACCLBDAS_MIN, ZACCLBDAR_MIN, &
              ZFDINFTY, XKER_Q_RACCS                                      )
!
CALL RSCOLRG (IND, ZALPHAS, ZNUS, ZALPHAR, ZNUR,                          &
              ZESR, XFS, ZCS, ZDS, 0., ZCR, ZDR,                          &
              ZACCLBDAS_MAX, ZACCLBDAR_MAX, ZACCLBDAS_MIN, ZACCLBDAR_MIN, &
              ZFDINFTY, XKER_Q_SACCRG, ZAG, ZBS, ZAS                      )
!
!-------------------------------------------------------------------------------
!
!*	11.	DRY GROWTH OF GRAUPELN BY CAPTURE OF SNOW OR ICE
!		------------------------------------------------
!
!*      11.1    charge transfer associated to mass transfer
!
IF( .NOT.ALLOCATED(XKER_Q_SDRYG)) ALLOCATE( XKER_Q_SDRYG(IDRYLBDAG,IDRYLBDAS) )
!
IF (HELEC == 'ELE4') THEN
  XFQSDRYG = (XPI / 4.0) * (PRHO00**ZCEXVT)
ELSE
  XFQSDRYG = (XPI / 4.0) * ZCCS * ZCCG * (PRHO00**ZCEXVT)
END IF
!
XLBQSDRYG1 =      MOMG(ZALPHAS,ZNUS,2.+XFS)
XLBQSDRYG2 = 2. * MOMG(ZALPHAS,ZNUS,1.+XFS) * MOMG(ZALPHAG,ZNUG,1.)
XLBQSDRYG3 =      MOMG(ZALPHAS,ZNUS,XFS)    * MOMG(ZALPHAG,ZNUG,2.)
!
ZEGS = 1. ! also initialized in ini_rain_ice_elec
!
CALL RZCOLX (IND, ZALPHAG, ZNUG, ZALPHAS, ZNUS,                          &
             ZEGS, XFS, ZCG, ZDG, 0., ZCS, ZDS, 0.,                      &
             ZDRYLBDAG_MAX, ZDRYLBDAS_MAX, ZDRYLBDAG_MIN, ZDRYLBDAS_MIN, &
             ZFDINFTY, XKER_Q_SDRYG                                      )
!
!
!*      11.2    NI process: Heldson et Farley (1987) parameterization
!
IF (CNI_CHARGING == 'HELFA') THEN
  XHIDRYG = 2.E-15   ! Charge exchanged per collision between ice and graupel
  XHSDRYG = 2.E-14
!
  IF (HELEC == 'ELE4') THEN
    XFQSDRYGBH   = (XPI / 4.0) * PRHO00**(ZCEXVT) * XHSDRYG 
  ELSE
    XFQSDRYGBH   = (XPI / 4.0) * ZCCG * ZCCS * (PRHO00**(ZCEXVT)) * XHSDRYG
  END IF
! 
  XLBQSDRYGB4H =      MOMG(ZALPHAS,ZNUS,2.)  
  XLBQSDRYGB5H = 2. * MOMG(ZALPHAS,ZNUS,1.) * MOMG(ZALPHAG,ZNUG,1.) 
  XLBQSDRYGB6H =      MOMG(ZALPHAG,ZNUG,2.) 
!
  IF( .NOT.ALLOCATED(XKER_Q_SDRYGB)) ALLOCATE( XKER_Q_SDRYGB(IDRYLBDAG,IDRYLBDAS) )
  CALL RZCOLX (IND, ZALPHAG, ZNUG, ZALPHAS, ZNUS,                          &
               ZEGS, 0., ZCG, ZDG, 0., ZCS, ZDS, 0.,                       &
               ZDRYLBDAG_MAX, ZDRYLBDAS_MAX, ZDRYLBDAG_MIN, ZDRYLBDAS_MIN, &
               ZFDINFTY, XKER_Q_SDRYGB                                     )
! Delta vqb1_sg
ENDIF
!
!
!*      11.3    NI process: Gardiner et al. (1985) parameterization
!
IF (CNI_CHARGING == 'GARDI') THEN
  IF (HELEC == 'ELE4') THEN
    XFQIDRYGBG   = (XPI / 4.0) * PRHO00**(4.*ZCEXVT) * ZCG**4. * 7.3
    XFQSDRYGBG   = (XPI / 4.0) * PRHO00**(4.*ZCEXVT) * 7.3
  ELSE
    XFQIDRYGBG   = (XPI / 4.0) * ZCCG * (PRHO00**(4.*ZCEXVT)) * ZCG**4. * 7.3
    XFQSDRYGBG   = (XPI / 4.0) * ZCCS * ZCCG * (PRHO00**(4.*ZCEXVT)) * 7.3
  END IF
  XLBQIDRYGBG  = MOMG(ZALPHAI,ZNUI,4.) * MOMG(ZALPHAG,ZNUG,2.+4.*ZDG)
!
  XLBQSDRYGB4G =      MOMG(ZALPHAS,ZNUS,4.) * MOMG(ZALPHAG,ZNUG,2.)
  XLBQSDRYGB5G = 2. * MOMG(ZALPHAS,ZNUS,5.) * MOMG(ZALPHAG,ZNUG,1.)
  XLBQSDRYGB6G =      MOMG(ZALPHAS,ZNUS,6.)    
!
  IF( .NOT.ALLOCATED(XKER_Q_SDRYGB)) ALLOCATE( XKER_Q_SDRYGB(IDRYLBDAG,IDRYLBDAS) )
  CALL VQZCOLX (IND, ZALPHAG, ZNUG, ZALPHAS, ZNUS,                          &
                ZEGS, 4., ZCG, ZDG, ZCS, ZDS, 4.,                           &
                ZDRYLBDAG_MAX, ZDRYLBDAS_MAX, ZDRYLBDAG_MIN, ZDRYLBDAS_MIN, &
                ZFDINFTY, XKER_Q_SDRYGB                                     )
END IF
!
!
!*      11.4    NI process: Saunders et al. (1991) and 
!*                          Saunders and Peck (1998) parameterizations
!
IF (CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
    CNI_CHARGING == 'SAP98' .OR.                              &
    CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
    CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
  IF (HELEC == 'ELE4') THEN
    XFQIDRYGBS   = XPI / 4.0
    XFQSDRYGBS   = XPI / 4.0
  ELSE
    XFQIDRYGBS   = (XPI / 4.0) * ZCCG
    XFQSDRYGBS   = (XPI / 4.0) * ZCCS * ZCCG
  END IF
  XLBQSDRYGB1S = MOMG(ZALPHAG,ZNUG,2.)
  XLBQSDRYGB2S = 2. * MOMG(ZALPHAG,ZNUG,1.)
!
  IF( .NOT.ALLOCATED(XKER_Q_SDRYGB1)) ALLOCATE( XKER_Q_SDRYGB1(IDRYLBDAG,IDRYLBDAS) )
  IF( .NOT.ALLOCATED(XKER_Q_SDRYGB2)) ALLOCATE( XKER_Q_SDRYGB2(IDRYLBDAG,IDRYLBDAS) )
!
! Positive charging region
  CALL VQZCOLX (IND, ZALPHAG, ZNUG, ZALPHAS, ZNUS,                          &
                ZEGS, XSMP, ZCG, ZDG, ZCS, ZDS, (1.+XSNP),                  &
                ZDRYLBDAG_MAX, ZDRYLBDAS_MAX, ZDRYLBDAG_MIN, ZDRYLBDAS_MIN, &
                ZFDINFTY, XKER_Q_SDRYGB1                                    )
!
! Negative charging region
  CALL VQZCOLX (IND, ZALPHAG, ZNUG, ZALPHAS, ZNUS,                          &
                ZEGS, XSMN, ZCG, ZDG, ZCS, ZDS, (1.+XSNN),                  &
                ZDRYLBDAG_MAX, ZDRYLBDAS_MAX, ZDRYLBDAG_MIN, ZDRYLBDAS_MIN, &
                ZFDINFTY, XKER_Q_SDRYGB2                                    )
ENDIF
!
!
!*      11.5    NI process: Takahashi (1978) parameterization
!
IF (CNI_CHARGING == 'TAKAH') THEN
!
! IDRYG_boun
  IF (HELEC == 'ELE4') THEN
    XFQIDRYGBT1 = (XPI / 4.0) * ZCG
  ELSE
    XFQIDRYGBT1 = (XPI / 4.0) * ZCCG * ZCG
  END IF
  XFQIDRYGBT2 = 10.0 * MOMG(ZALPHAG,ZNUG,2.+ZDG) 
  XFQIDRYGBT3 = 5.0 * ZCG * MOMG(ZALPHAI,ZNUI,2.) *               &
                MOMG(ZALPHAG,ZNUG,2.+2.*ZDG) / ((2.E-4)**2 * 8. * & 
               (ZAI * MOMG(ZALPHAI,ZNUI,ZBI))**(2./ZBI))  
!
! SDRYG_boun
  IF (HELEC == 'ELE4') THEN
    XFQSDRYGBT1  = XPI / 4.0
  ELSE
    XFQSDRYGBT1  = (XPI / 4.0) * ZCCG * ZCCS
  END IF
  XFQSDRYGBT2  = ZCG * MOMG(ZALPHAG,ZNUG,ZDG) * MOMG(ZALPHAS,ZNUS,2.)
  XFQSDRYGBT3  = ZCS * MOMG(ZALPHAS,ZNUS,2.+ZDS)
  XFQSDRYGBT4  = ZCG * MOMG(ZALPHAG,ZNUG,2.+ZDG)
  XFQSDRYGBT5  = ZCS * MOMG(ZALPHAG,ZNUG,2.) * MOMG(ZALPHAS,ZNUS,ZDS)
  XFQSDRYGBT6  = 2. * ZCG * MOMG(ZALPHAG,ZNUG,1.+ZDG) * MOMG(ZALPHAS,ZNUS,1.)
  XFQSDRYGBT7  = 2. * ZCS * MOMG(ZALPHAG,ZNUG,1.) * MOMG(ZALPHAS,ZNUS,1.+ZDS)
  XFQSDRYGBT8  = 5. / ((1.E-4)**2 * 8.)
  XFQSDRYGBT9  = MOMG(ZALPHAG,ZNUG,2.) * MOMG(ZALPHAS,ZNUS,2.)
  XFQSDRYGBT10 = MOMG(ZALPHAS,ZNUS,4.)
  XFQSDRYGBT11 = 2. * MOMG(ZALPHAG,ZNUG,1.) * MOMG(ZALPHAS,ZNUS,3.)
!
  IF( .NOT.ALLOCATED(XKER_Q_SDRYGB)) ALLOCATE( XKER_Q_SDRYGB(IDRYLBDAG,IDRYLBDAS) )
  CALL VQZCOLX (IND, ZALPHAG, ZNUG, ZALPHAS, ZNUS,                          &
                ZEGS, 2., ZCG, ZDG, ZCS, ZDS, 2.,                           &
                ZDRYLBDAG_MAX, ZDRYLBDAS_MAX, ZDRYLBDAG_MIN, ZDRYLBDAS_MIN, &
                ZFDINFTY, XKER_Q_SDRYGB                                     )
END IF
!
!
!*      11.6    NI process: limit the charge exchanged during QSDRYG_boun
!
IF (CNI_CHARGING == 'TAKAH' .OR. CNI_CHARGING == 'SAP98' .OR. &
    CNI_CHARGING == 'SAUN1' .OR. CNI_CHARGING == 'SAUN2' .OR. &
    CNI_CHARGING == 'GARDI' .OR.                              &
    CNI_CHARGING == 'BSMP1' .OR. CNI_CHARGING == 'BSMP2' .OR. &
    CNI_CHARGING == 'TEEWC' .OR. CNI_CHARGING == 'TERAR') THEN
  IF (HELEC == 'ELE4') THEN
    XAUX_LIM  = XPI / 4.0
  ELSE
    XAUX_LIM  = (XPI / 4.0) * ZCCG * ZCCS
  END IF
  XAUX_LIM1 =      MOMG(ZALPHAS,ZNUS,2.)  
  XAUX_LIM2 = 2. * MOMG(ZALPHAS,ZNUS,1.) * MOMG(ZALPHAG,ZNUG,1.) 
  XAUX_LIM3 =      MOMG(ZALPHAG,ZNUG,2.)
  IF( .NOT.ALLOCATED(XKER_Q_LIMSG)) ALLOCATE( XKER_Q_LIMSG(IDRYLBDAG,IDRYLBDAS) )
  CALL RZCOLX (IND, ZALPHAG, ZNUG, ZALPHAS, ZNUS,                          &
               ZEGS, 0., ZCG, ZDG, 0., ZCS, ZDS, 0.,                       &
               ZDRYLBDAG_MAX, ZDRYLBDAS_MAX, ZDRYLBDAG_MIN, ZDRYLBDAS_MIN, &
               ZFDINFTY, XKER_Q_LIMSG)
ENDIF
!
!
!-------------------------------------------------------------------------------
!
!*	12.	DRY GROWTH OF GRAUPELN BY CAPTURE OF RAINDROP
!		---------------------------------------------
!
IF( .NOT.ALLOCATED(XKER_Q_RDRYG)) ALLOCATE( XKER_Q_RDRYG(IDRYLBDAG,IDRYLBDAR) ) 
!
IF (HELEC == 'ELE4') THEN
  XFQRDRYG = (XPI / 4.0) * PRHO00**ZCEXVT
ELSE
  XFQRDRYG = (XPI / 4.0) * ZCCG * ZCCR * (PRHO00**ZCEXVT)
END IF
!
XLBQRDRYG1 =      MOMG(ZALPHAR,ZNUR,2.+XFR)
XLBQRDRYG2 = 2. * MOMG(ZALPHAR,ZNUR,1.+XFR) * MOMG(ZALPHAG,ZNUG,1.)
XLBQRDRYG3 =      MOMG(ZALPHAR,ZNUR,XFR)    * MOMG(ZALPHAG,ZNUG,2.)
!
ZEGR = 1.0 
!
CALL RZCOLX (IND, ZALPHAG, ZNUG, ZALPHAR, ZNUR,                            & 
             ZEGR, XFR, ZCG, ZDG, 0., ZCR, ZDR, 0.,                        &
             ZDRYLBDAG_MAX, ZDRYLBDAR_MAX, ZDRYLBDAG_MIN, ZDRYLBDAR_MIN,   &
             ZFDINFTY, XKER_Q_RDRYG                                        )
!
!
!-------------------------------------------------------------------------------
!
!*	13.	UPDATE THE Q=f(D) RELATION
!		--------------------------
!
IF (HCLOUD(1:3) == 'ICE') THEN
  XFQUPDC   = 400.E6 * MOMG(XALPHACQ,XNUCQ,XFC) / XLBDACQ**XFC ! Nc~400E6 m-3 as 
                                                          ! proposed for RCHONI
  !
  XFQUPDR   = ZCCR * MOMG(ZALPHAR,ZNUR,XFR)
  XEXFQUPDI = XFI / ZBI
  XFQUPDI   = MOMG(ZALPHAI,ZNUI,XFI) * (ZAI * MOMG(ZALPHAI,ZNUI,ZBI))**(-XEXFQUPDI)
END IF
XFQUPDS   = ZCCS * MOMG(ZALPHAS,ZNUS,XFS)
XFQUPDG   = ZCCG * MOMG(ZALPHAG,ZNUG,XFG)
XFQUPDH   = ZCCH * MOMG(ZALPHAH,ZNUH,XFH)
!
!
!------------------------------------------------------------------------------
!
!*      14.     INDUCTIVE PROCESS
!               -----------------
!
! d = 15 microns and N_c = 400 cm**(-3)
!
XCOLCG_IND = 0.8
XEBOUND    = 0.1
XALPHA_IND = 0.07    ! moderate inductive charging
XCOS_THETA = 0.2  
!
IF (HELEC == 'ELE4') THEN
  XIND1 = (XPI**3 / 8.) * (15.E-6)**2 * &
           ZCG * 400.E6 *               &
          XCOLCG_IND * XEBOUND * XALPHA_IND
ELSE
  XIND1 = (XPI**3 / 8.) * (15.E-6)**2 * &
           ZCG * 400.E6 * ZCCG *        &
           XCOLCG_IND * XEBOUND * XALPHA_IND
END IF
XIND2 = XPI * XEPSILON * XCOS_THETA * MOMG(ZALPHAG,ZNUG,2.+ZDG)
XIND3 = MOMG(ZALPHAG,ZNUG,ZDG+XFG) / 3.
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_PARAM_ELEC
