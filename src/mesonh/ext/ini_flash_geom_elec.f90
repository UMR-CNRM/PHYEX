!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ###############################
      MODULE MODI_INI_FLASH_GEOM_ELEC
!     ###############################
!
INTERFACE
!
      SUBROUTINE INI_FLASH_GEOM_ELEC (HCLOUD)
!
CHARACTER(LEN=4), INTENT(IN) :: HCLOUD   ! microphysics scheme
!
END SUBROUTINE INI_FLASH_GEOM_ELEC
END INTERFACE
END MODULE MODI_INI_FLASH_GEOM_ELEC
!
!	#######################################
        SUBROUTINE INI_FLASH_GEOM_ELEC (HCLOUD)
!	#######################################
!
!!****  *INI_FLASH_GEOM_ELEC* - routine to initialize the lightning flashes
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize the variables
!     of the lightning flashes routine 
!
!!**  METHOD
!!    ------
!!      The initialization of the scheme is performed as follows :
!!   
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!    REFERENCE
!!    ---------
!!      
!!    AUTHOR
!!    ------
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     29/11/02
!!
!!      Modifications
!!        J.-P. Pinty  jan 2015 : add LMA simulator
!!        J.Escobar 20/06/2018 : truly set NBRANCH_MAX = 5000 !
!!        C. Barthe 30/11/2022 : add parameters for LIMA
!!        C. Barthe 11/09/2023 : modify some parameters to use with LIMA2 
!!
!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
USE MODD_CST, ONLY : XPI
USE MODD_RAIN_ICE_DESCR_n,ONLY : XALPHAR_I=>XALPHAR, XNUR_I=>XNUR, XCCR,                       &
                                 XALPHAI_I=>XALPHAI, XNUI_I=>XNUI, XAI_I=>XAI, XBI_I=>XBI,     &
                                 XALPHAS_I=>XALPHAS, XNUS_I=>XNUS, XCCS_I=>XCCS, XCXS_I=>XCXS, &
                                 XALPHAG_I=>XALPHAG, XNUG_I=>XNUG, XCCG_I=>XCCG, XCXG_I=>XCXG, &
                                 XALPHAH_I=>XALPHAH, XNUH_I=>XNUH, XCCH_I=>XCCH, XCXH_I=>XCXH
USE MODD_PARAM_LIMA,      ONLY : XALPHAC, XNUC,                                                      &
                                 XALPHAR_L=>XALPHAR, XNUR_L=>XNUR, XALPHAI_L=>XALPHAI, XNUI_L=>XNUI, & 
                                 XALPHAS_L=>XALPHAS, XNUS_L=>XNUS, XALPHAG_L=>XALPHAG, XNUG_L=>XNUG, &
                                 NMOM_S, NMOM_G, NMOM_H
USE MODD_PARAM_LIMA_COLD, ONLY : XAI_L=>XAI, XBI_L=>XBI, XCCS_L=>XCCS, XCXS_L=>XCXS
USE MODD_PARAM_LIMA_MIXED,ONLY : XCCG_L=>XCCG, XCXG_L=>XCXG, &
                                 XCCH_L=>XCCH, XCXH_L=>XCXH, XALPHAH_L=>XALPHAH, XNUH_L=>XNUH
USE MODD_ELEC_DESCR
USE MODD_ELEC_PARAM
USE MODD_DIM_n, ONLY : NKMAX
USE MODD_TIME_n, ONLY : TDTCUR
USE MODD_LMA_SIMULATOR, ONLY : LLMA, TDTLMA, LWRITE_LMA, XDTLMA, CLMA_FILE
!
USE MODI_MOMG
!
IMPLICIT NONE
!
!*	0.1	Declaration of dummy arguments
!
CHARACTER(LEN=4), INTENT(IN) :: HCLOUD   ! microphysics scheme
!
!
!*	0.2	Declaration of local variables
!
! variables used to cope with the module variables common to icex and lima
REAL :: ZALPHAR, ZNUR,             &
        ZAI, ZBI, ZALPHAI, ZNUI,   &
        ZCCS, ZCXS, ZALPHAS, ZNUS, &
        ZCCG, ZCXG, ZALPHAG, ZNUG, &
        ZCCH, ZCXH, ZALPHAH, ZNUH
!
!-------------------------------------------------------------------------------
!
!*	1.	PRELIMINARIES
!		-------------
!
!*      1.1     Address module variables common to ICEx and LIMA
!
IF (HCLOUD(1:3) == 'ICE') THEN
  ZALPHAR = XALPHAR_I
  ZNUR    = XNUR_I
  !
  ZAI     = XAI_I
  ZBI     = XBI_I
  ZALPHAI = XALPHAI_I
  ZNUI    = XNUI_I
  !
  ZCCS    = XCCS_I
  ZCXS    = XCXS_I
  ZALPHAS = XALPHAS_I
  ZNUS    = XNUS_I
  !
  ZCCG    = XCCG_I
  ZCXG    = XCXG_I
  ZALPHAG = XALPHAG_I
  ZNUG    = XNUG_I
  !
  ZCCH    = XCCH_I
  ZCXH    = XCXH_I
  ZALPHAH = XALPHAH_I
  ZNUH    = XNUH_I
  !
ELSE IF (HCLOUD == 'LIMA') THEN
  ZALPHAR = XALPHAR_L
  ZNUR    = XNUR_L
  !
  ZAI     = XAI_L
  ZBI     = XBI_L
  ZALPHAI = XALPHAI_L
  ZNUI    = XNUI_L
  !
  ZCCS    = XCCS_L
  ZCXS    = XCXS_L
  ZALPHAS = XALPHAS_L
  ZNUS    = XNUS_L
  !
  ZCCG    = XCCG_L
  ZCXG    = XCXG_L
  ZALPHAG = XALPHAG_L
  ZNUG    = XNUG_L
  !
  ZCCH    = XCCH_L
  ZCXH    = XCXH_L
  ZALPHAH = XALPHAH_L
  ZNUH    = XNUH_L
END IF  
!
!-------------------------------------------------------------------------------
!
!*      2.     SOME CONSTANTS FOR NEUTRALIZATION
!              ---------------------------------
!
IF (HCLOUD(1:3) == 'ICE') THEN
  XFQLIGHTC = 660. * MOMG(3.,3.,2.) / MOMG(3.,3.,3.)   ! PI/A*lbda^(b-2) = 660.
ELSE IF (HCLOUD == 'LIMA') THEN
  XFQLIGHTC = XPI * MOMG(XALPHAC,XNUC,2.)
END IF
!
IF (HCLOUD(1:3) == 'ICE') THEN
  XFQLIGHTR  = XPI * XCCR * MOMG(ZALPHAR,ZNUR,2.)
  XEXQLIGHTR = XCXR - 2.
ELSE IF (HCLOUD == 'LIMA') THEN
  XFQLIGHTR  = XPI * MOMG(ZALPHAR,ZNUR,2.)
  XEXQLIGHTR = -2.
END IF
!
XEXQLIGHTI = 2. / ZBI
XFQLIGHTI  = XPI / 4. * MOMG(ZALPHAI,ZNUI,2.) * &
            (ZAI * MOMG(ZALPHAI,ZNUI,ZBI))**(-XEXQLIGHTI)
!
IF (HCLOUD(1:3) == 'ICE' .OR. &
   (HCLOUD == 'LIMA' .AND. NMOM_S == 1)) THEN
  XFQLIGHTS  = XPI * ZCCS * MOMG(ZALPHAS,ZNUS,2.)
  XEXQLIGHTS = ZCXS - 2.
ELSE IF (HCLOUD == 'LIMA' .AND. NMOM_S == 2) THEN
  XFQLIGHTS  = XPI * MOMG(ZALPHAS,ZNUS,2.)
  XEXQLIGHTS = -2.
END IF
!
IF (HCLOUD(1:3) == 'ICE' .OR. &
   (HCLOUD == 'LIMA' .AND. NMOM_G == 1)) THEN
  XFQLIGHTG  = XPI * ZCCG * MOMG(ZALPHAG,ZNUG,2.)
  XEXQLIGHTG = ZCXG - 2.
ELSE IF (HCLOUD == 'LIMA' .AND. NMOM_G == 2) THEN
  XFQLIGHTG  = XPI * MOMG(ZALPHAG,ZNUG,2.)
  XEXQLIGHTG = -2.
END IF
!
IF (HCLOUD(1:3) == 'ICE' .OR. &
   (HCLOUD == 'LIMA' .AND. NMOM_H == 1)) THEN
  XFQLIGHTH  = XPI * ZCCH * MOMG(ZALPHAH,ZNUH,2.)
  XEXQLIGHTH = ZCXH - 2.
ELSE  IF (HCLOUD == 'LIMA' .AND. NMOM_H == 2) THEN
  XFQLIGHTH  = XPI * MOMG(ZALPHAH,ZNUH,2.)
  XEXQLIGHTH = -2.
END IF
!
IF( .NOT.ALLOCATED(XNEUT_POS)) ALLOCATE( XNEUT_POS(NLGHTMAX) )
IF( .NOT.ALLOCATED(XNEUT_NEG)) ALLOCATE( XNEUT_NEG(NLGHTMAX) )
XNEUT_POS(:) = 0.
XNEUT_NEG(:) = 0.
!
!----------------------------------------------------------------------------
!
!*      3.      INITIALIZE SOME THRESHOLDS
!               --------------------------
!
! electric field threshold for cell detection
! from Marshall et al. (1995) JGR, the breakeven electric field is 
! 200 kV/m at the ground, ~ 33 kV/m at 15 km, and ~ 18 kV/m at 20 km height.
! To be sure all the electrified cells are detected, this threshold is set to 
! 20 kV/m
XE_THRESH =  35.E3 ! (V/m)
!
! the maximum of segments in the bi-leader corresponds to the number of 
! altitude levels in the domain since the bi-leader is hypothesized to 
! propagate only along the vertical 
NLEADER_MAX = NKMAX
!
! the maximum number of branches is arbitriraly set to 5000
NBRANCH_MAX = 5000
!
! the maximum number of electrified cells in the domain is arbitrarily 
! set to 10
NMAX_CELL = 10
!
! the altitude for CG to be prolongated to the ground is set to 2 km
! this threshold could be modified once ions will be taken into account
XALT_CG = 2000.  ! m
!
!
!----------------------------------------------------------------------------
!
!*	3.	INITIALIZATIONS
!		---------------
!
NNBLIGHT   = 0
NNB_CG     = 0
NNB_CG_POS = 0
!
!
!----------------------------------------------------------------------------
!
!*      4.      INITIALIZE LMA RECORDS
!               ----------------------
!
! needs LLMA = .TRUE. to operate
XDTLMA = 600.
TDTLMA = TDTCUR
LWRITE_LMA = .FALSE.
CLMA_FILE(1:5) = "BEGIN"
!
!----------------------------------------------------------------------------
!
END SUBROUTINE INI_FLASH_GEOM_ELEC
