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
      SUBROUTINE INI_FLASH_GEOM_ELEC
!
END SUBROUTINE INI_FLASH_GEOM_ELEC
END INTERFACE
END MODULE MODI_INI_FLASH_GEOM_ELEC
!
!	##############################
        SUBROUTINE INI_FLASH_GEOM_ELEC
!	##############################
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
!!
!-------------------------------------------------------------------------------
!
!*	0.	DECLARATIONS
!		------------
!
USE MODD_CST, ONLY : XPI
USE MODD_RAIN_ICE_DESCR
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
!
!*	0.2	Declaration of local variables
!
!
!----------------------------------------------------------------------------
!
!*      1.     SOME CONSTANTS FOR NEUTRALIZATION
!              ---------------------------------
!
XFQLIGHTC  = 660. * MOMG(3.,3.,2.) / MOMG(3.,3.,3.)   ! PI/A*lbda^(b-2) = 660.
!
XFQLIGHTR  = XPI * XCCR * MOMG(XALPHAR,XNUR,2.)
XEXQLIGHTR = XCXR - 2.
!
XEXQLIGHTI = 2. / XBI
XFQLIGHTI  = XPI / 4. * MOMG(XALPHAI,XNUI,2.) *                   &
            (XAI * MOMG(XALPHAI,XNUI,XBI))**(-XEXQLIGHTI)
!
XFQLIGHTS  = XPI * XCCS * MOMG(XALPHAS,XNUS,2.)
XEXQLIGHTS = XCXS - 2.
!
XFQLIGHTG  = XPI * XCCG * MOMG(XALPHAG,XNUG,2.)
XEXQLIGHTG = XCXG - 2.
!
!
!----------------------------------------------------------------------------
!
!*      2.      INITIALIZE SOME THRESHOLDS
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
