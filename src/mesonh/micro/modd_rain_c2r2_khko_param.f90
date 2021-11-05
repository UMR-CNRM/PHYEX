!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/10/16 14:23:23
!-----------------------------------------------------------------
!     ###########################
      MODULE MODD_RAIN_C2R2_KHKO_PARAM
!     ###########################
!
!!****  *MODD_RAIN_C2R2_KHKO_PARAM* - declaration of some microphysical factors
!!                               extensively used in the warm scheme.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare some precomputed
!     microphysical paramters directly used in routine RAIN_C2R2_KHKO
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!          
!!    AUTHOR
!!    ------
!!	J.-P. Pinty  *Laboratoire d'Aerologie*
!!	O.Geoffroy (GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/12/95                      
!!       J.-P. Pinty   29/11/02 add cloud droplet cond/eva parameters for C3R5
!!       G.Delautier 2014 : fusion MODD_RAIN_KHKO_PARAM et MODD_RAIN_C2R2_PARAM
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
IMPLICIT NONE 
!
REAL,SAVE :: XFSEDRR,XFSEDCR,                  & ! Constants for sedimentation
             XFSEDRC,XFSEDCC                     ! fluxes of R, C
!
!
REAL,SAVE :: XDIVA,                            & ! Diffusivity of water vapor
	     XTHCO                               ! Thermal conductivity
REAL,SAVE :: XWMIN                               ! Min value of updraft velocity
				                 ! to enable nucleation process
REAL,SAVE :: XTMIN                               ! Min value of
                                                 ! temperature evolution
				                 ! to enable nucleation process
REAL,SAVE :: XCSTHEN,XCSTDCRIT                   ! Cst for HEN precalculations
INTEGER, SAVE :: NHYP                            ! Number of value of the HYP
						 !    functions
REAL,SAVE :: XHYPINTP1, XHYPINTP2                ! Factors defining the
						 ! supersaturation log scale
REAL, DIMENSION(:), SAVE, ALLOCATABLE          & ! Tabulated HYPgeometric
	  :: XHYPF12, XHYPF32                    !   functions used in HEN
INTEGER, SAVE :: NAHEN                           ! Number of value of the AHEN
		                        	 !    functions
REAL,SAVE :: XAHENINTP1, XAHENINTP2              ! Factors defining the
						 ! temperatures in lin scale
REAL, DIMENSION(:), SAVE, ALLOCATABLE          & ! 
          :: XAHENG,XPSI1, XPSI3,              & ! Twomey-CPB98 and
	     XAHENF,XAHENY                       ! Feingold-Heymsfield
	                                         ! parameterization to compute Smax
REAL,SAVE :: XWCOEF_F1, XWCOEF_F2, XWCOEF_F3,  & ! COEF_F of the polynomial temp.
             XWCOEF_Y1, XWCOEF_Y2, XWCOEF_Y3     ! COEF_Y of the polynomial temp.
						 ! function powering W
!
!
REAL,SAVE :: XKERA1, XKERA2                      ! Constants to define the lin
						 ! and parabolic kernel param. 
REAL,SAVE :: XSELFC                              ! Constants for cloud droplet
                                                 ! selfcollection : SELF
!
REAL,SAVE :: XAUTO1, XAUTO2, XCAUTR,           & ! Constants for cloud droplet
    	     XLAUTR,   XLAUTR_THRESHOLD,       & ! autoconversion : AUT
    	     XITAUTR, XITAUTR_THRESHOLD
!
REAL,SAVE :: XACCR1, XACCR2, XACCR3,           & ! Constants for the accretion
	     XACCR4, XACCR5, XACCR6,           & ! process
             XACCR_CLARGE1, XACCR_CLARGE2, XACCR_RLARGE1, XACCR_RLARGE2, &
             XACCR_CSMALL1, XACCR_CSMALL2, XACCR_RSMALL1, XACCR_RSMALL2
!
REAL,SAVE :: XSCBU2, XSCBU3,                   & ! Constants for the raindrop
             XSCBU_EFF1, XSCBU_EFF2, XSCBUEXP1   ! breakup-selfcollection: SCBU
!
REAL,SAVE :: XSPONBUD1,XSPONBUD2,XSPONBUD3,    & ! Spontaneous Break-up
             XSPONCOEF2                          ! (drop size limiter)
!
REAL,SAVE :: X0EVAR, X1EVAR,                   & ! Constants for raindrop
	     XEX0EVAR, XEX1EVAR, XEX2EVAR        ! evaporation: EVA 
!
REAL,SAVE :: XCONCC_INI, XCONCR_PARAM_INI        ! Used to initialize the 
                                      ! concentrations from mixing ratios
                                      ! (init and grid-nesting from Kessler)
!
REAL,SAVE :: X0CNDC, X2CNDC                   ! Constants for cloud droplet
                                              ! condensation/evaporation
REAL,SAVE :: XFREFFC  ! Factor to compute the cloud droplet effective radius
REAL,SAVE :: XFREFFR  ! Factor to compute the rain drop     effective radius
REAL,SAVE :: XCREC, XCRER
                      ! Factors to compute reff when cloud and rain are present
!
REAL,SAVE ::   XR0                               ! new drizzle drops radius
    	                                         ! autoconversion
!
REAL,SAVE :: XCEVAP                              ! Constants for raindrop
                                                 ! evaporation 

END MODULE MODD_RAIN_C2R2_KHKO_PARAM 
