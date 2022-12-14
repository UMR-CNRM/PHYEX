!MNH_LIC Copyright 2013-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###########################
      MODULE MODD_PARAM_LIMA_WARM
!     ###########################
!
!!****  *MODD_PARAM_LIMA_WARM* - declaration of some descriptive parameters and
!!                               microphysical factors extensively used in 
!!                               the LIMA warm scheme.
!!    AUTHOR
!!    ------
!!  	J.-P. Pinty  *Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!
!-------------------------------------------------------------------------------
USE MODD_PARAMETERS, ONLY: JPSVNAMELGTMAX
!
IMPLICIT NONE 
!
!*       1.   DESCRIPTIVE PARAMETERS
!             ----------------------
!
REAL,SAVE ::  XLBC, XLBEXC,          & ! shape parameters of the cloud droplets
              XLBR, XLBEXR, XNR        ! shape parameters of the raindrops
!
REAL,SAVE :: XAR,XBR,XCR,XDR,XF0R,XF1R,     & ! Raindrop       charact.
                             XCCR,XCXR,     & !For diagnostics
             XAC,XBC,XCC,XDC,XF0C,XF2C,XC1C   ! Cloud droplet  charact.
!
!
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(5),PARAMETER &
                     :: CLIMA_WARM_NAMES=(/'CCLOUD  ','CRAIN   ','CCCNFREE','CCCNACTI','SPRO    '/)
                                       ! basenames of the SV articles stored
                                       ! in the binary files
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(5),PARAMETER &
                     :: CLIMA_WARM_CONC=(/'NC   ','NR   ','NFREE','NCCN ','SS   '/)
!                                       ! basenames of the SV articles stored
!                                       ! in the binary files for DIAG
!
!* Special issue for Below-Cloud SCAVenging of Aerosol particles 
CHARACTER(LEN=JPSVNAMELGTMAX),DIMENSION(2) :: CAERO_MASS =(/'MASSAP', 'MAP   '/)
!
!-------------------------------------------------------------------------------
!
!*       2.   MICROPHYSICAL FACTORS
!             ---------------------
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
REAL, DIMENSION(:,:), SAVE, ALLOCATABLE        & ! Tabulated HYPgeometric
	  :: XHYPF12, XHYPF32                    !   functions used in HEN
INTEGER, SAVE :: NAHEN                           ! Number of value of the AHEN
		                        	 !    functions
REAL,SAVE :: XAHENINTP1, XAHENINTP2              ! Factors defining the
						 ! temperatures in lin scale
REAL, DIMENSION(:), SAVE, ALLOCATABLE          & ! 
          :: XAHENG,XAHENG2,XAHENG3,XPSI1, XPSI3,      & ! Twomey-CPB98 and
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
    	     XLAUTR, XLAUTR_THRESHOLD,         & ! autoconversion : AUT
    	     XITAUTR, XITAUTR_THRESHOLD, XR0     ! XR0 for KHKO autoconversion
!
REAL,SAVE :: XACCR1, XACCR2, XACCR3,           & ! Constants for the accretion
	     XACCR4, XACCR5, XACCR6,           & ! process
             XACCR_CLARGE1, XACCR_CLARGE2, XACCR_RLARGE1, XACCR_RLARGE2, &
             XACCR_CSMALL1, XACCR_CSMALL2, XACCR_RSMALL1, XACCR_RSMALL2, &
             XFCACCR, XEXCACCR
!
REAL,SAVE :: XSCBU2, XSCBU3,                   & ! Constants for the raindrop
             XSCBU_EFF1, XSCBU_EFF2, XSCBUEXP1   ! breakup-selfcollection: SCBU
!
REAL,SAVE :: XSPONBUD1,XSPONBUD2,XSPONBUD3,    & ! Spontaneous Break-up
             XSPONCOEF2                          ! (drop size limiter)
!
REAL,SAVE :: X0EVAR, X1EVAR,                   & ! Constants for raindrop
	     XEX0EVAR, XEX1EVAR, XEX2EVAR,     & ! evaporation: EVA 
             XCEVAP                              ! for KHKO
!
REAL,DIMENSION(:,:,:,:), SAVE, ALLOCATABLE :: XCONCC_INI
REAL,SAVE                                  :: XCONCR_PARAM_INI        
                                      ! Used to initialize the 
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
!-------------------------------------------------------------------------------
!
END MODULE MODD_PARAM_LIMA_WARM
