!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!      #######################
          MODULE MODD_CTURB
!      #######################
!
!!****   *MODD_CTURB*  - declaration of the turbulent scheme constants
!!
!!     PURPOSE
!!     -------
!        The purpose of this declarative module is to declare the 
!      turbulence scheme constants.
!
!!
!!**   IMPLICIT ARGUMENTS
!!     ------------------
!!       NONE
!!
!!     REFERENCE
!!     ---------
!!       Book 2 of Meso-NH documentation (MODD_CTURB)
!!       Book 1 of Meso-NH documentation (Chapter Turbulence)
!!
!!     AUTHOR
!!     ------
!1       Joan Cuxart         * INM  and Meteo-France *
!!
!!     MODIFICATIONS
!!     -------------
!!       Original            08/08/94
!!     Nov 06, 2002 (V. Masson)  add XALPSBL and XASBL
!!     May 06                    Remove EPS
!!     Jan 2019     (Q. Rodier)  Remove XASBL
!!     Jan 2022     (Q. Rodier)  introduction of a strucuture
!----------------------------------------------------------------------------
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE
TYPE CSTURB_t
!
REAL :: XCMFS        ! constant for the momentum flux due to shear   
REAL :: XCMFB        ! constant for the momentum flux due to buoyancy
REAL :: XCSHF        ! constant for the sensible heat flux 
REAL :: XCHF         ! constant for the humidity flux 
REAL :: XCTV         ! constant for the temperature variance
REAL :: XCHV         ! constant for the humidity variance
REAL :: XCHT1        ! first ct. for the humidity-temperature correlation
REAL :: XCHT2        ! second ct. for the humidity-temperature correlation
!
REAL :: XCPR1        ! first ct. for the turbulent Prandtl numbers
REAL :: XCPR2        ! second ct. for the turbulent Prandtl numbers
REAL :: XCPR3        ! third ct. for the turbulent Prandtl numbers
REAL :: XCPR4        ! fourth ct. for the turbulent Prandtl numbers
REAL :: XCPR5        ! fifth ct. for the turbulent Prandtl numbers
!
REAL :: XCET         ! constant into the transport term of the TKE eq.
REAL :: XCED         ! constant into the dissipation term of the TKE eq.
!
REAL :: XCDP         ! ct. for the production term in the dissipation eq.
REAL :: XCDD         ! ct. for the destruction term in the dissipation eq.
REAL :: XCDT         ! ct. for the transport term in the dissipation eq.
!
REAL :: XTKEMIN      ! mimimum value for the TKE
REAL :: XRM17        ! Rodier et al 2017 constant in shear term for mixing length
!
REAL :: XLINI        ! initial value for BL mixing length
REAL :: XLINF        ! to prevent division by zero in the BL algorithm
!
REAL :: XALPSBL      ! constant linking TKE and friction velocity in the SBL
!
REAL :: XCEP         ! Constant for wind pressure-correlations
REAL :: XA0          ! Constant a0 for wind pressure-correlations
REAL :: XA2          ! Constant a2 for wind pressure-correlations
REAL :: XA3          ! Constant a3 for wind pressure-correlations
REAL :: XA5          ! Constant a5 for temperature pressure-correlations
REAL :: XCTD         ! Constant for temperature and vapor dissipation
REAL :: XCTP         ! Constant for temperature and vapor pressure-correlations
!
REAL :: XPHI_LIM     ! Threshold value for Phi3 and Psi3
REAL :: XSBL_O_BL    ! SBL height / BL height ratio
REAL :: XFTOP_O_FSURF! Fraction of surface (heat or momentum) flux used to define top of BL
!
END TYPE CSTURB_t
!
TYPE(CSTURB_t), TARGET, SAVE :: CSTURB
!
REAL,POINTER :: XCMFS=>CSTURB%XCMFS
REAL,POINTER :: XCMFB=>CSTURB%XCMFB
REAL,POINTER :: XCSHF=>CSTURB%XCSHF
REAL,POINTER :: XCHF=>CSTURB%XCHF
REAL,POINTER :: XCTV=>CSTURB%XCTV
REAL,POINTER :: XCHV=>CSTURB%XCHV
REAL,POINTER :: XCHT1=>CSTURB%XCHT1
REAL,POINTER :: XCHT2=>CSTURB%XCHT2
!
REAL,POINTER :: XCPR1=>CSTURB%XCPR1
REAL,POINTER :: XCPR2=>CSTURB%XCPR2
REAL,POINTER :: XCPR3=>CSTURB%XCPR3
REAL,POINTER :: XCPR4=>CSTURB%XCPR4
REAL,POINTER :: XCPR5=>CSTURB%XCPR5
!
REAL,POINTER :: XCET=>CSTURB%XCET
REAL,POINTER :: XCED=>CSTURB%XCED
!
REAL,POINTER :: XCDP=>CSTURB%XCDP
REAL,POINTER :: XCDD=>CSTURB%XCDD
REAL,POINTER :: XCDT=>CSTURB%XCDT
!
REAL,POINTER :: XTKEMIN=>CSTURB%XTKEMIN
REAL,POINTER :: XRM17=>CSTURB%XRM17
!
REAL,POINTER :: XLINI=>CSTURB%XLINI
REAL,POINTER :: XLINF=>CSTURB%XLINF
!
REAL,POINTER :: XALPSBL=>CSTURB%XALPSBL
!
REAL,POINTER :: XCEP=>CSTURB%XCEP
REAL,POINTER :: XA0=>CSTURB%XA0
REAL,POINTER :: XA2=>CSTURB%XA2
REAL,POINTER :: XA3=>CSTURB%XA3
REAL,POINTER :: XA5=>CSTURB%XA5
REAL,POINTER :: XCTD=>CSTURB%XCTD
REAL,POINTER :: XCTP=>CSTURB%XCTP
!
REAL,POINTER :: XPHI_LIM=>CSTURB%XPHI_LIM
REAL,POINTER :: XSBL_O_BL=>CSTURB%XSBL_O_BL
REAL,POINTER :: XFTOP_O_FSURF=>CSTURB%XFTOP_O_FSURF
END MODULE MODD_CTURB
