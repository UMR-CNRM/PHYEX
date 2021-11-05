!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$ $Date$
!-----------------------------------------------------------------
!-----------------------------------------------------------------
!       #####################
        MODULE MODI_INI_CTURB
!       #####################
!
INTERFACE
!
SUBROUTINE INI_CTURB
END SUBROUTINE INI_CTURB
!
END INTERFACE
!
END MODULE MODI_INI_CTURB
!
!
!
!       ####################
        SUBROUTINE INI_CTURB
!       ####################
!
!!****     *INI_CTURB*  - routine to initialize the turbulence scheme 
!!                        constants.
!!
!!      PURPOSE
!!      -------
!         The purpose of this routine is to initialize the turbulence 
!       scheme constants that are stored in module MODD_CTURB
!
!!      METHOD
!!      ------
!!        The constants are set to their numerical values
!!
!!      EXTERNAL
!!      --------
!!        NONE
!!
!!      IMPLICIT ARGUMENTS
!!      ------------------
!!        Module MODD_CTURB
!!
!!      REFERENCE
!!      ---------
!!        Book 2 of Meso-NH documentation (module INI_CTURB)
!!        Book 1 of Meso-NH documentation (Chapter Turbulence)
!!
!!      AUTHOR
!!      ------
!!        Joan Cuxart       * INM and Meteo-France *
!!
!!      MODIFICATIONS
!!      -------------
!!        Original          08/08/94
!!        J.Cuxart          15/06/95   document more precisely the Shuman cts
!!        P.Jabouille       20/10/99   XCET=0.4
!!        V.Masson          13/11/02   XALPSBL and XASBL
!!                             05/06   Remove KEPS
!!        Q.Rodier             01/19   XCED replaced by XCEDIS in read_exsegn.f90 and ini_modeln.f90
!!                                     Remove XASBL (not used)
!! --------------------------------------------------------------------------
!
!*        0. DECLARATIONS
!            ------------
!
USE MODD_CST
USE MODD_CTURB
!
IMPLICIT NONE
!
!  ---------------------------------------------------------------------------
!
!         1. SETTING THE NUMERICAL VALUES
!            ----------------------------
!
!         1.1 Constant for dissipation of Tke
!
!XCED is now replaced by XCEDIS
!XCED  = 0.70
!XCED  = 0.84
!
!       Redelsperger-Sommeria (1981) = 0.70
!       Schmidt-Schumann      (1989) = 0.845
!       Cheng-Canuto-Howard   (2002) = 0.845
!       Rodier, Masson, Couvreux, Paci (2017) = 0.34
!
!
!         1.2 Constant for wind pressure-correlations
!
!XCEP  = 4.
XCEP  = 2.11 
!       Redelsperger-Sommeria (1981) = 4.
!       Schmidt-Schumann      (1989) = 3.5
!       Cheng-Canuto-Howard   (2002) = 2.11
!
!
!         1.3 Constant a0 for wind pressure-correlations
!
XA0   = 0.6
!       Redelsperger-Sommeria (1981) = 0.6
!       Schmidt-Schumann      (1989) = 0.55
!       Cheng-Canuto-Howard   (2002) = 0.6
!
!
!         1.4 Constant a2 for wind pressure-correlations
!
XA2   = 1.
!       Redelsperger-Sommeria (1981) = 1.
!       Schmidt-Schumann      (1989) = 1.
!       Cheng-Canuto-Howard   (2002) = 0.57
!
!
!         1.5 Constant a3 for wind pressure-correlations
!
XA3   = 0.
!       Redelsperger-Sommeria (1981) = 0.
!       Schmidt-Schumann      (1989) = 0.45
!       Cheng-Canuto-Howard   (2002) = 0.5
!
!
!         1.6 Constant for dissipation of th'2, r'2, th'r'
!
XCTD  = 1.2
!       Redelsperger-Sommeria (1981) = 1.2
!       Schmidt-Schumann      (1989) = 1.01
!       Cheng-Canuto-Howard   (2002) = 0.98
!
!
!         1.7 Constant for temperature and vapor pressure-correlations
!
XCTP  = 4.65
!       Redelsperger-Sommeria (1981) = 4.
!       Schmidt-Schumann      (1989) = 3.25
!       Cheng-Canuto-Howard   (2002) = 4.65
!
!
!         1.8 Constant a5 for temperature pressure-correlations
!
XA5   = 1./3.
!       Redelsperger-Sommeria (1981) = 1./3.
!       Schmidt-Schumann      (1989) = 0.
!       Cheng-Canuto-Howard   (2002) = 1./3.
!
!
!         1.9 Values in the evolution equation of the TKE
!
XCET  = 0.40
!
!       Redelsperger-Sommeria (1981) = 0.20
!       Schmidt-Schumann      (1989) = 0.33
!       Krettenauer-Schumann  (1992) = 0.33
!       Bougeault and Lacarrere(1989)= 0.40
!
!
!         1.10  Value related to the TKE universal function within SBL
!
XALPSBL = 4.63
!       Redelsperger et al 2001     = 4.63
!       Wyngaard et al. 1974        = 3.75
!       Stull 1988                  = 4.75
!
!
!         1.11  Value related to the shear term in mixing length computation
!
XRM17 = 0.5  ! Rodier et al 2017
!
!
!         2. Derivated constants
!            -------------------
!
!         2.1 Constant in fluxes equations
!
XCMFS= 2./3./XCEP*(1.-XA0)   !Constant for the momentum flux due to shear (RS)
!
! Redelsperger-Sommeria (1981) ......... 0.066
! Schmidt-Schumann      (1989) ......... 0.086
!
!
XCSHF= 2./3./XCTP            !Constant for the sensible heat flux(RS)
!
! Redelsperger-Sommeria (1981) ......... 0.167
! Schmidt-Schumann      (1989) ......... 0.204
!
!
XCHF= XCSHF                  !Constant for the humidity flux(RS)
!
!         2.2 Constant in variances and covariances equations
!
XCTV= 2./3./XCTP/XCTD        !Constant for the temperature variance(RS)
!
! Redelsperger-Sommeria (1981) ......... 0.139
! Schmidt-Schumann      (1989) ......... 0.202
!
XCHV=  XCTV                  !Constant for the humidity variance(RS)
!
! Redelsperger-Sommeria (1981) ......... 0.139
!
!
XCHT1= XCTV/2.      !Constants for the temperature-humidity correlation(RS)
XCHT2= XCTV/2.
!
!         2.3 Constant in Prandtl numbers
!
XCPR1= XCTV         !Constants for the turbulent Prandtl and Schmidt numbers
XCPR2= XCHT1
XCPR3= XCPR2        ! used only for the Schmidt number for scalar variables
XCPR4= XCPR2
XCPR5= XCPR2
!
!         3. MINIMUM VALUES 
!            --------------
!
XTKEMIN=0.01
!
!XLINI=10.   ! BL mixing length
XLINI=0.1   ! BL mixing length
XLINF=1.E-10! to prevent division by zero
!
!
!         4. MAXIMUM VALUES 
!            --------------
!
XPHI_LIM = 3.
!
!
!         5. Constants in K-eps scheme
!            -------------------------
!
!         1.3 Values in the evolution equation of the dissipation of TKE
XCDP  =  1.46
!       Duynkerke (1988)             = 1.46
!
XCDD  =  1.83
!       Duynkerke (1988)             = 1.83
!
XCDT  =  0.42
!       Duynkerke (1988)             = 1./(2.38)
!
!
!         6. Constants in RMC01
!            ------------------
!
XSBL_O_BL     = 0.05 ! SBL height / BL height ratio
XFTOP_O_FSURF = 0.05 ! Fraction of surface (heat or momentum) flux used to define top of BL
!
!
END SUBROUTINE INI_CTURB
