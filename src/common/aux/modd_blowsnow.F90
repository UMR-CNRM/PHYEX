!MNH_LIC Copyright 1994-2018 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!!     ######################
       MODULE MODD_BLOWSNOW
!!     ######################
!!
!!     PURPOSE
!!     -------
!!
!!     Declaration of variables and types for the blowing snow scheme
!!
!!     METHOD
!!     ------
!!
!!
!!     REFERENCE
!!     ---------
!!     Etudes du transport de la neige par le vent en conditions alpines : 
!!     Observations et simulations à l'aide d'un modèle couplé atmosphère/
!!     manteau neigeux (Thèse, Uni. Paris Est, 2012)
!!
!!
!!     AUTHOR
!!     ------
!!     Vincent Vionnet (CNRM)
!!
!!
!!     MODIFICATIONS
!!     -------------
!!
!!--------------------------------------------------------------------
!!     DECLARATIONS
!!     ------------
IMPLICIT NONE

LOGICAL      :: LBLOWSNOW  = .FALSE.   ! switch to active pronostic blowing snow
!
INTEGER      :: NBLOWSNOW3D = 2 ! Number of blowing snow variables
! as scalar in Meso-NH. The curent version of the model use two scalars:
!				- Number concentration (#/kg)
!				- Mass concentration (kg/kg)

INTEGER     :: NBLOWSNOW_2D = 3  ! Number of 2D blowing snow variables
! adected in Meso-NH. The curent version of the model advectes three variables:
!             - total number concentration in Canopy
!             - total mass concentration in Canopy
!             - equivalent concentration in the saltation layer
!
REAL            :: XALPHA_SNOW ! Gamma distribution shape factor
!
REAL            :: XRSNOW  ! Ratio between diffusion coefficient for scalar
                           ! variables and blowing snow variables
                           ! RSNOW = KSCA/KSNOW = 4. (if Redelsperger-Sommeria (1981) used in ini_cturb)
                           ! RSNOW = KSCA/KSNOW = 2.5 ( if Cheng-Canuto-Howard (2002) used in ini_cturb)
                           ! Cheng-Canuto-Howard (2002) is the default in MNH V5.3
                           ! See Vionnet (PhD, 2012, In French) and Vionnet et al (TC, 2014)
                           ! for a complete dicsussion
CHARACTER(LEN=6),DIMENSION(:),ALLOCATABLE  :: CSNOWNAMES

CHARACTER(LEN=6),DIMENSION(2), PARAMETER  :: YPSNOW_INI = &
     (/'SNWM01','SNWM02'/)
!
CHARACTER(LEN=6),DIMENSION(3), PARAMETER  :: YPBLOWSNOW_2D = &
     (/'SNWCNU','SNWCMA','SNWCSA' /)

CHARACTER(LEN=4)  :: CSNOWSEDIM ! type of formulation for snow
!              sedimentation : MITC : Mitchell (1996)
!                              CARR : Carrier's drag coefficient (cf PIEKTUK)
!                              TABC : Tabulated values from Carrier's drag coefficient
!                              NONE : no seidmentation
!Minimal mean radius (um) 
REAL           :: XINIRADIUS_SNW = 5.e-6
!Minimum allowed number concentration (#/m3)
REAL           :: XN0MIN_SNW    =  1
!
!
END MODULE MODD_BLOWSNOW
