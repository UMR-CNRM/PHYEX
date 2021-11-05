!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 modd 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     ######################
      MODULE MODD_PARAM_C1R3
!     ######################
!
!!****  *MODD_PARAM_C1R3* - declaration of the control parameters
!!                               for use in the cold scheme.
!!
!!    PURPOSE
!!    -------
!       The purpose of this declarative module is to declare the microphysical
!     constants. This includes the descriptive parameters for the raindrop 
!     and the parameters relevant of the dimensional distributions.
!
!
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!      None 
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (MODD_PARAM_C1R3)
!!          
!!    AUTHOR
!!    ------
!!	J.-P. Pinty  *Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/04/2001
!!      Jean-Pierre PINTY    29/ 6/01  Add RHHONI process (freezing haze part.)
!!
!-------------------------------------------------------------------------------
!
!*       0.   DECLARATIONS
!             ------------
!
REAL,SAVE :: XALPHAI,XNUI,            & ! Pristine ice   distribution parameters
	     XALPHAS,XNUS,            & ! Snow/aggregate distribution parameters
	     XALPHAG,XNUG               ! Graupel        distribution parameters
REAL,SAVE :: XFACTNUC_DEP,XFACTNUC_CON  ! Amplification factor for IN conc.
                                        !  DEP refers to DEPosition mode
                                        !  CON refers to CONtact    mode
!
LOGICAL, SAVE :: LSEDI                  ! TRUE to enable the pristine ice 
                                        !                sedimentation
LOGICAL, SAVE :: LHHONI                 ! TRUE to enable the freezing of haze
                                        !                particules
!
CHARACTER(LEN=4), SAVE :: CPRISTINE_ICE_C1R3 ! Pristine type PLAT, COLU or BURO
CHARACTER(LEN=4), SAVE :: CHEVRIMED_ICE_C1R3 ! Heavily rimed type GRAU or HAIL
!
END MODULE MODD_PARAM_C1R3
!
!
