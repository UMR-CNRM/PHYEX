!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 operators 2006/05/18 13:07:25
!-----------------------------------------------------------------
!########################
MODULE MODI_GENERAL_GAMMA
!########################
!
INTERFACE
!
FUNCTION GENERAL_GAMMA(PALPHA,PNU,PLBDA,PX)  RESULT(PGENERAL_GAMMA)
REAL, INTENT(IN)                                  :: PALPHA
REAL, INTENT(IN)                                  :: PNU
REAL, INTENT(IN)                                  :: PLBDA
REAL, INTENT(IN)                                  :: PX
REAL                                              :: PGENERAL_GAMMA
END FUNCTION GENERAL_GAMMA
!
END INTERFACE
!
END MODULE MODI_GENERAL_GAMMA
!     ###################################################################
      FUNCTION GENERAL_GAMMA(PALPHA,PNU,PLBDA,PX)  RESULT(PGENERAL_GAMMA)
!     ###################################################################
!     
!
!!****  *GENERAL_GAMMA * -  Generalized gamma  function  
!!                   
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the Generalized gamma
!    function of its argument.
!    
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODULE MODI_GAMMA : computation of the Gamma function
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine CONDENS)
!!
!!
!!    AUTHOR
!!    ------
!!  	Jean-Pierre Pinty *LA/OMP*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     7/11/95
!
!*       0. DECLARATIONS
!           ------------
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                     :: PALPHA
REAL, INTENT(IN)                     :: PNU
REAL, INTENT(IN)                     :: PLBDA
REAL, INTENT(IN)                     :: PX
REAL                                 :: PGENERAL_GAMMA
!
!*       0.2 declarations of local variables
!
REAL                                 :: ZARG,ZPOWER
!
ZARG   = PLBDA*PX
ZPOWER = PALPHA*PNU - 1.0
!
PGENERAL_GAMMA = (PALPHA/GAMMA(PNU))*(ZARG**ZPOWER)*PLBDA*EXP(-(ZARG**PALPHA))
RETURN
!
END FUNCTION GENERAL_GAMMA
