!MNH_LIC Copyright 2006-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODI_COMPUTE_FRAC_ICE  
!    ############################ 
!
INTERFACE COMPUTE_FRAC_ICE
!
    SUBROUTINE COMPUTE_FRAC_ICE3D(HFRAC_ICE,PFRAC_ICE,PT)
!
CHARACTER(len=1),       INTENT(IN) :: HFRAC_ICE
REAL, DIMENSION(:,:,:), INTENT(IN) :: PT
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFRAC_ICE                                             
!
   END SUBROUTINE COMPUTE_FRAC_ICE3D
!
   SUBROUTINE COMPUTE_FRAC_ICE2D(HFRAC_ICE,PFRAC_ICE,PT)
!
CHARACTER(len=1),     INTENT(IN) :: HFRAC_ICE
REAL, DIMENSION(:,:), INTENT(IN) :: PT
REAL, DIMENSION(:,:), INTENT(INOUT) :: PFRAC_ICE                                             
!
   END SUBROUTINE COMPUTE_FRAC_ICE2D

   SUBROUTINE COMPUTE_FRAC_ICE1D(HFRAC_ICE,PFRAC_ICE,PT)
!
CHARACTER(len=1),   INTENT(IN) :: HFRAC_ICE
REAL, DIMENSION(:), INTENT(IN) :: PT
REAL, DIMENSION(:), INTENT(INOUT) :: PFRAC_ICE                                             

   END SUBROUTINE COMPUTE_FRAC_ICE1D

END INTERFACE
!
END MODULE MODI_COMPUTE_FRAC_ICE
!
!     ##############################
      MODULE MODI_COMPUTE_FRAC_ICE3D  
!     ############################## 

INTERFACE 
!
    SUBROUTINE COMPUTE_FRAC_ICE3D(HFRAC_ICE,PFRAC_ICE,PT)
!
CHARACTER(len=1),       INTENT(IN) :: HFRAC_ICE
REAL, DIMENSION(:,:,:), INTENT(IN) :: PT
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFRAC_ICE                                             
!
   END SUBROUTINE COMPUTE_FRAC_ICE3D
END INTERFACE
END MODULE MODI_COMPUTE_FRAC_ICE3D
!
!     ##############################
       MODULE MODI_COMPUTE_FRAC_ICE1D  
!     ############################## 

INTERFACE 
!
    SUBROUTINE COMPUTE_FRAC_ICE1D(HFRAC_ICE,PFRAC_ICE,PT)
!
CHARACTER(len=1),   INTENT(IN)    :: HFRAC_ICE
REAL, DIMENSION(:), INTENT(IN)    :: PT
REAL, DIMENSION(:), INTENT(INOUT) :: PFRAC_ICE
!
   END SUBROUTINE COMPUTE_FRAC_ICE1D
END INTERFACE
END MODULE MODI_COMPUTE_FRAC_ICE1D
!    ##########################################################
      SUBROUTINE COMPUTE_FRAC_ICE3D(HFRAC_ICE,PFRAC_ICE,PT)
!     #################################################################
!
!
!!****  *COMPUTE_FRAC_ICE* - computes ice fraction
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Julien PERGAUD      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         13/03/06
!!      S. Riette        April 2011 optimisation
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODI_COMPUTE_FRAC_ICE1D
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER(len=1),       INTENT(IN)    :: HFRAC_ICE ! scheme to use
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PT        ! Temperature
REAL, DIMENSION(:,:,:), INTENT(INOUT) :: PFRAC_ICE ! Ice fraction (1 for ice only, 0 for liquid only)
!-------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER :: JJ, JK
!-------------------------------------------------------------------------
!
!       0.3  Initialisation
!
!
!----------------------------------------------------------------------------
!
!       1 Compute FRAC_ICE
!         ----------------
!
DO JK=1, SIZE(PT,3)
  DO JJ=1, SIZE(PT,2)
    CALL COMPUTE_FRAC_ICE1D(HFRAC_ICE,PFRAC_ICE(:,JJ,JK),PT(:,JJ,JK))
  ENDDO
ENDDO


END SUBROUTINE COMPUTE_FRAC_ICE3D
!    ##########################################################
      SUBROUTINE COMPUTE_FRAC_ICE2D(HFRAC_ICE,PFRAC_ICE,PT)
!    ##########################################################
!
!
!!****  *COMPUTE_FRAC_ICE* - computes ice fraction
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Julien PERGAUD      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         13/03/06
!!      S. Riette        April 2011 optimisation
!!
!! --------------------------------------------------------------------------
!       0. DECLARATIONS
!          ------------
!
USE MODI_COMPUTE_FRAC_ICE1D
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER(len=1),     INTENT(IN)    :: HFRAC_ICE ! scheme to use
REAL, DIMENSION(:,:), INTENT(IN)    :: PT        ! Temperature
REAL, DIMENSION(:,:), INTENT(INOUT) :: PFRAC_ICE ! Ice fraction (1 for ice only, 0 for liquid only)
!-------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
INTEGER :: JK
!-------------------------------------------------------------------------
!
!       0.3  Initialisation
!
!
!----------------------------------------------------------------------------
!
!       1 Compute FRAC_ICE
!         ----------------
!
DO JK=1, SIZE(PT,2)
  CALL COMPUTE_FRAC_ICE1D(HFRAC_ICE,PFRAC_ICE(:,JK),PT(:,JK))
ENDDO


END SUBROUTINE COMPUTE_FRAC_ICE2D
!    ##########################################################
      SUBROUTINE COMPUTE_FRAC_ICE1D(HFRAC_ICE,PFRAC_ICE,PT)
!    ##########################################################
!
!
!!****  *COMPUTE_FRAC_ICE* - computes ice fraction
!!
!!    PURPOSE
!!    -------
!!
!!**  METHOD
!!    ------
!!    
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      Julien PERGAUD      * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         13/03/06
!!      S. Riette        April 2011 optimisation
!!      S. Riette        08/2016 add option O
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!!
!! --------------------------------------------------------------------------
!       0. DECLARATIONS
!          ------------
!
USE MODD_NEB, ONLY : XTMINMIX, XTMAXMIX
USE MODD_CST, ONLY : XTT
USE MODE_MSG
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
CHARACTER(len=1),   INTENT(IN)    :: HFRAC_ICE  ! scheme to use
REAL, DIMENSION(:), INTENT(IN)    :: PT         ! temperature
REAL, DIMENSION(:), INTENT(INOUT) :: PFRAC_ICE  ! Ice fraction (1 for ice only, 0 for liquid only)
!
!               0.2  declaration of local variables
! 
!
!               0.2  initialisation
!
!
!------------------------------------------------------------------------
!                1. Compute FRAC_ICE
!
IF (HFRAC_ICE=='T') THEN !using Temperature
  PFRAC_ICE(:) = ( XTMAXMIX - PT(:) ) / ( XTMAXMIX - XTMINMIX ) ! freezing interval
ELSEIF (HFRAC_ICE=='O') THEN !using Temperature with old formulae
  PFRAC_ICE(:) = ( XTT - PT(:) ) / 40. ! freezing interval
ELSEIF (HFRAC_ICE=='N') THEN !No ice
  PFRAC_ICE(:) = 0.
ELSEIF (HFRAC_ICE=='S') THEN !Same as previous
  !nothing to do
ELSE
  call Print_msg(NVERB_FATAL,'GEN','COMPUTE_FRAC_ICE','invalid option for HFRAC_ICE='//HFRAC_ICE)
ENDIF

PFRAC_ICE(:) = MAX( 0., MIN(1., PFRAC_ICE(:) ) )


END SUBROUTINE COMPUTE_FRAC_ICE1D
