!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ################
      MODULE MODI_MOMG
!     ################
!
INTERFACE MOMG
!
FUNCTION MOMG_X0D(PALPHA, PNU, PP)  RESULT(PMOMG)
REAL, INTENT(IN) :: PALPHA, PNU
REAL, INTENT(IN) :: PP
REAL             :: PMOMG
END FUNCTION MOMG_X0D
!
FUNCTION MOMG_X1D(PALPHA, PNU, PP)  RESULT(PMOMG)
REAL,               INTENT(IN) :: PALPHA, PNU
REAL, DIMENSION(:), INTENT(IN) :: PP
REAL, DIMENSION(SIZE(PP))      :: PMOMG
END FUNCTION MOMG_X1D
!
END INTERFACE
END MODULE MODI_MOMG
!
!--------------------------------------------------------------------------
!
!
!!****  *MOMG* -   
!!
!!    PURPOSE
!!    -------
!!      Compute: G(p) = Gamma(nu + p/alpha) / Gamma(nu) 
!!                    = M(p) * lambda^p
!!
!!**  IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      C. Barthe    * Laboratoire de l'Atmosphere et des Cyclones *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original   26 Nov. 2009
!!
!--------------------------------------------------------------------------------
!
!*       1.   FUNCTION MOMG FOR SCALAR VARIABLE
!             ---------------------------------
!
!     ##############################################
      FUNCTION MOMG_X0D(PALPHA, PNU, PP)  RESULT(PMOMG)
!     ##############################################
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
REAL, INTENT(IN) :: PALPHA, PNU
REAL, INTENT(IN) :: PP
REAL             :: PMOMG
!
!
PMOMG = GAMMA(PNU+PP/PALPHA) / GAMMA(PNU)
RETURN
!
END FUNCTION MOMG_X0D
!
!-------------------------------------------------------------------------------
!
!*       2.   FUNCTION MOMG FOR 1D ARRAY
!             --------------------------
!
!     ##############################################
      FUNCTION MOMG_X1D(PALPHA, PNU, PP)  RESULT(PMOMG)
!     ##############################################
!
USE MODI_GAMMA
!
IMPLICIT NONE
!
REAL,               INTENT(IN) :: PALPHA, PNU
REAL, DIMENSION(:), INTENT(IN) :: PP
REAL, DIMENSION(SIZE(PP))      :: PMOMG
!
!
PMOMG(:) = GAMMA(PNU+PP(:)/PALPHA) / GAMMA(PNU)
RETURN
!
END FUNCTION MOMG_X1D
!
!------------------------------------------------------------------------------
