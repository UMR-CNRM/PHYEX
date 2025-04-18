!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!########################
        MODULE MODI_GAMMA
!       #################
!
IMPLICIT NONE
INTERFACE GAMMA
!
PURE FUNCTION GAMMA_X0D(PX)  RESULT(PGAMMA)
IMPLICIT NONE
!$acc routine seq
REAL, INTENT(IN)                                  :: PX
REAL                                              :: PGAMMA
END FUNCTION GAMMA_X0D
!
PURE FUNCTION GAMMA_X1D(PX)  RESULT(PGAMMA)
IMPLICIT NONE
REAL, DIMENSION(:), INTENT(IN)        :: PX
REAL, DIMENSION(SIZE(PX))             :: PGAMMA
END FUNCTION GAMMA_X1D
!
END INTERFACE GAMMA
!
END MODULE MODI_GAMMA
