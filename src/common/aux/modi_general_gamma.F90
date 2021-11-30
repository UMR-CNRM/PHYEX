!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
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
