!MNH_LIC Copyright 1995-2019 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!####################
MODULE MODI_GAMMA_INC
!####################
!
INTERFACE
!
FUNCTION GAMMA_INC(PA,PX)  RESULT(PGAMMA_INC)
REAL, INTENT(IN)                                  :: PA
REAL, INTENT(IN)                                  :: PX
REAL                                              :: PGAMMA_INC
END FUNCTION GAMMA_INC
!
END INTERFACE
!
END MODULE MODI_GAMMA_INC
