!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
      MODULE MODD_CONVPAREXT

!$ACDC methods 

!     ######################
!
IMPLICIT NONE
!
TYPE CONVPAREXT
INTEGER :: JCVEXB ! start vertical computations at
                        ! 1 + JCVEXB = 1 + ( KBDIA - 1 )
INTEGER :: JCVEXT ! limit vertical computations to
                        ! KLEV - JCVEXT = KLEV - ( KTDIA - 1 )
END TYPE CONVPAREXT
!Keep global variables for parts of the code not ported to the type ye
INTEGER :: JCVEXB ! start vertical computations at
                        ! 1 + JCVEXB = 1 + ( KBDIA - 1 )
INTEGER :: JCVEXT ! limit vertical computations to
                        ! KLEV - JCVEXT = KLEV - ( KTDIA - 1 )
!
END MODULE MODD_CONVPAREXT
