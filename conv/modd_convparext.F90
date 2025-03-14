!     ######spl
      MODULE MODD_CONVPAREXT
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
