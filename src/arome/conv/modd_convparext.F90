!     ######spl
      MODULE MODD_CONVPAREXT
!     ######################
!
IMPLICIT NONE
!
INTEGER, SAVE :: JCVEXB ! start vertical computations at
                        ! 1 + JCVEXB = 1 + ( KBDIA - 1 )
INTEGER, SAVE :: JCVEXT ! limit vertical computations to
                        ! KLEV - JCVEXT = KLEV - ( KTDIA - 1 )
!
END MODULE MODD_CONVPAREXT
