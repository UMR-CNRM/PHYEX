!     ######spl
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
