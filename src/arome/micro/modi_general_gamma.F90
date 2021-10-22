!     ######spl
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
