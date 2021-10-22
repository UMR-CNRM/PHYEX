!     ######spl
        MODULE MODI_GAMMA
!       #################
!
INTERFACE GAMMA
!
FUNCTION GAMMA_X0D(PX)  RESULT(PGAMMA)
REAL, INTENT(IN)                                  :: PX
REAL                                              :: PGAMMA
END FUNCTION GAMMA_X0D
!
FUNCTION GAMMA_X1D(PX)  RESULT(PGAMMA)
REAL, DIMENSION(:), INTENT(IN)        :: PX
REAL, DIMENSION(SIZE(PX))             :: PGAMMA
END FUNCTION GAMMA_X1D
!
END INTERFACE GAMMA
!
END MODULE MODI_GAMMA
