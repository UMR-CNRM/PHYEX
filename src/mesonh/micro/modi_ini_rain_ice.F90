!     ######spl
       MODULE MODI_INI_RAIN_ICE
!      ########################
!
INTERFACE
      SUBROUTINE INI_RAIN_ICE ( KLUOUT, PTSTEP, PDZMIN, KSPLITR, HCLOUD )
!
INTEGER,                 INTENT(IN) :: KLUOUT    ! Logical unit number for prints
INTEGER,                 INTENT(OUT):: KSPLITR   ! Number of small time step
                                                 ! integration for  rain
                                                 ! sedimendation
!
REAL,                    INTENT(IN) :: PTSTEP    ! Effective Time step
!
REAL,                    INTENT(IN) :: PDZMIN    ! minimun vertical mesh size
!
CHARACTER (LEN=4), INTENT(IN)       :: HCLOUD    ! Indicator of the cloud scheme
!
END SUBROUTINE INI_RAIN_ICE
!
END INTERFACE
!
END MODULE MODI_INI_RAIN_ICE
