!     ######spl
       MODULE MODI_INI_RAIN_ICE
!      ########################
!
INTERFACE
      SUBROUTINE INI_RAIN_ICE ( KLUOUT, HCLOUD )
!
INTEGER,                 INTENT(IN) :: KLUOUT    ! Logical unit number for prints
!
CHARACTER (LEN=4), INTENT(IN)       :: HCLOUD    ! Indicator of the cloud scheme
!
END SUBROUTINE INI_RAIN_ICE
!
END INTERFACE
!
END MODULE MODI_INI_RAIN_ICE
