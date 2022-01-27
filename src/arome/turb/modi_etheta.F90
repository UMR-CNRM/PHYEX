!     ######spl
MODULE MODI_ETHETA
!#################
!
INTERFACE
!
FUNCTION ETHETA(KRR,KRRI,PTHLM,PRM,PLOCPEXNM,PATHETA,PSRCM) RESULT(PETHETA)
!
INTEGER                              :: KRR          ! number of moist var.
INTEGER                              :: KRRI         ! number of ice var.
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PTHLM     ! Conservative pot. temperature
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::   PRM       ! Mixing ratios, where
!                                      PRM(:,:,:,1) = conservative mixing ratio
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PLOCPEXNM ! Lv(T)/Cp/Exner at time t-1
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PATHETA   ! Atheta
!                                                    
REAL, DIMENSION(:,:,:),  INTENT(IN)  ::   PSRCM     ! Normalized 2dn_order
                                                    ! moment s'r'c/2Sigma_s2
!
REAL,DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2),SIZE(PTHLM,3)):: PETHETA ! result
!
!
END FUNCTION ETHETA
!
END INTERFACE
!
END MODULE MODI_ETHETA
