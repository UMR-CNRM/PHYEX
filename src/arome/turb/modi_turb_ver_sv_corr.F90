!     ######spl
     MODULE MODI_TURB_VER_SV_CORR
!    ####################
!
INTERFACE 
!
      SUBROUTINE TURB_VER_SV_CORR(KKA,KKU,KKL,KRR,KRRL,KRRI,                    &
                      PDZZ,                                         &
                      PTHLM,PRM,PTHVREF,                            &
                      PLOCPEXNM,PATHETA,PAMOIST,PSRCM,PPHI3,PPSI3,  &
                      PWM,PSVM,                                     &
                      PTKEM,PLM,PLEPS,PPSI_SV                       )
!
INTEGER,                 INTENT(IN)  ::  KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                 INTENT(IN)  ::  KKL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
INTEGER,                INTENT(IN)   ::  KRRL         ! number of liquid var.
INTEGER,                INTENT(IN)   ::  KRRI         ! number of ice var.
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDZZ 
                                                      ! Metric coefficients
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHLM        ! potential temperature at t-Delta t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM          ! Mixing ratios at t-Delta t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF      ! reference Thv
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM    ! Lv(T)/Cp/Exnref at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM        ! normalized 
                  ! 2nd-order flux s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PPHI3        ! Inv.Turb.Sch.for temperature
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PPSI3        ! Inv.Turb.Sch.for humidity
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PWM          ! w at time t
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM         ! scalar var. at t-Delta t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTKEM        ! TKE at time t
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM          ! Turb. mixing length   
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS        ! dissipative length   
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PPSI_SV      ! Inv.Turb.Sch.for scalars
                            ! cumulated sources for the prognostic variables
!
!
END SUBROUTINE TURB_VER_SV_CORR
!
END INTERFACE
!
END MODULE MODI_TURB_VER_SV_CORR
