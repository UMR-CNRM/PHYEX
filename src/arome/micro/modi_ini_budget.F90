!     ######spl
      MODULE MODI_INI_BUDGET
!     ###################### 
INTERFACE
      SUBROUTINE INI_BUDGET(KLUOUT, HLUOUT,PTSTEP,KSV,KRR,            &
      ONUMDIFU,ONUMDIFTH,ONUMDIFSV,                                   &
      OHORELAX_UVWTH,OHORELAX_RV,OHORELAX_RC,OHORELAX_RR,             &
      OHORELAX_RI,OHORELAX_RS, OHORELAX_RG, OHORELAX_RH,OHORELAX_TKE, & 
      OHORELAX_SV,OVE_RELAX,OCHTRANS,ONUDGING,ODRAGTREE,              &
      HRAD,HDCONV,HSCONV,HTURB,HTURBDIM,HCLOUD,                       &
      HMET_ADV_SCHEME,HSV_ADV_SCHEME)
!
INTEGER,         INTENT(IN) :: KLUOUT   ! Logical unit number for prints
CHARACTER (LEN=*), INTENT(IN) :: HLUOUT ! name of output listing
REAL, INTENT(IN) :: PTSTEP              ! time step
INTEGER, INTENT(IN) :: KSV              ! number of scalar variables
INTEGER, INTENT(IN) :: KRR              ! number of moist variables
LOGICAL, INTENT(IN) :: ONUMDIFU         ! switch to activate the numerical
                                        ! diffusion for momentum
LOGICAL, INTENT(IN) :: ONUMDIFTH        ! for meteorological scalar variables
LOGICAL, INTENT(IN) :: ONUMDIFSV        ! for tracer scalar variables
LOGICAL, INTENT(IN) :: OHORELAX_UVWTH  ! switch for the
                       ! horizontal relaxation for U,V,W,TH
LOGICAL, INTENT(IN) :: OHORELAX_RV     ! switch for the
                       ! horizontal relaxation for Rv
LOGICAL, INTENT(IN) :: OHORELAX_RC     ! switch for the
                       ! horizontal relaxation for Rc
LOGICAL, INTENT(IN) :: OHORELAX_RR     ! switch for the
                       ! horizontal relaxation for Rr
LOGICAL, INTENT(IN) :: OHORELAX_RI     ! switch for the
                       ! horizontal relaxation for Ri
LOGICAL, INTENT(IN) :: OHORELAX_RS     ! switch for the
                       ! horizontal relaxation for Rs
LOGICAL, INTENT(IN) :: OHORELAX_RG     ! switch for the
                       ! horizontal relaxation for Rg
LOGICAL, INTENT(IN) :: OHORELAX_RH     ! switch for the
                       ! horizontal relaxation for Rh
LOGICAL, INTENT(IN) :: OHORELAX_TKE    ! switch for the
                       ! horizontal relaxation for tke
LOGICAL,DIMENSION(:),INTENT(IN):: OHORELAX_SV     ! switch for the
                       ! horizontal relaxation for scalar variables
LOGICAL, INTENT(IN) :: OVE_RELAX        ! switch to activate the vertical 
                                        ! relaxation
LOGICAL, INTENT(IN) :: OCHTRANS         ! switch to activate convective 
                                        !transport for SV
LOGICAL, INTENT(IN) :: ONUDGING         ! switch to activate nudging
LOGICAL, INTENT(IN) :: ODRAGTREE        ! switch to activate vegetation drag
CHARACTER (LEN=*), INTENT(IN) :: HRAD   ! type of the radiation scheme
CHARACTER (LEN=*), INTENT(IN) :: HDCONV ! type of the deep convection scheme
CHARACTER (LEN=*), INTENT(IN) :: HSCONV ! type of the deep convection scheme
CHARACTER (LEN=*), INTENT(IN) :: HTURB  ! type of the turbulence scheme
CHARACTER (LEN=*), INTENT(IN) :: HTURBDIM! dimensionnality of the turbulence 
                                        ! scheme
CHARACTER (LEN=*), INTENT(IN) :: HCLOUD ! type of microphysical scheme
CHARACTER (LEN=*), INTENT(IN) :: HMET_ADV_SCHEME ! type of advection scheme    
                                        ! for meteorological scalar variables
CHARACTER (LEN=*), INTENT(IN) :: HSV_ADV_SCHEME ! type of advection scheme    
                                        ! for tracer scalar variables
!
      END SUBROUTINE INI_BUDGET
!
END INTERFACE
!
END MODULE MODI_INI_BUDGET
