!     ######spl
     MODULE MODI_PRANDTL
!    ###################
!
INTERFACE
!
      SUBROUTINE PRANDTL(KKA,KKU,KKL,KRR,KRRI,OCLOSE_OUT,OTURB_DIAG,&
                         HTURBDIM,                             &
                         HFMFILE,HLUOUT,                       &
                         PDXX,PDYY,PDZZ,PDZX,PDZY,             &
                         PTHVREF,PLOCPEXNM,PATHETA,PAMOIST,    &
                         PLM,PLEPS,PTKEM,PTHLM,PRM,PSVM,PSRCM, &
                         PREDTH1,PREDR1,                       &
                         PRED2TH3, PRED2R3, PRED2THR3,         &
                         PREDS1,PRED2THS3, PRED2RS3,           &
                         PBLL_O_E,                             &
                         PETHETA, PEMOIST                      )
!
!
INTEGER,                INTENT(IN)   :: KKA           !near ground array index  
INTEGER,                INTENT(IN)   :: KKU           !uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL           !vert. levels type 1=MNH -1=ARO
INTEGER,                INTENT(IN)   :: KRR           ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRI          ! number of ice var.
!
LOGICAL,                INTENT(IN)   ::  OCLOSE_OUT   ! switch for syncronous
                                                      ! file opening
LOGICAL,                INTENT(IN)   ::  OTURB_DIAG   ! switch to write some
                                 ! diagnostic fields in the syncronous FM-file
CHARACTER*4           , INTENT(IN)   ::  HTURBDIM     ! Kind of turbulence param.
CHARACTER(LEN=*),       INTENT(IN)   ::  HFMFILE      ! Name of the output
                                                      ! FM-file
CHARACTER(LEN=*),       INTENT(IN)   ::  HLUOUT       ! Output-listing name for
                                                      ! model n
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PDXX,PDYY,PDZZ,PDZX,PDZY
                                                  ! metric coefficients
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHVREF  ! Virtual Potential Temp.
                                                  ! of the reference state
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLOCPEXNM ! Lv(T)/Cp/Exner at t-1 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PATHETA      ! coefficients between 
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PAMOIST      ! s and Thetal and Rnp
!
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLM      ! Turbulent Mixing length
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PLEPS    ! Dissipative length
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PTHLM,PTKEM! Conservative Potential 
                                                  ! Temperature and TKE at t-1
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PRM      ! Mixing ratios at  t-1
                                                  ! with PRM(:,:,:,1) = cons.
                                                  ! mixing ratio
REAL, DIMENSION(:,:,:,:), INTENT(IN) ::  PSVM     ! Scalars at t-1      
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PSRCM
                                  ! s'r'c/2Sigma_s2 at t-1 multiplied by Lambda_3
!
!
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PREDTH1 ! Redelsperger number R_theta
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PREDR1  ! Redelsperger number R_q
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PRED2TH3 ! Redelsperger number R*2_theta
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PRED2R3  ! Redelsperger number R*2_q
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PRED2THR3! Redelsperger number R*2_thq
REAL, DIMENSION(:,:,:,:), INTENT(OUT)::  PREDS1   ! Redelsperger number R_sv
REAL, DIMENSION(:,:,:,:), INTENT(OUT)::  PRED2THS3! Redelsperger number R*2_thsv
REAL, DIMENSION(:,:,:,:), INTENT(OUT)::  PRED2RS3 ! Redelsperger number R*2_qsv
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PBLL_O_E! beta*Lk*Leps/tke
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PETHETA ! coefficient E_theta
REAL, DIMENSION(:,:,:), INTENT(OUT)  ::  PEMOIST ! coefficient E_moist
!
END SUBROUTINE PRANDTL
!
END INTERFACE
!
END MODULE MODI_PRANDTL
