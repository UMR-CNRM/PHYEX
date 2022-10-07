!     ######spl
      MODULE MODI_CONDENSATION
!     ########################
!
INTERFACE
!
    SUBROUTINE CONDENSATION(D, CST, ICEP, NEB, TURBN, &
                           &HFRAC_ICE, HCONDENS, HLAMBDA3,                                                  &
                           &PPABS, PZZ, PRHODREF, PT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT,    &
                           &PRR, PRS, PRG, PSIGS, LMFCONV, PMFCONV, PCLDFR, PSIGRC, OUSERI,                 &
                           &OSIGMAS, OCND2, LHGT_QS,                                                        &
                           &PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR, PSIGQSAT,                                 &
                           &PLV, PLS, PCPH,                                                                 &
                           &PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,                                         &
                           &PICE_CLD_WGT)
!
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_CST,        ONLY: CST_t
USE MODD_NEB,        ONLY: NEB_t
USE MODD_TURB_n,     ONLY: TURB_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
!
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)    :: ICEP
TYPE(NEB_t),                  INTENT(IN)    :: NEB
TYPE(TURB_t),                 INTENT(IN)    :: TURBN
CHARACTER(LEN=1),             INTENT(IN)    :: HFRAC_ICE
CHARACTER(LEN=4),             INTENT(IN)    :: HCONDENS
CHARACTER(LEN=*),             INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODREF
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRV_IN ! grid scale water vapor mixing ratio (kg/kg) in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PRV_OUT! grid scale water vapor mixing ratio (kg/kg) in output
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRC_IN ! grid scale r_c mixing ratio (kg/kg) in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PRC_OUT! grid scale r_c mixing ratio (kg/kg) in output
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRI_IN ! grid scale r_i (kg/kg) in input
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PRI_OUT! grid scale r_i (kg/kg) in output
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRR    ! grid scale mixing ration of rain (kg/kg)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRS    ! grid scale mixing ration of snow (kg/kg)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRG    ! grid scale mixing ration of graupel (kg/kg)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PSIGS  ! Sigma_s from turbulence scheme
LOGICAL,                                                       INTENT(IN)    ::  LMFCONV ! =SIZE(PMFCONV)!=0
REAL, DIMENSION(MERGE(D%NIJT,0,LMFCONV),&
                MERGE(D%NKT,0,LMFCONV)),              INTENT(IN)    :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PCLDFR ! cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PSIGRC ! s r_c / sig_s^2

LOGICAL, INTENT(IN)                         :: OUSERI ! logical switch to compute both
                                                      ! liquid and solid condensate (OUSERI=.TRUE.)
                                                      ! or only solid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                         :: OSIGMAS! use present global Sigma_s values
                                                      ! or that from turbulence scheme
LOGICAL, INTENT(IN)                         :: OCND2  ! logical switch to sparate liquid and ice
                                                      ! more rigid (DEFALT value : .FALSE.)
LOGICAL, INTENT(IN)                         :: LHGT_QS! logical switch for height dependent VQSIGSAT
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PICLDFR  ! ice cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PWCLDFR  ! water or mixed-phase cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PSSIO    ! Super-saturation with respect to ice in the  
                                                              ! supersaturated fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PSSIU    ! Sub-saturation with respect to ice in the  
                                                              ! subsaturated fraction
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PIFR     ! Ratio cloud ice moist part
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)
                                                              ! multiplied by PSIGQSAT

REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PLV    ! Latent heat L_v
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PLS    ! Latent heat L_s
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(IN)    :: PCPH   ! Specific heat C_ph
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLC_HRC
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLC_HCF ! cloud fraction
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLI_HRI
REAL, DIMENSION(D%NIJT,D%NKT), OPTIONAL, INTENT(OUT)   :: PHLI_HCF
REAL, DIMENSION(D%NIJT),       OPTIONAL, INTENT(IN)    :: PICE_CLD_WGT
END SUBROUTINE CONDENSATION
!
END INTERFACE
!
END MODULE MODI_CONDENSATION
