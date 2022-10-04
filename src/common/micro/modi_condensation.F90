!     ######spl
      MODULE MODI_CONDENSATION
!     ########################
!
INTERFACE
!
       SUBROUTINE CONDENSATION( KIU, KJU, KKU, KIB, KIE, KJB, KJE, KKB, KKE, KKL,&
          HFRAC_ICE, HCONDENS, HLAMBDA3, &
          PPABS, PZZ, PRHODREF, PT, PRV_IN, PRV_OUT, PRC_IN, PRC_OUT, PRI_IN, PRI_OUT, &
          PRR, PRS, PRG, PSIGS, PMFCONV, PCLDFR, &
          PSIGRC, OUSERI, OSIGMAS, OCND2, LHGT_QS, &
          PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR, PSIGQSAT, &
          PLV, PLS, PCPH, &
          PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF, &
          PICE_CLD_WGT)
!
!
INTEGER,                      INTENT(IN)    :: KIU    ! horizontal dimension in x
INTEGER,                      INTENT(IN)    :: KJU    ! horizontal dimension in y
INTEGER,                      INTENT(IN)    :: KKU    ! vertical dimension
INTEGER,                      INTENT(IN)    :: KIB    ! value of the first point in x
INTEGER,                      INTENT(IN)    :: KIE    ! value of the last  point in x
INTEGER,                      INTENT(IN)    :: KJB    ! value of the first point in y
INTEGER,                      INTENT(IN)    :: KJE    ! value of the last  point in y
INTEGER,                      INTENT(IN)    :: KKB    ! value of the first point in z
INTEGER,                      INTENT(IN)    :: KKE    ! value of the last  point in z
INTEGER,                      INTENT(IN)    :: KKL    ! +1 if grid goes from ground to atmosphere top, -1 otherwise
CHARACTER(LEN=1),             INTENT(IN)    :: HFRAC_ICE
CHARACTER(LEN=4),             INTENT(IN)    :: HCONDENS
CHARACTER(LEN=*),             INTENT(IN)    :: HLAMBDA3 ! formulation for lambda3 coeff
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PPABS  ! pressure (Pa)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PZZ    ! height of model levels (m)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRHODREF
REAL, DIMENSION(KIU,KJU,KKU), INTENT(INOUT) :: PT     ! grid scale T  (K)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRV_IN ! grid scale water vapor mixing ratio (kg/kg) in input
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PRV_OUT! grid scale water vapor mixing ratio (kg/kg) in output
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRC_IN ! grid scale r_c mixing ratio (kg/kg) in input
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PRC_OUT! grid scale r_c mixing ratio (kg/kg) in output
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRI_IN ! grid scale r_i (kg/kg) in input
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PRI_OUT! grid scale r_i (kg/kg) in output
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRR    ! grid scale mixing ration of rain (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRS    ! grid scale mixing ration of snow (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PRG    ! grid scale mixing ration of graupel (kg/kg)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(IN)    :: PSIGS  ! Sigma_s from turbulence scheme
REAL, DIMENSION(:,:,:),       INTENT(IN)    :: PMFCONV! convective mass flux (kg /s m^2)
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PCLDFR ! cloud fraction
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PSIGRC ! s r_c / sig_s^2

LOGICAL, INTENT(IN)                         :: OUSERI ! logical switch to compute both
                                                      ! liquid and solid condensate (OUSERI=.TRUE.)
                                                      ! or only solid condensate (OUSERI=.FALSE.)
LOGICAL, INTENT(IN)                         :: OSIGMAS! use present global Sigma_s values
                                                      ! or that from turbulence scheme
LOGICAL, INTENT(IN)                         :: OCND2  ! logical switch to sparate liquid and ice
                                                      ! more rigid (DEFAULT value : .FALSE.)
LOGICAL, INTENT(IN)                         :: LHGT_QS! logical switch for height dependent VQSIGSAT
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PICLDFR! ice cloud fraction
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PWCLDFR! water or mixed-phase cloud fraction
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PSSIO  ! Super-saturation with respect to ice in the  
                                                      ! supersaturated fraction
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PSSIU  ! Sub-saturation with respect to ice in the  
                                                      ! subsaturated fraction
REAL, DIMENSION(KIU,KJU,KKU), INTENT(OUT)   :: PIFR   ! Ratio cloud ice moist part
REAL, DIMENSION(KIU,KJU),     INTENT(IN)    :: PSIGQSAT ! use an extra "qsat" variance contribution (OSIGMAS case)
                                                        ! multiplied by PSIGQSAT

REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PLV    ! Latent heat L_v
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PLS    ! Latent heat L_s
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(IN)    :: PCPH   ! Specific heat C_ph
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLC_HRC
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLC_HCF ! cloud fraction
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLI_HRI
REAL, DIMENSION(KIU,KJU,KKU), OPTIONAL, INTENT(OUT)   :: PHLI_HCF
REAL, DIMENSION(KIU,KJU),   OPTIONAL, INTENT(IN)   :: PICE_CLD_WGT
END SUBROUTINE CONDENSATION
!
END INTERFACE
!
END MODULE MODI_CONDENSATION
