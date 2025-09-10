!MNH_LIC Copyright 2009-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODE_COMPUTE_MF_CLOUD
!    ############################
!
IMPLICIT NONE
CONTAINS
!
!     ######spl
      SUBROUTINE COMPUTE_MF_CLOUD(D, CST, TURBN, PARAMMF, ICEP, OSTATNW,    &
                                  KRR, KRRL, KRRI,                          &
                                  PFRAC_ICE,                                &
                                  PRV_UP, PRC_UP, PRI_UP, PEMF,             &
                                  PTHL_UP, PRT_UP, PFRAC_UP, PTH_UP,        &
                                  PEXNM, PTHLM, PRTM, PTHM, PRM,            &
                                  PDZZ, KKLCL,                              &
                                  PPABSM, PRHODREF,                         &
                                  PRC_MF, PRI_MF, PCF_MF, PSIGMF,           &
                                  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF,   &
                                  PWEIGHT_MF_CLOUD)

!     #################################################################
!!
!!****  *COMPUTE_MF_CLOUD* - 
!!       compute diagnostic subgrid cumulus cloud caracteristics
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to compute the cloud fraction and 
!!      the mean cloud content associated with clouds described by the 
!!      mass flux scheme
!!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     REFERENCE
!!     ---------
!!
!!
!!     AUTHOR
!!     ------
!!
!!    MODIFICATIONS
!!    -------------
!!      Original
!!      S. Riette Dec 2010 BIGA case
!!      S. Riette Aug 2011 code is split into subroutines
!!      S. Riette Jan 2012: support for both order of vertical levels
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!!      A. Marcel Jan 2025: bi-Gaussian PDF and associated subgrid precipitation
!!      A. Marcel Jan 2025: relaxation of the small fraction assumption
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,                ONLY: DIMPHYEX_t
USE MODD_CST,                     ONLY: CST_t
USE MODD_TURB_n,                  ONLY: TURB_t
USE MODD_PARAM_MFSHALL_n,         ONLY: PARAM_MFSHALL_t
USE MODD_RAIN_ICE_PARAM_n,        ONLY: RAIN_ICE_PARAM_t
USE MODE_MSG,                     ONLY: PRINT_MSG, NVERB_FATAL
USE MODE_COMPUTE_MF_CLOUD_DIRECT, ONLY: COMPUTE_MF_CLOUD_DIRECT
USE MODE_COMPUTE_MF_CLOUD_STAT,   ONLY: COMPUTE_MF_CLOUD_STAT
USE MODE_COMPUTE_MF_CLOUD_BIGAUS, ONLY: COMPUTE_MF_CLOUD_BIGAUS
USE YOMHOOK,                      ONLY: LHOOK, DR_HOOK, JPHOOK

IMPLICIT NONE

!*                    1.1  Declaration of Arguments
!
!
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(TURB_t),           INTENT(IN)   :: TURBN
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
TYPE(RAIN_ICE_PARAM_t), INTENT(IN)   :: ICEP
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
INTEGER,                INTENT(IN)   ::  KRRL         ! number of liquid water var.
INTEGER,                INTENT(IN)   ::  KRRI         ! number of ice water var.
LOGICAL,                INTENT(IN)   :: OSTATNW      ! cloud scheme inclues convect. covar. contrib
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PFRAC_ICE    ! liquid/ice fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PRV_UP, PRC_UP, PRI_UP, PEMF! updraft characteritics
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHL_UP, PRT_UP   ! rc,w,Mass Flux,Thetal,rt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PFRAC_UP          ! Updraft Fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTH_UP            ! updraft theta
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PEXNM             ! exner function
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHLM, PRTM       ! cons. var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PTHM              ! theta
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN)   ::  PRM               ! water var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PDZZ
INTEGER, DIMENSION(D%NIJT),  INTENT(IN)   ::  KKLCL             ! index of updraft condensation level
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   ::  PPABSM, PRHODREF  ! environement
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PRC_MF, PRI_MF    ! cloud content (INPUT=environment, OUTPUT=conv. cloud)
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PCF_MF            ! and cloud fraction for MF scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PSIGMF            ! SQRT(variance) for statistical cloud scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF ! low/high cloud diagnostics
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  ::  PWEIGHT_MF_CLOUD ! weight coefficient for the mass-flux cloud
!
!                       1.2  Declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!------------------------------------------------------------------------

!                     1. INITIALISATION
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD',0,ZHOOK_HANDLE)
!
!                     2.1 Internal domain
PRC_MF(:,:) = 0.
PRI_MF(:,:) = 0.
PCF_MF(:,:) = 0.
PSIGMF(:,:) = 0.
PHLC_HCF(:,:) = 0.
PHLC_HRC(:,:) = 0.
PHLI_HCF(:,:) = 0.
PHLI_HRI(:,:) = 0.
PWEIGHT_MF_CLOUD(:,:) = 0.

IF (PARAMMF%CMF_CLOUD == 'DIRE') THEN
  !Direct cloud scheme
  CALL COMPUTE_MF_CLOUD_DIRECT(D, PARAMMF, &
                              &KKLCL(:), PFRAC_UP(:,:), PRC_UP(:,:), PRI_UP(:,:),&
                              &PRC_MF(:,:), PRI_MF(:,:), PCF_MF(:,:), PWEIGHT_MF_CLOUD(:,:))
  !
ELSEIF (PARAMMF%CMF_CLOUD == 'STAT') THEN
  !Statistical scheme using the PDF proposed by Bougeault (81, 82) and
  !Bechtold et al (95).
  CALL COMPUTE_MF_CLOUD_STAT(D, CST, TURBN, PARAMMF, &
                            &KRR, KRRL, KRRI, OSTATNW, &
                            &PFRAC_ICE,&
                            &PTHLM, PRTM, PPABSM, PRM,&
                            &PDZZ, PTHM, PEXNM,&
                            &PEMF, PTHL_UP, PRT_UP,&
                            &PSIGMF)
ELSEIF (PARAMMF%CMF_CLOUD == 'BIGA') THEN
  !Statistical scheme using the bi-gaussian PDF
  CALL COMPUTE_MF_CLOUD_BIGAUS(D, CST, PARAMMF, ICEP, KRR, &
                              &PRT_UP, PRV_UP, PRC_UP, PRI_UP, PTH_UP, PFRAC_UP, &
                              &PRTM, PTHM, PRM, &
                              &PRHODREF, PEXNM, PPABSM, &
                              &PRC_MF, PRI_MF, PCF_MF, PSIGMF, &
                              &PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF, PWEIGHT_MF_CLOUD(:,:))
  !
ELSEIF  (PARAMMF%CMF_CLOUD == 'NONE') THEN
  ! No CONVECTIVE CLOUD SCHEME
  ! Nothing to do: PRC_MF, PRI_MF, PCF_MF, PSIGMF are already filled with zero
ELSE
  CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'COMPUTE_MF_CLOUD', &
                 'Shallow convection cloud scheme not valid: PARAMMF%CMF_CLOUD='//PARAMMF%CMF_CLOUD)
ENDIF

IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_MF_CLOUD
END MODULE MODE_COMPUTE_MF_CLOUD
