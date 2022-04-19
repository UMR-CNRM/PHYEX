!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODE_COMPUTE_MF_CLOUD_STAT
!    ############################
!
IMPLICIT NONE
CONTAINS
!     ######spl
      SUBROUTINE COMPUTE_MF_CLOUD_STAT(D, CST, PARAMMF, KRR, KRRL, KRRI,&
                            &PFRAC_ICE,&
                            &PTHLM, PRTM, PPABSM, PRM,&
                            &PDZZ, PTHM, PEXNM, &
                            &PEMF, PTHL_UP, PRT_UP,&
                            &PSIGMF)
!     #################################################################
!!
!!****  *COMPUTE_MF_CLOUD_STAT* -
!!       compute diagnostic subgrid cumulus cloud caracteristics with a statistical scheme
!!
!!    PURPOSE
!!    -------
!!****  With this option, a formulation for the computation of the variance of the departure
!!      to saturation is proposed.
!!
!
!!**  METHOD
!!    ------
!!      Updraft variables are used to diagnose the variance
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
!!     S. Riette moving of code previously in compute_mf_cloud code
!!
!!    MODIFICATIONS
!!    -------------
!!      Original 25 Aug 2011
!!      S. Riette Jan 2012: support for both order of vertical levels
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
!
USE MODI_SHUMAN_MF, ONLY: MZF_MF, MZM_MF, GZ_M_W_MF
USE MODE_COMPUTE_FUNCTION_THERMO_MF, ONLY: COMPUTE_FUNCTION_THERMO_MF
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*                    0.1  Declaration of Arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(CST_t),            INTENT(IN)   :: CST
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
INTEGER,                INTENT(IN)   :: KRR                     ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL                    ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI                    ! number of ice water var.
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PFRAC_ICE               ! liquid/ice fraction
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PTHLM, PRTM             ! cons. var. at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PPABSM                  ! Pressure at time t-1
REAL, DIMENSION(D%NIT,D%NKT,KRR), INTENT(IN)   :: PRM                     ! water var. at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PDZZ
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PTHM                    ! environement
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PEXNM
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PEMF                    ! updraft characteritics
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PTHL_UP, PRT_UP         ! rc,w,Mass Flux,Thetal,rt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT)  :: PSIGMF                  ! SQRT(variance) for statistical cloud scheme
!
!*                    0.1  Declaration of local variables
!
!
REAL, DIMENSION(D%NIT,D%NKT) :: ZFLXZ
REAL, DIMENSION(D%NIT,D%NKT) :: ZT
REAL, DIMENSION(D%NIT,D%NKT) :: ZAMOIST, ZATHETA
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*                    0.2 initialisation
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_STAT',0,ZHOOK_HANDLE)
!
!----------------------------------------------------------------------------
!
!*      1. COMPUTE SIGMA_MF (saturation deviation variance)
!          Soares et al (2004) formulation
!          ------------------------------------------------
!
! Thermodynamics functions
CALL COMPUTE_FUNCTION_THERMO_MF( D, CST, KRR,KRRL,KRRI,                   &
                                 PTHM,PRM,PEXNM,PFRAC_ICE,PPABSM, &
                                 ZT,ZAMOIST,ZATHETA               )
!
IF (KRRL > 0)  THEN
!
!*       1.1 contribution from <THl THl>
!

!
    ZFLXZ(:,:) = -2 * PARAMMF%XTAUSIGMF * PEMF(:,:)*(PTHL_UP(:,:)-MZM_MF(PTHLM(:,:), D%NKA, D%NKU, D%NKL)) * &
                      GZ_M_W_MF(PTHLM(:,:),PDZZ(:,:), D%NKA, D%NKU, D%NKL)
!
!   Avoid negative values
    ZFLXZ(:,:) = MAX(0.,ZFLXZ(:,:))


    PSIGMF(:,:) = MZF_MF(ZFLXZ(:,:), D%NKA, D%NKU, D%NKL) * ZATHETA(:,:)**2

!
!
!*       1.2  contribution from <Rnp Rnp>
!
!
!
!
    ZFLXZ(:,:) = -2 * PARAMMF%XTAUSIGMF * PEMF(:,:)*(PRT_UP(:,:)-MZM_MF(PRTM(:,:), D%NKA, D%NKU, D%NKL)) * &
                      GZ_M_W_MF(PRTM(:,:),PDZZ(:,:), D%NKA, D%NKU, D%NKL)
!
!   Avoid negative values
    ZFLXZ(:,:) = MAX(0.,ZFLXZ(:,:))
!

    PSIGMF(:,:) = PSIGMF(:,:) + ZAMOIST(:,:) **2 * MZF_MF(ZFLXZ(:,:), D%NKA, D%NKU, D%NKL)
!
!        1.3  Vertical part of Sigma_s
!
  PSIGMF(:,:) =  SQRT( MAX (PSIGMF(:,:) , 0.) )
ELSE
  PSIGMF(:,:) = 0.
END IF
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_STAT',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_MF_CLOUD_STAT
END MODULE MODE_COMPUTE_MF_CLOUD_STAT
