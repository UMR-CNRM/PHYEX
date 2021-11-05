!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODI_COMPUTE_MF_CLOUD_STAT
!    ############################
!
INTERFACE
!     #################################################################
      SUBROUTINE COMPUTE_MF_CLOUD_STAT(KKA, KKB, KKE, KKU, KKL, KRR, KRRL, KRRI,&
                            &PFRAC_ICE,&
                            &PTHLM, PRTM, PPABSM, PRM,&
                            &PDZZ, PTHM, PEXNM,&
                            &PEMF, PTHL_UP, PRT_UP,&
                            &PSIGMF)
!     #################################################################
!!
!
!*               1.1  Declaration of Arguments
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL                     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   :: KRR                     ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL                    ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI                    ! number of ice water var.
REAL, DIMENSION(:,:),   INTENT(IN)   :: PFRAC_ICE               ! liquid/ice fraction
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHLM, PRTM             ! cons. var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   :: PPABSM                  ! Pressure at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PRM                     ! water var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   :: PDZZ
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHM                    ! environement
REAL, DIMENSION(:,:),   INTENT(IN)   :: PEXNM
REAL, DIMENSION(:,:),   INTENT(IN)   :: PEMF                    ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHL_UP, PRT_UP         ! rc,w,Mass Flux,Thetal,rt
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PSIGMF                  ! SQRT(variance) for statistical cloud scheme


END SUBROUTINE COMPUTE_MF_CLOUD_STAT

END INTERFACE
!
END MODULE MODI_COMPUTE_MF_CLOUD_STAT
!     ######spl
      SUBROUTINE COMPUTE_MF_CLOUD_STAT(KKA, KKB, KKE, KKU, KKL, KRR, KRRL, KRRI,&
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
USE MODD_PARAM_MFSHALL_n, ONLY :  XTAUSIGMF
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT
!
USE MODI_SHUMAN_MF
USE MODI_COMPUTE_FUNCTION_THERMO_MF
!
!
IMPLICIT NONE
!
!*                    0.1  Declaration of Arguments
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   :: KKL                     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   :: KRR                     ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL                    ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI                    ! number of ice water var.
REAL, DIMENSION(:,:),   INTENT(IN)   :: PFRAC_ICE               ! liquid/ice fraction
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHLM, PRTM             ! cons. var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   :: PPABSM                  ! Pressure at time t-1
REAL, DIMENSION(:,:,:), INTENT(IN)   :: PRM                     ! water var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   :: PDZZ
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHM                    ! environement
REAL, DIMENSION(:,:),   INTENT(IN)   :: PEXNM
REAL, DIMENSION(:,:),   INTENT(IN)   :: PEMF                    ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHL_UP, PRT_UP         ! rc,w,Mass Flux,Thetal,rt
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PSIGMF                  ! SQRT(variance) for statistical cloud scheme
!
!*                    0.1  Declaration of local variables
!
!
REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2)) :: ZFLXZ
REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2)) :: ZT
REAL, DIMENSION(SIZE(PTHLM,1),SIZE(PTHLM,2)) :: ZAMOIST, ZATHETA
!
!*                    0.2 initialisation
!
!
!----------------------------------------------------------------------------
!
!*      1. COMPUTE SIGMA_MF (saturation deviation variance)
!          Soares et al (2004) formulation
!          ------------------------------------------------
!
! Thermodynamics functions
CALL COMPUTE_FUNCTION_THERMO_MF( KRR,KRRL,KRRI,                   &
                                 PTHM,PRM,PEXNM,PFRAC_ICE,PPABSM, &
                                 ZT,ZAMOIST,ZATHETA               )
!
IF (KRRL > 0)  THEN
!
!*       1.1 contribution from <THl THl>
!

!
    ZFLXZ(:,:) = -2 * XTAUSIGMF * PEMF(:,:)*(PTHL_UP(:,:)-MZM_MF(KKA,KKU,KKL,PTHLM(:,:))) * &
                      GZ_M_W_MF(KKA,KKU,KKL,PTHLM(:,:),PDZZ(:,:))
!
!   Avoid negative values
    ZFLXZ(:,:) = MAX(0.,ZFLXZ(:,:))


    PSIGMF(:,:) = MZF_MF(KKA,KKU,KKL,ZFLXZ(:,:)) * ZATHETA(:,:)**2

!
!
!*       1.2  contribution from <Rnp Rnp>
!
!
!
!
    ZFLXZ(:,:) = -2 * XTAUSIGMF * PEMF(:,:)*(PRT_UP(:,:)-MZM_MF(KKA,KKU,KKL,PRTM(:,:))) * &
                      GZ_M_W_MF(KKA,KKU,KKL,PRTM(:,:),PDZZ(:,:))
!
!   Avoid negative values
    ZFLXZ(:,:) = MAX(0.,ZFLXZ(:,:))
!

    PSIGMF(:,:) = PSIGMF(:,:) + ZAMOIST(:,:) **2 * MZF_MF(KKA,KKU,KKL,ZFLXZ(:,:))
!
!        1.3  Vertical part of Sigma_s
!
  PSIGMF(:,:) =  SQRT( MAX (PSIGMF(:,:) , 0.) )
ELSE
  PSIGMF(:,:) = 0.
END IF
!
!
END SUBROUTINE COMPUTE_MF_CLOUD_STAT
