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
      SUBROUTINE COMPUTE_MF_CLOUD_STAT(D, CST, CSTURB, PARAMMF, &
                            &KRR, KRRL, KRRI, OSTATNW,     &
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
!!      Wim de Rooy June 2019: update statistical cloud scheme (now including
!!                             covariance term for MF contribution)
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
USE MODD_CTURB,           ONLY: CSTURB_t

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
TYPE(CSTURB_t),         INTENT(IN)   :: CSTURB
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
LOGICAL,                INTENT(IN)   :: OSTATNW      ! cloud scheme inclues convect. covar. contrib
INTEGER,                INTENT(IN)   :: KRR                     ! number of moist var.
INTEGER,                INTENT(IN)   :: KRRL                    ! number of liquid water var.
INTEGER,                INTENT(IN)   :: KRRI                    ! number of ice water var.
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PFRAC_ICE               ! liquid/ice fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PTHLM, PRTM             ! cons. var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PPABSM                  ! Pressure at time t-1
REAL, DIMENSION(D%NIJT,D%NKT,KRR), INTENT(IN)   :: PRM                     ! water var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PDZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PTHM                    ! environement
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PEXNM
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PEMF                    ! updraft characteritics
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PTHL_UP, PRT_UP         ! rc,w,Mass Flux,Thetal,rt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PSIGMF                  ! SQRT(variance) for statistical cloud scheme
!
!*                    0.1  Declaration of local variables
!
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZFLXZ,ZFLXZ2,ZFLXZ3
REAL, DIMENSION(D%NIJT,D%NKT) :: ZT
REAL, DIMENSION(D%NIJT,D%NKT) :: ZAMOIST, ZATHETA
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWK,ZWK2
INTEGER :: JI, JK
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
CALL COMPUTE_FUNCTION_THERMO_MF( D, CST, KRR,KRRL,KRRI,OSTATNW,   &
                                 PTHM,PRM,PEXNM,PFRAC_ICE,PPABSM, &
                                 ZT,ZAMOIST,ZATHETA               )
!
IF (KRRL > 0)  THEN
!
!*       1.1 contribution from <THl THl>
!

!
    CALL MZM_MF(D, PTHLM(:,:), ZFLXZ(:,:))
    CALL GZ_M_W_MF(D, PTHLM(:,:), PDZZ(:,:), ZWK(:,:))
    IF (OSTATNW) THEN
    DO JK=1,D%NKT 
      DO JI=D%NIJB,D%NIJE 
        ZFLXZ(JI,JK) = -2 * CSTURB%XCTV* PARAMMF%XTAUSIGMF * PEMF(JI,JK)* &
        & (PTHL_UP(JI,JK)-ZFLXZ(JI,JK)) * ZWK(JI,JK)
      ENDDO
    ENDDO
    ELSE
    DO JK=1,D%NKT 
      DO JI=D%NIJB,D%NIJE 
        ZFLXZ(JI,JK) = -2 * PARAMMF%XTAUSIGMF * PEMF(JI,JK)* &
        & (PTHL_UP(JI,JK)-ZFLXZ(JI,JK)) * ZWK(JI,JK)
      ENDDO
    ENDDO
    END IF
    !
    !   Avoid negative values
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      ZFLXZ(JI,JK) = MAX(0.,ZFLXZ(JI,JK))
    ENDDO
  ENDDO


    CALL MZF_MF(D, ZFLXZ(:,:), PSIGMF(:,:))
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      PSIGMF(JI,JK) = PSIGMF(JI,JK) * ZATHETA(JI,JK)**2
    ENDDO
  ENDDO

!
!
!*       1.2  contribution from <Rnp Rnp>
!
!
!
!
    CALL MZM_MF(D, PRTM(:,:), ZFLXZ2(:,:))
    CALL GZ_M_W_MF(D, PRTM(:,:), PDZZ(:,:), ZWK2(:,:))
    IF (OSTATNW) THEN
    DO JK=1,D%NKT 
      DO JI=D%NIJB,D%NIJE 
        ZFLXZ2(JI,JK) = -2 * CSTURB%XCTV * PARAMMF%XTAUSIGMF * PEMF(JI,JK)* &
        & (PRT_UP(JI,JK)-ZFLXZ2(JI,JK)) * ZWK2(JI,JK)
      ENDDO
    ENDDO
    ELSE
    DO JK=1,D%NKT 
      DO JI=D%NIJB,D%NIJE 
        ZFLXZ2(JI,JK) = -2 * PARAMMF%XTAUSIGMF * PEMF(JI,JK)* &
        & (PRT_UP(JI,JK)-ZFLXZ2(JI,JK)) * ZWK2(JI,JK) 
      ENDDO
    ENDDO
    END IF
    !
    !   Avoid negative values
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      ZFLXZ2(JI,JK) = MAX(0.,ZFLXZ2(JI,JK))
    ENDDO
  ENDDO

    CALL MZF_MF(D, ZFLXZ2(:,:), ZWK2(:,:))
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      PSIGMF(JI,JK) = PSIGMF(JI,JK) + ZAMOIST(JI,JK) **2 *ZWK2(JI,JK)
    ENDDO
  ENDDO
    IF (OSTATNW) THEN
      !wc Now including convection covariance contribution in case of OSTATNW=TRUE
      !
      !       1.2.2 contribution from <Rnp Thl>
    DO JK=1,D%NKT 
      DO JI=D%NIJB,D%NIJE 
        ZFLXZ3(JI,JK) = - CSTURB%XCTV * PARAMMF%XTAUSIGMF * &
        (PEMF(JI,JK)*(PRT_UP(JI,JK)-ZFLXZ2(JI,JK)) * &
        ZWK(JI,JK) + &
        PEMF(JI,JK)*(PTHL_UP(JI,JK)-ZFLXZ(JI,JK)) * &
        ZWK2(JI,JK))
      ENDDO
    ENDDO
      CALL MZF_MF(D, ZFLXZ3, ZFLXZ)
    DO JK=1,D%NKT 
      DO JI=D%NIJB,D%NIJE 
        PSIGMF(JI,JK) = PSIGMF(JI,JK) - &
        MIN(0.,2.*ZAMOIST(JI,JK)*ZATHETA(JI,JK)*&
        &ZFLXZ(JI,JK))
      ENDDO
    ENDDO
    ENDIF
!
!        1.3  Vertical part of Sigma_s
!
  DO JK=1,D%NKT 
    DO JI=D%NIJB,D%NIJE 
      PSIGMF(JI,JK) =  SQRT( MAX (PSIGMF(JI,JK) , 0.) )
    ENDDO
  ENDDO
ELSE
  PSIGMF(:,:) = 0.
END IF
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_STAT',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_MF_CLOUD_STAT
END MODULE MODE_COMPUTE_MF_CLOUD_STAT
