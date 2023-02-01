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
INTEGER :: JIJ, JK
INTEGER :: IIJB,IIJE ! physical horizontal domain indices
INTEGER :: IKT
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!*                    0.2 initialisation
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_STAT',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
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
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZFLXZ(IIJB:IIJE,1:IKT) = -2 * CSTURB%XCTV* PARAMMF%XTAUSIGMF * PEMF(IIJB:IIJE,1:IKT)* &
                           & (PTHL_UP(IIJB:IIJE,1:IKT)-ZFLXZ(IIJB:IIJE,1:IKT)) * ZWK(IIJB:IIJE,1:IKT)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZFLXZ(IIJB:IIJE,1:IKT) = -2 * PARAMMF%XTAUSIGMF * PEMF(IIJB:IIJE,1:IKT)* &
                           & (PTHL_UP(IIJB:IIJE,1:IKT)-ZFLXZ(IIJB:IIJE,1:IKT)) * ZWK(IIJB:IIJE,1:IKT)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    END IF
    !
    !   Avoid negative values
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZFLXZ(IIJB:IIJE,1:IKT) = MAX(0.,ZFLXZ(IIJB:IIJE,1:IKT))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)


    CALL MZF_MF(D, ZFLXZ(:,:), PSIGMF(:,:))
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PSIGMF(IIJB:IIJE,1:IKT) = PSIGMF(IIJB:IIJE,1:IKT) * ZATHETA(IIJB:IIJE,1:IKT)**2
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

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
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZFLXZ2(IIJB:IIJE,1:IKT) = -2 * CSTURB%XCTV * PARAMMF%XTAUSIGMF * PEMF(IIJB:IIJE,1:IKT)* &
                           & (PRT_UP(IIJB:IIJE,1:IKT)-ZFLXZ2(IIJB:IIJE,1:IKT)) * ZWK2(IIJB:IIJE,1:IKT)
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ELSE
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZFLXZ2(IIJB:IIJE,1:IKT) = -2 * PARAMMF%XTAUSIGMF * PEMF(IIJB:IIJE,1:IKT)* &
                           & (PRT_UP(IIJB:IIJE,1:IKT)-ZFLXZ2(IIJB:IIJE,1:IKT)) * ZWK2(IIJB:IIJE,1:IKT) 
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    END IF
    !
    !   Avoid negative values
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ZFLXZ2(IIJB:IIJE,1:IKT) = MAX(0.,ZFLXZ2(IIJB:IIJE,1:IKT))
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)

    CALL MZF_MF(D, ZFLXZ2(:,:), ZWK2(:,:))
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PSIGMF(IIJB:IIJE,1:IKT) = PSIGMF(IIJB:IIJE,1:IKT) + ZAMOIST(IIJB:IIJE,1:IKT) **2 *ZWK2(IIJB:IIJE,1:IKT)
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    IF (OSTATNW) THEN
      !wc Now including convection covariance contribution in case of OSTATNW=TRUE
      !
      !       1.2.2 contribution from <Rnp Thl>
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      ZFLXZ3(IIJB:IIJE,1:IKT) = - CSTURB%XCTV * PARAMMF%XTAUSIGMF * &
                    (PEMF(IIJB:IIJE,1:IKT)*(PRT_UP(IIJB:IIJE,1:IKT)-ZFLXZ2(IIJB:IIJE,1:IKT)) * &
                                   ZWK(IIJB:IIJE,1:IKT) + &
                                   PEMF(IIJB:IIJE,1:IKT)*(PTHL_UP(IIJB:IIJE,1:IKT)-ZFLXZ(IIJB:IIJE,1:IKT)) * &
                                   ZWK2(IIJB:IIJE,1:IKT))
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      CALL MZF_MF(D, ZFLXZ3, ZFLXZ)
      !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
      PSIGMF(IIJB:IIJE,1:IKT) = PSIGMF(IIJB:IIJE,1:IKT) - &
                                MIN(0.,2.*ZAMOIST(IIJB:IIJE,1:IKT)*ZATHETA(IIJB:IIJE,1:IKT)*&
                                      &ZFLXZ(IIJB:IIJE,1:IKT))
      !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    ENDIF
!
!        1.3  Vertical part of Sigma_s
!
    !$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
    PSIGMF(IIJB:IIJE,1:IKT) =  SQRT( MAX (PSIGMF(IIJB:IIJE,1:IKT) , 0.) )
    !$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ELSE
  PSIGMF(:,:) = 0.
END IF
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_STAT',1,ZHOOK_HANDLE)
!
END SUBROUTINE COMPUTE_MF_CLOUD_STAT
END MODULE MODE_COMPUTE_MF_CLOUD_STAT
