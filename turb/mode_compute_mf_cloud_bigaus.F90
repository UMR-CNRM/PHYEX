!MNH_LIC Copyright 2011-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODE_COMPUTE_MF_CLOUD_BIGAUS
!    ###################################
!
IMPLICIT NONE
CONTAINS
      SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS(D, CST, PARAMMF,&
                                  PEMF, PDEPTH,&
                                  PRT_UP, PTHV_UP, PFRAC_ICE_UP, PRSAT_UP,&
                                  PRTM, PTHM, PTHVM,&
                                  PDZZ, PZZ, PRHODREF,&
                                  PRC_MF, PRI_MF, PCF_MF)
!     #################################################################
!!
!!****  *COMPUTE_MF_CLOUD_BIGAUS* -
!!       compute diagnostic subgrid cumulus cloud caracteristics with a statistical scheme
!!       based on a bi-gaussian PDF. In this routine, we only compute the shallow convection
!!       part of this bi-gaussian
!!
!!    PURPOSE
!!    -------
!!****  With this option, a formulation for the computation of the variance of the departure
!!      to saturation is proposed. This variance is used to compute the cloud fraction and
!!      the mean convective cloud content from the bi-gaussian PDF proposed by E. Perraud
!!
!
!!**  METHOD
!!    ------
!!      Updraft variables are used to diagnose the variance
!!      Perraud et al (2011)
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
!!      S. Riette Jun 2019: remove unused PRC_UP and PRI_UP, use SIGN in ERFC computation
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_CST,             ONLY: CST_t
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALL_t
!
USE MODI_SHUMAN_MF, ONLY: MZF_MF, GZ_M_W_MF
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
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PEMF                    ! updraft characteritics
REAL, DIMENSION(D%NIJT),     INTENT(IN)   :: PDEPTH                  ! Deepness of cloud
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PTHV_UP, PRSAT_UP, PRT_UP ! updraft characteritics
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PFRAC_ICE_UP            ! liquid/ice fraction in updraft
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PTHM, PRTM, PTHVM       ! env. var. at t-dt
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PDZZ, PZZ
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PRHODREF
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PRC_MF, PRI_MF          ! cloud content
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PCF_MF                  ! and cloud fraction for MF scheme
!
!*                    0.1  Declaration of local variables
!
!
REAL, DIMENSION(D%NIJT,D%NKT) :: ZGRAD_Z_RT, &            !
                                            & ZALPHA_UP_M, &           ! Variables used to compute variance
                                            & ZSIGMF                   ! and sqrt(variance)
REAL, DIMENSION(D%NIJT)              :: ZOMEGA_UP_M              !
REAL, DIMENSION(D%NIJT,D%NKT) :: ZW1 ! working array
INTEGER                                    :: JI, JK  ! loop control
REAL, DIMENSION(D%NIJT,D%NKT) :: ZEMF_M, ZTHV_UP_M, &   !
                                            & ZRSAT_UP_M, ZRT_UP_M,& ! Interpolation on mass points
                                            & ZFRAC_ICE_UP_M         !
REAL, DIMENSION(D%NIJT,D%NKT) :: ZCOND ! condensate
REAL, DIMENSION(D%NIJT,D%NKT) :: ZA, ZGAM ! used for integration
REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_BIGAUS',0,ZHOOK_HANDLE)

!Computation is done on mass points
!----------------------------------------------------------------------------
!
!*      1. Computation of the variance
!          ------------------------------------------------
!
!
!Vertical gradient of RT, result on mass points
CALL GZ_M_W_MF(D, PRTM(:,:), PDZZ(:,:), ZW1(:,:))
CALL MZF_MF(D, ZW1(:,:), ZGRAD_Z_RT(:,:))

!Interpolation on mass points
CALL MZF_MF(D, PTHV_UP(:,:), ZTHV_UP_M(:,:))
CALL MZF_MF(D, PRSAT_UP(:,:), ZRSAT_UP_M(:,:))
CALL MZF_MF(D, PRT_UP(:,:), ZRT_UP_M(:,:))
CALL MZF_MF(D, PEMF(:,:), ZEMF_M(:,:))
CALL MZF_MF(D, PFRAC_ICE_UP(:,:), ZFRAC_ICE_UP_M(:,:))

!computation of omega star up
ZOMEGA_UP_M(:)=0.
DO JK=D%NKB,D%NKE-D%NKL,D%NKL
  DO JI=D%NIJB,D%NIJE 
  !Vertical integration over the entire column but only buoyant points are used
  !ZOMEGA_UP_M(JI)=ZOMEGA_UP_M(JI) + &
  !                ZEMF_M(JI,JK) * &
  !                MAX(0.,(ZTHV_UP_M(JI,JK)-PTHVM(JI,JK))) * &
  !                (PZZ(JI,JK+KKL)-PZZ(JI,JK)) / &
  !                (PTHM(JI,JK) * PRHODREF(JI,JK))
    
  !Vertical integration over the entire column
    ZOMEGA_UP_M(JI)=ZOMEGA_UP_M(JI) + &
    ZEMF_M(JI,JK) * &
    (ZTHV_UP_M(JI,JK)-PTHVM(JI,JK)) * &
    (PZZ(JI,JK+D%NKL)-PZZ(JI,JK)) / &
    (PTHM(JI,JK) * PRHODREF(JI,JK))
  ENDDO
ENDDO
DO JI=D%NIJB,D%NIJE 
  ZOMEGA_UP_M(JI)=MAX(ZOMEGA_UP_M(JI), 1.E-20)
  ZOMEGA_UP_M(JI)=(CST%XG*ZOMEGA_UP_M(JI))**(1./3.)
ENDDO

!computation of alpha up
DO JK=D%NKA,D%NKU,D%NKL
  DO JI=D%NIJB,D%NIJE 
    ZALPHA_UP_M(JI,JK)=ZEMF_M(JI,JK)/(PARAMMF%XALPHA_MF*PRHODREF(JI,JK)*ZOMEGA_UP_M(JI))
    ZALPHA_UP_M(JI,JK)=MAX(0., MIN(ZALPHA_UP_M(JI,JK), 1.))
  ENDDO
ENDDO

!computation of sigma of the distribution
DO JK=D%NKA,D%NKU,D%NKL
  DO JI=D%NIJB,D%NIJE 
    ZSIGMF(JI,JK)=ZEMF_M(JI,JK) * &
    (ZRT_UP_M(JI,JK) - PRTM(JI,JK)) * &
    PDEPTH(JI) * ZGRAD_Z_RT(JI,JK) / &
    (PARAMMF%XSIGMA_MF * ZOMEGA_UP_M(JI) * PRHODREF(JI,JK))
  ENDDO
ENDDO
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
    ZSIGMF(JI,JK)=SQRT(MAX(ABS(ZSIGMF(JI,JK)), 1.E-40))
  ENDDO
ENDDO
!
!*      2. PDF integration
!          ------------------------------------------------
!
DO JK=1,D%NKT 
  DO JI=D%NIJB,D%NIJE 
!The mean of the distribution is ZRT_UP
!Computation of ZA and ZGAM (=efrc(ZA)) coefficient
    ZA(JI,JK)=(ZRSAT_UP_M(JI,JK)-ZRT_UP_M(JI,JK))/&
    &(sqrt(2.)*ZSIGMF(JI,JK))
    
!Approximation of erf function
    ZGAM(JI,JK)=1-SIGN(1., ZA(JI,JK))*SQRT(1-EXP(-4*ZA(JI,JK)**2/CST%XPI))
    
!computation of cloud fraction
    PCF_MF(JI,JK)=MAX( 0., MIN(1.,0.5*ZGAM(JI,JK) * ZALPHA_UP_M(JI,JK)))
    
!computation of condensate, then PRC and PRI
    ZCOND(JI,JK)=(EXP(-ZA(JI,JK)**2)-&
    &ZA(JI,JK)*SQRT(CST%XPI)*ZGAM(JI,JK))* &
    &ZSIGMF(JI,JK)/SQRT(2.*CST%XPI) * ZALPHA_UP_M(JI,JK)
    ZCOND(JI,JK)=MAX(ZCOND(JI,JK), 0.) !due to approximation of ZGAM value, ZCOND could be slightly negative
    PRC_MF(JI,JK)=(1.-ZFRAC_ICE_UP_M(JI,JK)) * ZCOND(JI,JK)
    PRI_MF(JI,JK)=(   ZFRAC_ICE_UP_M(JI,JK)) * ZCOND(JI,JK)
  ENDDO
ENDDO
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_BIGAUS',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS
END MODULE MODE_COMPUTE_MF_CLOUD_BIGAUS
