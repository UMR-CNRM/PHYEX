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
!USE MODI_GAMMA_INC
!
!USE MODE_THERMO
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
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PEMF                    ! updraft characteritics
REAL, DIMENSION(D%NIT),     INTENT(IN)   :: PDEPTH                  ! Deepness of cloud
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PTHV_UP, PRSAT_UP, PRT_UP ! updraft characteritics
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PFRAC_ICE_UP            ! liquid/ice fraction in updraft
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PTHM, PRTM, PTHVM       ! env. var. at t-dt
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PDZZ, PZZ
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(IN)   :: PRHODREF
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT)  :: PRC_MF, PRI_MF          ! cloud content
REAL, DIMENSION(D%NIT,D%NKT),   INTENT(OUT)  :: PCF_MF                  ! and cloud fraction for MF scheme
!
!*                    0.1  Declaration of local variables
!
!
REAL, DIMENSION(D%NIT,D%NKT) :: ZGRAD_Z_RT, &            !
                                            & ZALPHA_UP_M, &           ! Variables used to compute variance
                                            & ZSIGMF                   ! and sqrt(variance)
REAL, DIMENSION(D%NIT)              :: ZOMEGA_UP_M              !
REAL, DIMENSION(D%NIT,D%NKT) :: ZW1 ! working array
INTEGER                                    :: JK  ! vertical loop control
REAL, DIMENSION(D%NIT,D%NKT) :: ZEMF_M, ZTHV_UP_M, &   !
                                            & ZRSAT_UP_M, ZRT_UP_M,& ! Interpolation on mass points
                                            & ZFRAC_ICE_UP_M         !
REAL, DIMENSION(D%NIT,D%NKT) :: ZCOND ! condensate
REAL, DIMENSION(D%NIT,D%NKT) :: ZA, ZGAM ! used for integration
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
ZW1(:,:)=GZ_M_W_MF(PRTM(:,:), PDZZ(:,:), D%NKA, D%NKU, D%NKL)
ZGRAD_Z_RT(:,:)=MZF_MF(ZW1(:,:), D%NKA, D%NKU, D%NKL)

!Interpolation on mass points
ZTHV_UP_M(:,:) = MZF_MF(PTHV_UP(:,:), D%NKA, D%NKU, D%NKL)
ZRSAT_UP_M(:,:)= MZF_MF(PRSAT_UP(:,:), D%NKA, D%NKU, D%NKL)
ZRT_UP_M(:,:)  = MZF_MF(PRT_UP(:,:), D%NKA, D%NKU, D%NKL)
ZEMF_M(:,:)    = MZF_MF(PEMF(:,:), D%NKA, D%NKU, D%NKL)
ZFRAC_ICE_UP_M(:,:) = MZF_MF(PFRAC_ICE_UP(:,:), D%NKA, D%NKU, D%NKL)

!computation of omega star up
ZOMEGA_UP_M(:)=0.
DO JK=D%NKB,D%NKE,D%NKL
  !Vertical integration over the entire column but only buoyant points are used
  !ZOMEGA_UP_M(:)=ZOMEGA_UP_M(:) + &
  !                ZEMF_M(:,JK) * &
  !                MAX(0.,(ZTHV_UP_M(:,JK)-PTHVM(:,JK))) * &
  !                (PZZ(:,JK+KKL)-PZZ(:,JK)) / &
  !                (PTHM(:,JK) * PRHODREF(:,JK))

  !Vertical integration over the entire column
  ZOMEGA_UP_M(:)=ZOMEGA_UP_M(:) + &
                 ZEMF_M(:,JK) * &
                 (ZTHV_UP_M(:,JK)-PTHVM(:,JK)) * &
                 (PZZ(:,JK+D%NKL)-PZZ(:,JK)) / &
                 (PTHM(:,JK) * PRHODREF(:,JK))
ENDDO
ZOMEGA_UP_M(:)=MAX(ZOMEGA_UP_M(:), 1.E-20)
ZOMEGA_UP_M(:)=(CST%XG*ZOMEGA_UP_M(:))**(1./3.)

!computation of alpha up
DO JK=D%NKA,D%NKU,D%NKL
  ZALPHA_UP_M(:,JK)=ZEMF_M(:,JK)/(PARAMMF%XALPHA_MF*PRHODREF(:,JK)*ZOMEGA_UP_M(:))
ENDDO
ZALPHA_UP_M(:,:)=MAX(0., MIN(ZALPHA_UP_M(:,:), 1.))

!computation of sigma of the distribution
DO JK=D%NKA,D%NKU,D%NKL
  ZSIGMF(:,JK)=ZEMF_M(:,JK) * &
               (ZRT_UP_M(:,JK) - PRTM(:,JK)) * &
               PDEPTH(:) * ZGRAD_Z_RT(:,JK) / &
               (PARAMMF%XSIGMA_MF * ZOMEGA_UP_M(:) * PRHODREF(:,JK))
ENDDO
ZSIGMF(:,:)=SQRT(MAX(ABS(ZSIGMF(:,:)), 1.E-40))
!
!*      2. PDF integration
!          ------------------------------------------------
!
!The mean of the distribution is ZRT_UP
!Computation of ZA and ZGAM (=efrc(ZA)) coefficient
ZA(:,:)=(ZRSAT_UP_M(:,:)-ZRT_UP_M(:,:))/(sqrt(2.)*ZSIGMF(:,:))

!Approximation of erf function
ZGAM(:,:)=1-SIGN(1., ZA(:,:))*SQRT(1-EXP(-4*ZA(:,:)**2/CST%XPI))

!computation of cloud fraction
PCF_MF(:,:)=MAX( 0., MIN(1.,0.5*ZGAM(:,:) * ZALPHA_UP_M(:,:)))

!computation of condensate, then PRC and PRI
ZCOND(:,:)=(EXP(-ZA(:,:)**2)-ZA(:,:)*SQRT(CST%XPI)*ZGAM(:,:))*ZSIGMF(:,:)/SQRT(2.*CST%XPI) * ZALPHA_UP_M(:,:)
ZCOND(:,:)=MAX(ZCOND(:,:), 0.) !due to approximation of ZGAM value, ZCOND could be slightly negative
PRC_MF(:,:)=(1.-ZFRAC_ICE_UP_M(:,:)) * ZCOND(:,:)
PRI_MF(:,:)=(   ZFRAC_ICE_UP_M(:,:)) * ZCOND(:,:)

IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_BIGAUS',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS
END MODULE MODE_COMPUTE_MF_CLOUD_BIGAUS
