!     ######spl
      SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS(KKA, KKB, KKE, KKU, KKL,&
                                  PRC_UP, PRI_UP, PEMF, PDEPTH,&
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
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
USE MODD_CMFSHALL, ONLY : XALPHA_MF, XSIGMA_MF
USE MODD_CST, ONLY  : XPI, XG
!
USE MODI_SHUMAN_MF
USE MODI_GAMMA_INC
!
USE MODE_THERMO
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
REAL, DIMENSION(:,:),   INTENT(IN)   :: PRC_UP,PRI_UP,PEMF      ! updraft characteritics
REAL, DIMENSION(:),     INTENT(IN)   :: PDEPTH                  ! Deepness of cloud
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHV_UP, PRSAT_UP, PRT_UP ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(IN)   :: PFRAC_ICE_UP            ! liquid/ice fraction in updraft
REAL, DIMENSION(:,:),   INTENT(IN)   :: PTHM, PRTM, PTHVM       ! env. var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   :: PDZZ, PZZ
REAL, DIMENSION(:,:),   INTENT(IN)   :: PRHODREF
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PRC_MF, PRI_MF          ! cloud content
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PCF_MF                  ! and cloud fraction for MF scheme
!
!*                    0.1  Declaration of local variables
!
!
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZGRAD_Z_RT, &            !
                                            & ZALPHA_UP_M, &           ! Variables used to compute variance
                                            & ZSIGMF                   ! and sqrt(variance)
REAL, DIMENSION(SIZE(PTHM,1))              :: ZOMEGA_UP_M              !
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZW1 ! working array
INTEGER                                    :: JK  ! vertical loop control
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZEMF_M, ZTHV_UP_M, &   !
                                            & ZRSAT_UP_M, ZRC_UP_M,& ! Interpolation on mass points
                                            & ZRI_UP_M, ZRT_UP_M,&   !
                                            & ZFRAC_ICE_UP_M         !
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZCOND ! condensate
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2)) :: ZA, ZGAM ! used for integration
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
ZW1(:,:)=GZ_M_W_MF(KKA,KKU,KKL, PRTM(:,:), PDZZ(:,:))
ZGRAD_Z_RT(:,:)=MZF_MF(KKA,KKU,KKL, ZW1(:,:))

!Interpolation on mass points
ZTHV_UP_M(:,:) = MZF_MF(KKA,KKU,KKL, PTHV_UP(:,:))
ZRSAT_UP_M(:,:)= MZF_MF(KKA,KKU,KKL, PRSAT_UP(:,:))
ZRC_UP_M(:,:)  = MZF_MF(KKA,KKU,KKL, PRC_UP(:,:))
ZRI_UP_M(:,:)  = MZF_MF(KKA,KKU,KKL, PRI_UP(:,:))
ZRT_UP_M(:,:)  = MZF_MF(KKA,KKU,KKL, PRT_UP(:,:))
ZEMF_M(:,:)    = MZF_MF(KKA,KKU,KKL, PEMF(:,:))
ZFRAC_ICE_UP_M(:,:) = MZF_MF(KKA,KKU,KKL, PFRAC_ICE_UP(:,:))

!computation of omega star up
ZOMEGA_UP_M(:)=0.
DO JK=KKB,KKE,KKL
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
                 (PZZ(:,JK+KKL)-PZZ(:,JK)) / &
                 (PTHM(:,JK) * PRHODREF(:,JK))
ENDDO
ZOMEGA_UP_M(:)=MAX(ZOMEGA_UP_M(:), 1.E-20)
ZOMEGA_UP_M(:)=(XG*ZOMEGA_UP_M(:))**(1./3.)

!computation of alpha up
DO JK=KKA,KKU,KKL
  ZALPHA_UP_M(:,JK)=ZEMF_M(:,JK)/(XALPHA_MF*PRHODREF(:,JK)*ZOMEGA_UP_M(:))
ENDDO
ZALPHA_UP_M(:,:)=MAX(0., MIN(ZALPHA_UP_M(:,:), 1.))

!computation of sigma of the distribution
DO JK=KKA,KKU,KKL
  ZSIGMF(:,JK)=ZEMF_M(:,JK) * &
               (ZRT_UP_M(:,JK) - PRTM(:,JK)) * &
               PDEPTH(:) * ZGRAD_Z_RT(:,JK) / &
               (XSIGMA_MF * ZOMEGA_UP_M(:) * PRHODREF(:,JK))
ENDDO
ZSIGMF(:,:)=SQRT(MAX(ABS(ZSIGMF(:,:)), 1.E-40))
!
!*      2. PDF integration
!          ------------------------------------------------
!
!The mean of the distribution is ZRT_UP
!Computation of ZA and ZGAM (=efrc(ZA)) coefficient
ZA(:,:)=(ZRSAT_UP_M(:,:)-ZRT_UP_M(:,:))/(sqrt(2.)*ZSIGMF(:,:))

!erf computed by an incomplete gamma function approximation
!DO JK=KKA,KKU,KKL
!  DO JI=1, SIZE(PCF_MF,1)
!    IF(ZA(JI,JK)>1E-20) THEN
!      ZGAM(JI,JK)=1-GAMMA_INC(0.5,ZA(JI,JK)**2)
!    ELSEIF(ZA(JI,JK)<-1E-20) THEN
!      ZGAM(JI,JK)=1+GAMMA_INC(0.5,ZA(JI,JK)**2)
!    ELSE
!      ZGAM(JI,JK)=1
!    ENDIF
!  ENDDO
!ENDDO

!alternative approximation of erf function (better for vectorisation)
WHERE(ZA(:,:)>0)
  ZGAM(:,:)=1-SQRT(1-EXP(-4*ZA(:,:)**2/XPI))
ELSEWHERE
  ZGAM(:,:)=1+SQRT(1-EXP(-4*ZA(:,:)**2/XPI))
ENDWHERE

!computation of cloud fraction
PCF_MF(:,:)=MAX( 0., MIN(1.,0.5*ZGAM(:,:) * ZALPHA_UP_M(:,:)))

!computation of condensate, then PRC and PRI
ZCOND(:,:)=(EXP(-ZA(:,:)**2)-ZA(:,:)*SQRT(XPI)*ZGAM(:,:))*ZSIGMF(:,:)/SQRT(2.*XPI) * ZALPHA_UP_M(:,:)
ZCOND(:,:)=MAX(ZCOND(:,:), 0.) !due to approximation of ZGAM value, ZCOND could be slightly negative
PRC_MF(:,:)=(1.-ZFRAC_ICE_UP_M(:,:)) * ZCOND(:,:)
PRI_MF(:,:)=(   ZFRAC_ICE_UP_M(:,:)) * ZCOND(:,:)

IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_BIGAUS',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_MF_CLOUD_BIGAUS
