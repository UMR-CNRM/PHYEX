!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODE_COMPUTE_MF_CLOUD_DIRECT
!    ###################################
!
IMPLICIT NONE
CONTAINS
      SUBROUTINE COMPUTE_MF_CLOUD_DIRECT(D, PARAMMF, &
                                        &KKLCL, PFRAC_UP, PRC_UP, PRI_UP,&
                                        &PRC_MF, PRI_MF, PCF_MF, PWEIGHT_MF_CLOUD)
!     #################################################################
!!
!!****  *COMPUTE_MF_CLOUD_DIRECT* -
!!       compute diagnostic subgrid cumulus cloud caracteristics with a direct scheme
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
!!      Updraft variables are used directly to diagnose subgrid clouds
!!      This scheme may be activated only if the selected updraft model
!!      gives the updraft fraction as an output
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
!!      S. Riette Apr 2013: computation begins one level lower (to be able to have a cloud
!!                          on mass level just below the first saturated flux level)
!!      A. Marcel Jan 2025: relaxation of the small fraction assumption
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_PARAM_MFSHALL_n, ONLY : PARAM_MFSHALL_t
USE MODI_SHUMAN_MF, ONLY: MZF_MF
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*                    0.1  Declaration of Arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
TYPE(PARAM_MFSHALL_t),  INTENT(IN)   :: PARAMMF
INTEGER, DIMENSION(D%NIJT),  INTENT(IN)   :: KKLCL          ! index of updraft condensation level
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PFRAC_UP       ! Updraft Fraction
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(IN)   :: PRC_UP,PRI_UP  ! updraft characteritics
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PRC_MF, PRI_MF ! cloud content
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PCF_MF         ! and cloud fraction for MF scheme
REAL, DIMENSION(D%NIJT,D%NKT),   INTENT(OUT)  :: PWEIGHT_MF_CLOUD ! weight coefficient for the mass-flux cloud
!
!*                    0.1  Declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER  :: JIJ,JK, JK0, IKB,IKE,IKL,IKT,IIJB,IIJE
REAL, DIMENSION(D%NIJT,D%NKT) :: ZFRAC_UP_M
!
!*                    0.2 Initialisation
!
IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_DIRECT',0,ZHOOK_HANDLE)
!
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
IKT=D%NKT
IIJB=D%NIJB
IIJE=D%NIJE
!*      1. COMPUTATION OF SUBGRID CLOUD
!          ----------------------------

!
! Warning: updraft variables are on flux levels
! and PRC_MF, PRI_MF and PCF_MF are on mass levels
PRC_MF(:,:)=0.
PRI_MF(:,:)=0.
PCF_MF(:,:)=0.

IF(PARAMMF%LRELAX_ALPHA_MF) THEN
  CALL MZF_MF(D, PFRAC_UP(:,:), ZFRAC_UP_M(:,:))
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PWEIGHT_MF_CLOUD(JIJ, JK) = ZFRAC_UP_M(JIJ, JK)*PARAMMF%XALPHA_MF
      PWEIGHT_MF_CLOUD(JIJ, JK) = MAX(0., MIN(PWEIGHT_MF_CLOUD(JIJ, JK), 1.))
    END DO
  END DO
ENDIF

DO JIJ=IIJB,IIJE
  JK0=KKLCL(JIJ)-IKL ! first mass level with cloud
  JK0=MAX(JK0, MIN(IKB,IKE)) !protection if KKL=1
  JK0=MIN(JK0, MAX(IKB,IKE)) !protection if KKL=-1
  DO JK=JK0,IKE-IKL,IKL
    PCF_MF(JIJ,JK ) = MAX( 0., MIN(1.,PARAMMF%XKCF_MF *0.5* (       &
                &    PFRAC_UP(JIJ,JK) +  PFRAC_UP(JIJ,JK+IKL) ) ))
    PRC_MF(JIJ,JK)  = 0.5* PARAMMF%XKCF_MF * ( PFRAC_UP(JIJ,JK)*PRC_UP(JIJ,JK)  &
                         + PFRAC_UP(JIJ,JK+IKL)*PRC_UP(JIJ,JK+IKL) )
    PRI_MF(JIJ,JK)  = 0.5* PARAMMF%XKCF_MF * ( PFRAC_UP(JIJ,JK)*PRI_UP(JIJ,JK)  &
                         + PFRAC_UP(JIJ,JK+IKL)*PRI_UP(JIJ,JK+IKL) )
  END DO
END DO

IF (LHOOK) CALL DR_HOOK('COMPUTE_MF_CLOUD_DIRECT',1,ZHOOK_HANDLE)

END SUBROUTINE COMPUTE_MF_CLOUD_DIRECT
END MODULE MODE_COMPUTE_MF_CLOUD_DIRECT
