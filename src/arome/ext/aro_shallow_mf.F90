!     ######spl
      SUBROUTINE  ARO_SHALLOW_MF(PHYEX,&
                KKL, KLON, KLEV, KFDIA, KRR, KRRL, KRRI,KSV,&
                KSV_LGBEG,KSV_LGEND, &
                PTSTEP, PDX, PDY,                                     &
                PZZ, PZZF, PDZZF,                                     &
                PRHODJ, PRHODREF,                                     &
                PPABSM, PEXNM,                                        &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,                                             &
                PUM,PVM,PTKEM,PSVM,                                   &
                PDUDT_MF,PDVDT_MF,PDTKEDT_MF,                         &
                PDTHLDT_MF,PDRTDT_MF,PDSVDT_MF,                       &
                PSIGMF,PRC_MF,PRI_MF,PCF_MF,                          &
                PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF,   &
                PWEIGHT_MF_CLOUD, PFLXZTHVMF, PFLXZUMF, PFLXZVMF,     &
                PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,                  &
                PU_UP, PV_UP, PTHV_UP, PW_UP, PFRAC_UP, PEMF,         &
                YDDDH,YDLDDH,YDMDDH                                   )

      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ##########################################################################
!
!!****  * -  interface to call SHALLOW_MF :
!!             computation of turbulence "mass flux" fluxes and their divergence
!!
!!
!!
!!    PURPOSE
!!    -------
!!
!!
!!
!!
!!**  METHOD
!!    ------
!!
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine SHALLOW_MF (routine de MesoNH)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!      Documentation AROME
!!
!!    AUTHOR
!!    ------
!!    S.Malardel
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    07/2006
!!      Y. Seity : new arguments for EDMF scheme 04/2009
!!      S. Riette 18 May 2010: aro_shallow_mf and shallow_mf interfaces changed
!!      S. Riette Jan 2012: support for both order of vertical levels
!!      A. Marcel Jan 2025: EDMF contribution to dynamic TKE production
!!      A. Marcel Jan 2025: TKE mixing
!!      A. Marcel Jan 2025: bi-Gaussian PDF and associated subgrid precipitation
!!      A. Marcel Jan 2025: relaxation of the small fraction assumption
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_BUDGET, ONLY: NBUDGET_SV1, TBUDGETDATA
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODI_SHALLOW_MF
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
USE YOMLSFORC, ONLY: LMUSCLFA, NMUSCLFA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
TYPE(PHYEX_t),            INTENT(IN)   :: PHYEX
INTEGER,                  INTENT(IN)   :: KKL      ! +1 if grid goes from ground to
                                                   ! atmosphere top, -1 otherwise
INTEGER,                  INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KFDIA
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KRRL     ! Number of liquide water variables
INTEGER,                  INTENT(IN)   :: KRRI     ! Number of ice variables
INTEGER,                  INTENT(IN)   :: KSV      ! Number of passive scalar variables
!
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer

REAL,                     INTENT(IN)   :: PTSTEP   ! Time step
REAL,                     INTENT(IN)   :: PDX      ! grid size along x-axis
REAL,                     INTENT(IN)   :: PDY      ! grid size along y-axis
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PZZ     ! Height of layer boundaries
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PZZF    ! Height of level
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PDZZF    !thikness between layers

REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PRHODJ  ! Dry density * Jacobian
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PRHODREF  ! Dry density
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PPABSM      ! Pressure at time t-1
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PEXNM   ! Exner function
!
! normal surface fluxes of theta and Rv
REAL, DIMENSION(KLON),          INTENT(IN)   ::  PSFTH,PSFRV
!   prognostic variables at t- deltat
!
!   thermodynamical variables which are transformed in conservative var.
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   ::  PTHM       ! pot. temp.
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(IN) ::  PRM         ! mixing ratio
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   ::  PUM,PVM       ! momentum
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   ::  PTKEM
REAL, DIMENSION(KLON,KLEV,KSV), INTENT(IN) ::  PSVM         ! passive scalar
                                                             ! variables for EDMF scheme
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::  PDUDT_MF     ! tendency of U   by massflux scheme
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::  PDVDT_MF     ! tendency of V   by massflux scheme
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::  PDTKEDT_MF   ! tendency of TKE by massflux scheme
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::  PDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::  PDRTDT_MF    ! tendency of rt  by massflux scheme
REAL, DIMENSION(KLON,KLEV,KSV), INTENT(OUT)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme

REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   ::  PHLC_HRC_MF, PHLC_HCF_MF, PHLI_HRI_MF, PHLI_HCF_MF ! high/low cloud diagnostics
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   ::  PWEIGHT_MF_CLOUD ! weight coefficient for the mass-flux cloud
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   ::  PFLXZTHVMF           ! Thermal production for TKE scheme
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   ::  PFLXZUMF             ! Dynamic production for TKE scheme
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   ::  PFLXZVMF             ! Dynamic production for TKE scheme
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PTHL_UP   ! Thl updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PRT_UP    ! Rt  updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PRV_UP    ! Vapor updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PU_UP     ! U wind updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PV_UP     ! V wind updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PRC_UP    ! cloud content updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PRI_UP    ! ice content   updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PTHV_UP   ! Thv   updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PW_UP     ! vertical speed updraft characteristics
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PFRAC_UP  ! updraft fraction
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) ::  PEMF      ! updraft mass flux
!
TYPE(TYP_DDH), INTENT(INOUT), TARGET   :: YDDDH
TYPE(TLDDH),   INTENT(IN), TARGET      :: YDLDDH
TYPE(TMDDH),   INTENT(IN), TARGET      :: YDMDDH
!
!
!*       0.2   Declarations of local variables :
!
TYPE(TBUDGETDATA), DIMENSION(NBUDGET_SV1) :: YLBUDGET !NBUDGET_SV1 is the one with the highest number needed for shallow_mf
INTEGER, DIMENSION(KLON) :: IKLCL,IKETL,IKCTL
REAL,DIMENSION(KLON,KLEV) :: ZFLXZTHMF,ZFLXZRMF,ZFLXZTKEMF
REAL,DIMENSION(KLON,KLEV) :: ZDETR,ZENTR
REAL,DIMENSION(KLON,KLEV) :: ZTKE_UP
REAL,DIMENSION(size(PRHODJ,1),size(PRHODJ,2)) :: ZW
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JBU ! Loop index for budgets
!
#include "wrarom.intfb.h"
!
!
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_SHALLOW_MF',0,ZHOOK_HANDLE)

!Dimensions
CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, JPVEXT, KFDIA)
!
!
!------------------------------------------------------------------------------
!
!*       2.   INITIALISATION
!
!             ---------------


! tableau a recalculer a chaque pas de temps
! attention, ZDZZ est l'altitude entre deux niveaux (et pas l'�paisseur de la couche)

!DO JL = IIB,IIE
!   DO JK = 2, SIZE(PZZF,2)-1
!      ZDZZ(JL,JK)=PZZF(JL,JK)-PZZF(JL,JK-KKL)
!   ENDDO
!   ZDZZ(JL,IKA)=PZZF(JL,IKA)-(1.5*PZZ(JL,IKA)-0.5*PZZ(JL,IKA+KKL)) ! must work with JPVEXT=0 or 1
!   ZDZZ(JL,IKU)=PZZF(JL,IKU)-PZZF(JL,IKU-KKL) ! excluded from the loop because depending on KKL, IKU can be 1 or SIZE()
!ENDDO
!
!
!------------------------------------------------------------------------------
!
!
! Budgets
DO JBU=1, NBUDGET_SV1
  YLBUDGET(JBU)%NBUDGET=JBU
  YLBUDGET(JBU)%YDDDH=>YDDDH
  YLBUDGET(JBU)%YDLDDH=>YDLDDH
  YLBUDGET(JBU)%YDMDDH=>YDMDDH
ENDDO
!
!------------------------------------------------------------------------------
!
!
!*       4.   APPEL DE LA CONVECTION PEU PROFONDE MESONH
!
!         ---------------------------------
!
  CALL SHALLOW_MF(YLDIMPHYEX, PHYEX%CST, PHYEX%NEBN, PHYEX%PARAM_MFSHALLN, PHYEX%TURBN, PHYEX%CSTURB, PHYEX%RAIN_ICE_PARAMN, &
     &KRR=KRR, KRRL=KRRL, KRRI=KRRI, KSV=KSV,                                             &
     &ONOMIXLG=PHYEX%MISC%ONOMIXLG,KSV_LGBEG=KSV_LGBEG,KSV_LGEND=KSV_LGEND, &
     &PTSTEP=PTSTEP,                                                                      &
     &PDZZ=PDZZF,PZZ=PZZ,                                                                 &
     &PRHODJ=PRHODJ,PRHODREF=PRHODREF,                                                    &
     &PPABSM=PPABSM,PEXNM=PEXNM,                                                          &
     &PSFTH=PSFTH,PSFRV=PSFRV,                                                            &
     &PTHM=PTHM,PRM=PRM,PUM=PUM,PVM=PVM,PTKEM=PTKEM,PSVM=PSVM,                            &
     &PDUDT_MF=PDUDT_MF,PDVDT_MF=PDVDT_MF,PDTKEDT_MF=PDTKEDT_MF,                          &
     &PDTHLDT_MF=PDTHLDT_MF,PDRTDT_MF=PDRTDT_MF,PDSVDT_MF=PDSVDT_MF,                      &
     &PSIGMF=PSIGMF,PRC_MF=PRC_MF,PRI_MF=PRI_MF,PCF_MF=PCF_MF,                            &
     &PHLC_HRC=PHLC_HRC_MF, PHLC_HCF=PHLC_HCF_MF, PHLI_HRI=PHLI_HRI_MF, PHLI_HCF=PHLI_HCF_MF,&
     &PWEIGHT_MF_CLOUD=PWEIGHT_MF_CLOUD, PFLXZTHVMF=PFLXZTHVMF,      &
     &PFLXZTHMF=ZFLXZTHMF,PFLXZRMF=ZFLXZRMF,PFLXZUMF=PFLXZUMF,PFLXZVMF=PFLXZVMF,PFLXZTKEMF=ZFLXZTKEMF, &
     &PTHL_UP=PTHL_UP,PRT_UP=PRT_UP,PRV_UP=PRV_UP,PRC_UP=PRC_UP,PRI_UP=PRI_UP,            &
     &PU_UP=PU_UP, PV_UP=PV_UP, PTKE_UP=ZTKE_UP, PTHV_UP=PTHV_UP, PW_UP=PW_UP,            &
     &PFRAC_UP=PFRAC_UP,PEMF=PEMF,PDETR=ZDETR,PENTR=ZENTR,                                &
     &KKLCL=IKLCL,KKETL=IKETL,KKCTL=IKCTL,PDX=PDX,PDY=PDY,                                &
     &BUCONF=PHYEX%MISC%TBUCONF, TBUDGETS=YLBUDGET, KBUDGETS=SIZE(YLBUDGET)               )
!
!
!------------------------------------------------------------------------------
!
!
!*       5.   DIVISION PAR RHODJ DES TERMES SOURCES DE MESONH
!             (on obtient des termes homog�nes � des tendances)
!
!             -----------------------------------------------
!
IF (LMUSCLFA) THEN
  !Some updraft variables are initialised with environmental values.
  !ZW mask must be used to suppress those values, for variables on flux levels
  ZW=1.
  WHERE(PFRAC_UP==0.) ZW=0.

  !Variables on flux levels
  CALL WRAROM(NMUSCLFA, 'MF',         PEMF*ZW,               &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'CF_MF',      PCF_MF,                &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .FALSE.)
  CALL WRAROM(NMUSCLFA, 'ZWU',        PW_UP*ZW,              &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'PUDAL',      PFRAC_UP*ZW,           &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'ZEPSILON_U', ZENTR,                 &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .FALSE.)
  CALL WRAROM(NMUSCLFA, 'ZDELTA_U',   ZDETR,                 &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .FALSE.)
  CALL WRAROM(NMUSCLFA, 'PQLSHCONV',  PRC_UP/(1.+PRT_UP)*ZW, &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'PQISHCONV',  PRI_UP/(1.+PRT_UP)*ZW, &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'PQU',        PRV_UP/(1.+PRT_UP)*ZW, &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'THV_UP',     PTHV_UP*ZW,            &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'RT_UP',      PRT_UP*ZW,             &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'TKE_UP',     ZTKE_UP*ZW,            &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'WE_UP',      ZFLXZTKEMF*ZW,         &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'WRT_UP',     ZFLXZRMF*ZW,           &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
  CALL WRAROM(NMUSCLFA, 'WTHL_UP',    ZFLXZTHMF*ZW,          &
              KLON, YLDIMPHYEX%NKT, YLDIMPHYEX%NKB, YLDIMPHYEX%NKE, YLDIMPHYEX%NKL, .TRUE.)
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ARO_SHALLOW_MF',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_SHALLOW_MF
