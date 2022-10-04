!     ######spl
      SUBROUTINE  ARO_SHALLOW_MF(KKL, KLON, KLEV, KFDIA, KRR, KRRL, KRRI,KSV,     &
                HMF_UPDRAFT, HMF_CLOUD, HFRAC_ICE, OMIXUV,            &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                KTCOUNT, PTSTEP, PDX, PDY,                            &
                PZZ, PZZF, PDZZF,                                            &
                PRHODJ, PRHODREF,                                     &
                PPABSM, PEXNM,                                        &
                PSFTH,PSFRV,                                          &
                PTHM,PRM,                                             &
                PUM,PVM,PTKEM,PSVM,                                   &
                PDUDT_MF,PDVDT_MF,                                    &
                PDTHLDT_MF,PDRTDT_MF,PDSVDT_MF,                       &
                PSIGMF,PRC_MF,PRI_MF,PCF_MF,PFLXZTHVMF,                      &
                PTHL_UP,PRT_UP,PRV_UP,PRC_UP,PRI_UP,                  &
                PU_UP, PV_UP, PTHV_UP, PW_UP, PFRAC_UP, PEMF)

      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPVEXT
USE MODD_CST, ONLY: CST
USE MODD_NEB, ONLY: NEB
USE MODD_TURB_n, ONLY: TURBN
USE MODD_CTURB, ONLY: CSTURB
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALLN
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODI_SHALLOW_MF
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!
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
CHARACTER (LEN=4), INTENT(IN)   :: HMF_UPDRAFT  ! Type of Mass Flux Scheme
CHARACTER (LEN=4), INTENT(IN)   :: HMF_CLOUD    ! Type of statistical cloud scheme
CHARACTER*1,       INTENT(IN)   :: HFRAC_ICE    ! partition liquid/ice scheme
LOGICAL,                        INTENT(IN)   :: OMIXUV    ! True if mixing of momentum
!
LOGICAL,                INTENT(IN)   :: ONOMIXLG  ! False if mixing of lagrangian tracer
INTEGER,                INTENT(IN)   :: KSV_LGBEG ! first index of lag. tracer
INTEGER,                INTENT(IN)   :: KSV_LGEND ! last  index of lag. tracer

INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
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
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::  PDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(KLON,KLEV),   INTENT(OUT)::  PDRTDT_MF    ! tendency of rt  by massflux scheme
REAL, DIMENSION(KLON,KLEV,KSV), INTENT(OUT)::  PDSVDT_MF    ! tendency of Sv  by massflux scheme

REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   ::  PSIGMF,PRC_MF,PRI_MF,PCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(KLON,KLEV), INTENT(OUT)   ::  PFLXZTHVMF           ! Thermal production for TKE scheme
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
!
!*       0.2   Declarations of local variables :
!
INTEGER, DIMENSION(size(PRHODJ,1)) :: IKLCL,IKETL,IKCTL
REAL,DIMENSION(size(PRHODJ,1),size(PRHODJ,2)) :: ZFLXZTHMF,ZFLXZRMF,ZFLXZUMF,ZFLXZVMF
REAL,DIMENSION(size(PRHODJ,1),size(PRHODJ,2)) :: ZDETR,ZENTR
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
REAL          ::  ZIMPL        ! degree of implicitness
REAL(KIND=JPRB) :: ZHOOK_HANDLE

! vars for incomplete NPROMA test
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PDZZF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PZZ
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PRHODJ
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PRHODREF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PPABSM
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PEXNM
REAL, DIMENSION(KLON*2) :: ZP_PSFTH
REAL, DIMENSION(KLON*2) :: ZP_PSFRV
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PTHM
REAL, DIMENSION(KLON*2, KLEV, KRR) :: ZP_PRM
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PUM
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PVM
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PTKEM
REAL, DIMENSION(KLON*2, KLEV, KSV) :: ZP_PSVM
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PDUDT_MF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PDVDT_MF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PDTHLDT_MF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PDRTDT_MF
REAL, DIMENSION(KLON*2, KLEV, KSV) :: ZP_PDSVDT_MF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PSIGMF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PRC_MF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PRI_MF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PCF_MF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PFLXZTHVMF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_ZFLXZTHMF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_ZFLXZRMF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_ZFLXZUMF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_ZFLXZVMF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PTHL_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PRT_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PRV_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PRC_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PRI_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PU_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PV_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PTHV_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PW_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PFRAC_UP
REAL, DIMENSION(KLON*2, KLEV) :: ZP_PEMF
REAL, DIMENSION(KLON*2, KLEV) :: ZP_ZDETR
REAL, DIMENSION(KLON*2, KLEV) :: ZP_ZENTR
INTEGER, DIMENSION(KLON*2, KLEV) :: IP_IKLCL
INTEGER, DIMENSION(KLON*2, KLEV) :: IP_IKETL
INTEGER, DIMENSION(KLON*2, KLEV) :: IP_IKCTL

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
!CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, JPVEXT, KFDIA)
CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON*2, 1, KLEV, JPVEXT, KFDIA)
!
!
!------------------------------------------------------------------------------
!
!*       2.   INITIALISATION
!
!             ---------------


ZIMPL=1.
!ZIMPL=0.
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
!*       3.   MULTIPLICATION PAR RHODJ
!             POUR OBTENIR LES TERMES SOURCES DE MESONH
!
!         -----------------------------------------------

!
!------------------------------------------------------------------------------
!
!
!*       4.   APPEL DE LA TURBULENCE MESONH
!
!         ---------------------------------
!
ZP_PDZZF(:,:)=HUGE(PDZZF)
ZP_PZZ(:,:)=HUGE(PDZZF)
ZP_PRHODJ(:,:)=HUGE(PDZZF)
ZP_PRHODREF(:,:)=HUGE(PDZZF)
ZP_PPABSM(:,:)=HUGE(PDZZF)
ZP_PEXNM(:,:)=HUGE(PDZZF)
ZP_PSFTH(:)=HUGE(PDZZF)
ZP_PSFRV(:)=HUGE(PDZZF)
ZP_PTHM(:,:)=HUGE(PDZZF)
ZP_PRM(:,:,:)=HUGE(PDZZF)
ZP_PUM(:,:)=HUGE(PDZZF)
ZP_PVM(:,:)=HUGE(PDZZF)
ZP_PTKEM(:,:)=HUGE(PDZZF)
ZP_PSVM(:,:,:)=HUGE(PDZZF)
ZP_PDUDT_MF(:,:)=HUGE(PDZZF)
ZP_PDVDT_MF(:,:)=HUGE(PDZZF)
ZP_PDTHLDT_MF(:,:)=HUGE(PDZZF)
ZP_PDRTDT_MF(:,:)=HUGE(PDZZF)
ZP_PDSVDT_MF(:,:,:)=HUGE(PDZZF)
ZP_PSIGMF(:,:)=HUGE(PDZZF)
ZP_PRC_MF(:,:)=HUGE(PDZZF)
ZP_PRI_MF(:,:)=HUGE(PDZZF)
ZP_PCF_MF(:,:)=HUGE(PDZZF)
ZP_PFLXZTHVMF(:,:)=HUGE(PDZZF)
ZP_ZFLXZTHMF(:,:)=HUGE(PDZZF)
ZP_ZFLXZRMF(:,:)=HUGE(PDZZF)
ZP_ZFLXZUMF(:,:)=HUGE(PDZZF)
ZP_ZFLXZVMF(:,:)=HUGE(PDZZF)
ZP_PTHL_UP(:,:)=HUGE(PDZZF)
ZP_PRT_UP(:,:)=HUGE(PDZZF)
ZP_PRV_UP(:,:)=HUGE(PDZZF)
ZP_PRC_UP(:,:)=HUGE(PDZZF)
ZP_PRI_UP(:,:)=HUGE(PDZZF)
ZP_PU_UP(:,:)=HUGE(PDZZF)
ZP_PV_UP(:,:)=HUGE(PDZZF)
ZP_PTHV_UP(:,:)=HUGE(PDZZF)
ZP_PW_UP(:,:)=HUGE(PDZZF)
ZP_PFRAC_UP(:,:)=HUGE(PDZZF)
ZP_PEMF(:,:)=HUGE(PDZZF)
ZP_ZDETR(:,:)=HUGE(PDZZF)
ZP_ZENTR(:,:)=HUGE(PDZZF)
IP_IKLCL(:,:)=HUGE(IKLCL)
IP_IKETL(:,:)=HUGE(IKLCL)
IP_IKCTL(:,:)=HUGE(IKLCL)

ZP_PDZZF(:KLON,:)=PDZZF
ZP_PZZ(:KLON,:)=PZZ
ZP_PRHODJ(:KLON,:)=PRHODJ
ZP_PRHODREF(:KLON,:)=PRHODREF
ZP_PPABSM(:KLON,:)=PPABSM
ZP_PEXNM(:KLON,:)=PEXNM
ZP_PSFTH(:KLON)=PSFTH
ZP_PSFRV(:KLON)=PSFRV
ZP_PTHM(:KLON,:)=PTHM
ZP_PRM(:KLON,:,:)=PRM
ZP_PUM(:KLON,:)=PUM
ZP_PVM(:KLON,:)=PVM
ZP_PTKEM(:KLON,:)=PTKEM
ZP_PSVM(:KLON,:,:)=PSVM
!ZP_PDUDT_MF(:KLON,:)=PDUDT_MF
!ZP_PDVDT_MF(:KLON,:)=PDVDT_MF
!ZP_PDTHLDT_MF(:KLON,:)=PDTHLDT_MF
!ZP_PDRTDT_MF(:KLON,:)=PDRTDT_MF
!ZP_PDSVDT_MF(:KLON,:,:)=PDSVDT_MF
!ZP_PSIGMF(:KLON,:)=PSIGMF
!ZP_PRC_MF(:KLON,:)=PRC_MF
!ZP_PRI_MF(:KLON,:)=PRI_MF
!ZP_PCF_MF(:KLON,:)=PCF_MF
!ZP_PFLXZTHVMF(:KLON,:)=PFLXZTHVMF
!ZP_ZFLXZTHMF(:KLON,:)=ZFLXZTHMF
!ZP_ZFLXZRMF(:KLON,:)=ZFLXZRMF
!ZP_ZFLXZUMF(:KLON,:)=ZFLXZUMF
!ZP_ZFLXZVMF(:KLON,:)=ZFLXZVMF
ZP_PTHL_UP(:KLON,:)=PTHL_UP
ZP_PRT_UP(:KLON,:)=PRT_UP
!ZP_PRV_UP(:KLON,:)=PRV_UP
ZP_PRC_UP(:KLON,:)=PRC_UP
ZP_PRI_UP(:KLON,:)=PRI_UP
ZP_PU_UP(:KLON,:)=PU_UP
ZP_PV_UP(:KLON,:)=PV_UP
ZP_PTHV_UP(:KLON,:)=PTHV_UP
ZP_PW_UP(:KLON,:)=PW_UP
ZP_PFRAC_UP(:KLON,:)=PFRAC_UP
ZP_PEMF(:KLON,:)=PEMF
!ZP_ZDETR(:KLON,:)=ZDETR
!ZP_ZENTR(:KLON,:)=ZENTR
!IP_IKLCL(:KLON,:)=IKLCL
!IP_IKETL(:KLON,:)=IKETL
!IP_IKCTL(:KLON,:)=IKCTL

  CALL SHALLOW_MF(YLDIMPHYEX, CST, NEB, PARAM_MFSHALLN, TURBN, CSTURB,                    &
     &KRR=KRR, KRRL=KRRL, KRRI=KRRI, KSV=KSV,                                             &
     &HMF_UPDRAFT=HMF_UPDRAFT, HMF_CLOUD=HMF_CLOUD,HFRAC_ICE=HFRAC_ICE,OMIXUV=OMIXUV,     &
     &ONOMIXLG=ONOMIXLG,KSV_LGBEG=KSV_LGBEG,KSV_LGEND=KSV_LGEND,                          &
     &PIMPL_MF=ZIMPL, PTSTEP=PTSTEP,                                                      &
     &PDZZ=ZP_PDZZF,PZZ=ZP_PZZ,                                                                 &
     &PRHODJ=ZP_PRHODJ,PRHODREF=ZP_PRHODREF,                                                    &
     &PPABSM=ZP_PPABSM,PEXNM=ZP_PEXNM,                                                          &
     &PSFTH=ZP_PSFTH,PSFRV=ZP_PSFRV,                                                            &
     &PTHM=ZP_PTHM,PRM=ZP_PRM,PUM=ZP_PUM,PVM=ZP_PVM,PTKEM=ZP_PTKEM,PSVM=ZP_PSVM,                            &
     &PDUDT_MF=ZP_PDUDT_MF,PDVDT_MF=ZP_PDVDT_MF,                                                &
     &PDTHLDT_MF=ZP_PDTHLDT_MF,PDRTDT_MF=ZP_PDRTDT_MF,PDSVDT_MF=ZP_PDSVDT_MF,                      &
     &PSIGMF=ZP_PSIGMF,PRC_MF=ZP_PRC_MF,PRI_MF=ZP_PRI_MF,PCF_MF=ZP_PCF_MF,PFLXZTHVMF=ZP_PFLXZTHVMF,      &
     &PFLXZTHMF=ZP_ZFLXZTHMF,PFLXZRMF=ZP_ZFLXZRMF,PFLXZUMF=ZP_ZFLXZUMF,PFLXZVMF=ZP_ZFLXZVMF,          &
     &PTHL_UP=ZP_PTHL_UP,PRT_UP=ZP_PRT_UP,PRV_UP=ZP_PRV_UP,PRC_UP=ZP_PRC_UP,PRI_UP=ZP_PRI_UP,            &
     &PU_UP=ZP_PU_UP, PV_UP=ZP_PV_UP, PTHV_UP=ZP_PTHV_UP, PW_UP=ZP_PW_UP,                             &
     &PFRAC_UP=ZP_PFRAC_UP,PEMF=ZP_PEMF,PDETR=ZP_ZDETR,PENTR=ZP_ZENTR,                                &
     &KKLCL=IP_IKLCL,KKETL=IP_IKETL,KKCTL=IP_IKCTL,PDX=PDX,PDY=PDY                                 )
!PDZZF=ZP_PDZZF(:KLON,:)
!PZZ=ZP_PZZ(:KLON,:)
!PRHODJ=ZP_PRHODJ(:KLON,:)
!PRHODREF=ZP_PRHODREF(:KLON,:)
!PPABSM=ZP_PPABSM(:KLON,:)
!PEXNM=ZP_PEXNM(:KLON,:)
!PSFTH=ZP_PSFTH(:KLON)
!PSFRV=ZP_PSFRV(:KLON)
!PTHM=ZP_PTHM(:KLON,:)
!PRM=ZP_PRM(:KLON,:,:)
!PUM=ZP_PUM(:KLON,:)
!PVM=ZP_PVM(:KLON,:)
!PTKEM=ZP_PTKEM(:KLON,:)
!PSVM=ZP_PSVM(:KLON,:,:)
PDUDT_MF=ZP_PDUDT_MF(:KLON,:)
PDVDT_MF=ZP_PDVDT_MF(:KLON,:)
PDTHLDT_MF=ZP_PDTHLDT_MF(:KLON,:)
PDRTDT_MF=ZP_PDRTDT_MF(:KLON,:)
PDSVDT_MF=ZP_PDSVDT_MF(:KLON,:,:)
PSIGMF=ZP_PSIGMF(:KLON,:)
PRC_MF=ZP_PRC_MF(:KLON,:)
PRI_MF=ZP_PRI_MF(:KLON,:)
PCF_MF=ZP_PCF_MF(:KLON,:)
PFLXZTHVMF=ZP_PFLXZTHVMF(:KLON,:)
!ZFLXZTHMF=ZP_ZFLXZTHMF(:KLON,:)
!ZFLXZRMF=ZP_ZFLXZRMF(:KLON,:)
!ZFLXZUMF=ZP_ZFLXZUMF(:KLON,:)
!ZFLXZVMF=ZP_ZFLXZVMF(:KLON,:)
PTHL_UP=ZP_PTHL_UP(:KLON,:)
PRT_UP=ZP_PRT_UP(:KLON,:)
PRV_UP=ZP_PRV_UP(:KLON,:)
PRC_UP=ZP_PRC_UP(:KLON,:)
PRI_UP=ZP_PRI_UP(:KLON,:)
PU_UP=ZP_PU_UP(:KLON,:)
PV_UP=ZP_PV_UP(:KLON,:)
PTHV_UP=ZP_PTHV_UP(:KLON,:)
PW_UP=ZP_PW_UP(:KLON,:)
PFRAC_UP=ZP_PFRAC_UP(:KLON,:)
PEMF=ZP_PEMF(:KLON,:)
!ZDETR=ZP_ZDETR(:KLON,:)
!ZENTR=ZP_ZENTR(:KLON,:)
!IKLCL=IP_IKLCL(:KLON,:)
!IKETL=IP_IKETL(:KLON,:)
!IKCTL=IP_IKCTL(:KLON,:)

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
IF (LHOOK) CALL DR_HOOK('ARO_SHALLOW_MF',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_SHALLOW_MF
