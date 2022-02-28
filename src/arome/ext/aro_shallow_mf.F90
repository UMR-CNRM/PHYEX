!     ######spl
      SUBROUTINE  ARO_SHALLOW_MF(KKL, KLON,KLEV, KRR, KRRL, KRRI,KSV,     &
                HMF_UPDRAFT, HMF_CLOUD, HFRAC_ICE, OMIXUV,            &
                ONOMIXLG,KSV_LGBEG,KSV_LGEND,                         &
                KTCOUNT, PTSTEP,                                      &
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
USE MODD_PARAMETERS, ONLY: JPVEXT, JPHEXT
!
USE MODI_SHALLOW_MF
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
INTEGER :: JRR           ! Loop index for the moist
INTEGER :: IIB           ! Define the physical domain
INTEGER :: IIE           !
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKB           !
INTEGER :: IKE           !
INTEGER :: IKA, IKU
INTEGER :: JI, JJ, JL, JK !
INTEGER ::II
INTEGER, DIMENSION(size(PRHODJ,1)) :: IKLCL,IKETL,IKCTL
REAL,DIMENSION(size(PRHODJ,1),size(PRHODJ,2)) :: ZFLXZTHMF,ZFLXZRMF,ZFLXZUMF,ZFLXZVMF
REAL,DIMENSION(size(PRHODJ,1),size(PRHODJ,2)) :: ZDETR,ZENTR
!
!

REAL          ::  ZIMPL        ! degree of implicitness
!
!
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_SHALLOW_MF',0,ZHOOK_HANDLE)


IIB=1+JPHEXT
IIE=SIZE(PZZ,1) - JPHEXT
IJB=1+JPHEXT
IJE=1 - JPHEXT
IF(KKL==1)THEN
  IKA=1
  IKU=SIZE(PZZ,2)
ELSE
  IKA=SIZE(PZZ,2)
  IKU=1
ENDIF
IKB=IKA+KKL*JPVEXT
IKE=IKU-KKL*JPVEXT
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
  CALL SHALLOW_MF(KKA=IKA,KKU=IKU,KKL=KKL,KRR=KRR,KRRL=KRRL,KRRI=KRRI,                    &
     &HMF_UPDRAFT=HMF_UPDRAFT, HMF_CLOUD=HMF_CLOUD,HFRAC_ICE=HFRAC_ICE,OMIXUV=OMIXUV,     &
     &ONOMIXLG=ONOMIXLG,KSV_LGBEG=KSV_LGBEG,KSV_LGEND=KSV_LGEND,                          &
     &PIMPL_MF=ZIMPL, PTSTEP=PTSTEP,                                                      &
     &PDZZ=PDZZF,PZZ=PZZ,                                                                  &
     &PRHODJ=PRHODJ,PRHODREF=PRHODREF,                                                    &
     &PPABSM=PPABSM,PEXNM=PEXNM,                                                          &
     &PSFTH=PSFTH,PSFRV=PSFRV,                                                            &
     &PTHM=PTHM,PRM=PRM,PUM=PUM,PVM=PVM,PTKEM=PTKEM,PSVM=PSVM,                            &
     &PDUDT_MF=PDUDT_MF,PDVDT_MF=PDVDT_MF,                                                &
     &PDTHLDT_MF=PDTHLDT_MF,PDRTDT_MF=PDRTDT_MF,PDSVDT_MF=PDSVDT_MF,                      &
     &PSIGMF=PSIGMF,PRC_MF=PRC_MF,PRI_MF=PRI_MF,PCF_MF=PCF_MF,PFLXZTHVMF=PFLXZTHVMF,      &
     &PFLXZTHMF=ZFLXZTHMF,PFLXZRMF=ZFLXZRMF,PFLXZUMF=ZFLXZUMF,PFLXZVMF=ZFLXZVMF,          &
     &PTHL_UP=PTHL_UP,PRT_UP=PRT_UP,PRV_UP=PRV_UP,PRC_UP=PRC_UP,PRI_UP=PRI_UP,            &
     &PU_UP=PU_UP, PV_UP=PV_UP, PTHV_UP=PTHV_UP, PW_UP=PW_UP,                             &
     &PFRAC_UP=PFRAC_UP,PEMF=PEMF,PDETR=ZDETR,PENTR=ZENTR,                                &
     &KKLCL=IKLCL,KKETL=IKETL,KKCTL=IKCTL,PDX=0.,PDY=0.                                   )
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
