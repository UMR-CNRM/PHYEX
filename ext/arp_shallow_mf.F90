!     ######spl
      SUBROUTINE ARP_SHALLOW_MF(KIDIA,KFDIA,KLON,KTDIA,KLEV,PIMPL,TSPHY,PZZ,PZZF,PR,PCP, &
                               & CMF_UPDRAFT,CMF_CLOUD,LMIXUV, &
                               & PU, PV, PT,PQV,PQL,PQI,PQR,PQS,PTKE,PAPRSF, &
                               & PDELP,PDIFTQ,PDIFTS,PSTRTU,PSTRTV,PSFTH,PSFRV,&
                               & PRODTH_CVPP,PQLI,PNEB,KNLAB,PMF_UP)


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
!!      Documentation ARPEGE
!!
!!    AUTHOR
!!    ------
!!    Y.Bouteloup from aro_shallow_mf
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    11/2010
!!      S. Riette shallow_mf now outputs ice cloud
!!      S. Riette Jan 2012: support for both order of vertical levels
!!      S. Riette April 2022: call abort, waiting for an update from an arpege developper...
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------

USE YOMCST   , ONLY : RG, RATM, RKAPPA, RD, RCPD, RCPV

!USE MODD_PARAMETERS
!
USE MODD_CST, ONLY: CST
USE MODD_NEB, ONLY: NEB
USE MODD_TURB_n, ONLY: TURBN
USE MODD_CTURB, ONLY: CSTURB
USE MODD_PARAM_MFSHALL_n, ONLY: PARAM_MFSHALLN
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODI_SHALLOW_MF
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
USE MODD_CST
USE YOMCT3
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
!

INTEGER,                  INTENT(IN)   :: KIDIA
INTEGER,                  INTENT(IN)   :: KFDIA
INTEGER,                  INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels (bottom of atmosphere in ARP)
INTEGER,                  INTENT(IN)   :: KTDIA    !Top of atmosphere in ARPEGE
REAL,                     INTENT(IN)   :: TSPHY   ! Time step
REAL,                     INTENT(IN)   :: PIMPL

CHARACTER (LEN=4), INTENT(IN)   :: CMF_UPDRAFT  ! Type of Mass Flux Scheme
CHARACTER (LEN=4), INTENT(IN)   :: CMF_CLOUD    ! Type of statistical cloud scheme
LOGICAL,                        INTENT(IN)   :: LMIXUV    ! True if mixing of momentum

!REAL, DIMENSION(KLON,KLEV+2),   INTENT(IN)   :: PZZ     ! Height of layer boundaries 
REAL, DIMENSION(KLON,0:KLEV), INTENT(IN)   :: PZZ     ! Height of layer boundaries 
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PZZF    ! Height of level
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PR   ! Air gaz constant
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PCP  ! Cp
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PU
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PV
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PT
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PQV
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PQL
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PQI
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PQR
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PQS
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PTKE
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PAPRSF
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)   :: PDELP

REAL, DIMENSION(KLON,0:KLEV),   INTENT(INOUT)   :: PDIFTQ
REAL, DIMENSION(KLON,0:KLEV),   INTENT(INOUT)   :: PDIFTS
REAL, DIMENSION(KLON,0:KLEV),   INTENT(INOUT)   :: PSTRTU
REAL, DIMENSION(KLON,0:KLEV),   INTENT(INOUT)   :: PSTRTV

REAL, DIMENSION(KLON,0:KLEV),   INTENT(INOUT)   :: PRODTH_CVPP
REAL, DIMENSION(KLON,KLEV)  ,   INTENT(INOUT)   :: PQLI
REAL, DIMENSION(KLON,KLEV)  ,   INTENT(INOUT)   :: PNEB
INTEGER, DIMENSION(KLON,KLEV),  INTENT(INOUT)   :: KNLAB

REAL, DIMENSION(KLON,0:KLEV),   INTENT(OUT)   :: PMF_UP   

!
! normal surface fluxes of theta and Rv
REAL, DIMENSION(KLON),          INTENT(IN)   ::  PSFTH,PSFRV
!   prognostic variables at t- deltat
!
!CHARACTER (LEN=14),  INTENT(IN)  :: CPNAME
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JRR           ! Loop index for the moist
INTEGER :: IIB           ! Define the physical domain
INTEGER :: IIE           !
INTEGER :: IJB           !
INTEGER :: IJE           !
INTEGER :: IKA
INTEGER :: IKB           !
INTEGER :: IKE           !
INTEGER :: IKU
INTEGER :: IKL           !
INTEGER :: IKR           !
INTEGER :: IKRL           !
INTEGER :: IKRI           !
INTEGER :: JI, JJ, JL, JK, JLON, JLEV !
INTEGER ::II, IUSCM, IKK, ILEV
INTEGER :: ISV_LGBEG, ISV_LGEND, ITCOUNT
INTEGER, DIMENSION(KIDIA:KFDIA)    :: IKLCL,IKETL,IKCTL
REAL,DIMENSION(KIDIA:KFDIA,KLEV+2) :: ZFLXZTHMF,ZFLXZRMF,ZFLXZUMF,ZFLXZVMF
REAL,DIMENSION(KIDIA:KFDIA,KLEV+2) :: ZEMF,ZDETR,ZENTR


REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)   :: ZDZZ,ZZZ,ZTHETA,ZEXNER,ZHRO,ZHRODJ,ZHRODREF
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)   :: ZDUDT_TURB,ZDVDT_TURB,ZDRTDT_TURB,ZDTHLDT_TURB
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2,5) :: ZRM
REAL, DIMENSION(KIDIA:KFDIA)          ::  ZSFTH,ZSFRV

REAL          ::  ZINVG, ZDT, ZEMF_MAX, ZTDCP, ZVMD

REAL, DIMENSION(KIDIA:KFDIA,KLEV+2,1) ::  ZSVM, ZDSVDT_TURB, ZSVDT_MF

CHARACTER (LEN=4) :: HMF_UPDRAFT, HMF_CLOUD


LOGICAL LLOMIXUV, LLONOMIXLG
!
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)   ::  ZDUDT_MF     ! tendency of U   by massflux scheme
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)   ::  ZDVDT_MF     ! tendency of V   by massflux scheme
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)   ::  ZDTHLDT_MF   ! tendency of thl by massflux scheme
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)   ::  ZDRTDT_MF    ! tendency of rt  by massflux scheme

REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZSIGMF,ZRC_MF,ZRI_MF,ZCF_MF ! cloud info for the cloud scheme
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZFLXZTHVMF           ! Thermal production for TKE scheme
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZTHL_UP   ! Thl updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZRT_UP    ! Rt  updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZRV_UP    ! Vapor updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZU_UP     ! U wind updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZV_UP     ! V wind updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZRC_UP    ! cloud content updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZRI_UP    ! ice content   updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZTHV_UP   ! Thv   updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZW_UP     ! vertical speed updraft characteristics
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZFRAC_UP  ! updraft fraction
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZDPSG     ! Delta P / g
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZFQ_MF    ! Flux de qv by massflux scheme
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZFH_MF    ! Flux d'hentalpy by massflux scheme
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZFU_MF
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZFV_MF
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZQDM

REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZAPRSF
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZTKE
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZU
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZV
REAL, DIMENSION(KIDIA:KFDIA,KLEV+2)  ::  ZZZF
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
LOGICAL :: OSTATNW
#include "abor1.intfb.h"

!------------------------------------------------------------------------------

!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------


!  Controle :

!shallow_mf code is now ready to deal with KIDIA/KFDIA
!Array copies can be suppressed (no need to limit the horizontal domain nor to add the two extra levels)
!CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, 0, KFDIA)

!For now, copies are done
CALL FILL_DIMPHYEX(YLDIMPHYEX, KFDIA, 1, KLEV, 1, KFDIA)

CALL ABOR1('ARP_SHALLOW_MF: code must be checked before being activated again')

! Avec inversion des boucles
IKA=1        ! <== Bottom index of array
IKB=2        ! <== Physical bottom
IKE=KLEV+1   ! <== Physical top
IKU=KLEV+2   ! <== Top index of array
IKL=1        ! <== Loop increment to go from top to bottom

IIB=KIDIA
IIE=KFDIA
IKR=5         ! <== Number of water species
IKRL=2
IKRI=2
ZINVG = 1./RG

!------------------------------------------------------------------------------

!*       2.   INITIALISATION

!             ---------------

! tableau a recalculer a chaque pas de temps
! attention, ZDZZ est l'altitude entre deux niveaux (et pas l'�paisseur de la couche)

! Inversion des niveaux

DO JK = IKB , IKE, IKL
   DO JL = IIB,IIE
      IKK = IKE + 1 - JK
      ZAPRSF(JL,JK) = PAPRSF(JL,IKK)
      ZTKE  (JL,JK) = PTKE  (JL,IKK)
      ZU    (JL,JK) = PU    (JL,IKK)
      ZV    (JL,JK) = PV    (JL,IKK)
   ENDDO
ENDDO


! AROME type initialisation
  !initialisation de ZZZ
DO JK = IKB , IKE+1
  DO JL = IIB,IIE
    IKK = IKE + 1 - JK
    ZZZ(JL,JK) = PZZ(JL,IKK)*ZINVG
  ENDDO
ENDDO


DO JL = IIB,IIE
  ZZZ(JL,1) = 2*ZZZ(JL,2)-ZZZ(JL,3)
ENDDO
!initialisation de ZZZF
DO JK = IKB , IKE
  DO JL = IIB,IIE
    IKK = IKE + 1 - JK
    ZZZF(JL,JK) = PZZF(JL,IKK)*ZINVG
  ENDDO
ENDDO
DO JL = IIB,IIE
  ZZZF(JL,1)=1.5*ZZZ(JL,2)-0.5*ZZZ(JL,3)
  ZZZF(JL,IKE+1)=ZZZF(JL,IKE)+ (ZZZ(JL,IKE+1)-ZZZ(JL,IKE))
  ZDZZ(JL,1)=-999.
ENDDO

DO JK = IKB , IKE+1
   DO JL = IIB,IIE
      ZDZZ(JL,JK)=ZZZF(JL,JK)-ZZZF(JL,JK-1)
   ENDDO
ENDDO

! Inversion des niveaux
DO JL = IIB,IIE
   DO JK = IKB , IKE, IKL
      IKK = IKE + 1 - JK
      ZEXNER(JL,JK)=(ZAPRSF(JL,JK)/RATM)**RKAPPA
      ZTHETA(JL,JK)=PT(JL,IKK)/ZEXNER(JL,JK)   
      ZHRO(JL,JK)=ZAPRSF(JL,JK)/(PT(JL,IKK)*PR(JL,IKK))
      ZQDM(JL,JK)=(1.-PQV(JL,IKK)-PQL(JL,IKK)-PQI(JL,IKK)-PQR(JL,IKK)-PQS(JL,IKK))

      ZHRODREF(JL,JK)=ZHRO(JL,JK)*ZQDM(JL,JK)
      ZHRODJ(JL,JK)=PDELP(JL,IKK)*ZINVG
      ZRM(JL,JK,1)=PQV(JL,IKK)/ZQDM(JL,JK)
      ZRM(JL,JK,2)=PQL(JL,IKK)/ZQDM(JL,JK)
      ZRM(JL,JK,3)=PQR(JL,IKK)/ZQDM(JL,JK)
      ZRM(JL,JK,4)=PQI(JL,IKK)/ZQDM(JL,JK)
      ZRM(JL,JK,5)=PQS(JL,IKK)/ZQDM(JL,JK)

      ZDPSG(JL,JK) = MAX(1.E-15,PDELP(JL,IKK)/RG)

!  Copy KLEV array into KLEV+2 array, or how to spend cpu time unnecessarily but that's the arome physics rule !
   ENDDO
ENDDO

ZDUDT_TURB(KIDIA:KFDIA,:)   = 0.
ZDVDT_TURB(KIDIA:KFDIA,:)   = 0.
ZDRTDT_TURB(KIDIA:KFDIA,:)  = 0.
ZDTHLDT_TURB(KIDIA:KFDIA,:) = 0.
ZDSVDT_TURB(KIDIA:KFDIA,:,:)= 0.

ZSVM(:,:,:)=0.

 ZAPRSF  (KIDIA:KFDIA,IKB-IKL)     = ZAPRSF  (KIDIA:KFDIA,IKB)
 ZAPRSF  (KIDIA:KFDIA,IKE+IKL)     = ZAPRSF  (KIDIA:KFDIA,IKE)
 ZTKE    (KIDIA:KFDIA,IKB-IKL)     = ZTKE    (KIDIA:KFDIA,IKB)
 ZTKE    (KIDIA:KFDIA,IKE+IKL)     = ZTKE    (KIDIA:KFDIA,IKE)
 ZU      (KIDIA:KFDIA,IKB-IKL)     = ZU      (KIDIA:KFDIA,IKB)
 ZU      (KIDIA:KFDIA,IKE+IKL)     = ZU      (KIDIA:KFDIA,IKE)
 ZV      (KIDIA:KFDIA,IKB-IKL)     = ZV      (KIDIA:KFDIA,IKB)
 ZV      (KIDIA:KFDIA,IKE+IKL)     = ZV      (KIDIA:KFDIA,IKE)
 ZEXNER  (KIDIA:KFDIA,IKB-IKL)     = ZEXNER  (KIDIA:KFDIA,IKB)
 ZEXNER  (KIDIA:KFDIA,IKE+IKL)     = ZEXNER  (KIDIA:KFDIA,IKE)
 ZTHETA  (KIDIA:KFDIA,IKB-IKL)     = ZTHETA  (KIDIA:KFDIA,IKB)
 ZTHETA  (KIDIA:KFDIA,IKE+IKL)     = ZTHETA  (KIDIA:KFDIA,IKE)
 ZHRO    (KIDIA:KFDIA,IKB-IKL)     = ZHRO    (KIDIA:KFDIA,IKB)
 ZHRO    (KIDIA:KFDIA,IKE+IKL)     = ZHRO    (KIDIA:KFDIA,IKE)
 ZQDM    (KIDIA:KFDIA,IKB-IKL)     = ZQDM    (KIDIA:KFDIA,IKB)
 ZQDM    (KIDIA:KFDIA,IKE+IKL)     = ZQDM    (KIDIA:KFDIA,IKE)
 ZHRODREF(KIDIA:KFDIA,IKB-IKL)     = ZHRODREF(KIDIA:KFDIA,IKB)
 ZHRODREF(KIDIA:KFDIA,IKE+IKL)     = ZHRODREF(KIDIA:KFDIA,IKE)
 ZHRODJ  (KIDIA:KFDIA,IKB-IKL)     = ZHRODJ  (KIDIA:KFDIA,IKB)
 ZHRODJ  (KIDIA:KFDIA,IKE+IKL)     = ZHRODJ  (KIDIA:KFDIA,IKE)
 ZRM     (KIDIA:KFDIA,IKB-IKL,1)   = ZRM     (KIDIA:KFDIA,IKB,1)
 ZRM     (KIDIA:KFDIA,IKE+IKL,1)   = ZRM     (KIDIA:KFDIA,IKE,1)
 ZRM     (KIDIA:KFDIA,IKB-IKL,2)   = ZRM     (KIDIA:KFDIA,IKB,2)
 ZRM     (KIDIA:KFDIA,IKE+IKL,2)   = ZRM     (KIDIA:KFDIA,IKE,2)
 ZRM     (KIDIA:KFDIA,IKB-IKL,3)   = ZRM     (KIDIA:KFDIA,IKB,3)
 ZRM     (KIDIA:KFDIA,IKE+IKL,3)   = ZRM     (KIDIA:KFDIA,IKE,3)
 ZRM     (KIDIA:KFDIA,IKB-IKL,4)   = ZRM     (KIDIA:KFDIA,IKB,4)
 ZRM     (KIDIA:KFDIA,IKE+IKL,4)   = ZRM     (KIDIA:KFDIA,IKE,4)
 ZRM     (KIDIA:KFDIA,IKB-IKL,5)   = ZRM     (KIDIA:KFDIA,IKB,5)
 ZRM     (KIDIA:KFDIA,IKE+IKL,5)   = ZRM     (KIDIA:KFDIA,IKE,5)
 ZSVM    (KIDIA:KFDIA,IKB-IKL,:)   = ZSVM    (KIDIA:KFDIA,IKB,:)
 ZSVM    (KIDIA:KFDIA,IKE+IKL,:)   = ZSVM    (KIDIA:KFDIA,IKE,:)
 ZDUDT_TURB(KIDIA:KFDIA,IKB-IKL)   = ZDUDT_TURB(KIDIA:KFDIA,IKB)   
 ZDUDT_TURB(KIDIA:KFDIA,IKE+IKL)   = ZDUDT_TURB(KIDIA:KFDIA,IKE)  
 ZDVDT_TURB(KIDIA:KFDIA,IKB-IKL)   = ZDVDT_TURB(KIDIA:KFDIA,IKB)    
 ZDVDT_TURB(KIDIA:KFDIA,IKE+IKL)   = ZDVDT_TURB(KIDIA:KFDIA,IKE)  
 ZDTHLDT_TURB(KIDIA:KFDIA,IKB-IKL) = ZDTHLDT_TURB(KIDIA:KFDIA,IKB)  
 ZDTHLDT_TURB(KIDIA:KFDIA,IKE+IKL) = ZDTHLDT_TURB(KIDIA:KFDIA,IKE)  
 ZDRTDT_TURB(KIDIA:KFDIA,IKB-IKL)  = ZDRTDT_TURB(KIDIA:KFDIA,IKB)  
 ZDRTDT_TURB(KIDIA:KFDIA,IKE+IKL)  = ZDRTDT_TURB(KIDIA:KFDIA,IKE)  
 ZDSVDT_TURB(KIDIA:KFDIA,IKB-IKL,:)= ZDSVDT_TURB(KIDIA:KFDIA,IKB,:)
 ZDSVDT_TURB(KIDIA:KFDIA,IKE+IKL,:)= ZDSVDT_TURB(KIDIA:KFDIA,IKE,:)

DO JL = IIB,IIE
  ZSFTH(JL) = -PSFTH(JL)/ZHRO(JL,IKB)/RCPD
  ZSFRV(JL) = -PSFRV(JL)/ZHRO(JL,IKB)
ENDDO

LLOMIXUV    = .TRUE.
HMF_UPDRAFT =  CMF_UPDRAFT
HMF_CLOUD   =  CMF_CLOUD
LLOMIXUV    =  LMIXUV
LLONOMIXLG  = .FALSE.
ISV_LGBEG   =  0
ISV_LGEND   =  0
ITCOUNT     =  1
ZDT         =  TSPHY


!  Mise � 0 des tendances 

ZDUDT_MF(:,:)   = 0.
ZDVDT_MF(:,:)   = 0.
ZDTHLDT_MF(:,:) = 0.
ZDRTDT_MF(:,:)  = 0.

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
OSTATNW = .FALSE.
  CALL SHALLOW_MF(YLDIMPHYEX, CST, NEB, PARAM_MFSHALLN, TURBN, CSTURB,   &
       KRR=IKR,KRRL=IKRL,KRRI=IKRI, KSV=1,                             &
       HMF_UPDRAFT=HMF_UPDRAFT, HMF_CLOUD=HMF_CLOUD,HFRAC_ICE='N',OMIXUV=LLOMIXUV,     &
       OSTATNW=OSTATNW,                                                  &
       ONOMIXLG=LLONOMIXLG,KSV_LGBEG=ISV_LGBEG,KSV_LGEND=ISV_LGEND,      &
      PIMPL_MF=PIMPL, PTSTEP=ZDT,                                        &
      PDZZ=ZDZZ,PZZ=ZZZ,                                                 &
      PRHODJ=ZHRODJ,PRHODREF=ZHRODREF,                                   &
      PPABSM=ZAPRSF,PEXNM=ZEXNER,                                        &
      PSFTH=ZSFTH,PSFRV=ZSFRV,                                           &
      PTHM=ZTHETA,PRM=ZRM,PUM=ZU,PVM=ZV,PTKEM=ZTKE,PSVM=ZSVM,            &
!  Output
      PDUDT_MF=ZDUDT_MF,PDVDT_MF=ZDVDT_MF,                                      &
      PDTHLDT_MF=ZDTHLDT_MF,PDRTDT_MF=ZDRTDT_MF,PDSVDT_MF=ZSVDT_MF,             &
      PSIGMF=ZSIGMF,PRC_MF=ZRC_MF,PRI_MF=ZRI_MF,PCF_MF=ZCF_MF,PFLXZTHVMF=ZFLXZTHVMF,          &
      PFLXZTHMF=ZFLXZTHMF,PFLXZRMF=ZFLXZRMF,PFLXZUMF=ZFLXZUMF,PFLXZVMF=ZFLXZVMF,&
      PTHL_UP=ZTHL_UP,PRT_UP=ZRT_UP,PRV_UP=ZRV_UP,PRC_UP=ZRC_UP,PRI_UP=ZRI_UP,  &
      PU_UP=ZU_UP, PV_UP=ZV_UP, PTHV_UP=ZTHV_UP, PW_UP=ZW_UP,                   &
      PFRAC_UP=ZFRAC_UP,PEMF=ZEMF,PDETR=ZDETR,PENTR=ZENTR,                      &
      KKLCL=IKLCL,KKETL=IKETL,KKCTL=IKCTL,                                      &
!
      PDX=0., PDY=0.)


!  Conversion des tendances de theta en tendance de cpT
!  et conversion en qi en multipliant par qd
! Puis calcul des flux

ZFQ_MF(:,:) = 0.
ZFH_MF(:,:) = 0.
ZFU_MF(:,:) = 0.
ZFV_MF(:,:) = 0.


ZVMD=RCPV-RCPD

DO JL = IIB,IIE
   DO JK = IKE , IKB, -IKL   ! Loop from top to bottom

      IKK = IKE + 1 - JK

      ZFQ_MF(JL,JK) = ZFQ_MF(JL,JK+IKL) - ZDPSG(JL,JK)*ZDRTDT_MF(JL,JK)*ZQDM(JL,JK)
      ZTDCP=ZVMD*ZDRTDT_MF(JL,JK)



      ZFH_MF(JL,JK) = ZFH_MF(JL,JK+IKL) - ZDPSG(JL,JK) &
      & * (ZDTHLDT_MF(JL,JK)*ZEXNER(JL,JK)*(PCP(JL,IKK)+TSPHY*ZTDCP)+PT(JL,IKK)*ZTDCP)

!      ZFH_MF(JL,JK) = ZFH_MF(JL,JK+IKL) - ZDPSG(JL,JK) &
!      & * (ZDTHLDT_MF(JL,JK)*ZEXNER(JL,JK)*PCP(JL,IKK)+PT(JL,IKK)*ZTDCP)

      ZFU_MF(JL,JK) = ZFU_MF(JL,JK+IKL) - ZDPSG(JL,JK)*ZDUDT_MF(JL,JK)
      ZFV_MF(JL,JK) = ZFV_MF(JL,JK+IKL) - ZDPSG(JL,JK)*ZDVDT_MF(JL,JK)
   ENDDO
ENDDO

ZRC_UP(:,:) = ZRC_UP(:,:)*ZFRAC_UP(:,:)
ZRI_UP(:,:) = ZRI_UP(:,:)*ZFRAC_UP(:,:)
PRODTH_CVPP(:,:) = 0.


! stockage dans les flux turbulents (Inversion des niveaux !)

DO JL = IIB,IIE
   DO JK = IKE , IKB, -IKL   ! Loop from top to bottom
      IKK = IKE + 1 - JK
      PDIFTQ(JL,IKK) = PDIFTQ(JL,IKK) + ZFQ_MF(JL,JK)
      PDIFTS(JL,IKK) = PDIFTS(JL,IKK) + ZFH_MF(JL,JK)
      PSTRTU(JL,IKK) = PSTRTU(JL,IKK) + ZFU_MF(JL,JK)
      PSTRTV(JL,IKK) = PSTRTV(JL,IKK) + ZFV_MF(JL,JK)
      PRODTH_CVPP(JL,IKK) = RG/ZTHETA(JL,JK)*ZFLXZTHVMF(JL,JK)

! Shallow cloud information     
      PQLI  (JL,IKK) =  (ZRC_MF(JL,JK)+ZRI_MF(JL,JK))/(1.+ZRT_UP(JL,JK)) ! with HFRAC_ICE='N', ZRI_MF=0
      KNLAB (JL,IKK) =  INT(MAX(0.,SIGN(1.,PQLI(JL,IKK)-1.E-8)))
      PNEB  (JL,IKK) =  ZCF_MF(JL,JK) 
      PMF_UP(JL,IKK) = -ZEMF(JL,JK)/ZHRODREF(JL,JK) ! <== On ne sort pas le flux de masse mais
!                                                    ! ce dont on aura besoin dans ACDIFV1 !!!!!!!
   ENDDO
ENDDO

END SUBROUTINE ARP_SHALLOW_MF
