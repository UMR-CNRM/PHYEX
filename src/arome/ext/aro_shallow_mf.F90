!     ######spl
      SUBROUTINE  ARO_SHALLOW_MF(KKL, KLON, KLEV, KFDIA, KRR, KRRL, KRRI,KSV,     &
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
USE YOMMP0, ONLY : MYPROC
USE OML_MOD, ONLY : OML_MY_THREAD
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
CHARACTER(LEN=20) :: CLOFILE
INTEGER           :: IFILE
LOGICAL           :: LFILEEXISTS
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
  IF(.TRUE.) THEN ! .FALSE. to activate the writing, .TRUE. to not activate the writing
    LFILEEXISTS=.FALSE.
  ELSE
    IFILE=0
    LFILEEXISTS=.TRUE.
    DO WHILE(LFILEEXISTS)
      IFILE=IFILE+1
      WRITE(CLOFILE, '(I4.4,"_",I2.2,"_",I8.8,".dat")') MYPROC, OML_MY_THREAD(), IFILE
      INQUIRE(FILE=CLOFILE, EXIST=LFILEEXISTS)
    ENDDO
    IF(IFILE==1 .AND. MYPROC==1) THEN
      PRINT*, "SHALLOWMF CST: KRR, KRRL, KRRI, KSV", KRR, KRRL, KRRI, KSV
      PRINT*, "SHALLOWMF CST: HMF_UPDRAFT, HMF_CLOUD, HFRAC_ICE", HMF_UPDRAFT, HMF_CLOUD, HFRAC_ICE
      PRINT*, "SHALLOWMF CST: OMIXUV, ONOMIXLG, KSV_LGBEG, KSV_LGEND", OMIXUV, ONOMIXLG, KSV_LGBEG, KSV_LGEND
      PRINT*, "SHALLOWMF CST: ZIMPL, PTSTEP", ZIMPL, PTSTEP
    ENDIF
    LFILEEXISTS=IFILE<500
    IF(LFILEEXISTS) THEN
      IFILE=7+OML_MY_THREAD()
      OPEN(IFILE, FILE=CLOFILE, FORM='unformatted') !,access='stream')
      WRITE(IFILE) SIZE(PRHODJ, 1), SIZE(PRHODJ, 2)
      WRITE(IFILE) PDZZF
      !IN
      WRITE(IFILE) PZZ
      WRITE(IFILE) PRHODJ
      WRITE(IFILE) PRHODREF
      WRITE(IFILE) PPABSM
      WRITE(IFILE) PEXNM
      WRITE(IFILE) PSFTH
      WRITE(IFILE) PSFRV
      WRITE(IFILE) PTHM
      WRITE(IFILE) PRM
      WRITE(IFILE) PUM
      WRITE(IFILE) PVM
      WRITE(IFILE) PTKEM
      WRITE(IFILE) PSVM
      !INOUT
      WRITE(IFILE) PTHL_UP
      WRITE(IFILE) PRT_UP
      WRITE(IFILE) PRV_UP
      WRITE(IFILE) PRC_UP
      WRITE(IFILE) PRI_UP
      WRITE(IFILE) PU_UP
      WRITE(IFILE) PV_UP
      WRITE(IFILE) PTHV_UP
      WRITE(IFILE) PW_UP
      WRITE(IFILE) PFRAC_UP
      WRITE(IFILE) PEMF
    ENDIF
  ENDIF
  CALL SHALLOW_MF(YLDIMPHYEX, CST, NEB, PARAM_MFSHALLN, TURBN, CSTURB,                    &
     &KRR=KRR, KRRL=KRRL, KRRI=KRRI, KSV=KSV,                                             &
     &HMF_UPDRAFT=HMF_UPDRAFT, HMF_CLOUD=HMF_CLOUD,HFRAC_ICE=HFRAC_ICE,OMIXUV=OMIXUV,     &
     &ONOMIXLG=ONOMIXLG,KSV_LGBEG=KSV_LGBEG,KSV_LGEND=KSV_LGEND,                          &
     &PIMPL_MF=ZIMPL, PTSTEP=PTSTEP,                                                      &
     &PDZZ=PDZZF,PZZ=PZZ,                                                                 &
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
  IF(LFILEEXISTS) THEN                                                                                                              
    WRITE(IFILE) 
    !OUT
    WRITE(IFILE) PDUDT_MF
    WRITE(IFILE) PDVDT_MF
    WRITE(IFILE) PDTHLDT_MF
    WRITE(IFILE) PDRTDT_MF
    WRITE(IFILE) PDSVDT_MF
    WRITE(IFILE) PSIGMF
    WRITE(IFILE) PRC_MF
    WRITE(IFILE) PRI_MF
    WRITE(IFILE) PCF_MF
    WRITE(IFILE) PFLXZTHVMF
    WRITE(IFILE) ZFLXZTHMF
    WRITE(IFILE) ZFLXZRMF
    WRITE(IFILE) ZFLXZUMF
    WRITE(IFILE) ZFLXZVMF
    !INOUT
    WRITE(IFILE) PTHL_UP
    WRITE(IFILE) PRT_UP
    WRITE(IFILE) PRV_UP
    WRITE(IFILE) PRC_UP
    WRITE(IFILE) PRI_UP
    WRITE(IFILE) PU_UP
    WRITE(IFILE) PV_UP
    WRITE(IFILE) PTHV_UP
    WRITE(IFILE) PW_UP
    WRITE(IFILE) PFRAC_UP
    WRITE(IFILE) PEMF
    !OUT
    WRITE(IFILE) ZDETR
    WRITE(IFILE) ZENTR
    WRITE(IFILE) IKLCL
    WRITE(IFILE) IKETL
    WRITE(IFILE) IKCTL
    CLOSE(IFILE)                                                                                                                    
  ENDIF  
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
