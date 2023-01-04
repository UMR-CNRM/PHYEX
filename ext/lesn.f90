!MNH_LIC Copyright 2000-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #################
      SUBROUTINE  LES_n
!     #################
!
!
!!****  *LES_n* computes the current time-step LES diagnostics for model _n
!!
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original  07/02/00
!!                01/02/01 (D. Gazen) add module MODD_NSV for NSV variable
!!                06/11/02 (V. Masson) add LES budgets and use of anomalies
!!                         in LES quantities computations
!!                01/04/03 (V. Masson and F. Couvreux) bug in BL height loop
!!                   10/07 (J.Pergaud) Add mass flux diagnostics
!!                06/08    (O.Thouron) Add radiative diagnostics
!!                12/10    (R.Honnert) Add EDKF mass flux in BL height
!!                10/09    (P. Aumond) Add possibility of user maskS
!!                10/14    (C.Lac) Correction on user masks
!!                10/16    (C.Lac) Add ground droplet deposition amount
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                     02/2019 (C. Lac) Add rain fraction as a LES diagnostic
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_CTURB, ONLY : XFTOP_O_FSURF
!
USE MODD_LES
USE MODD_LES_BUDGET
USE MODD_CONF
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_LES_n
USE MODD_RADIATIONS_n
USE MODD_GRID_n
USE MODD_REF_n
USE MODD_FIELD_n
USE MODD_CONF_n
USE MODD_PARAM_n
USE MODD_TURB_n
USE MODD_METRICS_n
USE MODD_LUNIT_n, ONLY: TLUOUT
USE MODD_PARAM_n, ONLY: CCLOUD
USE MODD_PRECIP_n, ONLY: XINPRR,XACPRR,XINPRR3D,XEVAP3D,XINPRC,XINDEP
USE MODD_NSV, ONLY : NSV, NSV_CS
USE MODD_PARAM_ICE, ONLY: LDEPOSC,LSEDIC
USE MODD_PARAM_C2R2, ONLY: LDEPOC,LSEDC
USE MODD_PARAM_LIMA, ONLY : MSEDC=>LSEDC
!
USE MODI_SHUMAN
USE MODI_GRADIENT_M
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_W
USE MODI_LES_VER_INT
USE MODI_SPEC_VER_INT
USE MODI_LES_MEAN_ll
USE MODI_THL_RT_FROM_TH_R
USE MODI_LES_RES_TR
USE MODI_BUDGET_FLAGS
USE MODI_LES_BUDGET_TEND_n
USE MODE_BL_DEPTH_DIAG
!
USE MODE_ll
USE MODE_MODELN_HANDLER
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZEXN      ! Exner function
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHL      ! liquid potential temperature
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHV      ! virtual potential temperature
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRHO      ! air density
!
REAL, DIMENSION(:,:,:),    ALLOCATABLE :: CHAMPXY1    !tableau intermediaire
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTEMP     ! Temperature
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZEW
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZINDCLD   !indice cloud si rc>0
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZINDCLD2  !indice cloud rc>1E-5
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZCLDFR_LES! CLDFR    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZICEFR_LES! ICEFR    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRAINFR_LES! RAINFR   on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZMASSF    ! massflux=rho*w
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZREHU     ! relative humidity


REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZZ_LES    ! alt.  on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE ::ZZZ_LES
REAL, DIMENSION(:,:,:),   ALLOCATABLE ::ZINPRR3D_LES   ! precipitation flux 3D
REAL, DIMENSION(:,:,:),   ALLOCATABLE ::ZEVAP3D_LES !evaporation 3D
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZP_LES    ! pres. on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDP_LES   ! dynamical production TKE   
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTP_LES   ! thermal production TKE    
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTR_LES   ! transport production TKE    
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDISS_LES ! dissipation TKE    
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLM_LES    ! mixing length

REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDPDZ_LES   ! dp/dz on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDTHLDZ_LES ! dThl/dz on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDTHDZ_LES ! dTh/dz on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDRTDZ_LES  ! dRt/dz  on LES vertical grid
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZDSvDZ_LES  ! dSv/dz  on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDUDZ_LES ! du/dz on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDVDZ_LES ! dv/dz on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDWDZ_LES ! dw/dz on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZEXN_LES  ! Exner on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRHO_LES  ! rho   on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZU_LES    ! U     on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZV_LES    ! V     on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZW_LES    ! W     on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZMF_LES   ! mass flux on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTH_LES   ! Theta on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHV_LES  ! thv   on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHL_LES  ! thl   on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTKE_LES  ! tke   on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZKE_LES   ! ke   on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRV_LES   ! Rv    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZREHU_LES ! Rehu  on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRC_LES   ! Rc    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRR_LES   ! Rr    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRI_LES   ! Ri    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRS_LES   ! Rs    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRG_LES   ! Rg    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRH_LES   ! Rh    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRT_LES   ! Rt    on LES vertical grid
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSV_LES   ! Sv    on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTH_ANOM  ! Theta anomaly on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHV_ANOM ! thv   anomaly on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRV_ANOM  ! Rv    anomaly on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRC_ANOM  ! Rc    anomaly on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRI_ANOM  ! Ri    anomaly on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRR_ANOM  ! Rr    anomaly on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZP_ANOM   ! p     anomaly on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRHO_ANOM ! rho   anomaly on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDPDZ_ANOM! dp/dz anomaly on LES vertical grid
REAL, DIMENSION(:),       ALLOCATABLE :: ZMEAN_DPDZ! dp/dz mean    on LES vertical grid
REAL, DIMENSION(:),       ALLOCATABLE :: ZLES_MEAN_DRtDZ! drt/dz mean    on LES vertical grid
REAL, DIMENSION(:),       ALLOCATABLE :: ZLES_MEAN_DTHDZ! dth/dz mean    on LES vertical grid
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZLES_MEAN_DSVDZ! drt/dz mean    on LES vertical grid
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZLWP_LES, ZRWP_LES, ZTKET_LES
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZIWP_LES, ZSWP_LES, ZGWP_LES, ZHWP_LES
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZINDCLD2D  !
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZINDCLD2D2 !
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZLWP_ANOM ! lwp anomaly
REAL, DIMENSION(:,:),     ALLOCATABLE :: ZMAXWRR2D ! maxwrr2D
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZU_SPEC   ! U     on SPEC vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZV_SPEC   ! V     on SPEC vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZW_SPEC   ! W     on SPEC vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTH_SPEC  ! Theta on SPEC vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHL_SPEC ! thl   on SPEC vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRV_SPEC  ! Rv    on SPEC vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRC_SPEC  ! Rc    on SPEC vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRI_SPEC  ! Ri    on SPEC vertical grid
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSV_SPEC  ! Sv    on SPEC vertical grid
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE  :: ZRT      ! rv+rc+rr+ri+rs+rg+rh
REAL, DIMENSION(:),       ALLOCATABLE  :: ZWORK1D,ZWORK1DT
REAL, DIMENSION(:,:),     ALLOCATABLE  :: ZWORK2D
REAL :: ZINPRRm,ZCOUNT
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRADEFF_LES   ! Re on LES vertical grid
!!fl sw, lw, dthrad on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZSWU_LES   ! SWU on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZSWD_LES   ! SWD on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLWU_LES   ! LWU on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZLWD_LES   ! LWD on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDTHRADSW_LES   ! DTHRADSW on LES vertical grid
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZDTHRADLW_LES   ! DTHRADLW on LES vertical grid
!
REAL, DIMENSION(:),   ALLOCATABLE :: ZWORK   ! 
!
INTEGER :: IRR      ! moist variables counter
INTEGER :: JSV      ! scalar variables counter
INTEGER :: IIU, IJU ! array sizes
INTEGER :: IKE,IKB
INTEGER :: JI, JJ, JK   ! loop counters
INTEGER :: IIU_ll, IJU_ll    ! total domain I size (fin)
INTEGER :: IIA_ll, IJA_ll    ! total domain I size (debut)
INTEGER :: IINFO_ll      ! return code of parallel routine
INTEGER :: IIMAX_ll, IJMAX_ll  ! total physical domain I size
INTEGER :: JLOOP
!
INTEGER :: IMASK    ! mask counter
INTEGER :: IMASKUSER! mask user number 
!
INTEGER :: IRESP, ILUOUT
INTEGER :: IMI      ! Current model index
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
!-------------------------------------------------------------------------------
!
IMI = GET_CURRENT_MODEL_INDEX()
!
IF (.NOT. LLES_CALL) RETURN
!
CALL GET_GLOBALDIMS_ll(IIMAX_ll,IJMAX_ll)
IIU_ll = IIMAX_ll+JPHEXT
IJU_ll = IJMAX_ll+JPHEXT
IIA_ll=JPHEXT+1
IJA_ll=JPHEXT+1
IKE=SIZE(XVT,3)-JPVEXT
IKB=1+JPVEXT
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL FILL_DIMPHYEX(YLDIMPHYEX, SIZE(XTHT,1), SIZE(XTHT,2), SIZE(XTHT,3),.TRUE.)
!
ILUOUT = TLUOUT%NLU
!
!-------------------------------------------------------------------------------
!
!* interpolation coefficients for Z type grid
!
IF (CSPECTRA_LEVEL_TYPE=='Z') THEN
  IF (ASSOCIATED(XCOEFLIN_CURRENT_SPEC)) CALL LES_DEALLOCATE('XCOEFLIN_CURRENT_SPEC')
  IF (ASSOCIATED(NKLIN_CURRENT_SPEC   )) CALL LES_DEALLOCATE('NKLIN_CURRENT_SPEC')
  !
  CALL LES_ALLOCATE('XCOEFLIN_CURRENT_SPEC',(/IIU,IJU,NSPECTRA_K/))
  CALL LES_ALLOCATE('NKLIN_CURRENT_SPEC',(/IIU,IJU,NSPECTRA_K/))
  !
  XCOEFLIN_CURRENT_SPEC(:,:,:) = XCOEFLIN_SPEC(:,:,:)
  NKLIN_CURRENT_SPEC   (:,:,:) = NKLIN_SPEC   (:,:,:)
END IF
!
!-------------------------------------------------------------------------------
!
!*      1.   Allocations
!            -----------
!
ALLOCATE(ZP_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZDP_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZTP_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZTR_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZDISS_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZLM_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZDTHLDZ_LES(IIU,IJU,NLES_K))
ALLOCATE(ZDTHDZ_LES(IIU,IJU,NLES_K))
ALLOCATE(ZDRTDZ_LES(IIU,IJU,NLES_K))
ALLOCATE(ZDUDZ_LES(IIU,IJU,NLES_K))
ALLOCATE(ZDVDZ_LES(IIU,IJU,NLES_K))
ALLOCATE(ZDWDZ_LES(IIU,IJU,NLES_K))
ALLOCATE(ZDSVDZ_LES(IIU,IJU,NLES_K,NSV))

ALLOCATE(ZDPDZ_LES(IIU,IJU,NLES_K))
ALLOCATE(ZEXN_LES (IIU,IJU,NLES_K))
ALLOCATE(ZRHO_LES (IIU,IJU,NLES_K))
ALLOCATE(ZU_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZV_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZW_LES   (IIU,IJU,NLES_K))
ALLOCATE(ZMF_LES  (IIU,IJU,NLES_K))
ALLOCATE(ZTH_LES  (IIU,IJU,NLES_K))
IF (CRAD /= 'NONE') THEN
  ALLOCATE(ZRADEFF_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZSWU_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZSWD_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZLWU_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZLWD_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZDTHRADSW_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZDTHRADLW_LES  (IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZRADEFF_LES  (0,0,0))
  ALLOCATE(ZSWU_LES  (0,0,0))
  ALLOCATE(ZSWD_LES  (0,0,0))
  ALLOCATE(ZLWU_LES  (0,0,0))
  ALLOCATE(ZLWD_LES  (0,0,0))
  ALLOCATE(ZDTHRADSW_LES  (0,0,0))
  ALLOCATE(ZDTHRADLW_LES  (0,0,0))
END IF
IF (LUSERV) THEN
  ALLOCATE(ZTHV_LES (IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZTHV_LES (0,0,0))
END IF
ALLOCATE(ZTHL_LES (IIU,IJU,NLES_K))
ALLOCATE(ZTKE_LES (IIU,IJU,NLES_K))
ALLOCATE(ZKE_LES(IIU,IJU,NLES_K))
ALLOCATE(ZTKET_LES(IIU,IJU))
ALLOCATE(ZWORK1D (NLES_K))
ALLOCATE(ZWORK1DT (NLES_K))
ALLOCATE(ZZZ_LES(IIU,IJU,NLES_K))
IF (LUSERV) THEN
  ALLOCATE(ZRV_LES    (IIU,IJU,NLES_K))
  ALLOCATE(ZRT_LES    (IIU,IJU,NLES_K))
  ALLOCATE(ZREHU_LES  (IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZRV_LES    (0,0,0))
  ALLOCATE(ZRT_LES    (0,0,0))
  ALLOCATE(ZREHU_LES  (0,0,0))
END IF
IF (LUSERC) THEN
  ALLOCATE(ZRC_LES    (IIU,IJU,NLES_K))
  ALLOCATE(ZLWP_LES(IIU,IJU))
  ALLOCATE(ZINDCLD2D(IIU,IJU))
  ALLOCATE(ZINDCLD2D2(IIU,IJU))
  ALLOCATE(ZCLDFR_LES(IIU,IJU,NLES_K))
  ALLOCATE(ZWORK2D(IIU,IJU))
  ALLOCATE(ZLWP_ANOM(IIU,IJU))
ELSE
  ALLOCATE(ZRC_LES    (0,0,0))
  ALLOCATE(ZLWP_LES(0,0))
  ALLOCATE(ZINDCLD2D(0,0))
  ALLOCATE(ZINDCLD2D2(0,0))
  ALLOCATE(ZCLDFR_LES(0,0,0))
  ALLOCATE(ZWORK2D(0,0))
  ALLOCATE(ZLWP_ANOM(0,0))
END IF
IF (LUSERR) THEN
  ALLOCATE(ZRR_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZMAXWRR2D(IIU,IJU))
  ALLOCATE(ZRWP_LES(IIU,IJU))
  ALLOCATE(ZINPRR3D_LES (IIU,IJU,NLES_K))
  ALLOCATE(ZEVAP3D_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZRAINFR_LES(IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZRR_LES  (0,0,0))
  ALLOCATE(ZMAXWRR2D(0,0))
  ALLOCATE(ZRWP_LES(0,0))
  ALLOCATE(ZINPRR3D_LES(0,0,0))
  ALLOCATE(ZEVAP3D_LES(0,0,0))
  ALLOCATE(ZRAINFR_LES(0,0,0))
END IF
IF (LUSERI) THEN
  ALLOCATE(ZRI_LES    (IIU,IJU,NLES_K))
  ALLOCATE(ZIWP_LES(IIU,IJU))
  ALLOCATE(ZICEFR_LES(IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZRI_LES    (0,0,0))
  ALLOCATE(ZIWP_LES(0,0))
  ALLOCATE(ZICEFR_LES(0,0,0))
END IF
IF (LUSERS) THEN
  ALLOCATE(ZRS_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZSWP_LES(IIU,IJU))
ELSE
  ALLOCATE(ZRS_LES  (0,0,0))
  ALLOCATE(ZSWP_LES(0,0))
END IF
IF (LUSERG) THEN
  ALLOCATE(ZRG_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZGWP_LES(IIU,IJU))
ELSE
  ALLOCATE(ZRG_LES  (0,0,0))
  ALLOCATE(ZGWP_LES(0,0))
END IF
IF (LUSERH) THEN
  ALLOCATE(ZRH_LES  (IIU,IJU,NLES_K))
  ALLOCATE(ZHWP_LES(IIU,IJU))
ELSE
  ALLOCATE(ZRH_LES  (0,0,0))
  ALLOCATE(ZHWP_LES(0,0))
END IF
IF (NSV>0) THEN
  ALLOCATE(ZSV_LES  (IIU,IJU,NLES_K,NSV))
ELSE
  ALLOCATE(ZSV_LES  (0,0,0,0))
END IF
!
ALLOCATE(ZP_ANOM   (IIU,IJU,NLES_K))
ALLOCATE(ZRHO_ANOM (IIU,IJU,NLES_K))
ALLOCATE(ZTH_ANOM  (IIU,IJU,NLES_K))
ALLOCATE(ZDPDZ_ANOM(IIU,IJU,NLES_K))
IF (LUSERV) THEN
  ALLOCATE(ZTHV_ANOM(IIU,IJU,NLES_K))
  ALLOCATE(ZRV_ANOM (IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZTHV_ANOM(0,0,0))
  ALLOCATE(ZRV_ANOM (0,0,0))
END IF
IF (LUSERC) THEN
  ALLOCATE(ZRC_ANOM (IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZRC_ANOM (0,0,0))
END IF
IF (LUSERI) THEN
  ALLOCATE(ZRI_ANOM (IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZRI_ANOM (0,0,0))
END IF
IF (LUSERR) THEN
  ALLOCATE(ZRR_ANOM (IIU,IJU,NLES_K))
ELSE
  ALLOCATE(ZRR_ANOM (0,0,0))
END IF
ALLOCATE(ZMEAN_DPDZ(NLES_K))
ALLOCATE(ZLES_MEAN_DTHDZ(NLES_K))
!
!
ALLOCATE(ZU_SPEC  (NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K))
ALLOCATE(ZV_SPEC  (NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K))
ALLOCATE(ZW_SPEC  (NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K))
ALLOCATE(ZTH_SPEC (NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K))
IF (LUSERC) THEN
  ALLOCATE(ZTHL_SPEC(NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K))
ELSE
  ALLOCATE(ZTHL_SPEC(0,0,0))
END IF
IF (LUSERV) THEN
  ALLOCATE(ZRV_SPEC  (NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K))
ELSE
  ALLOCATE(ZRV_SPEC  (0,0,0))
END IF
IF (LUSERC) THEN
  ALLOCATE(ZRC_SPEC  (NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K))
ELSE
  ALLOCATE(ZRC_SPEC  (0,0,0))
END IF
IF (LUSERI) THEN
  ALLOCATE(ZRI_SPEC  (NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K))
ELSE
  ALLOCATE(ZRI_SPEC  (0,0,0))
END IF
IF (NSV>0) THEN
  ALLOCATE(ZSV_SPEC  (NSPECTRA_NI,NSPECTRA_NJ,NSPECTRA_K,NSV))
ELSE
  ALLOCATE(ZSV_SPEC  (0,0,0,0))
END IF
!
!
ALLOCATE(ZEXN (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(ZRHO (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(ZRT  (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(ZTHV (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(ZTHL (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(ZEW  (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(ZMASSF (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(ZTEMP (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(ZREHU (IIU,IJU,SIZE(XTHT,3)))
ALLOCATE(CHAMPXY1 (IIU,IJU,1))
!
!-------------------------------------------------------------------------------
!
!*      1.2  preliminary calculations
!            ------------------------
!
ZEXN(:,:,:) = (XPABST/XP00)**(XRD/XCPD)
!
!
!* computation of relative humidity
ZTEMP=XTHT*ZEXN
ZEW=EXP (XALPW -XBETAW/ZTEMP-XGAMW*ALOG(ZTEMP))
IF (LUSERV) THEN
  ZREHU(:,:,:)=100.*XRT(:,:,:,1)*XPABST(:,:,:)/((XRD/XRV+XRT(:,:,:,1))*ZEW(:,:,:))
ELSE
  ZREHU(:,:,:)=0.
END IF
!
CALL THL_RT_FROM_TH_R(LUSERV, LUSERC, LUSERR,             &
                      LUSERI, LUSERS, LUSERG, LUSERH,     &
                      XCURRENT_L_O_EXN_CP,                &
                      XTHT, XRT,                          &
                      ZTHL, ZRT                           )
!
!* computation of density and virtual potential temperature
!
ZTHV=XTHT
IF (LUSERV) ZTHV=ZTHV*(1.+XRV/XRD*XRT(:,:,:,1))/(1.+ZRT(:,:,:))
!
IF (CEQNSYS=='DUR') THEN
  ZRHO=XPABST/(XRD*ZTHV*ZEXN)
ELSE
  ZRHO=XRHODREF*( 1. + (XCPD-XRD)/XRD*(ZEXN/XEXNREF - 1.) - (ZTHV/XTHVREF - 1.) )
END IF
!
! computation of mass flux
ZMASSF=MZM(ZRHO)*XWT
!
!-------------------------------------------------------------------------------
!
!*      2.   Vertical interpolations to LES vertical grid
!            --------------------------------------------
!
!* note that velocity fields are first localized on the MASS points
!
!
IF (CRAD /= 'NONE') THEN
  CALL LES_VER_INT(   XRADEFF, ZRADEFF_LES)
  CALL LES_VER_INT(   XSWU, ZSWU_LES)
  CALL LES_VER_INT(   XSWD, ZSWD_LES)
  CALL LES_VER_INT(   XLWU, ZLWU_LES)
  CALL LES_VER_INT(   XLWD, ZLWD_LES)
  CALL LES_VER_INT(   XDTHRADSW, ZDTHRADSW_LES)
  CALL LES_VER_INT(   XDTHRADLW, ZDTHRADLW_LES)
END IF
!
CALL LES_VER_INT( XZZ   , ZZZ_LES)
CALL LES_VER_INT( XPABST, ZP_LES )
CALL LES_VER_INT( XDYP, ZDP_LES )
CALL LES_VER_INT( XTHP, ZTP_LES )
CALL LES_VER_INT( XTR, ZTR_LES )
CALL LES_VER_INT( XDISS, ZDISS_LES )
CALL LES_VER_INT( XLEM, ZLM_LES )
CALL LES_VER_INT( GZ_M_M(XPABST,XDZZ), ZDPDZ_LES )
!
CALL LES_VER_INT( MXF(XUT)  ,ZU_LES  )
CALL LES_VER_INT( MYF(XVT)  ,ZV_LES  )
CALL LES_VER_INT( MZF(XWT)  ,ZW_LES  )
CALL LES_VER_INT( MZF(ZMASSF) ,ZMF_LES)
CALL LES_VER_INT(     XTHT  ,ZTH_LES )
CALL LES_VER_INT( MXF(MZF(GZ_U_UW(XUT,XDZZ))), ZDUDZ_LES )
CALL LES_VER_INT( MYF(MZF(GZ_V_VW(XVT,XDZZ))), ZDVDZ_LES )
CALL LES_VER_INT( GZ_W_M(XWT,XDZZ), ZDWDZ_LES )
CALL LES_VER_INT( ZEXN, ZEXN_LES)  
!
CALL LES_VER_INT( GZ_M_M(XTHT,XDZZ), ZDTHDZ_LES )
!
CALL LES_VER_INT(ZRHO, ZRHO_LES)
!
IF (LUSERV) CALL LES_VER_INT(ZTHV, ZTHV_LES)
CALL LES_VER_INT(ZTHL, ZTHL_LES)
CALL LES_VER_INT( GZ_M_M(ZTHL,XDZZ), ZDTHLDZ_LES )
!
CALL LES_VER_INT(     XTKET ,ZTKE_LES)
IRR = 0
IF (LUSERV) THEN
  IRR = IRR + 1
  CALL LES_VER_INT(     XRT(:,:,:,IRR)  ,ZRV_LES )
  CALL LES_VER_INT(     ZRT(:,:,:)      ,ZRT_LES )
  CALL LES_VER_INT( GZ_M_M(ZRT,XDZZ), ZDRTDZ_LES )
  CALL LES_VER_INT(   ZREHU(:,:,:)     ,ZREHU_LES)
END IF
IF (LUSERC) THEN
  IRR = IRR + 1
  CALL LES_VER_INT(     XRT(:,:,:,IRR)  ,ZRC_LES )
  ALLOCATE(ZINDCLD (IIU,IJU,NLES_K))
  ALLOCATE(ZINDCLD2(IIU,IJU,NLES_K))
  ZINDCLD = CEILING(ZRC_LES-1.E-6)
  ZINDCLD2 = CEILING(ZRC_LES-1.E-5)
  CALL LES_VER_INT( XCLDFR(:,:,:)  ,ZCLDFR_LES )
ELSE
  ALLOCATE(ZINDCLD (0,0,0))
  ALLOCATE(ZINDCLD2(0,0,0))
END IF
IF (LUSERR) THEN
  IRR = IRR + 1
  CALL LES_VER_INT(     XRT(:,:,:,IRR)  ,ZRR_LES )
  CALL LES_VER_INT(     XINPRR3D(:,:,:), ZINPRR3D_LES)
  CALL LES_VER_INT(    XEVAP3D(:,:,:), ZEVAP3D_LES)
  CALL LES_VER_INT( XRAINFR(:,:,:)  ,ZRAINFR_LES )
END IF
IF (LUSERC) THEN
    DO JJ=1,IJU
     DO JI=1,IIU
      ZINDCLD2D(JI,JJ) = maxval(ZINDCLD(JI,JJ,:))
      ZINDCLD2D2(JI,JJ)= maxval(ZINDCLD2(JI,JJ,:))
     END DO
    END DO
  !* integration of rho rc
  !!!ZLWP_LES only for cloud water           
    ZLWP_LES(:,:)  = 0.
      DO JK=1,NLES_K-1
    ZLWP_LES(:,:) = ZLWP_LES(:,:) + (ZZZ_LES(:,:,JK+1)-ZZZ_LES(:,:,JK))      &
                      * (ZRC_LES(:,:,JK)) * ZRHO_LES(:,:,JK)
      END DO
  CALL LES_MEAN_ll ( ZLWP_LES, LLES_CURRENT_CART_MASK(:,:,1),               &
                    XLES_LWP(NLES_CURRENT_TCOUNT)     )    
!
END IF

  !!!ZRWP_LES only for rain water   
IF (LUSERR) THEN
   ZRWP_LES(:,:)=0.
  DO JK=1,NLES_K-1
    ZRWP_LES(:,:) = ZRWP_LES(:,:) + (ZZZ_LES(:,:,JK+1)-ZZZ_LES(:,:,JK))      &
                    * (ZRR_LES(:,:,JK)) * ZRHO_LES(:,:,JK)
  END DO
  CALL LES_MEAN_ll ( ZRWP_LES, LLES_CURRENT_CART_MASK(:,:,1),               &
                    XLES_RWP(NLES_CURRENT_TCOUNT)     )
ENDIF
!                    
IF (LUSERI) THEN
  IRR = IRR + 1
  CALL LES_VER_INT(     XRT(:,:,:,IRR)  ,ZRI_LES )
  ZIWP_LES(:,:)=0.
  DO JK=1,NLES_K-1
    ZIWP_LES(:,:) = ZIWP_LES(:,:) + (ZZZ_LES(:,:,JK+1)-ZZZ_LES(:,:,JK))      &
                    * (ZRI_LES(:,:,JK)) * ZRHO_LES(:,:,JK)
  END DO
  CALL LES_MEAN_ll ( ZIWP_LES, LLES_CURRENT_CART_MASK(:,:,1),               &
                    XLES_IWP(NLES_CURRENT_TCOUNT)     )
  CALL LES_VER_INT( XICEFR(:,:,:)  ,ZICEFR_LES )
END IF
IF (LUSERS) THEN
  IRR = IRR + 1
  CALL LES_VER_INT(     XRT(:,:,:,IRR)  ,ZRS_LES )
  ZSWP_LES(:,:)=0.
  DO JK=1,NLES_K-1
    ZSWP_LES(:,:) = ZSWP_LES(:,:) + (ZZZ_LES(:,:,JK+1)-ZZZ_LES(:,:,JK))      &
                    * (ZRS_LES(:,:,JK)) * ZRHO_LES(:,:,JK)
  END DO
  CALL LES_MEAN_ll ( ZSWP_LES, LLES_CURRENT_CART_MASK(:,:,1),               &
                    XLES_SWP(NLES_CURRENT_TCOUNT)     )
END IF
IF (LUSERG) THEN
  IRR = IRR + 1
  CALL LES_VER_INT(     XRT(:,:,:,IRR)  ,ZRG_LES )
  ZGWP_LES(:,:)=0.
  DO JK=1,NLES_K-1
    ZGWP_LES(:,:) = ZGWP_LES(:,:) + (ZZZ_LES(:,:,JK+1)-ZZZ_LES(:,:,JK))      &
                    * (ZRG_LES(:,:,JK)) * ZRHO_LES(:,:,JK)
  END DO
  CALL LES_MEAN_ll ( ZGWP_LES, LLES_CURRENT_CART_MASK(:,:,1),               &
                    XLES_GWP(NLES_CURRENT_TCOUNT)     )
END IF
IF (LUSERH) THEN
  IRR = IRR + 1
  CALL LES_VER_INT(     XRT(:,:,:,IRR)  ,ZRH_LES )
  ZHWP_LES(:,:)=0.
  DO JK=1,NLES_K-1
    ZHWP_LES(:,:) = ZHWP_LES(:,:) + (ZZZ_LES(:,:,JK+1)-ZZZ_LES(:,:,JK))      &
                    * (ZRH_LES(:,:,JK)) * ZRHO_LES(:,:,JK)
  END DO
  CALL LES_MEAN_ll ( ZHWP_LES, LLES_CURRENT_CART_MASK(:,:,1),               &
                    XLES_HWP(NLES_CURRENT_TCOUNT)     )
END IF
IF (NSV>0) THEN
  DO JSV=1,NSV
    CALL LES_VER_INT(  XSVT(:,:,:,JSV), ZSV_LES(:,:,:,JSV) )
    CALL LES_VER_INT( GZ_M_M(XSVT(:,:,:,JSV),XDZZ), ZDSVDZ_LES(:,:,:,JSV) )
  END DO
END IF
!
!*mean sw and lw fluxes  
  CALL LES_MEAN_ll ( ZSWU_LES, LLES_CURRENT_CART_MASK,               &
                    XLES_SWU(:,NLES_CURRENT_TCOUNT)     )
  CALL LES_MEAN_ll ( ZSWD_LES, LLES_CURRENT_CART_MASK,               &
                    XLES_SWD(:,NLES_CURRENT_TCOUNT)     )
  CALL LES_MEAN_ll ( ZLWU_LES, LLES_CURRENT_CART_MASK,               &
                    XLES_LWU(:,NLES_CURRENT_TCOUNT)     )
  CALL LES_MEAN_ll ( ZLWD_LES, LLES_CURRENT_CART_MASK,               &
                    XLES_LWD(:,NLES_CURRENT_TCOUNT)     )
  CALL LES_MEAN_ll ( ZDTHRADSW_LES, LLES_CURRENT_CART_MASK,          &
                    XLES_DTHRADSW(:,NLES_CURRENT_TCOUNT)     )
  CALL LES_MEAN_ll ( ZDTHRADLW_LES, LLES_CURRENT_CART_MASK,          &
                    XLES_DTHRADLW(:,NLES_CURRENT_TCOUNT)     )
  CALL LES_MEAN_ll ( ZRADEFF_LES, LLES_CURRENT_CART_MASK,             &
                    XLES_RADEFF(:,NLES_CURRENT_TCOUNT)     )
!* mean vertical profiles on the LES grid
!
  CALL LES_MEAN_ll ( ZU_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_U(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZV_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_V(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZW_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_W(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZP_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_P(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZDP_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_DP(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZTP_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_TP(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZTR_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_TR(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZDISS_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_DISS(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZLM_LES, LLES_CURRENT_CART_MASK,                   &
                     XLES_MEAN_LM(:,NLES_CURRENT_TCOUNT,1)     )
!
  CALL LES_MEAN_ll ( ZRHO_LES, LLES_CURRENT_CART_MASK,                 &
                     XLES_MEAN_RHO(:,NLES_CURRENT_TCOUNT,1)   )
!
  CALL LES_MEAN_ll ( ZMF_LES, LLES_CURRENT_CART_MASK,                  &
                     XLES_MEAN_Mf(:,NLES_CURRENT_TCOUNT,1)     )                             
!
  CALL LES_MEAN_ll ( ZTH_LES*ZEXN_LES, LLES_CURRENT_CART_MASK,         &
                     ZWORK1DT(:)                                 )                    
!                     
!computation of es
  ZWORK1D(:)=EXP(XALPW -                                          &
                  XBETAW/ZWORK1DT(:)         &
                  -XGAMW*ALOG(ZWORK1DT(:)))
!computation of qs

  IF (LUSERV)                                                          &
  XLES_MEAN_Qs(:,NLES_CURRENT_TCOUNT,1)=XRD/XRV*ZWORK1D(:)/        &
  (XLES_MEAN_P(:,NLES_CURRENT_TCOUNT,1)-ZWORK1D(:)*(1-XRD/XRV))   
! qs is determined from the temperature average over the current_mask
!
  CALL LES_MEAN_ll ( ZTH_LES, LLES_CURRENT_CART_MASK,                  &
                     XLES_MEAN_Th(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERV)                                                          &
  CALL LES_MEAN_ll ( ZTHV_LES, LLES_CURRENT_CART_MASK,                 &
                     XLES_MEAN_Thv(:,NLES_CURRENT_TCOUNT,1)    )
!
  IF (LUSERC)                                                       &
  CALL LES_MEAN_ll ( ZTHL_LES, LLES_CURRENT_CART_MASK,              &
                     XLES_MEAN_Thl(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERC)                                                       &
  CALL LES_MEAN_ll ( ZRT_LES, LLES_CURRENT_CART_MASK,               &
                     XLES_MEAN_Rt(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERV) &
  CALL LES_MEAN_ll ( ZRV_LES, LLES_CURRENT_CART_MASK,               &
                     XLES_MEAN_Rv(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERV) &
  CALL LES_MEAN_ll ( ZREHU_LES, LLES_CURRENT_CART_MASK,             &
                   XLES_MEAN_Rehu(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERC) &
  CALL LES_MEAN_ll ( ZRC_LES, LLES_CURRENT_CART_MASK,               &
                    XLES_MEAN_Rc(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERC) THEN 
     CALL LES_MEAN_ll ( ZINDCLD, LLES_CURRENT_CART_MASK,               &
                    XLES_MEAN_INDCf(:,NLES_CURRENT_TCOUNT,1)     )
     CALL LES_MEAN_ll ( ZINDCLD2, LLES_CURRENT_CART_MASK,              &
                    XLES_MEAN_INDCf2(:,NLES_CURRENT_TCOUNT,1)    )
     CALL LES_MEAN_ll ( ZCLDFR_LES, LLES_CURRENT_CART_MASK,            &
                    XLES_MEAN_Cf(:,NLES_CURRENT_TCOUNT,1)        )
!
!* cf total
     CALL LES_MEAN_ll( ZINDCLD2D, LLES_CURRENT_CART_MASK(:,:,1) ,      & 
                   XLES_CFtot(NLES_CURRENT_TCOUNT)      ) 
     CALL LES_MEAN_ll( ZINDCLD2D2, LLES_CURRENT_CART_MASK(:,:,1),      &
                    XLES_CF2tot(NLES_CURRENT_TCOUNT)    )
  ENDIF
!
  IF (LUSERR) THEN

    CALL LES_MEAN_ll ( XINPRR, LLES_CURRENT_CART_MASK(:,:,1),    &
                     XLES_INPRR(NLES_CURRENT_TCOUNT)    )
    ZINPRRm=0.
    ZCOUNT=0.
    ZINDCLD2D(:,:)=0.
    DO JJ=1,IJU
     DO JI=1,IIU
       IF (ZRR_LES(JI,JJ,1) .GT. 1.E-6) ZINPRRm = ZINPRRm+XINPRR(JI,JJ)
       IF (ZRR_LES(JI,JJ,1) .GT. 1.E-6) ZINDCLD2D(JI,JJ)=1.
       IF (ZRR_LES(JI,JJ,1) .GT. 1.E-6) ZCOUNT=ZCOUNT+1.
     END DO
    END DO                 
    IF (ZCOUNT .GE. 1) ZINPRRm=ZINPRRm/ZCOUNT
    XLES_RAIN_INPRR(NLES_CURRENT_TCOUNT)=ZINPRRm
    CALL  LES_MEAN_ll ( ZINDCLD2D, LLES_CURRENT_CART_MASK(:,:,1),  &
                     XLES_PRECFR(NLES_CURRENT_TCOUNT)         )   
    CALL LES_MEAN_ll ( ZINPRR3D_LES, LLES_CURRENT_CART_MASK,  &
                     XLES_INPRR3D(:,NLES_CURRENT_TCOUNT,1)  )
    CALL LES_MEAN_ll ( ZEVAP3D_LES, LLES_CURRENT_CART_MASK,  &
                     XLES_EVAP3D(:,NLES_CURRENT_TCOUNT,1)   )         
    DO JK=1,NLES_K
     CHAMPXY1(:,:,1)=ZINPRR3D_LES(:,:,JK)                
     XLES_MAX_INPRR3D(JK,NLES_CURRENT_TCOUNT,1)=MAX_ll (CHAMPXY1,IINFO_ll, &
                         IIA_ll,IJA_ll,1,IIU_ll,IJU_ll,1)
    END DO
!

!   conversion de m/s en mm/day
 XLES_RAIN_INPRR(NLES_CURRENT_TCOUNT)=XLES_RAIN_INPRR(NLES_CURRENT_TCOUNT)*3.6E6*24.
    XLES_INPRR(NLES_CURRENT_TCOUNT)=XLES_INPRR(NLES_CURRENT_TCOUNT)*3.6E6*24.

    CALL LES_MEAN_ll ( XACPRR, LLES_CURRENT_CART_MASK(:,:,1),    &
                     XLES_ACPRR(NLES_CURRENT_TCOUNT)    )
!   conversion de m en mm
    XLES_ACPRR(NLES_CURRENT_TCOUNT)=XLES_ACPRR(NLES_CURRENT_TCOUNT)*1000.
    CALL LES_MEAN_ll ( ZRAINFR_LES, LLES_CURRENT_CART_MASK,            &
                    XLES_MEAN_RF(:,NLES_CURRENT_TCOUNT,1)        )

  ENDIF
!
  IF (LUSERC ) THEN
    IF (( CCLOUD(1:3) == 'ICE'                                   .AND.LSEDIC) .OR. &
    ((CCLOUD=='C2R2' .OR. CCLOUD=='C3R5' .OR. CCLOUD=='KHKO').AND.LSEDC)  .OR. &
    ( CCLOUD=='LIMA'                                         .AND.MSEDC))  THEN
      CALL LES_MEAN_ll ( XINPRC, LLES_CURRENT_CART_MASK(:,:,1),    &
                     XLES_INPRC(NLES_CURRENT_TCOUNT)    )
!     conversion from m/s to mm/day
      XLES_INPRC(NLES_CURRENT_TCOUNT)=XLES_INPRC(NLES_CURRENT_TCOUNT)*3.6E6*24.
    ENDIF
    IF ( (((CCLOUD == 'KHKO') .OR.(CCLOUD == 'C2R2')) .AND. LDEPOC) &
         .OR. ( (CCLOUD(1:3) == 'ICE') .AND. LDEPOSC) ) THEN
        CALL LES_MEAN_ll ( XINDEP, LLES_CURRENT_CART_MASK(:,:,1),    &
                     XLES_INDEP(NLES_CURRENT_TCOUNT)    )
!       conversion from m/s to mm/day
        XLES_INDEP(NLES_CURRENT_TCOUNT)=XLES_INDEP(NLES_CURRENT_TCOUNT)*3.6E6*24.
    ENDIF
  ENDIF
!
  IF (LUSERR) &
  CALL LES_MEAN_ll ( ZRR_LES, LLES_CURRENT_CART_MASK,               &
                     XLES_MEAN_Rr(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERI) &
  CALL LES_MEAN_ll ( ZRI_LES, LLES_CURRENT_CART_MASK,               &
                     XLES_MEAN_Ri(:,NLES_CURRENT_TCOUNT,1)     )
  CALL LES_MEAN_ll ( ZICEFR_LES, LLES_CURRENT_CART_MASK,            &
                    XLES_MEAN_If(:,NLES_CURRENT_TCOUNT,1)        )
!
  IF (LUSERS) &
  CALL LES_MEAN_ll ( ZRS_LES, LLES_CURRENT_CART_MASK,               &
                   XLES_MEAN_Rs(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERG) &
  CALL LES_MEAN_ll ( ZRG_LES, LLES_CURRENT_CART_MASK,               &
                   XLES_MEAN_Rg(:,NLES_CURRENT_TCOUNT,1)     )
!
  IF (LUSERH) &
  CALL LES_MEAN_ll ( ZRH_LES, LLES_CURRENT_CART_MASK,               &
                     XLES_MEAN_Rh(:,NLES_CURRENT_TCOUNT,1)   )
!
  DO JSV=1,NSV
    CALL LES_MEAN_ll ( ZSV_LES(:,:,:,JSV), LLES_CURRENT_CART_MASK,   &
                       XLES_MEAN_Sv(:,NLES_CURRENT_TCOUNT,1,JSV)     )
  END DO
!
  CALL LES_MEAN_ll ( ZDPDZ_LES, LLES_CURRENT_CART_MASK,      &
                     ZMEAN_DPDZ(:)                           )
  CALL LES_MEAN_ll ( ZDTHDZ_LES, LLES_CURRENT_CART_MASK,      &
                     ZLES_MEAN_DTHDZ(:)                      )
                     
!
!* build the 3D resolved turbulent fields by removing the mean field
!
DO JJ=1,IJU
  DO JI=1,IIU
    ZP_ANOM(JI,JJ,:) = ZP_LES(JI,JJ,:) - XLES_MEAN_P(:,NLES_CURRENT_TCOUNT,1)
    ZDPDZ_ANOM(JI,JJ,:) = ZDPDZ_LES(JI,JJ,:) - ZMEAN_DPDZ(:)
    ZTH_ANOM(JI,JJ,:) = ZTH_LES(JI,JJ,:) - XLES_MEAN_Th(:,NLES_CURRENT_TCOUNT,1)
    ZRHO_ANOM(JI,JJ,:) = ZRHO_LES(JI,JJ,:) - XLES_MEAN_Rho(:,NLES_CURRENT_TCOUNT,1)
    IF (LUSERV) THEN
      ZTHV_ANOM(JI,JJ,:) = ZTHV_LES(JI,JJ,:) - XLES_MEAN_Thv(:,NLES_CURRENT_TCOUNT,1)
      ZRV_ANOM(JI,JJ,:) = ZRV_LES(JI,JJ,:) - XLES_MEAN_Rv(:,NLES_CURRENT_TCOUNT,1)
    END IF
    IF (LUSERC) THEN
      ZRC_ANOM(JI,JJ,:) = ZRC_LES(JI,JJ,:) - XLES_MEAN_Rc(:,NLES_CURRENT_TCOUNT,1)
      ZLWP_ANOM(JI,JJ) =ZLWP_LES(JI,JJ)-XLES_LWP(NLES_CURRENT_TCOUNT)
    END IF
    IF (LUSERI) THEN
      ZRI_ANOM(JI,JJ,:) = ZRI_LES(JI,JJ,:) - XLES_MEAN_Ri(:,NLES_CURRENT_TCOUNT,1)
    END IF
    IF (LUSERR) THEN
      ZRR_ANOM(JI,JJ,:) = ZRR_LES(JI,JJ,:) - XLES_MEAN_Rr(:,NLES_CURRENT_TCOUNT,1)
    END IF
  END DO
END DO
!
!
!--------------------------------------------------------------------------------
!
!* vertical grid computed at first LES call for this model
!
IF (NLES_CURRENT_TCOUNT==1) THEN
  ALLOCATE(ZZ_LES    (IIU,IJU,NLES_K))
  CALL LES_VER_INT( MZF(XZZ)  ,ZZ_LES   )
  CALL LES_MEAN_ll ( ZZ_LES, LLES_CURRENT_CART_MASK, XLES_Z  )
  DEALLOCATE(ZZ_LES)
  CALL LES_MEAN_ll ( XZS,    LLES_CURRENT_CART_MASK(:,:,1), XLES_ZS )
END IF
!
!-------------------------------------------------------------------------------
!
!*      3.   Vertical interpolations to SECTRA computations vertical grid
!            ------------------------------------------------------------
!
!* note that velocity fields are previously localized on the MASS points
!
CALL SPEC_VER_INT(IMI, MXF(XUT)  ,ZU_SPEC  )
CALL SPEC_VER_INT(IMI, MYF(XVT)  ,ZV_SPEC  )
CALL SPEC_VER_INT(IMI, MZF(XWT)  ,ZW_SPEC  )
CALL SPEC_VER_INT(IMI,     XTHT  ,ZTH_SPEC )
IF (LUSERC) CALL SPEC_VER_INT(IMI,     ZTHL  ,ZTHL_SPEC)
IRR = 0
IF (LUSERV) THEN
  IRR = IRR + 1
  CALL SPEC_VER_INT(IMI,     XRT(:,:,:,IRR)  ,ZRV_SPEC )
END IF
IF (LUSERC) THEN
  IRR = IRR + 1
  CALL SPEC_VER_INT(IMI,     XRT(:,:,:,IRR)  ,ZRC_SPEC )
END IF
IF (LUSERR) THEN
  IRR = IRR + 1
END IF
IF (LUSERI) THEN
  IRR = IRR + 1
  CALL SPEC_VER_INT(IMI,     XRT(:,:,:,IRR)  ,ZRI_SPEC )
END IF
IF (NSV>0) THEN
  DO JSV=1,NSV
    CALL SPEC_VER_INT(IMI,  XSVT(:,:,:,JSV), ZSV_SPEC(:,:,:,JSV) )
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*      4.   Call to LES computations on cartesian (sub-)domain
!            --------------------------------------------------
!
IMASK=1
!
CALL LES(LLES_CURRENT_CART_MASK)
!
!-------------------------------------------------------------------------------
!
!*      5.   Call to LES computations on nebulosity mask
!            -------------------------------------------
!
IF (LLES_NEB_MASK) THEN
  IMASK=IMASK+1
  CALL LES(LLES_CURRENT_NEB_MASK .AND. LLES_CURRENT_CART_MASK)
!
  IMASK=IMASK+1
  CALL LES((.NOT. LLES_CURRENT_NEB_MASK) .AND. LLES_CURRENT_CART_MASK)
END IF
!
!-------------------------------------------------------------------------------
!
!*      6.   Call to LES computations on cloud core mask
!            -------------------------------------------
!
IF (LLES_CORE_MASK) THEN
  IMASK=IMASK+1
  CALL LES(LLES_CURRENT_CORE_MASK .AND. LLES_CURRENT_CART_MASK)
!
  IMASK=IMASK+1
  CALL LES((.NOT. LLES_CURRENT_CORE_MASK) .AND. LLES_CURRENT_CART_MASK)
END IF
!
!-------------------------------------------------------------------------------
!
!*      7.   Call to LES computations on user mask
!            -------------------------------------
!
IF (LLES_MY_MASK) THEN
 DO JI=1,NLES_MASKS_USER
  IMASK=IMASK+1
  CALL LES(LLES_CURRENT_MY_MASKS(:,:,:,JI))
 END DO
END IF
!
!-------------------------------------------------------------------------------
!
!*      7b.  Call to LES computations on conditional sampling mask
!            -----------------------------------------------------
!
IF (LLES_CS_MASK) THEN
  IMASK=IMASK+1
  CALL LES(LLES_CURRENT_CS1_MASK)
  IMASK=IMASK+1
  CALL LES(LLES_CURRENT_CS2_MASK)
  IMASK=IMASK+1
  CALL LES(LLES_CURRENT_CS3_MASK)
END IF
!
!-------------------------------------------------------------------------------
!
!*      8.   budgets
!            -------
!
!*      8.1  tendencies
!            ----------
!
!
!*      8.2  dynamical production, transport and mean advection
!            --------------------------------------------------
!
ALLOCATE(ZLES_MEAN_DRtDZ(NLES_K))
ALLOCATE(ZLES_MEAN_DSVDZ(NLES_K,NSV))
!
IF (LUSERV) THEN
  ZLES_MEAN_DRtDZ(:) = XLES_MEAN_DRtDZ(:,NLES_CURRENT_TCOUNT,1)
ELSE
  ZLES_MEAN_DRtDZ(:) = XUNDEF
END IF
!
ZLES_MEAN_DSVDZ = 0.
DO JSV=1,NSV
  ZLES_MEAN_DSvDZ(:,JSV) = XLES_MEAN_DSvDZ(:,NLES_CURRENT_TCOUNT,1,JSV)
END DO
!
CALL LES_RES_TR(LUSERV,                                    &
                XLES_MEAN_DUDZ(:,NLES_CURRENT_TCOUNT,1),   &
                XLES_MEAN_DVDZ(:,NLES_CURRENT_TCOUNT,1),   &
                XLES_MEAN_DWDZ(:,NLES_CURRENT_TCOUNT,1),   &
                XLES_MEAN_DThlDZ(:,NLES_CURRENT_TCOUNT,1), &
                ZLES_MEAN_DRtDZ(:),                        &
                ZLES_MEAN_DSvDZ(:,:)                       )
!
DEALLOCATE(ZLES_MEAN_DRtDZ)
DEALLOCATE(ZLES_MEAN_DSVDZ)
!
CALL LES_BUDGET_TEND_n
!*      8.3  end of LES budgets computations
!            -------------------------------
!
DO JLOOP=1,NLES_TOT
  XLES_BU_RES_KE   (:,NLES_CURRENT_TCOUNT,JLOOP) = X_LES_BU_RES_KE   (:,JLOOP)
  XLES_BU_RES_WThl (:,NLES_CURRENT_TCOUNT,JLOOP) = X_LES_BU_RES_WThl (:,JLOOP)
  XLES_BU_RES_Thl2 (:,NLES_CURRENT_TCOUNT,JLOOP) = X_LES_BU_RES_Thl2 (:,JLOOP)
  XLES_BU_SBG_Tke  (:,NLES_CURRENT_TCOUNT,JLOOP) = X_LES_BU_SBG_Tke  (:,JLOOP)
  IF (LUSERV) THEN
    XLES_BU_RES_WRt  (:,NLES_CURRENT_TCOUNT,JLOOP) = X_LES_BU_RES_WRt  (:,JLOOP)
    XLES_BU_RES_Rt2  (:,NLES_CURRENT_TCOUNT,JLOOP) = X_LES_BU_RES_Rt2  (:,JLOOP)
    XLES_BU_RES_ThlRt(:,NLES_CURRENT_TCOUNT,JLOOP) = X_LES_BU_RES_ThlRt(:,JLOOP)
  END IF
  DO JSV=1,NSV
    XLES_BU_RES_Sv2  (:,NLES_CURRENT_TCOUNT,JLOOP,JSV) = X_LES_BU_RES_Sv2  (:,JLOOP,JSV)
    XLES_BU_RES_WSv  (:,NLES_CURRENT_TCOUNT,JLOOP,JSV) = X_LES_BU_RES_WSv  (:,JLOOP,JSV)
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
!*      9.   Deallocations
!            -------------
!
!*      9.1  local variables
!            ---------------
!
DEALLOCATE(ZEXN      )
DEALLOCATE(ZTHL)
DEALLOCATE(ZRT       )
DEALLOCATE(ZTHV      )
DEALLOCATE(ZRHO      )
DEALLOCATE(ZEW       )

DEALLOCATE(ZINDCLD   )
DEALLOCATE(ZINDCLD2  )
DEALLOCATE(ZINDCLD2D )
DEALLOCATE(ZINDCLD2D2)
DEALLOCATE(ZCLDFR_LES)
DEALLOCATE(ZICEFR_LES)
DEALLOCATE(ZRAINFR_LES)
DEALLOCATE(ZMASSF    )  
DEALLOCATE(ZTEMP     )
DEALLOCATE(ZREHU     )
DEALLOCATE(CHAMPXY1  )
!
DEALLOCATE(ZU_LES)
DEALLOCATE(ZV_LES)
DEALLOCATE(ZW_LES)
DEALLOCATE(ZTHL_LES)
DEALLOCATE(ZRT_LES)
DEALLOCATE(ZSV_LES)
DEALLOCATE(ZP_LES   )
DEALLOCATE(ZDP_LES   )
DEALLOCATE(ZTP_LES   )
DEALLOCATE(ZTR_LES   )
DEALLOCATE(ZDISS_LES   )
DEALLOCATE(ZLM_LES   )
DEALLOCATE(ZDPDZ_LES)
DEALLOCATE(ZLWP_ANOM)
DEALLOCATE(ZWORK2D)
DEALLOCATE(ZWORK1D)
DEALLOCATE(ZWORK1DT)
DEALLOCATE(ZMAXWRR2D)
DEALLOCATE(ZDTHLDZ_LES)
DEALLOCATE(ZDTHDZ_LES)
DEALLOCATE(ZDRTDZ_LES)
DEALLOCATE(ZDSVDZ_LES)
DEALLOCATE(ZDUDZ_LES)
DEALLOCATE(ZDVDZ_LES)
DEALLOCATE(ZDWDZ_LES)
DEALLOCATE(ZRHO_LES )
DEALLOCATE(ZEXN_LES )
DEALLOCATE(ZTH_LES  )
DEALLOCATE(ZMF_LES  )
DEALLOCATE(ZTHV_LES )
DEALLOCATE(ZTKE_LES )
DEALLOCATE(ZKE_LES )
DEALLOCATE(ZTKET_LES)
DEALLOCATE(ZRV_LES  )
DEALLOCATE(ZREHU_LES  )
DEALLOCATE(ZRC_LES  )
DEALLOCATE(ZRR_LES  )
DEALLOCATE(ZZZ_LES)
DEALLOCATE(ZLWP_LES )
DEALLOCATE(ZRWP_LES )
DEALLOCATE(ZIWP_LES )
DEALLOCATE(ZSWP_LES )
DEALLOCATE(ZGWP_LES )
DEALLOCATE(ZHWP_LES )
DEALLOCATE(ZINPRR3D_LES)
DEALLOCATE(ZEVAP3D_LES)
DEALLOCATE(ZRI_LES  )
DEALLOCATE(ZRS_LES  )
DEALLOCATE(ZRG_LES  )
DEALLOCATE(ZRH_LES  )
DEALLOCATE(ZP_ANOM  )
DEALLOCATE(ZRHO_ANOM)
DEALLOCATE(ZTH_ANOM )
DEALLOCATE(ZTHV_ANOM)
DEALLOCATE(ZRV_ANOM )
DEALLOCATE(ZRC_ANOM )
DEALLOCATE(ZRI_ANOM )
DEALLOCATE(ZRR_ANOM )
DEALLOCATE(ZDPDZ_ANOM)
DEALLOCATE(ZMEAN_DPDZ)
DEALLOCATE(ZLES_MEAN_DTHDZ)
!
DEALLOCATE(ZU_SPEC   )
DEALLOCATE(ZV_SPEC   )
DEALLOCATE(ZW_SPEC   )
DEALLOCATE(ZTH_SPEC  )
DEALLOCATE(ZTHL_SPEC )
DEALLOCATE(ZRV_SPEC  )
DEALLOCATE(ZRC_SPEC  )
DEALLOCATE(ZRI_SPEC  )
DEALLOCATE(ZSV_SPEC  )
!
DEALLOCATE(ZRADEFF_LES  )
DEALLOCATE(ZSWU_LES  )
DEALLOCATE(ZSWD_LES  )
DEALLOCATE(ZLWD_LES  )
DEALLOCATE(ZLWU_LES  )
DEALLOCATE(ZDTHRADSW_LES  )
DEALLOCATE(ZDTHRADLW_LES  )
!
!*      9.2  current time-step LES masks (in MODD_LES)
!            ---------------------------
!
CALL LES_DEALLOCATE('LLES_CURRENT_CART_MASK')
IF (LLES_NEB_MASK)     CALL LES_DEALLOCATE('LLES_CURRENT_NEB_MASK')
IF (LLES_CORE_MASK)    CALL LES_DEALLOCATE('LLES_CURRENT_CORE_MASK')
IF (LLES_MY_MASK)   THEN
   CALL LES_DEALLOCATE('LLES_CURRENT_MY_MASKS')
END IF
IF (LLES_CS_MASK)   THEN
    CALL LES_DEALLOCATE('LLES_CURRENT_CS1_MASK')
    IF (NSV_CS >= 2) CALL LES_DEALLOCATE('LLES_CURRENT_CS2_MASK')
    IF (NSV_CS == 3) CALL LES_DEALLOCATE('LLES_CURRENT_CS3_MASK')
END IF
!
!
!*      9.3  variables in MODD_LES_BUDGET
!            ----------------------------
!

DEALLOCATE(XU_ANOM  )
DEALLOCATE(XV_ANOM  )
DEALLOCATE(XW_ANOM  )
DEALLOCATE(XTHL_ANOM)
DEALLOCATE(XRT_ANOM )
DEALLOCATE(XSV_ANOM )
!
DEALLOCATE(XCURRENT_L_O_EXN_CP)
DEALLOCATE(XCURRENT_RHODJ     )
!
DEALLOCATE(XCURRENT_RUS  )
DEALLOCATE(XCURRENT_RVS  )
DEALLOCATE(XCURRENT_RWS  )
DEALLOCATE(XCURRENT_RTHS )
DEALLOCATE(XCURRENT_RTKES)
DEALLOCATE(XCURRENT_RRS  )
DEALLOCATE(XCURRENT_RSVS )
DEALLOCATE(XCURRENT_RTHLS)
DEALLOCATE(XCURRENT_RRTS )

DEALLOCATE(X_LES_BU_RES_KE   )
DEALLOCATE(X_LES_BU_RES_WThl )
DEALLOCATE(X_LES_BU_RES_Thl2 )
DEALLOCATE(X_LES_BU_RES_WRt  )
DEALLOCATE(X_LES_BU_RES_Rt2  )
DEALLOCATE(X_LES_BU_RES_ThlRt)
DEALLOCATE(X_LES_BU_RES_Sv2  )
DEALLOCATE(X_LES_BU_RES_WSv  )
DEALLOCATE(X_LES_BU_SBG_TKE  )
!-------------------------------------------------------------------------------
!
!*      10.  end of LES computations for this time-step
!            ------------------------------------------
!
LLES_CALL=.FALSE.
CALL BUDGET_FLAGS(LUSERV, LUSERC, LUSERR,         &
                  LUSERI, LUSERS, LUSERG, LUSERH  )
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!     ##########################################################################
      SUBROUTINE  LES(OMASK)
!     ##########################################################################
!
!
!!****  *LES* computes the current time-step LES diagnostics for one mask.
!!
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_PARAMETERS
!
USE MODI_LES_FLUX_ll
USE MODI_LES_3RD_MOMENT_ll
USE MODI_LES_4TH_MOMENT_ll
USE MODI_LES_MEAN_1PROC
USE MODI_LES_MEAN_MPROC
USE MODI_LES_PDF_ll
!
USE MODI_LES_HOR_CORR
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
LOGICAL, DIMENSION(:,:,:),   INTENT(IN) :: OMASK     ! 2D mask for computations
!
!
!
!       0.2  declaration of local variables
!
INTEGER :: JSV      ! scalar variables counter
INTEGER :: JI
INTEGER :: JK       ! vertical loop counter
INTEGER :: JPDF     ! pdf counter
!
LOGICAL, DIMENSION(SIZE(ZW_LES,1),SIZE(ZW_LES,2),SIZE(ZW_LES,3)) :: GUPDRAFT_MASK
LOGICAL, DIMENSION(SIZE(ZW_LES,1),SIZE(ZW_LES,2),SIZE(ZW_LES,3)) :: GDOWNDRAFT_MASK
REAL,    DIMENSION(SIZE(ZW_LES,1),SIZE(ZW_LES,2),SIZE(ZW_LES,3)) :: ZUPDRAFT
REAL,    DIMENSION(SIZE(ZW_LES,1),SIZE(ZW_LES,2),SIZE(ZW_LES,3)) :: ZDOWNDRAFT
REAL,    DIMENSION(SIZE(ZW_LES,1),SIZE(ZW_LES,2),SIZE(ZW_LES,3)) :: ZW_UP
REAL,    DIMENSION(SIZE(ZW_LES,1),SIZE(ZW_LES,2),SIZE(ZW_LES,3)) :: ZWORK_LES
!
INTEGER, DIMENSION(SIZE(ZW_LES,3)) :: IAVG_PTS
INTEGER, DIMENSION(SIZE(ZW_LES,3)) :: IUND_PTS
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZAVG
!
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_RESOLVED_U3
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_RESOLVED_UV2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_RESOLVED_UW2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_RESOLVED_VU2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_RESOLVED_V3
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_RESOLVED_VW2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_RESOLVED_WU2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_RESOLVED_WV2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_UPDRAFT_U2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_UPDRAFT_V2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_UPDRAFT_W2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_DOWNDRAFT_U2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_DOWNDRAFT_V2
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZLES_DOWNDRAFT_W2
REAL,    DIMENSION(SIZE(ZW_LES,3),NPDF) :: ZPDF
!
INTEGER, DIMENSION(1)               :: IKMIN_FLUX ! vertical index of min. W'thl'
INTEGER, DIMENSION(1)               :: IKMAX_TH !vertical index maxdth
INTEGER, DIMENSION(1)               :: IKMAX_CF   ! vertical index of max. Cf
!
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZKE_TOT    ! total turbulent kinetic energy
REAL                               :: ZINT_KE_TOT! integral of KE_TOT
REAL                               :: ZINT_RHOKE! integral of RHO*KE
REAL                               :: ZFRIC_SURF ! surface friction
REAL,    DIMENSION(SIZE(ZW_LES,3)) :: ZFRIC_LES  ! friction at all LES levels
!
!-------------------------------------------------------------------------------
!
!       1.   local diagnostics (for any mask type)
!            -----------------
!
!
!       1.2  Number of points used for averaging on current processor
!            --------------------------------------------------------
!
!* to be sure to be coherent with other computations,
!  a field on LES vertical grid (and horizontal mass point grid) is used.
!  This information is necessary for the subgrid fluxes computations, because
!  half of the work is already done, but the number of averaging points was
!  not kept.
!
CALL LES_MEAN_1PROC ( XW_ANOM, OMASK,                              &
                      ZAVG(:),                                     &
                      IAVG_PTS(:),                                 &
                      IUND_PTS(:)                                  )
!
!
!       1.3  Number of points used for averaging on all processor
!            ----------------------------------------------------
!
CALL LES_MEAN_ll ( XW_ANOM, OMASK,                                  &
                   ZAVG(:),                                         &
                   NLES_AVG_PTS_ll(:,NLES_CURRENT_TCOUNT,IMASK),    &
                   NLES_UND_PTS_ll(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!
!       1.4  Mean quantities
!            ---------------
!
IF (LLES_MEAN .AND. IMASK > 1) THEN
!
!* horizontal wind velocities
!
  CALL LES_MEAN_ll ( ZU_LES, OMASK,                               &
                     XLES_MEAN_U(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
  CALL LES_MEAN_ll ( ZV_LES, OMASK,                               &
                     XLES_MEAN_V(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* vertical wind velocity
!
  CALL LES_MEAN_ll ( ZW_LES, OMASK,                               &
                     XLES_MEAN_W(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* pressure
!
  CALL LES_MEAN_ll ( ZP_LES, OMASK,                               &
                     XLES_MEAN_P(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* dynamical production TKE
!
  CALL LES_MEAN_ll ( ZDP_LES, OMASK,                               &
                     XLES_MEAN_DP(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* thermal production TKE
!
  CALL LES_MEAN_ll ( ZTP_LES, OMASK,                               &
                     XLES_MEAN_TP(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* transport TKE
!
  CALL LES_MEAN_ll ( ZTR_LES, OMASK,                               &
                     XLES_MEAN_TR(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* dissipation TKE
!
  CALL LES_MEAN_ll ( ZDISS_LES, OMASK,                               &
                     XLES_MEAN_DISS(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* mixing length            
!
  CALL LES_MEAN_ll ( ZLM_LES, OMASK,                               &
                     XLES_MEAN_LM(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* density
!
  CALL LES_MEAN_ll ( ZRHO_LES, OMASK,                             &
                     XLES_MEAN_RHO(:,NLES_CURRENT_TCOUNT,IMASK)   )
!
!
!* potential temperature
!
  CALL LES_MEAN_ll ( ZTH_LES, OMASK,                               &
                     XLES_MEAN_Th(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* mass flux
  CALL LES_MEAN_ll ( ZMF_LES, OMASK,                               &
                     XLES_MEAN_Mf(:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!
!* virtual potential temperature
!
  IF (LUSERV)                                                      &
  CALL LES_MEAN_ll ( ZTHV_LES, OMASK,                              &
                     XLES_MEAN_Thv(:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* liquid potential temperature
!
  IF (LUSERC) THEN
    CALL LES_MEAN_ll ( ZTHL_LES, OMASK,                               &
                       XLES_MEAN_Thl(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!* vapor mixing ratio
!
  IF (LUSERV) THEN
    CALL LES_MEAN_ll ( ZRV_LES, OMASK,                               &
                       XLES_MEAN_Rv(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!*relative humidity
!
  IF (LUSERV) THEN
    CALL LES_MEAN_ll ( ZREHU_LES, OMASK,                               &
                       XLES_MEAN_Rehu(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!  
!* cloud mixing ratio
!
  IF (LUSERC) THEN
    CALL LES_MEAN_ll ( ZRC_LES, OMASK,                              &
                      XLES_MEAN_Rc(:,NLES_CURRENT_TCOUNT,IMASK)     )
    CALL LES_MEAN_ll ( ZRT_LES, OMASK,                               &
                      XLES_MEAN_Rt(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!* rain mixing ratio
!
  IF (LUSERR) THEN
    CALL LES_MEAN_ll ( ZRR_LES, OMASK,                               &
                       XLES_MEAN_Rr(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!* ice mixing ratio
!
  IF (LUSERI) THEN
    CALL LES_MEAN_ll ( ZRI_LES, OMASK,                               &
                       XLES_MEAN_Ri(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!* snow mixing ratio
!
  IF (LUSERS) THEN
    CALL LES_MEAN_ll ( ZRS_LES, OMASK,                             &
                     XLES_MEAN_Rs(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!* graupel mixing ratio
!
  IF (LUSERG) THEN
    CALL LES_MEAN_ll ( ZRG_LES, OMASK,                             &
                     XLES_MEAN_Rg(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!* hail mixing ratio
!
  IF (LUSERH) THEN
    CALL LES_MEAN_ll ( ZRH_LES, OMASK,                             &
                       XLES_MEAN_Rh(:,NLES_CURRENT_TCOUNT,IMASK)   )
  END IF
!
!* scalar variables mixing ratio
!
  DO JSV=1,NSV
    CALL LES_MEAN_ll ( ZSV_LES(:,:,:,JSV), OMASK,                        &
                       XLES_MEAN_Sv(:,NLES_CURRENT_TCOUNT,IMASK,JSV)     )
  END DO
END IF
!
!* wind modulus
!
IF (LLES_MEAN) THEN
!
  ZWORK_LES =SQRT( ZU_LES**2 +ZV_LES**2 )
  CALL LES_MEAN_ll ( ZWORK_LES, OMASK,                           &
                     XLES_MEAN_WIND(:,NLES_CURRENT_TCOUNT,IMASK) )
!
!* vertical speed larger than mean vertical speed (updraft)
!
  DO JK=1,NLES_K
    ZW_UP(:,:,JK) =  MAX(ZW_LES(:,:,JK), XLES_MEAN_W(JK,NLES_CURRENT_TCOUNT,IMASK))
  END DO
!
!* upward mass flux
!
  ZWORK_LES = ZW_UP * ZRHO_LES
  CALL LES_MEAN_ll ( ZWORK_LES, OMASK,                                 &
                     XLES_RESOLVED_MASSFX(:,NLES_CURRENT_TCOUNT,IMASK) )
!
!* pdf calculation
!
 IF (LLES_PDF) THEN
  CALL LES_PDF_ll  ( ZTH_LES,OMASK,XTH_PDF_MIN,XTH_PDF_MAX,    &
                        ZPDF(:,:) )
  DO JSV=1,NPDF
       XLES_PDF_TH(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
  END DO
                        
  CALL LES_PDF_ll  ( ZW_LES,OMASK,XW_PDF_MIN,XW_PDF_MAX,    &
                        ZPDF(:,:) )
  DO JSV=1,NPDF
       XLES_PDF_W(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
  END DO
  CALL LES_PDF_ll  ( ZTHV_LES,OMASK,XTHV_PDF_MIN,XTHV_PDF_MAX,    &
                        ZPDF(:,:) )
  DO JSV=1,NPDF
       XLES_PDF_THV(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
  END DO
  IF (LUSERV) THEN
   CALL LES_PDF_ll  ( ZRV_LES,OMASK,XRV_PDF_MIN,XRV_PDF_MAX,    &
                        ZPDF(:,:) )
   DO JSV=1,NPDF
       XLES_PDF_RV(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
   END DO
  END IF
  IF (LUSERC) THEN
    CALL LES_PDF_ll  ( ZRC_LES,OMASK,XRC_PDF_MIN,XRC_PDF_MAX,    &
                        ZPDF(:,:) )
    DO JSV=1,NPDF
       XLES_PDF_RC(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
    END DO
    CALL LES_PDF_ll  ( ZRT_LES,OMASK,XRT_PDF_MIN,XRT_PDF_MAX,    &
                        ZPDF(:,:) )
    DO JSV=1,NPDF
       XLES_PDF_RT(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
    END DO
    CALL LES_PDF_ll  ( ZTHL_LES,OMASK,XTHL_PDF_MIN,XTHL_PDF_MAX,    &
                        ZPDF(:,:) )
    DO JSV=1,NPDF
       XLES_PDF_THL(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
    END DO
  END IF
  IF (LUSERR) THEN
    CALL LES_PDF_ll  ( ZRR_LES,OMASK,XRR_PDF_MIN,XRR_PDF_MAX,    &
                        ZPDF(:,:) )
    DO JSV=1,NPDF
       XLES_PDF_RR(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
    END DO
  END IF
  IF (LUSERI) THEN
    CALL LES_PDF_ll  ( ZRI_LES,OMASK,XRI_PDF_MIN,XRI_PDF_MAX,    &
                        ZPDF(:,:) )
    DO JSV=1,NPDF
       XLES_PDF_RI(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
    END DO
  END IF
  IF (LUSERS) THEN
    CALL LES_PDF_ll  ( ZRS_LES,OMASK,XRS_PDF_MIN,XRS_PDF_MAX,    &
                        ZPDF(:,:) )
    DO JSV=1,NPDF
       XLES_PDF_RS(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
    END DO
  END IF
  IF (LUSERG) THEN
    CALL LES_PDF_ll  ( ZRG_LES,OMASK,XRG_PDF_MIN,XRG_PDF_MAX,    &
                        ZPDF(:,:) )
    DO JSV=1,NPDF
       XLES_PDF_RG(:,NLES_CURRENT_TCOUNT,IMASK,JSV)=ZPDF(:,JSV)
    END DO
  END IF
 END IF
!
!* mean vertical gradients
!
  CALL LES_MEAN_ll ( ZDTHLDZ_LES, OMASK, XLES_MEAN_DTHLDZ(:,NLES_CURRENT_TCOUNT,IMASK) )
  CALL LES_MEAN_ll ( ZDUDZ_LES, OMASK, XLES_MEAN_DUDZ(:,NLES_CURRENT_TCOUNT,IMASK) )
  CALL LES_MEAN_ll ( ZDVDZ_LES, OMASK, XLES_MEAN_DVDZ(:,NLES_CURRENT_TCOUNT,IMASK) )
  CALL LES_MEAN_ll ( ZDWDZ_LES, OMASK, XLES_MEAN_DWDZ(:,NLES_CURRENT_TCOUNT,IMASK) )
  IF (LUSERV) CALL LES_MEAN_ll ( ZDRtDZ_LES, OMASK, XLES_MEAN_DRtDZ(:,NLES_CURRENT_TCOUNT,IMASK) )
  DO JSV=1,NSV
    CALL LES_MEAN_ll ( ZDSVDZ_LES(:,:,:,JSV), OMASK, XLES_MEAN_DSVDZ(:,NLES_CURRENT_TCOUNT,IMASK,JSV) )
  END DO

END IF
!-------------------------------------------------------------------------------
!
!       1.5  Resolved quantities
!            -------------------
!
!* horizontal wind variances
!
  CALL LES_FLUX_ll ( XU_ANOM, XU_ANOM,                                 &
                     OMASK,                                            &
                     XLES_RESOLVED_U2 (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
  CALL LES_FLUX_ll ( XV_ANOM, XV_ANOM,                                 &
                     OMASK,                                            &
                     XLES_RESOLVED_V2 (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* vertical wind variance
!
  CALL LES_FLUX_ll ( XW_ANOM, XW_ANOM,                                 &
                     OMASK,                                            &
                     XLES_RESOLVED_W2 (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* pressure variance
!
  CALL LES_FLUX_ll ( ZP_ANOM, ZP_ANOM,                                 &
                     OMASK,                                            &
                     XLES_RESOLVED_P2 (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* potential temperature variance
!
  CALL LES_FLUX_ll ( ZTH_ANOM, ZTH_ANOM,                               &
                     OMASK,                                            &
                     XLES_RESOLVED_TH2(:,NLES_CURRENT_TCOUNT,IMASK)    )

!
!* resolved turbulent kinetic energy
!
  XLES_RESOLVED_Ke(:,NLES_CURRENT_TCOUNT,IMASK) = XUNDEF
!
  WHERE(XLES_RESOLVED_U2 (:,NLES_CURRENT_TCOUNT,IMASK) /= XUNDEF) &
  XLES_RESOLVED_Ke(:,NLES_CURRENT_TCOUNT,IMASK) = 0.5*(           &
                   XLES_RESOLVED_U2 (:,NLES_CURRENT_TCOUNT,IMASK) &
                 + XLES_RESOLVED_V2 (:,NLES_CURRENT_TCOUNT,IMASK) &
                 + XLES_RESOLVED_W2 (:,NLES_CURRENT_TCOUNT,IMASK))
!
!* potential temperature - virtual potential temperature covariance
!
  IF (LUSERV) THEN
    CALL LES_FLUX_ll ( ZTH_ANOM, ZTHV_ANOM,                              &
                       OMASK,                                            &
                       XLES_RESOLVED_THTHV(:,NLES_CURRENT_TCOUNT,IMASK)  )

!
!* vapor mixing ratio variance
!
    CALL LES_FLUX_ll ( ZRV_ANOM, ZRV_ANOM,                               &
                       OMASK,                                            &
                       XLES_RESOLVED_Rv2(:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!
!* potential temperature - vapor mixing ratio correlation
!
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRV_ANOM,                               &
                       OMASK,                                            &
                       XLES_RESOLVED_ThRv(:,NLES_CURRENT_TCOUNT,IMASK)   )
!
!* virtual potential temperature - vapor mixing ratio correlation
!
    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRV_ANOM,                              &
                       OMASK,                                            &
                       XLES_RESOLVED_ThvRv(:,NLES_CURRENT_TCOUNT,IMASK)  )
  END IF
!
!
!* liquid potential temperature - virtual potential temperature covariance
!
  IF (LUSERC) THEN
    CALL LES_FLUX_ll ( XTHL_ANOM, ZTHV_ANOM,                             &
                       OMASK,                                            &
                       XLES_RESOLVED_THLTHV(:,NLES_CURRENT_TCOUNT,IMASK) )
!
!* liquid potential temperature variance
!
    CALL LES_FLUX_ll ( XTHL_ANOM, XTHL_ANOM,                             &
                       OMASK,                                            &
                       XLES_RESOLVED_THL2(:,NLES_CURRENT_TCOUNT,IMASK)   )
!
!* total water mixing ratio variance
!
    CALL LES_FLUX_ll ( XRT_ANOM, XRT_ANOM,                               &
                       OMASK,                                            &
                       XLES_RESOLVED_Rt2(:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* cloud mixing ratio variance
!
    CALL LES_FLUX_ll ( ZRC_ANOM, ZRC_ANOM,                               &
                       OMASK,                                            &
                       XLES_RESOLVED_Rc2(:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* potential temperature - cloud mixing ratio correlation
!
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRC_ANOM,                               &
                       OMASK,                                            &
                       XLES_RESOLVED_ThRc(:,NLES_CURRENT_TCOUNT,IMASK)   )
!
!* liquid potential temperature - vapor mixing ratio correlation
!
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRV_ANOM,                               &
                       OMASK,                                             &
                       XLES_RESOLVED_ThlRv(:,NLES_CURRENT_TCOUNT,IMASK)   )
!
!* liquid potential temperature - cloud mixing ratio correlation
!
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRC_ANOM,                               &
                       OMASK,                                             &
                       XLES_RESOLVED_ThlRc(:,NLES_CURRENT_TCOUNT,IMASK)   )
!
!* virtual potential temperature - cloud mixing ratio correlation
!
    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRC_ANOM,                              &
                       OMASK,                                            &
                       XLES_RESOLVED_ThvRc(:,NLES_CURRENT_TCOUNT,IMASK)  )
!
! variance of lwp 
!
    IF (IMASK .EQ. 1) THEN
    CALL LES_FLUX_ll (ZLWP_ANOM, ZLWP_ANOM,                           &
                      OMASK(:,:,1),                                   &
                       XLES_LWPVAR(NLES_CURRENT_TCOUNT)  ) 
    END IF
  END IF
!
!* ice mixing ratio variance
!
  IF (LUSERI) THEN
    CALL LES_FLUX_ll ( ZRI_ANOM, ZRI_ANOM,                               &
                       OMASK,                                            &
                       XLES_RESOLVED_Ri2(:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* potential temperature - ice mixing ratio correlation
!
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRI_ANOM,                               &
                       OMASK,                                            &
                       XLES_RESOLVED_ThRi(:,NLES_CURRENT_TCOUNT,IMASK)   )
!
!* liquid potential temperature - ice mixing ratio correlation
!
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRI_ANOM,                               &
                       OMASK,                                             &
                       XLES_RESOLVED_ThlRi(:,NLES_CURRENT_TCOUNT,IMASK)   )
!
!* virtual potential temperature - ice mixing ratio correlation
!
    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRI_ANOM,                              &
                       OMASK,                                            &
                       XLES_RESOLVED_ThvRi(:,NLES_CURRENT_TCOUNT,IMASK)  )
  END IF
!
!* scalar variable mixing ratio variances
!
  DO JSV=1,NSV
    CALL LES_FLUX_ll ( XSV_ANOM(:,:,:,JSV), XSV_ANOM(:,:,:,JSV),             &
                       OMASK,                                                &
                       XLES_RESOLVED_Sv2(:,NLES_CURRENT_TCOUNT,IMASK,JSV)    )
!
!* potential temperature - scalar variables ratio correlation
!
    CALL LES_FLUX_ll ( ZTH_ANOM, XSV_ANOM(:,:,:,JSV),                      &
                       OMASK,                                              &
                       XLES_RESOLVED_ThSv(:,NLES_CURRENT_TCOUNT,IMASK,JSV) )
!
!* liquid potential temperature - scalar variables ratio correlation
!
    IF (LUSERC) THEN
      CALL LES_FLUX_ll ( XTHL_ANOM, XSV_ANOM(:,:,:,JSV),                     &
                         OMASK,                                              &
                         XLES_RESOLVED_ThlSv(:,NLES_CURRENT_TCOUNT,IMASK,JSV) )
    END IF
!
!* virtual potential temperature - scalar variables ratio correlation
!
    IF (LUSERV) THEN
      CALL LES_FLUX_ll ( ZTHV_ANOM, XSV_ANOM(:,:,:,JSV),                      &
                         OMASK,                                               &
                         XLES_RESOLVED_ThvSv(:,NLES_CURRENT_TCOUNT,IMASK,JSV) )
    END IF
  END DO
!
!
!* wind fluxes
!
  CALL LES_FLUX_ll ( XU_ANOM, XV_ANOM,                                 &
                     OMASK,                                            &
                     XLES_RESOLVED_UV (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
  CALL LES_FLUX_ll ( XW_ANOM, XU_ANOM,                                 &
                     OMASK,                                            &
                     XLES_RESOLVED_WU (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
  CALL LES_FLUX_ll ( XW_ANOM, XV_ANOM,                                 &
                     OMASK,                                            &
                     XLES_RESOLVED_WV (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* pressure fluxes
!
  CALL LES_FLUX_ll ( XU_ANOM, ZDPDZ_ANOM,                               &
                     OMASK,                                             &
                     XLES_RESOLVED_UP  (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
  CALL LES_FLUX_ll ( XV_ANOM, ZDPDZ_ANOM,                               &
                     OMASK,                                             &
                     XLES_RESOLVED_VP  (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
  CALL LES_FLUX_ll ( XW_ANOM, ZDPDZ_ANOM,                               &
                     OMASK,                                             &
                     XLES_RESOLVED_WP  (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* theta fluxes
!
  CALL LES_FLUX_ll ( XU_ANOM, ZTH_ANOM,                                 &
                     OMASK,                                             &
                     XLES_RESOLVED_UTh (:,NLES_CURRENT_TCOUNT,IMASK)    )

  CALL LES_FLUX_ll ( XV_ANOM, ZTH_ANOM,                                 &
                     OMASK,                                             &
                     XLES_RESOLVED_VTh (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
  CALL LES_FLUX_ll ( XW_ANOM, ZTH_ANOM,                                 &
                     OMASK,                                             &
                     XLES_RESOLVED_WTh (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* virtual theta fluxes
!
  IF (LUSERV) THEN
    CALL LES_FLUX_ll ( XU_ANOM, ZTHV_ANOM,                                 &
                       OMASK,                                              &
                       XLES_RESOLVED_UThv (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
    CALL LES_FLUX_ll ( XV_ANOM, ZTHV_ANOM,                                 &
                       OMASK,                                              &
                       XLES_RESOLVED_VThv (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
    CALL LES_FLUX_ll ( XW_ANOM, ZTHV_ANOM,                                 &
                       OMASK,                                              &
                       XLES_RESOLVED_WThv (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* vapor mixing ratio fluxes
!
    CALL LES_FLUX_ll ( XU_ANOM, ZRV_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_URv (:,NLES_CURRENT_TCOUNT,IMASK)     )
!
    CALL LES_FLUX_ll ( XV_ANOM, ZRV_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_VRv (:,NLES_CURRENT_TCOUNT,IMASK)     )
!
    CALL LES_FLUX_ll ( XW_ANOM, ZRV_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_WRv (:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!* cloud water mixing ratio fluxes
!
  IF (LUSERC) THEN
    CALL LES_FLUX_ll ( XU_ANOM, ZRC_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_URc (:,NLES_CURRENT_TCOUNT,IMASK)     )
!
    CALL LES_FLUX_ll ( XV_ANOM, ZRC_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_VRc (:,NLES_CURRENT_TCOUNT,IMASK)     )
!
    CALL LES_FLUX_ll ( XW_ANOM, ZRC_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_WRc (:,NLES_CURRENT_TCOUNT,IMASK)     )
!
!* liquid theta fluxes
!
    CALL LES_FLUX_ll ( XU_ANOM, XTHL_ANOM,                                 &
                       OMASK,                                              &
                       XLES_RESOLVED_UThl (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
    CALL LES_FLUX_ll ( XV_ANOM, XTHL_ANOM,                                 &
                       OMASK,                                              &
                       XLES_RESOLVED_VThl (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
    CALL LES_FLUX_ll ( XW_ANOM, XTHL_ANOM,                                 &
                       OMASK,                                              &
                       XLES_RESOLVED_WThl (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* total water mixing ratio fluxes
!
    CALL LES_FLUX_ll ( XW_ANOM, XRT_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_WRt (:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!
!* cloud ice mixing ratio fluxes
!
  IF (LUSERI) THEN
    CALL LES_FLUX_ll ( XU_ANOM, ZRI_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_URi (:,NLES_CURRENT_TCOUNT,IMASK)     )
!
  CALL LES_FLUX_ll ( XV_ANOM, ZRI_ANOM,                                  &
                     OMASK,                                              &
                     XLES_RESOLVED_VRi (:,NLES_CURRENT_TCOUNT,IMASK)     )
!
    CALL LES_FLUX_ll ( XW_ANOM, ZRI_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_WRi (:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
  IF (LUSERR) THEN
     CALL LES_FLUX_ll ( XW_ANOM, ZRR_ANOM,                                  &
                       OMASK,                                              &
                       XLES_RESOLVED_WRr (:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF
!

!
!* scalar variables fluxes
!
  DO JSV=1,NSV
    CALL LES_FLUX_ll ( XU_ANOM, XSV_ANOM(:,:,:,JSV),                           &
                       OMASK,                                                  &
                       XLES_RESOLVED_USv (:,NLES_CURRENT_TCOUNT,IMASK,JSV)     )
!
   CALL LES_FLUX_ll ( XV_ANOM, XSV_ANOM(:,:,:,JSV),                            &
                       OMASK,                                                  &
                       XLES_RESOLVED_VSv (:,NLES_CURRENT_TCOUNT,IMASK,JSV)     )
!
    CALL LES_FLUX_ll ( XW_ANOM, XSV_ANOM(:,:,:,JSV),                           &
                       OMASK,                                                  &
                       XLES_RESOLVED_WSv (:,NLES_CURRENT_TCOUNT,IMASK,JSV)     )
  END DO
!
!* skewness
!
  CALL LES_3RD_MOMENT_ll ( XU_ANOM, XU_ANOM, XU_ANOM,                       &
                          OMASK,                                            &
                          XLES_RESOLVED_U3 (:,NLES_CURRENT_TCOUNT,IMASK)    )
 
  CALL LES_3RD_MOMENT_ll ( XV_ANOM, XV_ANOM, XV_ANOM,                       &
                          OMASK,                                            &
                          XLES_RESOLVED_V3 (:,NLES_CURRENT_TCOUNT,IMASK)    )
 
  CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, XW_ANOM,                       &
                          OMASK,                                            &
                          XLES_RESOLVED_W3 (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* kurtosis
!
  CALL LES_4TH_MOMENT_ll ( XU_ANOM, XU_ANOM, XU_ANOM,  XU_ANOM,              &
                          OMASK,                                            &
                          XLES_RESOLVED_U4 (:,NLES_CURRENT_TCOUNT,IMASK)    )
 
  CALL LES_4TH_MOMENT_ll ( XV_ANOM, XV_ANOM, XV_ANOM,  XV_ANOM,             &
                          OMASK,                                            &
                          XLES_RESOLVED_V4 (:,NLES_CURRENT_TCOUNT,IMASK)    )
 
  CALL LES_4TH_MOMENT_ll ( XW_ANOM, XW_ANOM, XW_ANOM,   XW_ANOM,            &
                          OMASK,                                            &
                          XLES_RESOLVED_W4 (:,NLES_CURRENT_TCOUNT,IMASK)    )
!
!* third moments of liquid potential temperature
!
  IF (LUSERC) THEN
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XTHL_ANOM, XTHL_ANOM,                    &
                            OMASK,                                             &
                            XLES_RESOLVED_WThl2(:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, XTHL_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_W2Thl(:,NLES_CURRENT_TCOUNT,IMASK)   )

  ELSE
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZTH_ANOM, ZTH_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WThl2(:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, ZTH_ANOM,                       &
                            OMASK,                                             &
                            XLES_RESOLVED_W2Thl(:,NLES_CURRENT_TCOUNT,IMASK)   )
  END IF
!
!* third moments of water vapor
!
  IF (LUSERV) THEN
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZRV_ANOM, ZRV_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WRv2 (:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, ZRV_ANOM,                       &
                            OMASK,                                             &
                            XLES_RESOLVED_W2Rv (:,NLES_CURRENT_TCOUNT,IMASK)   )
  END IF

  IF (LUSERC) THEN
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XTHL_ANOM, ZRV_ANOM,                     &
                            OMASK,                                             &
                            XLES_RESOLVED_WThlRv(:,NLES_CURRENT_TCOUNT,IMASK)  )
  ELSE IF (LUSERV) THEN
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZTH_ANOM, ZRV_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WThlRv(:,NLES_CURRENT_TCOUNT,IMASK)  )
  END IF
!
!* third moments of total water
! 
  IF (LUSERC) THEN
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XRT_ANOM, XRT_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WRt2 (:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, XRT_ANOM,                       &
                            OMASK,                                             &
                            XLES_RESOLVED_W2Rt (:,NLES_CURRENT_TCOUNT,IMASK)   )
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XTHL_ANOM, XRT_ANOM,                       &
                            OMASK,                                             &
                            XLES_RESOLVED_WThlRt (:,NLES_CURRENT_TCOUNT,IMASK)   )                        
  ELSE IF (LUSERV) THEN
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZRV_ANOM, ZRV_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WRt2 (:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, ZRV_ANOM,                       &
                            OMASK,                                             &
                            XLES_RESOLVED_W2Rt (:,NLES_CURRENT_TCOUNT,IMASK)   )
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZTH_ANOM, ZRV_ANOM,                       &
                            OMASK,                                             &
                            XLES_RESOLVED_WThlRt (:,NLES_CURRENT_TCOUNT,IMASK)   )                        
  END IF                            
!
!* third moments of cloud water
!
  IF (LUSERC) THEN
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZRC_ANOM, ZRC_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WRc2 (:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, ZRC_ANOM,                       &
                            OMASK,                                             &
                            XLES_RESOLVED_W2Rc (:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XTHL_ANOM, ZRC_ANOM,                     &
                            OMASK,                                             &
                            XLES_RESOLVED_WThlRc(:,NLES_CURRENT_TCOUNT,IMASK)  )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZRV_ANOM, ZRC_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WRvRc (:,NLES_CURRENT_TCOUNT,IMASK)  )
  END IF
!
!* third moments of cloud ice
!
  IF (LUSERI) THEN
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZRI_ANOM, ZRI_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WRi2 (:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, ZRI_ANOM,                      &
                           OMASK,                                             &
                           XLES_RESOLVED_W2Ri (:,NLES_CURRENT_TCOUNT,IMASK)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XTHL_ANOM, ZRI_ANOM,                     &
                            OMASK,                                             &
                            XLES_RESOLVED_WThlRi(:,NLES_CURRENT_TCOUNT,IMASK)  )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZRV_ANOM, ZRI_ANOM,                      &
                            OMASK,                                             &
                            XLES_RESOLVED_WRvRi (:,NLES_CURRENT_TCOUNT,IMASK)  )
  END IF
!
!* third moments of scalar variables
!
  DO JSV=1,NSV
    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XSV_ANOM(:,:,:,JSV), XSV_ANOM(:,:,:,JSV),    &
                            OMASK,                                                 &
                            XLES_RESOLVED_WSv2 (:,NLES_CURRENT_TCOUNT,IMASK,JSV)   )

    CALL LES_3RD_MOMENT_ll ( XW_ANOM, XW_ANOM, XSV_ANOM(:,:,:,JSV),                &
                            OMASK,                                                 &
                            XLES_RESOLVED_W2Sv (:,NLES_CURRENT_TCOUNT,IMASK,JSV)   )
    IF (LUSERC) THEN
      CALL LES_3RD_MOMENT_ll ( XW_ANOM, XTHL_ANOM, XSV_ANOM(:,:,:,JSV),              &
                              OMASK,                                                 &
                              XLES_RESOLVED_WThlSv(:,NLES_CURRENT_TCOUNT,IMASK,JSV)  )
    ELSE
      CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZTH_ANOM, XSV_ANOM(:,:,:,JSV),               &
                              OMASK,                                                 &
                              XLES_RESOLVED_WThlSv(:,NLES_CURRENT_TCOUNT,IMASK,JSV)  )
    END IF

    IF (LUSERV) THEN
      CALL LES_3RD_MOMENT_ll ( XW_ANOM, ZRV_ANOM, XSV_ANOM(:,:,:,JSV),               &
                              OMASK,                                                 &
                              XLES_RESOLVED_WRvSv (:,NLES_CURRENT_TCOUNT,IMASK,JSV)  )
    END IF
  END DO
!
!* presso-correlations
!
!
  CALL LES_FLUX_ll ( XTHL_ANOM, ZDPDZ_ANOM,                              &
                     OMASK,                                              &
                     XLES_RESOLVED_ThlPz(:,NLES_CURRENT_TCOUNT,IMASK)    )

  IF (LUSERV) &
  CALL LES_FLUX_ll ( ZRV_ANOM, ZDPDZ_ANOM,                               &
                     OMASK,                                              &
                     XLES_RESOLVED_RvPz(:,NLES_CURRENT_TCOUNT,IMASK)     )

  IF (LUSERC) THEN
    CALL LES_FLUX_ll ( XRT_ANOM, ZDPDZ_ANOM,                               &
                       OMASK,                                              &
                       XLES_RESOLVED_RtPz(:,NLES_CURRENT_TCOUNT,IMASK)     )

    CALL LES_FLUX_ll ( ZRC_ANOM, ZDPDZ_ANOM,                               &
                       OMASK,                                              &
                       XLES_RESOLVED_RcPz(:,NLES_CURRENT_TCOUNT,IMASK)     )
  END IF

  IF (LUSERI) &
  CALL LES_FLUX_ll ( ZRI_ANOM, ZDPDZ_ANOM,                               &
                     OMASK,                                              &
                     XLES_RESOLVED_RiPz(:,NLES_CURRENT_TCOUNT,IMASK)     )

!
!
!* resolved turbulent kinetic energy fluxes
!

  CALL LES_3RD_MOMENT_ll ( XU_ANOM, XU_ANOM, XU_ANOM,                        &
                          OMASK,                                             &
                          ZLES_RESOLVED_U3  (:)                              )

  CALL LES_3RD_MOMENT_ll ( XU_ANOM, XV_ANOM, XV_ANOM,                        &
                          OMASK,                                             &
                          ZLES_RESOLVED_UV2 (:)                              )

  CALL LES_3RD_MOMENT_ll ( XU_ANOM, XW_ANOM, XW_ANOM,                        &
                          OMASK,                                             &
                          ZLES_RESOLVED_UW2 (:)                              )

  XLES_RESOLVED_UKe(:,NLES_CURRENT_TCOUNT,IMASK) = 0.5*(  ZLES_RESOLVED_U3  &
                                                        + ZLES_RESOLVED_UV2 &
                                                        + ZLES_RESOLVED_UW2 )



  CALL LES_3RD_MOMENT_ll ( XV_ANOM, XU_ANOM, XU_ANOM,                        &
                          OMASK,                                             &
                          ZLES_RESOLVED_VU2 (:)                              )

  CALL LES_3RD_MOMENT_ll ( XV_ANOM, XV_ANOM, XV_ANOM,                        &
                          OMASK,                                             &
                          ZLES_RESOLVED_V3  (:)                              )

  CALL LES_3RD_MOMENT_ll ( XV_ANOM, XW_ANOM, XW_ANOM,                        &
                          OMASK,                                             &
                          ZLES_RESOLVED_VW2 (:)                              )

  XLES_RESOLVED_VKe(:,NLES_CURRENT_TCOUNT,IMASK) = 0.5*(  ZLES_RESOLVED_VU2 &
                                                        + ZLES_RESOLVED_V3  &
                                                        + ZLES_RESOLVED_VW2 )


  CALL LES_3RD_MOMENT_ll ( XW_ANOM, XU_ANOM, XU_ANOM,                        &
                          OMASK,                                             &
                          ZLES_RESOLVED_WU2 (:)                              )

  CALL LES_3RD_MOMENT_ll ( XW_ANOM, XV_ANOM, XV_ANOM,                        &
                          OMASK,                                             &
                          ZLES_RESOLVED_WV2 (:)                              )

  XLES_RESOLVED_WKe(:,NLES_CURRENT_TCOUNT,IMASK) = 0.5*(  ZLES_RESOLVED_WU2  &
                                                        + ZLES_RESOLVED_WV2  &
                           + XLES_RESOLVED_W3(:,NLES_CURRENT_TCOUNT,IMASK)   )

!
!
!-------------------------------------------------------------------------------
!
!       1.6  Subgrid quantities
!            ------------------
!
IF (LLES_SUBGRID) THEN
!
!* wind fluxes and variances
!
  CALL LES_MEAN_ll ( ZTKE_LES, OMASK,                                  &
                     XLES_SUBGRID_Tke(:,NLES_CURRENT_TCOUNT,IMASK)     )

  CALL LES_MEAN_MPROC ( XLES_SUBGRID_UV(:,NLES_CURRENT_TCOUNT,IMASK),   &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WU(:,NLES_CURRENT_TCOUNT,IMASK),   &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WV(:,NLES_CURRENT_TCOUNT,IMASK),   &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_U2(:,NLES_CURRENT_TCOUNT,IMASK),   &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_V2(:,NLES_CURRENT_TCOUNT,IMASK),   &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_W2(:,NLES_CURRENT_TCOUNT,IMASK),   &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
!
!* liquid potential temperature fluxes
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_UThl(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_VThl(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WThl(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!

!* liquid potential temperature variance
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_Thl2(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )
!
!* Mass flux scheme of shallow convection
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_THLUP_MF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )  
                        
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_RTUP_MF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )   
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_RVUP_MF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )   
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_RCUP_MF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )   
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_RIUP_MF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )   
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WUP_MF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )    
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_MASSFLUX(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )  
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_DETR(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )      
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_ENTR(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )      
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_FRACUP(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )     
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_THVUP_MF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )  
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WTHLMF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )    
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WRTMF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )  
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WTHVMF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        ) 
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WUMF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )   
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WVMF(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                        )   

                        
!* total water mixing ratio fluxes, correlation and variance
!
  IF (LUSERV) THEN
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_URt(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                        )
    !
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_VRt(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                        )
    !
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_WRt(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                        )
    !
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_ThlRt(:,NLES_CURRENT_TCOUNT,IMASK),&
                          IAVG_PTS(:), IUND_PTS(:)                        )
    !
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_Rt2(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                        )
  END IF
!
!* scalar variances
!
  DO JSV=1,NSV
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_Sv2(:,NLES_CURRENT_TCOUNT,IMASK,JSV), &
                          IAVG_PTS(:), IUND_PTS(:)                           )
  END DO
!
!* cloud water mixing ratio fluxes
!
  IF (LUSERC) THEN
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_URc(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                        )
    !
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_VRc(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                        )
    !
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_WRc(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                        )
  END IF
!
!* scalar variables fluxes
!
  DO JSV=1,NSV
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_USv(:,NLES_CURRENT_TCOUNT,IMASK,JSV),  &
                          IAVG_PTS(:), IUND_PTS(:)                            )
    !
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_VSv(:,NLES_CURRENT_TCOUNT,IMASK,JSV),  &
                          IAVG_PTS(:), IUND_PTS(:)                            )
    !
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_WSv(:,NLES_CURRENT_TCOUNT,IMASK,JSV),  &
                          IAVG_PTS(:), IUND_PTS(:)                            )
  END DO
!
!* subgrid turbulent kinetic energy fluxes
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_UTke(:,NLES_CURRENT_TCOUNT,IMASK),  &
                        IAVG_PTS(:), IUND_PTS(:)                         )
  !
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_VTke(:,NLES_CURRENT_TCOUNT,IMASK),  &
                        IAVG_PTS(:), IUND_PTS(:)                         )
  !
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WTke(:,NLES_CURRENT_TCOUNT,IMASK),  &
                        IAVG_PTS(:), IUND_PTS(:)                         )

  CALL LES_MEAN_MPROC ( XLES_SUBGRID_ddz_WTke(:,NLES_CURRENT_TCOUNT,IMASK),&
                        IAVG_PTS(:), IUND_PTS(:)                           )
!
!* fluxes and correlations with virtual potential temperature
!
  IF (LUSERV) THEN
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_WThv(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                         )

    CALL LES_MEAN_MPROC ( XLES_SUBGRID_ThlThv(:,NLES_CURRENT_TCOUNT,IMASK),&
                          IAVG_PTS(:), IUND_PTS(:)                         )

    CALL LES_MEAN_MPROC ( XLES_SUBGRID_RtThv(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                         )

    DO JSV=1,NSV
      CALL LES_MEAN_MPROC ( XLES_SUBGRID_SvThv(:,NLES_CURRENT_TCOUNT,IMASK,JSV),  &
                            IAVG_PTS(:), IUND_PTS(:)                              )
   END DO
  END IF
!
!* third order fluxes
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_W2Thl(:,NLES_CURRENT_TCOUNT,IMASK),   &
                          IAVG_PTS(:), IUND_PTS(:)                         )

  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WThl2(:,NLES_CURRENT_TCOUNT,IMASK),   &
                          IAVG_PTS(:), IUND_PTS(:)                         )

  IF (LUSERV) THEN
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_W2Rt(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                         )

    CALL LES_MEAN_MPROC ( XLES_SUBGRID_WThlRt(:,NLES_CURRENT_TCOUNT,IMASK),&
                          IAVG_PTS(:), IUND_PTS(:)                         )

    CALL LES_MEAN_MPROC ( XLES_SUBGRID_WRt2(:,NLES_CURRENT_TCOUNT,IMASK),  &
                          IAVG_PTS(:), IUND_PTS(:)                         )

  END IF
    DO JSV=1,NSV
      CALL LES_MEAN_MPROC ( XLES_SUBGRID_W2Sv(:,NLES_CURRENT_TCOUNT,IMASK,JSV),  &
                            IAVG_PTS(:), IUND_PTS(:)                             )

      CALL LES_MEAN_MPROC ( XLES_SUBGRID_WSv2(:,NLES_CURRENT_TCOUNT,IMASK,JSV),  &
                            IAVG_PTS(:), IUND_PTS(:)                             )
   END DO
!
!* dissipative terms
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_DISS_Tke(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                          )

  CALL LES_MEAN_MPROC ( XLES_SUBGRID_DISS_Thl2(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                           )

  IF (LUSERV) THEN
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_DISS_Rt2(:,NLES_CURRENT_TCOUNT,IMASK),&
                          IAVG_PTS(:), IUND_PTS(:)                           )

    CALL LES_MEAN_MPROC ( XLES_SUBGRID_DISS_ThlRt(:,NLES_CURRENT_TCOUNT,IMASK),&
                          IAVG_PTS(:), IUND_PTS(:)                             )
  END IF

  DO JSV=1,NSV
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_DISS_Sv2(:,NLES_CURRENT_TCOUNT,IMASK,JSV), &
                          IAVG_PTS(:), IUND_PTS(:)                                )
  END DO
!
!* presso-correlation terms
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_WP(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                    )

  CALL LES_MEAN_MPROC ( XLES_SUBGRID_ThlPz(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                       )

  IF (LUSERV) THEN
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_RtPz(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                        )

  END IF

  DO JSV=1,NSV
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_SvPz(:,NLES_CURRENT_TCOUNT,IMASK,JSV), &
                          IAVG_PTS(:), IUND_PTS(:)                            )
  END DO

!* phi3 and psi3 terms
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_PHI3(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                    )

  IF (LUSERV) THEN
    CALL LES_MEAN_MPROC ( XLES_SUBGRID_PSI3(:,NLES_CURRENT_TCOUNT,IMASK), &
                            IAVG_PTS(:), IUND_PTS(:)                    )
  END IF 
!
!* subgrid mixing length
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_LMix(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                      )
!
!* subgrid dissipative length
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_LDiss(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                       )
!
!* eddy diffusivities
!
  CALL LES_MEAN_MPROC ( XLES_SUBGRID_Km(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                    )

  CALL LES_MEAN_MPROC ( XLES_SUBGRID_Kh(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                    )
  
END IF
!
! computation of KHT and KHR depending on LLES
  IF (LUSERC) THEN 
   IF (LLES_RESOLVED) THEN 
    XLES_MEAN_KHt(:,NLES_CURRENT_TCOUNT,IMASK)=0.       
    WHERE(XLES_MEAN_DTHLDZ(:,NLES_CURRENT_TCOUNT,IMASK)/=0)              &
          XLES_MEAN_KHt(:,NLES_CURRENT_TCOUNT,IMASK)=-1.                 &
                      *XLES_RESOLVED_WThl (:,NLES_CURRENT_TCOUNT,IMASK)/ &
                       XLES_MEAN_DTHLDZ(:,NLES_CURRENT_TCOUNT,IMASK)
    XLES_MEAN_KHr(:,NLES_CURRENT_TCOUNT,IMASK)=0.                   
    WHERE(XLES_MEAN_DRtDZ(:,NLES_CURRENT_TCOUNT,IMASK)/=0)               &
          XLES_MEAN_KHr(:,NLES_CURRENT_TCOUNT,IMASK)=-1.*                &
                        XLES_RESOLVED_WRt (:,NLES_CURRENT_TCOUNT,IMASK)/ &
                       XLES_MEAN_DRtDZ(:,NLES_CURRENT_TCOUNT,IMASK)
   END IF                    
   IF (LLES_SUBGRID) THEN 
    XLES_MEAN_KHt(:,NLES_CURRENT_TCOUNT,IMASK)=0.       
    WHERE(XLES_MEAN_DTHLDZ(:,NLES_CURRENT_TCOUNT,IMASK)/=0)              & 
          XLES_MEAN_KHt(:,NLES_CURRENT_TCOUNT,IMASK)=-1.                 &
                      *XLES_SUBGRID_WThl (:,NLES_CURRENT_TCOUNT,IMASK) /  &
                       XLES_MEAN_DTHLDZ(:,NLES_CURRENT_TCOUNT,IMASK)
    XLES_MEAN_KHr(:,NLES_CURRENT_TCOUNT,IMASK)=0.                   
    WHERE(XLES_MEAN_DRtDZ(:,NLES_CURRENT_TCOUNT,IMASK)/=0)               &
          XLES_MEAN_KHr(:,NLES_CURRENT_TCOUNT,IMASK)=-1.*                &
                     XLES_SUBGRID_WRt (:,NLES_CURRENT_TCOUNT,IMASK) / &
                       XLES_MEAN_DRtDZ(:,NLES_CURRENT_TCOUNT,IMASK)
   END IF
   IF (LLES_RESOLVED .AND. LLES_SUBGRID) THEN 
    XLES_MEAN_KHt(:,NLES_CURRENT_TCOUNT,IMASK)=0.       
    WHERE(XLES_MEAN_DTHLDZ(:,NLES_CURRENT_TCOUNT,IMASK)/=0)              &
          XLES_MEAN_KHt(:,NLES_CURRENT_TCOUNT,IMASK)=-1.                 &
                      *(XLES_RESOLVED_WThl (:,NLES_CURRENT_TCOUNT,IMASK)+ &
                        XLES_SUBGRID_WThl (:,NLES_CURRENT_TCOUNT,IMASK))/ &
                       XLES_MEAN_DTHLDZ(:,NLES_CURRENT_TCOUNT,IMASK)
    XLES_MEAN_KHr(:,NLES_CURRENT_TCOUNT,IMASK)=0.                   
    WHERE(XLES_MEAN_DRtDZ(:,NLES_CURRENT_TCOUNT,IMASK)/=0)               &
          XLES_MEAN_KHr(:,NLES_CURRENT_TCOUNT,IMASK)=-1.*                &
                        (XLES_RESOLVED_WRt (:,NLES_CURRENT_TCOUNT,IMASK)+ &
                     XLES_SUBGRID_WRt (:,NLES_CURRENT_TCOUNT,IMASK)) / &
                       XLES_MEAN_DRtDZ(:,NLES_CURRENT_TCOUNT,IMASK)
   END IF
  END IF
!-------------------------------------------------------------------------------
!
!       1.7  Interaction of subgrid and resolved quantities
!            ----------------------------------------------
!
!* WARNING: these terms also contain the term due to the mean flow.
!           this mean flow contribution will be removed from them
!           when treated in write_les_budgetn.f90
!
!
!* subgrid turbulent kinetic energy fluxes
!
IF (LLES_RESOLVED) THEN
  CALL LES_FLUX_ll ( XU_ANOM, ZTKE_LES,                                 &
                     OMASK,                                             &
                     XLES_RES_U_SBG_Tke(:,NLES_CURRENT_TCOUNT,IMASK)    )
!
  CALL LES_FLUX_ll ( XV_ANOM, ZTKE_LES,                                 &
                     OMASK,                                             &
                     XLES_RES_V_SBG_Tke(:,NLES_CURRENT_TCOUNT,IMASK)    )
!
  CALL LES_FLUX_ll ( XW_ANOM, ZTKE_LES,                                 &
                     OMASK,                                             &
                     XLES_RES_W_SBG_Tke(:,NLES_CURRENT_TCOUNT,IMASK)    )
END IF
!
!* WARNING: these terms also contain the term due to the mean flow.
!           this mean flow contribution will be removed from them
!           when treated in write_les_budgetn.f90
!
!* production terms for subgrid quantities
!
IF (LLES_RESOLVED .AND. LLES_SUBGRID) THEN
  CALL LES_MEAN_MPROC ( XLES_RES_ddxa_U_SBG_UaU(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                              )

  CALL LES_MEAN_MPROC ( XLES_RES_ddxa_V_SBG_UaV(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                              )

  CALL LES_MEAN_MPROC ( XLES_RES_ddxa_W_SBG_UaW(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                              )

  CALL LES_MEAN_MPROC ( XLES_RES_ddxa_W_SBG_UaThl(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                                )

  CALL LES_MEAN_MPROC ( XLES_RES_ddxa_Thl_SBG_UaW(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                                )

  CALL LES_MEAN_MPROC ( XLES_RES_ddz_Thl_SBG_W2  (:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                                )

  CALL LES_MEAN_MPROC ( XLES_RES_ddxa_Thl_SBG_UaThl(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                                  )

  IF (LUSERV) THEN
    CALL LES_MEAN_MPROC ( XLES_RES_ddxa_W_SBG_UaRt(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                               )

    CALL LES_MEAN_MPROC ( XLES_RES_ddxa_Rt_SBG_UaW(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                               )

    CALL LES_MEAN_MPROC ( XLES_RES_ddz_Rt_SBG_W2  (:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                               )

    CALL LES_MEAN_MPROC ( XLES_RES_ddxa_Thl_SBG_UaRt(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                                 )

    CALL LES_MEAN_MPROC ( XLES_RES_ddxa_Rt_SBG_UaThl(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                                 )

    CALL LES_MEAN_MPROC ( XLES_RES_ddxa_Rt_SBG_UaRt(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                                )

  END IF
!
!* WARNING: these terms also contain the term due to the mean flow.
!           this mean flow contribution will be removed from them
!           when treated in write_les_budgetn.f90
!
!* turbulent transport and advection terms for subgrid quantities
!
  CALL LES_MEAN_MPROC ( XLES_RES_W_SBG_WThl(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                          )

  CALL LES_MEAN_MPROC ( XLES_RES_W_SBG_Thl2(:,NLES_CURRENT_TCOUNT,IMASK), &
                        IAVG_PTS(:), IUND_PTS(:)                          )

  IF (LUSERV) THEN
    CALL LES_MEAN_MPROC ( XLES_RES_W_SBG_WRt(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                         )

    CALL LES_MEAN_MPROC ( XLES_RES_W_SBG_Rt2(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                         )

    CALL LES_MEAN_MPROC ( XLES_RES_W_SBG_ThlRt(:,NLES_CURRENT_TCOUNT,IMASK), &
                          IAVG_PTS(:), IUND_PTS(:)                           )

  END IF

  DO JSV=1,NSV
    CALL LES_MEAN_MPROC ( XLES_RES_W_SBG_WSv(:,NLES_CURRENT_TCOUNT,IMASK,JSV), &
                          IAVG_PTS(:), IUND_PTS(:)                             )

    CALL LES_MEAN_MPROC ( XLES_RES_W_SBG_Sv2(:,NLES_CURRENT_TCOUNT,IMASK,JSV), &
                          IAVG_PTS(:), IUND_PTS(:)                             )
  END DO
END IF
!
!-------------------------------------------------------------------------------
!
!       2.   The following is for cartesian mask only
!            ----------------------------------------
!
IF (IMASK>1) RETURN
!
!-------------------------------------------------------------------------------
!
!       3.   Updraft diagnostics
!            -------------------
!
IF (LLES_UPDRAFT) THEN
!
  DO JK=1,NLES_K
    GUPDRAFT_MASK(:,:,JK) = (XW_ANOM(:,:,JK) > 0.) .AND. LLES_CURRENT_CART_MASK(:,:,JK)
  END DO
!
!
!       3.1  Updraft fraction
!            ----------------
!
  ZUPDRAFT(:,:,:) = 0.
  WHERE (GUPDRAFT_MASK(:,:,:))
    ZUPDRAFT(:,:,:) = 1.
  END WHERE
!
  CALL LES_MEAN_ll ( ZUPDRAFT, OMASK,                        &
                     XLES_UPDRAFT(:,NLES_CURRENT_TCOUNT)     )
!
!
!       3.2  Updraft mean quantities
!            -----------------------
!
!* vertical wind velocity
!
  CALL LES_MEAN_ll ( ZW_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_W(:,NLES_CURRENT_TCOUNT)     )
!
!* potential temperature
!
  CALL LES_MEAN_ll ( ZTH_LES, GUPDRAFT_MASK,                   &
                     XLES_UPDRAFT_Th(:,NLES_CURRENT_TCOUNT)    )
!
!* liquid potential temperature
!
  IF (LUSERC) &
  CALL LES_MEAN_ll ( ZTHL_LES, GUPDRAFT_MASK,                   &
                     XLES_UPDRAFT_Thl(:,NLES_CURRENT_TCOUNT)    )
!
!* virtual potential temperature
!
  IF (LUSERV) &
  CALL LES_MEAN_ll ( ZTHV_LES, GUPDRAFT_MASK,                   &
                     XLES_UPDRAFT_Thv(:,NLES_CURRENT_TCOUNT)    )
!
!* vapor mixing ratio
!
  IF (LUSERV) &
  CALL LES_MEAN_ll ( ZRV_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_Rv(:,NLES_CURRENT_TCOUNT)     )
!
!* cloud water mixing ratio
!
  IF (LUSERC) &
  CALL LES_MEAN_ll ( ZRC_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_Rc(:,NLES_CURRENT_TCOUNT)     )
!
!* rain mixing ratio
!
  IF (LUSERR) &
  CALL LES_MEAN_ll ( ZRR_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_Rr(:,NLES_CURRENT_TCOUNT)     )
!
!* cloud ice mixing ratio
!
  IF (LUSERI) &
  CALL LES_MEAN_ll ( ZRI_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_Ri(:,NLES_CURRENT_TCOUNT)     )
!
!* snow mixing ratio
!
  IF (LUSERS) &
  CALL LES_MEAN_ll ( ZRS_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_Rs(:,NLES_CURRENT_TCOUNT)     )
!
!* graupel mixing ratio
!
  IF (LUSERG) &
  CALL LES_MEAN_ll ( ZRG_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_Rg(:,NLES_CURRENT_TCOUNT)     )
!
!* hail mixing ratio
!
  IF (LUSERH) &
  CALL LES_MEAN_ll ( ZRG_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_Rh(:,NLES_CURRENT_TCOUNT)     )
!
!* scalar variables
!
  DO JSV=1,NSV
    CALL LES_MEAN_ll ( ZSV_LES(:,:,:,JSV), GUPDRAFT_MASK,             &
                       XLES_UPDRAFT_Sv(:,NLES_CURRENT_TCOUNT,JSV)     )
  END DO
!
!* subgrid turbulent kinetic energy
!
  CALL LES_MEAN_ll ( ZTKE_LES, GUPDRAFT_MASK,                    &
                     XLES_UPDRAFT_Tke(:,NLES_CURRENT_TCOUNT)     )
!
!
!       3.3  Updraft resolved quantities
!            ---------------------------
!
!
!* resolved turbulent kinetic energy
!
  CALL LES_FLUX_ll ( XU_ANOM, XU_ANOM,                           &
                     GUPDRAFT_MASK,                              &
                     ZLES_UPDRAFT_U2(:)                          )

  CALL LES_FLUX_ll ( XV_ANOM, XV_ANOM,                           &
                     GUPDRAFT_MASK,                              &
                     ZLES_UPDRAFT_V2(:)                          )

  CALL LES_FLUX_ll ( XW_ANOM, XW_ANOM,                           &
                     GUPDRAFT_MASK,                              &
                     ZLES_UPDRAFT_W2(:)                          )

  XLES_UPDRAFT_Ke(:,NLES_CURRENT_TCOUNT)  = 0.5 * (  ZLES_UPDRAFT_U2(:) &
                                                   + ZLES_UPDRAFT_V2(:) &
                                                   + ZLES_UPDRAFT_W2(:) )
!
!* vertical potential temperature flux
!
  CALL LES_FLUX_ll ( XW_ANOM, ZTH_ANOM,                          &
                     GUPDRAFT_MASK,                              &
                     XLES_UPDRAFT_WTh(:,NLES_CURRENT_TCOUNT)     )
!
!* vertical liquid potential temperature flux
!
  IF (LUSERC) &
  CALL LES_FLUX_ll ( XW_ANOM, XTHL_ANOM,                          &
                     GUPDRAFT_MASK,                               &
                     XLES_UPDRAFT_WThl(:,NLES_CURRENT_TCOUNT)     )
!
!* vertical virtual potential temperature flux
!
  IF (LUSERV) &
  CALL LES_FLUX_ll ( XW_ANOM, ZTHV_ANOM,                          &
                     GUPDRAFT_MASK,                               &
                     XLES_UPDRAFT_WThv(:,NLES_CURRENT_TCOUNT)     )
!
!* potential temperature variance
!
  CALL LES_FLUX_ll ( ZTH_ANOM, ZTH_ANOM,                          &
                     GUPDRAFT_MASK,                               &
                     XLES_UPDRAFT_Th2(:,NLES_CURRENT_TCOUNT)      )
!
!* liquid potential temperature variance
!
  IF (LUSERC) &
  CALL LES_FLUX_ll ( XTHL_ANOM, XTHL_ANOM,                         &
                     GUPDRAFT_MASK,                                &
                     XLES_UPDRAFT_Thl2(:,NLES_CURRENT_TCOUNT)      )
!
!* potential temperature - virtual potential temperature covariance
!
  IF (LUSERV) &
  CALL LES_FLUX_ll ( ZTH_ANOM, ZTHV_ANOM,                          &
                     GUPDRAFT_MASK,                                &
                     XLES_UPDRAFT_ThThv (:,NLES_CURRENT_TCOUNT)    )
!
!* liquid potential temperature - virtual potential temperature covariance
!
  IF (LUSERC) &
  CALL LES_FLUX_ll ( XTHL_ANOM, ZTHV_ANOM,                         &
                     GUPDRAFT_MASK,                                &
                     XLES_UPDRAFT_ThlThv(:,NLES_CURRENT_TCOUNT)    )
!
!* water vapor mixing ratio flux, variance and correlations
!
  IF (LUSERV) THEN
    CALL LES_FLUX_ll ( XW_ANOM, ZRV_ANOM,                           &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_WRv(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZRV_ANOM, ZRV_ANOM,                          &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_Rv2(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRV_ANOM,                          &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThRv (:,NLES_CURRENT_TCOUNT)    )
    !
    IF (LUSERC) &
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRV_ANOM,                         &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThlRv(:,NLES_CURRENT_TCOUNT)    )

    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRV_ANOM,                         &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThvRv(:,NLES_CURRENT_TCOUNT)    )
  END IF
!
!* cloud water mixing ratio flux
!
  IF (LUSERC) THEN
    CALL LES_FLUX_ll ( XW_ANOM, ZRC_ANOM,                           &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_WRc(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZRC_ANOM, ZRC_ANOM,                          &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_Rc2(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRC_ANOM,                          &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThRc (:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRC_ANOM,                         &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThlRc(:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRC_ANOM,                         &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThvRc(:,NLES_CURRENT_TCOUNT)    )
  END IF
!
!* cloud ice mixing ratio flux
!
  IF (LUSERI) THEN
    CALL LES_FLUX_ll ( XW_ANOM, ZRI_ANOM,                           &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_WRi(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZRI_ANOM, ZRI_ANOM,                          &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_Ri2(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRI_ANOM,                          &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThRi (:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRI_ANOM,                         &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThlRi(:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRI_ANOM,                         &
                       GUPDRAFT_MASK,                               &
                       XLES_UPDRAFT_ThvRi(:,NLES_CURRENT_TCOUNT)    )
  END IF
!
!* scalar variables flux
!
  DO JSV=1,NSV
    CALL LES_FLUX_ll ( XW_ANOM, XSV_ANOM(:,:,:,JSV),                    &
                       GUPDRAFT_MASK,                                   &
                       XLES_UPDRAFT_WSv(:,NLES_CURRENT_TCOUNT,JSV)      )
    !
    CALL LES_FLUX_ll ( XSV_ANOM(:,:,:,JSV), XSV_ANOM(:,:,:,JSV),        &
                       GUPDRAFT_MASK,                                   &
                       XLES_UPDRAFT_Sv2(:,NLES_CURRENT_TCOUNT,JSV)      )
    !
    CALL LES_FLUX_ll ( ZTH_ANOM, XSV_ANOM(:,:,:,JSV),                   &
                       GUPDRAFT_MASK,                                   &
                       XLES_UPDRAFT_ThSv(:,NLES_CURRENT_TCOUNT,JSV)     )
    !
    IF (LUSERC) &
    CALL LES_FLUX_ll ( XTHL_ANOM, XSV_ANOM(:,:,:,JSV),                  &
                       GUPDRAFT_MASK,                                   &
                       XLES_UPDRAFT_ThlSv(:,NLES_CURRENT_TCOUNT,JSV)    )
    !
    IF (LUSERV) &
    CALL LES_FLUX_ll ( ZTHV_ANOM, XSV_ANOM(:,:,:,JSV),                  &
                       GUPDRAFT_MASK,                                   &
                       XLES_UPDRAFT_ThvSv(:,NLES_CURRENT_TCOUNT,JSV)    )
  END DO
!
END IF
!
!-------------------------------------------------------------------------------
!
!       4.   Downdraft diagnostics
!            ---------------------
!
IF (LLES_DOWNDRAFT) THEN
!
  DO JK=1,NLES_K
    GDOWNDRAFT_MASK(:,:,JK) = (XW_ANOM(:,:,JK) <= 0.) .AND. LLES_CURRENT_CART_MASK(:,:,JK)
  END DO
!
!
!       4.1  Downdraft fraction
!            ------------------
!
  ZDOWNDRAFT(:,:,:) = 0.
  WHERE (GDOWNDRAFT_MASK(:,:,:))
    ZDOWNDRAFT(:,:,:) = 1.
  END WHERE
!
  CALL LES_MEAN_ll ( ZDOWNDRAFT, OMASK,                        &
                     XLES_DOWNDRAFT(:,NLES_CURRENT_TCOUNT)     )
!
!
!       4.2  Downdraft mean quantities
!            -------------------------
!
!* vertical wind velocity
!
  CALL LES_MEAN_ll ( ZW_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_W(:,NLES_CURRENT_TCOUNT)     )
!
!* potential temperature
!
  CALL LES_MEAN_ll ( ZTH_LES, GDOWNDRAFT_MASK,                   &
                     XLES_DOWNDRAFT_Th(:,NLES_CURRENT_TCOUNT)    )
!
!* liquid potential temperature
!
  IF (LUSERC) &
  CALL LES_MEAN_ll ( ZTHL_LES, GDOWNDRAFT_MASK,                   &
                     XLES_DOWNDRAFT_Thl(:,NLES_CURRENT_TCOUNT)    )
!
!* virtual potential temperature
!
  IF (LUSERV) &
  CALL LES_MEAN_ll ( ZTHV_LES, GDOWNDRAFT_MASK,                   &
                     XLES_DOWNDRAFT_Thv(:,NLES_CURRENT_TCOUNT)    )
!
!* vapor mixing ratio
!
  IF (LUSERV) &
  CALL LES_MEAN_ll ( ZRV_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_Rv(:,NLES_CURRENT_TCOUNT)     )
!
!* cloud water mixing ratio
!
  IF (LUSERC) &
  CALL LES_MEAN_ll ( ZRC_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_Rc(:,NLES_CURRENT_TCOUNT)     )
!
!* rain mixing ratio
!
  IF (LUSERR) &
  CALL LES_MEAN_ll ( ZRR_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_Rr(:,NLES_CURRENT_TCOUNT)     )
!
!* cloud ice mixing ratio
!
  IF (LUSERI) &
  CALL LES_MEAN_ll ( ZRI_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_Ri(:,NLES_CURRENT_TCOUNT)     )
!
!* snow mixing ratio
!
  IF (LUSERS) &
  CALL LES_MEAN_ll ( ZRS_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_Rs(:,NLES_CURRENT_TCOUNT)     )
!
!* graupel mixing ratio
!
  IF (LUSERG) &
  CALL LES_MEAN_ll ( ZRG_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_Rg(:,NLES_CURRENT_TCOUNT)     )
!
!* hail mixing ratio
!
  IF (LUSERH) &
  CALL LES_MEAN_ll ( ZRG_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_Rh(:,NLES_CURRENT_TCOUNT)     )
!
!* scalar variables
!
  DO JSV=1,NSV
    CALL LES_MEAN_ll ( ZSV_LES(:,:,:,JSV), GDOWNDRAFT_MASK,            &
                       XLES_DOWNDRAFT_Sv(:,NLES_CURRENT_TCOUNT,JSV)    )
  END DO
!
!* subgrid turbulent kinetic energy
!
  CALL LES_MEAN_ll ( ZTKE_LES, GDOWNDRAFT_MASK,                    &
                     XLES_DOWNDRAFT_Tke(:,NLES_CURRENT_TCOUNT)     )
!
!
!       4.3  Downdraft resolved quantities
!            -----------------------------
!
!* resolved turbulent kinetic energy
!
  CALL LES_FLUX_ll ( XU_ANOM, XU_ANOM,                             &
                     GDOWNDRAFT_MASK,                              &
                     ZLES_DOWNDRAFT_U2(:)                          )

  CALL LES_FLUX_ll ( XV_ANOM, XV_ANOM,                             &
                     GDOWNDRAFT_MASK,                              &
                     ZLES_DOWNDRAFT_V2(:)                          )

  CALL LES_FLUX_ll ( XW_ANOM, XW_ANOM,                             &
                     GDOWNDRAFT_MASK,                              &
                     ZLES_DOWNDRAFT_W2(:)                          )

  XLES_DOWNDRAFT_Ke(:,NLES_CURRENT_TCOUNT)  = 0.5 * (  ZLES_DOWNDRAFT_U2(:) &
                                                     + ZLES_DOWNDRAFT_V2(:) &
                                                     + ZLES_DOWNDRAFT_W2(:) )
!
!* vertical potential temperature flux
!
  CALL LES_FLUX_ll ( XW_ANOM, ZTH_ANOM,                            &
                     GDOWNDRAFT_MASK,                              &
                     XLES_DOWNDRAFT_WTh(:,NLES_CURRENT_TCOUNT)     )
!
!* vertical liquid potential temperature flux
!
  IF (LUSERC) &
  CALL LES_FLUX_ll ( XW_ANOM, XTHL_ANOM,                            &
                     GDOWNDRAFT_MASK,                               &
                     XLES_DOWNDRAFT_WThl(:,NLES_CURRENT_TCOUNT)     )
!
!* vertical virtual potential temperature flux
!
  IF (LUSERV) &
  CALL LES_FLUX_ll ( XW_ANOM, ZTHV_ANOM,                            &
                     GDOWNDRAFT_MASK,                               &
                     XLES_DOWNDRAFT_WThv(:,NLES_CURRENT_TCOUNT)     )
!
!* potential temperature variance
!
  CALL LES_FLUX_ll ( ZTH_ANOM, ZTH_ANOM,                            &
                     GDOWNDRAFT_MASK,                               &
                     XLES_DOWNDRAFT_Th2(:,NLES_CURRENT_TCOUNT)      )
!
!* liquid potential temperature variance
!
  IF (LUSERC) &
  CALL LES_FLUX_ll ( XTHL_ANOM, XTHL_ANOM,                           &
                     GDOWNDRAFT_MASK,                                &
                     XLES_DOWNDRAFT_Thl2(:,NLES_CURRENT_TCOUNT)      )
!
!* potential temperature - virtual potential temperature covariance
!
  IF (LUSERV) &
  CALL LES_FLUX_ll ( ZTH_ANOM, ZTHV_ANOM,                            &
                     GDOWNDRAFT_MASK,                                &
                     XLES_DOWNDRAFT_ThThv (:,NLES_CURRENT_TCOUNT)    )
!
!* liquid potential temperature - virtual potential temperature covariance
!
  IF (LUSERC) &
  CALL LES_FLUX_ll ( XTHL_ANOM, ZTHV_ANOM,                           &
                     GDOWNDRAFT_MASK,                                &
                     XLES_DOWNDRAFT_ThlThv(:,NLES_CURRENT_TCOUNT)    )
!
!
!* water vapor mixing ratio flux, variance and correlations
!
  IF (LUSERV) THEN
    CALL LES_FLUX_ll ( XW_ANOM, ZRV_ANOM,                             &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_WRv(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZRV_ANOM, ZRV_ANOM,                            &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_Rv2(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRV_ANOM,                            &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThRv (:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRV_ANOM,                           &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThvRv(:,NLES_CURRENT_TCOUNT)    )
    !
    IF (LUSERC) &
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRV_ANOM,                           &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThlRv(:,NLES_CURRENT_TCOUNT)    )
  END IF
!
!* cloud water mixing ratio flux
!
  IF (LUSERC) THEN
    CALL LES_FLUX_ll ( XW_ANOM, ZRC_ANOM,                             &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_WRc(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZRC_ANOM, ZRC_ANOM,                            &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_Rc2(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRC_ANOM,                            &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThRc (:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRC_ANOM,                           &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThvRc(:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRC_ANOM,                           &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThlRc(:,NLES_CURRENT_TCOUNT)    )
  END IF
!
!* cloud ice mixing ratio flux
!
  IF (LUSERI) THEN
    CALL LES_FLUX_ll ( XW_ANOM, ZRI_ANOM,                             &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_WRi(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZRI_ANOM, ZRI_ANOM,                            &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_Ri2(:,NLES_CURRENT_TCOUNT)      )
    !
    CALL LES_FLUX_ll ( ZTH_ANOM, ZRI_ANOM,                            &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThRi (:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( ZTHV_ANOM, ZRI_ANOM,                           &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThvRi(:,NLES_CURRENT_TCOUNT)    )
    !
    CALL LES_FLUX_ll ( XTHL_ANOM, ZRI_ANOM,                           &
                       GDOWNDRAFT_MASK,                               &
                       XLES_DOWNDRAFT_ThlRi(:,NLES_CURRENT_TCOUNT)    )
  END IF
!
!* scalar variables flux
!
  DO JSV=1,NSV
    CALL LES_FLUX_ll ( XW_ANOM, XSV_ANOM(:,:,:,JSV),                      &
                       GDOWNDRAFT_MASK,                                   &
                       XLES_DOWNDRAFT_WSv(:,NLES_CURRENT_TCOUNT,JSV)      )
    !
    CALL LES_FLUX_ll ( XSV_ANOM(:,:,:,JSV), XSV_ANOM(:,:,:,JSV),          &
                       GDOWNDRAFT_MASK,                                   &
                       XLES_DOWNDRAFT_Sv2(:,NLES_CURRENT_TCOUNT,JSV)      )
    !
    CALL LES_FLUX_ll ( ZTH_ANOM, XSV_ANOM(:,:,:,JSV),                     &
                       GDOWNDRAFT_MASK,                                   &
                       XLES_DOWNDRAFT_ThSv(:,NLES_CURRENT_TCOUNT,JSV)     )
    !
    IF (LUSERC) &
    CALL LES_FLUX_ll ( XTHL_ANOM, XSV_ANOM(:,:,:,JSV),                    &
                       GDOWNDRAFT_MASK,                                   &
                       XLES_DOWNDRAFT_ThlSv(:,NLES_CURRENT_TCOUNT,JSV)    )
    !
    IF (LUSERV) &
    CALL LES_FLUX_ll ( ZTHV_ANOM, XSV_ANOM(:,:,:,JSV),                    &
                       GDOWNDRAFT_MASK,                                   &
                       XLES_DOWNDRAFT_ThvSv(:,NLES_CURRENT_TCOUNT,JSV)    )
  END DO
!
END IF
!
!-------------------------------------------------------------------------------
!
!       5.   surface or 2D variables (only for the cartesian mask)
!            -----------------------
!
!* surface flux of temperature Qo
!
CALL LES_MEAN_MPROC ( XLES_Q0   (NLES_CURRENT_TCOUNT), IAVG_PTS(1), IUND_PTS(1) )
!
!* surface flux of water vapor Eo
!
CALL LES_MEAN_MPROC ( XLES_E0   (NLES_CURRENT_TCOUNT), IAVG_PTS(1), IUND_PTS(1) )
!
!* surface flux for scalar variables
!
DO JSV=1,NSV
  CALL LES_MEAN_MPROC ( XLES_SV0  (NLES_CURRENT_TCOUNT,JSV), IAVG_PTS(1), IUND_PTS(1) )
END DO
!
!* surface flux of U wind component
!
CALL LES_MEAN_MPROC ( XLES_UW0  (NLES_CURRENT_TCOUNT), IAVG_PTS(1), IUND_PTS(1) )
!
!* surface flux of V wind component
!
CALL LES_MEAN_MPROC ( XLES_VW0  (NLES_CURRENT_TCOUNT), IAVG_PTS(1), IUND_PTS(1) )
!
!* friction velocity u*
!
!* average of local u*
!!CALL LES_MEAN_MPROC ( XLES_USTAR(NLES_CURRENT_TCOUNT), IAVG_PTS(1), IUND_PTS(1) )
!* or true global u*
XLES_USTAR(NLES_CURRENT_TCOUNT) = SQRT(SQRT(XLES_UW0(NLES_CURRENT_TCOUNT)**2  &
                                           +XLES_VW0(NLES_CURRENT_TCOUNT)**2 ))
!
!* Boundary layer height
!
IF (CBL_HEIGHT_DEF=='WTV') THEN
!
!* level where temperature flux is minimum
!
ALLOCATE(ZWORK(SIZE(XLES_SUBGRID_WTHVMF(:,NLES_CURRENT_TCOUNT,IMASK),1)))
ZWORK=XLES_SUBGRID_WTHVMF(:,NLES_CURRENT_TCOUNT,IMASK)
WHERE(ZWORK==XUNDEF) ZWORK=0.

  IF (LUSERC) THEN
    IKMIN_FLUX = MINLOC(  XLES_RESOLVED_WThv(:,NLES_CURRENT_TCOUNT,1)  &
                        + XLES_SUBGRID_WThl (:,NLES_CURRENT_TCOUNT,1)  &
                        + ZWORK                                        & ! flux if EDKF   
      + (XRV/XRD - 1.) *( XLES_SUBGRID_WRt  (:,NLES_CURRENT_TCOUNT,1)  &
                         -XLES_SUBGRID_WRc  (:,NLES_CURRENT_TCOUNT,1)) )
  ELSE IF (LUSERV) THEN
    IKMIN_FLUX = MINLOC(  XLES_RESOLVED_WThv(:,NLES_CURRENT_TCOUNT,1) &
                        + ZWORK                                       & ! flux if EDKF 
                        + XLES_SUBGRID_WThl (:,NLES_CURRENT_TCOUNT,1) &
      + (XRV/XRD - 1.) *  XLES_SUBGRID_WRt  (:,NLES_CURRENT_TCOUNT,1) )
  ELSE
    IKMIN_FLUX = MINLOC(  XLES_RESOLVED_WTh(:,NLES_CURRENT_TCOUNT,1) &
                        + ZWORK                                      & ! flux if EDKF
                        + XLES_SUBGRID_WThl(:,NLES_CURRENT_TCOUNT,1) )
  END IF
DEALLOCATE(ZWORK)
!
!* boundary layer height
!
  XLES_BL_HEIGHT(NLES_CURRENT_TCOUNT) = XLES_Z(IKMIN_FLUX(1)) - XLES_ZS
!
ELSE IF (CBL_HEIGHT_DEF=='DTH') THEN
  IKMAX_TH=MAXLOC( ZLES_MEAN_DTHDZ(:)) 
  XLES_BL_HEIGHT(NLES_CURRENT_TCOUNT) = XLES_Z(IKMAX_TH(1)) - XLES_ZS
!  
ELSE IF (CBL_HEIGHT_DEF=='KE ') THEN

  XLES_BL_HEIGHT(NLES_CURRENT_TCOUNT) = XLES_Z(NLES_K) - XLES_ZS
!
!* total Turbulent Kinetic Energy
!
  ZKE_TOT(:) = 0.
!
  ZKE_TOT(:) = ZKE_TOT(:) + XLES_SUBGRID_TKE   (:,NLES_CURRENT_TCOUNT,1)
!
  IF (CTURBLEN/='BL89' .AND. CTURBLEN/='RM17' .AND. LLES_RESOLVED) &
  ZKE_TOT(:) = ZKE_TOT(:) + XLES_RESOLVED_KE(:,NLES_CURRENT_TCOUNT,1)
!
  ZINT_KE_TOT = 0.
!
!* integration of total kinetic energy on boundary layer depth
!
  ZINT_KE_TOT = ZINT_KE_TOT +XLES_Z(1)*ZKE_TOT(1)
  DO JK=1,NLES_K-1
    ZINT_KE_TOT = ZINT_KE_TOT + (XLES_Z(JK+1)-XLES_Z(JK))                     &
                              * 0.5 *( ZKE_TOT(JK+1) + ZKE_TOT(JK) )
!
!* test of total kinetic energy smaller than 5% of the averaged value below
!
    IF ( ZKE_TOT(JK+1) < 0.05 * ZINT_KE_TOT / (XLES_Z(JK+1)-XLES_Z(1)) ) THEN
      XLES_BL_HEIGHT(NLES_CURRENT_TCOUNT) = XLES_Z(JK) - XLES_ZS
      EXIT
    END IF
!
  END DO
!
ELSE IF (CBL_HEIGHT_DEF=='TKE') THEN

  XLES_BL_HEIGHT(NLES_CURRENT_TCOUNT) = XLES_Z(NLES_K) - XLES_ZS
!
!* subgrid Turbulent Kinetic Energy
!
  ZKE_TOT(:) = XLES_SUBGRID_TKE   (:,NLES_CURRENT_TCOUNT,1)
!
  ZINT_KE_TOT = 0.
!
!* integration of subgrid kinetic energy on boundary layer depth
!
  DO JK=1,NLES_K-1
    ZINT_KE_TOT = ZINT_KE_TOT + (XLES_Z(JK+1)-XLES_Z(JK))                     &
                              * 0.5 *( ZKE_TOT(JK+1) + ZKE_TOT(JK) )
!
!* test of subgrid kinetic energy smaller than 0.1% of the averaged value below
!
    IF ( ZKE_TOT(JK+1) < 0.001 * ZINT_KE_TOT / (XLES_Z(JK+1)-XLES_Z(1)) ) THEN
      XLES_BL_HEIGHT(NLES_CURRENT_TCOUNT) = XLES_Z(JK) - XLES_ZS
      EXIT
    END IF
  END DO
ELSE IF (CBL_HEIGHT_DEF=='FRI') THEN
  ZFRIC_LES = SQRT( ( XLES_SUBGRID_WU (:,NLES_CURRENT_TCOUNT,1)      &
                     +XLES_RESOLVED_WU(:,NLES_CURRENT_TCOUNT,1))**2  &
                   +( XLES_SUBGRID_WV (:,NLES_CURRENT_TCOUNT,1)      &
                     +XLES_RESOLVED_WV(:,NLES_CURRENT_TCOUNT,1))**2 )
  ZFRIC_SURF = XLES_USTAR(NLES_CURRENT_TCOUNT)**2
  CALL BL_DEPTH_DIAG(YLDIMPHYEX,ZFRIC_SURF, XLES_ZS, &
                     ZFRIC_LES,  XLES_Z,  &
                     XFTOP_O_FSURF,XLES_BL_HEIGHT(NLES_CURRENT_TCOUNT))
END IF
!
!
!* integration of total kinetic energy on boundary layer depth
!
XLES_INT_TKE(NLES_CURRENT_TCOUNT)=ZINT_KE_TOT
    !* integration of tke 
    ZTKET_LES(:,:)  = 0.
      DO JK=1,NLES_K-1
      ZKE_LES(:,:,JK)=0.5*(XU_ANOM(:,:,JK)*XU_ANOM(:,:,JK)+&
      XV_ANOM(:,:,JK)*XV_ANOM(:,:,JK)+XW_ANOM(:,:,JK)*XW_ANOM(:,:,JK))

    ZTKET_LES(:,:) = ZTKET_LES(:,:) + (ZZZ_LES(:,:,JK+1)-ZZZ_LES(:,:,JK))      &
                      * (ZTKE_LES(:,:,JK)+ZKE_LES(:,:,JK)) 
      END DO
  CALL LES_MEAN_ll ( ZTKET_LES, LLES_CURRENT_CART_MASK(:,:,1),               &
                    XLES_INT_TKE(NLES_CURRENT_TCOUNT)     )    
!
!* convective velocity
!
XLES_WSTAR(NLES_CURRENT_TCOUNT) = 0.
!
IF (               XLES_Q0(NLES_CURRENT_TCOUNT)                  &
    + (XRV/XRD-1.)*XLES_E0(NLES_CURRENT_TCOUNT) >0.) THEN
  IF (LUSERV) THEN
    XLES_WSTAR(NLES_CURRENT_TCOUNT) =                            &
     ( XG / XLES_MEAN_Thv (1,NLES_CURRENT_TCOUNT,1)              &
          * (                    XLES_Q0( NLES_CURRENT_TCOUNT )  &
              + (XRV/XRD - 1.) * XLES_E0( NLES_CURRENT_TCOUNT )) &
          * XLES_BL_HEIGHT(  NLES_CURRENT_TCOUNT  )              &
     ) ** (1./3.)
  ELSE
    XLES_WSTAR(NLES_CURRENT_TCOUNT) =                            &
     ( XG / XLES_MEAN_Th  (1,NLES_CURRENT_TCOUNT,1)              &
          * (                    XLES_Q0( NLES_CURRENT_TCOUNT )  &
              + (XRV/XRD - 1.) * XLES_E0( NLES_CURRENT_TCOUNT )) &
          * XLES_BL_HEIGHT( NLES_CURRENT_TCOUNT  )               &
     ) ** (1./3.)
  END IF
END IF
!
!* cloud base height
 IF (LUSERC) THEN
  ZINT_RHOKE =0.
  JJ=1
  DO JI=1,NLES_K
  IF ((ZINT_RHOKE .EQ. 0) .AND.  &
      (XLES_MEAN_RC(JI,NLES_CURRENT_TCOUNT,1) .GT. 1.E-6)) THEN
     ZINT_RHOKE=1.
     JJ=JI 
   END IF
  END DO
  XLES_ZCB(NLES_CURRENT_TCOUNT)= XLES_Z(JJ)-XLES_ZS
 ENDIF
!
!* height of max of cf
 IF (LUSERC) THEN
  IKMAX_CF= MAXLOC(  XLES_MEAN_INDCf(:,NLES_CURRENT_TCOUNT,1))
  XLES_ZMAXCF(NLES_CURRENT_TCOUNT) = XLES_Z(IKMAX_CF(1)) - XLES_ZS
  IKMAX_CF= MAXLOC(  XLES_MEAN_INDCf2(:,NLES_CURRENT_TCOUNT,1))
  XLES_ZMAXCF2(NLES_CURRENT_TCOUNT) = XLES_Z(IKMAX_CF(1)) - XLES_ZS
 ENDIF
!
!* Monin-Obukhov length
!
XLES_MO_LENGTH(NLES_CURRENT_TCOUNT) = 0.
!
IF (LUSERV) THEN
  IF ( XLES_Q0(NLES_CURRENT_TCOUNT)+(XRV/XRD-1.)*XLES_E0(NLES_CURRENT_TCOUNT) /=0. )&
  XLES_MO_LENGTH(NLES_CURRENT_TCOUNT) = (- (XLES_USTAR(NLES_CURRENT_TCOUNT))**3)  &
                    / (XKARMAN*(  XLES_Q0(NLES_CURRENT_TCOUNT)                    &
                                +(XRV/XRD-1.)*XLES_E0(NLES_CURRENT_TCOUNT))       &
                              *XG/XLES_MEAN_Thv(1,NLES_CURRENT_TCOUNT,1) )
ELSE
   IF ( XLES_Q0(NLES_CURRENT_TCOUNT) /=0. )                                       &
   XLES_MO_LENGTH(NLES_CURRENT_TCOUNT) = (- (XLES_USTAR(NLES_CURRENT_TCOUNT))**3) &
                    / (XKARMAN*XLES_Q0(NLES_CURRENT_TCOUNT)                       &
                              *XG/XLES_MEAN_Th(1,NLES_CURRENT_TCOUNT,1) )
END IF
!
!-------------------------------------------------------------------------------
!
!       6.   correlations along x and y axes
!            -------------------------------
!
!*  u * u
!
DO JK=1,NSPECTRA_K
  CALL LES_HOR_CORR( ZU_SPEC(:,:,JK), ZU_SPEC(:,:,JK),          &
                     CLES_LBCX , CLES_LBCY,                     &
                     XCORRi_UU(:,JK,NLES_CURRENT_TCOUNT),       &
                     XCORRj_UU(:,JK,NLES_CURRENT_TCOUNT)        )
END DO
!
!*  v * v
!
DO JK=1,NSPECTRA_K
  CALL LES_HOR_CORR( ZV_SPEC(:,:,JK), ZV_SPEC(:,:,JK),          &
                     CLES_LBCX , CLES_LBCY,                     &
                     XCORRi_VV(:,JK,NLES_CURRENT_TCOUNT),       &
                     XCORRj_VV(:,JK,NLES_CURRENT_TCOUNT)        )
END DO
!
!*  u * v
!
DO JK=1,NSPECTRA_K
  CALL LES_HOR_CORR( ZU_SPEC(:,:,JK), ZV_SPEC(:,:,JK),          &
                     CLES_LBCX , CLES_LBCY,                     &
                     XCORRi_UV(:,JK,NLES_CURRENT_TCOUNT),       &
                     XCORRj_UV(:,JK,NLES_CURRENT_TCOUNT)        )
END DO
!
!*  w * u
!
DO JK=1,NSPECTRA_K
  CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZU_SPEC(:,:,JK),          &
                     CLES_LBCX , CLES_LBCY,                     &
                     XCORRi_WU(:,JK,NLES_CURRENT_TCOUNT),       &
                     XCORRj_WU(:,JK,NLES_CURRENT_TCOUNT)        )
END DO
!
!*  w * v
!
DO JK=1,NSPECTRA_K
  CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZV_SPEC(:,:,JK),          &
                     CLES_LBCX , CLES_LBCY,                     &
                     XCORRi_WV(:,JK,NLES_CURRENT_TCOUNT),       &
                     XCORRj_WV(:,JK,NLES_CURRENT_TCOUNT)        )
END DO
!
!*  w * w
!
DO JK=1,NSPECTRA_K
  CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZW_SPEC(:,:,JK),          &
                     CLES_LBCX , CLES_LBCY,                     &
                     XCORRi_WW(:,JK,NLES_CURRENT_TCOUNT),       &
                     XCORRj_WW(:,JK,NLES_CURRENT_TCOUNT)        )
END DO
!
!*  w * th
!
DO JK=1,NSPECTRA_K
  CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZTH_SPEC(:,:,JK),         &
                     CLES_LBCX , CLES_LBCY,                     &
                     XCORRi_WTh(:,JK,NLES_CURRENT_TCOUNT),      &
                     XCORRj_WTh(:,JK,NLES_CURRENT_TCOUNT)       )
END DO
!
!*  w * thl
!
DO JK=1,NSPECTRA_K
  IF (LUSERC) &
  CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZTHL_SPEC(:,:,JK),         &
                     CLES_LBCX , CLES_LBCY,                      &
                     XCORRi_WThl(:,JK,NLES_CURRENT_TCOUNT),      &
                     XCORRj_WThl(:,JK,NLES_CURRENT_TCOUNT)       )
END DO
!
!*  th * th
!
DO JK=1,NSPECTRA_K
  CALL LES_HOR_CORR( ZTH_SPEC(:,:,JK), ZTH_SPEC(:,:,JK),          &
                     CLES_LBCX , CLES_LBCY,                       &
                     XCORRi_ThTh(:,JK,NLES_CURRENT_TCOUNT),       &
                     XCORRj_ThTh(:,JK,NLES_CURRENT_TCOUNT)        )
END DO
!
!*  thl * thl
!
DO JK=1,NSPECTRA_K
  IF (LUSERC) &
  CALL LES_HOR_CORR( ZTHL_SPEC(:,:,JK), ZTHL_SPEC(:,:,JK),        &
                     CLES_LBCX , CLES_LBCY,                       &
                     XCORRi_ThlThl(:,JK,NLES_CURRENT_TCOUNT),     &
                     XCORRj_ThlThl(:,JK,NLES_CURRENT_TCOUNT)      )
END DO
!
!* correlations with water vapor
!
IF (LUSERV) THEN
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZRV_SPEC(:,:,JK),        &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_WRv(:,JK,NLES_CURRENT_TCOUNT),     &
                       XCORRj_WRv(:,JK,NLES_CURRENT_TCOUNT)      )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZTH_SPEC(:,:,JK), ZRV_SPEC(:,:,JK),       &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_ThRv(:,JK,NLES_CURRENT_TCOUNT),    &
                       XCORRj_ThRv(:,JK,NLES_CURRENT_TCOUNT)     )
  END DO
  !
  DO JK=1,NSPECTRA_K
    IF (LUSERC) &
    CALL LES_HOR_CORR( ZTHL_SPEC(:,:,JK), ZRV_SPEC(:,:,JK),      &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_ThlRv(:,JK,NLES_CURRENT_TCOUNT),   &
                       XCORRj_ThlRv(:,JK,NLES_CURRENT_TCOUNT)    )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZRV_SPEC(:,:,JK), ZRV_SPEC(:,:,JK),       &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_RvRv(:,JK,NLES_CURRENT_TCOUNT),    &
                       XCORRj_RvRv(:,JK,NLES_CURRENT_TCOUNT)     )
  END DO
END IF
!
!
!* correlations with cloud water
!
IF (LUSERC) THEN
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZRC_SPEC(:,:,JK),        &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_WRc(:,JK,NLES_CURRENT_TCOUNT),     &
                       XCORRj_WRc(:,JK,NLES_CURRENT_TCOUNT)      )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZTH_SPEC(:,:,JK), ZRC_SPEC(:,:,JK),       &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_ThRc(:,JK,NLES_CURRENT_TCOUNT),    &
                       XCORRj_ThRc(:,JK,NLES_CURRENT_TCOUNT)     )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZTHL_SPEC(:,:,JK), ZRC_SPEC(:,:,JK),      &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_ThlRc(:,JK,NLES_CURRENT_TCOUNT),   &
                       XCORRj_ThlRc(:,JK,NLES_CURRENT_TCOUNT)    )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZRC_SPEC(:,:,JK), ZRC_SPEC(:,:,JK),       &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_RcRc(:,JK,NLES_CURRENT_TCOUNT),    &
                       XCORRj_RcRc(:,JK,NLES_CURRENT_TCOUNT)     )
  END DO
END IF
!
!* correlations with cloud ice
!
IF (LUSERI) THEN
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZRI_SPEC(:,:,JK),        &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_WRi(:,JK,NLES_CURRENT_TCOUNT),     &
                       XCORRj_WRi(:,JK,NLES_CURRENT_TCOUNT)      )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZTH_SPEC(:,:,JK), ZRI_SPEC(:,:,JK),       &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_ThRi(:,JK,NLES_CURRENT_TCOUNT),    &
                       XCORRj_ThRi(:,JK,NLES_CURRENT_TCOUNT)     )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZTHL_SPEC(:,:,JK), ZRI_SPEC(:,:,JK),      &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_ThlRi(:,JK,NLES_CURRENT_TCOUNT),   &
                       XCORRj_ThlRi(:,JK,NLES_CURRENT_TCOUNT)    )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZRI_SPEC(:,:,JK), ZRI_SPEC(:,:,JK),       &
                       CLES_LBCX , CLES_LBCY,                    &
                       XCORRi_RiRi(:,JK,NLES_CURRENT_TCOUNT),    &
                       XCORRj_RiRi(:,JK,NLES_CURRENT_TCOUNT)     )
  END DO
END IF
!
!* correlations with scalar variables
!
DO JSV=1,NSV
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZW_SPEC(:,:,JK), ZSV_SPEC(:,:,JK,JSV),       &
                       CLES_LBCX , CLES_LBCY,                       &
                       XCORRi_WSv(:,JK,NLES_CURRENT_TCOUNT,JSV),    &
                       XCORRj_WSv(:,JK,NLES_CURRENT_TCOUNT,JSV)     )
  END DO
  !
  DO JK=1,NSPECTRA_K
    CALL LES_HOR_CORR( ZSV_SPEC(:,:,JK,JSV), ZSV_SPEC(:,:,JK,JSV),  &
                       CLES_LBCX , CLES_LBCY,                       &
                       XCORRi_SvSv(:,JK,NLES_CURRENT_TCOUNT,JSV),   &
                       XCORRj_SvSv(:,JK,NLES_CURRENT_TCOUNT,JSV)    )
  END DO
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LES_n
