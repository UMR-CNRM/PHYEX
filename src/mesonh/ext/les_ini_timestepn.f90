!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      #######################
MODULE MODI_LES_INI_TIMESTEP_n
!      #######################
!
!
INTERFACE LES_INI_TIMESTEP_n
!
      SUBROUTINE LES_INI_TIMESTEP_n(KTCOUNT)
!
INTEGER, INTENT(IN) :: KTCOUNT ! current model time-step
!
END SUBROUTINE LES_INI_TIMESTEP_n
!
END INTERFACE
!
END MODULE MODI_LES_INI_TIMESTEP_n

!     ##############################
      SUBROUTINE  LES_INI_TIMESTEP_n(KTCOUNT)
!     ##############################
!
!
!!****  *LES_INI_TIMESTEP_n* initializes the LES variables for
!!                    the current time-step of model _n
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
!!      Original         06/11/02
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  P. Wautelet 30/03/2021: budgets: LES cartesian subdomain limits are defined in the physical domain
! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_NSV
USE MODD_LES
USE MODD_LES_n
USE MODD_FIELD_n
USE MODD_METRICS_n
USE MODD_REF_n
USE MODD_CONF_n
USE MODD_TIME_n
USE MODD_DYN_n
USE MODD_TIME
USE MODD_CONF
USE MODD_LES_BUDGET
!
use mode_datetime,       only: Datetime_distance
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODI_LES_VER_INT
USE MODI_THL_RT_FROM_TH_R
USE MODI_LES_MEAN_ll
USE MODI_SHUMAN
!
USE MODI_SECOND_MNH
USE MODI_LES_CLOUD_MASKS_N
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
INTEGER, INTENT(IN) :: KTCOUNT ! current model time-step
!
!
!       0.2  declaration of local variables
!
INTEGER :: IXOR_ll, IYOR_ll                    ! origine point coordinates
!                                              ! of current processor domain
!                                              ! on model domain on all
!                                              ! processors
INTEGER :: IIB_ll, IJB_ll                      ! SO point coordinates of
!                                              ! current processor phys. domain
!                                              ! on model domain on all
!                                              ! processors
INTEGER :: IIE_ll, IJE_ll                      ! NE point coordinates of
!                                              ! current processor phys. domain
!                                              ! on model domain on all
!                                              ! processors
INTEGER :: IIINF_MASK, IISUP_MASK              ! cart. mask local proc. limits
INTEGER :: IJINF_MASK, IJSUP_MASK              ! cart. mask local proc. limits
!
INTEGER :: JK                                  ! vertical loop counter
INTEGER :: IIB, IJB, IIE, IJE                  ! hor. indices
INTEGER :: IIU, IJU                            ! hor. indices
INTEGER :: IKU                                 ! ver. index
INTEGER :: IRR, IRRC, IRRR, IRRI, IRRS, IRRG   ! moist variables indices
!
INTEGER :: JSV                                 ! scalar variables counter
!
REAL    :: ZTIME1, ZTIME2                      ! CPU time counters
!
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZTHL ! theta_l
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZRT  ! total water
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZL   ! Latent heat of vaporization
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZCP  ! Cp
REAL, DIMENSION(:,:,:),   ALLOCATABLE :: ZEXN ! Exner function
INTEGER                               :: IMI  ! current model index
!-------------------------------------------------------------------------------
!
!*      1.   Does current time-step is a LES time-step?
!            -----------------------------------------
!
LLES_CALL= .FALSE.
!
CALL SECOND_MNH(ZTIME1)
!
IF (NLES_TCOUNT==NLES_TIMES) LLES_CALL=.FALSE.
!
IF ( KTCOUNT>1 .AND. MOD (KTCOUNT-1,NLES_DTCOUNT)==0) LLES_CALL=.TRUE.
!
IF (.NOT. LLES_CALL) RETURN
!
CALL BUDGET_FLAGS(LUSERV, LUSERC, LUSERR,         &
                  LUSERI, LUSERS, LUSERG, LUSERH  )
!
NLES_TCOUNT = NLES_TCOUNT + 1
!
NLES_CURRENT_TCOUNT = NLES_TCOUNT
!
tles_dates(nles_tcount) = tdtcur
call Datetime_distance( tdtseg, tdtcur, xles_times(nles_tcount) )
!
!* forward-in-time time-step
!
XCURRENT_TSTEP = XTSTEP
!
!-------------------------------------------------------------------------------
!
CALL GET_OR_ll     ('B',IXOR_ll,IYOR_ll)
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
IIB_ll=IXOR_ll+IIB-1
IJB_ll=IYOR_ll+IJB-1
IIE_ll=IXOR_ll+IIE-1
IJE_ll=IYOR_ll+IJE-1
!
IKU = SIZE(XTHT,3)
!
IMI = GET_CURRENT_MODEL_INDEX()
!
!-------------------------------------------------------------------------------
!
!*      2.   Definition of masks
!            -------------------
!
!*      2.1  Cartesian (sub-)domain (on local processor)
!            ----------------------
!
CALL LES_ALLOCATE('LLES_CURRENT_CART_MASK',(/IIU,IJU,NLES_K/))
!
IIINF_MASK = MAX(IIB, NLESn_IINF(IMI)+JPHEXT-(IIB_ll-1-JPHEXT))
IJINF_MASK = MAX(IJB, NLESn_JINF(IMI)+JPHEXT-(IJB_ll-1-JPHEXT))
IISUP_MASK = MIN(IIE, NLESn_ISUP(IMI)+JPHEXT-(IIB_ll-1-JPHEXT))
IJSUP_MASK = MIN(IJE, NLESn_JSUP(IMI)+JPHEXT-(IJB_ll-1-JPHEXT))
!
!
LLES_CURRENT_CART_MASK(:,:,:) = .FALSE.
LLES_CURRENT_CART_MASK(IIINF_MASK:IISUP_MASK,IJINF_MASK:IJSUP_MASK,:) = .TRUE.
!
CLES_CURRENT_LBCX(:) = CLES_LBCX(:,IMI)
CLES_CURRENT_LBCY(:) = CLES_LBCY(:,IMI)
!
!-------------------------------------------------------------------------------
!
!*      3.   Definition of LES vertical grid for this model
!            ----------------------------------------------
!
IF (CLES_LEVEL_TYPE=='Z') THEN
  IF (ASSOCIATED(XCOEFLIN_CURRENT_LES)) CALL LES_DEALLOCATE('XCOEFLIN_CURRENT_LES')
  IF (ASSOCIATED(NKLIN_CURRENT_LES   )) CALL LES_DEALLOCATE('NKLIN_CURRENT_LES')
  !
  CALL LES_ALLOCATE('XCOEFLIN_CURRENT_LES',(/IIU,IJU,NLES_K/))
  CALL LES_ALLOCATE('NKLIN_CURRENT_LES',(/IIU,IJU,NLES_K/))
  !
  XCOEFLIN_CURRENT_LES(:,:,:) = XCOEFLIN_LES(:,:,:)
  NKLIN_CURRENT_LES   (:,:,:) = NKLIN_LES   (:,:,:)
END IF
!
!-------------------------------------------------------------------------------
!
!*      4.   Definition of variables used in budgets for current model
!            ---------------------------------------------------------
!
IF (LUSERC) THEN
  ALLOCATE(XCURRENT_L_O_EXN_CP (IIU,IJU,IKU))
ELSE
  ALLOCATE(XCURRENT_L_O_EXN_CP (0,0,0))
END IF
ALLOCATE(XCURRENT_RHODJ      (IIU,IJU,IKU))
!
!* coefficients for Th to Thl conversion
!
IF (LUSERC) THEN
  ALLOCATE(ZL  (IIU,IJU,IKU))
  ALLOCATE(ZEXN(IIU,IJU,IKU))
  ALLOCATE(ZCP (IIU,IJU,IKU))
  !
  !* Exner function
  !
  ZEXN(:,:,:) = (XPABST/XP00)**(XRD/XCPD)
  !
  !* Latent heat of vaporization
  !
  ZL(:,:,:) = XLVTT + (XCPD-XCL) * (XTHT(:,:,:)*ZEXN(:,:,:)-XTT)
  !
  !* heat capacity at constant pressure of the humid air
  !
  ZCP(:,:,:) = XCPD
  IRR=2
  ZCP(:,:,:) = ZCP(:,:,:) + XCPV * XRT(:,:,:,1)
  ZCP(:,:,:) = ZCP(:,:,:) + XCL  * XRT(:,:,:,2)
  IF (LUSERR) THEN
    IRR=IRR+1
    ZCP(:,:,:) = ZCP(:,:,:) + XCL  * XRT(:,:,:,IRR)
  END IF
  IF (LUSERI) THEN
    IRR=IRR+1
    ZCP(:,:,:) = ZCP(:,:,:) + XCI  * XRT(:,:,:,IRR)
  END IF
  IF (LUSERS) THEN
    IRR=IRR+1
    ZCP(:,:,:) = ZCP(:,:,:) + XCI  * XRT(:,:,:,IRR)
  END IF
  IF (LUSERG) THEN
    IRR=IRR+1
    ZCP(:,:,:) = ZCP(:,:,:) + XCI  * XRT(:,:,:,IRR)
  END IF
  IF (LUSERH) THEN
    IRR=IRR+1
    ZCP(:,:,:) = ZCP(:,:,:) + XCI  * XRT(:,:,:,IRR)
  END IF
  !
  !* L / (Exn * Cp)
  !
  XCURRENT_L_O_EXN_CP(:,:,:) = ZL(:,:,:) / ZEXN(:,:,:) / ZCP(:,:,:)
  !
  DEALLOCATE(ZL  )
  DEALLOCATE(ZEXN)
  DEALLOCATE(ZCP )
END IF
!
!* other initializations
!
XCURRENT_RHODJ=XRHODJ
!
LCURRENT_USERV=LUSERV
LCURRENT_USERC=LUSERC
LCURRENT_USERR=LUSERR
LCURRENT_USERI=LUSERI
LCURRENT_USERS=LUSERS
LCURRENT_USERG=LUSERG
LCURRENT_USERH=LUSERH
!
NCURRENT_RR = NRR
!
ALLOCATE(XCURRENT_RUS  (IIU,IJU,IKU))
ALLOCATE(XCURRENT_RVS  (IIU,IJU,IKU))
ALLOCATE(XCURRENT_RWS  (IIU,IJU,IKU))
ALLOCATE(XCURRENT_RTHS (IIU,IJU,IKU))
ALLOCATE(XCURRENT_RTKES(IIU,IJU,IKU))
ALLOCATE(XCURRENT_RRS  (IIU,IJU,IKU,NRR))
ALLOCATE(XCURRENT_RSVS (IIU,IJU,IKU,NSV))
ALLOCATE(XCURRENT_RTHLS(IIU,IJU,IKU))
ALLOCATE(XCURRENT_RRTS (IIU,IJU,IKU))
!
XCURRENT_RUS  =XRUS
XCURRENT_RVS  =XRVS
XCURRENT_RWS  =XRWS
XCURRENT_RTHS =XRTHS
XCURRENT_RTKES=XRTKES
XCURRENT_RRS  =XRRS
XCURRENT_RSVS =XRSVS
CALL THL_RT_FROM_TH_R(LUSERV, LUSERC, LUSERR,             &
                      LUSERI, LUSERS, LUSERG, LUSERH,     &
                      XCURRENT_L_O_EXN_CP,                &
                      XCURRENT_RTHS, XCURRENT_RRS,        &
                      XCURRENT_RTHLS, XCURRENT_RRTS       )

ALLOCATE(X_LES_BU_RES_KE   (NLES_K,NLES_TOT))
ALLOCATE(X_LES_BU_RES_WThl (NLES_K,NLES_TOT))
ALLOCATE(X_LES_BU_RES_Thl2 (NLES_K,NLES_TOT))
ALLOCATE(X_LES_BU_SBG_Tke  (NLES_K,NLES_TOT))
ALLOCATE(X_LES_BU_RES_WRt  (NLES_K,NLES_TOT))
ALLOCATE(X_LES_BU_RES_Rt2  (NLES_K,NLES_TOT))
ALLOCATE(X_LES_BU_RES_ThlRt(NLES_K,NLES_TOT))
ALLOCATE(X_LES_BU_RES_Sv2  (NLES_K,NLES_TOT,NSV))
ALLOCATE(X_LES_BU_RES_WSv  (NLES_K,NLES_TOT,NSV))

X_LES_BU_RES_KE   = 0.
X_LES_BU_RES_WThl = 0.
X_LES_BU_RES_Thl2 = 0.
X_LES_BU_SBG_Tke  = 0.
X_LES_BU_RES_WRt  = 0.
X_LES_BU_RES_Rt2  = 0.
X_LES_BU_RES_ThlRt= 0.
X_LES_BU_RES_Sv2  = 0.
X_LES_BU_RES_WSv  = 0.
!
!-------------------------------------------------------------------------------
!
!*      4.   Definition of anomaly fields
!            ----------------------------
!
ALLOCATE (XU_ANOM  (IIU,IJU,NLES_K))
ALLOCATE (XV_ANOM  (IIU,IJU,NLES_K))
ALLOCATE (XW_ANOM  (IIU,IJU,NLES_K))
ALLOCATE (XTHL_ANOM(IIU,IJU,NLES_K))
IF (LUSERV) THEN
  ALLOCATE (XRT_ANOM (IIU,IJU,NLES_K))
ELSE
  ALLOCATE (XRT_ANOM (0,0,0))
END IF
ALLOCATE (XSV_ANOM (IIU,IJU,NLES_K,NSV))
!
!*      4.1  conservative variables
!            ----------------------
!
ALLOCATE(ZTHL(IIU,IJU,IKU))
ALLOCATE(ZRT (IIU,IJU,IKU))
CALL THL_RT_FROM_TH_R(LUSERV, LUSERC, LUSERR,             &
                      LUSERI, LUSERS, LUSERG, LUSERH,     &
                      XCURRENT_L_O_EXN_CP,                &
                      XTHT, XRT,                          &
                      ZTHL, ZRT                           )
!
!*      4.2  anomaly fields on the LES grid
!            ------------------------------
!
CALL LES_ANOMALY_FIELD(MXF(XUT),XU_ANOM)
CALL LES_ANOMALY_FIELD(MYF(XVT),XV_ANOM)
CALL LES_ANOMALY_FIELD(MZF(XWT),XW_ANOM)
CALL LES_ANOMALY_FIELD(ZTHL,XTHL_ANOM)
IF (LUSERV) CALL LES_ANOMALY_FIELD(ZRT,XRT_ANOM)
DO JSV=1,NSV
  CALL LES_ANOMALY_FIELD(XSVT(:,:,:,JSV),XSV_ANOM(:,:,:,JSV))
END DO
!
!-------------------------------------------------------------------------------
!
DEALLOCATE(ZTHL)
DEALLOCATE(ZRT )
!-------------------------------------------------------------------------------
!
!*      6.0  Nebulosity masks
!            ----------------
!
CALL LES_CLOUD_MASKS_n
!
!-------------------------------------------------------------------------------
CALL SECOND_MNH(ZTIME2)
XTIME_LES_BU = XTIME_LES_BU + ZTIME2 - ZTIME1
!--------------------------------------------------------------------------------
!
CONTAINS
!
!--------------------------------------------------------------------------------
!
SUBROUTINE LES_ANOMALY_FIELD(PF,PF_ANOM)
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PF
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PF_ANOM

REAL, DIMENSION(SIZE(PF_ANOM,3)) :: ZMEAN
INTEGER :: JI, JJ

CALL LES_VER_INT(PF, PF_ANOM)
CALL LES_MEAN_ll(PF_ANOM, LLES_CURRENT_CART_MASK, ZMEAN  )
DO JJ=1,SIZE(PF_ANOM,2)
  DO JI=1,SIZE(PF_ANOM,1)
    PF_ANOM(JI,JJ,:) = PF_ANOM(JI,JJ,:) - ZMEAN(:)
  END DO
END DO

END SUBROUTINE LES_ANOMALY_FIELD
!--------------------------------------------------------------------------------
!
END SUBROUTINE LES_INI_TIMESTEP_n   

