!     ######spl
      SUBROUTINE  ARO_LIMA(PHYEX,KKA,KKU,KKL,KLON,KLEV,KIDIA,KFDIA,KRR, KSV, &
                                  PTSTEP, PDZZ, PRHODJ, PRHODREF, PEXNREF,&
                                  PPABSM, PW_NU, PDTHRAD, PTHT, PRT, PSVT, PCIT, &
                                  PTHS, PRS, PSVS, PEVAP,  &
                                  PINPRR,PINPRS,                 &
                                  PINPRG,PINPRH,PFPR,     &
                                  PCLDFR,PICEFR,PPRCFR,         &
                                  PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF, &
                                  YDDDH, YDLDDH, YDMDDH    )

      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ##########################################################################
!
!!****  * -  compute the  resolved clouds and precipitation
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the  microphysical sources
!!    related to the resolved clouds and precipitation in LIMA
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    Vie et al., 2015 GMD
!!
!!    AUTHOR
!!    ------
!!    B. Vie 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    17/09/13
!!
!! 02/2025 - S. Antoine : Correction of negative values dependent on the number of moments of each species
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
!USE MODD_CONF
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_PARAMETERS
!
USE MODD_PARAM_LIMA
USE MODD_NSV
USE MODD_CH_AEROSOL, ONLY: NSP, NCARB, NSOA, LORILAM
USE MODD_DUST, ONLY: LDUST
USE MODD_SALT, ONLY: LSALT
!
USE MODD_BUDGET
USE MODE_BUDGET_IAL, ONLY: BUDGET_DDH, TBUDGETDATA_IAL
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX 
!
USE MODI_LIMA
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
TYPE(PHYEX_t),            INTENT(IN)   :: PHYEX
INTEGER,                  INTENT(IN)   :: KKA  !near ground array index
INTEGER,                  INTENT(IN)   :: KKU  !uppest atmosphere array index
INTEGER,                  INTENT(IN)   :: KKL  !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KIDIA    !
INTEGER,                  INTENT(IN)   :: KFDIA    !
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSV      ! Number of LIMA variables
REAL,                     INTENT(IN)   :: PTSTEP   ! Time step
!
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)      :: PDZZ     ! Height (z)
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)      :: PRHODJ   ! Dry density * Jacobian
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)      :: PRHODREF ! Reference dry air density
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)      :: PEXNREF  ! Reference Exner function
!
!
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)      :: PPABSM   ! abs. pressure at time t-dt
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)      :: PW_NU    ! w for CCN activation
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)      :: PDTHRAD  ! radiative Theta tendency for CCN act.
REAL, DIMENSION(KLON,KLEV),   INTENT(IN)      :: PTHT     ! Theta at time t
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(INOUT) :: PRT      ! Moist variables at time t
REAL, DIMENSION(KLON,KLEV,KSV), INTENT(INOUT) :: PSVT     ! LIMA variables at time t
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT)   :: PCIT  ! Pristine ice number
!
!
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT)   :: PTHS     ! Theta source
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(INOUT) :: PRS      ! Moist  variable sources
REAL, DIMENSION(KLON,KLEV,KSV), INTENT(INOUT) :: PSVS     ! LIMA variable sources
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT)     :: PEVAP    ! Rain evap profile
!
REAL, DIMENSION(KLON), INTENT(INOUT)          :: PINPRR   ! Rain instant precip
REAL, DIMENSION(KLON), INTENT(INOUT)          :: PINPRS   ! Snow instant precip
REAL, DIMENSION(KLON), INTENT(INOUT)          :: PINPRG   ! Graupel instant precip
REAL, DIMENSION(KLON), INTENT(INOUT)          :: PINPRH   ! Hail instant precip
REAL, DIMENSION(KLON,KLEV,KRR), INTENT(INOUT) :: PFPR     ! upper-air precip
!
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT)   :: PCLDFR   ! liquid cloud fraction
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT)   :: PICEFR   ! ice cloud fraction
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT)   :: PPRCFR   ! precipitation fraction
!
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT):: PHLC_HRC
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT):: PHLC_HCF
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT):: PHLI_HRI
REAL, DIMENSION(KLON,KLEV),   INTENT(INOUT):: PHLI_HCF
!
TYPE(TYP_DDH), INTENT(INOUT), TARGET :: YDDDH
TYPE(TLDDH), INTENT(IN), TARGET :: YDLDDH
TYPE(TMDDH), INTENT(IN), TARGET :: YDMDDH
!
!
!*       0.2   Declarations of local variables :

CHARACTER(LEN=4)   :: HACTCCN  ! kind of CCN activation
REAL, DIMENSION(KLON,KLEV,10)    :: ZSOLORG ![%] solubility fraction of soa
REAL, DIMENSION(KLON,KLEV, NSP+NCARB+NSOA) :: ZMI
!
!
INTEGER :: JRR, JSV, JL           ! Loop index for the moist and scalar variables
INTEGER :: JLON, JLEV
!
REAL, DIMENSION(KLON,KLEV) :: ZT
REAL, DIMENSION(KLON,KLEV) :: ZLV
REAL, DIMENSION(KLON,KLEV) :: ZLS
REAL, DIMENSION(KLON,KLEV) :: ZCPH

REAL, DIMENSION(KLON,KLEV) :: ZCOR
REAL, DIMENSION(KLON,KLEV) :: ZDUM3DC
REAL, DIMENSION(KLON,KLEV) :: ZDUM3DR
REAL, DIMENSION(KLON,KLEV) :: ZDUM3DS
REAL, DIMENSION(KLON,KLEV) :: ZDUM3DG
REAL, DIMENSION(KLON,KLEV) :: ZDUM3DH

REAL, DIMENSION(KLON,KLEV) :: ZRAINFR
REAL, DIMENSION(KLON,KLEV) :: ZHLC_HCF
REAL, DIMENSION(KLON,KLEV) :: ZHLC_LCF
REAL, DIMENSION(KLON,KLEV) :: ZHLC_HRC
REAL, DIMENSION(KLON,KLEV) :: ZHLC_LRC

REAL, DIMENSION(KLON):: ZINPRC    ! surf cloud sedimentation
                                    ! for the correction of negative rv
REAL, DIMENSION(KLON):: ZINPRI    ! surf cloud ice sedimentation
REAL, DIMENSION(KLON):: ZINDEP    ! surf cloud ice sedimentation
REAL  :: ZMASSTOT                   ! total mass  for one water category
                                    ! including the negative values
REAL  :: ZMASSPOS                   ! total mass  for one water category
                                    ! after removing the negative values
REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR

REAL, DIMENSION(KLON,KLEV) :: ZWORK

REAL  :: ZTHVREFZIKB
LOGICAL :: LL_RRR_BUDGET
!
TYPE(TBUDGETDATA_IAL), DIMENSION(NBUDGET_SV1+NSV_LIMA-1), TARGET :: YLBUDGET
TYPE(TBUDGETDATA_PTR), DIMENSION(NBUDGET_SV1+NSV_LIMA-1) :: YLBUDGET_PTR
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_LIMA',0,ZHOOK_HANDLE)

!Dimensions
CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, 0, KIDIA, KFDIA)

ZINPRC(:)    = 0.
ZDUM3DC(:,:) = 0.
ZDUM3DR(:,:) = 0.
ZDUM3DS(:,:) = 0.
ZDUM3DG(:,:) = 0.
ZDUM3DH(:,:) = 0.
PINPRH(:)    = 0.
HACTCCN='    '
ZSOLORG=0.
ZMI=0.
IF (PHYEX%MISC%CELEC /='NONE') THEN
  CALL ABOR1('ARO_LIMA : CELEC ELECTRICITY SCHEME NOT YET CORRECLY PLUGGED HERE')
ELSE
  ! The following value of ZTHVREFZIKB must be removed from the electricity scheme or computed correctly here
  ZTHVREFZIKB = 0. ! for electricity use only
END IF

!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
!
!  complete the vertical boundaries
!
!
! personal comment:  tranfering these variables to the
!                    microphysical routines would save
!                    computing time
!
DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    ZT(JLON,JLEV) = PTHT(JLON,JLEV)*PEXNREF(JLON,JLEV)
  ENDDO
ENDDO

DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    ZLV(JLON,JLEV) = PHYEX%CST%XLVTT + (PHYEX%CST%XCPV-PHYEX%CST%XCL) * (ZT(JLON,JLEV)-PHYEX%CST%XTT)
  ENDDO
ENDDO

DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    ZLS(JLON,JLEV) = PHYEX%CST%XLSTT + (PHYEX%CST%XCPV-PHYEX%CST%XCI) * (ZT(JLON,JLEV)-PHYEX%CST%XTT)
  ENDDO
ENDDO

DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    ZCPH(JLON,JLEV) = PHYEX%CST%XCPD + PHYEX%CST%XCPV*2.*PTSTEP*PRS(JLON,JLEV,1)
  ENDDO
ENDDO
!

!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for 1-moment precipitating species (Rood 87)
!
DO JRR = 5, MIN(KRR,7) ! snow, graupel and hail
  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      IF (PRS(JLON,JLEV,JRR) < 1.E-15)  PRS(JLON,JLEV,JRR) = 0.
    ENDDO
  ENDDO
ENDDO

!
!*       3.2    Correct negative values
!
! Correction where rc<0
DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (NMOM_C.EQ.1 .AND. PRS(JLON,JLEV,2) < 1.E-15) THEN
      PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,2)
      PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                      & - PRS(JLON,JLEV,2) * ZLV(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
      PRS(JLON,JLEV,2)  = 0.0
    ELSEIF (NMOM_C .GE. 2 .AND. (PRS(JLON,JLEV,2) < 1.E-15 .OR. PSVS(JLON,JLEV,NSV_LIMA_NC) < 1.E-15)) THEN
      PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,2)
      PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                      & - PRS(JLON,JLEV,2) * ZLV(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
      PRS(JLON,JLEV,2)  = 0.0
      PSVS(JLON,JLEV,NSV_LIMA_NC) = 0.0
    ENDIF
  ENDDO
ENDDO

! Correction where rr<0
DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (NMOM_R .EQ. 1 .AND. PRS(JLON,JLEV,3) < 1.E-15) THEN
      PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,3)
      PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                      & - PRS(JLON,JLEV,3) * ZLV(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
      PRS(JLON,JLEV,3) = 0.0
    ELSEIF (NMOM_R .GE. 2 .AND. (PRS(JLON,JLEV,3) < 1.E-15 .OR. PSVS(JLON,JLEV,NSV_LIMA_NR) < 1.E-15)) THEN
      PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,3)
      PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                      & - PRS(JLON,JLEV,3) * ZLV(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
      PRS(JLON,JLEV,3) = 0.0
      PSVS(JLON,JLEV,NSV_LIMA_NR) = 0.0
    ENDIF
  ENDDO
ENDDO

! Correction where ri<0
DO JLEV = 1, KLEV
  DO JLON = KIDIA, KFDIA
    IF (NMOM_I .EQ. 1 .AND. PRS(JLON,JLEV,4) < 1.E-15) THEN
      PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,4)
      PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                      & - PRS(JLON,JLEV,4) * ZLS(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
      PRS(JLON,JLEV,4)  = 0.0
    ELSEIF (NMOM_I .GE. 2 .AND. (PRS(JLON,JLEV,4) < 1.E-15 .OR. PSVS(JLON,JLEV,NSV_LIMA_NI) < 1.E-15)) THEN
      PRS(JLON,JLEV,1) = PRS(JLON,JLEV,1) + PRS(JLON,JLEV,4)
      PTHS(JLON,JLEV) = PTHS(JLON,JLEV) &
                      & - PRS(JLON,JLEV,4) * ZLS(JLON,JLEV) / ZCPH(JLON,JLEV) / PEXNREF(JLON,JLEV)
      PRS(JLON,JLEV,4)  = 0.0
      PSVS(JLON,JLEV,NSV_LIMA_NI) = 0.0
    ENDIF
  ENDDO
ENDDO

DO JSV = 1, KSV
  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      PSVS(JLON,JLEV,JSV) = MAX(0.0, PSVS(JLON,JLEV,JSV))
    ENDDO
  ENDDO
ENDDO
!
!
!*       3.3  STORE THE BUDGET TERMS
!            ----------------------

LL_RRR_BUDGET = (TBUCONF%LBUDGET_RV).OR.(TBUCONF%LBUDGET_RC).OR.(TBUCONF%LBUDGET_RR).OR.(TBUCONF%LBUDGET_RI) &
             & .OR.(TBUCONF%LBUDGET_RS).OR.(TBUCONF%LBUDGET_RG).OR.(TBUCONF%LBUDGET_RH)       

IF (LL_RRR_BUDGET) THEN      
  DO JRR=1,KRR
    DO JLEV = 1, KLEV
      DO JLON = KIDIA, KFDIA
        ZWORK(JLON,JLEV) = PRS(JLON,JLEV,JRR)*PRHODJ(JLON,JLEV)
      ENDDO
    ENDDO
    CALL BUDGET_DDH(ZWORK, JRR+5,'NEGA_BU_RRR',YDDDH,YDLDDH, YDMDDH)
  END DO 
END IF

IF (TBUCONF%LBUDGET_TH) THEN
  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      ZWORK(JLON,JLEV) = PTHS(JLON,JLEV)*PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK,4,'NEGA_BU_RTH',YDDDH, YDLDDH, YDMDDH)
END IF

IF (TBUCONF%LBUDGET_SV) THEN

  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      ZWORK(JLON,JLEV) = PSVS(JLON,JLEV,NSV_LIMA_NC)*PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK,12+NSV_LIMA_NC,'NEGA_BU_RSV',YDDDH, YDLDDH, YDMDDH)

  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      ZWORK(JLON,JLEV) = PSVS(JLON,JLEV,NSV_LIMA_NR)*PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK,12+NSV_LIMA_NR,'NEGA_BU_RSV',YDDDH, YDLDDH, YDMDDH)

  DO JLEV = 1, KLEV
    DO JLON = KIDIA, KFDIA
      ZWORK(JLON,JLEV) = PSVS(JLON,JLEV,NSV_LIMA_NI)*PRHODJ(JLON,JLEV)
    ENDDO
  ENDDO
  CALL BUDGET_DDH (ZWORK,12+NSV_LIMA_NI,'NEGA_BU_RSV',YDDDH, YDLDDH, YDMDDH)

  IF (NMOD_CCN.GE.1) THEN
    DO JL=1, NMOD_CCN
      DO JLEV = 1, KLEV
        DO JLON = KIDIA, KFDIA
          ZWORK(JLON,JLEV) = PSVS(JLON,JLEV,NSV_LIMA_CCN_FREE+JL-1)*PRHODJ(JLON,JLEV)
        ENDDO
      ENDDO
      CALL BUDGET_DDH(ZWORK,12+NSV_LIMA_CCN_FREE+JL-1,'NEGA_BU_RSV',YDDDH,YDLDDH, YDMDDH) 
    END DO
  END IF

  IF (NMOD_IFN.GE.1) THEN
    DO JL=1, NMOD_IFN
      DO JLEV = 1, KLEV
        DO JLON = KIDIA, KFDIA
          ZWORK(JLON,JLEV) = PSVS(JLON,JLEV,NSV_LIMA_IFN_FREE+JL-1)*PRHODJ(JLON,JLEV)
        ENDDO
      ENDDO
      CALL BUDGET_DDH(ZWORK,12+NSV_LIMA_IFN_FREE+JL-1,'NEGA_BU_RSV',YDDDH,YDLDDH, YDMDDH) 
    END DO
  END IF
END IF

DO JRR=1, NBUDGET_SV1+NSV_LIMA-1
   YLBUDGET(JRR)%NBUDGET=JRR
   YLBUDGET(JRR)%YDDDH=>YDDDH
   YLBUDGET(JRR)%YDLDDH=>YDLDDH
   YLBUDGET(JRR)%YDMDDH=>YDMDDH
   YLBUDGET_PTR(JRR)%PTR=>YLBUDGET(JRR)
ENDDO
!
!
!-------------------------------------------------------------------------------
!

!*       9.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
!*                   Compute the explicit microphysical sources
!
!
!
CALL LIMA (LIMAP=PHYEX%PARAM_LIMA, LIMAW=PHYEX%PARAM_LIMA_WARM, LIMAC=PHYEX%PARAM_LIMA_COLD, &
           LIMAM=PHYEX%PARAM_LIMA_MIXED, TNSV=PHYEX%TNSV, &
           D=YLDIMPHYEX, CST=PHYEX%CST, NEBN=PHYEX%NEBN, ICED=PHYEX%RAIN_ICE_DESCRN, ICEP=PHYEX%RAIN_ICE_PARAMN,   &
           ELECD=PHYEX%ELEC_DESCR, ELECP=PHYEX%ELEC_PARAM, &
           BUCONF=TBUCONF, TBUDGETS=YLBUDGET_PTR,HACTCCN=HACTCCN, KBUDGETS=SIZE(YLBUDGET_PTR), KRR=KRR, &
           PTSTEP=2*PTSTEP, OELEC=PHYEX%MISC%OELEC,                 &
           PRHODREF=PRHODREF, PEXNREF=PEXNREF, PDZZ=PDZZ, PTHVREFZIKB=ZTHVREFZIKB,       &
           PRHODJ=PRHODJ, PPABST=PPABSM,                                 &
           KCARB=NCARB, KSOA=NSOA, KSP=NSP, ODUST=LDUST, OSALT=LSALT, OORILAM=LORILAM,  &
           ODTHRAD=.TRUE., PDTHRAD=PDTHRAD, PTHT=PTHT, PRT=PRT, PSVT=PSVT, PCIT=PCIT, PW_NU=PW_NU,  &
           PAERO=PSVT, PSOLORG=ZSOLORG, PMI=ZMI, &
           PTHS=PTHS, PRS=PRS, PSVS=PSVS,                                &
           PINPRC=ZINPRC, PINDEP=ZINDEP, PINPRR=PINPRR, PINPRI=ZINPRI, PINPRS=PINPRS, PINPRG=PINPRG, PINPRH=PINPRH, &
           PEVAP3D=PEVAP, PCLDFR=PCLDFR, PICEFR=PICEFR, PPRCFR=PPRCFR, PFPR=PFPR, &
           PHLC_HCF=PHLC_HCF, PHLC_HRC=PHLC_HRC, PHLI_HCF=PHLI_HCF, PHLI_HRI=PHLI_HRI )
!add ZINPRC in PINPRR
DO JLON = KIDIA, KFDIA
  PINPRR(JLON) = PINPRR(JLON) + ZINPRC(JLON)
ENDDO
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_LIMA',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_LIMA
