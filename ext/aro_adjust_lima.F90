!     ######spl
      SUBROUTINE  ARO_ADJUST_LIMA(PHYEX, &
                                  KKA,KKU,KKL,KLON,KLEV,KFDIA,  KRR, KSV, KTCOUNT,  &
                                  OSUBG_COND, OSIGMAS, &
                                  PTSTEP, PSIGQSAT, &
                                  PZZF, PRHODJ, PRHODREF, PEXNREF,&
                                  PPABSM, PTHT, PRT, PSVT, PSIGS, &
                                  PW_NU, PDTHRAD, &
                                  PMFCONV, PRC_MF, PRI_MF, PCF_MF, &
                                  PTHS, PRS,  PSVS, PSRCS, PCLDFR, PICEFR, PPRCFR, &
                                  YDDDH, YDLDDH, YDMDDH, LLIMAINIT )
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     ##########################################################################
!
!!****  * -  compute the  resolved clouds and precipitation
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the  microphysical sources
!!    related to the resolved clouds and precipitation
!!
!!
!!**  METHOD
!!    ------
!!      The main actions of this routine is to call the routines computing the
!!    microphysical sources. Before that:
!!        - it computes the real absolute pressure,
!!        - negative values of the current guess of all mixing ratio are removed.
!!          This is done by a global filling algorithm based on a multiplicative
!!          method (Rood, 1987), in order to conserved the total mass in the
!!          simulation domain.
!!        - Sources are transformed in physical tendencies, by removing the
!!          multiplicative term Rhod*J.
!!        - External points values are filled owing to the use of cyclic
!!          l.b.c., in order to performe computations on the full domain.
!!      After calling to microphysical routines, the physical tendencies are
!!    switched back to prognostic variables.
!!
!!
!!    EXTERNAL
!!    --------
!!      Subroutine FMLOOK: to recover the logical unit number linked to a FMfile
!!      Subroutine SLOW_TERMS: Computes the explicit microphysical sources
!!      Subroutine FAST_TERMS: Performs the saturation adjustment for l
!!      Subroutine RAIN_ICE  : Computes the explicit microphysical sources for i
!!      Subroutine ICE_ADJUST: Performs the saturation adjustment for i+l
!!      MIN_ll,SUM3D_ll : distributed functions equivalent to MIN and SUM
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : contains declarations of parameter variables
!!         JPHEXT       : Horizontal external points number
!!         JPVEXT       : Vertical external points number
!!      Module MODD_CST
!!          XP00               ! Reference pressure
!!          XRD                ! Gaz  constant for dry air
!!          XCPD               ! Cpd (dry air)
!!
!!    REFERENCE
!!    ---------
!!
!!      Documentation AROME
!!
!!    AUTHOR
!!    ------
!!    S.Malardel and Y.Seity
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    10/03/03
!!      T. Kovacic  11-05-05, Call to budgets for NEGA1_
!!      S. Riette ice for EDKF
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!USE MODD_CONF
USE MODD_PHYEX, ONLY: PHYEX_t
USE MODD_PARAMETERS
USE MODD_BUDGET, ONLY: TBUDGETDATA, NBUDGET_SV1, TBUCONF
!
USE MODD_PARAM_LIMA
USE MODD_NSV
!
USE MODI_LIMA_ADJUST_SPLIT
USE MODE_SET_CONC_LIMA
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!USE MODE_BUDGET_PHY, ONLY: BUDGET_DDH
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
TYPE(PHYEX_t),            INTENT(IN)   :: PHYEX
INTEGER,                  INTENT(IN)   :: KKA    !near ground array index
INTEGER,                  INTENT(IN)   :: KKU    !uppest atmosphere array index
INTEGER,                  INTENT(IN)   :: KKL    !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KFDIA    !
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KSV      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
LOGICAL,                  INTENT(IN)   :: OSIGMAS  ! Switch for Sigma_s:
                                        ! use values computed in CONDENSATION
                                        ! or that from turbulence scheme
REAL,                     INTENT(IN)   :: PTSTEP   ! Time step
REAL,                     INTENT(IN)   :: PSIGQSAT ! coeff applied to qsat variance contribution
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PZZF     ! Height (z)
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRHODREF
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PPABSM  ! abs. pressure at time t-dt
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) :: PRT     ! Moist variables at time t
REAL, DIMENSION(KLON,1,KLEV,KSV), INTENT(INOUT) :: PSVT     ! Moist variables at time t
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PW_NU   ! w for CCN activation
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PDTHRAD ! rad theta tendency for CCN activation
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PMFCONV ! convective mass flux
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRC_MF, PRI_MF, PCF_MF
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(INOUT) :: PTHS  ! Theta source
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(KLON,1,KLEV,KSV), INTENT(INOUT) :: PSVS   ! Moist  variable sources
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(OUT)   :: PSRCS ! Second-order flux
                                                 ! s'rc'/2Sigma_s2 at time t+1
                                                 ! multiplied by Lambda_3
REAL, DIMENSION(KLON,1,KLEV), INTENT(INOUT)   :: PCLDFR! Cloud fraction
REAL, DIMENSION(KLON,1,KLEV), INTENT(INOUT)   :: PICEFR! Cloud fraction
REAL, DIMENSION(KLON,1,KLEV), INTENT(INOUT)   :: PPRCFR! Cloud fraction
!
!
TYPE(TYP_DDH), INTENT(INOUT), TARGET :: YDDDH
TYPE(TLDDH), INTENT(IN), TARGET :: YDLDDH
TYPE(TMDDH), INTENT(IN), TARGET :: YDMDDH
!
LOGICAL,                  INTENT(IN)    :: LLIMAINIT
!
!*       0.2   Declarations of local variables :

!
INTEGER :: JRR           ! Loop index for the moist and scalar variables
!
REAL, DIMENSION(SIZE(PZZF,1),SIZE(PZZF,2),SIZE(PZZF,3)):: ZT,ZLV,ZLS,ZCPH
REAL, DIMENSION(SIZE(PZZF,1),SIZE(PZZF,2),SIZE(PZZF,3)):: ZCOR
                                    ! for the correction of negative rv
REAL, DIMENSION(SIZE(PZZF,1),SIZE(PZZF,2),SIZE(PZZF,3)):: ZZZ
                                    ! model layer height
REAL  :: ZMASSTOT                   ! total mass  for one water category
                                    ! including the negative values
REAL  :: ZMASSPOS                   ! total mass  for one water category
                                    ! after removing the negative values
REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR
!
TYPE(TBUDGETDATA), DIMENSION(NBUDGET_SV1+NSV_LIMA-1) :: YLBUDGET
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
!
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ARO_ADJUST_LIMA',0,ZHOOK_HANDLE)

CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, 0, KFDIA)

!
!*       2.     TRANSFORMATION INTO PHYSICAL TENDENCIES
!               ---------------------------------------
!
!
! personal comment:  tranfering these variables to the
!                    microphysical routines would save
!                    computing time
!
ZT(:,:,:)= PTHT(:,:,:)*PEXNREF(:,:,:)
ZLV(:,:,:)=PHYEX%CST%XLVTT +(PHYEX%CST%XCPV-PHYEX%CST%XCL) *(ZT(:,:,:)-PHYEX%CST%XTT)
ZLS(:,:,:)=PHYEX%CST%XLSTT +(PHYEX%CST%XCPV-PHYEX%CST%XCI) *(ZT(:,:,:)-PHYEX%CST%XTT)
ZCPH(:,:,:)=PHYEX%CST%XCPD +PHYEX%CST%XCPV*2.*PTSTEP*PRS(:,:,:,1)

!set concentration for LIMA
PRS = PRS * 2.*PTSTEP
PSVS = PSVS * 2.*PTSTEP
IF (LLIMAINIT) THEN
   CALL SET_CONC_LIMA (1,'ICE3',PRHODREF,PRT,PSVT)
   CALL SET_CONC_LIMA (1,'ICE3',PRHODREF,PRS,PSVS)
ELSE
   CALL SET_CONC_LIMA (1,'ICE3',PRHODREF,PRT,PSVT, .TRUE.)
   CALL SET_CONC_LIMA (1,'ICE3',PRHODREF,PRS,PSVS, .TRUE.)
END IF
PRS = PRS / (2.*PTSTEP)
PSVS = PSVS / (2.*PTSTEP)

!print *, "aro_adjust_lima 2"
!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for precipitating species (Rood 87)
!
DO JRR = 3,KRR
   SELECT CASE (JRR)
   CASE(3,5,6,7) ! rain, snow, graupel and hail

      IF ( MINVAL( PRS(:,:,:,JRR)) < 0.0 ) THEN
! For AROME, we cannot use MAX_ll so that according to JPP's advises
!  we only correct negative values but not the total mass
! compute the total water mass computation
!
!          ZMASSTOT = MAX( 0. , SUM( PRS(:,:,:,JRR) ))
!
! remove the negative values
!
         PRS(:,:,:,JRR) = MAX( 0., PRS(:,:,:,JRR) )
!
! compute the new total mass
!
!          ZMASSPOS = MAX(1.E-60,SUM( PRS(:,:,:,JRR) ))
!
! correct again in such a way to conserve the total mass
!
!          ZRATIO = ZMASSTOT / ZMASSPOS
!          PRS(:,:,:,JRR) = PRS(:,:,:,JRR) * ZRATIO

      END IF
   END SELECT
END DO
!
!*       3.2    Correct negative values
!
! Correction where rc<0
IF (NMOM_C.GE.1) THEN
!        WHERE (PRS(:,:,:,2) < 0. .OR. PSVS(:,:,:,NSV_LIMA_NC) < 0.)
   WHERE (PRS(:,:,:,2) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,2)  = 0.0
      PSVS(:,:,:,NSV_LIMA_NC) = 0.0
   END WHERE
END IF
! Correction where rr<0
IF (NMOM_R.GE.1) THEN
!        WHERE (PRS(:,:,:,3) < 0. .OR. PSVS(:,:,:,NSV_LIMA_NR) < 0.)
   WHERE (PRS(:,:,:,3) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,3)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,3) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,3)  = 0.0
      PSVS(:,:,:,NSV_LIMA_NR) = 0.0
   END WHERE
END IF
! Correction of IFN concentrations where ri<0 or Ni<0
!     IF (LCOLD_LIMA) THEN
!        DO JMOD = 1, NMOD_IFN 
!           WHERE (PRS(:,:,:,4) < 0. .OR. PSVS(:,:,:,NSV_LIMA_NI) < 0.) ! ri or Ni < 0.
!              PSVS(:,:,:,NSV_LIMA_IFN_FREE+JMOD-1) =               &
!                   PSVS(:,:,:,NSV_LIMA_IFN_FREE+JMOD-1) + &
!                   PSVS(:,:,:,NSV_LIMA_IFN_NUCL+JMOD-1)     ! N_IF =N_IF+N_IN
!              PSVS(:,:,:,NSV_LIMA_IFN_NUCL+JMOD-1) = 0.0             ! N_IN =0.
!           END WHERE
!        ENDDO
!     END IF
! Correction where ri<0
IF (NMOM_I.GE.1) THEN
!        WHERE (PRS(:,:,:,4) < 0. .OR. PSVS(:,:,:,NSV_LIMA_NI) < 0.)
   WHERE (PRS(:,:,:,4) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,4)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,4) * ZLS(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,4)  = 0.0
      PSVS(:,:,:,NSV_LIMA_NI) = 0.0
   END WHERE
END IF
!
PSVS(:,:,:,:) = MAX( 0.0,PSVS(:,:,:,:) )
!
!
!*       3.3  STORE THE BUDGET TERMS
!            ----------------------
!
!IF (LBUDGET_RV) CALL BUDGET (PRS(:,:,:,1) * PRHODJ(:,:,:), 6,'NEGA_BU_RRV')
!IF (LBUDGET_RC) CALL BUDGET (PRS(:,:,:,2) * PRHODJ(:,:,:), 7,'NEGA_BU_RRC')
!IF (LBUDGET_RR) CALL BUDGET (PRS(:,:,:,3) * PRHODJ(:,:,:), 8,'NEGA_BU_RRR')
!IF (LBUDGET_RI) CALL BUDGET (PRS(:,:,:,4) * PRHODJ(:,:,:) ,9,'NEGA_BU_RRI')
!IF (LBUDGET_RS) CALL BUDGET (PRS(:,:,:,5) * PRHODJ(:,:,:),10,'NEGA_BU_RRS')
!IF (LBUDGET_RG) CALL BUDGET (PRS(:,:,:,6) * PRHODJ(:,:,:),11,'NEGA_BU_RRG')
!IF (LBUDGET_RH) CALL BUDGET (PRS(:,:,:,7) * PRHODJ(:,:,:),12,'NEGA_BU_RRH')
!IF (LBUDGET_TH) CALL BUDGET (PTHS(:,:,:)  * PRHODJ(:,:,:), 4,'NEGA_BU_RTH')

DO JRR = 1, NBUDGET_SV1+NSV_LIMA-1
   YLBUDGET(JRR)%NBUDGET=JRR
   YLBUDGET(JRR)%YDDDH=>YDDDH
   YLBUDGET(JRR)%YDLDDH=>YDLDDH
   YLBUDGET(JRR)%YDMDDH=>YDMDDH
ENDDO
!
!-------------------------------------------------------------------------------
!

!*       9.     MIXED-PHASE MICROPHYSICAL SCHEME (WITH 3 ICE SPECIES)
!               -----------------------------------------------------
!
!
!*       9.2    Perform the saturation adjustment over cloud ice and cloud water
!
    ZZZ =  PZZF

    CALL LIMA_ADJUST_SPLIT(D=YLDIMPHYEX, CST=PHYEX%CST, BUCONF=TBUCONF, TBUDGETS=YLBUDGET, KBUDGETS=SIZE(YLBUDGET), &
         KRR=KRR, KMI=PHYEX%MISC%KMI, HCONDENS=PHYEX%NEBN%CCONDENS, HLAMBDA3=PHYEX%NEBN%CLAMBDA3, &
         OSUBG_COND=OSUBG_COND, OSIGMAS=OSIGMAS, PTSTEP=2*PTSTEP, PSIGQSAT=PSIGQSAT, &
         PRHODREF=PRHODREF, PRHODJ=PRHODJ, PEXNREF=PEXNREF, PSIGS=PSIGS, PMFCONV=PMFCONV, &
         PPABST=PPABSM, PPABSTT=PPABSM, PZZ=ZZZ, PDTHRAD=PDTHRAD, PW_NU=PW_NU, &
         PRT=PRT, PRS=PRS, PSVT=PSVT, PSVS=PSVS, &
         PTHS=PTHS, PSRCS=PSRCS, PCLDFR=PCLDFR, PICEFR=PICEFR, PRC_MF=PRC_MF, PRI_MF=PRI_MF, PCF_MF=PCF_MF )
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_ADJUST_LIMA',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_ADJUST_LIMA
