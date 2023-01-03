!     ######spl
      SUBROUTINE  ARO_RAIN_ICE(KPROMA,KKA,KKU,KKL,KLON,KLEV, KFDIA, KRR, KTCOUNT, KSPLITR,&
                                  OSUBG_COND, CSUBG_AUCV_RC, CSUBG_AUCV_RI,OSEDIC, CSEDIM, CMICRO, &
                                  PTSTEP, PDZZ, PRHODJ, PRHODREF, PEXNREF,&
                                  PPABSM, PHLC_HRC, PHLC_HCF, PHLI_HRI, PHLI_HCF, PTHT, PRT, PSIGS,PCLDFR, &
                                  PTHS, PRS, PEVAP,  &
                                  PCIT, OWARM, PSEA, PTOWN,   &
                                  PICLDFR, PWCLDFR, PSSIO, PSSIU, PIFR,  &
                                  OCND2, LKOGAN, LMODICEDEP,&
                                  PINPRR,PINPRS,PINPRG,PINPRH,PFPR,     &
                                  YDDDH, YDLDDH, YDMDDH, &
                                  YSPP_ICENU,YSPP_KGN_ACON,YSPP_KGN_SBGR)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
!!      15/05/05 T. Kovacic, budgets for negative correction
!!      29/09/08 Y. Seity, add PEVAP for chemistry
!!      2012-02 Y. Seity,  add possibility to run with reversed vertical levels
!!      2013-11, D. Degrauwe: Export upper-air precipitation fluxes PFPR.
!!      2013-11 S. Riette, subgrid precipitation
!!      2014-11 S. Riette, ICE3/ICE4 modified, old versions under OLD3/OLD4
!!      2014-11 S. Riette, ICE3/ICE4 modified, old versions under OLD3/OLD4
!!      2020-12 U. Andrae : Introduce SPP for HARMONIE-AROME
!!      2018-02 K.I: Ivarsson: More inputs to OCND2-option for saving computing time.
!!     R. El Khatib 24-Aug-2021 Specific cache-blocking factor for microphysics
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF
USE MODD_CST, ONLY: CST
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM
USE MODD_PARAM_ICE,      ONLY: PARAM_ICE
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
!
USE MODD_BUDGET, ONLY: TBUDGETDATA, NBUDGET_RH, TBUCONF
USE MODE_BUDGET, ONLY: BUDGET_DDH
USE MODE_FILL_DIMPHYEX, ONLY: FILL_DIMPHYEX
!
USE MODI_RAIN_ICE
!
USE SPP_MOD_TYPE, ONLY : TSPP_CONFIG_TYPE, APPLY_SPP
!
USE MODI_RAIN_ICE_OLD
!
USE DDH_MIX , ONLY : TYP_DDH
USE YOMLDDH , ONLY : TLDDH
USE YOMMDDH , ONLY : TMDDH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!

!
INTEGER,                  INTENT(IN)   :: KPROMA  !internal cache-blocking factor for microphysic loop in rain_ice
INTEGER,                  INTENT(IN)   :: KKA  !near ground array index
INTEGER,                  INTENT(IN)   :: KKU  !uppest atmosphere array index
INTEGER,                  INTENT(IN)   :: KKL  !vert. levels type 1=MNH -1=ARO
INTEGER,                  INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER,                  INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER,                  INTENT(IN)   :: KFDIA
INTEGER,                  INTENT(IN)   :: KRR      ! Number of moist variables
INTEGER,                  INTENT(IN)   :: KTCOUNT  ! Temporal loop counter
INTEGER,                  INTENT(IN)   :: KSPLITR  ! Number of small time step
                                       ! integrations for  rain sedimendation
LOGICAL,                  INTENT(IN)   :: OSUBG_COND ! Switch for Subgrid Cond.
CHARACTER (LEN=4),        INTENT(IN)   :: CSUBG_AUCV_RC
                                        ! type of subgrid rc->rr autoconvertion scheme
CHARACTER (LEN=80),       INTENT(IN)   :: CSUBG_AUCV_RI
                                        ! type of subgrid ri->rs autoconvertion scheme
LOGICAL,                  INTENT(IN)   :: OSEDIC ! Switch for cloud sedim.
CHARACTER (LEN=4),        INTENT(IN)   :: CSEDIM  ! Sedimentation scheme
                                                  ! (STAT or EULE)
CHARACTER (LEN=4),        INTENT(IN)   :: CMICRO  ! Microphysics scheme
REAL,                     INTENT(IN)   :: PTSTEP   ! Time step
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PDZZ     ! Height (z)
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRHODJ  !Dry density * Jacobian
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PRHODREF! Reference dry air density
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PEXNREF ! Reference Exner function
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PPABSM  ! abs. pressure at time t-dt
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PHLC_HRC
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PHLC_HCF
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PHLI_HRI
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PHLI_HCF
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PTHT    ! Theta at time t
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT):: PRT   ! Moist variables at time t
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PSIGS   ! Sigma_s at time t
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PCLDFR  ! Cloud fraction
! input from aro_adjust / condensation with OCND2, dummy if OCND2 = F
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PICLDFR ! ice cloud fraction
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PWCLDFR ! water or mixed-phase cloud fraction
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PSSIO   ! Super-saturation with respect to ice in the  
                                                        ! supersaturated fraction
REAL, DIMENSION(KLON,1,KLEV),   INTENT(IN)   :: PSSIU   ! Sub-saturation with respect to ice in the  
                                                        ! subsaturated fraction 
REAL, DIMENSION(KLON,1,KLEV),   INTENT(INOUT):: PIFR    ! Ratio cloud ice moist part to dry part 
!REAL, DIMENSION (KLON,1),       INTENT(IN)  :: PPBL    ! PBL top above ground (m)
! input from aro_adjust / condensation with OCND2 END.
!
!
REAL, DIMENSION(KLON,1,KLEV),   INTENT(INOUT) :: PTHS    ! Theta source
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) :: PRS   ! Moist  variable sources
REAL, DIMENSION(KLON,1,KLEV), INTENT(INOUT) :: PEVAP ! Rain evap profile
!
!

REAL, DIMENSION(KLON,1,KLEV), INTENT(INOUT)   :: PCIT  ! Pristine ice number
                                                 ! concentration at time t
LOGICAL,                  INTENT(IN)    :: OWARM ! Control of the rain formation
                                                 !  by slow warm microphysical
                                                 !         processes
LOGICAL,                  INTENT(IN)    :: OCND2 ! Logical switch to separate liquid and ice
LOGICAL,                  INTENT(IN)    :: LKOGAN! Logical switch for using Kogan autoconversion of liquid
LOGICAL,                  INTENT(IN)    :: LMODICEDEP ! Logical switch for alternative dep/evap of ice
REAL, DIMENSION(KLON,1), INTENT(IN)        :: PSEA  ! Land sea mask
REAL, DIMENSION(KLON,1), INTENT(IN)        :: PTOWN  ! Town mask
REAL, DIMENSION(KLON,1), INTENT(OUT)       :: PINPRR! Rain instant precip
REAL, DIMENSION(KLON,1), INTENT(OUT)       :: PINPRS! Snow instant precip
REAL, DIMENSION(KLON,1), INTENT(OUT)       :: PINPRG! Graupel instant precip
REAL, DIMENSION(KLON,1), INTENT(OUT)       :: PINPRH! Hail instant precip
REAL, DIMENSION(KLON,1,KLEV,KRR), INTENT(INOUT) :: PFPR ! upper-air precip
!
TYPE(TYP_DDH), INTENT(INOUT), TARGET :: YDDDH
TYPE(TLDDH), INTENT(IN), TARGET :: YDLDDH
TYPE(TMDDH), INTENT(IN), TARGET :: YDMDDH
!
TYPE(TSPP_CONFIG_TYPE), INTENT(INOUT) :: YSPP_ICENU,YSPP_KGN_ACON,YSPP_KGN_SBGR 
!
!
!*       0.2   Declarations of local variables :
INTEGER :: JRR           ! Loop index for the moist and scalar variables
!
!
!
REAL, DIMENSION(KLON,1,KLEV):: ZT,ZLV,ZLS,ZCPH
REAL, DIMENSION(KLON,1,KLEV):: ZCOR
REAL, DIMENSION(KLON,1):: ZINDEP     ! surf cloud deposition (already contained in sedimentation)
REAL, DIMENSION(KLON,1,KLEV):: ZRAINFR
REAL, DIMENSION(KLON,1) :: ZICENU, ZKGN_ACON, ZKGN_SBGR
REAL, DIMENSION(KLON,1):: ZINPRC    ! surf cloud sedimentation
                                    ! for the correction of negative rv
REAL  :: ZMASSTOT                   ! total mass  for one water category
                                    ! including the negative values
REAL  :: ZMASSPOS                   ! total mass  for one water category
                                    ! after removing the negative values
REAL  :: ZRATIO                     ! ZMASSTOT / ZMASSCOR

LOGICAL, DIMENSION(KLON, 1, KLEV) :: LLMICRO !mask to limit computation

TYPE(TBUDGETDATA), DIMENSION(NBUDGET_RH) :: YLBUDGET !NBUDGET_RH is the one with the highest number
TYPE(DIMPHYEX_t) :: YLDIMPHYEX
!
INTEGER :: IPROMA, ISIZE, IGPBLKS ! cache-blocking management
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "abor1.intfb.h"
!------------------------------------------------------------------------------
!
!*       1.     PRELIMINARY COMPUTATIONS
!               ------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_RAIN_ICE',0,ZHOOK_HANDLE)

!Dimensions
CALL FILL_DIMPHYEX(YLDIMPHYEX, KLON, 1, KLEV, 0, KFDIA)

ZINPRC=0.
PINPRH=0.

!Mask to limit computation
IF ( KRR == 7 ) THEN
  IF (CMICRO /= 'ICE4' .AND. CMICRO /= 'OLD4') THEN
    CALL ABOR1('ARO_RAIN_ICE : KRR==7 NOT COMPATIBLE WITH CMICRO /= ICE4 OR OLD4')
  ENDIF
  LLMICRO(:,:,:)=                          &
                PRT(:,:,:,2)>RAIN_ICE_DESCR%XRTMIN(2) .OR. &
                PRT(:,:,:,3)>RAIN_ICE_DESCR%XRTMIN(3) .OR. &
                PRT(:,:,:,4)>RAIN_ICE_DESCR%XRTMIN(4) .OR. &
                PRT(:,:,:,5)>RAIN_ICE_DESCR%XRTMIN(5) .OR. &
                PRT(:,:,:,6)>RAIN_ICE_DESCR%XRTMIN(6) .OR. &
                PRT(:,:,:,7)>RAIN_ICE_DESCR%XRTMIN(7)
ELSE IF( KRR == 6 ) THEN
  LLMICRO(:,:,:)=                          &
                PRT(:,:,:,2)>RAIN_ICE_DESCR%XRTMIN(2) .OR. &
                PRT(:,:,:,3)>RAIN_ICE_DESCR%XRTMIN(3) .OR. &
                PRT(:,:,:,4)>RAIN_ICE_DESCR%XRTMIN(4) .OR. &
                PRT(:,:,:,5)>RAIN_ICE_DESCR%XRTMIN(5) .OR. &
                PRT(:,:,:,6)>RAIN_ICE_DESCR%XRTMIN(6)
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
ZT(:,:,:)= PTHT(:,:,:)*PEXNREF(:,:,:)
ZLV(:,:,:)=CST%XLVTT +(CST%XCPV-CST%XCL) *(ZT(:,:,:)-CST%XTT)
ZLS(:,:,:)=CST%XLSTT +(CST%XCPV-CST%XCI) *(ZT(:,:,:)-CST%XTT)
ZCPH(:,:,:)=CST%XCPD +CST%XCPV*2.*PTSTEP*PRS(:,:,:,1)
!
!
!*       3.     REMOVE NEGATIVE VALUES
!               ----------------------
!
!*       3.1    Non local correction for precipitating species (Rood 87)
!
IF (CMICRO == 'KESS' .OR. CMICRO == 'ICE3' .OR. CMICRO == 'ICE2' &
    .OR.  CMICRO == 'C2R2' .OR. CMICRO == 'C3R5'.OR. CMICRO == 'ICE4' .OR. &
    CMICRO == 'OLD3' .OR. CMICRO == 'OLD4') THEN

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
END IF
!
!*       3.2    Adjustement for liquid and solid cloud
!
SELECT CASE ( CMICRO )
!
!
  CASE('ICE2','ICE3','ICE4', 'OLD3', 'OLD4')
    WHERE (PRS(:,:,:,4) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,4)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,4) * ZLS(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,4) = 0.
    END WHERE
!
!   cloud
    WHERE (PRS(:,:,:,2) < 0.)
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,2) = 0.
    END WHERE
!
! if rc or ri are positive, we can correct negative rv
!   cloud
    WHERE ((PRS(:,:,:,1) <0.) .AND. (PRS(:,:,:,2)> 0.) )
      PRS(:,:,:,1) = PRS(:,:,:,1) + PRS(:,:,:,2)
      PTHS(:,:,:) = PTHS(:,:,:) - PRS(:,:,:,2) * ZLV(:,:,:) /  &
           ZCPH(:,:,:) / PEXNREF(:,:,:)
      PRS(:,:,:,2) = 0.
    END WHERE
!   ice
    IF(KRR > 3) THEN
      WHERE ((PRS(:,:,:,1) < 0.).AND.(PRS(:,:,:,4) > 0.))
        ZCOR(:,:,:)=MIN(-PRS(:,:,:,1),PRS(:,:,:,4))
        PRS(:,:,:,1) = PRS(:,:,:,1) + ZCOR(:,:,:)
        PTHS(:,:,:) = PTHS(:,:,:) - ZCOR(:,:,:) * ZLS(:,:,:) /  &
             ZCPH(:,:,:) / PEXNREF(:,:,:)
        PRS(:,:,:,4) = PRS(:,:,:,4) -ZCOR(:,:,:)
      END WHERE
    END IF
!
END SELECT
!
!
!*       3.3  STORE THE BUDGET TERMS
!            ----------------------
IF (TBUCONF%LBUDGET_RV) CALL BUDGET_DDH (PRS(:,:,:,1) * PRHODJ(:,:,:), 6,'NEGA_BU_RRV',YDDDH, YDLDDH, YDMDDH)
IF (TBUCONF%LBUDGET_RC) CALL BUDGET_DDH (PRS(:,:,:,2) * PRHODJ(:,:,:), 7,'NEGA_BU_RRC',YDDDH, YDLDDH, YDMDDH)
IF (TBUCONF%LBUDGET_RR) CALL BUDGET_DDH (PRS(:,:,:,3) * PRHODJ(:,:,:), 8,'NEGA_BU_RRR',YDDDH, YDLDDH, YDMDDH)
IF (TBUCONF%LBUDGET_RI) CALL BUDGET_DDH (PRS(:,:,:,4) * PRHODJ(:,:,:) ,9,'NEGA_BU_RRI',YDDDH, YDLDDH, YDMDDH)
IF (TBUCONF%LBUDGET_RS) CALL BUDGET_DDH (PRS(:,:,:,5) * PRHODJ(:,:,:),10,'NEGA_BU_RRS',YDDDH, YDLDDH, YDMDDH)
IF (TBUCONF%LBUDGET_RG) CALL BUDGET_DDH (PRS(:,:,:,6) * PRHODJ(:,:,:),11,'NEGA_BU_RRG',YDDDH, YDLDDH, YDMDDH)
IF (TBUCONF%LBUDGET_RH .AND. KRR==7) CALL BUDGET_DDH (PRS(:,:,:,7) * PRHODJ(:,:,:),12,'NEGA_BU_RRH',YDDDH, YDLDDH, YDMDDH)
IF (TBUCONF%LBUDGET_TH) CALL BUDGET_DDH (PTHS(:,:,:)  * PRHODJ(:,:,:), 4,'NEGA_BU_RTH',YDDDH, YDLDDH, YDMDDH)

DO JRR=1, NBUDGET_RH
  YLBUDGET(JRR)%NBUDGET=JRR
  YLBUDGET(JRR)%YDDDH=>YDDDH
  YLBUDGET(JRR)%YDLDDH=>YDLDDH
  YLBUDGET(JRR)%YDMDDH=>YDMDDH
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
ISIZE=COUNT(LLMICRO)
IF (KPROMA > 0 .AND. ISIZE > 0) THEN
  ! Cache-blocking is active
  ! number of chunks :
  IGPBLKS = (ISIZE-1)/MIN(KPROMA,ISIZE)+1
  ! Adjust IPROMA to limit the number of small chunks
  IPROMA=(ISIZE-1)/IGPBLKS+1
ELSE
  IPROMA=ISIZE ! no cache-blocking
ENDIF
IF (CMICRO=='ICE4') THEN
    CALL RAIN_ICE(  YLDIMPHYEX, CST, PARAM_ICE, RAIN_ICE_PARAM, &
                 &  RAIN_ICE_DESCR, TBUCONF, &
                 &  IPROMA, ISIZE, &
                 &  OSEDIC=OSEDIC, OCND2=OCND2, HSEDIM=CSEDIM, &
                 &  HSUBG_AUCV_RC=CSUBG_AUCV_RC, HSUBG_AUCV_RI=CSUBG_AUCV_RI,&
                 &  OWARM=OWARM, &
                 &  PTSTEP=2*PTSTEP, &
                 &  KRR=KRR, ODMICRO=LLMICRO, PEXN=PEXNREF,            &
                 &  PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF, PEXNREF=PEXNREF,&
                 &  PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
                 &  PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, &
                 &  PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF, &
                 &  PTHT=PTHT,PRVT= PRT(:,:,:,1),PRCT= PRT(:,:,:,2), &
                 &  PRRT=PRT(:,:,:,3), &
                 &  PRIT=PRT(:,:,:,4), PRST=PRT(:,:,:,5), &
                 &  PRGT=PRT(:,:,:,6),       &
                 &  PTHS=PTHS, PRVS=PRS(:,:,:,1),PRCS=PRS(:,:,:,2),&
                 &  PRRS=PRS(:,:,:,3),&
                 &  PRIS=PRS(:,:,:,4),PRSS= PRS(:,:,:,5),PRGS= PRS(:,:,:,6),&
                 &  PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
                 &  PINPRS=PINPRS, PINPRG=PINPRG, PINDEP=ZINDEP, PRAINFR=ZRAINFR, &
                 &  PSIGS=PSIGS, &
                 &  TBUDGETS=YLBUDGET, KBUDGETS=SIZE(YLBUDGET), &
                 &  PSEA=PSEA, PTOWN=PTOWN, &
                 &  PRHT=PRT(:,:,:,7), PRHS=PRS(:,:,:,7), PINPRH=PINPRH, PFPR=PFPR)
ELSEIF (CMICRO=='ICE3') THEN
    CALL RAIN_ICE(  YLDIMPHYEX, CST, PARAM_ICE, RAIN_ICE_PARAM, &
                 &  RAIN_ICE_DESCR, TBUCONF, &
                 &  IPROMA, ISIZE, &
                 &  OSEDIC=OSEDIC, OCND2=OCND2, HSEDIM=CSEDIM, &
                 &  HSUBG_AUCV_RC=CSUBG_AUCV_RC, HSUBG_AUCV_RI=CSUBG_AUCV_RI,&
                 &  OWARM=OWARM, &
                 &  PTSTEP=2*PTSTEP, &
                 &  KRR=KRR, ODMICRO=LLMICRO, PEXN=PEXNREF,            &
                 &  PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF,PEXNREF=PEXNREF,&
                 &  PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
                 &  PHLC_HRC=PHLC_HRC, PHLC_HCF=PHLC_HCF, &
                 &  PHLI_HRI=PHLI_HRI, PHLI_HCF=PHLI_HCF, &
                 &  PTHT=PTHT,PRVT=PRT(:,:,:,1),PRCT=PRT(:,:,:,2), &
                 &  PRRT=PRT(:,:,:,3), &
                 &  PRIT=PRT(:,:,:,4), PRST=PRT(:,:,:,5), &
                 &  PRGT=PRT(:,:,:,6),       &
                 &  PTHS=PTHS, PRVS=PRS(:,:,:,1),PRCS=PRS(:,:,:,2),&
                 &  PRRS=PRS(:,:,:,3),&
                 &  PRIS=PRS(:,:,:,4),PRSS= PRS(:,:,:,5),PRGS= PRS(:,:,:,6),&
                 &  PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
                 &  PINPRS=PINPRS, PINPRG=PINPRG, PINDEP=ZINDEP, PRAINFR=ZRAINFR, &
                 &  PSIGS=PSIGS, &
                 &  TBUDGETS=YLBUDGET, KBUDGETS=SIZE(YLBUDGET), &
                 &  PSEA=PSEA, PTOWN=PTOWN, PFPR=PFPR)
ELSEIF (CMICRO=='OLD4') THEN
    IF (YSPP_ICENU%LPERT) THEN
     CALL APPLY_SPP(YSPP_ICENU,KLON,1,KLON,RAIN_ICE_PARAM%XFRMIN(9),ZICENU)
    ELSE
     ZICENU(:,:) = RAIN_ICE_PARAM%XFRMIN(9)
    ENDIF
   
    IF (YSPP_KGN_ACON%LPERT) THEN
     CALL APPLY_SPP(YSPP_KGN_ACON,KLON,1,KLON,RAIN_ICE_PARAM%XFRMIN(10),ZKGN_ACON)
    ELSE
     ZKGN_ACON(:,:) = RAIN_ICE_PARAM%XFRMIN(10)
    ENDIF
   
    IF (YSPP_KGN_SBGR%LPERT) THEN
     CALL APPLY_SPP(YSPP_KGN_SBGR,KLON,1,KLON,RAIN_ICE_PARAM%XFRMIN(11),ZKGN_SBGR)
    ELSE
     ZKGN_SBGR(:,:) = RAIN_ICE_PARAM%XFRMIN(11)
    ENDIF
    CALL RAIN_ICE_OLD( OSEDIC=OSEDIC, OCND2=OCND2, LKOGAN=LKOGAN, LMODICEDEP=LMODICEDEP, &
                 &  HSEDIM=CSEDIM, HSUBG_AUCV_RC=CSUBG_AUCV_RC, &
                 &  OWARM=OWARM,KKA=KKA,KKU=KKU,KKL=KKL,KSPLITR=KSPLITR, &
                 &  PTSTEP=2*PTSTEP, KRR=KRR,                              &
                 &  PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF, PEXNREF=PEXNREF,&
                 &  PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
                 &  PICLDFR=PICLDFR, PWCLDFR=PWCLDFR, &
                 &  PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR, &
                 &  PTHT=PTHT,PRVT= PRT(:,:,:,1),PRCT= PRT(:,:,:,2), &
                 &  PRRT=PRT(:,:,:,3), &
                 &  PRIT=PRT(:,:,:,4), PRST=PRT(:,:,:,5), &
                 &  PRGT=PRT(:,:,:,6),       &
                 &  PTHS=PTHS, PRVS=PRS(:,:,:,1),PRCS=PRS(:,:,:,2),&
                 &  PRRS=PRS(:,:,:,3),&
                 &  PRIS=PRS(:,:,:,4),PRSS= PRS(:,:,:,5),PRGS= PRS(:,:,:,6),&
                 &  PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
                 &  PINPRS=PINPRS, PINPRG=PINPRG, &
                 &  PSIGS=PSIGS, PSEA=PSEA, PTOWN=PTOWN, &
                 &  YDDDH=YDDDH,YDLDDH=YDLDDH,YDMDDH=YDMDDH, &
                 &  PRHT=PRT(:,:,:,7),&
                 &  PRHS=PRS(:,:,:,7), PINPRH=PINPRH, PFPR=PFPR, &
                 &  PICENU=ZICENU, &
                 &  PKGN_ACON=ZKGN_ACON, &
                 &  PKGN_SBGR=ZKGN_SBGR)
ELSE
    IF (YSPP_ICENU%LPERT) THEN
     CALL APPLY_SPP(YSPP_ICENU,KLON,1,KLON,RAIN_ICE_PARAM%XFRMIN(9),ZICENU)
    ELSE
     ZICENU(:,:) = RAIN_ICE_PARAM%XFRMIN(9)
    ENDIF
   
    IF (YSPP_KGN_ACON%LPERT) THEN
     CALL APPLY_SPP(YSPP_KGN_ACON,KLON,1,KLON,RAIN_ICE_PARAM%XFRMIN(10),ZKGN_ACON)
    ELSE
     ZKGN_ACON(:,:) = RAIN_ICE_PARAM%XFRMIN(10)
    ENDIF
   
    IF (YSPP_KGN_SBGR%LPERT) THEN
     CALL APPLY_SPP(YSPP_KGN_SBGR,KLON,1,KLON,RAIN_ICE_PARAM%XFRMIN(11),ZKGN_SBGR)
    ELSE
     ZKGN_SBGR(:,:) = RAIN_ICE_PARAM%XFRMIN(11)
    ENDIF
    CALL RAIN_ICE_OLD( OSEDIC=OSEDIC, OCND2=OCND2, LKOGAN=LKOGAN, LMODICEDEP=LMODICEDEP, &
                 &  HSEDIM=CSEDIM, HSUBG_AUCV_RC=CSUBG_AUCV_RC, &
                 &  OWARM=OWARM,KKA=KKA,KKU=KKU,KKL=KKL,KSPLITR=KSPLITR, &
                 &  PTSTEP=2*PTSTEP, KRR=KRR,                              &
                 &  PDZZ=PDZZ, PRHODJ=PRHODJ, PRHODREF=PRHODREF, PEXNREF=PEXNREF,&
                 &  PPABST=PPABSM, PCIT=PCIT, PCLDFR=PCLDFR,  &
                 &  PICLDFR=PICLDFR, PWCLDFR=PWCLDFR, &
                 &  PSSIO=PSSIO, PSSIU=PSSIU, PIFR=PIFR, &
                 &  PTHT=PTHT,PRVT= PRT(:,:,:,1),PRCT= PRT(:,:,:,2), &
                 &  PRRT=PRT(:,:,:,3), &
                 &  PRIT=PRT(:,:,:,4), PRST=PRT(:,:,:,5), &
                 &  PRGT=PRT(:,:,:,6),       &
                 &  PTHS=PTHS, PRVS=PRS(:,:,:,1),PRCS=PRS(:,:,:,2),&
                 &  PRRS=PRS(:,:,:,3),&
                 &  PRIS=PRS(:,:,:,4),PRSS= PRS(:,:,:,5),PRGS= PRS(:,:,:,6),&
                 &  PINPRC=ZINPRC,PINPRR=PINPRR,PEVAP3D=PEVAP,&
                 &  PINPRS=PINPRS, PINPRG=PINPRG, &
                 &  PSIGS=PSIGS, PSEA=PSEA, PTOWN=PTOWN, &
                 &  YDDDH=YDDDH,YDLDDH=YDLDDH,YDMDDH=YDMDDH, &
                 &  PFPR=PFPR, &
                 &  PICENU=ZICENU, &
                 &  PKGN_ACON=ZKGN_ACON, &
                 &  PKGN_SBGR=ZKGN_SBGR)
ENDIF
!add ZINPRC in PINPRR
PINPRR=PINPRR+ZINPRC
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ARO_RAIN_ICE',1,ZHOOK_HANDLE)
END SUBROUTINE ARO_RAIN_ICE
