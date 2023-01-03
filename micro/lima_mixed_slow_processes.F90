!      #####################################
       MODULE MODI_LIMA_MIXED_SLOW_PROCESSES
!      #####################################
!
INTERFACE
      SUBROUTINE LIMA_MIXED_SLOW_PROCESSES(ZRHODREF, ZZT, ZSSI, PTSTEP,  &
                                           ZLSFACT, ZLVFACT, ZAI, ZCJ,   &
                                           ZRGT, ZCIT,                   &
                                           ZRVS, ZRCS, ZRIS, ZRGS, ZTHS, &
                                           ZCCS, ZCIS, ZIFS, ZINS,       &
                                           ZLBDAI, ZLBDAG,               &
                                           ZRHODJ, GMICRO, PRHODJ, KMI,  &
                                           PTHS, PRVS, PRCS, PRIS, PRGS, &
                                           PCCS, PCIS, &
                            YDDDH, YDLDDH, YDMDDH                    )
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: ZZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: ZSSI      ! Supersaturation over ice
REAL,                 INTENT(IN)    :: PTSTEP    ! Time step          
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZAI       ! Thermodynamical function
REAL, DIMENSION(:),   INTENT(IN)    :: ZCJ       ! for the ventilation coefficient
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRGT      ! Graupel/hail m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZCIT      ! Pristine ice conc. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRVS      ! Water vapor m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRCS      ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRIS      ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRGS      ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZTHS      ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCCS      ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCIS      ! Pristine ice conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: ZIFS      ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: ZINS      ! Nucleated Ice nuclei conc. source 
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAI  ! Slope parameter of the ice crystal distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAG  ! Slope parameter of the graupel distr.
!
! used for budget storage
REAL,    DIMENSION(:),     INTENT(IN) :: ZRHODJ
LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: GMICRO 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
INTEGER,                   INTENT(IN) :: KMI 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PTHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRVS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRIS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRGS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCIS
!
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
END SUBROUTINE LIMA_MIXED_SLOW_PROCESSES
END INTERFACE
END MODULE MODI_LIMA_MIXED_SLOW_PROCESSES
!
!     #######################################################################
      SUBROUTINE LIMA_MIXED_SLOW_PROCESSES(ZRHODREF, ZZT, ZSSI, PTSTEP,  &
                                           ZLSFACT, ZLVFACT, ZAI, ZCJ,   &
                                           ZRGT, ZCIT,                   &
                                           ZRVS, ZRCS, ZRIS, ZRGS, ZTHS, &
                                           ZCCS, ZCIS, ZIFS, ZINS,       &
                                           ZLBDAI, ZLBDAG,               &
                                           ZRHODJ, GMICRO, PRHODJ, KMI,  &
                                           PTHS, PRVS, PRCS, PRIS, PRGS, &
                                           PCCS, PCIS, &
                            YDDDH, YDLDDH, YDMDDH                    )
!     #######################################################################
!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the mixed-phase 
!!    slow processes : 
!!
!!      Deposition of water vapor on graupeln
!!      Cloud ice Melting
!!      Bergeron-Findeisen effect
!!
!!**  METHOD
!!    ------
!!
!!    REFERENCE
!!    ---------
!!
!!      Most of the parameterizations come from the ICE3 scheme, described in
!!    the MESO-NH scientific documentation.
!!
!!      Cohard, J.-M. and J.-P. Pinty, 2000: A comprehensive two-moment warm 
!!      microphysical bulk scheme. 
!!        Part I: Description and tests
!!        Part II: 2D experiments with a non-hydrostatic model
!!      Accepted for publication in Quart. J. Roy. Meteor. Soc. 
!!
!!    AUTHOR
!!    ------
!!      J.-M. Cohard     * Laboratoire d'Aerologie*
!!      J.-P. Pinty      * Laboratoire d'Aerologie*
!!      S.    Berthet    * Laboratoire d'Aerologie*
!!      B.    Vié        * Laboratoire d'Aerologie*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             ??/??/13 
!!      C. Barthe  * LACy *  jan. 2014   add budgets
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,              ONLY : XTT, XALPI, XBETAI, XGAMI,          &
                                       XALPW, XBETAW, XGAMW
USE MODD_PARAM_LIMA,       ONLY : XRTMIN, XCTMIN, NMOD_IFN
USE MODD_PARAM_LIMA_COLD,  ONLY : XDI, X0DEPI, X2DEPI, XSCFAC  
USE MODD_PARAM_LIMA_MIXED, ONLY : XLBG, XLBEXG, XLBDAG_MAX,           &
                                  X0DEPG, XEX0DEPG, X1DEPG, XEX1DEPG 
!
USE MODD_NSV
USE MODD_BUDGET
USE MODE_BUDGET, ONLY: BUDGET_DDH
!
USE DDH_MIX, ONLY  : TYP_DDH
USE YOMLDDH, ONLY  : TLDDH
USE YOMMDDH, ONLY  : TMDDH
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRHODREF  ! RHO Dry REFerence
REAL, DIMENSION(:),   INTENT(IN)    :: ZZT       ! Temperature
REAL, DIMENSION(:),   INTENT(IN)    :: ZSSI      ! Supersaturation over ice
REAL,                 INTENT(IN)    :: PTSTEP    ! Time step          
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLSFACT   ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZLVFACT   ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(:),   INTENT(IN)    :: ZAI       ! Thermodynamical function
REAL, DIMENSION(:),   INTENT(IN)    :: ZCJ       ! for the ventilation coefficient
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZRGT      ! Graupel/hail m.r. at t
REAL, DIMENSION(:),   INTENT(IN)    :: ZCIT      ! Pristine ice conc. at t
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRVS      ! Water vapor m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRCS      ! Cloud water m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRIS      ! Pristine ice m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZRGS      ! Graupel/hail m.r. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZTHS      ! Theta source
!
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCCS      ! Cloud water conc. source
REAL, DIMENSION(:),   INTENT(INOUT) :: ZCIS      ! Pristine ice conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: ZIFS      ! Free Ice nuclei conc. source
REAL, DIMENSION(:,:), INTENT(INOUT) :: ZINS      ! Nucleated Ice nuclei conc. source 
!
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAI  ! Slope parameter of the ice crystal distr.
REAL, DIMENSION(:),   INTENT(IN)    :: ZLBDAG  ! Slope parameter of the graupel distr.
!
! used for budget storage
REAL,    DIMENSION(:),     INTENT(IN) :: ZRHODJ
LOGICAL, DIMENSION(:,:,:), INTENT(IN) :: GMICRO 
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRHODJ
INTEGER,                   INTENT(IN) :: KMI
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PTHS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRVS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRIS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PRGS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCCS
REAL,    DIMENSION(:,:,:), INTENT(IN) :: PCIS
!
TYPE(TYP_DDH), INTENT(INOUT) :: YDDDH
TYPE(TLDDH), INTENT(IN) :: YDLDDH
TYPE(TMDDH), INTENT(IN) :: YDMDDH
!
!*       0.2   Declarations of local variables :
!
REAL, DIMENSION(SIZE(ZZT)) :: ZZW, ZMASK    ! Work vectors
!
INTEGER :: JMOD_IFN
!
!-------------------------------------------------------------------------------
!
!*       1    Deposition of water vapor on r_g: RVDEPG
!        ---------------------------------------------
!
!
   ZZW(:) = 0.0
   WHERE ( (ZRGT(:)>XRTMIN(6)) .AND. (ZRGS(:)>XRTMIN(6)/PTSTEP) )
!Correction BVIE RHODREF
!      ZZW(:) = ( ZSSI(:)/(ZRHODREF(:)*ZAI(:)) ) *                               &
      ZZW(:) = ( ZSSI(:)/(ZAI(:)) ) *                               &
               ( X0DEPG*ZLBDAG(:)**XEX0DEPG + X1DEPG*ZCJ(:)*ZLBDAG(:)**XEX1DEPG )
      ZZW(:) =         MIN( ZRVS(:),ZZW(:)      )*(0.5+SIGN(0.5,ZZW(:))) &
                     - MIN( ZRGS(:),ABS(ZZW(:)) )*(0.5-SIGN(0.5,ZZW(:)))
      ZRGS(:) = ZRGS(:) + ZZW(:)
      ZRVS(:) = ZRVS(:) - ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:)*ZLSFACT(:)
   END WHERE
!
! Budget storage
   IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
     IF (LBUDGET_TH) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS)*PRHODJ(:,:,:),&
                                                                4,'DEPG_BU_RTH',YDDDH, YDLDDH, YDMDDH)
     IF (LBUDGET_RV) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZRVS(:),MASK=GMICRO(:,:,:),FIELD=PRVS)*PRHODJ(:,:,:),&
                                                                6,'DEPG_BU_RRV',YDDDH, YDLDDH, YDMDDH)
     IF (LBUDGET_RG) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZRGS(:),MASK=GMICRO(:,:,:),FIELD=PRGS)*PRHODJ(:,:,:), &
                                                               11,'DEPG_BU_RRG',YDDDH, YDLDDH, YDMDDH)
   END IF
!
!
!*       2    cloud ice Melting: RIMLTC and CIMLTC
!        -----------------------------------------
!
!
   ZMASK(:) = 1.0
   WHERE( (ZRIS(:)>XRTMIN(4)/PTSTEP) .AND. (ZZT(:)>XTT) )
      ZRCS(:) = ZRCS(:) + ZRIS(:)
      ZTHS(:) = ZTHS(:) - ZRIS(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(-RIMLTC))
      ZRIS(:) = 0.0
!
      ZCCS(:) = ZCCS(:) + ZCIS(:)
      ZCIS(:) = 0.0
      ZMASK(:)= 0.0
   END WHERE
   DO JMOD_IFN = 1,NMOD_IFN
! Correction BVIE aerosols not released but in droplets
!      ZIFS(:,JMOD_IFN) = ZIFS(:,JMOD_IFN) + ZINS(:,JMOD_IFN)*(1.-ZMASK(:)) 
      ZINS(:,JMOD_IFN) = ZINS(:,JMOD_IFN) * ZMASK(:)
   ENDDO
!
! Budget storage
   IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
     IF (LBUDGET_TH) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS)*PRHODJ(:,:,:),&
                                                                4,'IMLT_BU_RTH',YDDDH, YDLDDH, YDMDDH)
     IF (LBUDGET_RC) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZRCS(:),MASK=GMICRO(:,:,:),FIELD=PRCS)*PRHODJ(:,:,:), &
                                                                7,'IMLT_BU_RRC',YDDDH, YDLDDH, YDMDDH)
     IF (LBUDGET_RI) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZRIS(:),MASK=GMICRO(:,:,:),FIELD=PRIS)*PRHODJ(:,:,:), &
                                                                9,'IMLT_BU_RRI',YDDDH, YDLDDH, YDMDDH)
     IF (LBUDGET_SV) THEN
       CALL BUDGET_DDH (UNPACK(ZCCS(:),MASK=GMICRO(:,:,:),FIELD=PCCS)*PRHODJ(:,:,:), &
                                                               12+NSV_LIMA_NC,'IMLT_BU_RSV',YDDDH, YDLDDH, YDMDDH)
       CALL BUDGET_DDH (UNPACK(ZCIS(:),MASK=GMICRO(:,:,:),FIELD=PCIS)*PRHODJ(:,:,:), &
                                                               12+NSV_LIMA_NI,'IMLT_BU_RSV',YDDDH, YDLDDH, YDMDDH)
     END IF
   END IF
!
!
!*       3    Bergeron-Findeisen effect: RCBERI
!        --------------------------------------
!
!
   ZZW(:) = 0.0
   WHERE( (ZRCS(:)>XRTMIN(2)/PTSTEP) .AND. (ZRIS(:)>XRTMIN(4)/PTSTEP) .AND. (ZCIT(:)>XCTMIN(4)) )
      ZZW(:) = EXP( (XALPW-XALPI) - (XBETAW-XBETAI)/ZZT(:)          &
                                  - (XGAMW-XGAMI)*ALOG(ZZT(:)) ) -1.0 
                                      ! supersaturation of saturated water over ice
      ZZW(:) = MIN( ZRCS(:),( ZZW(:) / ZAI(:) ) * ZCIT(:) *        &
                    ( X0DEPI/ZLBDAI(:)+X2DEPI*ZCJ(:)*ZCJ(:)/ZLBDAI(:)**(XDI+2.0) ) )
      ZRCS(:) = ZRCS(:) - ZZW(:)
      ZRIS(:) = ZRIS(:) + ZZW(:)
      ZTHS(:) = ZTHS(:) + ZZW(:)*(ZLSFACT(:)-ZLVFACT(:)) ! f(L_f*(RCBERI))
   END WHERE
!
! Budget storage
   IF (NBUMOD==KMI .AND. LBU_ENABLE) THEN
     IF (LBUDGET_TH) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZTHS(:),MASK=GMICRO(:,:,:),FIELD=PTHS)*PRHODJ(:,:,:),&
                                                               4,'BERFI_BU_RTH',YDDDH, YDLDDH, YDMDDH)
     IF (LBUDGET_RC) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZRCS(:),MASK=GMICRO(:,:,:),FIELD=PRCS)*PRHODJ(:,:,:), &
                                                               7,'BERFI_BU_RRC',YDDDH, YDLDDH, YDMDDH)
     IF (LBUDGET_RI) CALL BUDGET_DDH (                                                 &
                   UNPACK(ZRIS(:),MASK=GMICRO(:,:,:),FIELD=PRIS)*PRHODJ(:,:,:), &
                                                               9,'BERFI_BU_RRI',YDDDH, YDLDDH, YDMDDH)
   END IF
!
!------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_MIXED_SLOW_PROCESSES
