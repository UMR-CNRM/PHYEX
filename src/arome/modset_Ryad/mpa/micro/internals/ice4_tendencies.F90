SUBROUTINE ICE4_TENDENCIES(KPROMA,KSIZE, KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL, &
                          &KRR, LDSOFT, PCOMPUTE, &
                          &OWARM, HSUBG_RC_RR_ACCR, HSUBG_RR_EVAP, HSUBG_AUCV_RC, HSUBG_PR_PDF, &
                          &PEXN, PRHODREF, PLVFACT, PLSFACT, K1, K2, K3,&
                          &PPRES, PCF, PSIGMA_RC, &
                          &PCIT, &
                          &PT, PVART, &
                          &PRRT3D, &
                          &PRVHENI_MR, PRRHONG_MR, PRIMLTC_MR, PRSRIMCG_MR, &
                          &PRCHONI, PRVDEPS, PRIAGGS, PRIAUTS, PRVDEPG, &
                          &PRCAUTR, PRCACCR, PRREVAV, &
                          &PRCRIMSS, PRCRIMSG, PRSRIMCG, PRRACCSS, PRRACCSG, PRSACCRG, PRSMLTG, PRCMLTSR, &
                          &PRICFRRG, PRRCFRIG, PRICFRR, PRCWETG, PRIWETG, PRRWETG, PRSWETG, &
                          &PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, PRWETGH, PRWETGH_MR, PRGMLTR, &
                          &PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, &
                          &PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, PRDRYHG, PRHMLTR, &
                          &PRCBERI, &
                          &PRS_TEND, PRG_TEND, PRH_TEND, PSSI, &
                          &PA, PB, &
                          &PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC)
!!
!!**  PURPOSE
!!    -------
!!      Computes the tendencies
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!!     R. El Khatib 24-Aug-2021 Optimizations
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_RAIN_ICE_PARAM
USE MODD_RAIN_ICE_DESCR
USE MODD_PARAM_ICE, ONLY : CSNOWRIMING, LFEEDBACKT
USE MODD_BUDGET,    ONLY : LBU_ENABLE
USE MODD_FIELDS_ADDRESS, ONLY : & ! common fields adress
      & ITH,     & ! Potential temperature
      & IRV,     & ! Water vapor
      & IRC,     & ! Cloud water
      & IRR,     & ! Rain water
      & IRI,     & ! Pristine ice
      & IRS,     & ! Snow/aggregate
      & IRG,     & ! Graupel
      & IRH        ! Hail
USE MODI_ICE4_NUCLEATION
USE MODI_ICE4_RSRIMCG_OLD
USE MODI_ICE4_COMPUTE_PDF
USE MODI_ICE4_RAINFR_VERT
USE MODI_ICE4_SLOW
USE MODI_ICE4_WARM
USE MODI_ICE4_FAST_RS
USE MODI_ICE4_FAST_RG
USE MODI_ICE4_FAST_RH
USE MODI_ICE4_FAST_RI

USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KPROMA,KSIZE, KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL
INTEGER,                      INTENT(IN)    :: KRR
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PCOMPUTE
LOGICAL,                      INTENT(IN)    :: OWARM
CHARACTER*80,                 INTENT(IN)    :: HSUBG_RC_RR_ACCR
CHARACTER*80,                 INTENT(IN)    :: HSUBG_RR_EVAP
CHARACTER(len=4),             INTENT(IN)    :: HSUBG_AUCV_RC
CHARACTER*80,                 INTENT(IN)    :: HSUBG_PR_PDF ! pdf for subgrid precipitation
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PEXN
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PRHODREF
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PLSFACT

INTEGER, DIMENSION(KPROMA),    INTENT(IN)    :: K1
INTEGER, DIMENSION(KPROMA),    INTENT(IN)    :: K2
INTEGER, DIMENSION(KPROMA),    INTENT(IN)    :: K3
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PPRES
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PCF
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PSIGMA_RC
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PCIT
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PT
REAL, DIMENSION(KPROMA,0:KRR), INTENT(IN)    :: PVART







REAL, DIMENSION(KIT,KJT,KKT),  INTENT(IN)    :: PRRT3D
REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PRVHENI_MR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRRHONG_MR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRIMLTC_MR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRSRIMCG_MR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCHONI
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRVDEPS
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRIAGGS
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRIAUTS
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRVDEPG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCAUTR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCACCR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRREVAV
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCRIMSS
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCRIMSG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRSRIMCG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRRACCSS
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRRACCSG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRSACCRG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRSMLTG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCMLTSR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRICFRRG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRRCFRIG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRICFRR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCWETG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRIWETG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRRWETG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRSWETG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCDRYG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRIDRYG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRRDRYG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRSDRYG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRWETGH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRWETGH_MR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRGMLTR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCWETH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRIWETH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRSWETH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRGWETH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRRWETH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCDRYH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRIDRYH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRSDRYH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRRDRYH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRGDRYH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRDRYHG
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRHMLTR
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRCBERI
REAL, DIMENSION(KPROMA, 8),    INTENT(INOUT) :: PRS_TEND
REAL, DIMENSION(KPROMA, 8),    INTENT(INOUT) :: PRG_TEND
REAL, DIMENSION(KPROMA, 10),   INTENT(INOUT) :: PRH_TEND
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PSSI
REAL, DIMENSION(KPROMA,0:KRR), INTENT(OUT)   :: PA







REAL, DIMENSION(KPROMA,0:KRR), INTENT(OUT)   :: PB







REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PHLC_HCF
REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PHLC_LCF
REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PHLC_HRC
REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PHLC_LRC
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KPROMA,0:KRR) :: ZVART
REAL, DIMENSION(KPROMA) :: ZT, ZRAINFR, &
                        & ZKA, ZDV, ZAI, ZCJ, &
                        & ZLBDAR, ZLBDAS, ZLBDAG, ZLBDAH, ZLBDAR_RF, &
                        & ZRGSI, ZRGSI_MR
REAL, DIMENSION(KIT,KJT,KKT) :: ZRRT3D, ZRAINFR3D
INTEGER :: JL, JV
REAL, DIMENSION(KPROMA) :: ZWETG ! 1. if graupel growths in wet mode, 0. otherwise
REAL(KIND=JPRB) :: ZZW
LOGICAL :: LLRFR
LOGICAL, DIMENSION(KPROMA) :: LLCOMPUTE

REAL(KIND=JPRB) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ICE4_TENDENCIES', 0, ZHOOK_HANDLE)

LLCOMPUTE(1:KSIZE)=PCOMPUTE(1:KSIZE)==1.

!
ZT(:)=PT(:)
DO JV=0,KRR
  ZVART(:,JV)=PVART(:,JV)
! better fill PA/PB arrays always inside this subroutine. REK.
  PA(:,JV)=0.
ENDDO

IF(LDSOFT) THEN

  PRVHENI_MR(:)=0.
  PRRHONG_MR(:)=0.
  PRIMLTC_MR(:)=0.
  PRSRIMCG_MR(:)=0.
  DO JV=0,KRR
    PB(:,JV)=0.
  ENDDO

ELSE
  !
  !*       2.     COMPUTES THE SLOW COLD PROCESS SOURCES
  !               --------------------------------------

  CALL ICE4_NUCLEATION(KSIZE, LLCOMPUTE, &
                       ZVART(:,ITH), PPRES, PRHODREF, PEXN, PLSFACT, ZT, &
                       ZVART(:,IRV), &
                       PCIT, PRVHENI_MR)

    DO JL=1, KSIZE
      ZVART(JL,ITH)=ZVART(JL,ITH) + PRVHENI_MR(JL)*PLSFACT(JL)
      ZT(JL) = ZVART(JL,ITH) * PEXN(JL)
      ZVART(JL,IRV)=ZVART(JL,IRV) - PRVHENI_MR(JL)
      ZVART(JL,IRI)=ZVART(JL,IRI) + PRVHENI_MR(JL)
    ENDDO
    !
    !*       3.3     compute the spontaneous freezing source: RRHONG
    !
    DO JL=1, KSIZE
      IF (ZT(JL)<XTT-35.0 .AND. ZVART(JL,IRR)>XRTMIN(3) .AND. LLCOMPUTE(JL)) THEN
        PRRHONG_MR(JL)=ZVART(JL,IRR)
      ELSE
        PRRHONG_MR(JL) = 0.
      ENDIF
    ENDDO
    IF (LFEEDBACKT) THEN
      DO JL=1, KSIZE
        !Limitation due to -35 crossing of temperature
        PRRHONG_MR(JL)=MIN(PRRHONG_MR(JL), MAX(0., ((XTT-35.)/PEXN(JL)-ZVART(JL,ITH))/(PLSFACT(JL)-PLVFACT(JL))))
      ENDDO
    ENDIF
    DO JL=1, KSIZE
      ZVART(JL,ITH) = ZVART(JL,ITH) + PRRHONG_MR(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RRHONG))
      ZT  (JL) = ZVART(JL,ITH) * PEXN(JL)
      ZVART(JL,IRR) = ZVART(JL,IRR) - PRRHONG_MR(JL)
      ZVART(JL,IRG) = ZVART(JL,IRG) + PRRHONG_MR(JL)
    ENDDO
    !
    !*       7.1    cloud ice melting
    !
    DO JL=1, KSIZE
      IF (ZVART(JL,IRI)>0. .AND. ZT(JL)>XTT .AND. LLCOMPUTE(JL)) THEN
        PRIMLTC_MR(JL)=ZVART(JL,IRI)
      ELSE
        PRIMLTC_MR(JL)=0.
      ENDIF
    ENDDO
    IF(LFEEDBACKT) THEN
      DO JL=1, KSIZE
        !Limitation due to 0 crossing of temperature
        PRIMLTC_MR(JL)=MIN(PRIMLTC_MR(JL), MAX(0., (ZVART(JL,ITH)-XTT/PEXN(JL)) / (PLSFACT(JL)-PLVFACT(JL))))
      ENDDO
    ENDIF
    DO JL=1, KSIZE
      ZVART(JL,ITH) = ZVART(JL,ITH) - PRIMLTC_MR(JL)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(-RIMLTC))
      ZT  (JL) = ZVART(JL,ITH) * PEXN(JL)
      ZVART(JL,IRC) = ZVART(JL,IRC) + PRIMLTC_MR(JL)
      ZVART(JL,IRI) = ZVART(JL,IRI) - PRIMLTC_MR(JL)
    ENDDO

    DO JL=1, KSIZE
      PB(JL,ITH) = PRVHENI_MR(JL)*PLSFACT(JL)
      PB(JL,ITH) = PB(JL,ITH) + PRRHONG_MR(JL)*(PLSFACT(JL)-PLVFACT(JL))
      PB(JL,ITH) = PB(JL,ITH) - PRIMLTC_MR(JL)*(PLSFACT(JL)-PLVFACT(JL))

      PB(JL,IRV) = - PRVHENI_MR(JL)

      PB(JL,IRC) = + PRIMLTC_MR(JL)

      PB(JL,IRR) = - PRRHONG_MR(JL)

      PB(JL,IRI) = + PRVHENI_MR(JL)
      PB(JL,IRI) = PB(JL,IRI) - PRIMLTC_MR(JL)

      PB(JL,IRS) = 0.

      PB(JL,IRG) = + PRRHONG_MR(JL)
    ENDDO
    IF (KRR==7) PB(:,IRH)=0.
  !
  !        5.1.6  riming-conversion of the large sized aggregates into graupel (old parametrisation)
  !
  IF(CSNOWRIMING=='OLD ') THEN
    ZLBDAS(1:KSIZE)=0.
    WHERE(ZVART(1:KSIZE,IRS)>0.)
      ZLBDAS(1:KSIZE)  = MIN(XLBDAS_MAX, XLBS*(PRHODREF(1:KSIZE)*MAX(ZVART(1:KSIZE,IRS), XRTMIN(5)))**XLBEXS)
    END WHERE
    CALL ICE4_RSRIMCG_OLD(KSIZE, LDSOFT, LLCOMPUTE, &
                         &PRHODREF, &
                         &ZLBDAS, &
                         &ZT, ZVART(:,IRC), ZVART(:,IRS), &
                         &PRSRIMCG_MR, PB(:,IRS), PB(:,IRG))
    DO JL=1, KSIZE
      ZVART(JL,IRS) = ZVART(JL,IRS) - PRSRIMCG_MR(JL)
      ZVART(JL,IRG) = ZVART(JL,IRG) + PRSRIMCG_MR(JL)
    ENDDO
  ELSE
    PRSRIMCG_MR(:) = 0.
  ENDIF

!
!* Derived fields
!
  DO JL=1, KSIZE
    ZZW = EXP(XALPI-XBETAI/ZT(JL)-XGAMI*ALOG(ZT(JL)))
    PSSI(JL) = ZVART(JL,IRV)*( PPRES(JL)-ZZW ) / ( XEPSILO * ZZW ) - 1.0
                                                      ! Supersaturation over ice
    ZKA(JL) = 2.38E-2 + 0.0071E-2*(ZT(JL)-XTT) ! k_a
    ZDV(JL) = 0.211E-4*(ZT(JL)/XTT)**1.94 * (XP00/PPRES(JL)) ! D_v
    ZAI(JL) = (XLSTT+(XCPV-XCI)*(ZT(JL)-XTT))**2 / (ZKA(JL)*XRV*ZT(JL)**2) &
                                 + ( XRV*ZT(JL) ) / (ZDV(JL)*ZZW)
    ZCJ(JL) = XSCFAC*PRHODREF(JL)**0.3 / SQRT(1.718E-5+0.0049E-5*(ZT(JL)-XTT))
  ENDDO

ENDIF ! LDSOFT
!
!Cloud water split between high and low content part is done here
CALL ICE4_COMPUTE_PDF(KSIZE, HSUBG_AUCV_RC, HSUBG_PR_PDF,&
                      PRHODREF, ZVART(:,IRC), PCF, PSIGMA_RC,&
                      PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC, ZRAINFR)
! ZRAINFR is computed above, then overwritten below :-( .REK.


j ai des doutes sur le codage de Ryad, normalement zrainfr sert à remplir prainfr qui remplit à son tour zrainfr

LLRFR=HSUBG_RC_RR_ACCR=='PRFR'
IF (LLRFR) THEN
  !Diagnostic of precipitation fraction
  ZRRT3D(:,:,:)=PRRT3D(:,:,:)
  DO JL=1,KSIZE
    ZRRT3D (K1(JL),K2(JL),K3(JL))=PVART(JL,IRR)-PRRHONG_MR(JL)
  ENDDO
  CALL ICE4_RAINFR_VERT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL, ZRAINFR3D(:,:,:), ZRRT3D(:,:,:))
  DO JL=1,KSIZE
    ZRAINFR(JL)=ZRAINFR3D(K1(JL), K2(JL), K3(JL))
  ENDDO
ELSE
  ZRAINFR(:)=1.
ENDIF
!
!*  compute the slope parameters
!
DO JL=1, KSIZE
  !ZLBDAR will be used when we consider rain diluted over the grid box
  IF (ZVART(JL,IRR)>0.) THEN
    ZLBDAR(JL)  = XLBR*( PRHODREF(JL)*MAX(ZVART(JL,IRR), XRTMIN(3)))**XLBEXR
  ELSE
    ZLBDAR(JL)=0.
  ENDIF
  !ZLBDAR_RF is used when we consider rain concentrated in its fraction
  IF (LLRFR) THEN
    IF (ZVART(JL,IRR)>0. .AND. ZRAINFR(JL)>0.) THEN
      ZLBDAR_RF(JL) = XLBR*( PRHODREF(JL)*MAX(ZVART(JL,IRR)/ZRAINFR(JL), XRTMIN(3)))**XLBEXR
    ELSE
      ZLBDAR_RF(JL)=0.
    ENDIF
  ELSE
    ZLBDAR_RF(JL) = ZLBDAR(JL)
  ENDIF
  IF (ZVART(JL,IRS)>0.) THEN
    ZLBDAS(JL)  = MIN(XLBDAS_MAX, XLBS*(PRHODREF(JL)*MAX(ZVART(JL,IRS), XRTMIN(5)))**XLBEXS)
  ELSE
    ZLBDAS(JL)=0.
  ENDIF
  IF (ZVART(JL,IRG)>0.) THEN
    ZLBDAG(JL)  = XLBG*(PRHODREF(JL)*MAX(ZVART(JL,IRG), XRTMIN(6)))**XLBEXG
  ELSE
    ZLBDAG(JL)=0.
  ENDIF
  IF (KRR==7) THEN
    IF (PVART(JL,IRH)>0.) THEN
      ZLBDAH(JL) = XLBH*(PRHODREF(JL)*MAX(PVART(JL,IRH), XRTMIN(7)))**XLBEXH
    ELSE
      ZLBDAH(JL)=0.
    ENDIF
  ENDIF
ENDDO
!
!
CALL ICE4_SLOW(KSIZE, LDSOFT, PCOMPUTE, PRHODREF, ZT, &
              &PSSI, PLVFACT, PLSFACT, &
              &ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), &
              &ZLBDAS, ZLBDAG, &
              &ZAI, ZCJ, &
              &PRCHONI, PRVDEPS, PRIAGGS, PRIAUTS, PRVDEPG, &
              &PA(:,ITH), PA(:,IRV), PA(:,IRC), PA(:,IRI), PA(:,IRS), PA(:,IRG))
!
!-------------------------------------------------------------------------------
!
!
!*       3.     COMPUTES THE SLOW WARM PROCESS SOURCES
!               --------------------------------------
!
!
IF(OWARM) THEN    !  Check if the formation of the raindrops by the slow
                  !  warm processes is allowed
  CALL ICE4_WARM(KSIZE, LDSOFT, PCOMPUTE, HSUBG_RC_RR_ACCR, HSUBG_RR_EVAP, &
                &PRHODREF, PLVFACT, ZT, PPRES, ZVART(:,ITH),&
                &ZLBDAR, ZLBDAR_RF, ZKA, ZDV, ZCJ, &
                &PHLC_LCF, PHLC_HCF, PHLC_LRC, PHLC_HRC, &
                &PCF, ZRAINFR, &
                &ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), &
                &PRCAUTR, PRCACCR, PRREVAV, &
                &PA(:,ITH), PA(:,IRV), PA(:,IRC), PA(:,IRR))
ELSE
  PRCAUTR(:)=0.
  PRCACCR(:)=0.
  PRREVAV(:)=0.
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_s
!               ----------------------------------------------
!
CALL ICE4_FAST_RS(KPROMA,KSIZE, LDSOFT, PCOMPUTE, &
                 &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                 &ZDV, ZKA, ZCJ, &
                 &ZLBDAR, ZLBDAS, &
                 &ZT, ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), ZVART(:,IRS), &
                 &PRIAGGS, &
                 &PRCRIMSS, PRCRIMSG, PRSRIMCG, &
                 &PRRACCSS, PRRACCSG, PRSACCRG, PRSMLTG, &
                 &PRCMLTSR, &
                 &PRS_TEND, &
                 &PA(:,ITH), PA(:,IRC), PA(:,IRR), PA(:,IRS), PA(:,IRG))
!
!-------------------------------------------------------------------------------
!
!
!*       5.        COMPUTES THE FAST COLD PROCESS SOURCES FOR r_g
!                  ------------------------------------------------------
!
! Move these inits inside. REK
DO JL=1, KSIZE
  ZRGSI(JL) = PRVDEPG(JL) + PRSMLTG(JL) + PRRACCSG(JL) + &
           & PRSACCRG(JL) + PRCRIMSG(JL) + PRSRIMCG(JL)
  ZRGSI_MR(JL) = PRRHONG_MR(JL) + PRSRIMCG_MR(JL)
ENDDO
CALL ICE4_FAST_RG(KPROMA,KSIZE, LDSOFT, PCOMPUTE, KRR, &
                 &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                 &ZDV, ZKA, ZCJ, PCIT, &
                 &ZLBDAR, ZLBDAS, ZLBDAG, &
                 &ZT, ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), &
                 &ZRGSI, ZRGSI_MR(:), &
                 &ZWETG, &
                 &PRICFRRG, PRRCFRIG, PRICFRR, PRCWETG, PRIWETG, PRRWETG, PRSWETG, &
                 &PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, PRWETGH, PRWETGH_MR, PRGMLTR, &
                 &PRG_TEND, &
                 &PA(:,ITH), PA(:,IRC), PA(:,IRR), PA(:,IRI), PA(:,IRS), PA(:,IRG), PA(:,IRH), PB(:,IRG), PB(:,IRH))
!
!-------------------------------------------------------------------------------
!
!
!*       6.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_h
!               ----------------------------------------------
!
IF (KRR==7) THEN
  CALL ICE4_FAST_RH(KPROMA,KSIZE, LDSOFT, PCOMPUTE, ZWETG, &
                   &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                   &ZDV, ZKA, ZCJ, &
                   &ZLBDAS, ZLBDAG, ZLBDAR, ZLBDAH, &
                   &ZT,  ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), PVART(:,IRH), &
                   &PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, &
                   &PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, PRDRYHG, PRHMLTR, &
                   &PRH_TEND, &
                   &PA(:,ITH), PA(:,IRC), PA(:,IRR), PA(:,IRI), PA(:,IRS), PA(:,IRG), PA(:,IRH))
ELSEIF (LBU_ENABLE) THEN
  PRCWETH(:)=0.
  PRIWETH(:)=0.
  PRSWETH(:)=0.
  PRGWETH(:)=0.
  PRRWETH(:)=0.
  PRCDRYH(:)=0.
  PRIDRYH(:)=0.
  PRSDRYH(:)=0.
  PRRDRYH(:)=0.
  PRGDRYH(:)=0.
  PRDRYHG(:)=0.
  PRHMLTR(:)=0.
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       7.     COMPUTES SPECIFIC SOURCES OF THE WARM AND COLD CLOUDY SPECIES
!               -------------------------------------------------------------
!
CALL ICE4_FAST_RI(KSIZE, LDSOFT, PCOMPUTE, &
                 &PRHODREF, PLVFACT, PLSFACT, &
                 &ZAI, ZCJ, PCIT, &
                 &PSSI, &
                 &ZVART(:,IRC), ZVART(:,IRI), &
                 &PRCBERI, PA(:,ITH), PA(:,IRC), PA(:,IRI))
!
IF (LHOOK) CALL DR_HOOK('ICE4_TENDENCIES', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_TENDENCIES
