!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_TENDENCIES
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_TENDENCIES(CST, PARAMI, ICEP, ICED, BUCONF, KPROMA, KSIZE, &
                          &KRR, ODSOFT, LDCOMPUTE, &
                          &OSAVE_MICRO, OELEC, &
                          &PEXN, PRHODREF, PLVFACT, PLSFACT, &
                          &PPRES, PCF, PSIGMA_RC, &
                          &PCIT, &
                          &PT, PTH, PVART, &
                          &PLATHAM_IAGGS, &
                          &PBU_INST, &
                          &PRS_TEND, PRG_TEND, PRH_TEND, PSSI, &
                          &PA, PB, PATH, PBTH, &
                          &PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC, &
                          &PHLI_HCF, PHLI_LCF, PHLI_HRI, PHLI_LRI, &
                          &PRAINFR)
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
!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!!     R. El Khatib 24-Aug-2021 Optimizations
!!     C. Barthe    06/2023: Add retroaction of electric field on IAGGS
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_BUDGET,         ONLY: TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
!
USE MODD_FIELDS_ADDRESS
USE MODE_ICE4_RRHONG, ONLY: ICE4_RRHONG
USE MODE_ICE4_RIMLTC, ONLY: ICE4_RIMLTC
USE MODE_ICE4_RSRIMCG_OLD, ONLY: ICE4_RSRIMCG_OLD
USE MODE_ICE4_COMPUTE_PDF, ONLY: ICE4_COMPUTE_PDF
USE MODE_ICE4_SLOW, ONLY: ICE4_SLOW
USE MODE_ICE4_WARM, ONLY: ICE4_WARM
USE MODE_ICE4_FAST_RS, ONLY: ICE4_FAST_RS
USE MODE_ICE4_FAST_RG, ONLY: ICE4_FAST_RG
USE MODE_ICE4_FAST_RH, ONLY: ICE4_FAST_RH
USE MODE_ICE4_FAST_RI, ONLY: ICE4_FAST_RI
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(TBUDGETCONF_t),      INTENT(IN)    :: BUCONF
INTEGER,                      INTENT(IN)    :: KPROMA, KSIZE
INTEGER,                      INTENT(IN)    :: KRR
LOGICAL,                      INTENT(IN)    :: ODSOFT
LOGICAL, DIMENSION(KPROMA),   INTENT(IN)    :: LDCOMPUTE
LOGICAL,                      INTENT(IN)    :: OSAVE_MICRO
LOGICAL,                      INTENT(IN)    :: OELEC
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PEXN
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PRHODREF
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PPRES
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PCF
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PSIGMA_RC
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PCIT
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PT
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PTH
REAL, DIMENSION(KPROMA,KRR),   INTENT(IN)    :: PVART
REAL, DIMENSION(MERGE(KPROMA,0,OELEC)), INTENT(IN) :: PLATHAM_IAGGS
REAL, DIMENSION(KPROMA, IBUNUM),INTENT(INOUT):: PBU_INST
REAL, DIMENSION(KPROMA, 8),    INTENT(INOUT) :: PRS_TEND
REAL, DIMENSION(KPROMA, 8),    INTENT(INOUT) :: PRG_TEND
REAL, DIMENSION(KPROMA, 10),   INTENT(INOUT) :: PRH_TEND
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PSSI
REAL, DIMENSION(KPROMA,KRR),   INTENT(OUT)   :: PA
REAL, DIMENSION(KPROMA,KRR),   INTENT(OUT)   :: PB
REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PATH
REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PBTH
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PHLC_LCF
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PHLC_LRC
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PHLI_HCF
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PHLI_LCF
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PHLI_HRI
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PHLI_LRI
REAL, DIMENSION(KPROMA),       INTENT(INOUT) :: PRAINFR   ! Rain fraction
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KPROMA,KRR) :: ZVART
REAL, DIMENSION(KPROMA) :: ZT, ZTH, &
                        & ZKA, ZDV, ZAI, ZCJ, &
                        & ZLBDAR, ZLBDAS, ZLBDAG, ZLBDAH, ZLBDAR_RF, &
                        & ZRGSI, ZRGSI_MR, ZRAINFR
INTEGER :: JL, JV
LOGICAL, DIMENSION(KPROMA) :: LLWETG ! .TRUE. if graupel growths in wet mode
LOGICAL :: LLMASK
REAL :: ZZW
LOGICAL :: LLRFR
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ICE4_TENDENCIES', 0, ZHOOK_HANDLE)


!
!$acc kernels
ZT(:)=PT(:)
ZTH(:)=PTH(:)
PATH(:)=0.
PBTH(:)=0.
DO JV=1,KRR
  ZVART(:,JV)=PVART(:,JV)
  PA(:,JV)=0.
  PB(:,JV)=0.
ENDDO
!$acc end kernels
!
IF(ODSOFT) THEN
!$acc kernels
  PBU_INST(:, IRVHENI_MR)=0.
  PBU_INST(:, IRRHONG_MR)=0.
  PBU_INST(:, IRIMLTC_MR)=0.
  PBU_INST(:, IRSRIMCG_MR)=0.
!$acc end kernels
ELSE
  !
  !*       2.     COMPUTES THE SLOW COLD PROCESS SOURCES
  !               --------------------------------------
  DO JL=1, KSIZE
!NEC$ noinline
    CALL ICE4_NUCLEATION(CST, PARAMI, ICEP, ICED, LDCOMPUTE(JL), &
                     ZTH(JL), PPRES(JL), PRHODREF(JL), PEXN(JL), PLSFACT(JL), ZT(JL), &
                     ZVART(JL,IRV), &
                     PCIT(JL), PBU_INST(JL, IRVHENI_MR))
  ENDDO
  DO JL=1, KSIZE
    ZTH(JL)=ZTH(JL) + PBU_INST(JL, IRVHENI_MR)*PLSFACT(JL)
    ZT(JL) = ZTH(JL) * PEXN(JL)
    ZVART(JL,IRV)=ZVART(JL,IRV) - PBU_INST(JL, IRVHENI_MR)
    ZVART(JL,IRI)=ZVART(JL,IRI) + PBU_INST(JL, IRVHENI_MR)
  ENDDO
!$acc end kernels
  !
  !*       3.3     compute the spontaneous freezing source: RRHONG
  !
  CALL ICE4_RRHONG(CST, PARAMI, ICED, KPROMA, KSIZE, LDCOMPUTE, &
                  &PEXN, PLVFACT, PLSFACT, &
                  &ZT, ZVART(:,IRR), &
                  &ZTH(:), &
                  &PBU_INST(:, IRRHONG_MR))
!$acc kernels
!$acc loop independent
  DO JL=1, KSIZE
    ZTH(JL) = ZTH(JL) + PBU_INST(JL, IRRHONG_MR)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(RRHONG))
    ZT(JL) = ZTH(JL) * PEXN(JL)
    ZVART(JL,IRR) = ZVART(JL,IRR) - PBU_INST(JL, IRRHONG_MR)
    ZVART(JL,IRG) = ZVART(JL,IRG) + PBU_INST(JL, IRRHONG_MR)
  ENDDO
!$acc end kernels
  !
  !*       7.1    cloud ice melting
  !
  CALL ICE4_RIMLTC(CST, PARAMI, KPROMA, KSIZE, LDCOMPUTE, &
                  &PEXN, PLVFACT, PLSFACT, &
                  &ZT, &
                  &ZTH(:), ZVART(:,IRI), &
                  &PBU_INST(:, IRIMLTC_MR))
!$acc kernels
!$acc loop independent
  DO JL=1, KSIZE
    ZTH(JL) = ZTH(JL) - PBU_INST(JL, IRIMLTC_MR)*(PLSFACT(JL)-PLVFACT(JL)) ! f(L_f*(-RIMLTC))
    ZT(JL) = ZTH(JL) * PEXN(JL)
    ZVART(JL,IRC) = ZVART(JL,IRC) + PBU_INST(JL, IRIMLTC_MR)
    ZVART(JL,IRI) = ZVART(JL,IRI) - PBU_INST(JL, IRIMLTC_MR)
  ENDDO
!$acc end kernels
  !
  !        5.1.6  riming-conversion of the large sized aggregates into graupel (old parametrisation)
  !
  IF(PARAMI%CSNOWRIMING=='OLD ') THEN
!$acc kernels
    !$mnh_expand_where(JL=1:KSIZE)
     IF (PARAMI%LSNOW_T) THEN 
        ZLBDAS(1:KSIZE)=0.
        WHERE (ZVART(1:KSIZE,IRS)>0. .AND. ZT(1:KSIZE)>263.15)
           ZLBDAS(1:KSIZE) = MAX(MIN(ICED%XLBDAS_MAX, 10**(14.554-0.0423*ZT(1:KSIZE))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
        END WHERE
        WHERE (ZVART(1:KSIZE,IRS)>0. .AND. ZT(1:KSIZE)<=263.15)
           ZLBDAS(1:KSIZE) = MAX(MIN(ICED%XLBDAS_MAX, 10**(6.226-0.0106*ZT(1:KSIZE))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
        END WHERE
     ELSE
      WHERE(ZVART(1:KSIZE,IRS)>0.)
        ZLBDAS(1:KSIZE)  = MIN(ICED%XLBDAS_MAX, ICED%XLBS*(PRHODREF(1:KSIZE)*MAX(ZVART(1:KSIZE,IRS), ICED%XRTMIN(5)))**ICED%XLBEXS)
      ELSEWHERE
        ZLBDAS(1:KSIZE)=0.
      END WHERE
     END IF
    !$mnh_end_expand_where(JL=1:KSIZE)
!$acc end kernels
    CALL ICE4_RSRIMCG_OLD(CST, PARAMI, ICEP, ICED, KPROMA, KSIZE, ODSOFT, LDCOMPUTE, &
                         &PRHODREF, &
                         &ZLBDAS, &
                         &ZT, ZVART(:,IRC), ZVART(:,IRS), &
                         &PBU_INST(:, IRSRIMCG_MR))
!$acc kernels
!$acc loop independent
    DO JL=1, KSIZE
      ZVART(JL,IRS) = ZVART(JL,IRS) - PBU_INST(JL, IRSRIMCG_MR)
      ZVART(JL,IRG) = ZVART(JL,IRG) + PBU_INST(JL, IRSRIMCG_MR)
    ENDDO
!$acc end kernels
  ELSE
!$acc kernels
    PBU_INST(:, IRSRIMCG_MR) = 0.
!$acc end kernels
  ENDIF
!
!$acc kernels
!$acc loop independent
  DO JL=1, KSIZE
    PBTH(JL)=PBTH(JL) + PBU_INST(JL, IRVHENI_MR)*PLSFACT(JL)
    PBTH(JL)=PBTH(JL) + PBU_INST(JL, IRRHONG_MR)*(PLSFACT(JL)-PLVFACT(JL))
    PBTH(JL)=PBTH(JL) - PBU_INST(JL, IRIMLTC_MR)*(PLSFACT(JL)-PLVFACT(JL))

    PB(JL, IRV)=PB(JL, IRV) - PBU_INST(JL, IRVHENI_MR)

    PB(JL, IRC)=PB(JL, IRC) + PBU_INST(JL, IRIMLTC_MR)

    PB(JL, IRR)=PB(JL, IRR) - PBU_INST(JL, IRRHONG_MR)

    PB(JL, IRI)=PB(JL, IRI) + PBU_INST(JL, IRVHENI_MR)
    PB(JL, IRI)=PB(JL, IRI) - PBU_INST(JL, IRIMLTC_MR)

    PB(JL, IRS)=PB(JL, IRS) - PBU_INST(JL, IRSRIMCG_MR)

    PB(JL, IRG)=PB(JL, IRG) + PBU_INST(JL, IRRHONG_MR)
    PB(JL, IRG)=PB(JL, IRG) + PBU_INST(JL, IRSRIMCG_MR)
  ENDDO
!$acc end kernels
 !
  !* Derived fields
  !
!$acc kernels
!$acc loop independent
  DO JL=1, KSIZE
    ZZW = EXP(CST%XALPI-CST%XBETAI/ZT(JL)-CST%XGAMI*LOG(ZT(JL)))
    PSSI(JL) = ZVART(JL,IRV)*( PPRES(JL)-ZZW ) / ( CST%XEPSILO * ZZW ) - 1.0
                                                      ! Supersaturation over ice
    ZKA(JL) = 2.38E-2 + 0.0071E-2*(ZT(JL)-CST%XTT) ! k_a
    ZDV(JL) = 0.211E-4*(ZT(JL)/CST%XTT)**1.94 * (CST%XP00/PPRES(JL)) ! D_v
    ZAI(JL) = (CST%XLSTT+(CST%XCPV-CST%XCI)*(ZT(JL)-CST%XTT))**2 / (ZKA(JL)*CST%XRV*ZT(JL)**2) &
                                 + ( CST%XRV*ZT(JL) ) / (ZDV(JL)*ZZW)
    ZCJ(JL) = ICEP%XSCFAC*PRHODREF(JL)**0.3 / SQRT(1.718E-5+0.0049E-5*(ZT(JL)-CST%XTT))
  ENDDO
!$acc end kernels
ENDIF ! ODSOFT
!
!Cloud water split between high and low content part is done here
CALL ICE4_COMPUTE_PDF(CST, ICEP, ICED, KSIZE, PARAMI%CSUBG_AUCV_RC, PARAMI%CSUBG_AUCV_RI, PARAMI%CSUBG_PR_PDF,&
                      LDCOMPUTE, PRHODREF, ZVART(:,IRC), ZVART(:,IRI), PCF, ZT, PSIGMA_RC, &
                      PHLC_HCF, PHLC_LCF, PHLC_HRC, PHLC_LRC, &
                      PHLI_HCF, PHLI_LCF, PHLI_HRI, PHLI_LRI, ZRAINFR)
LLRFR=PARAMI%CSUBG_RC_RR_ACCR=='PRFR' .OR. PARAMI%CSUBG_RR_EVAP=='PRFR'
IF (LLRFR) THEN
  !To be exact, ICE4_RAINFR_VERT should be called here with the updated PRAINFR
  !But this call would require the full 3D arrays for rain, snow and graupel which
  !are not available here (due to separation between 1D and 3D computation for GPU).
  !
  !We replace the full computation by a small update to ensure consistency
!$acc kernels
!$acc loop independent
  DO JL=1, KSIZE
    PRAINFR(JL)=MAX(PRAINFR(JL), ZRAINFR(JL))
    IF(KRR==7) THEN
      LLMASK=ZVART(JL,IRR) .GT. ICED%XRTMIN(3) .OR. ZVART(JL,IRS) .GT. ICED%XRTMIN(5) .OR. &
            &ZVART(JL,IRG) .GT. ICED%XRTMIN(6) .OR. ZVART(JL,IRH) .GT. ICED%XRTMIN(7)
    ELSE
      LLMASK=ZVART(JL,IRR) .GT. ICED%XRTMIN(3) .OR. ZVART(JL,IRS) .GT. ICED%XRTMIN(5) .OR. &
            &ZVART(JL,IRG) .GT. ICED%XRTMIN(6)
    ENDIF
    IF(LLMASK .AND. PRAINFR(JL)==0.) THEN
      PRAINFR(JL)=1.
    ENDIF
  ENDDO
!$acc end kernels
ELSE
!$acc kernels
  PRAINFR(:)=1.
!$acc end kernels
ENDIF
!
!*  compute the slope parameters
!
!$acc kernels
!$acc loop independent
DO JL=1, KSIZE
  !ZLBDAR will be used when we consider rain diluted over the grid box
  IF(ZVART(JL,IRR)>0.) THEN
    ZLBDAR(JL)=ICED%XLBR*(PRHODREF(JL)*MAX(ZVART(JL,IRR), ICED%XRTMIN(3)))**ICED%XLBEXR
  ELSE
    ZLBDAR(JL)=0.
  ENDIF
  !ZLBDAR_RF is used when we consider rain concentrated in its fraction
  IF(LLRFR) THEN
    IF(ZVART(JL,IRR)>0. .AND. PRAINFR(JL)>0.) THEN
      ZLBDAR_RF(JL)=ICED%XLBR*(PRHODREF(JL)*MAX(ZVART(JL,IRR)/PRAINFR(JL), ICED%XRTMIN(3)))**ICED%XLBEXR
    ELSE
      ZLBDAR_RF(JL)=0.
    ENDIF
  ELSE
    ZLBDAR_RF(JL)=ZLBDAR(JL)
  ENDIF
  IF (PARAMI%LSNOW_T) THEN 
   IF (ZVART(JL,IRS)>0. .AND. ZT(JL)>263.15) THEN
      ZLBDAS(JL) = MAX(MIN(ICED%XLBDAS_MAX, 10**(14.554-0.0423*ZT(JL))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
   ELSE IF (ZVART(JL,IRS)>0. .AND. ZT(JL)<=263.15) THEN
      ZLBDAS(JL) = MAX(MIN(ICED%XLBDAS_MAX, 10**(6.226-0.0106*ZT(JL))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
   ELSE
      ZLBDAS(JL)=0.
   END IF
  ELSE
   IF(ZVART(JL,IRS)>0.) THEN
    ZLBDAS(JL)=MIN(ICED%XLBDAS_MAX, ICED%XLBS*(PRHODREF(JL)*MAX(ZVART(JL,IRS), ICED%XRTMIN(5)))**ICED%XLBEXS)
   ELSE
    ZLBDAS(JL)=0.
   ENDIF
  END IF
  IF(ZVART(JL,IRG)>0.) THEN
    ZLBDAG(JL)=ICED%XLBG*(PRHODREF(JL)*MAX(ZVART(JL,IRG), ICED%XRTMIN(6)))**ICED%XLBEXG
  ELSE
    ZLBDAG(JL)=0.
  ENDIF
  IF(KRR==7) THEN
    IF(ZVART(JL,IRH)>0.) THEN
      ZLBDAH(JL)=ICED%XLBH*(PRHODREF(JL)*MAX(ZVART(JL,IRH), ICED%XRTMIN(7)))**ICED%XLBEXH
    ELSE
      ZLBDAH(JL)=0.
    ENDIF
  ENDIF
ENDDO
!$acc end kernels
!
!
CALL ICE4_SLOW(CST, ICEP, ICED, KPROMA, KSIZE, ODSOFT, OELEC, LDCOMPUTE, PRHODREF, ZT, &
              &PSSI, PLVFACT, PLSFACT, &
              &ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), &
              &ZLBDAS, ZLBDAG, &
              &ZAI, ZCJ, PHLI_HCF, PHLI_HRI, &
              &PLATHAM_IAGGS, &
              &PBU_INST(:, IRCHONI), PBU_INST(:, IRVDEPS), PBU_INST(:, IRIAGGS), PBU_INST(:, IRIAUTS), PBU_INST(:, IRVDEPG))
!
!-------------------------------------------------------------------------------
!
!
!*       3.     COMPUTES THE SLOW WARM PROCESS SOURCES
!               --------------------------------------
!
!
IF(PARAMI%LWARM) THEN    !  Check if the formation of the raindrops by the slow
                  !  warm processes is allowed
  CALL ICE4_WARM(CST, ICEP, ICED, KPROMA, KSIZE, ODSOFT,LDCOMPUTE, &
                &PARAMI%CSUBG_RC_RR_ACCR, PARAMI%CSUBG_RR_EVAP, &
                &PRHODREF, PLVFACT, ZT, PPRES, ZTH(:),&
                &ZLBDAR, ZLBDAR_RF, ZKA, ZDV, ZCJ, &
                &PHLC_LCF, PHLC_HCF, PHLC_LRC, PHLC_HRC, &
                &PCF, PRAINFR, &
                &ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), &
                &PBU_INST(:, IRCAUTR), PBU_INST(:, IRCACCR), PBU_INST(:, IRREVAV))
ELSE
!$acc kernels
  PBU_INST(:, IRCAUTR)=0.
  PBU_INST(:, IRCACCR)=0.
  PBU_INST(:, IRREVAV)=0.
!$acc end kernels
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_s
!               ----------------------------------------------
!
CALL ICE4_FAST_RS(CST, PARAMI, ICEP, ICED, KPROMA, KSIZE, ODSOFT, LDCOMPUTE, &
                 &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                 &ZDV, ZKA, ZCJ, &
                 &ZLBDAR, ZLBDAS, &
                 &ZT, ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), ZVART(:,IRS), &
                 &PBU_INST(:, IRIAGGS), &
                 &PBU_INST(:, IRCRIMSS), PBU_INST(:, IRCRIMSG), PBU_INST(:, IRSRIMCG), &
                 &PBU_INST(:, IRRACCSS), PBU_INST(:, IRRACCSG), PBU_INST(:, IRSACCRG), PBU_INST(:, IRSMLTG), &
                 &PBU_INST(:, IRCMLTSR), &
                 &PRS_TEND)
!
!-------------------------------------------------------------------------------
!
!
!*       5.        COMPUTES THE FAST COLD PROCESS SOURCES FOR r_g
!                  ------------------------------------------------------
!
!$acc kernels
!$acc loop independent
DO JL=1, KSIZE
  ZRGSI(JL) = PBU_INST(JL, IRVDEPG) + PBU_INST(JL, IRSMLTG) + PBU_INST(JL, IRRACCSG) + &
            & PBU_INST(JL, IRSACCRG) + PBU_INST(JL, IRCRIMSG) + PBU_INST(JL, IRSRIMCG)
  ZRGSI_MR(JL) = PBU_INST(JL, IRRHONG_MR) + PBU_INST(JL, IRSRIMCG_MR)
ENDDO
!$acc end kernels
CALL ICE4_FAST_RG(CST, PARAMI, ICEP, ICED, KPROMA, KSIZE, ODSOFT, LDCOMPUTE, KRR, &
                 &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                 &ZDV, ZKA, ZCJ, PCIT, &
                 &ZLBDAR, ZLBDAS, ZLBDAG, &
                 &ZT, ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), &
                 &ZRGSI, ZRGSI_MR(:), &
                 &LLWETG, &
                 &PBU_INST(:, IRICFRRG), PBU_INST(:, IRRCFRIG), PBU_INST(:, IRICFRR), PBU_INST(:, IRCWETG), &
                 &PBU_INST(:, IRIWETG), PBU_INST(:, IRRWETG), PBU_INST(:, IRSWETG), &
                 &PBU_INST(:, IRCDRYG), PBU_INST(:, IRIDRYG), PBU_INST(:, IRRDRYG), PBU_INST(:, IRSDRYG), &
                 &PBU_INST(:, IRWETGH), PBU_INST(:, IRWETGH_MR), PBU_INST(:, IRGMLTR), &
                 &PRG_TEND)
!
!-------------------------------------------------------------------------------
!
!
!*       6.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_h
!               ----------------------------------------------
!
IF (KRR==7) THEN
  CALL ICE4_FAST_RH(CST, PARAMI, ICEP, ICED, KPROMA, KSIZE, ODSOFT, LDCOMPUTE, LLWETG, &
                   &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                   &ZDV, ZKA, ZCJ, &
                   &ZLBDAS, ZLBDAG, ZLBDAR, ZLBDAH, &
                   &ZT,  ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), ZVART(:,IRH), &
                   &PBU_INST(:, IRCWETH), PBU_INST(:, IRIWETH), PBU_INST(:, IRSWETH), PBU_INST(:, IRGWETH), PBU_INST(:, IRRWETH), &
                   &PBU_INST(:, IRCDRYH), PBU_INST(:, IRIDRYH), PBU_INST(:, IRSDRYH), PBU_INST(:, IRRDRYH), &
                   &PBU_INST(:, IRGDRYH), PBU_INST(:, IRDRYHG), PBU_INST(:, IRHMLTR), &
                   &PRH_TEND)
ELSEIF (BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN
!$acc kernels
  PBU_INST(:, IRCWETH)=0.
  PBU_INST(:, IRIWETH)=0.
  PBU_INST(:, IRSWETH)=0.
  PBU_INST(:, IRGWETH)=0.
  PBU_INST(:, IRRWETH)=0.
  PBU_INST(:, IRCDRYH)=0.
  PBU_INST(:, IRIDRYH)=0.
  PBU_INST(:, IRSDRYH)=0.
  PBU_INST(:, IRRDRYH)=0.
  PBU_INST(:, IRGDRYH)=0.
  PBU_INST(:, IRDRYHG)=0.
  PBU_INST(:, IRHMLTR)=0.
!$acc end kernels
END IF
!
!-------------------------------------------------------------------------------
!
!
!*       7.     COMPUTES SPECIFIC SOURCES OF THE WARM AND COLD CLOUDY SPECIES
!               -------------------------------------------------------------
!
CALL ICE4_FAST_RI(ICEP, ICED, KPROMA, KSIZE, ODSOFT, LDCOMPUTE, &
                 &PRHODREF, PLVFACT, PLSFACT, &
                 &ZAI, ZCJ, PCIT, &
                 &PSSI, &
                 &ZVART(:,IRC), ZVART(:,IRI), &
                 &PBU_INST(:, IRCBERI))
!
!-------------------------------------------------------------------------------
!
!
!*       8.     COMPUTES TOTAL TENDENCIES
!               -------------------------
!
!$acc kernels
!$acc loop independent
DO JL=1, KSIZE
  PATH(JL) = PATH(JL) + PBU_INST(JL, IRVDEPG)*PLSFACT(JL)
  PATH(JL) = PATH(JL) + PBU_INST(JL, IRCHONI)*(PLSFACT(JL)-PLVFACT(JL))
  PATH(JL) = PATH(JL) + PBU_INST(JL, IRVDEPS)*PLSFACT(JL)
  PATH(JL) = PATH(JL) - PBU_INST(JL, IRREVAV)*PLVFACT(JL)
  PATH(JL) = PATH(JL) + PBU_INST(JL, IRCRIMSS)*(PLSFACT(JL)-PLVFACT(JL))
  PATH(JL) = PATH(JL) + PBU_INST(JL, IRCRIMSG)*(PLSFACT(JL)-PLVFACT(JL))
  PATH(JL) = PATH(JL) + PBU_INST(JL, IRRACCSS)*(PLSFACT(JL)-PLVFACT(JL))
  PATH(JL) = PATH(JL) + PBU_INST(JL, IRRACCSG)*(PLSFACT(JL)-PLVFACT(JL))
  PATH(JL) = PATH(JL) + (PBU_INST(JL, IRRCFRIG) - PBU_INST(JL, IRICFRR))*(PLSFACT(JL)-PLVFACT(JL))
  PATH(JL) = PATH(JL) + (PBU_INST(JL, IRCWETG) + PBU_INST(JL, IRRWETG))*(PLSFACT(JL)-PLVFACT(JL))
  PATH(JL) = PATH(JL) + (PBU_INST(JL, IRCDRYG)+PBU_INST(JL, IRRDRYG))*(PLSFACT(JL)-PLVFACT(JL))
  PATH(JL) = PATH(JL) - PBU_INST(JL, IRGMLTR)*(PLSFACT(JL)-PLVFACT(JL))
  IF (KRR==7) THEN
    PATH(JL) = PATH(JL) + (PBU_INST(JL, IRRWETH)+PBU_INST(JL, IRCWETH))*(PLSFACT(JL)-PLVFACT(JL))
    PATH(JL) = PATH(JL) + (PBU_INST(JL, IRCDRYH)+PBU_INST(JL, IRRDRYH))*(PLSFACT(JL)-PLVFACT(JL))
    PATH(JL) = PATH(JL) - PBU_INST(JL, IRHMLTR)*(PLSFACT(JL)-PLVFACT(JL))
  ENDIF
  PATH(JL) = PATH(JL) + PBU_INST(JL, IRCBERI)*(PLSFACT(JL)-PLVFACT(JL))

  PA(JL, IRV) = PA(JL, IRV) - PBU_INST(JL, IRVDEPG)
  PA(JL, IRV) = PA(JL, IRV) - PBU_INST(JL, IRVDEPS)
  PA(JL, IRV) = PA(JL, IRV) + PBU_INST(JL, IRREVAV)

  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCHONI)
  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCAUTR)
  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCACCR)
  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCRIMSS)
  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCRIMSG)
  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCMLTSR)
  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCWETG)
  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCDRYG)
  IF (KRR==7) THEN
    PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCWETH)
    PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCDRYH)
  ENDIF
  PA(JL, IRC) = PA(JL, IRC) - PBU_INST(JL, IRCBERI)

  PA(JL, IRR) = PA(JL, IRR) + PBU_INST(JL, IRCAUTR)
  PA(JL, IRR) = PA(JL, IRR) + PBU_INST(JL, IRCACCR)
  PA(JL, IRR) = PA(JL, IRR) - PBU_INST(JL, IRREVAV)
  PA(JL, IRR) = PA(JL, IRR) - PBU_INST(JL, IRRACCSS)
  PA(JL, IRR) = PA(JL, IRR) - PBU_INST(JL, IRRACCSG)
  PA(JL, IRR) = PA(JL, IRR) + PBU_INST(JL, IRCMLTSR)
  PA(JL, IRR) = PA(JL, IRR) - PBU_INST(JL, IRRCFRIG) + PBU_INST(JL, IRICFRR)
  PA(JL, IRR) = PA(JL, IRR) - PBU_INST(JL, IRRWETG)
  PA(JL, IRR) = PA(JL, IRR) - PBU_INST(JL, IRRDRYG)
  PA(JL, IRR) = PA(JL, IRR) + PBU_INST(JL, IRGMLTR)
  IF(KRR==7) THEN
    PA(JL, IRR) = PA(JL, IRR) - PBU_INST(JL, IRRWETH)
    PA(JL, IRR) = PA(JL, IRR) - PBU_INST(JL, IRRDRYH)
    PA(JL, IRR) = PA(JL, IRR) + PBU_INST(JL, IRHMLTR)
  ENDIF

  PA(JL, IRI) = PA(JL, IRI) + PBU_INST(JL, IRCHONI)
  PA(JL, IRI) = PA(JL, IRI) - PBU_INST(JL, IRIAGGS)
  PA(JL, IRI) = PA(JL, IRI) - PBU_INST(JL, IRIAUTS)
  PA(JL, IRI) = PA(JL, IRI) - PBU_INST(JL, IRICFRRG) - PBU_INST(JL, IRICFRR)
  PA(JL, IRI) = PA(JL, IRI) - PBU_INST(JL, IRIWETG)
  PA(JL, IRI) = PA(JL, IRI) - PBU_INST(JL, IRIDRYG)
  IF (KRR==7) THEN
    PA(JL, IRI) = PA(JL, IRI) - PBU_INST(JL, IRIWETH)
    PA(JL, IRI) = PA(JL, IRI) - PBU_INST(JL, IRIDRYH)
  ENDIF
  PA(JL, IRI) = PA(JL, IRI) + PBU_INST(JL, IRCBERI)

  PA(JL, IRS) = PA(JL, IRS) + PBU_INST(JL, IRVDEPS)
  PA(JL, IRS) = PA(JL, IRS) + PBU_INST(JL, IRIAGGS)
  PA(JL, IRS) = PA(JL, IRS) + PBU_INST(JL, IRIAUTS)
  PA(JL, IRS) = PA(JL, IRS) + PBU_INST(JL, IRCRIMSS)
  PA(JL, IRS) = PA(JL, IRS) - PBU_INST(JL, IRSRIMCG)
  PA(JL, IRS) = PA(JL, IRS) + PBU_INST(JL, IRRACCSS)
  PA(JL, IRS) = PA(JL, IRS) - PBU_INST(JL, IRSACCRG)
  PA(JL, IRS) = PA(JL, IRS) - PBU_INST(JL, IRSMLTG)
  PA(JL, IRS) = PA(JL, IRS) - PBU_INST(JL, IRSWETG)
  PA(JL, IRS) = PA(JL, IRS) - PBU_INST(JL, IRSDRYG)
  IF (KRR==7) THEN
    PA(JL, IRS) = PA(JL, IRS) - PBU_INST(JL, IRSWETH)
    PA(JL, IRS) = PA(JL, IRS) - PBU_INST(JL, IRSDRYH)
  ENDIF

  PA(JL, IRG) = PA(JL, IRG) + PBU_INST(JL, IRVDEPG)
  PA(JL, IRG) = PA(JL, IRG) + PBU_INST(JL, IRCRIMSG)+PBU_INST(JL, IRSRIMCG)
  PA(JL, IRG) = PA(JL, IRG) + PBU_INST(JL, IRRACCSG)+PBU_INST(JL, IRSACCRG)
  PA(JL, IRG) = PA(JL, IRG) + PBU_INST(JL, IRSMLTG)
  PA(JL, IRG) = PA(JL, IRG) + PBU_INST(JL, IRICFRRG) + PBU_INST(JL, IRRCFRIG)
  PA(JL, IRG) = PA(JL, IRG) + PBU_INST(JL, IRCWETG) + PBU_INST(JL, IRIWETG) + PBU_INST(JL, IRSWETG) + PBU_INST(JL, IRRWETG)
  PA(JL, IRG) = PA(JL, IRG) - PBU_INST(JL, IRWETGH)
  PB(JL, IRG) = PB(JL, IRG) - PBU_INST(JL, IRWETGH_MR)
  PA(JL, IRG) = PA(JL, IRG) + PBU_INST(JL, IRCDRYG) + PBU_INST(JL, IRIDRYG) + PBU_INST(JL, IRSDRYG) + PBU_INST(JL, IRRDRYG)
  PA(JL, IRG) = PA(JL, IRG) - PBU_INST(JL, IRGMLTR)
  IF (KRR==7) THEN
    PA(JL, IRG) = PA(JL, IRG) - PBU_INST(JL, IRGWETH)
    PA(JL, IRG) = PA(JL, IRG) - PBU_INST(JL, IRGDRYH) + PBU_INST(JL, IRDRYHG)
  ENDIF

  IF (KRR==7) THEN
    PA(JL, IRH) = PA(JL, IRH) + PBU_INST(JL, IRWETGH)
    PB(JL, IRH) = PB(JL, IRH) + PBU_INST(JL, IRWETGH_MR)
    PA(JL, IRH) = PA(JL, IRH) + PBU_INST(JL, IRCWETH)+PBU_INST(JL, IRIWETH)+PBU_INST(JL, IRSWETH)+&
                              & PBU_INST(JL, IRGWETH)+PBU_INST(JL, IRRWETH)
    PA(JL, IRH) = PA(JL, IRH) + PBU_INST(JL, IRCDRYH)+PBU_INST(JL, IRIDRYH)+PBU_INST(JL, IRSDRYH)+&
                              & PBU_INST(JL, IRRDRYH)+PBU_INST(JL, IRGDRYH) - PBU_INST(JL, IRDRYHG)
    PA(JL, IRH) = PA(JL, IRH) - PBU_INST(JL, IRHMLTR)
  ENDIF
ENDDO
!$acc end kernels
!
IF (LHOOK) CALL DR_HOOK('ICE4_TENDENCIES', 1, ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "ice4_nucleation.func.h"
END SUBROUTINE ICE4_TENDENCIES
END MODULE MODE_ICE4_TENDENCIES
