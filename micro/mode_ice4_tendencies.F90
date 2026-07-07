!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_TENDENCIES
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_TENDENCIES(CST, PARAMI, ICEP, ICED, BUCONF, D, &
                          &KRR, ODSOFT, LDCOMPUTE, &
                          &OSAVE_MICRO, OELEC, &
                          &PEXN, PRHODREF, PLVFACT, PLSFACT, &
                          &PPRES, PCF, PSIGMA_RC, &
                          &PCIT, &
                          &PT, PICLDFR, PZZZ, PCONC, &
                          &PSSIO, PSSIU, PIFR, PTH, PVART, &
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
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_BUDGET,         ONLY: TBUDGETCONF_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
!
USE MODD_FIELDS_ADDRESS, ONLY: IBUNUM, IRC, IRCACCR, IRCAUTR, IRCBERI, IRCDRYG, IRCDRYH, IRCHONI, &
                               IRCMLTSR, IRCRIMSG, IRCRIMSS, IRCWETG, IRCWETH, IRDRYHG, IRG, IRGDRYH, &
                               IRGMLTR, IRGWETH, IRH, IRHMLTR, IRI, IRIAGGS, IRIAUTS, IRICFRR, &
                               IRICFRRG, IRIDRYG, IRIDRYH, IRILARS, IRIMLTC_MR, IRIWETG, IRIWETH, &
                               IRR, IRRACCSG, IRRACCSS, IRRCFRIG, IRRDRYG, IRRDRYH, IRREVAV, &
                               IRRHONG_MR, IRRWETG, IRRWETH, IRS, IRSACCRG, IRSDRYG, IRSDRYH, &
                               IRSMLTG, IRSRIMCG, IRSRIMCG_MR, IRSWETG, IRSWETH, IRV, IRVDEPG, &
                               IRVDEPI, IRVDEPS, IRVHENI_MR, IRWETGH, IRWETGH_MR
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
USE MODE_ICE4_FAST_RI_RS, ONLY: ICE4_FAST_RI_RS
USE MODE_TIWMX,      ONLY: ESATI, AA2, BB3 
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
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
INTEGER,                      INTENT(IN)    :: KRR
LOGICAL,                      INTENT(IN)    :: ODSOFT
LOGICAL, DIMENSION(D%NIJT),   INTENT(IN)    :: LDCOMPUTE
LOGICAL,                      INTENT(IN)    :: OSAVE_MICRO
LOGICAL,                      INTENT(IN)    :: OELEC
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PEXN
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PRHODREF
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PPRES
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PCF
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PSIGMA_RC
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PCIT
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PT
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PICLDFR
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PZZZ
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PCONC
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PSSIO
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PSSIU
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PIFR
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PTH
REAL, DIMENSION(D%NIJT,KRR),   INTENT(IN)    :: PVART
REAL, DIMENSION(MERGE(D%NIJT,0,OELEC)), INTENT(IN) :: PLATHAM_IAGGS
REAL, DIMENSION(D%NIJT, IBUNUM),INTENT(INOUT):: PBU_INST
REAL, DIMENSION(D%NIJT, 8),    INTENT(INOUT) :: PRS_TEND
REAL, DIMENSION(D%NIJT, 8),    INTENT(INOUT) :: PRG_TEND
REAL, DIMENSION(D%NIJT, 10),   INTENT(INOUT) :: PRH_TEND
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PSSI
REAL, DIMENSION(D%NIJT,KRR),   INTENT(OUT)   :: PA
REAL, DIMENSION(D%NIJT,KRR),   INTENT(OUT)   :: PB
REAL, DIMENSION(D%NIJT),       INTENT(OUT)   :: PATH
REAL, DIMENSION(D%NIJT),       INTENT(OUT)   :: PBTH
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PHLC_HCF
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PHLC_LCF
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PHLC_HRC
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PHLC_LRC
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PHLI_HCF
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PHLI_LCF
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PHLI_HRI
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PHLI_LRI
REAL, DIMENSION(D%NIJT),       INTENT(INOUT) :: PRAINFR   ! Rain fraction
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,KRR) :: ZVART
REAL, DIMENSION(D%NIJT) :: ZT, ZTH, &
                        & ZKA, ZDV, ZAI, ZCJ, &
                        & ZLBDAR, ZLBDAS, ZLBDAG, ZLBDAH, ZLBDAR_RF, &
                        & ZRGSI, ZRGSI_MR, &
                        & ZRAINFR, &
                        & ZCOLF, ZACRF, ZESI, ZRCW, ZVT, ZST, ZR20, ZCONC
INTEGER :: JIJ, JV
LOGICAL, DIMENSION(D%NIJT) :: LLWETG ! .TRUE. if graupel growths in wet mode
REAL, DIMENSION(D%NIJT) :: ZZW, ZDEPG_S, ZNODEPG_S
LOGICAL :: LLMASK
LOGICAL :: LLRFR
LOGICAL :: LLTAB   ! Use saturation vapour pressure tables for optimization 
LOGICAL :: LLCOL   ! Collision factors OCND2-style
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE

IF (LHOOK) CALL DR_HOOK('ICE4_TENDENCIES', 0, ZHOOK_HANDLE)
!
LLTAB = .TRUE.
LLCOL = .TRUE.
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!

ZT(:)=PT(:)
ZTH(:)=PTH(:)
PATH(:)=0.
PBTH(:)=0.
DO JV=1,KRR
  ZVART(:,JV)=PVART(:,JV)
  PA(:,JV)=0.
  PB(:,JV)=0.
ENDDO

!
IF(ODSOFT) THEN

  PBU_INST(:, IRVHENI_MR)=0.
  PBU_INST(:, IRRHONG_MR)=0.
  PBU_INST(:, IRIMLTC_MR)=0.
  PBU_INST(:, IRSRIMCG_MR)=0.

ELSE
  !
  !*       2.     COMPUTES THE SLOW COLD PROCESS SOURCES
  !               --------------------------------------

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
!NEC$ noinline
    CALL ICE4_NUCLEATION(CST, PARAMI, ICEP, ICED, LDCOMPUTE(JIJ), &
                     ZTH(JIJ), PPRES(JIJ), PRHODREF(JIJ), PEXN(JIJ), PLSFACT(JIJ), ZT(JIJ), &
                     ZVART(JIJ,IRV), PICLDFR(JIJ), PZZZ(JIJ), &
                     PCIT(JIJ), PBU_INST(JIJ, IRVHENI_MR))
  ENDDO
  !$mnh_end_do()
  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    ZTH(JIJ)=ZTH(JIJ) + PBU_INST(JIJ, IRVHENI_MR)*PLSFACT(JIJ)
    ZT(JIJ) = ZTH(JIJ) * PEXN(JIJ)
    ZVART(JIJ,IRV)=ZVART(JIJ,IRV) - PBU_INST(JIJ, IRVHENI_MR)
    ZVART(JIJ,IRI)=ZVART(JIJ,IRI) + PBU_INST(JIJ, IRVHENI_MR)
  ENDDO
  !$mnh_end_do()

  !
  !*       3.3     compute the spontaneous freezing source: RRHONG
  !
  CALL ICE4_RRHONG(CST, PARAMI, ICED, D, LDCOMPUTE, &
                  &PEXN, PLVFACT, PLSFACT, &
                  &ZT, ZVART(:,IRR), &
                  &ZTH, &
                  &PBU_INST(:, IRRHONG_MR))

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    ZTH(JIJ) = ZTH(JIJ) + PBU_INST(JIJ, IRRHONG_MR)*(PLSFACT(JIJ)-PLVFACT(JIJ)) ! f(L_f*(RRHONG))
    ZT(JIJ) = ZTH(JIJ) * PEXN(JIJ)
    ZVART(JIJ,IRR) = ZVART(JIJ,IRR) - PBU_INST(JIJ, IRRHONG_MR)
    ZVART(JIJ,IRG) = ZVART(JIJ,IRG) + PBU_INST(JIJ, IRRHONG_MR)
  ENDDO
  !$mnh_end_do()

  !
  !*       7.1    cloud ice melting
  !
  CALL ICE4_RIMLTC(CST, PARAMI, D, LDCOMPUTE, &
                  &PEXN, PLVFACT, PLSFACT, &
                  &ZT, &
                  &ZTH, ZVART(:,IRI), &
                  &PBU_INST(:, IRIMLTC_MR))

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    ZTH(JIJ) = ZTH(JIJ) - PBU_INST(JIJ, IRIMLTC_MR)*(PLSFACT(JIJ)-PLVFACT(JIJ)) ! f(L_f*(-RIMLTC))
    ZT(JIJ) = ZTH(JIJ) * PEXN(JIJ)
    ZVART(JIJ,IRC) = ZVART(JIJ,IRC) + PBU_INST(JIJ, IRIMLTC_MR)
    ZVART(JIJ,IRI) = ZVART(JIJ,IRI) - PBU_INST(JIJ, IRIMLTC_MR)
  ENDDO
  !$mnh_end_do()

  !
  !        5.1.6  riming-conversion of the large sized aggregates into graupel (old parametrisation)
  !
  IF(PARAMI%CSNOWRIMING=='OLD ') THEN

    DO JIJ=D%NIJB, D%NIJE
       IF (PARAMI%LSNOW_T) THEN 
          ZLBDAS(JIJ)=0.
          IF (ZVART(JIJ, IRS)>0. .AND. ZT(JIJ)>263.15) THEN
             ZLBDAS(JIJ) = MAX(MIN(ICED%XLBDAS_MAX, 10**(14.554-0.0423*ZT(JIJ))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
          END IF
          IF (ZVART(JIJ, IRS)>0. .AND. ZT(JIJ)<=263.15) THEN
             ZLBDAS(JIJ) = MAX(MIN(ICED%XLBDAS_MAX, 10**(6.226-0.0106*ZT(JIJ))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
          END IF
       ELSE
        IF (ZVART(JIJ, IRS)>0.) THEN
          ZLBDAS(JIJ)  = MIN(ICED%XLBDAS_MAX, ICED%XLBS*(PRHODREF(JIJ)*MAX(ZVART(JIJ, IRS), ICED%XRTMIN(5)))**ICED%XLBEXS)
        ELSE
          ZLBDAS(JIJ)=0.
        END IF
       END IF
    END DO

    CALL ICE4_RSRIMCG_OLD(CST, PARAMI, ICEP, ICED, D, ODSOFT, LDCOMPUTE, &
                         &PRHODREF, &
                         &ZLBDAS, &
                         &ZT, ZVART(:,IRC), ZVART(:,IRS), &
                         &PBU_INST(:, IRSRIMCG_MR))

    !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
    DO JIJ=D%NIJB, D%NIJE
      ZVART(JIJ,IRS) = ZVART(JIJ,IRS) - PBU_INST(JIJ, IRSRIMCG_MR)
      ZVART(JIJ,IRG) = ZVART(JIJ,IRG) + PBU_INST(JIJ, IRSRIMCG_MR)
    ENDDO
    !$mnh_end_do()

  ELSE

    PBU_INST(:, IRSRIMCG_MR) = 0.

  ENDIF
!

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    PBTH(JIJ)=PBTH(JIJ) + PBU_INST(JIJ, IRVHENI_MR)*PLSFACT(JIJ)
    PBTH(JIJ)=PBTH(JIJ) + PBU_INST(JIJ, IRRHONG_MR)*(PLSFACT(JIJ)-PLVFACT(JIJ))
    PBTH(JIJ)=PBTH(JIJ) - PBU_INST(JIJ, IRIMLTC_MR)*(PLSFACT(JIJ)-PLVFACT(JIJ))

    PB(JIJ, IRV)=PB(JIJ, IRV) - PBU_INST(JIJ, IRVHENI_MR)

    PB(JIJ, IRC)=PB(JIJ, IRC) + PBU_INST(JIJ, IRIMLTC_MR)

    PB(JIJ, IRR)=PB(JIJ, IRR) - PBU_INST(JIJ, IRRHONG_MR)

    PB(JIJ, IRI)=PB(JIJ, IRI) + PBU_INST(JIJ, IRVHENI_MR)
    PB(JIJ, IRI)=PB(JIJ, IRI) - PBU_INST(JIJ, IRIMLTC_MR)

    PB(JIJ, IRS)=PB(JIJ, IRS) - PBU_INST(JIJ, IRSRIMCG_MR)

    PB(JIJ, IRG)=PB(JIJ, IRG) + PBU_INST(JIJ, IRRHONG_MR)
    PB(JIJ, IRG)=PB(JIJ, IRG) + PBU_INST(JIJ, IRSRIMCG_MR)
  ENDDO
!$mnh_end_do()

  !
  !* Derived fields
  !

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    IF(LLTAB.AND.PARAMI%LOCND2) THEN ! Use tabulated values, May go faster for some machines. (Is now separared from OCND2)
      ZESI(JIJ) = ESATI(ICEP%TIWMX, ZT(JIJ))
      PSSI(JIJ) = ZVART(JIJ,IRV)*( PPRES(JIJ)-ZESI(JIJ) ) / ( CST%XEPSILO * ZESI(JIJ) ) - 1.0
                                                    ! Supersaturation over ice0
      ZKA(JIJ) = 2.38E-2 + 0.0071E-2*(ZT(JIJ)-CST%XTT) ! k_a
      ZDV(JIJ) = 0.211E-4*(ZT(JIJ)/CST%XTT)**1.94 * (CST%XP00/PPRES(JIJ)) ! D_v
      ZAI(JIJ) =  AA2(ICEP%TIWMX, ZT(JIJ)) + BB3(ICEP%TIWMX, ZT(JIJ))*PPRES(JIJ)
    ELSE
      ZZW(JIJ) = EXP(CST%XALPI-CST%XBETAI/ZT(JIJ)-CST%XGAMI*LOG(ZT(JIJ)))
      PSSI(JIJ) = ZVART(JIJ,IRV)*( PPRES(JIJ)-ZZW(JIJ) ) / ( CST%XEPSILO * ZZW(JIJ) ) - 1.0
      ZKA(JIJ) = 2.38E-2 + 0.0071E-2*(ZT(JIJ)-CST%XTT) ! k_a
      ZDV(JIJ) = 0.211E-4*(ZT(JIJ)/CST%XTT)**1.94 * (CST%XP00/PPRES(JIJ)) ! D_v
      ZAI(JIJ) = (CST%XLSTT+(CST%XCPV-CST%XCI)*(ZT(JIJ)-CST%XTT))**2 / (ZKA(JIJ)*CST%XRV*ZT(JIJ)**2) &
                                   + ( CST%XRV*ZT(JIJ) ) / (ZDV(JIJ)*ZZW(JIJ))
    ENDIF

    ZCJ(JIJ) = ICEP%XSCFAC*PRHODREF(JIJ)**0.3 / SQRT(1.718E-5+0.0049E-5*(ZT(JIJ)-CST%XTT))

    ! Compute collison factors for cloud water - rain (ZACRF) and cloud water snow/graupel (ZCOLF)
    ZCOLF(JIJ)=1.0
    ZACRF(JIJ)=1.0
    ZCONC(JIJ) = PCONC(JIJ)*0.000001
    IF(LLCOL .AND. PARAMI%LOCND2)THEN
      ZCOLF(JIJ)=0.00001
      ZACRF(JIJ)=0.00001
      IF(ZVART(JIJ,IRC)>1.0E-10)THEN
        ! mean cloud droplet radius in cm
        ZRCW(JIJ) =  0.1*(0.75*ZVART(JIJ,IRC)*PRHODREF(JIJ)/(CST%XPI*ZCONC(JIJ)))**0.333 
        ! fall speed for mean cloud droplet with cloud droplet radius in cm/s
        IF(ZRCW(JIJ) < 0.0065 )THEN
          ZVT(JIJ)   =  1.19E6*ZRCW(JIJ)**2
        ELSE
          ZVT(JIJ)   =  8000.*ZRCW(JIJ)
        ENDIF
        ZVT(JIJ) = MIN(10.,ZVT(JIJ))
        ZST(JIJ) = MAX(0.01,2.*(100.-ZVT(JIJ))*ZVT(JIJ)/(CST%XG*10.))
        IF(ZST(JIJ) > 0.1) ZCOLF(JIJ) =  MAX(0.01,MIN(1.,0.939*ZST(JIJ)**2.657))
        IF( ZVART(JIJ,IRR) > 1.0E-10 .AND. ZRCW(JIJ) >1.0E-5)THEN
          ZR20(JIJ) = EXP(ZRCW(JIJ)*2000.)  ! This ZRCW is in cm . To convert to micro meter : x 10000
          ZACRF(JIJ)  = (ZR20(JIJ) -1.)/(ZR20(JIJ) +1.)             ! ZRCW is then multiplied with 0.2
        ENDIF
      ENDIF
    ENDIF
  ENDDO
  !$mnh_end_do()

ENDIF ! ODSOFT
!
!Cloud water split between high and low content part is done here
CALL ICE4_COMPUTE_PDF(CST, ICEP, ICED, D, PARAMI%CSUBG_AUCV_RC, PARAMI%CSUBG_AUCV_RI, PARAMI%CSUBG_PR_PDF,&
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

  !$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
  DO JIJ=D%NIJB, D%NIJE
    PRAINFR(JIJ)=MAX(PRAINFR(JIJ), ZRAINFR(JIJ))
    IF(KRR==7) THEN
      LLMASK=ZVART(JIJ,IRR) .GT. ICED%XRTMIN(3) .OR. ZVART(JIJ,IRS) .GT. ICED%XRTMIN(5) .OR. &
            &ZVART(JIJ,IRG) .GT. ICED%XRTMIN(6) .OR. ZVART(JIJ,IRH) .GT. ICED%XRTMIN(7)
    ELSE
      LLMASK=ZVART(JIJ,IRR) .GT. ICED%XRTMIN(3) .OR. ZVART(JIJ,IRS) .GT. ICED%XRTMIN(5) .OR. &
            &ZVART(JIJ,IRG) .GT. ICED%XRTMIN(6)
    ENDIF
    IF(LLMASK .AND. PRAINFR(JIJ)==0.) THEN
      PRAINFR(JIJ)=1.
    ENDIF
  ENDDO
  !$mnh_end_do()

ELSE

  PRAINFR(:)=1.

ENDIF
!
!*  compute the slope parameters
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  !ZLBDAR will be used when we consider rain diluted over the grid box
  IF(ZVART(JIJ,IRR)>0.) THEN
    ZLBDAR(JIJ)=ICED%XLBR*(PRHODREF(JIJ)*MAX(ZVART(JIJ,IRR), ICED%XRTMIN(3)))**ICED%XLBEXR
  ELSE
    ZLBDAR(JIJ)=0.
  ENDIF
  !ZLBDAR_RF is used when we consider rain concentrated in its fraction
  IF(LLRFR) THEN
    IF(ZVART(JIJ,IRR)>0. .AND. PRAINFR(JIJ)>0.) THEN
      ZLBDAR_RF(JIJ)=ICED%XLBR*(PRHODREF(JIJ)*MAX(ZVART(JIJ,IRR)/PRAINFR(JIJ), ICED%XRTMIN(3)))**ICED%XLBEXR
    ELSE
      ZLBDAR_RF(JIJ)=0.
    ENDIF
  ELSE
    ZLBDAR_RF(JIJ)=ZLBDAR(JIJ)
  ENDIF
  IF (PARAMI%LSNOW_T) THEN 
   IF (ZVART(JIJ,IRS)>0. .AND. ZT(JIJ)>263.15) THEN
      ZLBDAS(JIJ) = MAX(MIN(ICED%XLBDAS_MAX, 10**(14.554-0.0423*ZT(JIJ))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
   ELSE IF (ZVART(JIJ,IRS)>0. .AND. ZT(JIJ)<=263.15) THEN
      ZLBDAS(JIJ) = MAX(MIN(ICED%XLBDAS_MAX, 10**(6.226-0.0106*ZT(JIJ))),ICED%XLBDAS_MIN)*ICED%XTRANS_MP_GAMMAS
   ELSE
      ZLBDAS(JIJ)=0.
   END IF
  ELSE
   IF(ZVART(JIJ,IRS)>0.) THEN
    ZLBDAS(JIJ)=MIN(ICED%XLBDAS_MAX, ICED%XLBS*(PRHODREF(JIJ)*MAX(ZVART(JIJ,IRS), ICED%XRTMIN(5)))**ICED%XLBEXS)
   ELSE
    ZLBDAS(JIJ)=0.
   ENDIF
  END IF
  IF(ZVART(JIJ,IRG)>0.) THEN
    ZLBDAG(JIJ)=ICED%XLBG*(PRHODREF(JIJ)*MAX(ZVART(JIJ,IRG), ICED%XRTMIN(6)))**ICED%XLBEXG
  ELSE
    ZLBDAG(JIJ)=0.
  ENDIF
  IF(KRR==7) THEN
    IF(ZVART(JIJ,IRH)>0.) THEN
      ZLBDAH(JIJ)=ICED%XLBH*(PRHODREF(JIJ)*MAX(ZVART(JIJ,IRH), ICED%XRTMIN(7)))**ICED%XLBEXH
    ELSE
      ZLBDAH(JIJ)=0.
    ENDIF
  ENDIF
ENDDO
!$mnh_end_do()

!
!
CALL ICE4_SLOW(CST, PARAMI, ICEP, ICED, D, ODSOFT, OELEC, LDCOMPUTE, PRHODREF, ZT, &
              &PSSI, PLVFACT, PLSFACT, &
              &ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), &
              &ZLBDAS, ZLBDAG, &
              &ZAI, ZCJ, PHLI_HCF, PHLI_HRI, &
              &PLATHAM_IAGGS, &
              &PBU_INST(:, IRCHONI), PBU_INST(:, IRVDEPS), PBU_INST(:, IRIAGGS), PBU_INST(:, IRIAUTS), PBU_INST(:, IRVDEPG))
!
ZDEPG_S(:)=0.0
ZNODEPG_S(:)=1.0
IF (PARAMI%LOCND2) THEN
! Possibility to turn some growing graupels to snow in case of high ice supersaturation
! and low graupel content (both conditions should be present simultaneously)
  DO JIJ=D%NIJB, D%NIJE
    IF(ICEP%XFRMIN(5)>1.0E-12 .AND. ICEP%XFRMIN(6)>0.01) THEN
      ZDEPG_S(JIJ)=MAX(0., MIN(1., (ICEP%XFRMIN(5)-ZVART(JIJ,IRG)/ICEP%XFRMIN(5))))*&
                  MAX(0., MIN(1., (PSSI(JIJ)/ICEP%XFRMIN(6))))
      ZNODEPG_S(JIJ)=1.-ZDEPG_S(JIJ)
    ENDIF
  ENDDO
ENDIF
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
  CALL ICE4_WARM(CST, PARAMI, ICEP, ICED, D, ODSOFT,LDCOMPUTE, &
                &PARAMI%CSUBG_RC_RR_ACCR, PARAMI%CSUBG_RR_EVAP, &
                &PRHODREF, PLVFACT, ZT, PPRES, ZTH,&
                &ZLBDAR, ZLBDAR_RF, ZKA, ZDV, ZCJ, &
                &PHLC_LCF, PHLC_HCF, PHLC_LRC, PHLC_HRC, &
                &PCF, PRAINFR, &
                &ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), &
                &ZCONC, ZACRF, &
                &PBU_INST(:, IRCAUTR), PBU_INST(:, IRCACCR), PBU_INST(:, IRREVAV))
ELSE

  PBU_INST(:, IRCAUTR)=0.
  PBU_INST(:, IRCACCR)=0.
  PBU_INST(:, IRREVAV)=0.

END IF
!
!-------------------------------------------------------------------------------
!
!
!*       4.     COMPUTES THE FAST COLD PROCESS SOURCES FOR r_s
!               ----------------------------------------------
!
CALL ICE4_FAST_RS(CST, PARAMI, ICEP, ICED, D, ODSOFT, LDCOMPUTE, &
                 &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                 &ZDV, ZKA, ZCJ, ZCOLF, &
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

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  ZRGSI(JIJ) = PBU_INST(JIJ, IRVDEPG) + PBU_INST(JIJ, IRSMLTG) + PBU_INST(JIJ, IRRACCSG) + &
            & PBU_INST(JIJ, IRSACCRG) + PBU_INST(JIJ, IRCRIMSG) + PBU_INST(JIJ, IRSRIMCG)
  ZRGSI_MR(JIJ) = PBU_INST(JIJ, IRRHONG_MR) + PBU_INST(JIJ, IRSRIMCG_MR)
ENDDO
!$mnh_end_do()

CALL ICE4_FAST_RG(CST, PARAMI, ICEP, ICED, D, ODSOFT, LDCOMPUTE, KRR, &
                 &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                 &ZDV, ZKA, ZCJ, PCIT, &
                 &ZLBDAR, ZLBDAS, ZLBDAG, &
                 &ZT, ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), &
                 &ZRGSI, ZRGSI_MR, &
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
  CALL ICE4_FAST_RH(CST, PARAMI, ICEP, ICED, D, ODSOFT, LDCOMPUTE, LLWETG, &
                   &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                   &ZDV, ZKA, ZCJ, &
                   &ZLBDAS, ZLBDAG, ZLBDAR, ZLBDAH, &
                   &ZT,  ZVART(:,IRV), ZVART(:,IRC), ZVART(:,IRR), ZVART(:,IRI), ZVART(:,IRS), ZVART(:,IRG), ZVART(:,IRH), &
                   &PBU_INST(:, IRCWETH), PBU_INST(:, IRIWETH), PBU_INST(:, IRSWETH), PBU_INST(:, IRGWETH), PBU_INST(:, IRRWETH), &
                   &PBU_INST(:, IRCDRYH), PBU_INST(:, IRIDRYH), PBU_INST(:, IRSDRYH), PBU_INST(:, IRRDRYH), &
                   &PBU_INST(:, IRGDRYH), PBU_INST(:, IRDRYHG), PBU_INST(:, IRHMLTR), &
                   &PRH_TEND)
ELSEIF (BUCONF%LBU_ENABLE .OR. OSAVE_MICRO) THEN

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

END IF
!
!-------------------------------------------------------------------------------
!
!
!*       7.     COMPUTES SPECIFIC SOURCES OF THE WARM AND COLD CLOUDY SPECIES
!               -------------------------------------------------------------
!
IF (.NOT.PARAMI%LOCND2) THEN
   CALL ICE4_FAST_RI(ICEP, ICED, D, &
                 &ODSOFT, LDCOMPUTE, &
                 &PRHODREF, PLVFACT, PLSFACT, &
                 &ZAI, ZCJ, PCIT, &
                 &PSSI, &
                 &ZVART(:,IRC), ZVART(:,IRI), &
                 &PBU_INST(:, IRCBERI))
ELSE
   PBU_INST(:, IRCBERI)=0.0
ENDIF

!
!*       8.     COMPUTES DEPOSITION OF ICE AND CONVERSION OF LARGE ICE CRYSTALS TO SNOW
!               -------------------------------------------------------------
!

IF (PARAMI%LOCND2) THEN
  CALL ICE4_FAST_RI_RS(CST, PARAMI, ICEP, ICED, D,&
                 &ODSOFT, LDCOMPUTE, &
                 &PRHODREF, PLSFACT, &
                 &ZAI, PCIT, ZESI, PPRES, &
                 &PSSI, PSSIO, PSSIU, PICLDFR, PIFR, &
                 &ZVART(:,IRV),ZVART(:,IRI),ZVART(:,IRS), ZT, &
                 &PBU_INST(:, IRILARS),PBU_INST(:, IRVDEPI))
ELSE
  PBU_INST(:, IRILARS)=0.
  PBU_INST(:, IRVDEPI)=0.
ENDIF
!
!-------------------------------------------------------------------------------
!
!
!*       8.     COMPUTES TOTAL TENDENCIES
!               -------------------------
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRVDEPG)*PLSFACT(JIJ)
  PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRCHONI)*(PLSFACT(JIJ)-PLVFACT(JIJ))
  PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRVDEPS)*PLSFACT(JIJ)
  PATH(JIJ) = PATH(JIJ) - PBU_INST(JIJ, IRREVAV)*PLVFACT(JIJ)
  PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRCRIMSS)*(PLSFACT(JIJ)-PLVFACT(JIJ))
  PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRCRIMSG)*(PLSFACT(JIJ)-PLVFACT(JIJ))
  PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRRACCSS)*(PLSFACT(JIJ)-PLVFACT(JIJ))
  PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRRACCSG)*(PLSFACT(JIJ)-PLVFACT(JIJ))
  PATH(JIJ) = PATH(JIJ) + (PBU_INST(JIJ, IRRCFRIG) - PBU_INST(JIJ, IRICFRR))*(PLSFACT(JIJ)-PLVFACT(JIJ))
  PATH(JIJ) = PATH(JIJ) + (PBU_INST(JIJ, IRCWETG) + PBU_INST(JIJ, IRRWETG))*(PLSFACT(JIJ)-PLVFACT(JIJ))
  PATH(JIJ) = PATH(JIJ) + (PBU_INST(JIJ, IRCDRYG)+PBU_INST(JIJ, IRRDRYG))*(PLSFACT(JIJ)-PLVFACT(JIJ))
  PATH(JIJ) = PATH(JIJ) - PBU_INST(JIJ, IRGMLTR)*(PLSFACT(JIJ)-PLVFACT(JIJ))
  IF (KRR==7) THEN
    PATH(JIJ) = PATH(JIJ) + (PBU_INST(JIJ, IRRWETH)+PBU_INST(JIJ, IRCWETH))*(PLSFACT(JIJ)-PLVFACT(JIJ))
    PATH(JIJ) = PATH(JIJ) + (PBU_INST(JIJ, IRCDRYH)+PBU_INST(JIJ, IRRDRYH))*(PLSFACT(JIJ)-PLVFACT(JIJ))
    PATH(JIJ) = PATH(JIJ) - PBU_INST(JIJ, IRHMLTR)*(PLSFACT(JIJ)-PLVFACT(JIJ))
  ENDIF
  PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRCBERI)*(PLSFACT(JIJ)-PLVFACT(JIJ))
  IF (PARAMI%LOCND2) PATH(JIJ) = PATH(JIJ) + PBU_INST(JIJ, IRVDEPI)*PLSFACT(JIJ)

  PA(JIJ, IRV) = PA(JIJ, IRV) - PBU_INST(JIJ, IRVDEPG)
  PA(JIJ, IRV) = PA(JIJ, IRV) - PBU_INST(JIJ, IRVDEPS)
  PA(JIJ, IRV) = PA(JIJ, IRV) + PBU_INST(JIJ, IRREVAV)
  IF (PARAMI%LOCND2) PA(JIJ, IRV) = PA(JIJ, IRV) - PBU_INST(JIJ, IRVDEPI)

  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCHONI)
  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCAUTR)
  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCACCR)
  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCRIMSS)
  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCRIMSG)
  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCMLTSR)
  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCWETG)
  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCDRYG)
  IF (KRR==7) THEN
    PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCWETH)
    PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCDRYH)
  ENDIF
  PA(JIJ, IRC) = PA(JIJ, IRC) - PBU_INST(JIJ, IRCBERI)

  PA(JIJ, IRR) = PA(JIJ, IRR) + PBU_INST(JIJ, IRCAUTR)
  PA(JIJ, IRR) = PA(JIJ, IRR) + PBU_INST(JIJ, IRCACCR)
  PA(JIJ, IRR) = PA(JIJ, IRR) - PBU_INST(JIJ, IRREVAV)
  PA(JIJ, IRR) = PA(JIJ, IRR) - PBU_INST(JIJ, IRRACCSS)
  PA(JIJ, IRR) = PA(JIJ, IRR) - PBU_INST(JIJ, IRRACCSG)
  PA(JIJ, IRR) = PA(JIJ, IRR) + PBU_INST(JIJ, IRCMLTSR)
  PA(JIJ, IRR) = PA(JIJ, IRR) - PBU_INST(JIJ, IRRCFRIG) + PBU_INST(JIJ, IRICFRR)
  PA(JIJ, IRR) = PA(JIJ, IRR) - PBU_INST(JIJ, IRRWETG)
  PA(JIJ, IRR) = PA(JIJ, IRR) - PBU_INST(JIJ, IRRDRYG)
  PA(JIJ, IRR) = PA(JIJ, IRR) + PBU_INST(JIJ, IRGMLTR)
  IF(KRR==7) THEN
    PA(JIJ, IRR) = PA(JIJ, IRR) - PBU_INST(JIJ, IRRWETH)
    PA(JIJ, IRR) = PA(JIJ, IRR) - PBU_INST(JIJ, IRRDRYH)
    PA(JIJ, IRR) = PA(JIJ, IRR) + PBU_INST(JIJ, IRHMLTR)
  ENDIF

  PA(JIJ, IRI) = PA(JIJ, IRI) + PBU_INST(JIJ, IRCHONI)
  PA(JIJ, IRI) = PA(JIJ, IRI) - PBU_INST(JIJ, IRIAGGS)
  PA(JIJ, IRI) = PA(JIJ, IRI) - PBU_INST(JIJ, IRIAUTS)
  PA(JIJ, IRI) = PA(JIJ, IRI) - PBU_INST(JIJ, IRICFRRG) - PBU_INST(JIJ, IRICFRR)
  PA(JIJ, IRI) = PA(JIJ, IRI) - PBU_INST(JIJ, IRIWETG)
  PA(JIJ, IRI) = PA(JIJ, IRI) - PBU_INST(JIJ, IRIDRYG)
  IF (KRR==7) THEN
    PA(JIJ, IRI) = PA(JIJ, IRI) - PBU_INST(JIJ, IRIWETH)
    PA(JIJ, IRI) = PA(JIJ, IRI) - PBU_INST(JIJ, IRIDRYH)
  ENDIF
  PA(JIJ, IRI) = PA(JIJ, IRI) + PBU_INST(JIJ, IRCBERI)
  IF (PARAMI%LOCND2) PA(JIJ, IRI) = PA(JIJ, IRI) + PBU_INST(JIJ, IRVDEPI)
  IF (PARAMI%LOCND2) PA(JIJ, IRI) = PA(JIJ, IRI) - PBU_INST(JIJ, IRILARS)

  PA(JIJ, IRS) = PA(JIJ, IRS) + PBU_INST(JIJ, IRVDEPS)
  IF (PARAMI%LOCND2) PA(JIJ, IRS) = PA(JIJ, IRS) + PBU_INST(JIJ, IRVDEPG) * ZDEPG_S(JIJ)
  PA(JIJ, IRS) = PA(JIJ, IRS) + PBU_INST(JIJ, IRIAGGS)
  PA(JIJ, IRS) = PA(JIJ, IRS) + PBU_INST(JIJ, IRIAUTS)
  PA(JIJ, IRS) = PA(JIJ, IRS) + PBU_INST(JIJ, IRCRIMSS)
  PA(JIJ, IRS) = PA(JIJ, IRS) - PBU_INST(JIJ, IRSRIMCG)
  PA(JIJ, IRS) = PA(JIJ, IRS) + PBU_INST(JIJ, IRRACCSS)
  PA(JIJ, IRS) = PA(JIJ, IRS) - PBU_INST(JIJ, IRSACCRG)
  PA(JIJ, IRS) = PA(JIJ, IRS) - PBU_INST(JIJ, IRSMLTG)
  PA(JIJ, IRS) = PA(JIJ, IRS) - PBU_INST(JIJ, IRSWETG)
  PA(JIJ, IRS) = PA(JIJ, IRS) - PBU_INST(JIJ, IRSDRYG)
  IF (KRR==7) THEN
    PA(JIJ, IRS) = PA(JIJ, IRS) - PBU_INST(JIJ, IRSWETH)
    PA(JIJ, IRS) = PA(JIJ, IRS) - PBU_INST(JIJ, IRSDRYH)
  ENDIF
  IF (PARAMI%LOCND2) PA(JIJ, IRS) = PA(JIJ, IRS) + PBU_INST(JIJ, IRILARS)

  PA(JIJ, IRG) = PA(JIJ, IRG) + PBU_INST(JIJ, IRVDEPG) * ZNODEPG_S(JIJ)
  PA(JIJ, IRG) = PA(JIJ, IRG) + PBU_INST(JIJ, IRCRIMSG)+PBU_INST(JIJ, IRSRIMCG)
  PA(JIJ, IRG) = PA(JIJ, IRG) + PBU_INST(JIJ, IRRACCSG)+PBU_INST(JIJ, IRSACCRG)
  PA(JIJ, IRG) = PA(JIJ, IRG) + PBU_INST(JIJ, IRSMLTG)
  PA(JIJ, IRG) = PA(JIJ, IRG) + PBU_INST(JIJ, IRICFRRG) + PBU_INST(JIJ, IRRCFRIG)
  PA(JIJ, IRG) = PA(JIJ, IRG) + PBU_INST(JIJ, IRCWETG) + PBU_INST(JIJ, IRIWETG) + PBU_INST(JIJ, IRSWETG) + PBU_INST(JIJ, IRRWETG)
  PA(JIJ, IRG) = PA(JIJ, IRG) - PBU_INST(JIJ, IRWETGH)
  PB(JIJ, IRG) = PB(JIJ, IRG) - PBU_INST(JIJ, IRWETGH_MR)
  PA(JIJ, IRG) = PA(JIJ, IRG) + PBU_INST(JIJ, IRCDRYG) + PBU_INST(JIJ, IRIDRYG) + PBU_INST(JIJ, IRSDRYG) + PBU_INST(JIJ, IRRDRYG)
  PA(JIJ, IRG) = PA(JIJ, IRG) - PBU_INST(JIJ, IRGMLTR)
  IF (KRR==7) THEN
    PA(JIJ, IRG) = PA(JIJ, IRG) - PBU_INST(JIJ, IRGWETH)
    PA(JIJ, IRG) = PA(JIJ, IRG) - PBU_INST(JIJ, IRGDRYH) + PBU_INST(JIJ, IRDRYHG)
  ENDIF

  IF (KRR==7) THEN
    PA(JIJ, IRH) = PA(JIJ, IRH) + PBU_INST(JIJ, IRWETGH)
    PB(JIJ, IRH) = PB(JIJ, IRH) + PBU_INST(JIJ, IRWETGH_MR)
    PA(JIJ, IRH) = PA(JIJ, IRH) + PBU_INST(JIJ, IRCWETH)+PBU_INST(JIJ, IRIWETH)+PBU_INST(JIJ, IRSWETH)+&
                              & PBU_INST(JIJ, IRGWETH)+PBU_INST(JIJ, IRRWETH)
    PA(JIJ, IRH) = PA(JIJ, IRH) + PBU_INST(JIJ, IRCDRYH)+PBU_INST(JIJ, IRIDRYH)+PBU_INST(JIJ, IRSDRYH)+&
                              & PBU_INST(JIJ, IRRDRYH)+PBU_INST(JIJ, IRGDRYH) - PBU_INST(JIJ, IRDRYHG)
    PA(JIJ, IRH) = PA(JIJ, IRH) - PBU_INST(JIJ, IRHMLTR)
  ENDIF
ENDDO
!$mnh_end_do()

!
IF (LHOOK) CALL DR_HOOK('ICE4_TENDENCIES', 1, ZHOOK_HANDLE)
!
CONTAINS
INCLUDE "ice4_nucleation.func.h"
END SUBROUTINE ICE4_TENDENCIES
END MODULE MODE_ICE4_TENDENCIES
