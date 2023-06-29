!MNH_LIC Copyright 1994-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_FAST_RS
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RS(CST, PARAMI, ICEP, ICED, KPROMA, KSIZE, LDSOFT, LDCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAR, PLBDAS, &
                       &PT,  PRVT, PRCT, PRRT, PRST, &
                       &PRIAGGS, &
                       &PRCRIMSS, PRCRIMSG, PRSRIMCG, &
                       &PRRACCSS, PRRACCSG, PRSACCRG, PRSMLTG, &
                       &PRCMLTSR, &
                       &PRS_TEND)
!!
!!**  PURPOSE
!!    -------
!!      Computes the fast rs processes
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!!     R. El Khatib 24-Aug-2021 Optimizations
!  J. Wurtz       03/2022: New snow characteristics with LSNOW_T
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
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
INTEGER,                      INTENT(IN)    :: KPROMA, KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KPROMA),   INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRIAGGS  ! r_i aggregation on r_s
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRCRIMSS ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRCRIMSG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRSRIMCG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRRACCSS ! Rain accretion onto the aggregates
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRRACCSG ! Rain accretion onto the aggregates
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRSACCRG ! Rain accretion onto the aggregates
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRSMLTG  ! Conversion-Melting of the aggregates
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRCMLTSR ! Cloud droplet collection onto aggregates by positive temperature
REAL, DIMENSION(KPROMA, 8),   INTENT(INOUT) :: PRS_TEND ! Individual tendencies
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCRIMS=1, IRCRIMSS=2, IRSRIMCG=3, IRRACCS=4, IRRACCSS=5, IRSACCRG=6, &
                    & IFREEZ1=7, IFREEZ2=8
LOGICAL, DIMENSION(KPROMA) :: GRIM, GACC
INTEGER :: IGRIM, IGACC
INTEGER, DIMENSION(KPROMA) :: IBUF1, IBUF2, IBUF3
REAL, DIMENSION(KPROMA) :: ZBUF1, ZBUF2, ZBUF3
REAL, DIMENSION(KPROMA) :: ZZW, ZZW1, ZZW2, ZZW3, ZFREEZ_RATE
INTEGER :: JL
REAL :: ZZW0D
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RS', 0, ZHOOK_HANDLE)
!
!
!-------------------------------------------------------------------------------
!
!
!*       5.0    maximum freezing rate
!
DO JL=1, KSIZE
  IF(PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRS_TEND(JL, IFREEZ1)=PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRS_TEND(JL, IFREEZ1)=MIN(PRS_TEND(JL, IFREEZ1), EXP(CST%XALPI-CST%XBETAI/PT(JL)-CST%XGAMI*ALOG(PT(JL)))) ! min(ev, es_i(T))
      ENDIF
      PRS_TEND(JL, IFREEZ1)=PKA(JL)*(CST%XTT-PT(JL)) +                              &
                           &(PDV(JL)*(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JL)-CST%XTT)) &
                           &*(CST%XESTT-PRS_TEND(JL, IFREEZ1))/(CST%XRV*PT(JL))           )
#ifdef REPRO48
      PRS_TEND(JL, IFREEZ1)=PRS_TEND(JL, IFREEZ1)* (ICEP%X0DEPS*       PLBDAS(JL)**ICEP%XEX0DEPS +     &
                           &                        ICEP%X1DEPS*PCJ(JL)*PLBDAS(JL)**ICEP%XEX1DEPS )/ &
                           &(PRHODREF(JL)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JL))))
#else
      PRS_TEND(JL, IFREEZ1)=PRS_TEND(JL, IFREEZ1)* PRST(JL) *(ICEP%X0DEPS*       PLBDAS(JL)**ICEP%XEX0DEPS +     &
                           &                        ICEP%X1DEPS*PCJ(JL)*PLBDAS(JL)**(ICED%XBS+ICEP%XEX1DEPS )*   &
         (1+0.5*(ICED%XFVELOS/PLBDAS(JL))**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEX1DEPS/ICED%XALPHAS))/ &
                           &(PRHODREF(JL)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JL))))
#endif
      PRS_TEND(JL, IFREEZ2)=(PRHODREF(JL)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PT(JL)))   ) / &
                           &(PRHODREF(JL)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JL))))
    ENDIF
    !We must agregate, at least, the cold species
    !And we are only interested by the freezing rate of liquid species
    ZFREEZ_RATE(JL)=MAX(0., MAX(0., PRS_TEND(JL, IFREEZ1) + &
                                    &PRS_TEND(JL, IFREEZ2) * PRIAGGS(JL)) - &
                            PRIAGGS(JL))
  ELSE
    PRS_TEND(JL, IFREEZ1)=0.
    PRS_TEND(JL, IFREEZ2)=0.
    ZFREEZ_RATE(JL)=0.
  ENDIF
ENDDO
!
!*       5.1    cloud droplet riming of the aggregates
!
DO JL=1, KSIZE
  IF (PRCT(JL)>ICED%XRTMIN(2) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL)) THEN
#ifdef REPRO48
    ZZW(JL) = PLBDAS(JL)
#else
    ZZW(JL) = (PLBDAS(JL)**ICED%XALPHAS + ICED%XFVELOS**ICED%XALPHAS)**(1./ICED%XALPHAS)
#endif
    GRIM(JL) = .TRUE.
  ELSE
    GRIM(JL) = .FALSE.
    PRS_TEND(JL, IRCRIMS)=0.
    PRS_TEND(JL, IRCRIMSS)=0.
    PRS_TEND(JL, IRSRIMCG)=0.
  ENDIF
ENDDO
!
! Collection of cloud droplets by snow: this rate is used for riming (T<0) and for conversion/melting (T>0)
IF(.NOT. LDSOFT) THEN
  CALL INTERP_MICRO_1D(KPROMA, KSIZE, ZZW, ICEP%NGAMINC, ICEP%XRIMINTP1, ICEP%XRIMINTP2, &
                           PARAMI%LPACK_INTERP, GRIM(:), IBUF1, IBUF2, ZBUF1, ZBUF2, &
                           IGRIM, &
                           ICEP%XGAMINC_RIM1(:), ZZW1(:), ICEP%XGAMINC_RIM2(:), ZZW2(:), ICEP%XGAMINC_RIM4(:), ZZW3(:))
  IF(IGRIM>0) THEN
    !
    !        5.1.4  riming of the small sized aggregates
    !
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE (GRIM(1:KSIZE))
#ifdef REPRO48
      PRS_TEND(1:KSIZE, IRCRIMSS) = ICEP%XCRIMSS * ZZW1(1:KSIZE) * PRCT(1:KSIZE) & ! RCRIMSS
                                      * PLBDAS(1:KSIZE)**ICEP%XEXCRIMSS &
                                      * PRHODREF(1:KSIZE)**(-ICED%XCEXVT)
#else
      PRS_TEND(1:KSIZE, IRCRIMSS) = ICEP%XCRIMSS * ZZW1(1:KSIZE) * PRCT(1:KSIZE) & ! RCRIMSS
                                      * PRST(1:KSIZE)*(1+(ICED%XFVELOS/PLBDAS(1:KSIZE))**ICED%XALPHAS) &
                                        **(-ICED%XNUS+ICEP%XEXCRIMSS/ICED%XALPHAS) &
                                      * PRHODREF(1:KSIZE)**(-ICED%XCEXVT+1.) &
                                      * (PLBDAS(1:KSIZE)) ** (ICEP%XEXCRIMSS+ICED%XBS)
#endif
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
    !
    !        5.1.6  riming-conversion of the large sized aggregates into graupeln
    !
    !
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GRIM(1:KSIZE))
#ifdef REPRO48
      PRS_TEND(1:KSIZE, IRCRIMS)=ICEP%XCRIMSG * PRCT(1:KSIZE)               & ! RCRIMS
                                   * PLBDAS(1:KSIZE)**ICEP%XEXCRIMSG  &
                                   * PRHODREF(1:KSIZE)**(-ICED%XCEXVT)
#else
      PRS_TEND(1:KSIZE, IRCRIMS)=ICEP%XCRIMSG * PRCT(1:KSIZE)               & ! RCRIMS
                                   * PRST(1:KSIZE)*(1+(ICED%XFVELOS/PLBDAS(1:KSIZE))**(ICED%XALPHAS)) &
                                     **(-ICED%XNUS+ICEP%XEXCRIMSG/ICED%XALPHAS) &
                                   * PRHODREF(1:KSIZE)**(-ICED%XCEXVT+1.) &
                                   * PLBDAS(1:KSIZE)**(ICED%XBS+ICEP%XEXCRIMSG)
#endif
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)

    IF(PARAMI%CSNOWRIMING=='M90 ')THEN
      !Murakami 1990
      !$mnh_expand_where(JL=1:KSIZE)
      WHERE(GRIM(1:KSIZE))
        ZZW(1:KSIZE) = PRS_TEND(1:KSIZE, IRCRIMS) - PRS_TEND(1:KSIZE, IRCRIMSS) ! RCRIMSG
#ifdef REPRO48
        PRS_TEND(1:KSIZE, IRSRIMCG)=ICEP%XSRIMCG * PLBDAS(1:KSIZE)**ICEP%XEXSRIMCG*(1.0-ZZW2(1:KSIZE))
#else
        PRS_TEND(1:KSIZE, IRSRIMCG)=ICEP%XSRIMCG * PRST(1:KSIZE)*PRHODREF(1:KSIZE) &
                                                 * PLBDAS(1:KSIZE)**(ICEP%XEXSRIMCG+ICED%XBS)*(1.0-ZZW2(1:KSIZE))
#endif
#ifdef REPRO48
        PRS_TEND(1:KSIZE, IRSRIMCG)=ZZW(1:KSIZE)*PRS_TEND(1:KSIZE, IRSRIMCG)/ &
                       MAX(1.E-20, &
                           ICEP%XSRIMCG3*ICEP%XSRIMCG2*PLBDAS(1:KSIZE)**ICEP%XEXSRIMCG2*(1.-ZZW3(1:KSIZE)) - &
                           ICEP%XSRIMCG3*PRS_TEND(1:KSIZE, IRSRIMCG))
#else
        PRS_TEND(1:KSIZE, IRSRIMCG)=ZZW(1:KSIZE)*PRS_TEND(1:KSIZE, IRSRIMCG)/ &
                       MAX(1.E-20, &
                           ICEP%XSRIMCG3*ICEP%XSRIMCG2*PRST(1:KSIZE)*PRHODREF(1:KSIZE) &
                                                      *PLBDAS(1:KSIZE)**ICEP%XEXSRIMCG2*(1.-ZZW3(1:KSIZE)) - &
                           ICEP%XSRIMCG3*PRS_TEND(1:KSIZE, IRSRIMCG))
#endif
      END WHERE
      !$mnh_end_expand_where(JL=1:KSIZE)
    ELSE
      PRS_TEND(:, IRSRIMCG)=0.
    END IF
  ENDIF
ENDIF
!
DO JL=1, KSIZE
  ! More restrictive RIM mask to be used for riming by negative temperature only
  IF(GRIM(JL) .AND. PT(JL)<CST%XTT) THEN
    PRCRIMSS(JL)=MIN(ZFREEZ_RATE(JL), PRS_TEND(JL, IRCRIMSS))
    ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRCRIMSS(JL))
    ZZW0D = MIN(1., ZFREEZ_RATE(JL) / MAX(1.E-20, PRS_TEND(JL, IRCRIMS) - PRCRIMSS(JL))) ! proportion we are able to freeze
    PRCRIMSG(JL) = ZZW0D * MAX(0., PRS_TEND(JL, IRCRIMS) - PRCRIMSS(JL)) ! RCRIMSG
    ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRCRIMSG(JL))
    PRSRIMCG(JL) = ZZW0D * PRS_TEND(JL, IRSRIMCG)

    PRSRIMCG(JL) = PRSRIMCG(JL) * MAX(0., -SIGN(1., -PRCRIMSG(JL)))
    PRCRIMSG(JL)=MAX(0., PRCRIMSG(JL))
  ELSE
    PRCRIMSS(JL)=0.
    PRCRIMSG(JL)=0.
    PRSRIMCG(JL)=0.
  ENDIF
ENDDO
!
!*       5.2    rain accretion onto the aggregates
!
DO JL = 1, KSIZE
  IF (PRRT(JL)>ICED%XRTMIN(3) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL)) THEN
    GACC(JL) = .TRUE.
  ELSE
    GACC(JL) = .FALSE.
    PRS_TEND(JL, IRRACCS)=0.
    PRS_TEND(JL, IRRACCSS)=0.
    PRS_TEND(JL, IRSACCRG)=0.
  END IF
ENDDO
IF(.NOT. LDSOFT) THEN
  PRS_TEND(:, IRRACCS)=0.
  PRS_TEND(:, IRRACCSS)=0.
  PRS_TEND(:, IRSACCRG)=0.
  CALL INTERP_MICRO_2D(KPROMA, KSIZE, PLBDAS, PLBDAR, ICEP%NACCLBDAS, ICEP%NACCLBDAR, &
                      &ICEP%XACCINTP1S, ICEP%XACCINTP2S, ICEP%XACCINTP1R, ICEP%XACCINTP2R,&
                      &PARAMI%LPACK_INTERP, GACC(:), IBUF1(:), IBUF2(:), IBUF3(:), ZBUF1(:), ZBUF2(:), ZBUF3(:), &
                      &IGACC, &
                      &ICEP%XKER_RACCSS(:,:), ZZW1(:), ICEP%XKER_RACCS(:,:), ZZW2(:), ICEP%XKER_SACCRG(:,:), ZZW3(:))
  IF(IGACC>0)THEN
    !        5.2.4  raindrop accretion on the small sized aggregates
    !
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GACC(1:KSIZE))
#ifdef REPRO48
      ZZW(1:KSIZE) =                                                        & !! coef of RRACCS
            ICEP%XFRACCSS*( PLBDAS(1:KSIZE)**ICED%XCXS )*( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) ) &
       *( ICEP%XLBRACCS1/((PLBDAS(1:KSIZE)**2)               ) +                  &
          ICEP%XLBRACCS2/( PLBDAS(1:KSIZE)    * PLBDAR(1:KSIZE)    ) +                  &
          ICEP%XLBRACCS3/(               (PLBDAR(1:KSIZE)**2)) )/PLBDAR(1:KSIZE)**4
#else
      ZZW(1:KSIZE) =                                                        & !! coef of RRACCS
            ICEP%XFRACCSS*( PRST(1:KSIZE)*PLBDAS(1:KSIZE)**ICED%XBS )*( PRHODREF(1:KSIZE)**(-ICED%XCEXVT) ) &
       *( ICEP%XLBRACCS1/((PLBDAS(1:KSIZE)**2)               ) +                  &
          ICEP%XLBRACCS2/( PLBDAS(1:KSIZE)    * PLBDAR(1:KSIZE)    ) +                  &
          ICEP%XLBRACCS3/(               (PLBDAR(1:KSIZE)**2)) )/PLBDAR(1:KSIZE)**4
#endif
      PRS_TEND(1:KSIZE, IRRACCSS) =ZZW1(1:KSIZE)*ZZW(1:KSIZE)
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
    !
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GACC(1:KSIZE))
      PRS_TEND(1:KSIZE, IRRACCS) = ZZW2(1:KSIZE)*ZZW(1:KSIZE)
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
    !
    !        5.2.6  raindrop accretion-conversion of the large sized aggregates
    !               into graupeln
    !
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GACC(1:KSIZE))
#ifdef REPRO48
      PRS_TEND(1:KSIZE, IRSACCRG) = ICEP%XFSACCRG*ZZW3(1:KSIZE)*                    & ! RSACCRG
          ( PLBDAS(1:KSIZE)**(ICED%XCXS-ICED%XBS) )*( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) ) &
         *( ICEP%XLBSACCR1/((PLBDAR(1:KSIZE)**2)               ) +           &
            ICEP%XLBSACCR2/( PLBDAR(1:KSIZE)    * PLBDAS(1:KSIZE)    ) +           &
            ICEP%XLBSACCR3/(               (PLBDAS(1:KSIZE)**2)) )/PLBDAR(1:KSIZE)
#else
      PRS_TEND(1:KSIZE, IRSACCRG) = ICEP%XFSACCRG*ZZW3(1:KSIZE)*                    & ! RSACCRG
          ( PRST(1:KSIZE))*( PRHODREF(1:KSIZE)**(-ICED%XCEXVT) ) &
         *( ICEP%XLBSACCR1/((PLBDAR(1:KSIZE)**2)               ) +           &
            ICEP%XLBSACCR2/( PLBDAR(1:KSIZE)    * PLBDAS(1:KSIZE)    ) +           &
            ICEP%XLBSACCR3/(               (PLBDAS(1:KSIZE)**2)) )/PLBDAR(1:KSIZE)
#endif
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
  ENDIF
ENDIF
!
DO JL=1, KSIZE
  ! More restrictive ACC mask to be used for accretion by negative temperature only
  IF(GACC(JL) .AND. PT(JL)<CST%XTT) THEN
    PRRACCSS(JL)=MIN(ZFREEZ_RATE(JL), PRS_TEND(JL, IRRACCSS))
    ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRRACCSS(JL))
    ZZW(JL) = MIN(1., ZFREEZ_RATE(JL) / MAX(1.E-20, PRS_TEND(JL, IRRACCS)-PRRACCSS(JL))) ! proportion we are able to freeze
    PRRACCSG(JL)=ZZW(JL) * MAX(0., PRS_TEND(JL, IRRACCS)-PRRACCSS(JL))
    ZFREEZ_RATE(JL) = MAX(0., ZFREEZ_RATE(JL)-PRRACCSG(JL))
    PRSACCRG(JL)=ZZW(JL) * PRS_TEND(JL, IRSACCRG)

    PRSACCRG(JL) = PRSACCRG(JL) * MAX(0., -SIGN(1., -PRRACCSG(JL)))
    PRRACCSG(JL)=MAX(0., PRRACCSG(JL))
  ELSE
    PRRACCSS(JL)=0.
    PRRACCSG(JL)=0.
    PRSACCRG(JL)=0.
  ENDIF
ENDDO
!
!
!*       5.3    Conversion-Melting of the aggregates
!
DO JL=1, KSIZE
  IF(PRST(JL)>ICED%XRTMIN(5) .AND. PT(JL)>CST%XTT .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRSMLTG(JL)=PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRSMLTG(JL)=MIN(PRSMLTG(JL), EXP(CST%XALPW-CST%XBETAW/PT(JL)-CST%XGAMW*ALOG(PT(JL)))) ! min(ev, es_w(T))
      ENDIF
      PRSMLTG(JL)= PKA(JL)*(CST%XTT-PT(JL)) +                                 &
                  &(PDV(JL)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(JL) - CST%XTT )) &
                  & *(CST%XESTT-PRSMLTG(JL))/(CST%XRV*PT(JL))             )
      !
      ! compute RSMLT
      !
#ifdef REPRO48
      PRSMLTG(JL)  = ICEP%XFSCVMG*MAX(0., (-PRSMLTG(JL) * &
                 (ICEP%X0DEPS*       PLBDAS(JL)**ICEP%XEX0DEPS +     &
                 ICEP%X1DEPS*PCJ(JL)*PLBDAS(JL)**ICEP%XEX1DEPS)    &
                 -(PRS_TEND(JL, IRCRIMS) + PRS_TEND(JL, IRRACCS)) *       &
                 (PRHODREF(JL)*CST%XCL*(CST%XTT-PT(JL))) &
                 ) / (PRHODREF(JL)*CST%XLMTT))
#else
      PRSMLTG(JL)  = ICEP%XFSCVMG*MAX(0., (-PRSMLTG(JL) * &
                 PRST(JL)*PRHODREF(JL) *    &
                 (ICEP%X0DEPS*       PLBDAS(JL)**(ICED%XBS+ICEP%XEX0DEPS) + &
                 ICEP%X1DEPS*PCJ(JL)*(1+0.5*(ICED%XFVELOS/PLBDAS(JL))**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEX1DEPS/ICED%XALPHAS) &
                 *PLBDAS(JL)**(ICED%XBS+ICEP%XEX1DEPS)) &
                 -(PRS_TEND(JL, IRCRIMS) + PRS_TEND(JL, IRRACCS)) *       &
                 (PRHODREF(JL)*CST%XCL*(CST%XTT-PT(JL))) &
                 ) / (PRHODREF(JL)*CST%XLMTT))
#endif
      !
      ! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
      ! because the graupeln produced by this process are still icy!!!
      !
      ! When T < XTT, rc is collected by snow (riming) to produce snow and graupel
      ! When T > XTT, if riming was still enabled, rc would produce snow and graupel with snow becomming graupel (conversion/melting) and graupel becomming rain (melting)
      ! To insure consistency when crossing T=XTT, rc collected with T>XTT must be transformed in rain.
      ! rc cannot produce iced species with a positive temperature but is still collected with a good efficiency by snow
      PRCMLTSR(JL) = PRS_TEND(JL, IRCRIMS) ! both species are liquid, no heat is exchanged
    ENDIF
  ELSE
    PRSMLTG(JL)=0.
    PRCMLTSR(JL)=0.
  ENDIF
ENDDO

IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RS', 1, ZHOOK_HANDLE)
!
CONTAINS
!
INCLUDE "interp_micro.func.h"
!
END SUBROUTINE ICE4_FAST_RS
END MODULE MODE_ICE4_FAST_RS
