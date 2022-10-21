!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
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
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
INTEGER, DIMENSION(KPROMA) :: I1
REAL, DIMENSION(KPROMA) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(KPROMA) :: IVEC1, IVEC2
REAL, DIMENSION(KPROMA) :: ZZW, ZZW2, ZZW6, ZFREEZ_RATE
INTEGER :: JJ, JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
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
      PRS_TEND(JL, IFREEZ1)=PRS_TEND(JL, IFREEZ1)* (ICEP%X0DEPS*       PLBDAS(JL)**ICEP%XEX0DEPS +     &
                           &                        ICEP%X1DEPS*PCJ(JL)*PLBDAS(JL)**ICEP%XEX1DEPS )/ &
                           &(PRHODREF(JL)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JL))))
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
IGRIM = 0
DO JL=1, KSIZE
  IF (PRCT(JL)>ICED%XRTMIN(2) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL)) THEN
    IGRIM = IGRIM + 1
    I1(IGRIM) = JL
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
  IF(IGRIM>0) THEN
    !
    !        5.1.1  select the PLBDAS
    !
    DO JJ = 1, IGRIM
      ZVEC1(JJ) = PLBDAS(I1(JJ))
    END DO
    !
    !        5.1.2  find the next lower indice for the PLBDAS in the geometrical
    !               set of Lbda_s used to tabulate some moments of the incomplete
    !               gamma function
    !
      DO JJ=1,IGRIM 
        ZVEC2(JJ) = MAX( 1.00001, MIN( REAL(ICEP%NGAMINC)-0.00001,           &
        ICEP%XRIMINTP1 * LOG( ZVEC1(JJ) ) + ICEP%XRIMINTP2 ) )
        IVEC2(JJ) = INT( ZVEC2(JJ) )
        ZVEC2(JJ) = ZVEC2(JJ) - REAL( IVEC2(JJ) )
    !
    !        5.1.3  perform the linear interpolation of the normalized
    !               "2+XDS"-moment of the incomplete gamma function
    !
        ZVEC1(JJ) =   ICEP%XGAMINC_RIM1( IVEC2(JJ)+1 )* ZVEC2(JJ)      &
        - ICEP%XGAMINC_RIM1( IVEC2(JJ)   )*(ZVEC2(JJ) - 1.0)
      ENDDO
    ZZW(:) = 0.
    DO JJ = 1, IGRIM
      ZZW(I1(JJ)) = ZVEC1(JJ)
    END DO
    !
    !        5.1.4  riming of the small sized aggregates
    !
      DO JL=1,KSIZE 
        IF (GRIM(JL))THEN
          PRS_TEND(JL,IRCRIMSS) = ICEP%XCRIMSS * ZZW(JL) * PRCT(JL) & ! RCRIMSS
          *   PLBDAS(JL)**ICEP%XEXCRIMSS &
          * PRHODREF(JL)**(-ICED%XCEXVT)
        ENDIF
      ENDDO
    !
    !        5.1.5  perform the linear interpolation of the normalized
    !               "XBS"-moment of the incomplete gamma function (XGAMINC_RIM2) and
    !               "XBG"-moment of the incomplete gamma function (XGAMINC_RIM4)
    !
      DO JJ=1,IGRIM 
        ZVEC1(JJ) =  ICEP%XGAMINC_RIM2( IVEC2(JJ)+1 )* ZVEC2(JJ)      &
        - ICEP%XGAMINC_RIM2( IVEC2(JJ)   )*(ZVEC2(JJ) - 1.0)
      ENDDO
    ZZW(:) = 0.
    DO JJ = 1, IGRIM
      ZZW(I1(JJ)) = ZVEC1(JJ)
    END DO

      DO JJ=1,IGRIM 
        ZVEC1(JJ) =  ICEP%XGAMINC_RIM4( IVEC2(JJ)+1 )* ZVEC2(JJ)      &
        - ICEP%XGAMINC_RIM4( IVEC2(JJ)   )*(ZVEC2(JJ) - 1.0)
      ENDDO
    ZZW2(:) = 0.
    DO JJ = 1, IGRIM
      ZZW2(I1(JJ)) = ZVEC1(JJ)
    END DO
    !
    !        5.1.6  riming-conversion of the large sized aggregates into graupeln
    !
    !
      DO JL=1,KSIZE 
        IF(GRIM(JL))THEN
          PRS_TEND(JL,IRCRIMS)=ICEP%XCRIMSG * PRCT(JL)               & ! RCRIMS
          * PLBDAS(JL)**ICEP%XEXCRIMSG  &
          * PRHODREF(JL)**(-ICED%XCEXVT)
        ENDIF
      ENDDO

    IF(PARAMI%CSNOWRIMING=='M90 ')THEN
      !Murakami 1990
        DO JL=1,KSIZE 
          IF(GRIM(JL))THEN
            ZZW6(JL) = PRS_TEND(JL,IRCRIMS) - PRS_TEND(JL,IRCRIMSS) ! RCRIMSG
            PRS_TEND(JL,IRSRIMCG)=ICEP%XSRIMCG * PLBDAS(JL)**ICEP%XEXSRIMCG*(1.0-ZZW(JL))
            PRS_TEND(JL,IRSRIMCG)=ZZW6(JL)*PRS_TEND(JL,IRSRIMCG)/ &
            MAX(1.E-20, &
            ICEP%XSRIMCG3*ICEP%XSRIMCG2*PLBDAS(JL)**ICEP%XEXSRIMCG2*(1.-ZZW2(JL)) - &
            ICEP%XSRIMCG3*PRS_TEND(JL,IRSRIMCG))
          ENDIF
        ENDDO
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
    ZZW(JL) = MIN(1., ZFREEZ_RATE(JL) / MAX(1.E-20, PRS_TEND(JL, IRCRIMS) - PRCRIMSS(JL))) ! proportion we are able to freeze
    PRCRIMSG(JL) = ZZW(JL) * MAX(0., PRS_TEND(JL, IRCRIMS) - PRCRIMSS(JL)) ! RCRIMSG
    ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRCRIMSG(JL))
    PRSRIMCG(JL) = ZZW(JL) * PRS_TEND(JL, IRSRIMCG)

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
IGACC = 0
DO JL = 1, KSIZE
  IF (PRRT(JL)>ICED%XRTMIN(3) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL)) THEN
    IGACC = IGACC + 1
    I1(IGACC) = JL
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
  IF(IGACC>0)THEN
    !
    !
    !        5.2.1  select the (PLBDAS,PLBDAR) couplet
    !
    DO JJ = 1, IGACC
      ZVEC1(JJ) = PLBDAS(I1(JJ))
      ZVEC2(JJ) = PLBDAR(I1(JJ))
    ENDDO
    !
    !        5.2.2  find the next lower indice for the PLBDAS and for the PLBDAR
    !               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
    !               tabulate the RACCSS-kernel
    !
      DO JJ=1,IGACC 
        ZVEC1(JJ) = MAX( 1.00001, MIN( REAL(ICEP%NACCLBDAS)-0.00001,           &
        ICEP%XACCINTP1S * LOG( ZVEC1(JJ) ) + ICEP%XACCINTP2S ) )
        IVEC1(JJ) = INT( ZVEC1(JJ) )
        ZVEC1(JJ) = ZVEC1(JJ) - REAL( IVEC1(JJ) )
    !
        ZVEC2(JJ) = MAX( 1.00001, MIN( REAL(ICEP%NACCLBDAR)-0.00001,           &
        ICEP%XACCINTP1R * LOG( ZVEC2(JJ) ) + ICEP%XACCINTP2R ) )
        IVEC2(JJ) = INT( ZVEC2(JJ) )
        ZVEC2(JJ) = ZVEC2(JJ) - REAL( IVEC2(JJ) )
      ENDDO
    !
    !        5.2.3  perform the bilinear interpolation of the normalized
    !               RACCSS-kernel
    !
    DO JJ = 1, IGACC
      ZVEC3(JJ) =  (  ICEP%XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                 - (  ICEP%XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = 0.
    DO JJ = 1, IGACC
      ZZW(I1(JJ)) = ZVEC3(JJ)
    END DO
    !
    !        5.2.4  raindrop accretion on the small sized aggregates
    !
      DO JL=1,KSIZE 
        IF(GACC(JL))THEN
          ZZW6(JL) =                                                        & !! coef of RRACCS
          ICEP%XFRACCSS*( PLBDAS(JL)**ICED%XCXS )*( PRHODREF(JL)**(-ICED%XCEXVT-1.) ) &
          *( ICEP%XLBRACCS1/((PLBDAS(JL)**2)               ) +                  &
          ICEP%XLBRACCS2/( PLBDAS(JL)    * PLBDAR(JL)    ) +                  &
          ICEP%XLBRACCS3/(               (PLBDAR(JL)**2)) )/PLBDAR(JL)**4
          PRS_TEND(JL,IRRACCSS) =ZZW(JL)*ZZW6(JL)
        ENDIF
      ENDDO
    !
    !        5.2.4b perform the bilinear interpolation of the normalized
    !               RACCS-kernel
    !
    DO JJ = 1, IGACC
      ZVEC3(JJ) =  (   ICEP%XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  ICEP%XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                 - (   ICEP%XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  ICEP%XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = 0.
    DO JJ = 1, IGACC
      ZZW(I1(JJ)) = ZVEC3(JJ)
    END DO
      DO JL=1,KSIZE 
        IF(GACC(JL))THEN
          PRS_TEND(JL,IRRACCS) = ZZW(JL)*ZZW6(JL)
        ENDIF
      ENDDO
    !        5.2.5  perform the bilinear interpolation of the normalized
    !               SACCRG-kernel
    !
    DO JJ = 1, IGACC
        ZVEC3(JJ) =  (  ICEP%XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                      - ICEP%XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                            * ZVEC2(JJ) &
                   - (  ICEP%XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                      - ICEP%XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                            * (ZVEC2(JJ) - 1.0)
    END DO
    ZZW(:) = 0.
    DO JJ = 1, IGACC
      ZZW(I1(JJ)) = ZVEC3(JJ)
    END DO
    !
    !        5.2.6  raindrop accretion-conversion of the large sized aggregates
    !               into graupeln
    !
      DO JL=1,KSIZE 
        IF(GACC(JL))THEN
          PRS_TEND(JL,IRSACCRG) = ICEP%XFSACCRG*ZZW(JL)*                    & ! RSACCRG
          ( PLBDAS(JL)**(ICED%XCXS-ICED%XBS) )*( PRHODREF(JL)**(-ICED%XCEXVT-1.) ) &
          *( ICEP%XLBSACCR1/((PLBDAR(JL)**2)               ) +           &
          ICEP%XLBSACCR2/( PLBDAR(JL)    * PLBDAS(JL)    ) +           &
          ICEP%XLBSACCR3/(               (PLBDAS(JL)**2)) )/PLBDAR(JL)
        ENDIF
      ENDDO
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
      PRSMLTG(JL)  = ICEP%XFSCVMG*MAX(0., (-PRSMLTG(JL) * (ICEP%X0DEPS*       PLBDAS(JL)**ICEP%XEX0DEPS +     &
                                                           ICEP%X1DEPS*PCJ(JL)*PLBDAS(JL)**ICEP%XEX1DEPS)    &
                                           -(PRS_TEND(JL, IRCRIMS) + PRS_TEND(JL, IRRACCS)) *       &
                                            (PRHODREF(JL)*CST%XCL*(CST%XTT-PT(JL))) &
                                          ) / (PRHODREF(JL)*CST%XLMTT))
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
END SUBROUTINE ICE4_FAST_RS
END MODULE MODE_ICE4_FAST_RS
