!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ICE4_FAST_RS
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RS(CST, PARAMI, ICEP, ICED, KPROMA,KSIZE, LDSOFT, PCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAR, PLBDAS, &
                       &PT,  PRVT, PRCT, PRRT, PRST, &
                       &PRIAGGS, &
                       &PRCRIMSS, PRCRIMSG, PRSRIMCG, &
                       &PRRACCSS, PRRACCSG, PRSACCRG, PRSMLTG, &
                       &PRCMLTSR, &
                       &PRS_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RS, PA_RG)
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
INTEGER,                      INTENT(IN)    :: KPROMA,KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIAGGS  ! r_i aggregation on r_s
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCRIMSS ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCRIMSG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSRIMCG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRACCSS ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRACCSG ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSACCRG ! Rain accretion onto the aggregates
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSMLTG  ! Conversion-Melting of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCMLTSR ! Cloud droplet collection onto aggregates by positive temperature
REAL, DIMENSION(KPROMA, 8),   INTENT(INOUT) :: PRS_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCRIMS=1, IRCRIMSS=2, IRSRIMCG=3, IRRACCS=4, IRRACCSS=5, IRSACCRG=6, &
                    & IFREEZ1=7, IFREEZ2=8
REAL, DIMENSION(KSIZE) :: ZRIM, ZACC, ZMASK
LOGICAL, DIMENSION(KSIZE) :: GRIM, GACC
INTEGER :: IGRIM, IGACC
INTEGER, DIMENSION(KSIZE) :: I1
REAL, DIMENSION(KSIZE) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(KSIZE) :: IVEC1, IVEC2
REAL, DIMENSION(KSIZE) :: ZZW, ZZW2, ZZW6, ZFREEZ_RATE
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
  ZMASK(JL)=MAX(0., -SIGN(1., ICED%XRTMIN(5)-PRST(JL))) * & ! WHERE(PRST(:)>XRTMIN(5))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRS_TEND(JL, IFREEZ1)=ZMASK(JL) * PRS_TEND(JL, IFREEZ1)
    PRS_TEND(JL, IFREEZ2)=ZMASK(JL) * PRS_TEND(JL, IFREEZ2)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRS_TEND(JL, IFREEZ1)=ZMASK(JL) * PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(PARAMI%LEVLIMIT) THEN
    WHERE(ZMASK(1:KSIZE)==1.)
      PRS_TEND(1:KSIZE, IFREEZ1)=MIN(PRS_TEND(1:KSIZE, IFREEZ1), EXP(CST%XALPI-CST%XBETAI/PT(1:KSIZE)-CST%XGAMI*ALOG(PT(1:KSIZE)))) ! min(ev, es_i(T))
    END WHERE
  ENDIF
  PRS_TEND(:, IFREEZ2)=0.
  WHERE(ZMASK(1:KSIZE)==1.)
    PRS_TEND(1:KSIZE, IFREEZ1)=PKA(1:KSIZE)*(CST%XTT-PT(1:KSIZE)) +                              &
             (PDV(1:KSIZE)*(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(1:KSIZE)-CST%XTT)) &
                           *(CST%XESTT-PRS_TEND(1:KSIZE, IFREEZ1))/(CST%XRV*PT(1:KSIZE))           )
    PRS_TEND(1:KSIZE, IFREEZ1)=PRS_TEND(1:KSIZE, IFREEZ1)* ( ICEP%X0DEPS*       PLBDAS(1:KSIZE)**ICEP%XEX0DEPS +     &
                           ICEP%X1DEPS*PCJ(1:KSIZE)*PLBDAS(1:KSIZE)**ICEP%XEX1DEPS )/ &
                          ( PRHODREF(1:KSIZE)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(1:KSIZE))) )
    PRS_TEND(1:KSIZE, IFREEZ2)=(PRHODREF(1:KSIZE)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PT(1:KSIZE)))   ) / &
                          ( PRHODREF(1:KSIZE)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(1:KSIZE))) )
  END WHERE
ENDIF
DO JL=1, KSIZE
  !We must agregate, at least, the cold species
  !And we are only interested by the freezing rate of liquid species
  ZFREEZ_RATE(JL)=ZMASK(JL) * MAX(0., MAX(0., PRS_TEND(JL, IFREEZ1) + &
                                              &PRS_TEND(JL, IFREEZ2) * PRIAGGS(JL)) - &
                                      PRIAGGS(JL))
ENDDO
!
!*       5.1    cloud droplet riming of the aggregates
!
IGRIM = 0
DO JL=1, KSIZE
  ZRIM(JL)=MAX(0., -SIGN(1., ICED%XRTMIN(2)-PRCT(JL))) * & !WHERE(PRCT(:)>XRTMIN(2))
          &MAX(0., -SIGN(1., ICED%XRTMIN(5)-PRST(JL))) * & !WHERE(PRST(:)>XRTMIN(5))
          &PCOMPUTE(JL)
  IF (ZRIM(JL)>0) THEN
    IGRIM = IGRIM + 1
    I1(IGRIM) = JL
    GRIM(JL) = .TRUE.
  ELSE
    GRIM(JL) = .FALSE.
  ENDIF
ENDDO
!
! Collection of cloud droplets by snow: this rate is used for riming (T<0) and for conversion/melting (T>0)
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRS_TEND(JL, IRCRIMS)=ZRIM(JL) * PRS_TEND(JL, IRCRIMS)
    PRS_TEND(JL, IRCRIMSS)=ZRIM(JL) * PRS_TEND(JL, IRCRIMSS)
    PRS_TEND(JL, IRSRIMCG)=ZRIM(JL) * PRS_TEND(JL, IRSRIMCG)
  ENDDO
ELSE
  PRS_TEND(:, IRCRIMS)=0.
  PRS_TEND(:, IRCRIMSS)=0.
  PRS_TEND(:, IRSRIMCG)=0.
  !
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
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN( REAL(ICEP%NGAMINC)-0.00001,           &
                          ICEP%XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + ICEP%XRIMINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
    !
    !        5.1.3  perform the linear interpolation of the normalized
    !               "2+XDS"-moment of the incomplete gamma function
    !
    ZVEC1(1:IGRIM) =   ICEP%XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - ICEP%XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = 0.
    DO JJ = 1, IGRIM
      ZZW(I1(JJ)) = ZVEC1(JJ)
    END DO
    !
    !        5.1.4  riming of the small sized aggregates
    !
    WHERE (GRIM(1:KSIZE))
      PRS_TEND(1:KSIZE, IRCRIMSS) = ICEP%XCRIMSS * ZZW(1:KSIZE) * PRCT(1:KSIZE)                & ! RCRIMSS
                                      *   PLBDAS(1:KSIZE)**ICEP%XEXCRIMSS &
                                      * PRHODREF(1:KSIZE)**(-ICED%XCEXVT)
    END WHERE
    !
    !        5.1.5  perform the linear interpolation of the normalized
    !               "XBS"-moment of the incomplete gamma function (XGAMINC_RIM2) and
    !               "XBG"-moment of the incomplete gamma function (XGAMINC_RIM4)
    !
    ZVEC1(1:IGRIM) =  ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = 0.
    DO JJ = 1, IGRIM
      ZZW(I1(JJ)) = ZVEC1(JJ)
    END DO

    ZVEC1(1:IGRIM) =  ICEP%XGAMINC_RIM4( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - ICEP%XGAMINC_RIM4( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW2(:) = 0.
    DO JJ = 1, IGRIM
      ZZW2(I1(JJ)) = ZVEC1(JJ)
    END DO
    !
    !        5.1.6  riming-conversion of the large sized aggregates into graupeln
    !
    !
    WHERE(GRIM(1:KSIZE))
      PRS_TEND(1:KSIZE, IRCRIMS)=ICEP%XCRIMSG * PRCT(1:KSIZE)               & ! RCRIMS
                                   * PLBDAS(1:KSIZE)**ICEP%XEXCRIMSG  &
                                   * PRHODREF(1:KSIZE)**(-ICED%XCEXVT)
      ZZW6(1:KSIZE) = PRS_TEND(1:KSIZE, IRCRIMS) - PRS_TEND(1:KSIZE, IRCRIMSS) ! RCRIMSG
    END WHERE

    IF(PARAMI%CSNOWRIMING=='M90 ')THEN
      !Murakami 1990
      WHERE(GRIM(1:KSIZE))
        PRS_TEND(1:KSIZE, IRSRIMCG)=ICEP%XSRIMCG * PLBDAS(1:KSIZE)**ICEP%XEXSRIMCG*(1.0-ZZW(1:KSIZE))
        PRS_TEND(1:KSIZE, IRSRIMCG)=ZZW6(1:KSIZE)*PRS_TEND(1:KSIZE, IRSRIMCG)/ &
                       MAX(1.E-20, &
                           ICEP%XSRIMCG3*ICEP%XSRIMCG2*PLBDAS(1:KSIZE)**ICEP%XEXSRIMCG2*(1.-ZZW2(1:KSIZE)) - &
                           ICEP%XSRIMCG3*PRS_TEND(1:KSIZE, IRSRIMCG))
      END WHERE
    ELSE
      PRS_TEND(:, IRSRIMCG)=0.
    END IF
  ENDIF
ENDIF
!
DO JL=1, KSIZE
  ! More restrictive RIM mask to be used for riming by negative temperature only
  ZRIM(JL)=ZRIM(JL) * &
          &MAX(0., -SIGN(1., PT(JL)-CST%XTT)) ! WHERE(PT(:)<XTT)
  PRCRIMSS(JL)=ZRIM(JL)*MIN(ZFREEZ_RATE(JL), PRS_TEND(JL, IRCRIMSS))
  ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRCRIMSS(JL))
  ZZW(JL) = MIN(1., ZFREEZ_RATE(JL) / MAX(1.E-20, PRS_TEND(JL, IRCRIMS) - PRCRIMSS(JL))) ! proportion we are able to freeze
  PRCRIMSG(JL) = ZRIM(JL) * ZZW(JL) * MAX(0., PRS_TEND(JL, IRCRIMS) - PRCRIMSS(JL)) ! RCRIMSG
  ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRCRIMSG(JL))
  PRSRIMCG(JL) = ZRIM(JL) * ZZW(JL) * PRS_TEND(JL, IRSRIMCG)

  PRSRIMCG(JL) = PRSRIMCG(JL) * MAX(0., -SIGN(1., -PRCRIMSG(JL)))
  PRCRIMSG(JL)=MAX(0., PRCRIMSG(JL))

ENDDO
!
!*       5.2    rain accretion onto the aggregates
!
IGACC = 0
DO JJ = 1, KSIZE
  ZACC(JJ)=MAX(0., -SIGN(1., ICED%XRTMIN(3)-PRRT(JJ))) * & !WHERE(PRRT(:)>XRTMIN(3))
          &MAX(0., -SIGN(1., ICED%XRTMIN(5)-PRST(JJ))) * & !WHERE(PRST(:)>XRTMIN(5))
          &PCOMPUTE(JJ)
  IF (ZACC(JJ)>0) THEN
    IGACC = IGACC + 1
    I1(IGACC) = JJ
    GACC(JJ) = .TRUE.
  ELSE
    GACC(JJ) = .FALSE.
  END IF
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRS_TEND(JL, IRRACCS)=ZACC(JL) * PRS_TEND(JL, IRRACCS)
    PRS_TEND(JL, IRRACCSS)=ZACC(JL) * PRS_TEND(JL, IRRACCSS)
    PRS_TEND(JL, IRSACCRG)=ZACC(JL) * PRS_TEND(JL, IRSACCRG)
  ENDDO
ELSE
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
    ZVEC1(1:IGACC) = MAX( 1.00001, MIN( REAL(ICEP%NACCLBDAS)-0.00001,           &
                          ICEP%XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + ICEP%XACCINTP2S ) )
    IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
    ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - REAL( IVEC1(1:IGACC) )
    !
    ZVEC2(1:IGACC) = MAX( 1.00001, MIN( REAL(ICEP%NACCLBDAR)-0.00001,           &
                          ICEP%XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + ICEP%XACCINTP2R ) )
    IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
    ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - REAL( IVEC2(1:IGACC) )
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
    WHERE(GACC(1:KSIZE))
      ZZW6(1:KSIZE) =                                                        & !! coef of RRACCS
            ICEP%XFRACCSS*( PLBDAS(1:KSIZE)**ICED%XCXS )*( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) ) &
       *( ICEP%XLBRACCS1/((PLBDAS(1:KSIZE)**2)               ) +                  &
          ICEP%XLBRACCS2/( PLBDAS(1:KSIZE)    * PLBDAR(1:KSIZE)    ) +                  &
          ICEP%XLBRACCS3/(               (PLBDAR(1:KSIZE)**2)) )/PLBDAR(1:KSIZE)**4
      PRS_TEND(1:KSIZE, IRRACCSS) =ZZW(1:KSIZE)*ZZW6(1:KSIZE)
    END WHERE
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
    WHERE(GACC(1:KSIZE))
      PRS_TEND(1:KSIZE, IRRACCS) = ZZW(1:KSIZE)*ZZW6(1:KSIZE)
    END WHERE
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
    WHERE(GACC(1:KSIZE))
      PRS_TEND(1:KSIZE, IRSACCRG) = ICEP%XFSACCRG*ZZW(1:KSIZE)*                    & ! RSACCRG
          ( PLBDAS(1:KSIZE)**(ICED%XCXS-ICED%XBS) )*( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) ) &
         *( ICEP%XLBSACCR1/((PLBDAR(1:KSIZE)**2)               ) +           &
            ICEP%XLBSACCR2/( PLBDAR(1:KSIZE)    * PLBDAS(1:KSIZE)    ) +           &
            ICEP%XLBSACCR3/(               (PLBDAS(1:KSIZE)**2)) )/PLBDAR(1:KSIZE)
    END WHERE
  ENDIF
ENDIF
!
DO JL=1, KSIZE
  ! More restrictive ACC mask to be used for accretion by negative temperature only
  ZACC(JL) = ZACC(JL) * &
           &MAX(0., -SIGN(1., PT(JL)-CST%XTT)) ! WHERE(PT(:)<XTT)
  PRRACCSS(JL)=ZACC(JL)*MIN(ZFREEZ_RATE(JL), PRS_TEND(JL, IRRACCSS))
  ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRRACCSS(JL))
  ZZW(JL) = MIN(1., ZFREEZ_RATE(JL) / MAX(1.E-20, PRS_TEND(JL, IRRACCS)-PRRACCSS(JL))) ! proportion we are able to freeze
  PRRACCSG(JL)=ZACC(JL)*ZZW(JL) * MAX(0., PRS_TEND(JL, IRRACCS)-PRRACCSS(JL))
  ZFREEZ_RATE(JL) = MAX(0., ZFREEZ_RATE(JL)-PRRACCSG(JL))
  PRSACCRG(JL)=ZACC(JL)*ZZW(JL) * PRS_TEND(JL, IRSACCRG)

  PRSACCRG(JL) = PRSACCRG(JL) * MAX(0., -SIGN(1., -PRRACCSG(JL)))
  PRRACCSG(JL)=MAX(0., PRRACCSG(JL))

ENDDO
!
!
!*       5.3    Conversion-Melting of the aggregates
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., ICED%XRTMIN(5)-PRST(JL))) * & ! WHERE(PRST(:)>XRTMIN(5))
           &MAX(0., -SIGN(1., CST%XTT-PT(JL))) * & ! WHERE(PT(:)>XTT)
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRSMLTG(JL)=ZMASK(JL)*PRSMLTG(JL)
    PRCMLTSR(JL)=ZMASK(JL)*PRCMLTSR(JL)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRSMLTG(JL)=ZMASK(JL)*PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(PARAMI%LEVLIMIT) THEN
    WHERE(ZMASK(:)==1.)
      PRSMLTG(:)=MIN(PRSMLTG(:), EXP(CST%XALPW-CST%XBETAW/PT(:)-CST%XGAMW*ALOG(PT(:)))) ! min(ev, es_w(T))
    END WHERE
  ENDIF
  DO JL=1, KSIZE
    PRSMLTG(JL)=ZMASK(JL)*( &
                            & PKA(JL)*(CST%XTT-PT(JL)) +                                 &
                            & ( PDV(JL)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(JL) - CST%XTT )) &
                            & *(CST%XESTT-PRSMLTG(JL))/(CST%XRV*PT(JL))             ) &
                          &)
  ENDDO
  PRCMLTSR(:) = 0.
  WHERE(ZMASK(1:KSIZE)==1.)
    !
    ! compute RSMLT
    !
    PRSMLTG(1:KSIZE)  = ICEP%XFSCVMG*MAX( 0.0,( -PRSMLTG(1:KSIZE) *             &
                         ( ICEP%X0DEPS*       PLBDAS(1:KSIZE)**ICEP%XEX0DEPS +     &
                           ICEP%X1DEPS*PCJ(1:KSIZE)*PLBDAS(1:KSIZE)**ICEP%XEX1DEPS ) -   &
                                   ( PRS_TEND(1:KSIZE, IRCRIMS) + PRS_TEND(1:KSIZE, IRRACCS) ) *       &
                            ( PRHODREF(1:KSIZE)*CST%XCL*(CST%XTT-PT(1:KSIZE))) ) /    &
                                           ( PRHODREF(1:KSIZE)*CST%XLMTT ) )
    ! When T < XTT, rc is collected by snow (riming) to produce snow and graupel
    ! When T > XTT, if riming was still enabled, rc would produce snow and graupel with snow becomming graupel (conversion/melting) and graupel becomming rain (melting)
    ! To insure consistency when crossing T=XTT, rc collected with T>XTT must be transformed in rain.
    ! rc cannot produce iced species with a positive temperature but is still collected with a good efficiency by snow
    PRCMLTSR(1:KSIZE) = PRS_TEND(1:KSIZE, IRCRIMS) ! both species are liquid, no heat is exchanged
  END WHERE
ENDIF

DO JL=1, KSIZE
  PA_RC(JL) = PA_RC(JL) - PRCRIMSS(JL)
  PA_RS(JL) = PA_RS(JL) + PRCRIMSS(JL)
  PA_TH(JL) = PA_TH(JL) + PRCRIMSS(JL)*(PLSFACT(JL)-PLVFACT(JL))
  PA_RC(JL) = PA_RC(JL) - PRCRIMSG(JL)
  PA_RS(JL) = PA_RS(JL) - PRSRIMCG(JL)
  PA_RG(JL) = PA_RG(JL) + PRCRIMSG(JL)+PRSRIMCG(JL)
  PA_TH(JL) = PA_TH(JL) + PRCRIMSG(JL)*(PLSFACT(JL)-PLVFACT(JL))

  PA_RR(JL) = PA_RR(JL) - PRRACCSS(JL)
  PA_RS(JL) = PA_RS(JL) + PRRACCSS(JL)
  PA_TH(JL) = PA_TH(JL) + PRRACCSS(JL)*(PLSFACT(JL)-PLVFACT(JL))
  PA_RR(JL) = PA_RR(JL) - PRRACCSG(JL)
  PA_RS(JL) = PA_RS(JL) - PRSACCRG(JL)
  PA_RG(JL) = PA_RG(JL) + PRRACCSG(JL)+PRSACCRG(JL)
  PA_TH(JL) = PA_TH(JL) + PRRACCSG(JL)*(PLSFACT(JL)-PLVFACT(JL))

  ! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
  ! because the graupeln produced by this process are still icy!!!
  PA_RS(JL) = PA_RS(JL) - PRSMLTG(JL)
  PA_RG(JL) = PA_RG(JL) + PRSMLTG(JL)
  PA_RC(JL) = PA_RC(JL) - PRCMLTSR(JL)
  PA_RR(JL) = PA_RR(JL) + PRCMLTSR(JL)
ENDDO
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RS', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_FAST_RS
END MODULE MODE_ICE4_FAST_RS
