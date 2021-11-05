!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_FAST_RS
INTERFACE
SUBROUTINE ICE4_FAST_RS(KSIZE, LDSOFT, PCOMPUTE, &
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
IMPLICIT NONE
INTEGER,                      INTENT(IN)    :: KSIZE
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
REAL, DIMENSION(KSIZE, 8),    INTENT(INOUT) :: PRS_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
END SUBROUTINE ICE4_FAST_RS
END INTERFACE
END MODULE MODI_ICE4_FAST_RS
SUBROUTINE ICE4_FAST_RS(KSIZE, LDSOFT, PCOMPUTE, &
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
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: XALPI,XALPW,XBETAI,XBETAW,XCI,XCL,XCPV,XESTT,XGAMI,XGAMW,XLMTT,XLVTT,XMD,XMV,XRV,XTT,  &
                               XEPSILO
USE MODD_PARAM_ICE,      ONLY: LEVLIMIT, CSNOWRIMING
USE MODD_RAIN_ICE_DESCR, ONLY: XBS,XCEXVT,XCXS,XRTMIN
USE MODD_RAIN_ICE_PARAM, ONLY: NACCLBDAR,NACCLBDAS,NGAMINC,X0DEPS,X1DEPS,XACCINTP1R,XACCINTP1S,XACCINTP2R,XACCINTP2S, &
                               XCRIMSG,XCRIMSS,XEX0DEPS,XEX1DEPS,XEXCRIMSG,XEXCRIMSS,XEXSRIMCG,XEXSRIMCG2,XFRACCSS,   &
                               XFSACCRG,XFSCVMG,XGAMINC_RIM1,XGAMINC_RIM1,XGAMINC_RIM2,XGAMINC_RIM4,XKER_RACCS,       &
                               XKER_RACCSS,XKER_SACCRG,XLBRACCS1,XLBRACCS2,XLBRACCS3,XLBSACCR1,XLBSACCR2,XLBSACCR3,   &
                               XRIMINTP1,XRIMINTP2,XSRIMCG,XSRIMCG2,XSRIMCG3
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
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
REAL, DIMENSION(KSIZE, 8),    INTENT(INOUT) :: PRS_TEND ! Individual tendencies
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
!
REAL, DIMENSION(KSIZE) :: ZRIM, ZACC, ZMASK
LOGICAL, DIMENSION(KSIZE) :: GRIM, GACC
INTEGER :: IGRIM, IGACC
INTEGER, DIMENSION(KSIZE) :: I1
REAL, DIMENSION(KSIZE) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(KSIZE) :: IVEC1, IVEC2
REAL, DIMENSION(KSIZE) :: ZZW, ZZW2, ZZW6, ZFREEZ_RATE
INTEGER :: JJ, JL
!-------------------------------------------------------------------------------
!
!
!*       5.0    maximum freezing rate
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(5)-PRST(JL))) * & ! WHERE(PRST(:)>XRTMIN(5))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRS_TEND(JL, IFREEZ1)=ZMASK(JL) * PRS_TEND(JL, IFREEZ1)
    PRS_TEND(JL, IFREEZ2)=ZMASK(JL) * PRS_TEND(JL, IFREEZ2)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRS_TEND(JL, IFREEZ1)=ZMASK(JL) * PRVT(JL)*PPRES(JL)/(XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(LEVLIMIT) THEN
    WHERE(ZMASK(:)==1.)
      PRS_TEND(:, IFREEZ1)=MIN(PRS_TEND(:, IFREEZ1), EXP(XALPI-XBETAI/PT(:)-XGAMI*ALOG(PT(:)))) ! min(ev, es_i(T))
    END WHERE
  ENDIF
  PRS_TEND(:, IFREEZ2)=0.
  WHERE(ZMASK(:)==1.)
    PRS_TEND(:, IFREEZ1)=PKA(:)*(XTT-PT(:)) +                              &
             (PDV(:)*(XLVTT+(XCPV-XCL)*(PT(:)-XTT)) &
                           *(XESTT-PRS_TEND(:, IFREEZ1))/(XRV*PT(:))           )
    PRS_TEND(:, IFREEZ1)=PRS_TEND(:, IFREEZ1)* ( X0DEPS*       PLBDAS(:)**XEX0DEPS +     &
                           X1DEPS*PCJ(:)*PLBDAS(:)**XEX1DEPS )/ &
                          ( PRHODREF(:)*(XLMTT-XCL*(XTT-PT(:))) )
    PRS_TEND(:, IFREEZ2)=(PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PT(:)))   ) / &
                          ( PRHODREF(:)*(XLMTT-XCL*(XTT-PT(:))) )
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
DO JJ = 1, SIZE(GRIM)
  ZRIM(JJ)=MAX(0., -SIGN(1., XRTMIN(2)-PRCT(JJ))) * & !WHERE(PRCT(:)>XRTMIN(2))
          &MAX(0., -SIGN(1., XRTMIN(5)-PRST(JJ))) * & !WHERE(PRST(:)>XRTMIN(5))
          &PCOMPUTE(JJ)
  IF (ZRIM(JJ)>0) THEN
    IGRIM = IGRIM + 1
    I1(IGRIM) = JJ
    GRIM(JJ) = .TRUE.
  ELSE
    GRIM(JJ) = .FALSE.
  END IF
END DO
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
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN( REAL(NGAMINC)-0.00001,           &
                          XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + XRIMINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )
    !
    !        5.1.3  perform the linear interpolation of the normalized
    !               "2+XDS"-moment of the incomplete gamma function
    !
    ZVEC1(1:IGRIM) =   XGAMINC_RIM1( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                     - XGAMINC_RIM1( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = 0.
    DO JJ = 1, IGRIM
      ZZW(I1(JJ)) = ZVEC1(JJ)
    END DO
    !
    !        5.1.4  riming of the small sized aggregates
    !
    WHERE (GRIM(:))
      PRS_TEND(:, IRCRIMSS) = XCRIMSS * ZZW(:) * PRCT(:)                & ! RCRIMSS
                                      *   PLBDAS(:)**XEXCRIMSS &
                                      * PRHODREF(:)**(-XCEXVT)
    END WHERE
    !
    !        5.1.5  perform the linear interpolation of the normalized
    !               "XBS"-moment of the incomplete gamma function (XGAMINC_RIM2) and
    !               "XBG"-moment of the incomplete gamma function (XGAMINC_RIM4)
    !
    ZVEC1(1:IGRIM) =  XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = 0.
    DO JJ = 1, IGRIM
      ZZW(I1(JJ)) = ZVEC1(JJ)
    END DO

    ZVEC1(1:IGRIM) =  XGAMINC_RIM4( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM4( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW2(:) = 0.
    DO JJ = 1, IGRIM
      ZZW2(I1(JJ)) = ZVEC1(JJ)
    END DO
    !
    !        5.1.6  riming-conversion of the large sized aggregates into graupeln
    !
    !
    WHERE(GRIM(:))
      PRS_TEND(:, IRCRIMS)=XCRIMSG * PRCT(:)               & ! RCRIMS
                                   * PLBDAS(:)**XEXCRIMSG  &
                                   * PRHODREF(:)**(-XCEXVT)
      ZZW6(:) = PRS_TEND(:, IRCRIMS) - PRS_TEND(:, IRCRIMSS) ! RCRIMSG
    END WHERE

    IF(CSNOWRIMING=='M90 ')THEN
      !Murakami 1990
      WHERE(GRIM(:))
        PRS_TEND(:, IRSRIMCG)=XSRIMCG * PLBDAS(:)**XEXSRIMCG*(1.0-ZZW(:))
        PRS_TEND(:, IRSRIMCG)=ZZW6(:)*PRS_TEND(:, IRSRIMCG)/ &
                       MAX(1.E-20, &
                           XSRIMCG3*XSRIMCG2*PLBDAS(:)**XEXSRIMCG2*(1.-ZZW2(:)) - &
                           XSRIMCG3*PRS_TEND(:, IRSRIMCG))
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
          &MAX(0., -SIGN(1., PT(JL)-XTT)) ! WHERE(PT(:)<XTT)
  PRCRIMSS(JL)=ZRIM(JL)*MIN(ZFREEZ_RATE(JL), PRS_TEND(JL, IRCRIMSS))
  ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRCRIMSS(JL))
  ZZW(JL) = MIN(1., ZFREEZ_RATE(JL) / MAX(1.E-20, PRS_TEND(JL, IRCRIMS) - PRCRIMSS(JL))) ! proportion we are able to freeze
  PRCRIMSG(JL) = ZRIM(JL) * ZZW(JL) * MAX(0., PRS_TEND(JL, IRCRIMS) - PRCRIMSS(JL)) ! RCRIMSG
  ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRCRIMSG(JL))
  PRSRIMCG(JL) = ZRIM(JL) * ZZW(JL) * PRS_TEND(JL, IRSRIMCG)

  PRSRIMCG(JL) = PRSRIMCG(JL) * MAX(0., -SIGN(1., -PRCRIMSG(JL)))
  PRCRIMSG(JL)=MAX(0., PRCRIMSG(JL))

  PA_RC(JL) = PA_RC(JL) - PRCRIMSS(JL)
  PA_RS(JL) = PA_RS(JL) + PRCRIMSS(JL)
  PA_TH(JL) = PA_TH(JL) + PRCRIMSS(JL)*(PLSFACT(JL)-PLVFACT(JL))
  PA_RC(JL) = PA_RC(JL) - PRCRIMSG(JL)
  PA_RS(JL) = PA_RS(JL) - PRSRIMCG(JL)
  PA_RG(JL) = PA_RG(JL) + PRCRIMSG(JL)+PRSRIMCG(JL)
  PA_TH(JL) = PA_TH(JL) + PRCRIMSG(JL)*(PLSFACT(JL)-PLVFACT(JL))
ENDDO
!
!*       5.2    rain accretion onto the aggregates
!
IGACC = 0
DO JJ = 1, SIZE(GACC)
  ZACC(JJ)=MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JJ))) * & !WHERE(PRRT(:)>XRTMIN(3))
          &MAX(0., -SIGN(1., XRTMIN(5)-PRST(JJ))) * & !WHERE(PRST(:)>XRTMIN(5))
          &PCOMPUTE(JJ)
  IF (ZACC(JJ)>0) THEN
    IGACC = IGACC + 1
    I1(IGACC) = JJ
    GACC(JJ) = .TRUE.
  ELSE
    GACC(JJ) = .FALSE.
  END IF
END DO

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
    END DO
    !
    !        5.2.2  find the next lower indice for the PLBDAS and for the PLBDAR
    !               in the geometrical set of (Lbda_s,Lbda_r) couplet use to
    !               tabulate the RACCSS-kernel
    !
    ZVEC1(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAS)-0.00001,           &
                          XACCINTP1S * LOG( ZVEC1(1:IGACC) ) + XACCINTP2S ) )
    IVEC1(1:IGACC) = INT( ZVEC1(1:IGACC) )
    ZVEC1(1:IGACC) = ZVEC1(1:IGACC) - REAL( IVEC1(1:IGACC) )
    !
    ZVEC2(1:IGACC) = MAX( 1.00001, MIN( REAL(NACCLBDAR)-0.00001,           &
                          XACCINTP1R * LOG( ZVEC2(1:IGACC) ) + XACCINTP2R ) )
    IVEC2(1:IGACC) = INT( ZVEC2(1:IGACC) )
    ZVEC2(1:IGACC) = ZVEC2(1:IGACC) - REAL( IVEC2(1:IGACC) )
    !
    !        5.2.3  perform the bilinear interpolation of the normalized
    !               RACCSS-kernel
    !
    DO JJ = 1, IGACC
      ZVEC3(JJ) =  (  XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * ZVEC1(JJ) &
                 - (  XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RACCSS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                          * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = 0.
    DO JJ = 1, IGACC
      ZZW(I1(JJ)) = ZVEC3(JJ)
    END DO
    !
    !        5.2.4  raindrop accretion on the small sized aggregates
    !
    WHERE(GACC(:))
      ZZW6(:) =                                                        & !! coef of RRACCS
            XFRACCSS*( PLBDAS(:)**XCXS )*( PRHODREF(:)**(-XCEXVT-1.) ) &
       *( XLBRACCS1/((PLBDAS(:)**2)               ) +                  &
          XLBRACCS2/( PLBDAS(:)    * PLBDAR(:)    ) +                  &
          XLBRACCS3/(               (PLBDAR(:)**2)) )/PLBDAR(:)**4
      PRS_TEND(:, IRRACCSS) =ZZW(:)*ZZW6(:)
    END WHERE
    !
    !        5.2.4b perform the bilinear interpolation of the normalized
    !               RACCS-kernel
    !
    DO JJ = 1, IGACC
      ZVEC3(JJ) =  (   XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  XKER_RACCS(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                   * ZVEC1(JJ) &
                 - (   XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    -  XKER_RACCS(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                           * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = 0.
    DO JJ = 1, IGACC
      ZZW(I1(JJ)) = ZVEC3(JJ)
    END DO
    WHERE(GACC(:))
      PRS_TEND(:, IRRACCS) = ZZW(:)*ZZW6(:)
    END WHERE
    !        5.2.5  perform the bilinear interpolation of the normalized
    !               SACCRG-kernel
    !
    DO JJ = 1, IGACC
        ZVEC3(JJ) =  (  XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                      - XKER_SACCRG(IVEC2(JJ)+1,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
                                                            * ZVEC2(JJ) &
                   - (  XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)+1)* ZVEC1(JJ)          &
                      - XKER_SACCRG(IVEC2(JJ)  ,IVEC1(JJ)  )*(ZVEC1(JJ) - 1.0) ) &
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
    WHERE(GACC(:))
      PRS_TEND(:, IRSACCRG) = XFSACCRG*ZZW(:)*                    & ! RSACCRG
          ( PLBDAS(:)**(XCXS-XBS) )*( PRHODREF(:)**(-XCEXVT-1.) ) &
         *( XLBSACCR1/((PLBDAR(:)**2)               ) +           &
            XLBSACCR2/( PLBDAR(:)    * PLBDAS(:)    ) +           &
            XLBSACCR3/(               (PLBDAS(:)**2)) )/PLBDAR(:)
    END WHERE
  ENDIF
ENDIF
!
DO JL=1, KSIZE
  ! More restrictive ACC mask to be used for accretion by negative temperature only
  ZACC(JL) = ZACC(JL) * &
           &MAX(0., -SIGN(1., PT(JL)-XTT)) ! WHERE(PT(:)<XTT)
  PRRACCSS(JL)=ZACC(JL)*MIN(ZFREEZ_RATE(JL), PRS_TEND(JL, IRRACCSS))
  ZFREEZ_RATE(JL)=MAX(0., ZFREEZ_RATE(JL)-PRRACCSS(JL))
  ZZW(JL) = MIN(1., ZFREEZ_RATE(JL) / MAX(1.E-20, PRS_TEND(JL, IRRACCS)-PRRACCSS(JL))) ! proportion we are able to freeze
  PRRACCSG(JL)=ZACC(JL)*ZZW(JL) * MAX(0., PRS_TEND(JL, IRRACCS)-PRRACCSS(JL))
  ZFREEZ_RATE(JL) = MAX(0., ZFREEZ_RATE(JL)-PRRACCSG(JL))
  PRSACCRG(JL)=ZACC(JL)*ZZW(JL) * PRS_TEND(JL, IRSACCRG)

  PRSACCRG(JL) = PRSACCRG(JL) * MAX(0., -SIGN(1., -PRRACCSG(JL)))
  PRRACCSG(JL)=MAX(0., PRRACCSG(JL))

  PA_RR(JL) = PA_RR(JL) - PRRACCSS(JL)
  PA_RS(JL) = PA_RS(JL) + PRRACCSS(JL)
  PA_TH(JL) = PA_TH(JL) + PRRACCSS(JL)*(PLSFACT(JL)-PLVFACT(JL))
  PA_RR(JL) = PA_RR(JL) - PRRACCSG(JL)
  PA_RS(JL) = PA_RS(JL) - PRSACCRG(JL)
  PA_RG(JL) = PA_RG(JL) + PRRACCSG(JL)+PRSACCRG(JL)
  PA_TH(JL) = PA_TH(JL) + PRRACCSG(JL)*(PLSFACT(JL)-PLVFACT(JL))
ENDDO
!
!
!*       5.3    Conversion-Melting of the aggregates
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(5)-PRST(JL))) * & ! WHERE(PRST(:)>XRTMIN(5))
           &MAX(0., -SIGN(1., XTT-PT(JL))) * & ! WHERE(PT(:)>XTT)
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRSMLTG(JL)=ZMASK(JL)*PRSMLTG(JL)
    PRCMLTSR(JL)=ZMASK(JL)*PRCMLTSR(JL)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRSMLTG(JL)=ZMASK(JL)*PRVT(JL)*PPRES(JL)/(XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(LEVLIMIT) THEN
    WHERE(ZMASK(:)==1.)
      PRSMLTG(:)=MIN(PRSMLTG(:), EXP(XALPW-XBETAW/PT(:)-XGAMW*ALOG(PT(:)))) ! min(ev, es_w(T))
    END WHERE
  ENDIF
  DO JL=1, KSIZE
    PRSMLTG(JL)=ZMASK(JL)*( &
                            & PKA(JL)*(XTT-PT(JL)) +                                 &
                            & ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PT(JL) - XTT )) &
                            & *(XESTT-PRSMLTG(JL))/(XRV*PT(JL))             ) &
                          &)
  ENDDO
  PRCMLTSR(:) = 0.
  WHERE(ZMASK(:)==1.)
    !
    ! compute RSMLT
    !
    PRSMLTG(:)  = XFSCVMG*MAX( 0.0,( -PRSMLTG(:) *             &
                         ( X0DEPS*       PLBDAS(:)**XEX0DEPS +     &
                           X1DEPS*PCJ(:)*PLBDAS(:)**XEX1DEPS ) -   &
                                   ( PRS_TEND(:, IRCRIMS) + PRS_TEND(:, IRRACCS) ) *       &
                            ( PRHODREF(:)*XCL*(XTT-PT(:))) ) /    &
                                           ( PRHODREF(:)*XLMTT ) )
    ! When T < XTT, rc is collected by snow (riming) to produce snow and graupel
    ! When T > XTT, if riming was still enabled, rc would produce snow and graupel with snow becomming graupel (conversion/melting) and graupel becomming rain (melting)
    ! To insure consistency when crossing T=XTT, rc collected with T>XTT must be transformed in rain.
    ! rc cannot produce iced species with a positive temperature but is still collected with a good efficiency by snow
    PRCMLTSR(:) = PRS_TEND(:, IRCRIMS) ! both species are liquid, no heat is exchanged
  END WHERE
ENDIF
DO JL=1, KSIZE
  ! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
  ! because the graupeln produced by this process are still icy!!!
  PA_RS(JL) = PA_RS(JL) - PRSMLTG(JL)
  PA_RG(JL) = PA_RG(JL) + PRSMLTG(JL)
  PA_RC(JL) = PA_RC(JL) - PRCMLTSR(JL)
  PA_RR(JL) = PA_RR(JL) + PRCMLTSR(JL)
ENDDO

!
END SUBROUTINE ICE4_FAST_RS
