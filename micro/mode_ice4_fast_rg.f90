!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ICE4_FAST_RG
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RG(CST, PARAMI, ICEP, ICED, KPROMA, KSIZE, LDSOFT, LDCOMPUTE, KRR, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, PCIT, &
                       &PLBDAR, PLBDAS, PLBDAG, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                       &PRGSI, PRGSI_MR, &
                       &LDWETG, &
                       &PRICFRRG, PRRCFRIG, PRICFRR, PRCWETG, PRIWETG, PRRWETG, PRSWETG, &
                       &PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, PRWETGH, PRWETGH_MR, PRGMLTR, &
                       &PRG_TEND)
!!
!!**  PURPOSE
!!    -------
!!      Computes the fast rg processes
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
LOGICAL, DIMENSION(KPROMA),   INTENT(IN)    :: LDCOMPUTE
INTEGER,                      INTENT(IN)    :: KRR      ! Number of moist variable
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRGSI    ! Graupel tendency by other processes
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRGSI_MR ! Graupel mr change by other processes
LOGICAL, DIMENSION(KPROMA),   INTENT(OUT)   :: LDWETG   ! .TRUE. where graupel grows in wet mode
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRICFRRG ! Rain contact freezing
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRRCFRIG ! Rain contact freezing
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRICFRR  ! Rain contact freezing
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRCWETG  ! Graupel wet growth
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRIWETG  ! Graupel wet growth
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRRWETG  ! Graupel wet growth
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRSWETG  ! Graupel wet growth
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRCDRYG  ! Graupel dry growth
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRIDRYG  ! Graupel dry growth
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRRDRYG  ! Graupel dry growth
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRSDRYG  ! Graupel dry growth
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRWETGH  ! Conversion of graupel into hail
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRWETGH_MR ! Conversion of graupel into hail, mr change
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRGMLTR  ! Melting of the graupel
REAL, DIMENSION(KPROMA, 8),   INTENT(INOUT) :: PRG_TEND ! Individual tendencies
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCDRYG=1, IRIDRYG=2, IRIWETG=3, IRSDRYG=4, IRSWETG=5, IRRDRYG=6, &
                    & IFREEZ1=7, IFREEZ2=8
LOGICAL, DIMENSION(KPROMA) :: GDRY, LLDRYG
INTEGER :: IGDRY
REAL, DIMENSION(KPROMA) :: ZBUF1, ZBUF2, ZBUF3
INTEGER, DIMENSION(KPROMA) :: IBUF1, IBUF2, IBUF3
REAL, DIMENSION(KPROMA) :: ZZW, &
                           ZRDRYG_INIT, & !Initial dry growth rate of the graupeln
                           ZRWETG_INIT !Initial wet growth rate of the graupeln
REAL :: ZZW0D
INTEGER :: JL

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RG', 0, ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       6.1    rain contact freezing
!
DO JL=1, KSIZE
  IF(PRIT(JL)>ICED%XRTMIN(4) .AND. PRRT(JL)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRICFRRG(JL) = ICEP%XICFRR*PRIT(JL)                & ! RICFRRG
                                *PLBDAR(JL)**ICEP%XEXICFRR    &
                                *PRHODREF(JL)**(-ICED%XCEXVT)
      PRRCFRIG(JL) = ICEP%XRCFRI*PCIT(JL)                & ! RRCFRIG
                                * PLBDAR(JL)**ICEP%XEXRCFRI    &
                                * PRHODREF(JL)**(-ICED%XCEXVT-1.)
      IF(PARAMI%LCRFLIMIT) THEN
        !Comparison between heat to be released (to freeze rain) and heat sink (rain and ice temperature change)
        !ZZW0D is the proportion of process that can take place
        ZZW0D=MAX(0., MIN(1., (PRICFRRG(JL)*CST%XCI+PRRCFRIG(JL)*CST%XCL)*(CST%XTT-PT(JL)) / &
                              MAX(1.E-20, CST%XLVTT*PRRCFRIG(JL))))
        PRRCFRIG(JL) = ZZW0D * PRRCFRIG(JL) !Part of rain that can be freezed
        PRICFRR(JL) = (1.-ZZW0D) * PRICFRRG(JL) !Part of collected pristine ice converted to rain
        PRICFRRG(JL) = ZZW0D * PRICFRRG(JL) !Part of collected pristine ice that lead to graupel
      ELSE
        PRICFRR(JL) = 0.
      ENDIF
    ENDIF
  ELSE
    PRICFRRG(JL)=0.
    PRRCFRIG(JL)=0.
    PRICFRR(JL)=0.
  ENDIF
ENDDO
!
!
!*       6.3    compute the graupel growth
!
! Wet and dry collection of rc and ri on graupel
DO JL=1, KSIZE
  IF(PRGT(JL)>ICED%XRTMIN(6) .AND. PRCT(JL)>ICED%XRTMIN(2) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRG_TEND(JL, IRCDRYG)=PLBDAG(JL)**(ICED%XCXG-ICED%XDG-2.) * PRHODREF(JL)**(-ICED%XCEXVT)
      PRG_TEND(JL, IRCDRYG)=ICEP%XFCDRYG * PRCT(JL) * PRG_TEND(JL, IRCDRYG)
    ENDIF
  ELSE
    PRG_TEND(JL, IRCDRYG)=0.
  ENDIF

  IF(PRGT(JL)>ICED%XRTMIN(6) .AND. PRIT(JL)>ICED%XRTMIN(4) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRG_TEND(JL, IRIDRYG)=PLBDAG(JL)**(ICED%XCXG-ICED%XDG-2.) * PRHODREF(JL)**(-ICED%XCEXVT)
      PRG_TEND(JL, IRIDRYG)=ICEP%XFIDRYG*EXP(ICEP%XCOLEXIG*(PT(JL)-CST%XTT))*PRIT(JL)*PRG_TEND(JL, IRIDRYG)
      PRG_TEND(JL, IRIWETG)=PRG_TEND(JL, IRIDRYG) / (ICEP%XCOLIG*EXP(ICEP%XCOLEXIG*(PT(JL)-CST%XTT)))
    ENDIF
  ELSE
    PRG_TEND(JL, IRIDRYG)=0.
    PRG_TEND(JL, IRIWETG)=0.
  ENDIF
ENDDO

! Wet and dry collection of rs on graupel (6.2.1)
DO JL = 1, KSIZE
  IF (PRST(JL)>ICED%XRTMIN(5) .AND. PRGT(JL)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JL)) THEN
    GDRY(JL) = .TRUE.
  ELSE
    GDRY(JL) = .FALSE.
    PRG_TEND(JL, IRSDRYG)=0.
    PRG_TEND(JL, IRSWETG)=0.
  ENDIF
ENDDO

IF(.NOT. LDSOFT) THEN
  CALL INTERP_MICRO_2D(KPROMA, KSIZE, PLBDAG(:), PLBDAS(:), ICEP%NDRYLBDAG, ICEP%NDRYLBDAS, &
                       &ICEP%XDRYINTP1G, ICEP%XDRYINTP2G, ICEP%XDRYINTP1S, ICEP%XDRYINTP2S, &
                       &PARAMI%LPACK_INTERP, GDRY(:), IBUF1(:), IBUF2(:), IBUF3(:), ZBUF1(:), ZBUF2(:), ZBUF3(:), &
                       &IGDRY, &
                       &ICEP%XKER_SDRYG(:,:), ZZW(:))
  IF(IGDRY>0)THEN
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GDRY(1:KSIZE))
      PRG_TEND(1:KSIZE, IRSWETG)=ICEP%XFSDRYG*ZZW(1:KSIZE)                         & ! RSDRYG
                                    / ICEP%XCOLSG &
#if defined(REPRO48) || defined(REPRO55)
                  *(PLBDAS(1:KSIZE)**(ICED%XCXS-ICED%XBS))*( PLBDAG(1:KSIZE)**ICED%XCXG )    &
                  *(PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.))                    &
#else
                  *(PRST(1:KSIZE))*( PLBDAG(1:KSIZE)**ICED%XCXG )    &
                  *(PRHODREF(1:KSIZE)**(-ICED%XCEXVT))                    &
#endif
                       *( ICEP%XLBSDRYG1/( PLBDAG(1:KSIZE)**2              ) + &
                          ICEP%XLBSDRYG2/( PLBDAG(1:KSIZE)   * PLBDAS(1:KSIZE)   ) + &
                          ICEP%XLBSDRYG3/(               PLBDAS(1:KSIZE)**2))
      PRG_TEND(1:KSIZE, IRSDRYG)=PRG_TEND(1:KSIZE, IRSWETG)*ICEP%XCOLSG*EXP(ICEP%XCOLEXSG*(PT(1:KSIZE)-CST%XTT))
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
  ENDIF
ENDIF
!
!*       6.2.6  accretion of raindrops on the graupeln
!
DO JL = 1, KSIZE
  IF (PRRT(JL)>ICED%XRTMIN(3) .AND. PRGT(JL)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JL)) THEN
    GDRY(JL) = .TRUE.
  ELSE
    GDRY(JL) = .FALSE.
    PRG_TEND(JL, IRRDRYG)=0.
  ENDIF
ENDDO
IF(.NOT. LDSOFT) THEN
  !
  CALL INTERP_MICRO_2D(KPROMA, KSIZE, PLBDAG(:), PLBDAR(:), ICEP%NDRYLBDAG, ICEP%NDRYLBDAR, &
                       &ICEP%XDRYINTP1G, ICEP%XDRYINTP2G, ICEP%XDRYINTP1R, ICEP%XDRYINTP2R, &
                       &PARAMI%LPACK_INTERP, GDRY(:), IBUF1(:), IBUF2(:), IBUF3(:), ZBUF1(:), ZBUF2(:), ZBUF3(:), &
                       &IGDRY, &
                       &ICEP%XKER_RDRYG(:,:), ZZW(:))
  IF(IGDRY>0) THEN
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GDRY(1:KSIZE))
      PRG_TEND(1:KSIZE, IRRDRYG) = ICEP%XFRDRYG*ZZW(1:KSIZE)                    & ! RRDRYG
                        *( PLBDAR(1:KSIZE)**(-4) )*( PLBDAG(1:KSIZE)**ICED%XCXG ) &
                               *( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) )   &
                    *( ICEP%XLBRDRYG1/( PLBDAG(1:KSIZE)**2              ) + &
                       ICEP%XLBRDRYG2/( PLBDAG(1:KSIZE)   * PLBDAR(1:KSIZE)   ) + &
                       ICEP%XLBRDRYG3/(               PLBDAR(1:KSIZE)**2) )
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
  ENDIF
ENDIF

DO JL=1, KSIZE
  ZRDRYG_INIT(JL)=PRG_TEND(JL, IRCDRYG)+PRG_TEND(JL, IRIDRYG)+ &
                 &PRG_TEND(JL, IRSDRYG)+PRG_TEND(JL, IRRDRYG)
ENDDO

!Freezing rate and growth mode
DO JL=1, KSIZE
  IF(PRGT(JL)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JL)) THEN
    !Freezing rate
    IF(.NOT. LDSOFT) THEN
      PRG_TEND(JL, IFREEZ1)=PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRG_TEND(JL, IFREEZ1)=MIN(PRG_TEND(JL, IFREEZ1), EXP(CST%XALPI-CST%XBETAI/PT(JL)-CST%XGAMI*ALOG(PT(JL)))) ! min(ev, es_i(T))
      ENDIF
      PRG_TEND(JL, IFREEZ1)=PKA(JL)*(CST%XTT-PT(JL)) +                              &
               (PDV(JL)*(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JL)-CST%XTT)) &
                             *(CST%XESTT-PRG_TEND(JL, IFREEZ1))/(CST%XRV*PT(JL))           )
      PRG_TEND(JL, IFREEZ1)=PRG_TEND(JL, IFREEZ1)* ( ICEP%X0DEPG*       PLBDAG(JL)**ICEP%XEX0DEPG +     &
                             ICEP%X1DEPG*PCJ(JL)*PLBDAG(JL)**ICEP%XEX1DEPG )/ &
                            ( PRHODREF(JL)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JL))) )
      PRG_TEND(JL, IFREEZ2)=(PRHODREF(JL)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PT(JL)))   ) / &
                            ( PRHODREF(JL)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JL))) )
    ENDIF
    ZRWETG_INIT(JL)=MAX(PRG_TEND(JL, IRIWETG)+PRG_TEND(JL, IRSWETG), &
                       &MAX(0., PRG_TEND(JL, IFREEZ1) + &
                       &        PRG_TEND(JL, IFREEZ2) * (PRG_TEND(JL, IRIWETG)+PRG_TEND(JL, IRSWETG))))

    !Growth mode
    LDWETG(JL)=MAX(0., ZRWETG_INIT(JL)-PRG_TEND(JL, IRIWETG)-PRG_TEND(JL, IRSWETG))<= &
              &MAX(0., ZRDRYG_INIT(JL)-PRG_TEND(JL, IRIDRYG)-PRG_TEND(JL, IRSDRYG))

    IF(PARAMI%LNULLWETG) THEN
      LDWETG(JL) = LDWETG(JL) .AND. ZRDRYG_INIT(JL)>0.
    ELSE
      LDWETG(JL) = LDWETG(JL) .AND. ZRWETG_INIT(JL)>0.
    ENDIF
    IF(.NOT. PARAMI%LWETGPOST) THEN
      LDWETG(JL) = LDWETG(JL) .AND. PT(JL)<CST%XTT
    ENDIF

#ifdef REPRO48
    LLDRYG(JL)=PT(JL)<CST%XTT .AND. ZRDRYG_INIT(JL)>0. .AND. &
#else
    LLDRYG(JL)=PT(JL)<CST%XTT .AND. ZRDRYG_INIT(JL)>1.E-20 .AND. &
#endif
              &MAX(0., ZRWETG_INIT(JL)-PRG_TEND(JL, IRIWETG)-PRG_TEND(JL, IRSWETG))>&
              &MAX(0., ZRDRYG_INIT(JL)-PRG_TEND(JL, IRIDRYG)-PRG_TEND(JL, IRSDRYG))
  ELSE
    PRG_TEND(JL, IFREEZ1)=0.
    PRG_TEND(JL, IFREEZ2)=0.
    ZRWETG_INIT(JL)=0.
    LDWETG(JL)=.FALSE.
    LLDRYG(JL)=.FALSE.
  ENDIF
ENDDO

! Part of ZRWETG to be converted into hail
! Graupel can be produced by other processes instantaneously (inducing a mixing ratio change, PRGSI_MR) or
! as a tendency (PRWETGH)
IF(KRR==7) THEN
  !$mnh_expand_where(JL=1:KSIZE)
  WHERE(LDWETG(1:KSIZE))
    !assume a linear percent of conversion of produced graupel into hail
    PRWETGH(1:KSIZE)=(MAX(0., PRGSI(1:KSIZE)+PRICFRRG(1:KSIZE)+PRRCFRIG(1:KSIZE))+ZRWETG_INIT(1:KSIZE))*&
                    &ZRDRYG_INIT(1:KSIZE)/(ZRWETG_INIT(1:KSIZE)+ZRDRYG_INIT(1:KSIZE))
    PRWETGH_MR(1:KSIZE)=MAX(0., PRGSI_MR(1:KSIZE))*ZRDRYG_INIT(1:KSIZE)/(ZRWETG_INIT(1:KSIZE)+ZRDRYG_INIT(1:KSIZE))
  ELSEWHERE
    PRWETGH(1:KSIZE)=0.
    PRWETGH_MR(1:KSIZE)=0.
  END WHERE
  !$mnh_end_expand_where(JL=1:KSIZE)
ELSE
  PRWETGH(:)=0.
  PRWETGH_MR(:)=0.
ENDIF

DO JL=1, KSIZE
  !Aggregated minus collected
  IF(LDWETG(JL)) THEN
    PRRWETG(JL)=-(PRG_TEND(JL, IRIWETG)+PRG_TEND(JL, IRSWETG)+&
                 &PRG_TEND(JL, IRCDRYG)-ZRWETG_INIT(JL))
    PRCWETG(JL)=PRG_TEND(JL, IRCDRYG)
    PRIWETG(JL)=PRG_TEND(JL, IRIWETG)
    PRSWETG(JL)=PRG_TEND(JL, IRSWETG)
  ELSE
    PRRWETG(JL)=0.
    PRCWETG(JL)=0.
    PRIWETG(JL)=0.
    PRSWETG(JL)=0.
  ENDIF

  IF(LLDRYG(JL)) THEN
    PRCDRYG(JL)=PRG_TEND(JL, IRCDRYG)
    PRRDRYG(JL)=PRG_TEND(JL, IRRDRYG)
    PRIDRYG(JL)=PRG_TEND(JL, IRIDRYG)
    PRSDRYG(JL)=PRG_TEND(JL, IRSDRYG)
  ELSE
    PRCDRYG(JL)=0.
    PRRDRYG(JL)=0.
    PRIDRYG(JL)=0.
    PRSDRYG(JL)=0.
  ENDIF
ENDDO
!
!*       6.5    Melting of the graupeln
!
DO JL=1, KSIZE
  IF(PRGT(JL)>ICED%XRTMIN(6) .AND. PT(JL)>CST%XTT .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRGMLTR(JL)=PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRGMLTR(JL)=MIN(PRGMLTR(JL), EXP(CST%XALPW-CST%XBETAW/PT(JL)-CST%XGAMW*ALOG(PT(JL)))) ! min(ev, es_w(T))
      ENDIF
      PRGMLTR(JL)=PKA(JL)*(CST%XTT-PT(JL)) +                                 &
                  PDV(JL)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(JL) - CST%XTT )) &
                         *(CST%XESTT-PRGMLTR(JL))/(CST%XRV*PT(JL)) 
      PRGMLTR(JL)=MAX(0., (-PRGMLTR(JL)*                     &
                           (ICEP%X0DEPG*       PLBDAG(JL)**ICEP%XEX0DEPG +     &
                            ICEP%X1DEPG*PCJ(JL)*PLBDAG(JL)**ICEP%XEX1DEPG) -   &
                           (PRG_TEND(JL, IRCDRYG)+PRG_TEND(JL, IRRDRYG)) *     &
                           (PRHODREF(JL)*CST%XCL*(CST%XTT-PT(JL)))) /    &
                          ( PRHODREF(JL)*CST%XLMTT ) )
    ENDIF
  ELSE
    PRGMLTR(JL)=0.
  ENDIF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RG', 1, ZHOOK_HANDLE)
!
CONTAINS
!
INCLUDE "interp_micro.func.h"
!
END SUBROUTINE ICE4_FAST_RG
!
END MODULE MODE_ICE4_FAST_RG
