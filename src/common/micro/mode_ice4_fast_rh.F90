!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ICE4_FAST_RH
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RH(CST, PARAMI, ICEP, ICED, KPROMA, KSIZE, LDSOFT, LDCOMPUTE, LDWETG, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAS, PLBDAG, PLBDAR, PLBDAH, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                       &PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, &
                       &PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, PRDRYHG, PRHMLTR, &
                       &PRH_TEND)
!!
!!**  PURPOSE
!!    -------
!!      Computes the fast rh process
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
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),            INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)    :: ICED
INTEGER,                      INTENT(IN)    :: KPROMA, KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KPROMA),   INTENT(IN)    :: LDCOMPUTE
LOGICAL, DIMENSION(KPROMA),   INTENT(IN)    :: LDWETG    ! .TRUE. where graupel grows in wet mode
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAR   ! Slope parameter of the rain      distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PLBDAH   ! Slope parameter of the hail      distribution
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRRT     ! Rain m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(IN)    :: PRHT     ! Hail m.r. at t
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRCWETH  ! Dry growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRIWETH  ! Dry growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRSWETH  ! Dry growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRGWETH  ! Dry growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRRWETH  ! Dry growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRCDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRIDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRSDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRRDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRGDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KPROMA),      INTENT(OUT)   :: PRDRYHG  ! Conversion of hailstone into graupel
REAL, DIMENSION(KPROMA),      INTENT(INOUT) :: PRHMLTR  ! Melting of the hailstones
REAL, DIMENSION(KPROMA, 10),  INTENT(INOUT) :: PRH_TEND ! Individual tendencies
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCWETH=1, IRRWETH=2, IRIDRYH=3, IRIWETH=4, IRSDRYH=5, IRSWETH=6, IRGDRYH=7, IRGWETH=8, &
                    & IFREEZ1=9, IFREEZ2=10
LOGICAL, DIMENSION(KPROMA) :: GWET
INTEGER :: IGWET
REAL, DIMENSION(KPROMA) :: ZBUF1, ZBUF2, ZBUF3
INTEGER, DIMENSION(KPROMA) :: IBUF1, IBUF2, IBUF3
REAL, DIMENSION(KPROMA) :: ZZW, &
                           ZRDRYH_INIT, ZRWETH_INIT, &
                           ZRDRYHG
INTEGER :: JL
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL, DIMENSION(KPROMA) :: LLWETH, LLDRYH
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RH',0,ZHOOK_HANDLE)
!
!
!*       7.2    compute the Wet and Dry growth of hail
!
DO JL=1, KSIZE
  IF(PRHT(JL)>ICED%XRTMIN(7) .AND. PRCT(JL)>ICED%XRTMIN(2) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRH_TEND(JL, IRCWETH)=PLBDAH(JL)**(ICED%XCXH-ICED%XDH-2.0) * PRHODREF(JL)**(-ICED%XCEXVT)
      PRH_TEND(JL, IRCWETH)=ICEP%XFWETH * PRCT(JL) * PRH_TEND(JL, IRCWETH)
    ENDIF
  ELSE
    PRH_TEND(JL, IRCWETH)=0.
  ENDIF

  IF(PRHT(JL)>ICED%XRTMIN(7) .AND. PRIT(JL)>ICED%XRTMIN(4) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRH_TEND(JL, IRIWETH)=PLBDAH(JL)**(ICED%XCXH-ICED%XDH-2.0) * PRHODREF(JL)**(-ICED%XCEXVT)
      PRH_TEND(JL, IRIWETH)=ICEP%XFWETH * PRIT(JL) * PRH_TEND(JL, IRIWETH)   ! RIWETH
      PRH_TEND(JL, IRIDRYH)=PRH_TEND(JL, IRIWETH)*(ICEP%XCOLIH*EXP(ICEP%XCOLEXIH*(PT(JL)-CST%XTT)))   ! RIDRYH
    ENDIF
  ELSE
    PRH_TEND(JL, IRIWETH)=0.
    PRH_TEND(JL, IRIDRYH)=0.
  ENDIF
ENDDO

!
!*       7.2.1  accretion of aggregates on the hailstones
!
DO JL = 1, KSIZE
  IF (PRHT(JL)>ICED%XRTMIN(7) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL)) THEN
    GWET(JL) = .TRUE.
  ELSE
    GWET(JL) = .FALSE.
    PRH_TEND(JL, IRSWETH)=0.
    PRH_TEND(JL, IRSDRYH)=0.
  ENDIF
ENDDO
IF(.NOT. LDSOFT) THEN
   CALL INTERP_MICRO_2D(KPROMA, KSIZE, PLBDAH(:), PLBDAS(:), ICEP%NWETLBDAH, ICEP%NWETLBDAS, &
                       &ICEP%XWETINTP1H, ICEP%XWETINTP2H, ICEP%XWETINTP1S, ICEP%XWETINTP2S, &
                       &PARAMI%LPACK_INTERP, GWET(:), IBUF1(:), IBUF2(:), IBUF3(:), ZBUF1(:), ZBUF2(:), ZBUF3(:), &
                       &IGWET, &
                       &ICEP%XKER_SWETH(:,:), ZZW(:))
  IF(IGWET>0)THEN
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GWET(1:KSIZE))
#ifdef REPRO48
      PRH_TEND(1:KSIZE, IRSWETH)=ICEP%XFSWETH*ZZW(1:KSIZE)                       & ! RSWETH
                    *( PLBDAS(1:KSIZE)**(ICED%XCXS-ICED%XBS) )*( PLBDAH(1:KSIZE)**ICED%XCXH )  &
                       *( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) )               &
                       *( ICEP%XLBSWETH1/( PLBDAH(1:KSIZE)**2              ) + &                                  
                          ICEP%XLBSWETH2/( PLBDAH(1:KSIZE)   * PLBDAS(1:KSIZE)   ) + &
                          ICEP%XLBSWETH3/(               PLBDAS(1:KSIZE)**2) )
#else
      PRH_TEND(1:KSIZE, IRSWETH)=ICEP%XFSWETH*ZZW(1:KSIZE)                       & ! RSWETH
                    *( PRST(1:KSIZE))*( PLBDAH(1:KSIZE)**ICED%XCXH )  &
                       *( PRHODREF(1:KSIZE)**(-ICED%XCEXVT) )               &
                       *( ICEP%XLBSWETH1/( PLBDAH(1:KSIZE)**2              ) + &                                  
                          ICEP%XLBSWETH2/( PLBDAH(1:KSIZE)   * PLBDAS(1:KSIZE)   ) + &
                          ICEP%XLBSWETH3/(               PLBDAS(1:KSIZE)**2) )
#endif
      PRH_TEND(1:KSIZE, IRSDRYH)=PRH_TEND(1:KSIZE, IRSWETH)*(ICEP%XCOLSH*EXP(ICEP%XCOLEXSH*(PT(1:KSIZE)-CST%XTT)))
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
  ENDIF
ENDIF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
DO JL = 1, KSIZE
  IF (PRHT(JL)>ICED%XRTMIN(7) .AND. PRGT(JL)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JL)) THEN
    GWET(JL) = .TRUE.
  ELSE
    GWET(JL) = .FALSE.
    PRH_TEND(JL, IRGWETH)=0.
    PRH_TEND(JL, IRGDRYH)=0.
  END IF
ENDDO
IF(.NOT. LDSOFT) THEN
  CALL INTERP_MICRO_2D(KPROMA, KSIZE, PLBDAH(:), PLBDAG(:), ICEP%NWETLBDAH, ICEP%NWETLBDAG, &
                       &ICEP%XWETINTP1H, ICEP%XWETINTP2H, ICEP%XWETINTP1G, ICEP%XWETINTP2G, &
                       &PARAMI%LPACK_INTERP, GWET(:), IBUF1(:), IBUF2(:), IBUF3(:), ZBUF1(:), ZBUF2(:), ZBUF3(:), &
                       &IGWET, &
                       &ICEP%XKER_GWETH(:,:), ZZW(:))
  IF(IGWET>0)THEN
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GWET(1:KSIZE))
      PRH_TEND(1:KSIZE, IRGWETH)=ICEP%XFGWETH*ZZW(1:KSIZE)                       & ! RGWETH
                    *( PLBDAG(1:KSIZE)**(ICED%XCXG-ICED%XBG) )*( PLBDAH(1:KSIZE)**ICED%XCXH )  &
                       *( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) )               &
                       *( ICEP%XLBGWETH1/( PLBDAH(1:KSIZE)**2              ) + &
                          ICEP%XLBGWETH2/( PLBDAH(1:KSIZE)   * PLBDAG(1:KSIZE)   ) + &
                          ICEP%XLBGWETH3/(               PLBDAG(1:KSIZE)**2) )
      PRH_TEND(1:KSIZE, IRGDRYH)=PRH_TEND(1:KSIZE, IRGWETH)
    END WHERE
    !When graupel grows in wet mode, graupel is wet (!) and collection efficiency must remain the same
    WHERE(GWET(1:KSIZE) .AND. .NOT. LDWETG(1:KSIZE))
      PRH_TEND(1:KSIZE, IRGDRYH)=PRH_TEND(1:KSIZE, IRGDRYH)*(ICEP%XCOLGH*EXP(ICEP%XCOLEXGH*(PT(1:KSIZE)-CST%XTT)))
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
  END IF
ENDIF
!
!*       7.2.11  accretion of raindrops on the hailstones
!
DO JL = 1, KSIZE
  IF (PRHT(JL)>ICED%XRTMIN(7) .AND. PRRT(JL)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JL)) THEN
    GWET(JL) = .TRUE.
  ELSE
    GWET(JL) = .FALSE.
    PRH_TEND(JL, IRRWETH)=0.
  ENDIF
ENDDO
IF(.NOT. LDSOFT) THEN
  CALL INTERP_MICRO_2D(KPROMA, KSIZE, PLBDAH(:), PLBDAR(:), ICEP%NWETLBDAH, ICEP%NWETLBDAR, &
                       &ICEP%XWETINTP1H, ICEP%XWETINTP2H, ICEP%XWETINTP1R, ICEP%XWETINTP2R, &
                       &PARAMI%LPACK_INTERP, GWET(:), IBUF1(:), IBUF2(:), IBUF3(:), ZBUF1(:), ZBUF2(:), ZBUF3(:), &
                       &IGWET, &
                       &ICEP%XKER_RWETH(:,:), ZZW(:))
  IF(IGWET>0)THEN
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GWET(1:KSIZE))
      PRH_TEND(1:KSIZE, IRRWETH) = ICEP%XFRWETH*ZZW(1:KSIZE)                    & ! RRWETH
                        *( PLBDAR(1:KSIZE)**(-4) )*( PLBDAH(1:KSIZE)**ICED%XCXH ) &
                               *( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) )   &
                    *( ICEP%XLBRWETH1/( PLBDAH(1:KSIZE)**2              ) + &
                       ICEP%XLBRWETH2/( PLBDAH(1:KSIZE)   * PLBDAR(1:KSIZE)   ) + &
                       ICEP%XLBRWETH3/(               PLBDAR(1:KSIZE)**2) )
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
  ENDIF
ENDIF
!
DO JL=1, KSIZE
  ZRDRYH_INIT(JL)=PRH_TEND(JL, IRCWETH)+PRH_TEND(JL, IRIDRYH)+ &
                 &PRH_TEND(JL, IRSDRYH)+PRH_TEND(JL, IRRWETH)+PRH_TEND(JL, IRGDRYH)
ENDDO
!
!*       7.3    compute the Wet growth of hail
!    and
!*       7.4    Select Wet or Dry case
!
DO JL=1, KSIZE
  IF(PRHT(JL)>ICED%XRTMIN(7) .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRH_TEND(JL, IFREEZ1)=PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRH_TEND(JL, IFREEZ1)=MIN(PRH_TEND(JL, IFREEZ1), EXP(CST%XALPI-CST%XBETAI/PT(JL)-CST%XGAMI*ALOG(PT(JL)))) ! min(ev, es_i(T))
      ENDIF
      PRH_TEND(JL, IFREEZ1)=PKA(JL)*(CST%XTT-PT(JL)) +                              &
                            (PDV(JL)*(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JL)-CST%XTT)) &
                             *(CST%XESTT-PRH_TEND(JL, IFREEZ1))/(CST%XRV*PT(JL)))
      PRH_TEND(JL, IFREEZ1)=PRH_TEND(JL, IFREEZ1)* (ICEP%X0DEPH*        PLBDAH(JL)**ICEP%XEX0DEPH +     &
                                                    ICEP%X1DEPH*PCJ(JL)*PLBDAH(JL)**ICEP%XEX1DEPH)/ &
                            (PRHODREF(JL)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JL))))
      PRH_TEND(JL, IFREEZ2)=(PRHODREF(JL)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PT(JL)))) / &
                            (PRHODREF(JL)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JL))))
    ENDIF

    !We must agregate, at least, the cold species
    ZRWETH_INIT(JL)=MAX(PRH_TEND(JL, IRIWETH)+PRH_TEND(JL, IRSWETH)+PRH_TEND(JL, IRGWETH), &
                       &MAX(0., PRH_TEND(JL, IFREEZ1) + &
                               &PRH_TEND(JL, IFREEZ2)*(PRH_TEND(JL, IRIWETH)+PRH_TEND(JL, IRSWETH)+PRH_TEND(JL, IRGWETH))))

    !Wet case
    LLWETH(JL)=MAX(0., ZRWETH_INIT(JL)-PRH_TEND(JL, IRIWETH)-PRH_TEND(JL, IRSWETH)-PRH_TEND(JL, IRGWETH))<= &
               MAX(0., ZRDRYH_INIT(JL)-PRH_TEND(JL, IRIDRYH)-PRH_TEND(JL, IRSDRYH)-PRH_TEND(JL, IRGDRYH))
    IF(PARAMI%LNULLWETH) THEN
      LLWETH(JL) = LLWETH(JL) .AND. ZRDRYH_INIT(JL)>0.
    ELSE
      LLWETH(JL) = LLWETH(JL) .AND. ZRWETH_INIT(JL)>0.
    ENDIF
    IF(.NOT. PARAMI%LWETHPOST) THEN
      LLWETH(JL) = LLWETH(JL) .AND. PT(JL)<CST%XTT
    ENDIF

    !Dry case
    LLDRYH(JL)=PT(JL)<CST%XTT .AND. ZRDRYH_INIT(JL)>1.E-20 .AND. &
              &MAX(0., ZRWETH_INIT(JL)-PRH_TEND(JL, IRIWETH)-PRH_TEND(JL, IRSWETH))>&
              &MAX(0., ZRDRYH_INIT(JL)-PRH_TEND(JL, IRIDRYH)-PRH_TEND(JL, IRSDRYH))

  ELSE
    PRH_TEND(JL, IFREEZ1)=0.
    PRH_TEND(JL, IFREEZ2)=0.
    ZRWETH_INIT(JL)=0.
    LLWETH(JL)=.FALSE.
    LLDRYH(JL)=.FALSE.
  ENDIF
ENDDO

IF(PARAMI%LCONVHG)THEN
  !$mnh_expand_where(JL=1:KSIZE)
  WHERE(LLDRYH(1:KSIZE))
    ZRDRYHG(1:KSIZE)=ZRDRYH_INIT(1:KSIZE)*ZRWETH_INIT(1:KSIZE)/(ZRDRYH_INIT(1:KSIZE)+ZRWETH_INIT(1:KSIZE))
  ELSEWHERE
    ZRDRYHG(1:KSIZE)=0.
  END WHERE
  !$mnh_end_expand_where(JL=1:KSIZE)
ELSE
  ZRDRYHG(:)=0.
ENDIF

DO JL=1, KSIZE
  IF(LLWETH(JL)) THEN
    PRCWETH(JL) = PRH_TEND(JL, IRCWETH)
    PRIWETH(JL) = PRH_TEND(JL, IRIWETH)
    PRSWETH(JL) = PRH_TEND(JL, IRSWETH)
    PRGWETH(JL) = PRH_TEND(JL, IRGWETH)
    !Collected minus aggregated
    PRRWETH(JL) = (ZRWETH_INIT(JL) - PRH_TEND(JL, IRIWETH) - &
                   PRH_TEND(JL, IRSWETH) - PRH_TEND(JL, IRGWETH) - &
                   PRH_TEND(JL, IRCWETH))
  ELSE
    PRCWETH(JL) = 0.
    PRIWETH(JL) = 0.
    PRSWETH(JL) = 0.
    PRGWETH(JL) = 0.
    PRRWETH(JL) = 0.
  ENDIF

  IF(LLDRYH(JL)) THEN
    PRCDRYH(JL) = PRH_TEND(JL, IRCWETH)
    PRIDRYH(JL) = PRH_TEND(JL, IRIDRYH)
    PRSDRYH(JL) = PRH_TEND(JL, IRSDRYH)
    PRRDRYH(JL) = PRH_TEND(JL, IRRWETH)
    PRGDRYH(JL) = PRH_TEND(JL, IRGDRYH)
    PRDRYHG(JL) = ZRDRYHG(JL)
  ELSE
    PRCDRYH(JL) = 0.
    PRIDRYH(JL) = 0.
    PRSDRYH(JL) = 0.
    PRRDRYH(JL) = 0.
    PRGDRYH(JL) = 0.
    PRDRYHG(JL) = 0.
  ENDIF
ENDDO
!
!*       7.5    Melting of the hailstones
!
DO JL=1, KSIZE
  IF(PRHT(JL)>ICED%XRTMIN(7) .AND. PT(JL)>CST%XTT .AND. LDCOMPUTE(JL)) THEN
    IF(.NOT. LDSOFT) THEN
      PRHMLTR(JL) = PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRHMLTR(JL)=MIN(PRHMLTR(JL), EXP(CST%XALPW-CST%XBETAW/PT(JL)-CST%XGAMW*ALOG(PT(JL)))) ! min(ev, es_w(T))
      ENDIF
      PRHMLTR(JL) = PKA(JL)*(CST%XTT-PT(JL)) +                              &
                    PDV(JL)*(CST%XLVTT + (CST%XCPV - CST%XCL) * (PT(JL) - CST%XTT)) &
                           *(CST%XESTT-PRHMLTR(JL))/(CST%XRV*PT(JL))
      PRHMLTR(JL)  = MAX(0., (-PRHMLTR(JL) * (ICEP%X0DEPH*       PLBDAH(JL)**ICEP%XEX0DEPH +     &
                                              ICEP%X1DEPH*PCJ(JL)*PLBDAH(JL)**ICEP%XEX1DEPH) -   &
                              (PRH_TEND(JL, IRCWETH)+PRH_TEND(JL, IRRWETH)) *        &
                              (PRHODREF(JL)*CST%XCL*(CST%XTT-PT(JL)))) /    &
                             (PRHODREF(JL)*CST%XLMTT))
    ENDIF
  ELSE
    PRHMLTR(JL)=0.
  ENDIF
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RH', 1, ZHOOK_HANDLE)
!
CONTAINS
!
INCLUDE "interp_micro.func.h"
!
END SUBROUTINE ICE4_FAST_RH
END MODULE MODE_ICE4_FAST_RH
