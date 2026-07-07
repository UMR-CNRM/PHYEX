!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ICE4_FAST_RH
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RH(CST, PARAMI, ICEP, ICED, D, LDSOFT, LDCOMPUTE, LDWETG, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAS, PLBDAG, PLBDAR, PLBDAH, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                       &PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, &
                       &PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, PRDRYHG, PRHMLTR, &
                       &PRH_TEND)

!$ACDC singlecolumn

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
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
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
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(D%NIJT),   INTENT(IN)    :: LDCOMPUTE
LOGICAL, DIMENSION(D%NIJT),   INTENT(IN)    :: LDWETG    ! .TRUE. where graupel grows in wet mode
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAR   ! Slope parameter of the rain      distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAH   ! Slope parameter of the hail      distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRRT     ! Rain m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRHT     ! Hail m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRCWETH  ! Dry growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRIWETH  ! Dry growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRSWETH  ! Dry growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRGWETH  ! Dry growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRRWETH  ! Dry growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRCDRYH  ! Wet growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRIDRYH  ! Wet growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRSDRYH  ! Wet growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRRDRYH  ! Wet growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRGDRYH  ! Wet growth of hailstone
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRDRYHG  ! Conversion of hailstone into graupel
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRHMLTR  ! Melting of the hailstones
REAL, DIMENSION(D%NIJT, 10),  INTENT(INOUT) :: PRH_TEND ! Individual tendencies
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCWETH=1, IRRWETH=2, IRIDRYH=3, IRIWETH=4, IRSDRYH=5, IRSWETH=6, IRGDRYH=7, IRGWETH=8, &
                    & IFREEZ1=9, IFREEZ2=10
LOGICAL, DIMENSION(D%NIJT) :: GWET
INTEGER :: IGWET
REAL, DIMENSION(D%NIJT) :: ZBUF1, ZBUF2, ZBUF3
INTEGER, DIMENSION(D%NIJT) :: IBUF1, IBUF2, IBUF3
REAL, DIMENSION(D%NIJT) :: ZZW, &
                           ZRDRYH_INIT, ZRWETH_INIT, &
                           ZRDRYHG
INTEGER :: JIJ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
LOGICAL, DIMENSION(D%NIJT) :: LLWETH, LLDRYH
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RH',0,ZHOOK_HANDLE)
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!
!*       7.2    compute the Wet and Dry growth of hail
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRHT(JIJ)>ICED%XRTMIN(7) .AND. PRCT(JIJ)>ICED%XRTMIN(2) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRH_TEND(JIJ, IRCWETH)=PLBDAH(JIJ)**(ICED%XCXH-ICED%XDH-2.0) * PRHODREF(JIJ)**(-ICED%XCEXVT)
      PRH_TEND(JIJ, IRCWETH)=ICEP%XFWETH * PRCT(JIJ) * PRH_TEND(JIJ, IRCWETH)
    ENDIF
  ELSE
    PRH_TEND(JIJ, IRCWETH)=0.
  ENDIF

  IF(PRHT(JIJ)>ICED%XRTMIN(7) .AND. PRIT(JIJ)>ICED%XRTMIN(4) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRH_TEND(JIJ, IRIWETH)=PLBDAH(JIJ)**(ICED%XCXH-ICED%XDH-2.0) * PRHODREF(JIJ)**(-ICED%XCEXVT)
      PRH_TEND(JIJ, IRIWETH)=ICEP%XFWETH * PRIT(JIJ) * PRH_TEND(JIJ, IRIWETH)   ! RIWETH
      PRH_TEND(JIJ, IRIDRYH)=PRH_TEND(JIJ, IRIWETH)*(ICEP%XCOLIH*EXP(ICEP%XCOLEXIH*(PT(JIJ)-CST%XTT)))   ! RIDRYH
    ENDIF
  ELSE
    PRH_TEND(JIJ, IRIWETH)=0.
    PRH_TEND(JIJ, IRIDRYH)=0.
  ENDIF
ENDDO
!$mnh_end_do()


!
!*       7.2.1  accretion of aggregates on the hailstones
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF (PRHT(JIJ)>ICED%XRTMIN(7) .AND. PRST(JIJ)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JIJ)) THEN
    GWET(JIJ) = .TRUE.
  ELSE
    GWET(JIJ) = .FALSE.
    PRH_TEND(JIJ, IRSWETH)=0.
    PRH_TEND(JIJ, IRSDRYH)=0.
  ENDIF
ENDDO
!$mnh_end_do()

IF(.NOT. LDSOFT) THEN
   CALL INTERP_MICRO_2D(D, PLBDAH, PLBDAS, ICEP%NWETLBDAH, ICEP%NWETLBDAS, &
                       &ICEP%XWETINTP1H, ICEP%XWETINTP2H, ICEP%XWETINTP1S, ICEP%XWETINTP2S, &
                       &PARAMI%LPACK_INTERP, GWET, IBUF1, IBUF2, IBUF3, ZBUF1, ZBUF2, ZBUF3, &
                       &IGWET, &
                       &ICEP%XKER_SWETH, ZZW)
  IF(IGWET>0)THEN

    DO JIJ=D%NIJB, D%NIJE
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        IF (GWET(JIJ)) THEN
          PRH_TEND(JIJ, IRSWETH)=ICEP%XFSWETH*ZZW(JIJ)                       & ! RSWETH
                      *( PLBDAS(JIJ)**(ICED%XCXS-ICED%XBS) )*( PLBDAH(JIJ)**ICED%XCXH )  &
                           *( PRHODREF(JIJ)**(-ICED%XCEXVT-1.) )               &
                           *( ICEP%XLBSWETH1/( PLBDAH(JIJ)**2              ) + &
                              ICEP%XLBSWETH2/( PLBDAH(JIJ)   * PLBDAS(JIJ)   ) + &
                              ICEP%XLBSWETH3/(               PLBDAS(JIJ)**2) )
        END IF
      ELSE
        IF (GWET(JIJ)) THEN
          PRH_TEND(JIJ, IRSWETH)=ICEP%XFSWETH*ZZW(JIJ)                       & ! RSWETH
                      *( PRST(JIJ))*( PLBDAH(JIJ)**ICED%XCXH )  &
                           *( PRHODREF(JIJ)**(-ICED%XCEXVT) )               &
                           *( ICEP%XLBSWETH1/( PLBDAH(JIJ)**2              ) + &
                              ICEP%XLBSWETH2/( PLBDAH(JIJ)   * PLBDAS(JIJ)   ) + &
                              ICEP%XLBSWETH3/(               PLBDAS(JIJ)**2) )
        END IF
      ENDIF
      IF (GWET(JIJ)) THEN
        PRH_TEND(JIJ, IRSDRYH)=PRH_TEND(JIJ, IRSWETH)*(ICEP%XCOLSH*EXP(ICEP%XCOLEXSH*(PT(JIJ)-CST%XTT)))
      END IF
    END DO

  ENDIF
ENDIF
!
!*       7.2.6  accretion of graupeln on the hailstones
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF (PRHT(JIJ)>ICED%XRTMIN(7) .AND. PRGT(JIJ)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JIJ)) THEN
    GWET(JIJ) = .TRUE.
  ELSE
    GWET(JIJ) = .FALSE.
    PRH_TEND(JIJ, IRGWETH)=0.
    PRH_TEND(JIJ, IRGDRYH)=0.
  END IF
ENDDO
!$mnh_end_do()

IF(.NOT. LDSOFT) THEN
  CALL INTERP_MICRO_2D(D, PLBDAH, PLBDAG, ICEP%NWETLBDAH, ICEP%NWETLBDAG, &
                       &ICEP%XWETINTP1H, ICEP%XWETINTP2H, ICEP%XWETINTP1G, ICEP%XWETINTP2G, &
                       &PARAMI%LPACK_INTERP, GWET, IBUF1, IBUF2, IBUF3, ZBUF1, ZBUF2, ZBUF3, &
                       &IGWET, &
                       &ICEP%XKER_GWETH, ZZW)
  IF(IGWET>0)THEN

    DO JIJ=D%NIJB, D%NIJE
      IF (GWET(JIJ)) THEN
        PRH_TEND(JIJ, IRGWETH)=ICEP%XFGWETH*ZZW(JIJ)                       & ! RGWETH
                    *( PLBDAG(JIJ)**(ICED%XCXG-ICED%XBG) )*( PLBDAH(JIJ)**ICED%XCXH )  &
                         *( PRHODREF(JIJ)**(-ICED%XCEXVT-1.) )               &
                         *( ICEP%XLBGWETH1/( PLBDAH(JIJ)**2              ) + &
                            ICEP%XLBGWETH2/( PLBDAH(JIJ)   * PLBDAG(JIJ)   ) + &
                            ICEP%XLBGWETH3/(               PLBDAG(JIJ)**2) )
        PRH_TEND(JIJ, IRGDRYH)=PRH_TEND(JIJ, IRGWETH)
      END IF
      !When graupel grows in wet mode, graupel is wet (!) and collection efficiency must remain the same
      IF (GWET(JIJ) .AND. .NOT. LDWETG(JIJ)) THEN
        PRH_TEND(JIJ, IRGDRYH)=PRH_TEND(JIJ, IRGDRYH)*(ICEP%XCOLGH*EXP(ICEP%XCOLEXGH*(PT(JIJ)-CST%XTT)))
      END IF
    END DO

  END IF
ENDIF
!
!*       7.2.11  accretion of raindrops on the hailstones
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF (PRHT(JIJ)>ICED%XRTMIN(7) .AND. PRRT(JIJ)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JIJ)) THEN
    GWET(JIJ) = .TRUE.
  ELSE
    GWET(JIJ) = .FALSE.
    PRH_TEND(JIJ, IRRWETH)=0.
  ENDIF
ENDDO
!$mnh_end_do()

IF(.NOT. LDSOFT) THEN
  CALL INTERP_MICRO_2D(D, PLBDAH, PLBDAR, ICEP%NWETLBDAH, ICEP%NWETLBDAR, &
                       &ICEP%XWETINTP1H, ICEP%XWETINTP2H, ICEP%XWETINTP1R, ICEP%XWETINTP2R, &
                       &PARAMI%LPACK_INTERP, GWET, IBUF1, IBUF2, IBUF3, ZBUF1, ZBUF2, ZBUF3, &
                       &IGWET, &
                       &ICEP%XKER_RWETH, ZZW)
  IF(IGWET>0)THEN

    DO JIJ=D%NIJB, D%NIJE
      IF (GWET(JIJ)) THEN
        PRH_TEND(JIJ, IRRWETH) = ICEP%XFRWETH*ZZW(JIJ)                    & ! RRWETH
                        *( PLBDAR(JIJ)**(-4) )*( PLBDAH(JIJ)**ICED%XCXH ) &
                                 *( PRHODREF(JIJ)**(-ICED%XCEXVT-1.) )   &
                      *( ICEP%XLBRWETH1/( PLBDAH(JIJ)**2              ) + &
                         ICEP%XLBRWETH2/( PLBDAH(JIJ)   * PLBDAR(JIJ)   ) + &
                         ICEP%XLBRWETH3/(               PLBDAR(JIJ)**2) )
      END IF
    END DO

  ENDIF
ENDIF
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  ZRDRYH_INIT(JIJ)=PRH_TEND(JIJ, IRCWETH)+PRH_TEND(JIJ, IRIDRYH)+ &
                 &PRH_TEND(JIJ, IRSDRYH)+PRH_TEND(JIJ, IRRWETH)+PRH_TEND(JIJ, IRGDRYH)
ENDDO
!$mnh_end_do()

!
!*       7.3    compute the Wet growth of hail
!    and
!*       7.4    Select Wet or Dry case
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRHT(JIJ)>ICED%XRTMIN(7) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRH_TEND(JIJ, IFREEZ1)=PRVT(JIJ)*PPRES(JIJ)/(CST%XEPSILO+PRVT(JIJ)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRH_TEND(JIJ, IFREEZ1)=MIN(PRH_TEND(JIJ, IFREEZ1), EXP(CST%XALPI-CST%XBETAI/PT(JIJ)-CST%XGAMI*LOG(PT(JIJ)))) ! min(ev, es_i(T))
      ENDIF
      PRH_TEND(JIJ, IFREEZ1)=PKA(JIJ)*(CST%XTT-PT(JIJ)) +                              &
                            (PDV(JIJ)*(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JIJ)-CST%XTT)) &
                             *(CST%XESTT-PRH_TEND(JIJ, IFREEZ1))/(CST%XRV*PT(JIJ)))
      PRH_TEND(JIJ, IFREEZ1)=PRH_TEND(JIJ, IFREEZ1)* (ICEP%X0DEPH*        PLBDAH(JIJ)**ICEP%XEX0DEPH +     &
                                                    ICEP%X1DEPH*PCJ(JIJ)*PLBDAH(JIJ)**ICEP%XEX1DEPH)/ &
                            (PRHODREF(JIJ)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JIJ))))
      PRH_TEND(JIJ, IFREEZ2)=(PRHODREF(JIJ)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PT(JIJ)))) / &
                            (PRHODREF(JIJ)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JIJ))))
    ENDIF

    !We must agregate, at least, the cold species
    ZRWETH_INIT(JIJ)=MAX(PRH_TEND(JIJ, IRIWETH)+PRH_TEND(JIJ, IRSWETH)+PRH_TEND(JIJ, IRGWETH), &
                       &MAX(0., PRH_TEND(JIJ, IFREEZ1) + &
                               &PRH_TEND(JIJ, IFREEZ2)*(PRH_TEND(JIJ, IRIWETH)+PRH_TEND(JIJ, IRSWETH)+PRH_TEND(JIJ, IRGWETH))))

    !Wet case
    LLWETH(JIJ)=MAX(0., ZRWETH_INIT(JIJ)-PRH_TEND(JIJ, IRIWETH)-PRH_TEND(JIJ, IRSWETH)-PRH_TEND(JIJ, IRGWETH))<= &
               MAX(0., ZRDRYH_INIT(JIJ)-PRH_TEND(JIJ, IRIDRYH)-PRH_TEND(JIJ, IRSDRYH)-PRH_TEND(JIJ, IRGDRYH))
    IF(PARAMI%LNULLWETH) THEN
      LLWETH(JIJ) = LLWETH(JIJ) .AND. ZRDRYH_INIT(JIJ)>0.
    ELSE
      LLWETH(JIJ) = LLWETH(JIJ) .AND. ZRWETH_INIT(JIJ)>0.
    ENDIF
    IF(.NOT. PARAMI%LWETHPOST) THEN
      LLWETH(JIJ) = LLWETH(JIJ) .AND. PT(JIJ)<CST%XTT
    ENDIF

    !Dry case
    LLDRYH(JIJ)=PT(JIJ)<CST%XTT .AND. ZRDRYH_INIT(JIJ)>1.E-20 .AND. &
              &MAX(0., ZRWETH_INIT(JIJ)-PRH_TEND(JIJ, IRIWETH)-PRH_TEND(JIJ, IRSWETH))>&
              &MAX(0., ZRDRYH_INIT(JIJ)-PRH_TEND(JIJ, IRIDRYH)-PRH_TEND(JIJ, IRSDRYH))

  ELSE
    PRH_TEND(JIJ, IFREEZ1)=0.
    PRH_TEND(JIJ, IFREEZ2)=0.
    ZRWETH_INIT(JIJ)=0.
    LLWETH(JIJ)=.FALSE.
    LLDRYH(JIJ)=.FALSE.
  ENDIF
ENDDO
!$mnh_end_do()


IF(PARAMI%LCONVHG)THEN

  DO JIJ=D%NIJB, D%NIJE
    IF (LLDRYH(JIJ)) THEN
      ZRDRYHG(JIJ)=ZRDRYH_INIT(JIJ)*ZRWETH_INIT(JIJ)/(ZRDRYH_INIT(JIJ)+ZRWETH_INIT(JIJ))
    ELSE
      ZRDRYHG(JIJ)=0.
    END IF
  END DO

ELSE

  ZRDRYHG(:)=0.

ENDIF


!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(LLWETH(JIJ)) THEN
    PRCWETH(JIJ) = PRH_TEND(JIJ, IRCWETH)
    PRIWETH(JIJ) = PRH_TEND(JIJ, IRIWETH)
    PRSWETH(JIJ) = PRH_TEND(JIJ, IRSWETH)
    PRGWETH(JIJ) = PRH_TEND(JIJ, IRGWETH)
    !Collected minus aggregated
    PRRWETH(JIJ) = (ZRWETH_INIT(JIJ) - PRH_TEND(JIJ, IRIWETH) - &
                   PRH_TEND(JIJ, IRSWETH) - PRH_TEND(JIJ, IRGWETH) - &
                   PRH_TEND(JIJ, IRCWETH))
  ELSE
    PRCWETH(JIJ) = 0.
    PRIWETH(JIJ) = 0.
    PRSWETH(JIJ) = 0.
    PRGWETH(JIJ) = 0.
    PRRWETH(JIJ) = 0.
  ENDIF

  IF(LLDRYH(JIJ)) THEN
    PRCDRYH(JIJ) = PRH_TEND(JIJ, IRCWETH)
    PRIDRYH(JIJ) = PRH_TEND(JIJ, IRIDRYH)
    PRSDRYH(JIJ) = PRH_TEND(JIJ, IRSDRYH)
    PRRDRYH(JIJ) = PRH_TEND(JIJ, IRRWETH)
    PRGDRYH(JIJ) = PRH_TEND(JIJ, IRGDRYH)
    PRDRYHG(JIJ) = ZRDRYHG(JIJ)
  ELSE
    PRCDRYH(JIJ) = 0.
    PRIDRYH(JIJ) = 0.
    PRSDRYH(JIJ) = 0.
    PRRDRYH(JIJ) = 0.
    PRGDRYH(JIJ) = 0.
    PRDRYHG(JIJ) = 0.
  ENDIF
ENDDO
!$mnh_end_do()

!
!*       7.5    Melting of the hailstones
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRHT(JIJ)>ICED%XRTMIN(7) .AND. PT(JIJ)>CST%XTT .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRHMLTR(JIJ) = PRVT(JIJ)*PPRES(JIJ)/(CST%XEPSILO+PRVT(JIJ)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRHMLTR(JIJ)=MIN(PRHMLTR(JIJ), EXP(CST%XALPW-CST%XBETAW/PT(JIJ)-CST%XGAMW*LOG(PT(JIJ)))) ! min(ev, es_w(T))
      ENDIF
      PRHMLTR(JIJ) = PKA(JIJ)*(CST%XTT-PT(JIJ)) +                              &
                    PDV(JIJ)*(CST%XLVTT + (CST%XCPV - CST%XCL) * (PT(JIJ) - CST%XTT)) &
                           *(CST%XESTT-PRHMLTR(JIJ))/(CST%XRV*PT(JIJ))
      PRHMLTR(JIJ)  = MAX(0., (-PRHMLTR(JIJ) * (ICEP%X0DEPH*       PLBDAH(JIJ)**ICEP%XEX0DEPH +     &
                                              ICEP%X1DEPH*PCJ(JIJ)*PLBDAH(JIJ)**ICEP%XEX1DEPH) -   &
                              (PRH_TEND(JIJ, IRCWETH)+PRH_TEND(JIJ, IRRWETH)) *        &
                              (PRHODREF(JIJ)*CST%XCL*(CST%XTT-PT(JIJ)))) /    &
                             (PRHODREF(JIJ)*CST%XLMTT))
    ENDIF
  ELSE
    PRHMLTR(JIJ)=0.
  ENDIF
ENDDO
!$mnh_end_do()

!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RH', 1, ZHOOK_HANDLE)
!
CONTAINS
!
INCLUDE "interp_micro.func.h"
!
END SUBROUTINE ICE4_FAST_RH
END MODULE MODE_ICE4_FAST_RH
