!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ICE4_FAST_RG
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RG(CST, PARAMI, ICEP, ICED, D, LDSOFT, LDCOMPUTE, KRR, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, PCIT, &
                       &PLBDAR, PLBDAS, PLBDAG, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                       &PRGSI, PRGSI_MR, &
                       &LDWETG, &
                       &PRICFRRG, PRRCFRIG, PRICFRR, PRCWETG, PRIWETG, PRRWETG, PRSWETG, &
                       &PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, PRWETGH, PRWETGH_MR, PRGMLTR, &
                       &PRG_TEND)

!$ACDC singlecolumn

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
!  K.I Ivarsson   02/2023: Some modifications that can be activated from namelist,
!!                           e.g. for better forecasts of supercooled rain.
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
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(D%NIJT),   INTENT(IN)    :: LDCOMPUTE
INTEGER,                      INTENT(IN)    :: KRR      ! Number of moist variable
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRGSI    ! Graupel tendency by other processes
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRGSI_MR ! Graupel mr change by other processes
LOGICAL, DIMENSION(D%NIJT),   INTENT(OUT)   :: LDWETG   ! .TRUE. where graupel grows in wet mode
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRICFRRG ! Rain contact freezing
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRRCFRIG ! Rain contact freezing
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRICFRR  ! Rain contact freezing
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRCWETG  ! Graupel wet growth
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRIWETG  ! Graupel wet growth
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRRWETG  ! Graupel wet growth
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRSWETG  ! Graupel wet growth
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRCDRYG  ! Graupel dry growth
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRIDRYG  ! Graupel dry growth
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRRDRYG  ! Graupel dry growth
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRSDRYG  ! Graupel dry growth
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRWETGH  ! Conversion of graupel into hail
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRWETGH_MR ! Conversion of graupel into hail, mr change
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRGMLTR  ! Melting of the graupel
REAL, DIMENSION(D%NIJT, 8),   INTENT(INOUT) :: PRG_TEND ! Individual tendencies
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCDRYG=1, IRIDRYG=2, IRIWETG=3, IRSDRYG=4, IRSWETG=5, IRRDRYG=6, &
                    & IFREEZ1=7, IFREEZ2=8
LOGICAL, DIMENSION(D%NIJT) :: GDRY, LLDRYG
INTEGER :: IGDRY
REAL, DIMENSION(D%NIJT) :: ZBUF1, ZBUF2, ZBUF3
INTEGER, DIMENSION(D%NIJT) :: IBUF1, IBUF2, IBUF3
REAL, DIMENSION(D%NIJT) :: ZZW, &
                           ZRDRYG_INIT, & !Initial dry growth rate of the graupeln
                           ZRWETG_INIT !Initial wet growth rate of the graupeln
REAL :: ZZW0D
INTEGER :: JIJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RG', 0, ZHOOK_HANDLE)
!
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!
!-------------------------------------------------------------------------------
!
!*       6.1    rain contact freezing
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRIT(JIJ)>ICED%XRTMIN(4) .AND. PRRT(JIJ)>ICED%XRTMIN(3) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRICFRRG(JIJ) = ICEP%XICFRR*PRIT(JIJ)                & ! RICFRRG
                                *PLBDAR(JIJ)**ICEP%XEXICFRR    &
                                *PRHODREF(JIJ)**(-ICED%XCEXVT)
      PRRCFRIG(JIJ) = ICEP%XRCFRI*PCIT(JIJ)                & ! RRCFRIG
                                * PLBDAR(JIJ)**ICEP%XEXRCFRI    &
                                * PRHODREF(JIJ)**(-ICED%XCEXVT-1.)
      IF(PARAMI%LCRFLIMIT) THEN
        !Comparison between heat to be released (to freeze rain) and heat sink (rain and ice temperature change)
        !ZZW0D is the proportion of process that can take place
        ZZW0D=MAX(0., MIN(1., (PRICFRRG(JIJ)*CST%XCI+PRRCFRIG(JIJ)*CST%XCL)*(CST%XTT-PT(JIJ)) / &
                              MAX(1.E-20, CST%XLVTT*PRRCFRIG(JIJ))))
        PRRCFRIG(JIJ) = ZZW0D * PRRCFRIG(JIJ) !Part of rain that can be freezed
        PRICFRR(JIJ) = (1.-ZZW0D) * PRICFRRG(JIJ) !Part of collected pristine ice converted to rain
        PRICFRRG(JIJ) = ZZW0D * PRICFRRG(JIJ) !Part of collected pristine ice that lead to graupel
      ELSE
        PRICFRR(JIJ) = 0.
      ENDIF
    ENDIF
  ELSE
    PRICFRRG(JIJ)=0.
    PRRCFRIG(JIJ)=0.
    PRICFRR(JIJ)=0.
  ENDIF
ENDDO
!$mnh_end_do()

!
!
!*       6.3    compute the graupel growth
!
! Wet and dry collection of rc and ri on graupel

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRGT(JIJ)>ICED%XRTMIN(6) .AND. PRCT(JIJ)>ICED%XRTMIN(2) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRG_TEND(JIJ, IRCDRYG)=PLBDAG(JIJ)**(ICED%XCXG-ICED%XDG-2.) * PRHODREF(JIJ)**(-ICED%XCEXVT)
      PRG_TEND(JIJ, IRCDRYG)=ICEP%XFCDRYG * PRCT(JIJ) * PRG_TEND(JIJ, IRCDRYG)
    ENDIF
  ELSE
    PRG_TEND(JIJ, IRCDRYG)=0.
  ENDIF

  IF(PRGT(JIJ)>ICED%XRTMIN(6) .AND. PRIT(JIJ)>ICED%XRTMIN(4) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRG_TEND(JIJ, IRIDRYG)=PLBDAG(JIJ)**(ICED%XCXG-ICED%XDG-2.) * PRHODREF(JIJ)**(-ICED%XCEXVT)
      PRG_TEND(JIJ, IRIDRYG)=ICEP%XFIDRYG*EXP(ICEP%XCOLEXIG*(PT(JIJ)-CST%XTT))*PRIT(JIJ)*PRG_TEND(JIJ, IRIDRYG)
      PRG_TEND(JIJ, IRIWETG)=PRG_TEND(JIJ, IRIDRYG) / (ICEP%XCOLIG*EXP(ICEP%XCOLEXIG*(PT(JIJ)-CST%XTT)))
    ENDIF
  ELSE
    PRG_TEND(JIJ, IRIDRYG)=0.
    PRG_TEND(JIJ, IRIWETG)=0.
  ENDIF
ENDDO
!$mnh_end_do()


! Wet and dry collection of rs on graupel (6.2.1)

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF (PRST(JIJ)>ICED%XRTMIN(5) .AND. PRGT(JIJ)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JIJ)) THEN
    GDRY(JIJ) = .TRUE.
    IF(PARAMI%LOCND2)THEN
      PRG_TEND(JIJ, IRSDRYG)=0.
      PRG_TEND(JIJ, IRSWETG)=0.
    ENDIF
  ELSE
    GDRY(JIJ) = .FALSE.
    PRG_TEND(JIJ, IRSDRYG)=0.
    PRG_TEND(JIJ, IRSWETG)=0.
  ENDIF
ENDDO
!$mnh_end_do()


IF(.NOT. LDSOFT) THEN
  CALL INTERP_MICRO_2D(D, PLBDAG, PLBDAS, ICEP%NDRYLBDAG, ICEP%NDRYLBDAS, &
                       &ICEP%XDRYINTP1G, ICEP%XDRYINTP2G, ICEP%XDRYINTP1S, ICEP%XDRYINTP2S, &
                       &PARAMI%LPACK_INTERP, GDRY, IBUF1, IBUF2, IBUF3, ZBUF1, ZBUF2, ZBUF3, &
                       &IGDRY, &
                       &ICEP%XKER_SDRYG, ZZW)
  IF(IGDRY>0)THEN

    DO JIJ=D%NIJB, D%NIJE
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        IF(.NOT. PARAMI%LOCND2) THEN
          ! Here, OCND2 option is only used for disregarding collision snow-graupel
          ! which according to Greg Thompson can be neglected.
          IF (GDRY(JIJ)) THEN
            PRG_TEND(JIJ, IRSWETG)=ICEP%XFSDRYG*ZZW(JIJ)                         & ! RSDRYG
                                        / ICEP%XCOLSG &
                        *(PLBDAS(JIJ)**(ICED%XCXS-ICED%XBS))*( PLBDAG(JIJ)**ICED%XCXG )    &
                        *(PRHODREF(JIJ)**(-ICED%XCEXVT-1.))                    &
                             *( ICEP%XLBSDRYG1/( PLBDAG(JIJ)**2              ) + &
                                ICEP%XLBSDRYG2/( PLBDAG(JIJ)   * PLBDAS(JIJ)   ) + &
                                ICEP%XLBSDRYG3/(               PLBDAS(JIJ)**2))
          END IF
        ENDIF
      ELSE
        IF (GDRY(JIJ)) THEN
          PRG_TEND(JIJ, IRSWETG)=ICEP%XFSDRYG*ZZW(JIJ)                         & ! RSDRYG
                                      / ICEP%XCOLSG &
                      *(PRST(JIJ))*( PLBDAG(JIJ)**ICED%XCXG )    &
                      *(PRHODREF(JIJ)**(-ICED%XCEXVT))                    &
                           *( ICEP%XLBSDRYG1/( PLBDAG(JIJ)**2              ) + &
                              ICEP%XLBSDRYG2/( PLBDAG(JIJ)   * PLBDAS(JIJ)   ) + &
                              ICEP%XLBSDRYG3/(               PLBDAS(JIJ)**2))
        END IF
      ENDIF
      IF (GDRY(JIJ)) THEN
        PRG_TEND(JIJ, IRSDRYG)=PRG_TEND(JIJ, IRSWETG)*ICEP%XCOLSG*EXP(ICEP%XCOLEXSG*(PT(JIJ)-CST%XTT))
      END IF
    END DO

  ENDIF
ENDIF
!
!*       6.2.6  accretion of raindrops on the graupeln
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF (PRRT(JIJ)>ICED%XRTMIN(3) .AND. PRGT(JIJ)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JIJ)) THEN
    GDRY(JIJ) = .TRUE.
  ELSE
    GDRY(JIJ) = .FALSE.
    PRG_TEND(JIJ, IRRDRYG)=0.
  ENDIF
ENDDO
!$mnh_end_do()

IF(.NOT. LDSOFT) THEN
  !
  CALL INTERP_MICRO_2D(D, PLBDAG, PLBDAR, ICEP%NDRYLBDAG, ICEP%NDRYLBDAR, &
                       &ICEP%XDRYINTP1G, ICEP%XDRYINTP2G, ICEP%XDRYINTP1R, ICEP%XDRYINTP2R, &
                       &PARAMI%LPACK_INTERP, GDRY, IBUF1, IBUF2, IBUF3, ZBUF1, ZBUF2, ZBUF3, &
                       &IGDRY, &
                       &ICEP%XKER_RDRYG, ZZW)
  IF(IGDRY>0) THEN

    DO JIJ=D%NIJB, D%NIJE
      IF (GDRY(JIJ)) THEN
        PRG_TEND(JIJ, IRRDRYG) = ICEP%XFRDRYG*ZZW(JIJ)                    & ! RRDRYG
                        *( PLBDAR(JIJ)**(-4) )*( PLBDAG(JIJ)**ICED%XCXG ) &
                                 *( PRHODREF(JIJ)**(-ICED%XCEXVT-1.) )   &
                      *( ICEP%XLBRDRYG1/( PLBDAG(JIJ)**2              ) + &
                         ICEP%XLBRDRYG2/( PLBDAG(JIJ)   * PLBDAR(JIJ)   ) + &
                         ICEP%XLBRDRYG3/(               PLBDAR(JIJ)**2) )
      END IF
    END DO

  ENDIF
ENDIF


!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  ZRDRYG_INIT(JIJ)=PRG_TEND(JIJ, IRCDRYG)+PRG_TEND(JIJ, IRIDRYG)+ &
                 &PRG_TEND(JIJ, IRSDRYG)+PRG_TEND(JIJ, IRRDRYG)
ENDDO
!$mnh_end_do()


!Freezing rate and growth mode

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRGT(JIJ)>ICED%XRTMIN(6) .AND. LDCOMPUTE(JIJ)) THEN
    !Freezing rate
    IF(.NOT. LDSOFT) THEN
      PRG_TEND(JIJ, IFREEZ1)=PRVT(JIJ)*PPRES(JIJ)/(CST%XEPSILO+PRVT(JIJ)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRG_TEND(JIJ, IFREEZ1)=MIN(PRG_TEND(JIJ, IFREEZ1), EXP(CST%XALPI-CST%XBETAI/PT(JIJ)-CST%XGAMI*LOG(PT(JIJ)))) ! min(ev, es_i(T))
      ENDIF
      PRG_TEND(JIJ, IFREEZ1)=PKA(JIJ)*(CST%XTT-PT(JIJ)) +                              &
               (PDV(JIJ)*(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JIJ)-CST%XTT)) &
                             *(CST%XESTT-PRG_TEND(JIJ, IFREEZ1))/(CST%XRV*PT(JIJ))           )
      PRG_TEND(JIJ, IFREEZ1)=PRG_TEND(JIJ, IFREEZ1)* ( ICEP%X0DEPG*       PLBDAG(JIJ)**ICEP%XEX0DEPG +     &
                             ICEP%X1DEPG*PCJ(JIJ)*PLBDAG(JIJ)**ICEP%XEX1DEPG )/ &
                            ( PRHODREF(JIJ)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JIJ))) )
      PRG_TEND(JIJ, IFREEZ2)=(PRHODREF(JIJ)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PT(JIJ)))   ) / &
                            ( PRHODREF(JIJ)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JIJ))) )
    ENDIF
    ZRWETG_INIT(JIJ)=MAX(PRG_TEND(JIJ, IRIWETG)+PRG_TEND(JIJ, IRSWETG), &
                       &MAX(0., PRG_TEND(JIJ, IFREEZ1) + &
                       &        PRG_TEND(JIJ, IFREEZ2) * (PRG_TEND(JIJ, IRIWETG)+PRG_TEND(JIJ, IRSWETG))))

    !Growth mode
    LDWETG(JIJ)=MAX(0., ZRWETG_INIT(JIJ)-PRG_TEND(JIJ, IRIWETG)-PRG_TEND(JIJ, IRSWETG))<= &
              &MAX(0., ZRDRYG_INIT(JIJ)-PRG_TEND(JIJ, IRIDRYG)-PRG_TEND(JIJ, IRSDRYG))

    IF(PARAMI%LNULLWETG) THEN
      LDWETG(JIJ) = LDWETG(JIJ) .AND. ZRDRYG_INIT(JIJ)>0.
    ELSE
      LDWETG(JIJ) = LDWETG(JIJ) .AND. ZRWETG_INIT(JIJ)>0.
    ENDIF
    IF(.NOT. PARAMI%LWETGPOST) THEN
      LDWETG(JIJ) = LDWETG(JIJ) .AND. PT(JIJ)<CST%XTT
    ENDIF

    LLDRYG(JIJ)=PT(JIJ)<CST%XTT .AND. ZRDRYG_INIT(JIJ)>1.E-20 .AND. &
              &MAX(0., ZRWETG_INIT(JIJ)-PRG_TEND(JIJ, IRIWETG)-PRG_TEND(JIJ, IRSWETG))>&
              &MAX(0., ZRDRYG_INIT(JIJ)-PRG_TEND(JIJ, IRIDRYG)-PRG_TEND(JIJ, IRSDRYG))
  ELSE
    PRG_TEND(JIJ, IFREEZ1)=0.
    PRG_TEND(JIJ, IFREEZ2)=0.
    ZRWETG_INIT(JIJ)=0.
    LDWETG(JIJ)=.FALSE.
    LLDRYG(JIJ)=.FALSE.
  ENDIF
ENDDO
!$mnh_end_do()


! Part of ZRWETG to be converted into hail
! Graupel can be produced by other processes instantaneously (inducing a mixing ratio change, PRGSI_MR) or
! as a tendency (PRWETGH)
IF(KRR==7) THEN

  DO JIJ=D%NIJB, D%NIJE
    IF (LDWETG(JIJ)) THEN
      !assume a linear percent of conversion of produced graupel into hail
      PRWETGH(JIJ)=(MAX(0., PRGSI(JIJ)+PRICFRRG(JIJ)+PRRCFRIG(JIJ))+ZRWETG_INIT(JIJ))*&
                      &ZRDRYG_INIT(JIJ)/(ZRWETG_INIT(JIJ)+ZRDRYG_INIT(JIJ))
      PRWETGH_MR(JIJ)=MAX(0., PRGSI_MR(JIJ))*ZRDRYG_INIT(JIJ)/(ZRWETG_INIT(JIJ)+ZRDRYG_INIT(JIJ))
    ELSE
      PRWETGH(JIJ)=0.
      PRWETGH_MR(JIJ)=0.
    END IF
  END DO

ELSE

  PRWETGH(:)=0.
  PRWETGH_MR(:)=0.

ENDIF


!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  !Aggregated minus collected
  IF(LDWETG(JIJ)) THEN
    PRRWETG(JIJ)=-(PRG_TEND(JIJ, IRIWETG)+PRG_TEND(JIJ, IRSWETG)+&
                 &PRG_TEND(JIJ, IRCDRYG)-ZRWETG_INIT(JIJ))
    PRCWETG(JIJ)=PRG_TEND(JIJ, IRCDRYG)
    PRIWETG(JIJ)=PRG_TEND(JIJ, IRIWETG)
    PRSWETG(JIJ)=PRG_TEND(JIJ, IRSWETG)
  ELSE
    PRRWETG(JIJ)=0.
    PRCWETG(JIJ)=0.
    PRIWETG(JIJ)=0.
    PRSWETG(JIJ)=0.
  ENDIF

  IF(LLDRYG(JIJ)) THEN
    PRCDRYG(JIJ)=PRG_TEND(JIJ, IRCDRYG)
    PRRDRYG(JIJ)=PRG_TEND(JIJ, IRRDRYG)
    PRIDRYG(JIJ)=PRG_TEND(JIJ, IRIDRYG)
    PRSDRYG(JIJ)=PRG_TEND(JIJ, IRSDRYG)
  ELSE
    PRCDRYG(JIJ)=0.
    PRRDRYG(JIJ)=0.
    PRIDRYG(JIJ)=0.
    PRSDRYG(JIJ)=0.
  ENDIF
ENDDO
!$mnh_end_do()

!
!*       6.5    Melting of the graupeln
!

!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRGT(JIJ)>ICED%XRTMIN(6) .AND. PT(JIJ)>CST%XTT .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRGMLTR(JIJ)=PRVT(JIJ)*PPRES(JIJ)/(CST%XEPSILO+PRVT(JIJ)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRGMLTR(JIJ)=MIN(PRGMLTR(JIJ), EXP(CST%XALPW-CST%XBETAW/PT(JIJ)-CST%XGAMW*LOG(PT(JIJ)))) ! min(ev, es_w(T))
      ENDIF
      PRGMLTR(JIJ)=PKA(JIJ)*(CST%XTT-PT(JIJ)) +                                 &
                  PDV(JIJ)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(JIJ) - CST%XTT )) &
                         *(CST%XESTT-PRGMLTR(JIJ))/(CST%XRV*PT(JIJ)) 
      PRGMLTR(JIJ)=MAX(0., (-PRGMLTR(JIJ)*                     &
                           (ICEP%X0DEPG*       PLBDAG(JIJ)**ICEP%XEX0DEPG +     &
                            ICEP%X1DEPG*PCJ(JIJ)*PLBDAG(JIJ)**ICEP%XEX1DEPG) -   &
                           (PRG_TEND(JIJ, IRCDRYG)+PRG_TEND(JIJ, IRRDRYG)) *     &
                           (PRHODREF(JIJ)*CST%XCL*(CST%XTT-PT(JIJ)))) /    &
                          ( PRHODREF(JIJ)*CST%XLMTT ) )
    ENDIF
  ELSE
    PRGMLTR(JIJ)=0.
  ENDIF
ENDDO
!$mnh_end_do()

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
