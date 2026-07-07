!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_FAST_RS
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RS(CST, PARAMI, ICEP, ICED, D, LDSOFT, LDCOMPUTE, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, PCOLF, &
                       &PLBDAR, PLBDAS, &
                       &PT,  PRVT, PRCT, PRRT, PRST, &
                       &PRIAGGS, &
                       &PRCRIMSS, PRCRIMSG, PRSRIMCG, &
                       &PRRACCSS, PRRACCSG, PRSACCRG, PRSMLTG, &
                       &PRCMLTSR, &
                       &PRS_TEND)

!$ACDC singlecolumn --dummy

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
!  K.I Ivarsson   02/2023: Possibility to use alternative collision factor
!!                            and tunings for e.g. better forecasting supercooled rain
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
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLVFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCOLF    ! Collision factor cloud water - snow/graupel
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRIAGGS  ! r_i aggregation on r_s
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRCRIMSS ! Cloud droplet riming of the aggregates
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRCRIMSG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRSRIMCG ! Cloud droplet riming of the aggregates
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRRACCSS ! Rain accretion onto the aggregates
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRRACCSG ! Rain accretion onto the aggregates
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRSACCRG ! Rain accretion onto the aggregates
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRSMLTG  ! Conversion-Melting of the aggregates
REAL, DIMENSION(D%NIJT),      INTENT(INOUT) :: PRCMLTSR ! Cloud droplet collection onto aggregates by positive temperature
REAL, DIMENSION(D%NIJT, 8),   INTENT(INOUT) :: PRS_TEND ! Individual tendencies
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCRIMS=1, IRCRIMSS=2, IRSRIMCG=3, IRRACCS=4, IRRACCSS=5, IRSACCRG=6, &
                    & IFREEZ1=7, IFREEZ2=8
LOGICAL, DIMENSION(D%NIJT) :: GRIM, GACC
INTEGER :: IGRIM, IGACC
INTEGER, DIMENSION(D%NIJT) :: IBUF1, IBUF2, IBUF3
REAL, DIMENSION(D%NIJT) :: ZBUF1, ZBUF2, ZBUF3
REAL, DIMENSION(D%NIJT) :: ZZW, ZZW1, ZZW2, ZZW3, ZFREEZ_RATE
INTEGER :: JIJ
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
#ifdef MNH_COMPILER_CCE
!$mnh_undef(LOOP)
#endif
!
!$acc kernels
!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRST(JIJ)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRS_TEND(JIJ, IFREEZ1)=PRVT(JIJ)*PPRES(JIJ)/(CST%XEPSILO+PRVT(JIJ)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRS_TEND(JIJ, IFREEZ1)=MIN(PRS_TEND(JIJ, IFREEZ1), EXP(CST%XALPI-CST%XBETAI/PT(JIJ)-CST%XGAMI*LOG(PT(JIJ)))) ! min(ev, es_i(T))
      ENDIF
      PRS_TEND(JIJ, IFREEZ1)=PKA(JIJ)*(CST%XTT-PT(JIJ)) +                              &
                           &(PDV(JIJ)*(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(JIJ)-CST%XTT)) &
                           &*(CST%XESTT-PRS_TEND(JIJ, IFREEZ1))/(CST%XRV*PT(JIJ))           )
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        PRS_TEND(JIJ, IFREEZ1)=PRS_TEND(JIJ, IFREEZ1)* (ICEP%X0DEPS*       PLBDAS(JIJ)**ICEP%XEX0DEPS +     &
                             &                        ICEP%X1DEPS*PCJ(JIJ)*PLBDAS(JIJ)**ICEP%XEX1DEPS )/ &
                             &(PRHODREF(JIJ)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JIJ))))
      ELSE
        PRS_TEND(JIJ, IFREEZ1)=PRS_TEND(JIJ, IFREEZ1)* PRST(JIJ) *(ICEP%X0DEPS*       PLBDAS(JIJ)**ICEP%XEX0DEPS + &
                             &                        ICEP%X1DEPS*PCJ(JIJ)*PLBDAS(JIJ)**(ICED%XBS+ICEP%XEX1DEPS )* &
           (1+0.5*(ICED%XFVELOS/PLBDAS(JIJ))**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEX1DEPS/ICED%XALPHAS))/ &
                             &(PRHODREF(JIJ)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JIJ))))
      ENDIF
      PRS_TEND(JIJ, IFREEZ2)=(PRHODREF(JIJ)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PT(JIJ)))   ) / &
                           &(PRHODREF(JIJ)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(JIJ))))
    ENDIF
    !We must agregate, at least, the cold species
    !And we are only interested by the freezing rate of liquid species
    ZFREEZ_RATE(JIJ)=MAX(0., MAX(0., PRS_TEND(JIJ, IFREEZ1) + &
                                    &PRS_TEND(JIJ, IFREEZ2) * PRIAGGS(JIJ)) - &
                            PRIAGGS(JIJ))
  ELSE
    PRS_TEND(JIJ, IFREEZ1)=0.
    PRS_TEND(JIJ, IFREEZ2)=0.
    ZFREEZ_RATE(JIJ)=0.
  ENDIF
ENDDO
!$mnh_end_do()
!$acc end kernels
!
!*       5.1    cloud droplet riming of the aggregates
!
!$acc kernels
!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF (PRCT(JIJ)>ICED%XRTMIN(2) .AND. PRST(JIJ)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. ICEP%LNEWCOEFF) THEN
      ZZW(JIJ) = PLBDAS(JIJ)
    ELSE
      ZZW(JIJ) = (PLBDAS(JIJ)**ICED%XALPHAS + ICED%XFVELOS**ICED%XALPHAS)**(1./ICED%XALPHAS)
    ENDIF
    GRIM(JIJ) = .TRUE.
  ELSE
    GRIM(JIJ) = .FALSE.
    PRS_TEND(JIJ, IRCRIMS)=0.
    PRS_TEND(JIJ, IRCRIMSS)=0.
    PRS_TEND(JIJ, IRSRIMCG)=0.
  ENDIF
ENDDO
!$mnh_end_do()
!$acc end kernels
!
! Collection of cloud droplets by snow: this rate is used for riming (T<0) and for conversion/melting (T>0)
IF(.NOT. LDSOFT) THEN
  CALL INTERP_MICRO_1D(D, ZZW, ICEP%NGAMINC, ICEP%XRIMINTP1, ICEP%XRIMINTP2, &
                           PARAMI%LPACK_INTERP, GRIM, IBUF1, IBUF2, ZBUF1, ZBUF2, &
                           IGRIM, &
                           ICEP%XGAMINC_RIM1, ZZW1, ICEP%XGAMINC_RIM2, ZZW2, ICEP%XGAMINC_RIM4, ZZW3)
  IF(IGRIM>0) THEN
    !
    !        5.1.4  riming of the small sized aggregates
    !
!$acc kernels
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    IF(.NOT. ICEP%LNEWCOEFF) THEN
      WHERE (GRIM(D%NIJB:D%NIJE))
        PRS_TEND(D%NIJB:D%NIJE, IRCRIMSS) = ICEP%XCRIMSS * ZZW1(D%NIJB:D%NIJE) * PRCT(D%NIJB:D%NIJE) * PCOLF(D%NIJB:D%NIJE) & ! RCRIMSS
                                        * PLBDAS(D%NIJB:D%NIJE)**ICEP%XEXCRIMSS &
                                        * PRHODREF(D%NIJB:D%NIJE)**(-ICED%XCEXVT)
      END WHERE
    ELSE
      WHERE (GRIM(D%NIJB:D%NIJE))
        PRS_TEND(D%NIJB:D%NIJE, IRCRIMSS) = ICEP%XCRIMSS * ZZW1(D%NIJB:D%NIJE) * PRCT(D%NIJB:D%NIJE) * PCOLF(D%NIJB:D%NIJE) & ! RCRIMSS
                                        * PRST(D%NIJB:D%NIJE)*(1+(ICED%XFVELOS/PLBDAS(D%NIJB:D%NIJE))**ICED%XALPHAS) &
                                          **(-ICED%XNUS+ICEP%XEXCRIMSS/ICED%XALPHAS) &
                                        * PRHODREF(D%NIJB:D%NIJE)**(-ICED%XCEXVT+1.) &
                                        * (PLBDAS(D%NIJB:D%NIJE)) ** (ICEP%XEXCRIMSS+ICED%XBS)
      END WHERE
    ENDIF
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
    !
    !        5.1.6  riming-conversion of the large sized aggregates into graupeln
    !
    !
!$acc kernels
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    IF(.NOT. ICEP%LNEWCOEFF) THEN
      WHERE(GRIM(D%NIJB:D%NIJE))
        PRS_TEND(D%NIJB:D%NIJE, IRCRIMS)=ICEP%XCRIMSG * PRCT(D%NIJB:D%NIJE) * PCOLF(D%NIJB:D%NIJE)               & ! RCRIMS
                                     * PLBDAS(D%NIJB:D%NIJE)**ICEP%XEXCRIMSG  &
                                     * PRHODREF(D%NIJB:D%NIJE)**(-ICED%XCEXVT)
      END WHERE
    ELSE
      WHERE(GRIM(D%NIJB:D%NIJE))
        PRS_TEND(D%NIJB:D%NIJE, IRCRIMS)=ICEP%XCRIMSG * PRCT(D%NIJB:D%NIJE) * PCOLF(D%NIJB:D%NIJE)               & ! RCRIMS
                                     * PRST(D%NIJB:D%NIJE)*(1+(ICED%XFVELOS/PLBDAS(D%NIJB:D%NIJE))**(ICED%XALPHAS)) &
                                       **(-ICED%XNUS+ICEP%XEXCRIMSG/ICED%XALPHAS) &
                                     * PRHODREF(D%NIJB:D%NIJE)**(-ICED%XCEXVT+1.) &
                                     * PLBDAS(D%NIJB:D%NIJE)**(ICED%XBS+ICEP%XEXCRIMSG)
      END WHERE
    ENDIF
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
    IF(PARAMI%CSNOWRIMING=='M90 ')THEN
      !Murakami 1990
!$acc kernels
      !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
      WHERE(GRIM(D%NIJB:D%NIJE))
        ZZW(D%NIJB:D%NIJE) = PRS_TEND(D%NIJB:D%NIJE, IRCRIMS) - PRS_TEND(D%NIJB:D%NIJE, IRCRIMSS) ! RCRIMSG
      END WHERE
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        WHERE(GRIM(D%NIJB:D%NIJE))
          PRS_TEND(D%NIJB:D%NIJE, IRSRIMCG)=ICEP%XSRIMCG * PLBDAS(D%NIJB:D%NIJE)**ICEP%XEXSRIMCG*(1.0-ZZW2(D%NIJB:D%NIJE))
          PRS_TEND(D%NIJB:D%NIJE, IRSRIMCG)=ZZW(D%NIJB:D%NIJE)*PRS_TEND(D%NIJB:D%NIJE, IRSRIMCG)/ &
                         MAX(1.E-20, &
                             ICEP%XSRIMCG3*ICEP%XSRIMCG2*PLBDAS(D%NIJB:D%NIJE)**ICEP%XEXSRIMCG2*(1.-ZZW3(D%NIJB:D%NIJE)) - &
                             ICEP%XSRIMCG3*PRS_TEND(D%NIJB:D%NIJE, IRSRIMCG))
        END WHERE
      ELSE
        WHERE(GRIM(D%NIJB:D%NIJE))
          PRS_TEND(D%NIJB:D%NIJE, IRSRIMCG)=ICEP%XSRIMCG * PRST(D%NIJB:D%NIJE)*PRHODREF(D%NIJB:D%NIJE) &
                                                   * PLBDAS(D%NIJB:D%NIJE)**(ICEP%XEXSRIMCG+ICED%XBS)*(1.0-ZZW2(D%NIJB:D%NIJE))
          PRS_TEND(D%NIJB:D%NIJE, IRSRIMCG)=ZZW(D%NIJB:D%NIJE)*PRS_TEND(D%NIJB:D%NIJE, IRSRIMCG)/ &
                         MAX(1.E-20, &
                             ICEP%XSRIMCG3*ICEP%XSRIMCG2*PRST(D%NIJB:D%NIJE)*PRHODREF(D%NIJB:D%NIJE) &
                                                        *PLBDAS(D%NIJB:D%NIJE)**ICEP%XEXSRIMCG2*(1.-ZZW3(D%NIJB:D%NIJE)) - &
                             ICEP%XSRIMCG3*PRS_TEND(D%NIJB:D%NIJE, IRSRIMCG))

        END WHERE
      ENDIF
      !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
    ELSE
!$acc kernels
      PRS_TEND(:, IRSRIMCG)=0.
!$acc end kernels
    END IF
  ENDIF
ENDIF
!
!$acc kernels
!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  ! More restrictive RIM mask to be used for riming by negative temperature only
  IF(GRIM(JIJ) .AND. PT(JIJ)<CST%XTT) THEN
    PRCRIMSS(JIJ)=MIN(ZFREEZ_RATE(JIJ), PRS_TEND(JIJ, IRCRIMSS))
    ZFREEZ_RATE(JIJ)=MAX(0., ZFREEZ_RATE(JIJ)-PRCRIMSS(JIJ))
    ZZW0D = MIN(1., ZFREEZ_RATE(JIJ) / MAX(1.E-20, PRS_TEND(JIJ, IRCRIMS) - PRCRIMSS(JIJ))) ! proportion we are able to freeze
    PRCRIMSG(JIJ) = ZZW0D * MAX(0., PRS_TEND(JIJ, IRCRIMS) - PRCRIMSS(JIJ)) ! RCRIMSG
    ZFREEZ_RATE(JIJ)=MAX(0., ZFREEZ_RATE(JIJ)-PRCRIMSG(JIJ))
    PRSRIMCG(JIJ) = ZZW0D * PRS_TEND(JIJ, IRSRIMCG)

    PRSRIMCG(JIJ) = PRSRIMCG(JIJ) * MAX(0., -SIGN(1., -PRCRIMSG(JIJ)))
    PRCRIMSG(JIJ)=MAX(0., PRCRIMSG(JIJ))
  ELSE
    PRCRIMSS(JIJ)=0.
    PRCRIMSG(JIJ)=0.
    PRSRIMCG(JIJ)=0.
  ENDIF
ENDDO
!$mnh_end_do()
!$acc end kernels
!
!*       5.2    rain accretion onto the aggregates
!
!$acc kernels
!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF (PRRT(JIJ)>ICED%XRTMIN(3) .AND. PRST(JIJ)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JIJ)) THEN
    GACC(JIJ) = .TRUE.
  ELSE
    GACC(JIJ) = .FALSE.
    PRS_TEND(JIJ, IRRACCS)=0.
    PRS_TEND(JIJ, IRRACCSS)=0.
    PRS_TEND(JIJ, IRSACCRG)=0.
  END IF
ENDDO
!$mnh_end_do()
!$acc end kernels
IF(.NOT. LDSOFT) THEN
!$acc kernels
  PRS_TEND(:, IRRACCS)=0.
  PRS_TEND(:, IRRACCSS)=0.
  PRS_TEND(:, IRSACCRG)=0.
!$acc end kernels
  CALL INTERP_MICRO_2D(D, PLBDAS, PLBDAR, ICEP%NACCLBDAS, ICEP%NACCLBDAR, &
                      &ICEP%XACCINTP1S, ICEP%XACCINTP2S, ICEP%XACCINTP1R, ICEP%XACCINTP2R,&
                      &PARAMI%LPACK_INTERP, GACC, IBUF1, IBUF2, IBUF3, ZBUF1, ZBUF2, ZBUF3, &
                      &IGACC, &
                      &ICEP%XKER_RACCSS, ZZW1, ICEP%XKER_RACCS, ZZW2, ICEP%XKER_SACCRG, ZZW3)
  IF(IGACC>0)THEN
    !        5.2.4  raindrop accretion on the small sized aggregates
    !
!$acc kernels
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    IF(.NOT. ICEP%LNEWCOEFF) THEN
      WHERE(GACC(D%NIJB:D%NIJE))
        ZZW(D%NIJB:D%NIJE) =                                                        & !! coef of RRACCS
              ICEP%XFRACCSS*( PLBDAS(D%NIJB:D%NIJE)**ICED%XCXS )*( PRHODREF(D%NIJB:D%NIJE)**(-ICED%XCEXVT-1.) ) &
         *( ICEP%XLBRACCS1/((PLBDAS(D%NIJB:D%NIJE)**2)               ) +                  &
            ICEP%XLBRACCS2/( PLBDAS(D%NIJB:D%NIJE)    * PLBDAR(D%NIJB:D%NIJE)    ) +                  &
            ICEP%XLBRACCS3/(               (PLBDAR(D%NIJB:D%NIJE)**2)) )/PLBDAR(D%NIJB:D%NIJE)**4
      END WHERE
    ELSE
      WHERE(GACC(D%NIJB:D%NIJE))
        ZZW(D%NIJB:D%NIJE) =                                                        & !! coef of RRACCS
              ICEP%XFRACCSS*( PRST(D%NIJB:D%NIJE)*PLBDAS(D%NIJB:D%NIJE)**ICED%XBS )*( PRHODREF(D%NIJB:D%NIJE)**(-ICED%XCEXVT) ) &
         *( ICEP%XLBRACCS1/((PLBDAS(D%NIJB:D%NIJE)**2)               ) +                  &
            ICEP%XLBRACCS2/( PLBDAS(D%NIJB:D%NIJE)    * PLBDAR(D%NIJB:D%NIJE)    ) +                  &
            ICEP%XLBRACCS3/(               (PLBDAR(D%NIJB:D%NIJE)**2)) )/PLBDAR(D%NIJB:D%NIJE)**4
      END WHERE
    ENDIF
    WHERE(GACC(D%NIJB:D%NIJE))
      PRS_TEND(D%NIJB:D%NIJE, IRRACCSS) =ZZW1(D%NIJB:D%NIJE)*ZZW(D%NIJB:D%NIJE)
    END WHERE
    IF(PARAMI%LOCND2) THEN
      WHERE(GACC(D%NIJB:D%NIJE))
        PRS_TEND(D%NIJB:D%NIJE, IRRACCSS) = PRS_TEND(D%NIJB:D%NIJE, IRRACCSS) * ICEP%XFRMIN(7)
      END WHERE
    ENDIF
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
    !
!$acc kernels
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    WHERE(GACC(D%NIJB:D%NIJE))
      PRS_TEND(D%NIJB:D%NIJE, IRRACCS) = ZZW2(D%NIJB:D%NIJE)*ZZW(D%NIJB:D%NIJE)
    END WHERE
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
    !
    !        5.2.6  raindrop accretion-conversion of the large sized aggregates
    !               into graupeln
    !
!$acc kernels
    !$mnh_expand_where(JIJ=D%NIJB:D%NIJE)
    IF(.NOT. ICEP%LNEWCOEFF) THEN
      WHERE(GACC(D%NIJB:D%NIJE) .AND. (.NOT. PARAMI%LOCND2 .OR. PRST(:)>ICEP%XFRMIN(1) ))
        PRS_TEND(D%NIJB:D%NIJE, IRSACCRG) = ICEP%XFSACCRG*ZZW3(D%NIJB:D%NIJE)*                    & ! RSACCRG
            ( PLBDAS(D%NIJB:D%NIJE)**(ICED%XCXS-ICED%XBS) )*( PRHODREF(D%NIJB:D%NIJE)**(-ICED%XCEXVT-1.) ) &
           *( ICEP%XLBSACCR1/((PLBDAR(D%NIJB:D%NIJE)**2)               ) +           &
              ICEP%XLBSACCR2/( PLBDAR(D%NIJB:D%NIJE)    * PLBDAS(D%NIJB:D%NIJE)    ) +           &
              ICEP%XLBSACCR3/(               (PLBDAS(D%NIJB:D%NIJE)**2)) )/PLBDAR(D%NIJB:D%NIJE)
      END WHERE
    ELSE
      WHERE(GACC(D%NIJB:D%NIJE) .AND. (.NOT. PARAMI%LOCND2 .OR. PRST(:)>ICEP%XFRMIN(1) ))
        PRS_TEND(D%NIJB:D%NIJE, IRSACCRG) = ICEP%XFSACCRG*ZZW3(D%NIJB:D%NIJE)*                    & ! RSACCRG
            ( PRST(D%NIJB:D%NIJE))*( PRHODREF(D%NIJB:D%NIJE)**(-ICED%XCEXVT) ) &
           *( ICEP%XLBSACCR1/((PLBDAR(D%NIJB:D%NIJE)**2)               ) +           &
              ICEP%XLBSACCR2/( PLBDAR(D%NIJB:D%NIJE)    * PLBDAS(D%NIJB:D%NIJE)    ) +           &
              ICEP%XLBSACCR3/(               (PLBDAS(D%NIJB:D%NIJE)**2)) )/PLBDAR(D%NIJB:D%NIJE)
      END WHERE
    ENDIF
    !$mnh_end_expand_where(JIJ=D%NIJB:D%NIJE)
!$acc end kernels
  ENDIF
ENDIF
!
!$acc kernels
!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  ! More restrictive ACC mask to be used for accretion by negative temperature only
  IF(GACC(JIJ) .AND. PT(JIJ)<CST%XTT) THEN
    PRRACCSS(JIJ)=MIN(ZFREEZ_RATE(JIJ), PRS_TEND(JIJ, IRRACCSS))
    ZFREEZ_RATE(JIJ)=MAX(0., ZFREEZ_RATE(JIJ)-PRRACCSS(JIJ))
    ZZW(JIJ) = MIN(1., ZFREEZ_RATE(JIJ) / MAX(1.E-20, PRS_TEND(JIJ, IRRACCS)-PRRACCSS(JIJ))) ! proportion we are able to freeze
    PRRACCSG(JIJ)=ZZW(JIJ) * MAX(0., PRS_TEND(JIJ, IRRACCS)-PRRACCSS(JIJ))
    ZFREEZ_RATE(JIJ) = MAX(0., ZFREEZ_RATE(JIJ)-PRRACCSG(JIJ))
    PRSACCRG(JIJ)=ZZW(JIJ) * PRS_TEND(JIJ, IRSACCRG)

    PRSACCRG(JIJ) = PRSACCRG(JIJ) * MAX(0., -SIGN(1., -PRRACCSG(JIJ)))
    PRRACCSG(JIJ)=MAX(0., PRRACCSG(JIJ))
  ELSE
    PRRACCSS(JIJ)=0.
    PRRACCSG(JIJ)=0.
    PRSACCRG(JIJ)=0.
  ENDIF
ENDDO
!$mnh_end_do()
!$acc end kernels
!
!
!*       5.3    Conversion-Melting of the aggregates
!
!$acc kernels
!$mnh_do_concurrent( JIJ=D%NIJB:D%NIJE )
DO JIJ=D%NIJB, D%NIJE
  IF(PRST(JIJ)>ICED%XRTMIN(5) .AND. PT(JIJ)>CST%XTT .AND. LDCOMPUTE(JIJ)) THEN
    IF(.NOT. LDSOFT) THEN
      PRSMLTG(JIJ)=PRVT(JIJ)*PPRES(JIJ)/(CST%XEPSILO+PRVT(JIJ)) ! Vapor pressure
      IF(PARAMI%LEVLIMIT) THEN
        PRSMLTG(JIJ)=MIN(PRSMLTG(JIJ), EXP(CST%XALPW-CST%XBETAW/PT(JIJ)-CST%XGAMW*LOG(PT(JIJ)))) ! min(ev, es_w(T))
      ENDIF
      PRSMLTG(JIJ)= PKA(JIJ)*(CST%XTT-PT(JIJ)) +                                 &
                  &(PDV(JIJ)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(JIJ) - CST%XTT )) &
                  & *(CST%XESTT-PRSMLTG(JIJ))/(CST%XRV*PT(JIJ))             )
      !
      ! compute RSMLT
      !
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        PRSMLTG(JIJ)  = ICEP%XFSCVMG*MAX(0., (-PRSMLTG(JIJ) * &
                   (ICEP%X0DEPS*       PLBDAS(JIJ)**ICEP%XEX0DEPS +     &
                   ICEP%X1DEPS*PCJ(JIJ)*PLBDAS(JIJ)**ICEP%XEX1DEPS)    &
                   -(PRS_TEND(JIJ, IRCRIMS) + PRS_TEND(JIJ, IRRACCS)) *       &
                   (PRHODREF(JIJ)*CST%XCL*(CST%XTT-PT(JIJ))) &
                   ) / (PRHODREF(JIJ)*CST%XLMTT))
      ELSE
        PRSMLTG(JIJ)  = ICEP%XFSCVMG*MAX(0., (-PRSMLTG(JIJ) * &
                   PRST(JIJ)*PRHODREF(JIJ) *    &
                   (ICEP%X0DEPS*       PLBDAS(JIJ)**(ICED%XBS+ICEP%XEX0DEPS) + &
                   ICEP%X1DEPS*PCJ(JIJ)* &
                   (1+0.5*(ICED%XFVELOS/PLBDAS(JIJ))**ICED%XALPHAS)**(-ICED%XNUS+ICEP%XEX1DEPS/ICED%XALPHAS) &
                   *PLBDAS(JIJ)**(ICED%XBS+ICEP%XEX1DEPS)) &
                   -(PRS_TEND(JIJ, IRCRIMS) + PRS_TEND(JIJ, IRRACCS)) *       &
                   (PRHODREF(JIJ)*CST%XCL*(CST%XTT-PT(JIJ))) &
                   ) / (PRHODREF(JIJ)*CST%XLMTT))
      ENDIF
      !
      ! note that RSCVMG = RSMLT*XFSCVMG but no heat is exchanged (at the rate RSMLT)
      ! because the graupeln produced by this process are still icy!!!
      !
      ! When T < XTT, rc is collected by snow (riming) to produce snow and graupel
      ! When T > XTT, if riming was still enabled, rc would produce snow and graupel with snow becomming graupel (conversion/melting) and graupel becomming rain (melting)
      ! To insure consistency when crossing T=XTT, rc collected with T>XTT must be transformed in rain.
      ! rc cannot produce iced species with a positive temperature but is still collected with a good efficiency by snow
      PRCMLTSR(JIJ) = PRS_TEND(JIJ, IRCRIMS) ! both species are liquid, no heat is exchanged
    ENDIF
  ELSE
    PRSMLTG(JIJ)=0.
    PRCMLTSR(JIJ)=0.
  ENDIF
ENDDO
!$mnh_end_do()
!$acc end kernels
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RS', 1, ZHOOK_HANDLE)
!
CONTAINS
!
INCLUDE "interp_micro.func.h"
!
END SUBROUTINE ICE4_FAST_RS
END MODULE MODE_ICE4_FAST_RS
