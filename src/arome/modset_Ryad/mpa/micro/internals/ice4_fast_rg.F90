SUBROUTINE ICE4_FAST_RG(KPROMA,KSIZE, LDSOFT, PCOMPUTE, KRR, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, PCIT, &
                       &PLBDAR, PLBDAS, PLBDAG, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, &
                       &PRGSI, PRGSI_MR, &
                       &PWETG, &
                       &PRICFRRG, PRRCFRIG, PRICFRR, PRCWETG, PRIWETG, PRRWETG, PRSWETG, &
                       &PRCDRYG, PRIDRYG, PRRDRYG, PRSDRYG, PRWETGH, PRWETGH_MR, PRGMLTR, &
                       &PRG_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RI, PA_RS, PA_RG, PA_RH, PB_RG, PB_RH)
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
!!     R. El Khatib 24-Aug-2021 Optimizations
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_RAIN_ICE_PARAM
USE MODD_RAIN_ICE_DESCR
USE MODD_PARAM_ICE, ONLY : LEVLIMIT, LNULLWETG, LWETGPOST, LCRFLIMIT
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KPROMA,KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
INTEGER,                      INTENT(IN)    :: KRR      ! Number of moist variable
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the raindrop  distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGSI    ! Graupel tendency by other processes
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGSI_MR ! Graupel mr change by other processes
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PWETG    ! 1. where graupel grows in wet mode, 0. elsewhere
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRICFRRG ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRCFRIG ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRICFRR  ! Rain contact freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSWETG  ! Graupel wet growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSDRYG  ! Graupel dry growth
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRWETGH  ! Conversion of graupel into hail
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRWETGH_MR ! Conversion of graupel into hail, mr change
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRGMLTR  ! Melting of the graupel
REAL, DIMENSION(KPROMA, 8),   INTENT(INOUT) :: PRG_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RH
!
!*       0.2  declaration of local variables
!
LOGICAL, DIMENSION(KSIZE) :: GDRY
REAL, DIMENSION(KSIZE) :: ZDRY, ZDRYG, ZMASK
INTEGER :: IGDRY
REAL, DIMENSION(KSIZE) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(KSIZE) :: IVEC1, IVEC2
REAL, DIMENSION(KSIZE) :: ZZW, &
                                   ZRDRYG_INIT, & !Initial dry growth rate of the graupeln
                                   ZRWETG_INIT !Initial wet growth rate of the graupeln
INTEGER :: JJ, JL
INTEGER, PARAMETER :: IRCDRYG=1
INTEGER, PARAMETER :: IRIDRYG=2
INTEGER, PARAMETER :: IRIWETG=3
INTEGER, PARAMETER :: IRSDRYG=4
INTEGER, PARAMETER :: IRSWETG=5
INTEGER, PARAMETER :: IRRDRYG=6
INTEGER, PARAMETER :: IFREEZ1=7
INTEGER, PARAMETER :: IFREEZ2=8

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
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(4)-PRIT(JL))) * & ! WHERE(PRIT(:)>XRTMIN(4))
           &MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JL))) * & ! WHERE(PRRT(:)>XRTMIN(3))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRICFRRG(JL)=ZMASK(JL) * PRICFRRG(JL)
    PRRCFRIG(JL)=ZMASK(JL) * PRRCFRIG(JL)
    PRICFRR(JL)=ZMASK(JL) * PRICFRR(JL)
  ENDDO
ELSE
  PRICFRRG(:)=0.
  PRRCFRIG(:)=0.
  WHERE(ZMASK(:)==1.)
    PRICFRRG(:) = XICFRR*PRIT(:)                & ! RICFRRG
                                 *PLBDAR(:)**XEXICFRR    &
                                 *PRHODREF(:)**(-XCEXVT)
    PRRCFRIG(:) = XRCFRI*PCIT(:)                & ! RRCFRIG
                                 * PLBDAR(:)**XEXRCFRI    &
                                 * PRHODREF(:)**(-XCEXVT-1.)
  END WHERE

  IF(LCRFLIMIT) THEN
    DO JL=1, KSIZE
      !Comparison between heat to be released (to freeze rain) and heat sink (rain and ice temperature change)
      !ZZW is the proportion of process that can take place
      ZZW(JL)=(1.-ZMASK(JL)) + & ! 1. outside of mask
              ZMASK(JL) * MAX(0., MIN(1., (PRICFRRG(JL)*XCI+PRRCFRIG(JL)*XCL)*(XTT-PT(JL)) / &
                                          MAX(1.E-20, XLVTT*PRRCFRIG(JL))))
    ENDDO
  ELSE
    ZZW(:)=1.
  ENDIF
  DO JL=1, KSIZE
    PRRCFRIG(JL) = ZZW(JL) * PRRCFRIG(JL) !Part of rain that can be freezed
    PRICFRR(JL) = (1.-ZZW(JL)) * PRICFRRG(JL) !Part of collected pristine ice converted to rain
    PRICFRRG(JL) = ZZW(JL) * PRICFRRG(JL) !Part of collected pristine ice that lead to graupel
  ENDDO
ENDIF
!
!
!*       6.3    compute the graupel growth
!
! Wet and dry collection of rc and ri on graupel
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(6)-PRGT(JL))) * & ! WHERE(PRGT(:)>XRTMIN(6))
           &MAX(0., -SIGN(1., XRTMIN(2)-PRCT(JL))) * & ! WHERE(PRCT(:)>XRTMIN(2))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRG_TEND(JL, IRCDRYG)=ZMASK(JL)*PRG_TEND(JL, IRCDRYG)
  ENDDO
ELSE
  ZZW(:)=0.
  WHERE(ZMASK(:)==1.)
    ZZW(:)=PLBDAG(:)**(XCXG-XDG-2.) * PRHODREF(:)**(-XCEXVT)
  END WHERE
  DO JL=1, KSIZE
    PRG_TEND(JL, IRCDRYG)=ZMASK(JL)*XFCDRYG * PRCT(JL) * ZZW(JL)
  ENDDO
ENDIF

DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(6)-PRGT(JL))) * & ! WHERE(PRGT(:)>XRTMIN(6))
           &MAX(0., -SIGN(1., XRTMIN(4)-PRIT(JL))) * & ! WHERE(PRIT(:)>XRTMIN(4))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRG_TEND(JL, IRIDRYG)=ZMASK(JL) * PRG_TEND(JL, IRIDRYG)
    PRG_TEND(JL, IRIWETG)=ZMASK(JL) * PRG_TEND(JL, IRIWETG)
  ENDDO
ELSE
  PRG_TEND(:, IRIDRYG)=0.
  PRG_TEND(:, IRIWETG)=0.
  WHERE(ZMASK(1:KSIZE)==1.)
    ZZW(1:KSIZE)=PLBDAG(1:KSIZE)**(XCXG-XDG-2.) * PRHODREF(1:KSIZE)**(-XCEXVT)
    PRG_TEND(1:KSIZE, IRIDRYG)=XFIDRYG*EXP(XCOLEXIG*(PT(1:KSIZE)-XTT))*PRIT(1:KSIZE)*ZZW(1:KSIZE)
    PRG_TEND(1:KSIZE, IRIWETG)=PRG_TEND(1:KSIZE, IRIDRYG) / (XCOLIG*EXP(XCOLEXIG*(PT(1:KSIZE)-XTT)))
  END WHERE
ENDIF

! Wet and dry collection of rs on graupel (6.2.1)
DO JL=1, KSIZE
  ZDRY(JL)=MAX(0., -SIGN(1., XRTMIN(5)-PRST(JL))) * & ! WHERE(PRST(:)>XRTMIN(5))
          &MAX(0., -SIGN(1., XRTMIN(6)-PRGT(JL))) * & ! WHERE(PRGT(:)>XRTMIN(6))
          &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRG_TEND(JL, IRSDRYG)=ZDRY(JL) * PRG_TEND(JL, IRSDRYG)
    PRG_TEND(JL, IRSWETG)=ZDRY(JL) * PRG_TEND(JL, IRSWETG)
  ENDDO
ELSE
  PRG_TEND(:, IRSDRYG)=0.
  PRG_TEND(:, IRSWETG)=0.
  GDRY(:)=ZDRY(:)==1.
  IGDRY=COUNT(GDRY(:))
  IF(IGDRY>0)THEN
    !
    !*       6.2.3  select the (PLBDAG,PLBDAS) couplet
    !
    ZVEC1(1:IGDRY)=PACK(PLBDAG(:), MASK=GDRY(:))
    ZVEC2(1:IGDRY)=PACK(PLBDAS(:), MASK=GDRY(:))
    !
    !*       6.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
    !               in the geometrical set of (Lbda_g,Lbda_s) couplet use to
    !               tabulate the SDRYG-kernel
    !
    ZVEC1(1:IGDRY)=MAX(1.00001, MIN(FLOAT(NDRYLBDAG)-0.00001,           &
                          XDRYINTP1G*LOG(ZVEC1(1:IGDRY))+XDRYINTP2G))
    IVEC1(1:IGDRY)=INT(ZVEC1(1:IGDRY) )
    ZVEC1(1:IGDRY)=ZVEC1(1:IGDRY)-FLOAT(IVEC1(1:IGDRY))
    !
    ZVEC2(1:IGDRY)=MAX(1.00001, MIN( FLOAT(NDRYLBDAS)-0.00001,           &
                          XDRYINTP1S*LOG(ZVEC2(1:IGDRY))+XDRYINTP2S))
    IVEC2(1:IGDRY)=INT(ZVEC2(1:IGDRY))
    ZVEC2(1:IGDRY)=ZVEC2(1:IGDRY)-FLOAT(IVEC2(1:IGDRY))
    !
    !*       6.2.5  perform the bilinear interpolation of the normalized
    !               SDRYG-kernel
    !
    DO JJ=1, IGDRY
      ZVEC3(JJ) =  (  XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * ZVEC1(JJ) &
                 - (  XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_SDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         *(ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:)=UNPACK(VECTOR=ZVEC3(1:IGDRY), MASK=GDRY(:), FIELD=0.0)
    !
    WHERE(GDRY(1:KSIZE))
      PRG_TEND(1:KSIZE, IRSWETG)=XFSDRYG*ZZW(1:KSIZE)                         & ! RSDRYG
                                    / XCOLSG &
                  *(PLBDAS(1:KSIZE)**(XCXS-XBS))*( PLBDAG(1:KSIZE)**XCXG )    &
                  *(PRHODREF(1:KSIZE)**(-XCEXVT-1.))                    &
                       *( XLBSDRYG1/( PLBDAG(1:KSIZE)**2              ) + &
                          XLBSDRYG2/( PLBDAG(1:KSIZE)   * PLBDAS(1:KSIZE)   ) + &
                          XLBSDRYG3/(               PLBDAS(1:KSIZE)**2))
      PRG_TEND(1:KSIZE, IRSDRYG)=PRG_TEND(1:KSIZE, IRSWETG)*XCOLSG*EXP(XCOLEXSG*(PT(1:KSIZE)-XTT))
    END WHERE
  ENDIF
ENDIF
!
!*       6.2.6  accretion of raindrops on the graupeln
!
DO JL=1, KSIZE
  ZDRY(JL)=MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JL))) * & ! WHERE(PRRT(:)>XRTMIN(3))
          &MAX(0., -SIGN(1., XRTMIN(6)-PRGT(JL))) * & ! WHERE(PRGT(:)>XRTMIN(6))
          &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRG_TEND(JL, IRRDRYG)=ZDRY(JL) * PRG_TEND(JL, IRRDRYG)
  ENDDO
ELSE
  GDRY(:)=ZDRY(:)==1.
  PRG_TEND(:, IRRDRYG)=0.
  IGDRY=COUNT(GDRY(:))
  !
  IF(IGDRY>0) THEN
    !
    !*       6.2.8  select the (PLBDAG,PLBDAR) couplet
    !
    ZVEC1(1:IGDRY)=PACK(PLBDAG(:), MASK=GDRY(:))
    ZVEC2(1:IGDRY)=PACK(PLBDAR(:), MASK=GDRY(:))
    !
    !*       6.2.9  find the next lower indice for the PLBDAG and for the PLBDAR
    !               in the geometrical set of (Lbda_g,Lbda_r) couplet use to
    !               tabulate the RDRYG-kernel
    !
    ZVEC1(1:IGDRY)=MAX(1.00001, MIN( FLOAT(NDRYLBDAG)-0.00001,           &
                          XDRYINTP1G*LOG(ZVEC1(1:IGDRY))+XDRYINTP2G))
    IVEC1(1:IGDRY)=INT(ZVEC1(1:IGDRY))
    ZVEC1(1:IGDRY)=ZVEC1(1:IGDRY)-FLOAT(IVEC1(1:IGDRY))
    !
    ZVEC2(1:IGDRY)=MAX(1.00001, MIN( FLOAT(NDRYLBDAR)-0.00001,           &
                          XDRYINTP1R*LOG(ZVEC2(1:IGDRY))+XDRYINTP2R))
    IVEC2(1:IGDRY)=INT(ZVEC2(1:IGDRY))
    ZVEC2(1:IGDRY)=ZVEC2(1:IGDRY)-FLOAT(IVEC2(1:IGDRY))
    !
    !*       6.2.10 perform the bilinear interpolation of the normalized
    !               RDRYG-kernel
    !
    DO JJ=1, IGDRY
      ZVEC3(JJ)= (  XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                 - (  XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RDRYG(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         *(ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:)=UNPACK(VECTOR=ZVEC3(1:IGDRY), MASK=GDRY, FIELD=0.)
    !
    WHERE(GDRY(1:KSIZE))
      PRG_TEND(1:KSIZE, IRRDRYG) = XFRDRYG*ZZW(1:KSIZE)                    & ! RRDRYG
                        *( PLBDAR(1:KSIZE)**(-4) )*( PLBDAG(1:KSIZE)**XCXG ) &
                               *( PRHODREF(1:KSIZE)**(-XCEXVT-1.) )   &
                    *( XLBRDRYG1/( PLBDAG(1:KSIZE)**2              ) + &
                       XLBRDRYG2/( PLBDAG(1:KSIZE)   * PLBDAR(1:KSIZE)   ) + &
                       XLBRDRYG3/(               PLBDAR(1:KSIZE)**2) )
    END WHERE
  ENDIF
ENDIF

DO JL=1, KSIZE
  ZRDRYG_INIT(JL)=PRG_TEND(JL, IRCDRYG)+PRG_TEND(JL, IRIDRYG)+ &
                 &PRG_TEND(JL, IRSDRYG)+PRG_TEND(JL, IRRDRYG)
ENDDO

!Freezing rate
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(6)-PRGT(JL))) * & ! WHERE(PRGT(:)>XRTMIN(6))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRG_TEND(JL, IFREEZ1)=ZMASK(JL) * PRG_TEND(JL, IFREEZ1)
    PRG_TEND(JL, IFREEZ2)=ZMASK(JL) * PRG_TEND(JL, IFREEZ2)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRG_TEND(JL, IFREEZ1)=ZMASK(JL) * PRVT(JL)*PPRES(JL)/(XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(LEVLIMIT) THEN
    WHERE(ZMASK(1:KSIZE)==1.)
      PRG_TEND(1:KSIZE, IFREEZ1)=MIN(PRG_TEND(1:KSIZE, IFREEZ1), EXP(XALPI-XBETAI/PT(1:KSIZE)-XGAMI*ALOG(PT(1:KSIZE)))) ! min(ev, es_i(T))
    END WHERE
  ENDIF
  PRG_TEND(:, IFREEZ2)=0.
  WHERE(ZMASK(1:KSIZE)==1.)
    PRG_TEND(1:KSIZE, IFREEZ1)=PKA(1:KSIZE)*(XTT-PT(1:KSIZE)) +                              &
             (PDV(1:KSIZE)*(XLVTT+(XCPV-XCL)*(PT(1:KSIZE)-XTT)) &
                           *(XESTT-PRG_TEND(1:KSIZE, IFREEZ1))/(XRV*PT(1:KSIZE))           )
    PRG_TEND(1:KSIZE, IFREEZ1)=PRG_TEND(1:KSIZE, IFREEZ1)* ( X0DEPG*       PLBDAG(1:KSIZE)**XEX0DEPG +     &
                           X1DEPG*PCJ(1:KSIZE)*PLBDAG(1:KSIZE)**XEX1DEPG )/ &
                          ( PRHODREF(1:KSIZE)*(XLMTT-XCL*(XTT-PT(1:KSIZE))) )
    PRG_TEND(1:KSIZE, IFREEZ2)=(PRHODREF(1:KSIZE)*(XLMTT+(XCI-XCL)*(XTT-PT(1:KSIZE)))   ) / &
                          ( PRHODREF(1:KSIZE)*(XLMTT-XCL*(XTT-PT(1:KSIZE))) )
  END WHERE
ENDIF
DO JL=1, KSIZE
  !We must agregate, at least, the cold species
  ZRWETG_INIT(JL)=ZMASK(JL) * MAX(PRG_TEND(JL, IRIWETG)+PRG_TEND(JL, IRSWETG), &
                                 &MAX(0., PRG_TEND(JL, IFREEZ1) + &
                                         &PRG_TEND(JL, IFREEZ2) * ( &
                     &PRG_TEND(JL, IRIWETG)+PRG_TEND(JL, IRSWETG) )))
ENDDO

!Growth mode
DO JL=1, KSIZE
  PWETG(JL) = ZMASK(JL) * & !
            & MAX(0., SIGN(1., MAX(0., ZRDRYG_INIT(JL)-PRG_TEND(JL, IRIDRYG)-PRG_TEND(JL, IRSDRYG)) - &
                              &MAX(0., ZRWETG_INIT(JL)-PRG_TEND(JL, IRIWETG)-PRG_TEND(JL, IRSWETG))))
ENDDO
IF(LNULLWETG) THEN
  DO JL=1, KSIZE
    PWETG(JL) = PWETG(JL) * MAX(0., -SIGN(1., -ZRDRYG_INIT(JL)))
  ENDDO
ELSE
  DO JL=1, KSIZE
    PWETG(JL) = PWETG(JL) * MAX(0., -SIGN(1., -ZRWETG_INIT(JL)))
  ENDDO
ENDIF
IF(.NOT. LWETGPOST) THEN
  DO JL=1, KSIZE
    PWETG(JL) = PWETG(JL) * MAX(0., -SIGN(1., PT(JL)-XTT))
  ENDDO
ENDIF
DO JL=1, KSIZE
  ZDRYG(JL) = ZMASK(JL) * & !
            & MAX(0., -SIGN(1., PT(JL)-XTT)) * & ! WHERE(PT(:)<XTT)
            & MAX(0., -SIGN(1., -ZRDRYG_INIT(JL))) * & ! WHERE(ZRDRYG_INIT(:)>0.)
            & MAX(0., -SIGN(1., MAX(0., ZRDRYG_INIT(JL)-PRG_TEND(JL, IRIDRYG)-PRG_TEND(JL, IRSDRYG)) - &
                               &MAX(0., ZRWETG_INIT(JL)-PRG_TEND(JL, IRIWETG)-PRG_TEND(JL, IRSWETG))))
ENDDO

! Part of ZRWETG to be converted into hail
! Graupel can be produced by other processes instantaneously (inducing a mixing ratio change, PRGSI_MR) or
! as a tendency (PRWETGH)
PRWETGH(:)=0.
PRWETGH_MR(:)=0.
IF(KRR==7) THEN
  WHERE(PWETG(:)==1.)
    !assume a linear percent of conversion of produced graupel into hail
    PRWETGH(:)=(MAX(0., PRGSI(:)+PRICFRRG(:)+PRRCFRIG(:))+ZRWETG_INIT(:))*ZRDRYG_INIT(:)/(ZRWETG_INIT(:)+ZRDRYG_INIT(:))
    PRWETGH_MR(:)=MAX(0., PRGSI_MR(:))*ZRDRYG_INIT(:)/(ZRWETG_INIT(:)+ZRDRYG_INIT(:))
  END WHERE
ENDIF

DO JL=1, KSIZE
  !Aggregated minus collected
  PRRWETG(JL)=-PWETG(JL) * (PRG_TEND(JL, IRIWETG)+PRG_TEND(JL, IRSWETG)+&
                           &PRG_TEND(JL, IRCDRYG)-ZRWETG_INIT(JL))
  PRCWETG(JL)=PWETG(JL) * PRG_TEND(JL, IRCDRYG)
  PRIWETG(JL)=PWETG(JL) * PRG_TEND(JL, IRIWETG)
  PRSWETG(JL)=PWETG(JL) * PRG_TEND(JL, IRSWETG)

  PRCDRYG(JL)=ZDRYG(JL) * PRG_TEND(JL, IRCDRYG)
  PRRDRYG(JL)=ZDRYG(JL) * PRG_TEND(JL, IRRDRYG)
  PRIDRYG(JL)=ZDRYG(JL) * PRG_TEND(JL, IRIDRYG)
  PRSDRYG(JL)=ZDRYG(JL) * PRG_TEND(JL, IRSDRYG)

ENDDO
!
!*       6.5    Melting of the graupeln
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(6)-PRGT(JL))) * & ! WHERE(PRGT(:)>XRTMIN(6))
           &MAX(0., -SIGN(1., XTT-PT(JL))) * & ! WHERE(PT(:)>XTT)
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRGMLTR(JL)=ZMASK(JL) * PRGMLTR(JL)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRGMLTR(JL)=ZMASK(JL) * PRVT(JL)*PPRES(JL)/(XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(LEVLIMIT) THEN
    WHERE(ZMASK(:)==1.)
      PRGMLTR(:)=MIN(PRGMLTR(:), EXP(XALPW-XBETAW/PT(:)-XGAMW*ALOG(PT(:)))) ! min(ev, es_w(T))
    END WHERE
  ENDIF
  DO JL=1, KSIZE
    PRGMLTR(JL)=ZMASK(JL) * (PKA(JL)*(XTT-PT(JL)) +                                 &
               ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PT(JL) - XTT )) &
                           *(XESTT-PRGMLTR(JL))/(XRV*PT(JL))             ))
  ENDDO
  WHERE(ZMASK(1:KSIZE)==1.)
    !
    ! compute RGMLTR
    !
    PRGMLTR(1:KSIZE)  = MAX( 0.0,( -PRGMLTR(1:KSIZE) *                     &
                           ( X0DEPG*       PLBDAG(1:KSIZE)**XEX0DEPG +     &
                             X1DEPG*PCJ(1:KSIZE)*PLBDAG(1:KSIZE)**XEX1DEPG ) -   &
                         ( PRG_TEND(1:KSIZE, IRCDRYG)+PRG_TEND(1:KSIZE, IRRDRYG) ) *       &
                               ( PRHODREF(1:KSIZE)*XCL*(XTT-PT(1:KSIZE))) ) /    &
                                             ( PRHODREF(1:KSIZE)*XLMTT ) )
  END WHERE
ENDIF

DO JL=1, KSIZE
  PA_RI(JL) = PA_RI(JL) - PRICFRRG(JL) - PRICFRR(JL)
  PA_RR(JL) = PA_RR(JL) - PRRCFRIG(JL) + PRICFRR(JL)
  PA_RG(JL) = PA_RG(JL) + PRICFRRG(JL) + PRRCFRIG(JL)
  PA_TH(JL) = PA_TH(JL) + (PRRCFRIG(JL) - PRICFRR(JL))*(PLSFACT(JL)-PLVFACT(JL))

  PA_RC(JL) = PA_RC(JL) - PRCWETG(JL)
  PA_RI(JL) = PA_RI(JL) - PRIWETG(JL)
  PA_RS(JL) = PA_RS(JL) - PRSWETG(JL)
  PA_RG(JL) = PA_RG(JL) + PRCWETG(JL) + PRIWETG(JL) + PRSWETG(JL) + PRRWETG(JL)
  PA_RR(JL) = PA_RR(JL) - PRRWETG(JL)
  PA_TH(JL) = PA_TH(JL) + (PRCWETG(JL) + PRRWETG(JL))*(PLSFACT(JL)-PLVFACT(JL))
  PA_RG(JL) = PA_RG(JL) - PRWETGH(JL)
  PB_RG(JL) = PB_RG(JL) - PRWETGH_MR(JL)
  IF (KRR==7) THEN
    PA_RH(JL) = PA_RH(JL) + PRWETGH(JL)
    PB_RH(JL) = PB_RH(JL) + PRWETGH_MR(JL)
  ENDIF
  PA_RC(JL) = PA_RC(JL) - PRCDRYG(JL)
  PA_RI(JL) = PA_RI(JL) - PRIDRYG(JL)
  PA_RS(JL) = PA_RS(JL) - PRSDRYG(JL)
  PA_RR(JL) = PA_RR(JL) - PRRDRYG(JL)
  PA_RG(JL) = PA_RG(JL) + PRCDRYG(JL) + PRIDRYG(JL) + PRSDRYG(JL) + PRRDRYG(JL)
  PA_TH(JL) = PA_TH(JL) + (PRCDRYG(JL)+PRRDRYG(JL))*(PLSFACT(JL)-PLVFACT(JL))

  PA_RR(JL) = PA_RR(JL) + PRGMLTR(JL)
  PA_RG(JL) = PA_RG(JL) - PRGMLTR(JL)
  PA_TH(JL) = PA_TH(JL) - PRGMLTR(JL)*(PLSFACT(JL)-PLVFACT(JL))
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RG', 1, ZHOOK_HANDLE)

END SUBROUTINE ICE4_FAST_RG
