!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
MODULE MODE_ICE4_FAST_RH
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RH(CST, PARAMI, ICEP, ICED, KPROMA,KSIZE, LDSOFT, PCOMPUTE, PWETG, &
                       &PRHODREF, PLVFACT, PLSFACT, PPRES, &
                       &PDV, PKA, PCJ, &
                       &PLBDAS, PLBDAG, PLBDAR, PLBDAH, &
                       &PT,  PRVT, PRCT, PRRT, PRIT, PRST, PRGT, PRHT, &
                       &PRCWETH, PRIWETH, PRSWETH, PRGWETH, PRRWETH, &
                       &PRCDRYH, PRIDRYH, PRSDRYH, PRRDRYH, PRGDRYH, PRDRYHG, PRHMLTR, &
                       &PRH_TEND, &
                       &PA_TH, PA_RC, PA_RR, PA_RI, PA_RS, PA_RG, PA_RH)
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
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),                  INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),            INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),       INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)    :: ICED
INTEGER,                      INTENT(IN)    :: KPROMA,KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PWETG    ! 1. where graupel grows in wet mode, 0. elsewhere
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PPRES    ! absolute pressure at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PDV      ! Diffusivity of water vapor in the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PKA      ! Thermal conductivity of the air
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCJ      ! Function to compute the ventilation coefficient
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAG   ! Slope parameter of the graupel   distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAR   ! Slope parameter of the rain      distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAH   ! Slope parameter of the hail      distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRGT     ! Graupel m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHT     ! Hail m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRGWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRCDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRIDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRRDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRGDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRDRYHG  ! Conversion of hailstone into graupel
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRHMLTR  ! Melting of the hailstones
REAL, DIMENSION(KPROMA, 10),  INTENT(INOUT) :: PRH_TEND ! Individual tendencies
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RC
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RI
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RG
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PA_RH
!
!*       0.2  declaration of local variables
!
INTEGER, PARAMETER :: IRCWETH=1, IRRWETH=2, IRIDRYH=3, IRIWETH=4, IRSDRYH=5, IRSWETH=6, IRGDRYH=7, IRGWETH=8, &
                    & IFREEZ1=9, IFREEZ2=10
LOGICAL, DIMENSION(KSIZE) :: GWET
REAL, DIMENSION(KSIZE) :: ZHAIL, ZWET, ZMASK, ZWETH, ZDRYH
INTEGER :: IHAIL, IGWET
INTEGER, DIMENSION(KSIZE) :: I1
REAL, DIMENSION(KSIZE) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(KSIZE) :: IVEC1, IVEC2
REAL, DIMENSION(KSIZE) :: ZZW, &
                                   ZRDRYH_INIT, ZRWETH_INIT, &
                                   ZRDRYHG
INTEGER :: JJ, JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RH',0,ZHOOK_HANDLE)
!
!
!*       7.2    compute the Wet and Dry growth of hail
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., ICED%XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
           &MAX(0., -SIGN(1., ICED%XRTMIN(2)-PRCT(JL))) * & ! WHERE(PRCT(:)>XRTMIN(2))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRCWETH)=ZMASK(JL) * PRH_TEND(JL, IRCWETH)
  ENDDO
ELSE
  PRH_TEND(:, IRCWETH)=0.
  WHERE(ZMASK(1:KSIZE)==1.)
    ZZW(1:KSIZE) = PLBDAH(1:KSIZE)**(ICED%XCXH-ICED%XDH-2.0) * PRHODREF(1:KSIZE)**(-ICED%XCEXVT)
    PRH_TEND(1:KSIZE, IRCWETH)=ICEP%XFWETH * PRCT(1:KSIZE) * ZZW(1:KSIZE)    ! RCWETH
  END WHERE
ENDIF
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., ICED%XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
           &MAX(0., -SIGN(1., ICED%XRTMIN(4)-PRIT(JL))) * & ! WHERE(PRIT(:)>XRTMIN(4))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRIWETH)=ZMASK(JL) * PRH_TEND(JL, IRIWETH)
    PRH_TEND(JL, IRIDRYH)=ZMASK(JL) * PRH_TEND(JL, IRIDRYH)
  ENDDO
ELSE
  PRH_TEND(:, IRIWETH)=0.
  PRH_TEND(:, IRIDRYH)=0.
  WHERE(ZMASK(1:KSIZE)==1.)
    ZZW(1:KSIZE) = PLBDAH(1:KSIZE)**(ICED%XCXH-ICED%XDH-2.0) * PRHODREF(1:KSIZE)**(-ICED%XCEXVT)
    PRH_TEND(1:KSIZE, IRIWETH)=ICEP%XFWETH * PRIT(1:KSIZE) * ZZW(1:KSIZE)   ! RIWETH
    PRH_TEND(1:KSIZE, IRIDRYH)=PRH_TEND(1:KSIZE, IRIWETH)*(ICEP%XCOLIH*EXP(ICEP%XCOLEXIH*(PT(1:KSIZE)-CST%XTT)))   ! RIDRYH
  END WHERE
ENDIF

!
!*       7.2.1  accretion of aggregates on the hailstones
!
IGWET = 0
DO JJ = 1, KSIZE
  ZWET(JJ) = MAX(0., -SIGN(1., ICED%XRTMIN(7)-PRHT(JJ))) * & ! WHERE(PRHT(:)>XRTMIN(7))
            &MAX(0., -SIGN(1., ICED%XRTMIN(5)-PRST(JJ))) * & ! WHERE(PRST(:)>XRTMIN(5))
            &PCOMPUTE(JJ)
  IF (ZWET(JJ)>0) THEN
    IGWET = IGWET + 1
    I1(IGWET) = JJ
    GWET(JJ) = .TRUE.
  ELSE
    GWET(JJ) = .FALSE.
  ENDIF
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRSWETH)=ZWET(JL) * PRH_TEND(JL, IRSWETH)
    PRH_TEND(JL, IRSDRYH)=ZWET(JL) * PRH_TEND(JL, IRSDRYH)
  ENDDO
ELSE
  PRH_TEND(:, IRSWETH)=0.
  PRH_TEND(:, IRSDRYH)=0.
  IF(IGWET>0)THEN
    !
    !*       7.2.3  select the (PLBDAH,PLBDAS) couplet
    !
    DO JJ = 1, IGWET
      ZVEC1(JJ) = PLBDAH(I1(JJ))
      ZVEC2(JJ) = PLBDAS(I1(JJ))
    ENDDO
    !
    !*       7.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
    !               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
    !               tabulate the SWETH-kernel
    !
    ZVEC1(1:IGWET) = MAX( 1.00001, MIN( REAL(ICEP%NWETLBDAH)-0.00001,           &
                          ICEP%XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + ICEP%XWETINTP2H ) )
    IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
    ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - REAL( IVEC1(1:IGWET) )
    !
    ZVEC2(1:IGWET) = MAX( 1.00001, MIN( REAL(ICEP%NWETLBDAS)-0.00001,           &
                          ICEP%XWETINTP1S * LOG( ZVEC2(1:IGWET) ) + ICEP%XWETINTP2S ) )
    IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
    ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - REAL( IVEC2(1:IGWET) )
    !
    !*       7.2.5  perform the bilinear interpolation of the normalized
    !               SWETH-kernel
    !
    DO JJ = 1,IGWET
      ZVEC3(JJ) = (  ICEP%XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - ICEP%XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                  - ( ICEP%XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                  - ICEP%XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = 0.
    DO JJ = 1, IGWET
      ZZW(I1(JJ)) = ZVEC3(JJ)
    END DO
    !
    WHERE(GWET(1:KSIZE))
      PRH_TEND(1:KSIZE, IRSWETH)=ICEP%XFSWETH*ZZW(1:KSIZE)                       & ! RSWETH
                    *( PLBDAS(1:KSIZE)**(ICED%XCXS-ICED%XBS) )*( PLBDAH(1:KSIZE)**ICED%XCXH )  &
                       *( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) )               &
                       *( ICEP%XLBSWETH1/( PLBDAH(1:KSIZE)**2              ) + &
                          ICEP%XLBSWETH2/( PLBDAH(1:KSIZE)   * PLBDAS(1:KSIZE)   ) + &
                          ICEP%XLBSWETH3/(               PLBDAS(1:KSIZE)**2) )
      PRH_TEND(1:KSIZE, IRSDRYH)=PRH_TEND(1:KSIZE, IRSWETH)*(ICEP%XCOLSH*EXP(ICEP%XCOLEXSH*(PT(1:KSIZE)-CST%XTT)))
    END WHERE
  ENDIF
ENDIF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
IGWET = 0
DO JJ = 1, KSIZE
  ZWET(JJ)=MAX(0., -SIGN(1., ICED%XRTMIN(7)-PRHT(JJ))) * & ! WHERE(PRHT(:)>XRTMIN(7))
          &MAX(0., -SIGN(1., ICED%XRTMIN(6)-PRGT(JJ))) * & ! WHERE(PRGT(:)>XRTMIN(6))
          &PCOMPUTE(JJ)
  IF (ZWET(JJ)>0) THEN
    IGWET = IGWET + 1
    I1(IGWET) = JJ
    GWET(JJ) = .TRUE.
  ELSE
    GWET(JJ) = .FALSE.
  END IF
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRGWETH)=ZWET(JL) * PRH_TEND(JL, IRGWETH)
    PRH_TEND(JL, IRGDRYH)=ZWET(JL) * PRH_TEND(JL, IRGDRYH)
  ENDDO
ELSE
  PRH_TEND(:, IRGWETH)=0.
  PRH_TEND(:, IRGDRYH)=0.
  IF(IGWET>0)THEN
    !
    !*       7.2.8  select the (PLBDAH,PLBDAG) couplet
    !
    DO JJ = 1, IGWET
      ZVEC1(JJ) = PLBDAH(I1(JJ))
      ZVEC2(JJ) = PLBDAG(I1(JJ))
    END DO
    !
    !*       7.2.9  find the next lower indice for the PLBDAH and for the PLBDAG
    !               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
    !               tabulate the GWETH-kernel
    !
    ZVEC1(1:IGWET) = MAX( 1.00001, MIN( REAL(ICEP%NWETLBDAG)-0.00001,           &
                          ICEP%XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + ICEP%XWETINTP2H ) )
    IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
    ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - REAL( IVEC1(1:IGWET) )
    !
    ZVEC2(1:IGWET) = MAX( 1.00001, MIN( REAL(ICEP%NWETLBDAG)-0.00001,           &
                          ICEP%XWETINTP1G * LOG( ZVEC2(1:IGWET) ) + ICEP%XWETINTP2G ) )
    IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
    ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - REAL( IVEC2(1:IGWET) )
    !
    !*       7.2.10 perform the bilinear interpolation of the normalized
    !               GWETH-kernel
    !
    DO JJ = 1,IGWET
      ZVEC3(JJ) = (  ICEP%XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - ICEP%XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                - (  ICEP%XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - ICEP%XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                        * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = 0.
    DO JJ = 1, IGWET
      ZZW(I1(JJ)) = ZVEC3(JJ)
    END DO
    !
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
    WHERE(GWET(1:KSIZE) .AND. .NOT. PWETG(1:KSIZE)==1.)
      PRH_TEND(1:KSIZE, IRGDRYH)=PRH_TEND(1:KSIZE, IRGDRYH)*(ICEP%XCOLGH*EXP(ICEP%XCOLEXGH*(PT(1:KSIZE)-CST%XTT)))
    END WHERE
  END IF
ENDIF
!
!*       7.2.11  accretion of raindrops on the hailstones
!
IGWET = 0
DO JJ = 1, KSIZE
  ZWET(JJ)=MAX(0., -SIGN(1., ICED%XRTMIN(7)-PRHT(JJ))) * & ! WHERE(PRHT(:)>XRTMIN(7))
          &MAX(0., -SIGN(1., ICED%XRTMIN(3)-PRRT(JJ))) * & ! WHERE(PRRT(:)>XRTMIN(3))
          &PCOMPUTE(JJ)
  IF (ZWET(JJ)>0) THEN
    IGWET = IGWET + 1
    I1(IGWET) = JJ
    GWET(JJ) = .TRUE.
  ELSE
    GWET(JJ) = .FALSE.
  ENDIF
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRRWETH)=ZWET(JL) * PRH_TEND(JL, IRRWETH)
  ENDDO
ELSE
  PRH_TEND(:, IRRWETH)=0.
  IF(IGWET>0)THEN
    !
    !*       7.2.12  select the (PLBDAH,PLBDAR) couplet
    !
    DO JJ = 1, IGWET
      ZVEC1(JJ) = PLBDAH(I1(JJ))
      ZVEC2(JJ) = PLBDAR(I1(JJ))
    ENDDO
    !
    !*       7.2.13 find the next lower indice for the PLBDAH and for the PLBDAR
    !               in the geometrical set of (Lbda_h,Lbda_r) couplet use to
    !               tabulate the RWETH-kernel
    !
    ZVEC1(1:IGWET)=MAX(1.00001, MIN( REAL(ICEP%NWETLBDAH)-0.00001,           &
                          ICEP%XWETINTP1H*LOG(ZVEC1(1:IGWET))+ICEP%XWETINTP2H))
    IVEC1(1:IGWET)=INT(ZVEC1(1:IGWET))
    ZVEC1(1:IGWET)=ZVEC1(1:IGWET)-REAL(IVEC1(1:IGWET))
    !
    ZVEC2(1:IGWET)=MAX(1.00001, MIN( REAL(ICEP%NWETLBDAR)-0.00001,           &
                          ICEP%XWETINTP1R*LOG(ZVEC2(1:IGWET))+ICEP%XWETINTP2R))
    IVEC2(1:IGWET)=INT(ZVEC2(1:IGWET))
    ZVEC2(1:IGWET)=ZVEC2(1:IGWET)-REAL(IVEC2(1:IGWET))
    !
    !*       7.2.14 perform the bilinear interpolation of the normalized
    !               RWETH-kernel
    !
    DO JJ=1, IGWET
      ZVEC3(JJ)= (  ICEP%XKER_RWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_RWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                 - (  ICEP%XKER_RWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - ICEP%XKER_RWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         *(ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = 0.
    DO JJ = 1, IGWET
      ZZW(I1(JJ)) = ZVEC3(JJ)
    END DO
    !
    WHERE(GWET(1:KSIZE))
      PRH_TEND(1:KSIZE, IRRWETH) = ICEP%XFRWETH*ZZW(1:KSIZE)                    & ! RRWETH
                        *( PLBDAR(1:KSIZE)**(-4) )*( PLBDAH(1:KSIZE)**ICED%XCXH ) &
                               *( PRHODREF(1:KSIZE)**(-ICED%XCEXVT-1.) )   &
                    *( ICEP%XLBRWETH1/( PLBDAH(1:KSIZE)**2              ) + &
                       ICEP%XLBRWETH2/( PLBDAH(1:KSIZE)   * PLBDAR(1:KSIZE)   ) + &
                       ICEP%XLBRWETH3/(               PLBDAR(1:KSIZE)**2) )
    END WHERE
  ENDIF
ENDIF
!
DO JL=1, KSIZE
  ZRDRYH_INIT(JL)=PRH_TEND(JL, IRCWETH)+PRH_TEND(JL, IRIDRYH)+ &
                 &PRH_TEND(JL, IRSDRYH)+PRH_TEND(JL, IRRWETH)+PRH_TEND(JL, IRGDRYH)
ENDDO
!
!*       7.3    compute the Wet growth of hail
!
DO JL=1, KSIZE
  ZHAIL(JL)=MAX(0., -SIGN(1., ICED%XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IFREEZ1)=ZHAIL(JL) * PRH_TEND(JL, IFREEZ1)
    PRH_TEND(JL, IFREEZ2)=ZHAIL(JL) * PRH_TEND(JL, IFREEZ2)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRH_TEND(JL, IFREEZ1)=PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(PARAMI%LEVLIMIT) THEN
    WHERE(ZHAIL(1:KSIZE)==1.)
      PRH_TEND(1:KSIZE, IFREEZ1)=MIN(PRH_TEND(1:KSIZE, IFREEZ1), EXP(CST%XALPI-CST%XBETAI/PT(1:KSIZE)-CST%XGAMI*ALOG(PT(1:KSIZE)))) ! min(ev, es_i(T))
    END WHERE
  ENDIF
  PRH_TEND(:, IFREEZ2)=0.
  WHERE(ZHAIL(1:KSIZE)==1.)
    PRH_TEND(1:KSIZE, IFREEZ1)=PKA(1:KSIZE)*(CST%XTT-PT(1:KSIZE)) +                              &
             (PDV(1:KSIZE)*(CST%XLVTT+(CST%XCPV-CST%XCL)*(PT(1:KSIZE)-CST%XTT)) &
                           *(CST%XESTT-PRH_TEND(1:KSIZE, IFREEZ1))/(CST%XRV*PT(1:KSIZE))           )
    PRH_TEND(1:KSIZE, IFREEZ1)=PRH_TEND(1:KSIZE, IFREEZ1)* ( ICEP%X0DEPH*       PLBDAH(1:KSIZE)**ICEP%XEX0DEPH +     &
                           ICEP%X1DEPH*PCJ(1:KSIZE)*PLBDAH(1:KSIZE)**ICEP%XEX1DEPH )/ &
                          ( PRHODREF(1:KSIZE)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(1:KSIZE))) )
    PRH_TEND(1:KSIZE, IFREEZ2)=(PRHODREF(1:KSIZE)*(CST%XLMTT+(CST%XCI-CST%XCL)*(CST%XTT-PT(1:KSIZE)))   ) / &
                          ( PRHODREF(1:KSIZE)*(CST%XLMTT-CST%XCL*(CST%XTT-PT(1:KSIZE))) )
  END WHERE
ENDIF
DO JL=1, KSIZE
  !We must agregate, at least, the cold species
  ZRWETH_INIT(JL)=ZHAIL(JL) * MAX(PRH_TEND(JL, IRIWETH)+PRH_TEND(JL, IRSWETH)+PRH_TEND(JL, IRGWETH), &
                                 &MAX(0., PRH_TEND(JL, IFREEZ1) + &
                                         &PRH_TEND(JL, IFREEZ2) * ( &
                     &PRH_TEND(JL, IRIWETH)+PRH_TEND(JL, IRSWETH)+PRH_TEND(JL, IRGWETH) )))
ENDDO
!
!*       7.4    Select Wet or Dry case
!
!Wet case
DO JL=1, KSIZE
  ZWETH(JL) = ZHAIL(JL) * &
            & MAX(0., SIGN(1., MAX(0., ZRDRYH_INIT(JL)-PRH_TEND(JL, IRIDRYH)-PRH_TEND(JL, IRSDRYH)-PRH_TEND(JL, IRGDRYH)) - &
                              &MAX(0., ZRWETH_INIT(JL)-PRH_TEND(JL, IRIWETH)-PRH_TEND(JL, IRSWETH)-PRH_TEND(JL, IRGWETH))))
ENDDO
IF(PARAMI%LNULLWETH) THEN
  DO JL=1, KSIZE
    ZWETH(JL) = ZWETH(JL) * MAX(0., -SIGN(1., -ZRDRYH_INIT(JL))) ! WHERE(ZRDRYH_INIT(:)>0.)
  ENDDO
ELSE
  DO JL=1, KSIZE
    ZWETH(JL) = ZWETH(JL) * MAX(0., -SIGN(1., -ZRWETH_INIT(JL))) ! WHERE(ZRWETH_INIT(:)>0.)
  ENDDO
ENDIF
IF(.NOT. PARAMI%LWETHPOST) THEN
  DO JL=1, KSIZE
    ZWETH(JL) = ZWETH(JL) * MAX(0., -SIGN(1., PT(JL)-CST%XTT)) ! WHERE(PT(:)<XTT)
  ENDDO
ENDIF
DO JL=1, KSIZE
  ZDRYH(JL) = ZHAIL(JL) * &
            & MAX(0., -SIGN(1., PT(JL)-CST%XTT)) * & ! WHERE(PT(:)<XTT)
            & MAX(0., -SIGN(1., 1.E-20-ZRDRYH_INIT(JL))) * & !WHERE(ZRDRYH_INIT(:)>0.)
            & MAX(0., -SIGN(1., MAX(0., ZRDRYH_INIT(JL)-PRH_TEND(JL, IRIDRYH)-PRH_TEND(JL, IRSDRYH)) - &
                               &MAX(0., ZRWETH_INIT(JL)-PRH_TEND(JL, IRIWETH)-PRH_TEND(JL, IRSWETH))))
ENDDO
!
ZRDRYHG(:)=0.
IF(PARAMI%LCONVHG)THEN
  WHERE(ZDRYH(:)==1.)
    ZRDRYHG(:)=ZRDRYH_INIT(:)*ZRWETH_INIT(:)/(ZRDRYH_INIT(:)+ZRWETH_INIT(:))
  END WHERE
ENDIF
DO JL=1, KSIZE
  PRCWETH(JL) = ZWETH(JL) * PRH_TEND(JL, IRCWETH)
  PRIWETH(JL) = ZWETH(JL) * PRH_TEND(JL, IRIWETH)
  PRSWETH(JL) = ZWETH(JL) * PRH_TEND(JL, IRSWETH)
  PRGWETH(JL) = ZWETH(JL) * PRH_TEND(JL, IRGWETH)
  !Collected minus aggregated
  PRRWETH(JL) = ZWETH(JL) * (ZRWETH_INIT(JL) - PRH_TEND(JL, IRIWETH) - &
                             PRH_TEND(JL, IRSWETH) - PRH_TEND(JL, IRGWETH) - &
                             PRH_TEND(JL, IRCWETH))

  PRCDRYH(JL) = ZDRYH(JL) * PRH_TEND(JL, IRCWETH)
  PRIDRYH(JL) = ZDRYH(JL) * PRH_TEND(JL, IRIDRYH)
  PRSDRYH(JL) = ZDRYH(JL) * PRH_TEND(JL, IRSDRYH)
  PRRDRYH(JL) = ZDRYH(JL) * PRH_TEND(JL, IRRWETH)
  PRGDRYH(JL) = ZDRYH(JL) * PRH_TEND(JL, IRGDRYH)
  PRDRYHG(JL) = ZDRYH(JL) * ZRDRYHG(JL)

  PA_RC(JL) = PA_RC(JL) - PRCWETH(JL)
  PA_RI(JL) = PA_RI(JL) - PRIWETH(JL)
  PA_RS(JL) = PA_RS(JL) - PRSWETH(JL)
  PA_RG(JL) = PA_RG(JL) - PRGWETH(JL)
  PA_RH(JL) = PA_RH(JL) + PRCWETH(JL)+PRIWETH(JL)+PRSWETH(JL)+PRGWETH(JL)+PRRWETH(JL)
  PA_RR(JL) = PA_RR(JL) - PRRWETH(JL)
  PA_TH(JL) = PA_TH(JL) + (PRRWETH(JL)+PRCWETH(JL))*(PLSFACT(JL)-PLVFACT(JL))
  PA_RC(JL) = PA_RC(JL) - PRCDRYH(JL)
  PA_RI(JL) = PA_RI(JL) - PRIDRYH(JL)
  PA_RS(JL) = PA_RS(JL) - PRSDRYH(JL)
  PA_RR(JL) = PA_RR(JL) - PRRDRYH(JL)
  PA_RG(JL) = PA_RG(JL) - PRGDRYH(JL) + PRDRYHG(JL)
  PA_RH(JL) = PA_RH(JL) + PRCDRYH(JL)+PRIDRYH(JL)+PRSDRYH(JL)+&
                         &PRRDRYH(JL)+PRGDRYH(JL) - PRDRYHG(JL)
  PA_TH(JL) = PA_TH(JL) + (PRCDRYH(JL)+PRRDRYH(JL))*(PLSFACT(JL)-PLVFACT(JL))
ENDDO
!
!*       7.5    Melting of the hailstones
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., ICED%XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
           &MAX(0., -SIGN(1., CST%XTT-PT(JL))) * & ! WHERE(PT(:)>XTT)
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRHMLTR(JL)=ZMASK(JL)*PRHMLTR(JL)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRHMLTR(JL) = ZMASK(JL)* PRVT(JL)*PPRES(JL)/(CST%XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(PARAMI%LEVLIMIT) THEN
    WHERE(ZMASK(:)==1.)
      PRHMLTR(:)=MIN(PRHMLTR(:), EXP(CST%XALPW-CST%XBETAW/PT(:)-CST%XGAMW*ALOG(PT(:)))) ! min(ev, es_w(T))
    END WHERE
  ENDIF
  DO JL=1, KSIZE
    PRHMLTR(JL) = ZMASK(JL)* (PKA(JL)*(CST%XTT-PT(JL)) +                              &
           ( PDV(JL)*(CST%XLVTT + ( CST%XCPV - CST%XCL ) * ( PT(JL) - CST%XTT )) &
                           *(CST%XESTT-PRHMLTR(JL))/(CST%XRV*PT(JL))         ))
  ENDDO
  WHERE(ZMASK(1:KSIZE)==1.)
    !
    ! compute RHMLTR
    !
    PRHMLTR(1:KSIZE)  = MAX( 0.0,( -PRHMLTR(1:KSIZE) *                     &
                           ( ICEP%X0DEPH*       PLBDAH(1:KSIZE)**ICEP%XEX0DEPH +     &
                             ICEP%X1DEPH*PCJ(1:KSIZE)*PLBDAH(1:KSIZE)**ICEP%XEX1DEPH ) -   &
                         ( PRH_TEND(1:KSIZE, IRCWETH)+PRH_TEND(1:KSIZE, IRRWETH) )*        &
                               ( PRHODREF(1:KSIZE)*CST%XCL*(CST%XTT-PT(1:KSIZE))) ) /    &
                                             ( PRHODREF(1:KSIZE)*CST%XLMTT ) )
  END WHERE
END IF
DO JL=1, KSIZE
  PA_RR(JL) = PA_RR(JL) + PRHMLTR(JL)
  PA_RH(JL) = PA_RH(JL) - PRHMLTR(JL)
  PA_TH(JL) = PA_TH(JL) - PRHMLTR(JL)*(PLSFACT(JL)-PLVFACT(JL))
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RH', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_FAST_RH
END MODULE MODE_ICE4_FAST_RH
