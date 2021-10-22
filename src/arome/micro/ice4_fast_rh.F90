SUBROUTINE ICE4_FAST_RH(KSIZE, LDSOFT, PCOMPUTE, PWETG, &
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
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST
USE MODD_RAIN_ICE_PARAM
USE MODD_RAIN_ICE_DESCR
USE MODD_PARAM_ICE, ONLY : LEVLIMIT, LNULLWETH, LWETHPOST, LCONVHG
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN)    :: KSIZE
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
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRGWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRWETH  ! Dry growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRCDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRIDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRGDRYH  ! Wet growth of hailstone
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRDRYHG  ! Conversion of hailstone into graupel
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRHMLTR  ! Melting of the hailstones
REAL, DIMENSION(KSIZE, 10),   INTENT(INOUT) :: PRH_TEND ! Individual tendencies
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
LOGICAL, DIMENSION(KSIZE) :: GWET
REAL, DIMENSION(KSIZE) :: ZHAIL, ZWET, ZMASK, ZWETH, ZDRYH
INTEGER :: IHAIL, IGWET
REAL, DIMENSION(KSIZE) :: ZVEC1, ZVEC2, ZVEC3
INTEGER, DIMENSION(KSIZE) :: IVEC1, IVEC2
REAL, DIMENSION(KSIZE) :: ZZW, &
                                   ZRDRYH_INIT, ZRWETH_INIT, &
                                   ZRDRYHG
INTEGER :: JJ, JL
INTEGER :: IRCWETH, IRRWETH, IRIDRYH, IRIWETH, IRSDRYH, IRSWETH, IRGDRYH, IRGWETH, &
         & IFREEZ1, IFREEZ2
REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RH',0,ZHOOK_HANDLE)
!
IRCWETH=1
IRRWETH=2
IRIDRYH=3
IRIWETH=4
IRSDRYH=5
IRSWETH=6
IRGDRYH=7
IRGWETH=8
IFREEZ1=9
IFREEZ2=10
!
!
!
!*       7.2    compute the Wet and Dry growth of hail
!
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
           &MAX(0., -SIGN(1., XRTMIN(2)-PRCT(JL))) * & ! WHERE(PRCT(:)>XRTMIN(2))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRCWETH)=ZMASK(JL) * PRH_TEND(JL, IRCWETH)
  ENDDO
ELSE
  PRH_TEND(:, IRCWETH)=0.
  WHERE(ZMASK(:)==1.)
    ZZW(:) = PLBDAH(:)**(XCXH-XDH-2.0) * PRHODREF(:)**(-XCEXVT)
    PRH_TEND(:, IRCWETH)=XFWETH * PRCT(:) * ZZW(:)    ! RCWETH
  END WHERE
ENDIF
DO JL=1, KSIZE
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
           &MAX(0., -SIGN(1., XRTMIN(4)-PRIT(JL))) * & ! WHERE(PRIT(:)>XRTMIN(4))
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
  WHERE(ZMASK(:)==1.)
    ZZW(:) = PLBDAH(:)**(XCXH-XDH-2.0) * PRHODREF(:)**(-XCEXVT)
    PRH_TEND(:, IRIWETH)=XFWETH * PRIT(:) * ZZW(:)   ! RIWETH
    PRH_TEND(:, IRIDRYH)=PRH_TEND(:, IRIWETH)*(XCOLIH*EXP(XCOLEXIH*(PT(:)-XTT)))   ! RIDRYH
  END WHERE
ENDIF

!
!*       7.2.1  accretion of aggregates on the hailstones
!
DO JL=1, KSIZE
  ZWET(JL) = MAX(0., -SIGN(1., XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
            &MAX(0., -SIGN(1., XRTMIN(5)-PRST(JL))) * & ! WHERE(PRST(:)>XRTMIN(5))
            &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRSWETH)=ZWET(JL) * PRH_TEND(JL, IRSWETH)
    PRH_TEND(JL, IRSDRYH)=ZWET(JL) * PRH_TEND(JL, IRSDRYH)
  ENDDO
ELSE
  PRH_TEND(:, IRSWETH)=0.
  PRH_TEND(:, IRSDRYH)=0.
  GWET(:)=ZWET(:)==1.
  IGWET=COUNT(GWET(:))
  IF(IGWET>0)THEN
    !
    !*       7.2.3  select the (PLBDAH,PLBDAS) couplet
    !
    ZVEC1(1:IGWET) = PACK( PLBDAH(:),MASK=GWET(:) )
    ZVEC2(1:IGWET) = PACK( PLBDAS(:),MASK=GWET(:) )
    !
    !*       7.2.4  find the next lower indice for the PLBDAG and for the PLBDAS
    !               in the geometrical set of (Lbda_h,Lbda_s) couplet use to
    !               tabulate the SWETH-kernel
    !
    ZVEC1(1:IGWET) = MAX( 1.00001, MIN( FLOAT(NWETLBDAH)-0.00001,           &
                          XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
    IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
    ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )
    !
    ZVEC2(1:IGWET) = MAX( 1.00001, MIN( FLOAT(NWETLBDAS)-0.00001,           &
                          XWETINTP1S * LOG( ZVEC2(1:IGWET) ) + XWETINTP2S ) )
    IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
    ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
    !
    !*       7.2.5  perform the bilinear interpolation of the normalized
    !               SWETH-kernel
    !
    DO JJ = 1,IGWET
      ZVEC3(JJ) = (  XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_SWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                  - ( XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                  - XKER_SWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGWET),MASK=GWET,FIELD=0.0 )
    !
    WHERE(GWET(:))
      PRH_TEND(:, IRSWETH)=XFSWETH*ZZW(:)                       & ! RSWETH
                    *( PLBDAS(:)**(XCXS-XBS) )*( PLBDAH(:)**XCXH )  &
                       *( PRHODREF(:)**(-XCEXVT-1.) )               &
                       *( XLBSWETH1/( PLBDAH(:)**2              ) + &
                          XLBSWETH2/( PLBDAH(:)   * PLBDAS(:)   ) + &
                          XLBSWETH3/(               PLBDAS(:)**2) )
      PRH_TEND(:, IRSDRYH)=PRH_TEND(:, IRSWETH)*(XCOLSH*EXP(XCOLEXSH*(PT(:)-XTT)))
    END WHERE
  ENDIF
ENDIF
!
!*       7.2.6  accretion of graupeln on the hailstones
!
DO JL=1, KSIZE
  ZWET(JL)=MAX(0., -SIGN(1., XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
          &MAX(0., -SIGN(1., XRTMIN(6)-PRGT(JL))) * & ! WHERE(PRGT(:)>XRTMIN(6))
          &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRGWETH)=ZWET(JL) * PRH_TEND(JL, IRGWETH)
    PRH_TEND(JL, IRGDRYH)=ZWET(JL) * PRH_TEND(JL, IRGDRYH)
  ENDDO
ELSE
  PRH_TEND(:, IRGWETH)=0.
  PRH_TEND(:, IRGDRYH)=0.
  GWET(:)=ZWET(:)==1.
  IGWET=COUNT(GWET(:))
  IF(IGWET>0)THEN
    !
    !*       7.2.8  select the (PLBDAH,PLBDAG) couplet
    !
    ZVEC1(1:IGWET) = PACK( PLBDAH(:),MASK=GWET(:) )
    ZVEC2(1:IGWET) = PACK( PLBDAG(:),MASK=GWET(:) )
    !
    !*       7.2.9  find the next lower indice for the PLBDAH and for the PLBDAG
    !               in the geometrical set of (Lbda_h,Lbda_g) couplet use to
    !               tabulate the GWETH-kernel
    !
    ZVEC1(1:IGWET) = MAX( 1.00001, MIN( FLOAT(NWETLBDAG)-0.00001,           &
                          XWETINTP1H * LOG( ZVEC1(1:IGWET) ) + XWETINTP2H ) )
    IVEC1(1:IGWET) = INT( ZVEC1(1:IGWET) )
    ZVEC1(1:IGWET) = ZVEC1(1:IGWET) - FLOAT( IVEC1(1:IGWET) )
    !
    ZVEC2(1:IGWET) = MAX( 1.00001, MIN( FLOAT(NWETLBDAG)-0.00001,           &
                          XWETINTP1G * LOG( ZVEC2(1:IGWET) ) + XWETINTP2G ) )
    IVEC2(1:IGWET) = INT( ZVEC2(1:IGWET) )
    ZVEC2(1:IGWET) = ZVEC2(1:IGWET) - FLOAT( IVEC2(1:IGWET) )
    !
    !*       7.2.10 perform the bilinear interpolation of the normalized
    !               GWETH-kernel
    !
    DO JJ = 1,IGWET
      ZVEC3(JJ) = (  XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_GWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                 * ZVEC1(JJ) &
                - (  XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                   - XKER_GWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                        * (ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:) = UNPACK( VECTOR=ZVEC3(1:IGWET),MASK=GWET,FIELD=0.0 )
    !
    WHERE(GWET(:))
      PRH_TEND(:, IRGWETH)=XFGWETH*ZZW(:)                       & ! RGWETH
                    *( PLBDAG(:)**(XCXG-XBG) )*( PLBDAH(:)**XCXH )  &
                       *( PRHODREF(:)**(-XCEXVT-1.) )               &
                       *( XLBGWETH1/( PLBDAH(:)**2              ) + &
                          XLBGWETH2/( PLBDAH(:)   * PLBDAG(:)   ) + &
                          XLBGWETH3/(               PLBDAG(:)**2) )
      PRH_TEND(:, IRGDRYH)=PRH_TEND(:, IRGWETH)
    END WHERE
    !When graupel grows in wet mode, graupel is wet (!) and collection efficiency must remain the same
    WHERE(GWET(:) .AND. .NOT. PWETG(:)==1.)
      PRH_TEND(:, IRGDRYH)=PRH_TEND(:, IRGDRYH)*(XCOLGH*EXP(XCOLEXGH*(PT(:)-XTT)))
    END WHERE
  END IF
ENDIF
!
!*       7.2.11  accretion of raindrops on the hailstones
!
DO JL=1, KSIZE
  ZWET(JL)=MAX(0., -SIGN(1., XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
          &MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JL))) * & ! WHERE(PRRT(:)>XRTMIN(3))
          &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IRRWETH)=ZWET(JL) * PRH_TEND(JL, IRRWETH)
  ENDDO
ELSE
  PRH_TEND(:, IRRWETH)=0.
  GWET(:)=ZWET(:)==1.
  IGWET=COUNT(GWET(:))
  IF(IGWET>0)THEN
    !
    !*       7.2.12  select the (PLBDAH,PLBDAR) couplet
    !
    ZVEC1(1:IGWET)=PACK(PLBDAH(:), MASK=GWET(:))
    ZVEC2(1:IGWET)=PACK(PLBDAR(:), MASK=GWET(:))
    !
    !*       7.2.13 find the next lower indice for the PLBDAH and for the PLBDAR
    !               in the geometrical set of (Lbda_h,Lbda_r) couplet use to
    !               tabulate the RWETH-kernel
    !
    ZVEC1(1:IGWET)=MAX(1.00001, MIN( FLOAT(NWETLBDAH)-0.00001,           &
                          XWETINTP1H*LOG(ZVEC1(1:IGWET))+XWETINTP2H))
    IVEC1(1:IGWET)=INT(ZVEC1(1:IGWET))
    ZVEC1(1:IGWET)=ZVEC1(1:IGWET)-FLOAT(IVEC1(1:IGWET))
    !
    ZVEC2(1:IGWET)=MAX(1.00001, MIN( FLOAT(NWETLBDAR)-0.00001,           &
                          XWETINTP1R*LOG(ZVEC2(1:IGWET))+XWETINTP2R))
    IVEC2(1:IGWET)=INT(ZVEC2(1:IGWET))
    ZVEC2(1:IGWET)=ZVEC2(1:IGWET)-FLOAT(IVEC2(1:IGWET))
    !
    !*       7.2.14 perform the bilinear interpolation of the normalized
    !               RWETH-kernel
    !
    DO JJ=1, IGWET
      ZVEC3(JJ)= (  XKER_RWETH(IVEC1(JJ)+1,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RWETH(IVEC1(JJ)+1,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                                  * ZVEC1(JJ) &
                 - (  XKER_RWETH(IVEC1(JJ)  ,IVEC2(JJ)+1)* ZVEC2(JJ)          &
                    - XKER_RWETH(IVEC1(JJ)  ,IVEC2(JJ)  )*(ZVEC2(JJ) - 1.0) ) &
                                                         *(ZVEC1(JJ) - 1.0)
    END DO
    ZZW(:)=UNPACK(VECTOR=ZVEC3(1:IGWET), MASK=GWET, FIELD=0.)
    !
    WHERE(GWET(:))
      PRH_TEND(:, IRRWETH) = XFRWETH*ZZW(:)                    & ! RRWETH
                        *( PLBDAR(:)**(-4) )*( PLBDAH(:)**XCXH ) &
                               *( PRHODREF(:)**(-XCEXVT-1.) )   &
                    *( XLBRWETH1/( PLBDAH(:)**2              ) + &
                       XLBRWETH2/( PLBDAH(:)   * PLBDAR(:)   ) + &
                       XLBRWETH3/(               PLBDAR(:)**2) )
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
  ZHAIL(JL)=MAX(0., -SIGN(1., XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRH_TEND(JL, IFREEZ1)=ZHAIL(JL) * PRH_TEND(JL, IFREEZ1)
    PRH_TEND(JL, IFREEZ2)=ZHAIL(JL) * PRH_TEND(JL, IFREEZ2)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRH_TEND(JL, IFREEZ1)=PRVT(JL)*PPRES(JL)/(XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(LEVLIMIT) THEN
    WHERE(ZHAIL(:)==1.)
      PRH_TEND(:, IFREEZ1)=MIN(PRH_TEND(:, IFREEZ1), EXP(XALPI-XBETAI/PT(:)-XGAMI*ALOG(PT(:)))) ! min(ev, es_i(T))
    END WHERE
  ENDIF
  PRH_TEND(:, IFREEZ2)=0.
  WHERE(ZHAIL(:)==1.)
    PRH_TEND(:, IFREEZ1)=PKA(:)*(XTT-PT(:)) +                              &
             (PDV(:)*(XLVTT+(XCPV-XCL)*(PT(:)-XTT)) &
                           *(XESTT-PRH_TEND(:, IFREEZ1))/(XRV*PT(:))           )
    PRH_TEND(:, IFREEZ1)=PRH_TEND(:, IFREEZ1)* ( X0DEPH*       PLBDAH(:)**XEX0DEPH +     &
                           X1DEPH*PCJ(:)*PLBDAH(:)**XEX1DEPH )/ &
                          ( PRHODREF(:)*(XLMTT-XCL*(XTT-PT(:))) )
    PRH_TEND(:, IFREEZ2)=(PRHODREF(:)*(XLMTT+(XCI-XCL)*(XTT-PT(:)))   ) / &
                          ( PRHODREF(:)*(XLMTT-XCL*(XTT-PT(:))) )
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
IF(LNULLWETH) THEN
  DO JL=1, KSIZE
    ZWETH(JL) = ZWETH(JL) * MAX(0., -SIGN(1., -ZRDRYH_INIT(JL))) ! WHERE(ZRDRYH_INIT(:)>0.)
  ENDDO
ELSE
  DO JL=1, KSIZE
    ZWETH(JL) = ZWETH(JL) * MAX(0., -SIGN(1., -ZRWETH_INIT(JL))) ! WHERE(ZRWETH_INIT(:)>0.)
  ENDDO
ENDIF
IF(.NOT. LWETHPOST) THEN
  DO JL=1, KSIZE
    ZWETH(JL) = ZWETH(JL) * MAX(0., -SIGN(1., PT(JL)-XTT)) ! WHERE(PT(:)<XTT)
  ENDDO
ENDIF
DO JL=1, KSIZE
  ZDRYH(JL) = ZHAIL(JL) * &
            & MAX(0., -SIGN(1., PT(JL)-XTT)) * & ! WHERE(PT(:)<XTT)
            & MAX(0., -SIGN(1., -ZRDRYH_INIT(JL))) * & !WHERE(ZRDRYH_INIT(:)>0.)
            & MAX(0., -SIGN(1., MAX(0., ZRDRYH_INIT(JL)-PRH_TEND(JL, IRIDRYH)-PRH_TEND(JL, IRSDRYH)) - &
                               &MAX(0., ZRWETH_INIT(JL)-PRH_TEND(JL, IRIWETH)-PRH_TEND(JL, IRSWETH))))
ENDDO
!
ZRDRYHG(:)=0.
IF(LCONVHG)THEN
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
  ZMASK(JL)=MAX(0., -SIGN(1., XRTMIN(7)-PRHT(JL))) * & ! WHERE(PRHT(:)>XRTMIN(7))
           &MAX(0., -SIGN(1., XTT-PT(JL))) * & ! WHERE(PT(:)>XTT)
           &PCOMPUTE(JL)
ENDDO
IF(LDSOFT) THEN
  DO JL=1, KSIZE
    PRHMLTR(JL)=ZMASK(JL)*PRHMLTR(JL)
  ENDDO
ELSE
  DO JL=1, KSIZE
    PRHMLTR(JL) = ZMASK(JL)* PRVT(JL)*PPRES(JL)/(XEPSILO+PRVT(JL)) ! Vapor pressure
  ENDDO
  IF(LEVLIMIT) THEN
    WHERE(ZMASK(:)==1.)
      PRHMLTR(:)=MIN(PRHMLTR(:), EXP(XALPW-XBETAW/PT(:)-XGAMW*ALOG(PT(:)))) ! min(ev, es_w(T))
    END WHERE
  ENDIF
  DO JL=1, KSIZE
    PRHMLTR(JL) = ZMASK(JL)* (PKA(JL)*(XTT-PT(JL)) +                              &
           ( PDV(JL)*(XLVTT + ( XCPV - XCL ) * ( PT(JL) - XTT )) &
                           *(XESTT-PRHMLTR(JL))/(XRV*PT(JL))         ))
  ENDDO
  WHERE(ZMASK(:)==1.)
    !
    ! compute RHMLTR
    !
    PRHMLTR(:)  = MAX( 0.0,( -PRHMLTR(:) *                     &
                           ( X0DEPH*       PLBDAH(:)**XEX0DEPH +     &
                             X1DEPH*PCJ(:)*PLBDAH(:)**XEX1DEPH ) -   &
                         ( PRH_TEND(:, IRCWETH)+PRH_TEND(:, IRRWETH) )*        &
                               ( PRHODREF(:)*XCL*(XTT-PT(:))) ) /    &
                                             ( PRHODREF(:)*XLMTT ) )
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
