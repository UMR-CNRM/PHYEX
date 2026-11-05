MODULE MODE_ICE4_FAST_RI_RS
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_FAST_RI_RS(CST, PARAMI, ICEP, ICED, D, LDSOFT, LDCOMPUTE, &
                       &PRHODREF, PLSFACT, &
                       &PAI, PCIT, PESI, PPRES, &
                       &PSSI, PSSIO, PSSIU, PICLDFR, PIFR, &
                       &PRVT, PRIT, PRST, PT, &
                       &PRILARS, PRVDEPI)

!$ACDC singlecolumn --dummy

!!
!!**  PURPOSE
!!    -------
!!      Computes the deposition/evaporation of ice and conversion of large ice crystals into snow.
!!
!!    AUTHOR
!!    ------
!!    K. I Ivarsson  (Feb 2023)
!!
!!    MODIFICATIONS
!!    -------------
!!
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
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PLSFACT
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PAI      ! Thermodynamical function
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PCIT     ! Pristine ice conc. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRVT     ! Water vapor m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PSSI     ! Supersaturation over ice
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRIT     ! Pristine ice m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PESI     ! Saturation pressure over ice
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PPRES    ! Pressure
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PSSIO    ! Super-saturation with respect to ice in the 
                                                        ! supersaturated fraction
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PSSIU    ! Sub-saturation with respect to ice in the 
                                                        ! subsaturated fraction
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PICLDFR   ! Ice cloud fraction (fraction with supersat. with respect to ice)
REAL, DIMENSION(D%NIJT),      INTENT(IN)    :: PIFR    ! Ratio cloud ice moist part to dry part

REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRVDEPI  ! Conversion vapour to ice, non spherical effect
REAL, DIMENSION(D%NIJT),      INTENT(OUT)   :: PRILARS  ! Conversion of large ice crystals to snow

!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT) :: ZZW, &
                           ZZWC, &
                           ZCRYSHA, &
                           ZCI2S, &
                           ZW2D, &
                           ZXW2D13, &
                           ZWCITRED, &
                           ZWCITRED23
REAL :: ZDICRIT, ZTIMESC, ZKVO, ZTC, ZHU, ZQIMAX

INTEGER :: JIJ

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RI_RS',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       8.1   Turn large ice crystals to snow 
!
ZTIMESC = ICEP%XFRMIN(14) ! Time scale for few large ice crystal conversion to snow 

ZDICRIT = ICEP%XFRMIN(15)  ! Critical diameter of ice crystal to define what
                      ! is cloud ice and what is snow (m)


IF(PARAMI%LMODICEDEP .AND. ZTIMESC >0.)THEN

  ZDICRIT = (700.*CST%XPI/ICED%XAI/6.)**(1./ICED%XBI)*ZDICRIT**(3./ICED%XBI) ! from spherical diameter to maximum diameter
  ZKVO    = ((ICED%XALPHAI*ICED%XNUI + ICED%XBI -1.)/ICED%XALPHAI)**(1./ICED%XALPHAI)
  ZKVO    =  ZKVO/ZDICRIT/ZTIMESC
  !      PIFR = 1. ! Value of this parameter is very uncertain
  !      IF(XFRMIN(36)>0. )PIFR=  XFRMIN(36)
ENDIF

PRILARS(:) = 0.
ZWCITRED(:)= 0.1 ! ratio of ice crystal concentration wet to dry
                  ! part of a gridbox
DO JIJ=D%NIJB, D%NIJE
  IF( .NOT. LDSOFT) THEN
    ZZW(JIJ) = 0.
    IF(PARAMI%LMODICEDEP .AND. ZTIMESC > 0. )THEN
      ! Turn ice to snow if ice crystal distribution is such that
      ! the ice crystal diameter for the (mass x N_i) maximum 
      ! is larger than a prescribed size, 
      ! (ZDICRIT).
      ZZWC(JIJ) = (ICENUMBER3(PRIT(JIJ),PT(JIJ))+PCIT(JIJ))* &
            PRHODREF(JIJ)

      IF (PRIT(JIJ)>0.0 .AND. ZZWC(JIJ)>0.0) THEN   
        ZZW(JIJ) = MIN(1.E8,ICED%XLBI*( PRHODREF(JIJ)*PRIT(JIJ)/ZZWC(JIJ) )**ICED%XLBEXI) ! LAMBDA for ICE
        ZZW(JIJ) = MIN(0.5, 1. - 0.5**( ZKVO /ZZW(JIJ)))
        ZZW(JIJ) = MIN(0.9*PRIT(JIJ), ZZW(JIJ)*PRIT(JIJ))
        PRILARS(JIJ) = ZZW(JIJ)
      ENDIF
    ELSE
      ! Turn ice crystals lagrer than a precribed size into snow:
      ! (For the moment sperical ice crystals are assumed)
      IF (  (PRIT(JIJ)>0.0) .AND.(PSSI(JIJ)>0.001) ) THEN
        ZZW(JIJ) =   MIN(PRIT(JIJ),PRIT(JIJ)/(87.5*(ZDICRIT)**2*PAI(JIJ)/ PSSI(JIJ)))
        PRILARS(JIJ) =  ZZW(JIJ)
      ENDIF
    ENDIF
  ENDIF
ENDDO

PRVDEPI(:) = 0.
ZZW(:) = 0.
ZW2D(:) = 1./(PIFR(:)*PICLDFR(:) + 1. -PICLDFR(:))
ZZWC(:) = 0.
DO JIJ=D%NIJB, D%NIJE
  IF( .NOT. LDSOFT) THEN
    IF(PARAMI%LMODICEDEP)THEN
      ZZWC(JIJ) = (ICENUMBER3(PRIT(JIJ),PT(JIJ))+PCIT(JIJ))* &
                 PRHODREF(JIJ)
      IF( (PPRES(JIJ)*0.5 - PESI(JIJ))>0. .AND. ZZWC(JIJ)>0. .AND. LDCOMPUTE(JIJ)) THEN
        ZXW2D13(JIJ)= PIFR(JIJ)**(-ICED%XLBEXI)
        ZWCITRED23(JIJ)=ZWCITRED(JIJ)**(1.+ ICED%XLBEXI)     
        ZZW(JIJ)= ICEP%X0DEPI/(ICED%XLBI*PAI(JIJ)) *(ZZWC(JIJ)/PRHODREF(JIJ))**(1.+ICED%XLBEXI) * &
              & (MAX(ICED%XRTMIN(4),PRIT(JIJ))*ZW2D(JIJ) )**(-ICED%XLBEXI)
        ZZW(JIJ)= ZZW(JIJ)* ( PSSIO(JIJ)* PICLDFR(JIJ)* ZXW2D13(JIJ)  + ZWCITRED23(JIJ)*PSSIU(JIJ)* &
              & (1.-PICLDFR(JIJ)) )
      ENDIF
      PRVDEPI(JIJ) = PRVDEPI(JIJ) + ZZW(JIJ)
    ELSE
      IF ( (PPRES(JIJ)*0.5 - PESI(JIJ))>0. .AND. PCIT(JIJ)>0. .AND. LDCOMPUTE(JIJ) ) THEN
        ZXW2D13(JIJ)=PIFR(JIJ)**0.333
        ZWCITRED23(JIJ)=ZWCITRED(JIJ)**0.667
        ZTC =  MAX(-18.,MIN(-1.,PT(JIJ)-CST%XTT))
        ZHU =  MIN(0.15,MAX(0.,PSSI(JIJ)))
        ZCRYSHA(JIJ)=1.1+ 3.*ZHU*(1.+ SIN(0.64*ZTC -1.3))
        !       icedensity*4/3 *pi /8. =366.5 ; icedensity=700 kg/m3
        ZQIMAX = 366.5 * ZDICRIT**3 * PCIT(JIJ)*ZWCITRED(JIJ)/PRHODREF(JIJ)
        ZCI2S(JIJ) = 0.
        IF(PRIT(JIJ) > 1.0e-12 .AND. ZTIMESC > 0. )THEN
          ZCI2S(JIJ) = PRIT(JIJ)*(1. - MIN(1., 0.5*ZQIMAX /PRIT(JIJ)))* &
                   &  (1.-PICLDFR(JIJ))*ZW2D(JIJ)/(ZTIMESC*0.5)
        ENDIF

        ZZWC(JIJ)=ZCRYSHA(JIJ)*0.878/PAI(JIJ)*(PCIT(JIJ)/PRHODREF(JIJ))**0.667 &
              &*(MAX(ICED%XRTMIN(4),PRIT(JIJ))*ZW2D(JIJ))**0.333
        !     Ice supersaturated part of grid box:
        IF( PSSIO(JIJ)>0. .AND. PICLDFR(JIJ) > 0.02 ) THEN
          ZZW(JIJ) = ZZWC(JIJ)*ZXW2D13(JIJ)*PSSIO(JIJ)*PICLDFR(JIJ)
          PRVDEPI(JIJ) = PRVDEPI(JIJ) + ZZW(JIJ)
        ENDIF
        !    Ice subsaturated part of grid box:
        IF ( PSSIU(JIJ)<0. .AND. PICLDFR(JIJ) <0.98 ) THEN
          PRILARS(JIJ) = PRILARS(JIJ) + ZCI2S(JIJ)
          ZZW(JIJ) = ZZWC(JIJ)*ZWCITRED23(JIJ)*PSSIU(JIJ)*(1. - PICLDFR(JIJ))
          PRVDEPI(JIJ) = PRVDEPI(JIJ) + ZZW(JIJ)
        ENDIF
      ENDIF
    ENDIF
  ENDIF
ENDDO


IF (LHOOK) CALL DR_HOOK('ICE4_FAST_RI_RS', 1, ZHOOK_HANDLE)
!
CONTAINS
!
  FUNCTION ICENUMBER3 (Q_ICE, T3D) RESULT (RET)

    IMPLICIT NONE
    REAL, PARAMETER :: ICE_DENSITY = 890.0
    REAL, PARAMETER :: PI = 3.1415926536
    INTEGER :: IDX_REI
    REAL, INTENT(IN) :: Q_ICE, T3D
    REAL :: CORR, REICE, DEICE, RET
    DOUBLE PRECISION :: LAMBDA

    !+---+-----------------------------------------------------------------+
    !..Table of lookup values of radiative effective radius of ice crystals
    !.. as a function of Temperature from -94C to 0C.  Taken from WRF RRTMG
    !.. radiation code where it is attributed to Jon Egill Kristjansson
    !.. and coauthors.
    !+---+-----------------------------------------------------------------+

    REAL, DIMENSION(95), PARAMETER :: RETAB = (/                     &
      5.92779, 6.26422, 6.61973, 6.99539, 7.39234,                   &
      7.81177, 8.25496, 8.72323, 9.21800, 9.74075, 10.2930,          &
      10.8765, 11.4929, 12.1440, 12.8317, 13.5581, 14.2319,          &
      15.0351, 15.8799, 16.7674, 17.6986, 18.6744, 19.6955,          &
      20.7623, 21.8757, 23.0364, 24.2452, 25.5034, 26.8125,          &
      27.7895, 28.6450, 29.4167, 30.1088, 30.7306, 31.2943,          &
      31.8151, 32.3077, 32.7870, 33.2657, 33.7540, 34.2601,          &
      34.7892, 35.3442, 35.9255, 36.5316, 37.1602, 37.8078,          &
      38.4720, 39.1508, 39.8442, 40.5552, 41.2912, 42.0635,          &
      42.8876, 43.7863, 44.7853, 45.9170, 47.2165, 48.7221,          &
      50.4710, 52.4980, 54.8315, 57.4898, 60.4785, 63.7898,          &
      65.5604, 71.2885, 75.4113, 79.7368, 84.2351, 88.8833,          &
      93.6658, 98.5739, 103.603, 108.752, 114.025, 119.424,          &
      124.954, 130.630, 136.457, 142.446, 148.608, 154.956,          &
      161.503, 168.262, 175.248, 182.473, 189.952, 197.699,          &
      205.728, 214.055, 222.694, 231.661, 240.971, 250.639/)

    !+---+-----------------------------------------------------------------+
    !..From the model 3D temperature field, subtract 179K for which
    !.. index value of retab as a start.  Value of corr is for
    !.. interpolating between neighboring values in the table.
    !+---+-----------------------------------------------------------------+

    IDX_REI = INT(T3D-179.)
    IDX_REI = MIN(MAX(IDX_REI,1),95)
    CORR = T3D - INT(T3D)
    REICE = RETAB(IDX_REI)*(1.-CORR) + RETAB(MIN(95,IDX_REI+1))*CORR
    DEICE = 2.*REICE * 1.E-6

    !+---+-----------------------------------------------------------------+
    !..Now we have the final radiative effective size of ice (as function
    !.. of temperature only).  This size represents 3rd moment divided by
    !.. second moment of the ice size distribution, so we can compute a
    !.. number concentration from the mean size and mass mixing ratio.
    !.. The mean (radiative effective) diameter is 3./Slope for an inverse
    !.. exponential size distribution.  So, starting with slope, work
    !.. backwords to get number concentration.
    !+---+-----------------------------------------------------------------+

    LAMBDA = 3.0 / DEICE
    RET = Q_ICE * LAMBDA*LAMBDA*LAMBDA / (PI*ICE_DENSITY)

    !+---+-----------------------------------------------------------------+
    !..Example1:Common ice size coming from Thompson scheme is about 30 microns.
    !.. An example ice mixing ratio could be 0.001 g/kg for a temp. of -50C.
    !.. Remember to convert both into MKS units. This gives N_ice=357652 per kg.
    !..Example2: Lower in atmosphere at T=-10C matching ~162 microns in retab,
    !.. and assuming we have 0.1 g/kg mixing ratio, then N_ice=28122 per kg,
    !.. which is 28 crystals per liter of air if the air density is 1.0.
    !+---+-----------------------------------------------------------------+

  END FUNCTION ICENUMBER3
!
END SUBROUTINE ICE4_FAST_RI_RS
END MODULE MODE_ICE4_FAST_RI_RS
