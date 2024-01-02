!MNH_LIC Copyright 2013-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!
MODULE MODE_ELEC_BEARD_EFFECT
!
IMPLICIT NONE
CONTAINS
!
! ###################################################################
  SUBROUTINE ELEC_BEARD_EFFECT(D, CST, ICED, HCLOUD, KID, OSEDIM, PT, PRHODREF, PTHVREFZIKB, &
                               PRX, PQX, PEFIELDW, PLBDA, PBEARDCOEF)
! ####################################################################
!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the effect of the electric field
!!      on the terminal velocity of hydrometeors.
!!
!!    METHOD
!!    ------
!!      From Beard, K. V., 1980: The Effects of Altitude and Electrical Force on 
!!           the Terminal Velocity of Hydrometeors. J. Atmos. Sci., 37, 1363–1374, 
!!           https://doi.org/10.1175/1520-0469(1980)037<1363:TEOAAE>2.0.CO;2. 
!!
!!    AUTHOR
!!    ------
!!      J.-P. Pinty      * LAERO *
!!      C. Barthe        * LAERO *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             21/08/2013   first coded in rain_ice_elec
!!      C. Barthe   01/06/2023 : externalize the code to use it with ICE3 and LIMA
!!      C. Barthe   08/06/2023 : correction by 10-5 of the dynamic viscosity of air
!!                               (unecessary for eta0/eta but necessary for Re0)
!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CST,             ONLY: CST_t
USE MODD_DIMPHYEX,        ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t

USE MODD_ELEC_DESCR,      ONLY: XRTMIN_ELEC
USE MODD_PARAM_LIMA,      ONLY: XALPHAC_L=>XALPHAC, XNUC_L=>XNUC, XALPHAR_L=>XALPHAR, XNUR_L=>XNUR, &
                                XALPHAI_L=>XALPHAI, XNUI_L=>XNUI, XALPHAS_L=>XALPHAS, XNUS_L=>XNUS, &
                                XALPHAG_L=>XALPHAG, XNUG_L=>XNUG,                                   &
                                XCEXVT_L=>XCEXVT
USE MODD_PARAM_LIMA_COLD, ONLY: XBI_L=>XBI, XC_I_L=>XC_I, XDI_L=>XDI, &
                                XBS_L=>XBS, XCS_L=>XCS, XDS_L=>XDS
USE MODD_PARAM_LIMA_MIXED,ONLY: XBG_L=>XBG, XCG_L=>XCG, XDG_L=>XDG, &
                                XBH_L=>XBH, XCH_L=>XCH, XDH_L=>XDH, &
                                XALPHAH_L=>XALPHAH, XNUH_L=>XNUH
USE MODD_PARAM_LIMA_WARM, ONLY: XBR_L=>XBR, XCR_L=>XCR, XDR_L=>XDR, &
                                XBC_L=>XBC, XCC_L=>XCC, XDC_L=>XDC
!
USE MODI_MOMG
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments
!
TYPE(DIMPHYEX_t),                 INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
INTEGER,                          INTENT(IN)    :: KID         ! Hydrometeor ID
LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: OSEDIM      ! if T, compute the sedim. proc.
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODREF    ! Reference density
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PT          ! Temperature
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRX         ! m.r. source
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PQX         ! Elec. charge density source
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PEFIELDW    ! Vertical component of the electric field
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PLBDA       ! Slope param. of the distribution
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(INOUT) :: PBEARDCOEF  ! Beard coefficient
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
!
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
!*       0.2   Declarations of local variables
!
INTEGER :: JIJ, JK ! loop indexes
INTEGER :: IIJB, IIJE, IKTB, IKTE
REAL :: ZCEXVT, ZBX, ZCX, ZDX, ZALPHAX, ZNUX
REAL :: ZRE0
REAL :: ZETA0
REAL :: ZVX
REAL :: ZK
REAL :: ZCOR00, ZRHO00
REAL :: ZT   ! Temperature (C)
REAL :: ZCOR ! To remove the Foote-duToit correction
REAL :: ZF0, ZF1  ! Coef. in Beard's equation
real, dimension(D%NIJT,D%NKT) :: zreynolds
!
!-------------------------------------------------------------------------------
!

ASSOCIATE(XALPHAC_I=>ICED%XALPHAC, XNUC_I=>ICED%XNUC, XALPHAR_I=>ICED%XALPHAR, XNUR_I=>ICED%XNUR, &
                                XALPHAI_I=>ICED%XALPHAI, XNUI_I=>ICED%XNUI, XALPHAS_I=>ICED%XALPHAS, XNUS_I=>ICED%XNUS, &
                                XALPHAG_I=>ICED%XALPHAG, XNUG_I=>ICED%XNUG, XALPHAH_I=>ICED%XALPHAH, XNUH_I=>ICED%XNUH, &
                                XBC_I=>ICED%XBC, XCC_I=>ICED%XCC, XDC_I=>ICED%XDC,                                 &
                                XBR_I=>ICED%XBR, XCR_I=>ICED%XCR, XDR_I=>ICED%XDR,                                 &
                                XBI_I=>ICED%XBI, XC_I_I=>ICED%XC_I, XDI_I=>ICED%XDI,                               &
                                XBS_I=>ICED%XBS, XCS_I=>ICED%XCS, XDS_I=>ICED%XDS,                                 &
                                XBG_I=>ICED%XBG, XCG_I=>ICED%XCG, XDG_I=>ICED%XDG,                                 &
                                XBH_I=>ICED%XBH, XCH_I=>ICED%XCH, XDH_I=>ICED%XDH,                                 &
                                XCEXVT_I=>ICED%XCEXVT)


!*       1.    COMPUTE USEFULL PARAMETERS
!              --------------------------
!
IKTB = D%NKTB
IKTE = D%NKTE
IIJB = D%NIJB
IIJE = D%NIJE
!
!*       1.1   Select the right parameters
!              --> depend on the microphysics scheme and the hydrometeor species
!
IF (HCLOUD(1:3) == 'ICE') THEN
  ZCEXVT = XCEXVT_I
  !
  IF (KID == 2) THEN
    ZBX     = XBC_I
    ZCX     = XCC_I
    ZDX     = XDC_I
    ZALPHAX = XALPHAC_I
    ZNUX    = XNUC_I
  ELSE IF (KID == 3) THEN
    ZBX     = XBR_I
    ZCX     = XCR_I
    ZDX     = XDR_I
    ZALPHAX = XALPHAR_I
    ZNUX    = XNUR_I
  ELSE IF (KID == 4) THEN
    ! values for columns are used to be consistent with the McF&H formula
    ZBX     = 1.7
    ZCX     = 2.1E5
    ZDX     = 1.585
    ZALPHAX = XALPHAI_I
    ZNUX    = XNUI_I
  ELSE IF (KID == 5) THEN
    ZBX     = XBS_I
    ZCX     = XCS_I
    ZDX     = XDS_I
    ZALPHAX = XALPHAS_I
    ZNUX    = XNUS_I
  ELSE IF (KID == 6) THEN
    ZBX     = XBG_I
    ZCX     = XCG_I
    ZDX     = XDG_I
    ZALPHAX = XALPHAG_I
    ZNUX    = XNUG_I
  ELSE IF (KID == 7) THEN
    ZBX     = XBH_I
    ZCX     = XCH_I
    ZDX     = XDH_I
    ZALPHAX = XALPHAH_I
    ZNUX    = XNUH_I
  END IF
ELSE IF (HCLOUD == 'LIMA') THEN
  ZCEXVT = XCEXVT_L
  !
  IF (KID == 2) THEN
    ZBX     = XBC_L
    ZCX     = XCC_L
    ZDX     = XDC_L
    ZALPHAX = XALPHAC_L
    ZNUX    = XNUC_L
  ELSE IF (KID == 3) THEN
    ZBX     = XBR_L
    ZCX     = XCR_L
    ZDX     = XDR_L
    ZALPHAX = XALPHAR_L
    ZNUX    = XNUR_L
  ELSE IF (KID == 4) THEN
    ZBX     = 1.7
    ZCX     = 2.1E5
    ZDX     = 1.585
    ZALPHAX = XALPHAI_L
    ZNUX    = XNUI_L
  ELSE IF (KID == 5) THEN
    ZBX     = XBS_L
    ZCX     = XCS_L
    ZDX     = XDS_L
    ZALPHAX = XALPHAS_L
    ZNUX    = XNUS_L
  ELSE IF (KID == 6) THEN
    ZBX     = XBG_L
    ZCX     = XCG_L
    ZDX     = XDG_L
    ZALPHAX = XALPHAG_L
    ZNUX    = XNUG_L
  ELSE IF (KID == 7) THEN
    ZBX     = XBH_L
    ZCX     = XCH_L
    ZDX     = XDH_L
    ZALPHAX = XALPHAH_L
    ZNUX    = XNUH_L
  END IF
  !
END IF  
!
!*       1.2   Parameters from Table 1 in Beard (1980)
!
! Reference value of the dynamic viscosity of air
ZETA0 = (1.718E-5 + 0.0049E-5 * (PTHVREFZIKB - CST%XTT))
!
ZRHO00 = CST%XP00 / (CST%XRD * PTHVREFZIKB)
ZCOR00 = ZRHO00**ZCEXVT
!
! (rho_0 / eta_0) * (v * lambda^d)
ZVX = (ZRHO00 / ZETA0) * ZCX * MOMG(ZALPHAX,ZNUX,ZBX+ZDX) / MOMG(ZALPHAX,ZNUX,ZBX)
!
!-------------------------------------------------------------------------------
!
!*       2.    COMPUTE THE VELOCITY ADJUSTMENT FACTOR
!              --------------------------------------
!
zreynolds(:,:) = 0.
PBEARDCOEF(:,:) = 1.0
!
DO JK = IKTB, IKTE
  DO JIJ = IIJB, IIJE
!++cb++ 09/06/23 on n'applique l'effet Beard que pour les points ou le rapport de melange est
! suffisamment eleve pour eviter que qE >> mg => coef de Beard tres eleve !
! Ce pb intervient avec ICE3 pour lequel xrtmin est tres bas par rapport a LIMA.
    IF (OSEDIM(JIJ,JK) .AND. PRX(JIJ,JK) .GT. XRTMIN_ELEC(KID) .AND. PLBDA(JIJ,JK) .GT. 0.) THEN
!--cb--
      ! Temperature K --> C
      ZT = PT(JIJ,JK) - CST%XTT
      !
      ! Pre-factor of f_0
      IF (ZT >= 0.0) THEN
        ZF0 = ZETA0 / (1.718E-5 + 0.0049E-5 * ZT)
      ELSE
        ZF0 = ZETA0 / (1.718E-5 + 0.0049E-5 * ZT - 1.2E-10 * ZT * ZT)
      END IF
      !
      ! Pre-factor of f_infty
      ZF1 = SQRT(ZRHO00/PRHODREF(JIJ,JK))
      !
      ! compute (1 - K) = 1 - qE/mg
      ZK = 1. - PQX(JIJ,JK) * PEFIELDW(JIJ,JK) / (PRX(JIJ,JK) * CST%XG)
      !
      ! Hyp : K_0 ~ 0
      ! Hyp : si qE > mg, K > 1
      IF (ZK <= 0.0) THEN
        PBEARDCOEF(JIJ,JK) = 0.  ! levitation
      ELSE
        ! Reynolds number
        ZRE0 = ZVX / PLBDA(JIJ,JK)**(1.+ZDX)
        zreynolds(jij,jk) = zre0
        IF (ZRE0 <= 0.2) THEN
          PBEARDCOEF(JIJ,JK) = ZF0 * ZK
        ELSE IF (ZRE0 >= 1000.) THEN
          PBEARDCOEF(JIJ,JK) = ZF1 * SQRT(ZK)
        ELSE
          PBEARDCOEF(JIJ,JK) = ZF0 * ZK +                   &
                              (ZF1 * SQRT(ZK) - ZF0 * ZK) * &
                              (1.61 + LOG(ZRE0)) / 8.52
        END IF
        ! remove the Foote-duToit correction
        ZCOR = (PRHODREF(JIJ,JK) / ZRHO00)**ZCEXVT
        PBEARDCOEF(JIJ,JK) = PBEARDCOEF(JIJ,JK) * ZCOR
      END IF
    ELSE
      PBEARDCOEF(JIJ,JK) = 1.0 ! No "Beard" effect
    END IF
  END DO
END DO
!
END ASSOCIATE
!-------------------------------------------------------------------------------
!
END SUBROUTINE ELEC_BEARD_EFFECT
END MODULE MODE_ELEC_BEARD_EFFECT
