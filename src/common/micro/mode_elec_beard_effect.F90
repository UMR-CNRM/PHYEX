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
  SUBROUTINE ELEC_BEARD_EFFECT(D, CST, HCLOUD, KID, OSEDIM, PT, PRHODREF, PTHVREFZIKB, &
                               PRX, PQX, PEFIELDW, PLBDA, PBEARDCOEF, ICED, &
                               LIMAP, LIMAPC, LIMAPW, LIMAPM)
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
USE MODD_PARAM_LIMA,      ONLY: PARAM_LIMA_t
USE MODD_PARAM_LIMA_COLD, ONLY: PARAM_LIMA_COLD_t
USE MODD_PARAM_LIMA_WARM, ONLY: PARAM_LIMA_WARM_t
USE MODD_PARAM_LIMA_MIXED,ONLY: PARAM_LIMA_MIXED_t
!
USE MODI_MOMG
USE MODE_MSG,             ONLY: PRINT_MSG, NVERB_FATAL
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments
!
TYPE(DIMPHYEX_t),                 INTENT(IN)    :: D
TYPE(CST_t),              INTENT(IN)    :: CST
INTEGER,                          INTENT(IN)    :: KID         ! Hydrometeor ID
LOGICAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: OSEDIM      ! if T, compute the sedim. proc.
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRHODREF    ! Reference density
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PT          ! Temperature
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PRX         ! m.r. source
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PQX         ! Elec. charge density source
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PEFIELDW    ! Vertical component of the electric field
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PLBDA       ! Slope param. of the distribution
REAL,    DIMENSION(D%NIJT,D%NKT), INTENT(OUT)   :: PBEARDCOEF  ! Beard coefficient
CHARACTER (LEN=4),      INTENT(IN)   ::  HCLOUD       ! Kind of microphysical scheme
!
REAL, INTENT(IN)                :: PTHVREFZIKB ! Reference thv at IKB for electricity
TYPE(RAIN_ICE_DESCR_t),   OPTIONAL, INTENT(IN)    :: ICED
TYPE(PARAM_LIMA_t),       OPTIONAL, INTENT(IN)    :: LIMAP
TYPE(PARAM_LIMA_COLD_t),  OPTIONAL, INTENT(IN)    :: LIMAPC
TYPE(PARAM_LIMA_WARM_t),  OPTIONAL, INTENT(IN)    :: LIMAPW
TYPE(PARAM_LIMA_MIXED_t), OPTIONAL, INTENT(IN)    :: LIMAPM
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
  IF(.NOT. PRESENT(ICED)) THEN
    CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ELEC_BEARD_EFFECT', 'ICED mandatory with ICE3/ICE4')
  ENDIF
  ZCEXVT = ICED%XCEXVT
  !
  IF (KID == 2) THEN
    ZBX     = ICED%XBC
    ZCX     = ICED%XCC
    ZDX     = ICED%XDC
    ZALPHAX = ICED%XALPHAC
    ZNUX    = ICED%XNUC
  ELSE IF (KID == 3) THEN
    ZBX     = ICED%XBR
    ZCX     = ICED%XCR
    ZDX     = ICED%XDR
    ZALPHAX = ICED%XALPHAR
    ZNUX    = ICED%XNUR
  ELSE IF (KID == 4) THEN
    ! values for columns are used to be consistent with the McF&H formula
    ZBX     = 1.7
    ZCX     = 2.1E5
    ZDX     = 1.585
    ZALPHAX = ICED%XALPHAI
    ZNUX    = ICED%XNUI
  ELSE IF (KID == 5) THEN
    ZBX     = ICED%XBS
    ZCX     = ICED%XCS
    ZDX     = ICED%XDS
    ZALPHAX = ICED%XALPHAS
    ZNUX    = ICED%XNUS
  ELSE IF (KID == 6) THEN
    ZBX     = ICED%XBG
    ZCX     = ICED%XCG
    ZDX     = ICED%XDG
    ZALPHAX = ICED%XALPHAG
    ZNUX    = ICED%XNUG
  ELSE IF (KID == 7) THEN
    ZBX     = ICED%XBH
    ZCX     = ICED%XCH
    ZDX     = ICED%XDH
    ZALPHAX = ICED%XALPHAH
    ZNUX    = ICED%XNUH
  END IF
ELSE IF (HCLOUD == 'LIMA') THEN
  IF(.NOT. PRESENT(LIMAP)) THEN
    CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ELEC_BEARD_EFFECT', 'LIMAP mandatory with ICE3/ICE4')
  ENDIF
  IF(.NOT. PRESENT(LIMAPC)) THEN
    CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ELEC_BEARD_EFFECT', 'LIMAPC mandatory with ICE3/ICE4')
  ENDIF
  IF(.NOT. PRESENT(LIMAPW)) THEN
    CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ELEC_BEARD_EFFECT', 'LIMAPW mandatory with ICE3/ICE4')
  ENDIF
  IF(.NOT. PRESENT(LIMAPM)) THEN
    CALL PRINT_MSG( NVERB_FATAL, 'GEN', 'ELEC_BEARD_EFFECT', 'LIMAPM mandatory with ICE3/ICE4')
  ENDIF

  ZCEXVT = LIMAP%XCEXVT
  !
  IF (KID == 2) THEN
    ZBX     = LIMAPW%XBC
    ZCX     = LIMAPW%XCC
    ZDX     = LIMAPW%XDC
    ZALPHAX = LIMAP%XALPHAC
    ZNUX    = LIMAP%XNUC
  ELSE IF (KID == 3) THEN
    ZBX     = LIMAPW%XBR
    ZCX     = LIMAPW%XCR
    ZDX     = LIMAPW%XDR
    ZALPHAX = LIMAP%XALPHAR
    ZNUX    = LIMAP%XNUR
  ELSE IF (KID == 4) THEN
    ZBX     = 1.7
    ZCX     = 2.1E5
    ZDX     = 1.585
    ZALPHAX = LIMAP%XALPHAI
    ZNUX    = LIMAP%XNUI
  ELSE IF (KID == 5) THEN
    ZBX     = LIMAPC%XBS
    ZCX     = LIMAPC%XCS
    ZDX     = LIMAPC%XDS
    ZALPHAX = LIMAP%XALPHAS
    ZNUX    = LIMAP%XNUS
  ELSE IF (KID == 6) THEN
    ZBX     = LIMAPM%XBG
    ZCX     = LIMAPM%XCG
    ZDX     = LIMAPM%XDG
    ZALPHAX = LIMAP%XALPHAG
    ZNUX    = LIMAP%XNUG
  ELSE IF (KID == 7) THEN
    ZBX     = LIMAPM%XBH
    ZCX     = LIMAPM%XCH
    ZDX     = LIMAPM%XDH
    ZALPHAX = LIMAPM%XALPHAH
    ZNUX    = LIMAPM%XNUH
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
!-------------------------------------------------------------------------------
!
END SUBROUTINE ELEC_BEARD_EFFECT
END MODULE MODE_ELEC_BEARD_EFFECT
