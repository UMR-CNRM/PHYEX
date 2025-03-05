!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_CORRECT_NEGATIVITIES
IMPLICIT NONE
CONTAINS
       SUBROUTINE ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, PRV, PRC, PRR, &
                                    &PRI, PRS, PRG, &
                                    &PTH, PLVFACT, PLSFACT, PRH)
        !SUBROUTINE ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, OELEC, PRV, PRC, PRR, &
!                                    &PRI, PRS, PRG, &
!                                    &PTH, PLVFACT, PLSFACT, PRH, &
!                                    &PQPI, PQC, PQR, PQI, PQS, PQG, PQNI, &
!                                    &PTH, PLVFACT, PLSFACT, PRH, PQH)
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
!USE MODD_ELEC_DESCR, ONLY: XECHARGE
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),                         INTENT(IN)    :: D
TYPE(RAIN_ICE_DESCR_t),                   INTENT(IN)    :: ICED
INTEGER,                                  INTENT(IN)    :: KRR
!LOGICAL,                                  INTENT(IN)    :: OELEC
REAL, DIMENSION(D%NIJT, D%NKT),           INTENT(INOUT) :: PRV, PRC, PRR, PRI, PRS, PRG, PTH
!REAL, DIMENSION(D%NIJT, D%NKT),           INTENT(INOUT) :: PQPI, PQC, PQR, PQI, PQS, PQG, PQNI
REAL, DIMENSION(D%NIJT, D%NKT),           INTENT(IN)    :: PLVFACT, PLSFACT
REAL, DIMENSION(D%NIJT, D%NKT), OPTIONAL, INTENT(INOUT) :: PRH
!REAL, DIMENSION(D%NIJT, D%NKT), OPTIONAL, INTENT(INOUT) :: PQH
!
REAL :: ZW
!REAL :: ZION, ZADD
INTEGER :: JIJ, JK, IKTB, IKTE, IIJB, IIJE

REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ICE4_CORRECT_NEGATIVITIES', 0, ZHOOK_HANDLE)
!
IKTB=D%NKTB
IKTE=D%NKTE
IIJB=D%NIJB
IIJE=D%NIJE
!
!We correct negativities with conservation
!$acc kernels
!$acc loop independent
DO JK = IKTB, IKTE
  DO JIJ = IIJB, IIJE
    ! 1) deal with negative values for mixing ratio, except for vapor
    ZW         =PRC(JIJ,JK)-MAX(PRC(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLVFACT(JIJ,JK)
    PRC(JIJ,JK)=PRC(JIJ,JK)-ZW
!++cb-- pour l'elec, on peut eventuellement appeler une routine : ca evitera les pb avec xecharge ?
!    IF (OELEC .AND. ZW .LT. 0.) THEN
!      ZION = PQC(JIJ,JK) / XECHARGE
!      ZADD = 0.5 + SIGN(0.5, PQC(JIJ,JK))
!      PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
!      PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
!      PQC(JIJ,JK)  = 0.
!    END IF

    ZW         =PRR(JIJ,JK)-MAX(PRR(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLVFACT(JIJ,JK)
    PRR(JIJ,JK)=PRR(JIJ,JK)-ZW
!    IF (OELEC .AND. ZW .LT. 0.) THEN
!      ZION = PQR(JIJ,JK) / XECHARGE
!      ZADD = 0.5 + SIGN(0.5, PQR(JIJ,JK))
!      PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
!      PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
!      PQR(JIJ,JK)  = 0.
!    END IF

    ZW         =PRI(JIJ,JK)-MAX(PRI(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
    PRI(JIJ,JK)=PRI(JIJ,JK)-ZW
!    IF (OELEC .AND. ZW .LT. 0.) THEN
!      ZION = PQI(JIJ,JK) / XECHARGE
!      ZADD = 0.5 + SIGN(0.5, PQI(JIJ,JK))
!      PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
!      PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
!      PQI(JIJ,JK)  = 0.
!    END IF    

    ZW         =PRS(JIJ,JK)-MAX(PRS(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
    PRS(JIJ,JK)=PRS(JIJ,JK)-ZW
!    IF (OELEC .AND. ZW .LT. 0.) THEN
!      ZION = PQS(JIJ,JK) / XECHARGE
!      ZADD = 0.5 + SIGN(0.5, PQS(JIJ,JK))
!      PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
!      PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
!      PQS(JIJ,JK)  = 0.
!    END IF    

    ZW         =PRG(JIJ,JK)-MAX(PRG(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
    PRG(JIJ,JK)=PRG(JIJ,JK)-ZW
!    IF (OELEC .AND. ZW .LT. 0.) THEN
!      ZION = PQG(JIJ,JK) / XECHARGE
!      ZADD = 0.5 + SIGN(0.5, PQG(JIJ,JK))
!      PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
!      PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
!      PQG(JIJ,JK)  = 0.
!    END IF    

    IF(KRR==7) THEN
      ZW         =PRH(JIJ,JK)-MAX(PRH(JIJ,JK), 0.)
      PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
      PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
      PRH(JIJ,JK)=PRH(JIJ,JK)-ZW
!      IF (OELEC .AND. ZW .LT. 0.) THEN
!        ZION = PQH(JIJ,JK) / XECHARGE
!        ZADD = 0.5 + SIGN(0.5, PQH(JIJ,JK))
!        PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
!        PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
!        PQH(JIJ,JK)  = 0.
!      END IF      
    ENDIF

    ! 2) deal with negative vapor mixing ratio

    ! for rc and ri, we keep ice fraction constant
    ZW=MIN(1., MAX(ICED%XRTMIN(1)-PRV(JIJ,JK), 0.) / &
              &MAX(PRC(JIJ,JK)+PRI(JIJ,JK), 1.E-20)) ! Proportion of rc+ri to convert into rv
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW* &
                 &(PRC(JIJ,JK)*PLVFACT(JIJ,JK)+PRI(JIJ,JK)*PLSFACT(JIJ,JK))
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW*(PRC(JIJ,JK)+PRI(JIJ,JK))
    PRC(JIJ,JK)=(1.-ZW)*PRC(JIJ,JK)
    PRI(JIJ,JK)=(1.-ZW)*PRI(JIJ,JK)

    ZW=MIN(MAX(PRR(JIJ,JK), 0.), &
          &MAX(ICED%XRTMIN(1)-PRV(JIJ,JK), 0.)) ! Quantity of rr to convert into rv
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PRR(JIJ,JK)=PRR(JIJ,JK)-ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLVFACT(JIJ,JK)
!    IF (OELEC .AND. ZW .GT. 0.) THEN
!      ZION = PQR(JIJ,JK) / XECHARGE
!      ZADD = 0.5 + SIGN(0.5, PQG(JIJ,JK))
!      PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
!      PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
!      PQR(JIJ,JK)  = 0.
!    END IF    

    ZW=MIN(MAX(PRS(JIJ,JK), 0.), &
          &MAX(ICED%XRTMIN(1)-PRV(JIJ,JK), 0.)) ! Quantity of rs to convert into rv
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PRS(JIJ,JK)=PRS(JIJ,JK)-ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)

    ZW=MIN(MAX(PRG(JIJ,JK), 0.), &
          &MAX(ICED%XRTMIN(1)-PRV(JIJ,JK), 0.)) ! Quantity of rg to convert into rv
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PRG(JIJ,JK)=PRG(JIJ,JK)-ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)

    IF(KRR==7) THEN
      ZW=MIN(MAX(PRH(JIJ,JK), 0.), &
            &MAX(ICED%XRTMIN(1)-PRV(JIJ,JK), 0.)) ! Quantity of rh to convert into rv
      PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
      PRH(JIJ,JK)=PRH(JIJ,JK)-ZW
      PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
    ENDIF
  ENDDO
ENDDO
!$acc end kernels
!
IF (LHOOK) CALL DR_HOOK('ICE4_CORRECT_NEGATIVITIES', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_CORRECT_NEGATIVITIES
!
END MODULE MODE_ICE4_CORRECT_NEGATIVITIES
