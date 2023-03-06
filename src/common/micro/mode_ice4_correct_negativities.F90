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
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),                         INTENT(IN)    :: D
TYPE(RAIN_ICE_DESCR_t),                   INTENT(IN)    :: ICED
INTEGER,                                  INTENT(IN)    :: KRR
REAL, DIMENSION(D%NIJT, D%NKT),           INTENT(INOUT) :: PRV, PRC, PRR, PRI, PRS, PRG, PTH
REAL, DIMENSION(D%NIJT, D%NKT),           INTENT(IN)    :: PLVFACT, PLSFACT
REAL, DIMENSION(D%NIJT, D%NKT), OPTIONAL, INTENT(INOUT) :: PRH
!
REAL :: ZW
INTEGER :: JIJ, JK, IKTB, IKTE, IIJB, IIJE

REAL(KIND=JPRB) :: ZHOOK_HANDLE
!
IF (LHOOK) CALL DR_HOOK('ICE4_CORRECT_NEGATIVITIES', 0, ZHOOK_HANDLE)
!
IKTB=D%NKTB
IKTE=D%NKTE
IIJB=D%NIJB
IIJE=D%NIJE
!
!We correct negativities with conservation
DO JK = IKTB, IKTE
  DO JIJ = IIJB, IIJE
    ! 1) deal with negative values for mixing ratio, except for vapor
    ZW         =PRC(JIJ,JK)-MAX(PRC(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLVFACT(JIJ,JK)
    PRC(JIJ,JK)=PRC(JIJ,JK)-ZW

    ZW         =PRR(JIJ,JK)-MAX(PRR(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLVFACT(JIJ,JK)
    PRR(JIJ,JK)=PRR(JIJ,JK)-ZW

    ZW         =PRI(JIJ,JK)-MAX(PRI(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
    PRI(JIJ,JK)=PRI(JIJ,JK)-ZW

    ZW         =PRS(JIJ,JK)-MAX(PRS(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
    PRS(JIJ,JK)=PRS(JIJ,JK)-ZW

    ZW         =PRG(JIJ,JK)-MAX(PRG(JIJ,JK), 0.)
    PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
    PRG(JIJ,JK)=PRG(JIJ,JK)-ZW

    IF(KRR==7) THEN
      ZW         =PRH(JIJ,JK)-MAX(PRH(JIJ,JK), 0.)
      PRV(JIJ,JK)=PRV(JIJ,JK)+ZW
      PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
      PRH(JIJ,JK)=PRH(JIJ,JK)-ZW
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
!
IF (LHOOK) CALL DR_HOOK('ICE4_CORRECT_NEGATIVITIES', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_CORRECT_NEGATIVITIES
!
END MODULE MODE_ICE4_CORRECT_NEGATIVITIES
