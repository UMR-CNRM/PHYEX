!MNH_LIC Copyright 1995-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_CORRECT_NEGATIVITIES
IMPLICIT NONE
CONTAINS
       SUBROUTINE ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, PR, &
                                    &PTH, PLVFACT, PLSFACT)
        !SUBROUTINE ICE4_CORRECT_NEGATIVITIES(D, ICED, KRR, OELEC, PR, &
!                                    &PTH, PLVFACT, PLSFACT, &
!                                    &PQPI, PQC, PQR, PQI, PQS, PQG, PQNI, &
!                                    &PTH, PLVFACT, PLSFACT, PQH)
!
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
!USE MODD_ELEC_DESCR, ONLY: XECHARGE
USE MODD_FIELDS_ADDRESS
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),                         INTENT(IN)    :: D
TYPE(RAIN_ICE_DESCR_t),                   INTENT(IN)    :: ICED
INTEGER,                                  INTENT(IN)    :: KRR
!LOGICAL,                                  INTENT(IN)    :: OELEC
REAL, DIMENSION(D%NIJT, D%NKT, KRR),      INTENT(INOUT) :: PR
REAL, DIMENSION(D%NIJT, D%NKT),           INTENT(INOUT) :: PTH
!REAL, DIMENSION(D%NIJT, D%NKT),           INTENT(INOUT) :: PQPI, PQC, PQR, PQI, PQS, PQG, PQNI
REAL, DIMENSION(D%NIJT, D%NKT),           INTENT(IN)    :: PLVFACT, PLSFACT
!REAL, DIMENSION(D%NIJT, D%NKT), OPTIONAL, INTENT(INOUT) :: PQH
!
REAL :: ZW
!REAL :: ZION, ZADD
INTEGER :: JIJ, JK, IKTB, IKTE, IIJB, IIJE, JRR

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
    DO JRR=2, 3
      ZW         =PR(JIJ,JK,JRR)-MAX(PR(JIJ,JK,JRR), 0.)
      PR(JIJ,JK,IRV)=PR(JIJ,JK,IRV)+ZW
      PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLVFACT(JIJ,JK)
      PR(JIJ,JK,JRR)=PR(JIJ,JK,JRR)-ZW
      !++cb-- pour l'elec, on peut eventuellement appeler une routine : ca evitera les pb avec xecharge ?
      !IF (OELEC .AND. ZW .LT. 0.) THEN
      !  ZION = PQC(JIJ,JK) / XECHARGE
      !  ZADD = 0.5 + SIGN(0.5, PQC(JIJ,JK))
      !  PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
      !  PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
      !  PQC(JIJ,JK)  = 0.
      !END IF
    ENDDO
    DO JRR=4, KRR
      ZW         =PR(JIJ,JK,JRR)-MAX(PR(JIJ,JK,JRR), 0.)
      PR(JIJ,JK,IRV)=PR(JIJ,JK,IRV)+ZW
      PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
      PR(JIJ,JK,JRR)=PR(JIJ,JK,JRR)-ZW
      !IF (OELEC .AND. ZW .LT. 0.) THEN
      !  ZION = PQI(JIJ,JK) / XECHARGE
      !  ZADD = 0.5 + SIGN(0.5, PQI(JIJ,JK))
      !  PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
      !  PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
      !  PQI(JIJ,JK)  = 0.
      !END IF    
    ENDDO

    ! 2) deal with negative vapor mixing ratio

    ! for rc and ri, we keep ice fraction constant
    ZW=MIN(1., MAX(ICED%XRTMIN(IRV)-PR(JIJ,JK,IRV), 0.) / &
              &MAX(PR(JIJ,JK,IRC)+PR(JIJ,JK,IRI), 1.E-20)) ! Proportion of rc+ri to convert into rv
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW* &
                 &(PR(JIJ,JK,IRC)*PLVFACT(JIJ,JK)+PR(JIJ,JK,IRI)*PLSFACT(JIJ,JK))
    PR(JIJ,JK,IRV)=PR(JIJ,JK,IRV)+ZW*(PR(JIJ,JK,IRC)+PR(JIJ,JK,IRI))
    PR(JIJ,JK,IRC)=(1.-ZW)*PR(JIJ,JK,IRC)
    PR(JIJ,JK,IRI)=(1.-ZW)*PR(JIJ,JK,IRI)

    JRR=IRR
    ZW=MIN(MAX(PR(JIJ,JK,JRR), 0.), &
          &MAX(ICED%XRTMIN(IRV)-PR(JIJ,JK,IRV), 0.)) ! Quantity of rr to convert into rv
    PR(JIJ,JK,IRV)=PR(JIJ,JK,IRV)+ZW
    PR(JIJ,JK,JRR)=PR(JIJ,JK,JRR)-ZW
    PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLVFACT(JIJ,JK)
!    IF (OELEC .AND. ZW .GT. 0.) THEN
!      ZION = PQR(JIJ,JK) / XECHARGE
!      ZADD = 0.5 + SIGN(0.5, PQG(JIJ,JK))
!      PQPI(JIJ,JK) = PQPI(JIJ,JK) +       ZADD  * ZION
!      PQNI(JIJ,JK) = PQNI(JIJ,JK) + (1. - ZADD) * ZION
!      PQR(JIJ,JK)  = 0.
!    END IF    

    DO JRR=5, KRR
      ZW=MIN(MAX(PR(JIJ,JK,JRR), 0.), &
            &MAX(ICED%XRTMIN(IRV)-PR(JIJ,JK,IRV), 0.)) ! Quantity of rs to convert into rv
      PR(JIJ,JK,IRV)=PR(JIJ,JK,IRV)+ZW
      PR(JIJ,JK,JRR)=PR(JIJ,JK,JRR)-ZW
      PTH(JIJ,JK)=PTH(JIJ,JK)-ZW*PLSFACT(JIJ,JK)
    ENDDO
  ENDDO
ENDDO
!$acc end kernels
!
IF (LHOOK) CALL DR_HOOK('ICE4_CORRECT_NEGATIVITIES', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_CORRECT_NEGATIVITIES
!
END MODULE MODE_ICE4_CORRECT_NEGATIVITIES
