!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RAINFR_VERT
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_RAINFR_VERT(D, ICED, PPRFR, PRR, PRS, PRG, PRH)
!!
!!**  PURPOSE
!!    -------
!!      Computes the rain fraction
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the plitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!  P. Wautelet 13/02/2019: bugfix: intent of PPRFR OUT->INOUT
!
!
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_RAIN_ICE_DESCR, ONLY : RAIN_ICE_DESCR_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),             INTENT(IN)    :: D
TYPE(RAIN_ICE_DESCR_t),       INTENT(IN)    :: ICED
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(INOUT) :: PPRFR !Precipitation fraction
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRR !Rain field
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRS !Snow field
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PRG !Graupel field
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), OPTIONAL, INTENT(IN)    :: PRH !Hail field
!
INTEGER :: IKB, IKE, IKL, IIE, IIB, IJB, IJE
!*       0.2  declaration of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JI, JJ, JK
LOGICAL :: MASK
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_RAINFR_VERT',0,ZHOOK_HANDLE)
!
IKB=D%NKB
IKE=D%NKE
IKL=D%NKL
IIB=D%NIB
IIE=D%NIE
IJB=D%NJB
IJE=D%NJE
!
!-------------------------------------------------------------------------------
DO JI = IIB,IIE
   DO JJ = IJB, IJE
      PPRFR(JI,JJ,IKE)=0.
      DO JK=IKE-IKL, IKB, -IKL
         IF(PRESENT(PRH)) THEN
            MASK=PRR(JI,JJ,JK) .GT. ICED%XRTMIN(3) .OR. PRS(JI,JJ,JK) .GT. ICED%XRTMIN(5) &
            .OR. PRG(JI,JJ,JK) .GT. ICED%XRTMIN(6) .OR. PRH(JI,JJ,JK) .GT. ICED%XRTMIN(7)
         ELSE
            MASK=PRR(JI,JJ,JK) .GT. ICED%XRTMIN(3) .OR. PRS(JI,JJ,JK) .GT. ICED%XRTMIN(5) &
            .OR. PRG(JI,JJ,JK) .GT. ICED%XRTMIN(6) 
         END IF
         IF (MASK) THEN
            PPRFR(JI,JJ,JK)=MAX(PPRFR(JI,JJ,JK),PPRFR(JI,JJ,JK+IKL))
            IF (PPRFR(JI,JJ,JK)==0) THEN
               PPRFR(JI,JJ,JK)=1.
            END IF
         ELSE
            PPRFR(JI,JJ,JK)=0.
         END IF
      END DO
   END DO
END DO
!
IF (LHOOK) CALL DR_HOOK('ICE4_RAINFR_VERT',1,ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_RAINFR_VERT
END MODULE MODE_ICE4_RAINFR_VERT
