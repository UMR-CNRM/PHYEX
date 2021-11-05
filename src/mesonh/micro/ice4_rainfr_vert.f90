!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_RAINFR_VERT
INTERFACE
SUBROUTINE ICE4_RAINFR_VERT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL, PPRFR, PRR, PRS, PRG, PRH)
IMPLICIT NONE
INTEGER,                      INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT) :: PPRFR !Precipitation fraction
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)    :: PRR !Rain field
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)    :: PRS !Snow field
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)    :: PRG !Graupel field
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL,INTENT(IN)    :: PRH !Hail field
END SUBROUTINE ICE4_RAINFR_VERT
END INTERFACE
END MODULE MODI_ICE4_RAINFR_VERT
SUBROUTINE ICE4_RAINFR_VERT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL, PPRFR, PRR, PRS, PRG, PRH)
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
USE MODD_RAIN_ICE_DESCR, ONLY : XRTMIN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL
REAL, DIMENSION(KIT,KJT,KKT), INTENT(INOUT)    :: PPRFR !Precipitation fraction
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)    :: PRR !Rain field
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)    :: PRS !Snow field
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)    :: PRG !Graupel field
REAL, DIMENSION(KIT,KJT,KKT), OPTIONAL, INTENT(IN)    :: PRH !Hail field
!
!*       0.2  declaration of local variables
!
INTEGER :: JI, JJ, JK
LOGICAL :: MASK
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
DO JI = KIB,KIE
   DO JJ = KJB, KJE
      PPRFR(JI,JJ,KKE)=0.
      DO JK=KKE-KKL, KKB, -KKL
         IF(PRESENT(PRH)) THEN
            MASK=PRR(JI,JJ,JK) .GT. XRTMIN(3) .OR. PRS(JI,JJ,JK) .GT. XRTMIN(5) &
            .OR. PRG(JI,JJ,JK) .GT. XRTMIN(6) .OR. PRH(JI,JJ,JK) .GT. XRTMIN(7)
         ELSE
            MASK=PRR(JI,JJ,JK) .GT. XRTMIN(3) .OR. PRS(JI,JJ,JK) .GT. XRTMIN(5) &
            .OR. PRG(JI,JJ,JK) .GT. XRTMIN(6) 
         END IF
         IF (MASK) THEN
            PPRFR(JI,JJ,JK)=MAX(PPRFR(JI,JJ,JK),PPRFR(JI,JJ,JK+KKL))
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
!
END SUBROUTINE ICE4_RAINFR_VERT
