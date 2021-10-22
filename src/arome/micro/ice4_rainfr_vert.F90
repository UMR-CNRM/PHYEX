SUBROUTINE ICE4_RAINFR_VERT(KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL, PPRFR, PRR)
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
!
!
!*      0. DECLARATIONS
!          ------------
!
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE MODD_RAIN_ICE_DESCR, ONLY : XRTMIN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                      INTENT(IN) :: KIB, KIE, KIT, KJB, KJE, KJT, KKB, KKE, KKT, KKL
REAL, DIMENSION(KIT,KJT,KKT), INTENT(OUT)    :: PPRFR !Precipitation fraction
REAL, DIMENSION(KIT,KJT,KKT), INTENT(IN)    :: PRR !Rain field
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JI, JJ, JK
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_RAINFR_VERT',0,ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
DO JI = KIB,KIE
   DO JJ = KJB, KJE
      PPRFR(JI,JJ,KKE)=0.
      DO JK=KKE-KKL, KKB, -KKL
         IF (PRR(JI,JJ,JK) .GT. XRTMIN(3)) THEN
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
IF (LHOOK) CALL DR_HOOK('ICE4_RAINFR_VERT',1,ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_RAINFR_VERT
