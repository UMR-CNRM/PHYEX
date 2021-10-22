SUBROUTINE ICE4_RRHONG(KSIZE, LDSOFT, PCOMPUTE, &
                       &PEXN, PLVFACT, PLSFACT, &
                       &PT,   PRRT, &
                       &PTHT, &
                       &PRRHONG_MR, PB_TH, PB_RR, PB_RG)
!!
!!**  PURPOSE
!!    -------
!!      Computes the RRHONG process
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
USE MODD_PARAM_ICE, ONLY : LFEEDBACKT
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL,                  INTENT(IN)    :: LDSOFT
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PEXN     ! Exner function
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLVFACT  ! L_v/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLSFACT  ! L_s/(Pi_ref*C_ph)
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRRT     ! Rain water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PTHT     ! Theta at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRRHONG_MR ! Mixing ratio change due to spontaneous freezing
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RR
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KSIZE) :: ZMASK
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JL
!
!-------------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('ICE4_RRHONG',0,ZHOOK_HANDLE)
!
!*       3.3     compute the spontaneous freezing source: RRHONG
!
PRRHONG_MR(:) = 0.
IF(.NOT. LDSOFT) THEN
  DO JL=1, KSIZE
    ZMASK(JL)=MAX(0., -SIGN(1., PT(JL)-(XTT-35.0))) * & ! PT(:)<XTT-35.0
             &MAX(0., -SIGN(1., XRTMIN(3)-PRRT(JL))) * & ! PRRT(:)>XRTMIN(3)
             &PCOMPUTE(JL)
    PRRHONG_MR(JL)=PRRT(JL) * ZMASK(JL)
  ENDDO
  IF(LFEEDBACKT) THEN
    !Limitation due to -35 crossing of temperature
    DO JL=1, KSIZE
      PRRHONG_MR(JL)=MIN(PRRHONG_MR(JL), MAX(0., ((XTT-35.)/PEXN(JL)-PTHT(JL))/(PLSFACT(JL)-PLVFACT(JL))))
    ENDDO
  ENDIF
ENDIF
DO JL=1, KSIZE
  PB_RG(JL) = PB_RG(JL) + PRRHONG_MR(JL)
  PB_RR(JL) = PB_RR(JL) - PRRHONG_MR(JL)
  PB_TH(JL) = PB_TH(JL) + PRRHONG_MR(JL)*(PLSFACT(JL)-PLVFACT(JL))
ENDDO
!
IF (LHOOK) CALL DR_HOOK('ICE4_RRHONG', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_RRHONG
