!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_NUCLEATION
INTERFACE
SUBROUTINE ICE4_NUCLEATION(KSIZE, ODSOFT, ODCOMPUTE, &
                           PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                           PRVT, &
                           PCIT, PRVHENI_MR, PB_TH, PB_RV, PB_RI)
IMPLICIT NONE
INTEGER,                  INTENT(IN)    :: KSIZE
LOGICAL,                  INTENT(IN)    :: ODSOFT
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT      ! Temperature at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: PRVHENI_MR ! Mixing ratio change due to the heterogeneous nucleation
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_RV
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_RI
END SUBROUTINE ICE4_NUCLEATION
END INTERFACE
END MODULE MODI_ICE4_NUCLEATION
SUBROUTINE ICE4_NUCLEATION(KSIZE, ODSOFT, ODCOMPUTE, &
                           PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                           PRVT, &
                           PCIT, PRVHENI_MR, PB_TH, PB_RV, PB_RI)
!!
!!**  PURPOSE
!!    -------
!!      Computes the nucleation
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
USE MODD_CST,            ONLY: XALPI,XALPW,XBETAI,XBETAW,XGAMI,XGAMW,XMD,XMV,XTT,XEPSILO
USE MODD_PARAM_ICE,      ONLY: LFEEDBACKT
USE MODD_RAIN_ICE_PARAM, ONLY: XALPHA1,XALPHA2,XBETA1,XBETA2,XMNU0,XNU10,XNU20
USE MODD_RAIN_ICE_DESCR, ONLY: XRTMIN
!
USE MODE_MPPDB
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                  INTENT(IN)    :: KSIZE
LOGICAL,                  INTENT(IN)    :: ODSOFT
LOGICAL, DIMENSION(KSIZE),INTENT(IN)    :: ODCOMPUTE
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PT      ! Temperature at time t
REAL, DIMENSION(KSIZE),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KSIZE),   INTENT(OUT)   :: PRVHENI_MR ! Mixing ratio change due to the heterogeneous nucleation
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_TH
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_RV
REAL, DIMENSION(KSIZE),   INTENT(INOUT) :: PB_RI
!
!*       0.2  declaration of local variables
!
REAL, DIMENSION(KSIZE) :: ZW ! work array
LOGICAL, DIMENSION(KSIZE) :: GNEGT  ! Test where to compute the HEN process
REAL, DIMENSION(KSIZE)    :: ZZW,      & ! Work array
                             ZUSW,     & ! Undersaturation over water
                             ZSSI        ! Supersaturation over ice
!-------------------------------------------------------------------------------
!
!
PRVHENI_MR(:)=0.
IF(.NOT. ODSOFT) THEN
  GNEGT(:)=PT(:)<XTT .AND. PRVT>XRTMIN(1) .AND. ODCOMPUTE(:)
  PRVHENI_MR(:)=0.
  ZSSI(:)=0.
  ZUSW(:)=0.
  ZZW(:)=0.
  WHERE(GNEGT(:))
    ZZW(:)=ALOG(PT(:))
    ZUSW(:)=EXP(XALPW - XBETAW/PT(:) - XGAMW*ZZW(:))          ! es_w
    ZZW(:)=EXP(XALPI - XBETAI/PT(:) - XGAMI*ZZW(:))           ! es_i
  END WHERE
  WHERE(GNEGT(:))
    ZZW(:)=MIN(PPABST(:)/2., ZZW(:))             ! safety limitation
    ZSSI(:)=PRVT(:)*(PPABST(:)-ZZW(:)) / (XEPSILO*ZZW(:)) - 1.0
                                                 ! Supersaturation over ice
    ZUSW(:)=MIN(PPABST(:)/2., ZUSW(:))            ! safety limitation
    ZUSW(:)=(ZUSW(:)/ZZW(:))*((PPABST(:)-ZZW(:))/(PPABST(:)-ZUSW(:))) - 1.0
                               ! Supersaturation of saturated water vapor over ice
    !
    !*       3.1     compute the heterogeneous nucleation source RVHENI
    !
    !*       3.1.1   compute the cloud ice concentration
    !
    ZSSI(:)=MIN(ZSSI(:), ZUSW(:)) ! limitation of SSi according to SSw=0
  END WHERE
  ZZW(:)=0.
  WHERE(GNEGT(:) .AND. PT(:)<XTT-5.0 .AND. ZSSI(:)>0.0 )
    ZZW(:)=XNU20*EXP(XALPHA2*ZSSI(:)-XBETA2)
  ELSEWHERE(GNEGT(:) .AND. PT(:)<=XTT-2.0 .AND. PT(:)>=XTT-5.0 .AND. ZSSI(:)>0.0)
    ZZW(:)=MAX(XNU20*EXP(-XBETA2 ), &
               XNU10*EXP(-XBETA1*(PT(:)-XTT))*(ZSSI(:)/ZUSW(:))**XALPHA1)
  END WHERE
  WHERE(GNEGT(:))
    ZZW(:)=ZZW(:)-PCIT(:)
    ZZW(:)=MIN(ZZW(:), 50.E3) ! limitation provisoire a 50 l^-1
  END WHERE
  WHERE(GNEGT(:))
    !
    !*       3.1.2   update the r_i and r_v mixing ratios
    !
    PRVHENI_MR(:)=MAX(ZZW(:), 0.0)*XMNU0/PRHODREF(:)
    PRVHENI_MR(:)=MIN(PRVT(:), PRVHENI_MR(:))
  END WHERE
  !Limitation due to 0 crossing of temperature
  IF(LFEEDBACKT) THEN
    ZW(:)=0.
    WHERE(GNEGT(:))
      ZW(:)=MIN(PRVHENI_MR(:), &
                MAX(0., (XTT/PEXN(:)-PTHT(:))/PLSFACT(:))) / &
                MAX(PRVHENI_MR(:), 1.E-20)
    END WHERE
  ELSE
    ZW(:)=1.
  ENDIF
  PRVHENI_MR(:)=PRVHENI_MR(:)*ZW(:)
  PCIT(:)=MAX(ZZW(:)*ZW(:)+PCIT(:), PCIT(:))
  !
  PB_RI(:)=PB_RI(:) + PRVHENI_MR(:)
  PB_RV(:)=PB_RV(:) - PRVHENI_MR(:)
  PB_TH(:)=PB_TH(:) + PRVHENI_MR(:)*PLSFACT(:)
ENDIF
!
END SUBROUTINE ICE4_NUCLEATION
