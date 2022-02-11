!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_NUCLEATION_WRAPPER
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_NUCLEATION_WRAPPER(KIT, KJT, KKT, LDMASK, &
                                   PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                                   PRVT, &
                                   PCIT, PRVHENI_MR)
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
!  P. Wautelet 28/05/2019: move COUNTJV function to tools.f90
!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!!     R. El Khatib 24-Aug-2021 Optimizations
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,   ONLY: XTT
USE MODE_TOOLS, ONLY: COUNTJV
USE MODE_ICE4_NUCLEATION, ONLY: ICE4_NUCLEATION
USE PARKIND1,   ONLY : JPRB
USE YOMHOOK ,   ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,                        INTENT(IN)    :: KIT, KJT, KKT
LOGICAL, DIMENSION(KIT,KJT,KKT),INTENT(IN)    :: LDMASK
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PTHT    ! Theta at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PPABST  ! absolute pressure at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRHODREF! Reference density
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PEXN    ! Exner function
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PLSFACT
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PT      ! Temperature at time t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(IN)    :: PRVT    ! Water vapor m.r. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(INOUT) :: PCIT    ! Pristine ice n.c. at t
REAL, DIMENSION(KIT,KJT,KKT),   INTENT(OUT)   :: PRVHENI_MR ! Mixing ratio change due to the heterogeneous nucleation
!
!*       0.2  declaration of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER                            :: JL
INTEGER                            :: INEGT
INTEGER, DIMENSION(COUNT(PT<XTT .AND. LDMASK)) :: I1,I2,I3
LOGICAL, DIMENSION(COUNT(PT<XTT .AND. LDMASK))  :: GLDCOMPUTE ! computation criterium
LOGICAL, DIMENSION(KIT, KJT, KKT)  :: GNEGT  ! Test where to compute the HEN process
REAL, DIMENSION(COUNT(PT<XTT .AND. LDMASK))  :: ZZT,      & ! Temperature
                                      ZPRES,      & ! Pressure
                                      ZRVT,       & ! Water vapor m.r. at t
                                      ZCIT,       & ! Pristine ice conc. at t
                                      ZTHT,       & ! Theta at t
                                      ZRHODREF,   &
                                      ZEXN,       &
                                      ZLSFACT,    &
                                      ZRVHENI_MR
!! MNH version INTEGER, DIMENSION(:), ALLOCATABLE :: I1,I2,I3
!! MNH version LOGICAL, DIMENSION(:), ALLOCATABLE :: GLDCOMPUTE
!! MNH version LOGICAL, DIMENSION(KIT,KJT,KKT)    :: GNEGT  ! Test where to compute the HEN process
!! MNH version REAL, DIMENSION(:), ALLOCATABLE    :: ZZT,       & ! Temperature
!! MNH version                                       ZPRES,      & ! Pressure
!! MNH version                                       ZRVT,       & ! Water vapor m.r. at t
!! MNH version                                       ZCIT,       & ! Pristine ice conc. at t
!! MNH version                                       ZTHT,       & ! Theta at t
!! MNH version                                       ZRHODREF,   &
!! MNH version                                       ZEXN,       &
!! MNH version                                       ZLSFACT,    &
!! MNH version                                       ZRVHENI_MR
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_NUCLEATION_WRAPPER', 0, ZHOOK_HANDLE)!
!
!
!  optimization by looking for locations where
!  the temperature is negative only !!!
!
GNEGT(:,:,:)=PT(:,:,:)<XTT .AND. LDMASK
INEGT = COUNT(GNEGT(:,:,:))
!
!! MNH version ALLOCATE(GLDCOMPUTE(INEGT))
!! MNH version ALLOCATE(I1(INEGT),I2(INEGT),I3(INEGT))
!! MNH version ALLOCATE(ZZT(INEGT))
!! MNH version ALLOCATE(ZPRES(INEGT))
!! MNH version ALLOCATE(ZRVT(INEGT))
!! MNH version ALLOCATE(ZCIT(INEGT))
!! MNH version ALLOCATE(ZTHT(INEGT))
!! MNH version ALLOCATE(ZRHODREF(INEGT))
!! MNH version ALLOCATE(ZEXN(INEGT))
!! MNH version ALLOCATE(ZLSFACT(INEGT))
!! MNH version ALLOCATE(ZRVHENI_MR(INEGT))
!
IF(INEGT>0) INEGT=COUNTJV(GNEGT(:,:,:), I1(:), I2(:), I3(:))
!
PRVHENI_MR(:,:,:)=0.
IF(INEGT>0) THEN
  DO JL=1, INEGT
    ZRVT(JL)=PRVT(I1(JL), I2(JL), I3(JL))
    ZCIT(JL)=PCIT(I1(JL), I2(JL), I3(JL))
    ZPRES(JL)=PPABST(I1(JL), I2(JL), I3(JL))
    ZTHT(JL)=PTHT(I1(JL), I2(JL), I3(JL))
    ZRHODREF(JL)=PRHODREF(I1(JL), I2(JL), I3(JL))
    ZEXN(JL)=PEXN(I1(JL), I2(JL), I3(JL))
    ZLSFACT(JL)=PLSFACT(I1(JL), I2(JL), I3(JL)) / ZEXN(JL)
    ZZT(JL)=PT(I1(JL), I2(JL), I3(JL))
    GLDCOMPUTE(JL)=ZZT(JL)<XTT
  ENDDO
  CALL ICE4_NUCLEATION(INEGT, GLDCOMPUTE, &
                       ZTHT, ZPRES, ZRHODREF, ZEXN, ZLSFACT, ZZT, &
                       ZRVT, &
                       ZCIT, ZRVHENI_MR)
  DO JL=1, INEGT
    PRVHENI_MR(I1(JL), I2(JL), I3(JL)) = ZRVHENI_MR(JL)
    PCIT      (I1(JL), I2(JL), I3(JL)) = ZCIT      (JL)
  END DO
END IF
!
!! MNH versionDEALLOCATE(GLDCOMPUTE)
!! MNH versionDEALLOCATE(I1,I2,I3)
!! MNH versionDEALLOCATE(ZZT,ZPRES,ZRVT,ZCIT,ZTHT,ZRHODREF,ZEXN,ZLSFACT,ZRVHENI_MR)
!
IF (LHOOK) CALL DR_HOOK('ICE4_NUCLEATION_WRAPPER', 1, ZHOOK_HANDLE)

END SUBROUTINE ICE4_NUCLEATION_WRAPPER
END MODULE MODE_ICE4_NUCLEATION_WRAPPER
