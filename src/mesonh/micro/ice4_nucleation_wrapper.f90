!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODI_ICE4_NUCLEATION_WRAPPER
INTERFACE
SUBROUTINE ICE4_NUCLEATION_WRAPPER(KIT, KJT,KKT, LDMASK, &
                           PTHT, PPABST, PRHODREF, PEXN, PLSFACT, PT, &
                           PRVT, &
                           PCIT, PRVHENI_MR)
IMPLICIT NONE
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
END SUBROUTINE ICE4_NUCLEATION_WRAPPER
END INTERFACE
END MODULE MODI_ICE4_NUCLEATION_WRAPPER
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
!
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,   ONLY: XTT

use mode_tools, only: Countjv

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
INTEGER                            :: IDX, JI, JJ, JK
INTEGER                            :: JL
INTEGER                            :: INEGT, INEGT_TMP
INTEGER, DIMENSION(:), ALLOCATABLE :: I1,I2,I3
LOGICAL                            :: GDSOFT
LOGICAL, DIMENSION(:), ALLOCATABLE :: GLDCOMPUTE
LOGICAL, DIMENSION(KIT,KJT,KKT)    :: GNEGT  ! Test where to compute the HEN process
REAL, DIMENSION(:), ALLOCATABLE    :: ZZT,       & ! Temperature
                                      ZPRES,      & ! Pressure
                                      ZRVT,       & ! Water vapor m.r. at t
                                      ZCIT,       & ! Pristine ice conc. at t
                                      ZTHT,       & ! Theta at t
                                      ZRHODREF,   &
                                      ZEXN,       &
                                      ZLSFACT,    &
                                      ZRVHENI_MR, &
                                      ZB_TH, ZB_RV, ZB_RI
!
!-------------------------------------------------------------------------------
!
!
!
!  optimization by looking for locations where
!  the temperature is negative only !!!
!
GNEGT(:,:,:)=PT(:,:,:)<XTT .AND. LDMASK
INEGT = COUNT(GNEGT(:,:,:))
!
ALLOCATE(GLDCOMPUTE(INEGT))
ALLOCATE(I1(INEGT),I2(INEGT),I3(INEGT))
ALLOCATE(ZZT(INEGT))
ALLOCATE(ZPRES(INEGT))
ALLOCATE(ZRVT(INEGT))
ALLOCATE(ZCIT(INEGT))
ALLOCATE(ZTHT(INEGT))
ALLOCATE(ZRHODREF(INEGT))
ALLOCATE(ZEXN(INEGT))
ALLOCATE(ZLSFACT(INEGT))
ALLOCATE(ZRVHENI_MR(INEGT))
ALLOCATE(ZB_TH(INEGT))
ALLOCATE(ZB_RV(INEGT))
ALLOCATE(ZB_RI(INEGT))
!
IF(INEGT>0) INEGT_TMP=COUNTJV(GNEGT(:,:,:), I1(:), I2(:), I3(:))
!
PRVHENI_MR(:,:,:)=0.
IF(INEGT>0) THEN
  DO JL=1, INEGT
    ZRVT(JL)=PRVT(I1(JL), I2(JL), I3(JL))
    ZCIT(JL)=PCIT(I1(JL), I2(JL), I3(JL))
    ZZT(JL)=PT(I1(JL), I2(JL), I3(JL))
    ZPRES(JL)=PPABST(I1(JL), I2(JL), I3(JL))
    ZTHT(JL)=PTHT(I1(JL), I2(JL), I3(JL))
    ZRHODREF(JL)=PRHODREF(I1(JL), I2(JL), I3(JL))
    ZEXN(JL)=PEXN(I1(JL), I2(JL), I3(JL))
    ZLSFACT(JL)=PLSFACT(I1(JL), I2(JL), I3(JL))
  ENDDO
  GDSOFT = .FALSE.
  GLDCOMPUTE(:) = ZZT(:)<XTT
  ZB_TH(:) = 0.
  ZB_RV(:) = 0.
  ZB_RI(:) = 0.
  CALL ICE4_NUCLEATION(INEGT, GDSOFT, GLDCOMPUTE, &
                       ZTHT, ZPRES, ZRHODREF, ZEXN, ZLSFACT, ZZT, &
                       ZRVT, &
                       ZCIT, ZRVHENI_MR, ZB_TH, ZB_RV, ZB_RI)
  PRVHENI_MR(:,:,:)= 0.0
  DO JL=1, INEGT
    PRVHENI_MR(I1(JL), I2(JL), I3(JL)) = ZRVHENI_MR(JL)
    PCIT      (I1(JL), I2(JL), I3(JL)) = ZCIT      (JL)
  END DO
END IF
!
DEALLOCATE(GLDCOMPUTE)
DEALLOCATE(I1,I2,I3)
DEALLOCATE(ZZT,ZPRES,ZRVT,ZCIT,ZTHT,ZRHODREF,ZEXN,ZLSFACT,ZRVHENI_MR,ZB_TH,ZB_RV,ZB_RI)
!
END SUBROUTINE ICE4_NUCLEATION_WRAPPER
