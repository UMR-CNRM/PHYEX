!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RSRIMCG_OLD
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_RSRIMCG_OLD(CST, PARAMI, ICEP, ICED, D, LDSOFT, LDCOMPUTE, &
                           &PRHODREF, &
                           &PLBDAS, &
                           &PT, PRCT, PRST, &
                           &PRSRIMCG_MR)
!!
!!**  PURPOSE
!!    -------
!!      Computes the riming-conversion of the large sized aggregates into graupel
!!
!!    AUTHOR
!!    ------
!!      S. Riette from the splitting of rain_ice source code (nov. 2014)
!!
!!    MODIFICATIONS
!!    -------------
!!
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 29/05/2019: remove PACK/UNPACK intrinsics (to get more performance and better OpenACC support)
!
!
!*      0. DECLARATIONS
!          ------------
!
USE MODD_DIMPHYEX,       ONLY: DIMPHYEX_t
USE MODD_CST,            ONLY: CST_t
USE MODD_PARAM_ICE_n,      ONLY: PARAM_ICE_t
USE MODD_RAIN_ICE_DESCR_n, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM_n, ONLY: RAIN_ICE_PARAM_t
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
TYPE(PARAM_ICE_t),        INTENT(IN)    :: PARAMI
TYPE(RAIN_ICE_PARAM_t),   INTENT(IN)    :: ICEP
TYPE(RAIN_ICE_DESCR_t),   INTENT(IN)    :: ICED
TYPE(DIMPHYEX_t),INTENT(IN) :: D
LOGICAL,                       INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(D%NIJT),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(D%NIJT),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(D%NIJT),       INTENT(OUT)   :: PRSRIMCG_MR ! Mr change due to cloud droplet riming of the aggregates
!
!*       0.2  declaration of local variables
!
LOGICAL, DIMENSION(D%NIJT) :: GRIM
INTEGER :: IGRIM
REAL, DIMENSION(D%NIJT) :: ZBUF1, ZBUF2
INTEGER, DIMENSION(D%NIJT) :: IBUF1, IBUF2
REAL, DIMENSION(D%NIJT) :: ZZW
INTEGER :: JIJ
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('ICE4_RSRIMCG_OLD', 0, ZHOOK_HANDLE)
!
!-------------------------------------------------------------------------------
!
!*       5.1    cloud droplet riming of the aggregates
!

PRSRIMCG_MR(:)=0.

!
IF(.NOT. LDSOFT) THEN


  DO JIJ=D%NIJB, D%NIJE
    GRIM(JIJ)=PRCT(JIJ)>ICED%XRTMIN(2) .AND. PRST(JIJ)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JIJ) .AND. PT(JIJ)<CST%XTT
  ENDDO

  CALL INTERP_MICRO_1D(D, PLBDAS, ICEP%NGAMINC, ICEP%XRIMINTP1, ICEP%XRIMINTP2, &
                      &PARAMI%LPACK_INTERP, GRIM, IBUF1, IBUF2, ZBUF1, ZBUF2, &
                      &IGRIM, &
                      &ICEP%XGAMINC_RIM2, ZZW)
  !
  IF(IGRIM>0) THEN

    DO JIJ=D%NIJB, D%NIJE
      IF(.NOT. ICEP%LNEWCOEFF) THEN
        IF (GRIM(JIJ)) THEN
          PRSRIMCG_MR(JIJ) = ICEP%XSRIMCG * PLBDAS(JIJ)**ICEP%XEXSRIMCG   & ! RSRIMCG
                                 * (1.0 - ZZW(JIJ) )/PRHODREF(JIJ)
        END IF
      ELSE
        IF (GRIM(JIJ)) THEN
          PRSRIMCG_MR(JIJ) = ICEP%XSRIMCG * PLBDAS(JIJ)**ICEP%XEXSRIMCG   & ! RSRIMCG
                                 * (1.0 - ZZW(JIJ) )*PRST(JIJ)
        END IF
      ENDIF
      IF (GRIM(JIJ)) THEN
        PRSRIMCG_MR(JIJ)=MIN(PRST(JIJ), PRSRIMCG_MR(JIJ))
      END IF
    END DO

  END IF
ENDIF
!
IF (LHOOK) CALL DR_HOOK('ICE4_RSRIMCG_OLD', 1, ZHOOK_HANDLE)
!
CONTAINS
!
INCLUDE "interp_micro.func.h"
!
END SUBROUTINE ICE4_RSRIMCG_OLD
END MODULE MODE_ICE4_RSRIMCG_OLD
