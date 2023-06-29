!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RSRIMCG_OLD
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_RSRIMCG_OLD(CST, PARAMI, ICEP, ICED, KPROMA, KSIZE, LDSOFT, LDCOMPUTE, &
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
INTEGER, INTENT(IN) :: KPROMA, KSIZE
LOGICAL,                       INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KPROMA),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KPROMA),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KPROMA),       INTENT(OUT)   :: PRSRIMCG_MR ! Mr change due to cloud droplet riming of the aggregates
!
!*       0.2  declaration of local variables
!
LOGICAL, DIMENSION(KPROMA) :: GRIM
INTEGER :: IGRIM
REAL, DIMENSION(KPROMA) :: ZBUF1, ZBUF2
INTEGER, DIMENSION(KPROMA) :: IBUF1, IBUF2
REAL, DIMENSION(KPROMA) :: ZZW
INTEGER :: JL
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
  DO JL = 1, KSIZE
    GRIM(JL)=PRCT(JL)>ICED%XRTMIN(2) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL) .AND. PT(JL)<CST%XTT
  ENDDO
  CALL INTERP_MICRO_1D(KPROMA, KSIZE, PLBDAS(:), ICEP%NGAMINC, ICEP%XRIMINTP1, ICEP%XRIMINTP2, &
                      &PARAMI%LPACK_INTERP, GRIM(:), IBUF1, IBUF2, ZBUF1, ZBUF2, &
                      &IGRIM, &
                      &ICEP%XGAMINC_RIM2, ZZW)
  !
  IF(IGRIM>0) THEN
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GRIM(1:KSIZE))
#ifdef REPRO48
      PRSRIMCG_MR(1:KSIZE) = ICEP%XSRIMCG * PLBDAS(1:KSIZE)**ICEP%XEXSRIMCG   & ! RSRIMCG
                               * (1.0 - ZZW(1:KSIZE) )/PRHODREF(1:KSIZE)
#else
      PRSRIMCG_MR(1:KSIZE) = ICEP%XSRIMCG * PLBDAS(1:KSIZE)**ICEP%XEXSRIMCG   & ! RSRIMCG
                               * (1.0 - ZZW(1:KSIZE) )*PRST(1:KSIZE)
#endif
      PRSRIMCG_MR(1:KSIZE)=MIN(PRST(1:KSIZE), PRSRIMCG_MR(1:KSIZE))
    END WHERE
    !$mnh_end_expand_where(JL=1:KSIZE)
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
