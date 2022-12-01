!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RSRIMCG_OLD
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_RSRIMCG_OLD(CST, ICEP, ICED, KPROMA, KSIZE, LDSOFT, LDCOMPUTE, &
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
USE MODD_RAIN_ICE_DESCR, ONLY: RAIN_ICE_DESCR_t
USE MODD_RAIN_ICE_PARAM, ONLY: RAIN_ICE_PARAM_t
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(CST_t),              INTENT(IN)    :: CST
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
REAL, DIMENSION(KPROMA) :: ZVEC1, ZVEC2
INTEGER, DIMENSION(KPROMA) :: IVEC1, IVEC2
REAL, DIMENSION(KPROMA) :: ZZW
INTEGER :: JL
REAL(KIND=JPRB) :: ZHOOK_HANDLE
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
  IGRIM = 0
  DO JL = 1, KSIZE
    IF(PRCT(JL)>ICED%XRTMIN(2) .AND. PRST(JL)>ICED%XRTMIN(5) .AND. LDCOMPUTE(JL) .AND. PT(JL)<CST%XTT) THEN
      IGRIM = IGRIM + 1
      IVEC1(IGRIM) = JL
      GRIM(JL) = .TRUE.
    ELSE
      GRIM(JL) = .FALSE.
    ENDIF
  ENDDO
  !
  IF(IGRIM>0) THEN
    !
    !        5.1.1  select the PLBDAS
    !
    DO JL = 1, IGRIM
      ZVEC1(JL) = PLBDAS(IVEC1(JL))
    ENDDO
    !
    !        5.1.2  find the next lower indice for the PLBDAS in the geometrical
    !               set of Lbda_s used to tabulate some moments of the incomplete
    !               gamma function
    !
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN(REAL(ICEP%NGAMINC)-0.00001,           &
                          ICEP%XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + ICEP%XRIMINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )

    !
    !        5.1.5  perform the linear interpolation of the normalized
    !               "XBS"-moment of the incomplete gamma function (XGAMINC_RIM2)
    !
    ZVEC1(1:IGRIM) =  ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - ICEP%XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = 0.
    DO JL = 1, IGRIM
      ZZW(IVEC1(JL)) = ZVEC1(JL)
    ENDDO

    !
    !        5.1.6  riming-conversion of the large sized aggregates into graupeln
    !
    !
    !$mnh_expand_where(JL=1:KSIZE)
    WHERE(GRIM(1:KSIZE))
      PRSRIMCG_MR(1:KSIZE) = ICEP%XSRIMCG * PLBDAS(1:KSIZE)**ICEP%XEXSRIMCG   & ! RSRIMCG
#if defined(REPRO48) || defined(REPRO55)
                               * (1.0 - ZZW(1:KSIZE) )/PRHODREF(1:KSIZE)
#else
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
END SUBROUTINE ICE4_RSRIMCG_OLD
END MODULE MODE_ICE4_RSRIMCG_OLD
