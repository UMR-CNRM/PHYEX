!MNH_LIC Copyright 1994-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
MODULE MODE_ICE4_RSRIMCG_OLD
IMPLICIT NONE
CONTAINS
SUBROUTINE ICE4_RSRIMCG_OLD(KSIZE, LDSOFT, LDCOMPUTE, &
                           &PRHODREF, &
                           &PLBDAS, &
                           &PT, PRCT, PRST, &
                           &PRSRIMCG_MR, PB_RS, PB_RG)
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
USE MODD_CST,            ONLY: XTT
USE MODD_PARAM_ICE,      ONLY: CSNOWRIMING
USE MODD_RAIN_ICE_DESCR, ONLY: XRTMIN
USE MODD_RAIN_ICE_PARAM, ONLY: NGAMINC, XEXSRIMCG, XGAMINC_RIM2, XRIMINTP1, XRIMINTP2, XSRIMCG
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER, INTENT(IN) :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(OUT)   :: PRSRIMCG_MR ! Mr change due to cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
!
!*       0.2  declaration of local variables
!
LOGICAL, DIMENSION(KSIZE) :: GRIM
INTEGER :: IGRIM
REAL, DIMENSION(KSIZE) :: ZVEC1, ZVEC2
INTEGER, DIMENSION(KSIZE) :: IVEC1, IVEC2
REAL, DIMENSION(KSIZE) :: ZZW
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
  GRIM(:) = .FALSE.
  DO JL = 1, SIZE(GRIM)
    IF(PRCT(JL)>XRTMIN(2) .AND. PRST(JL)>XRTMIN(5) .AND. LDCOMPUTE(JL) .AND. PT(JL)<XTT) THEN
      IGRIM = IGRIM + 1
      IVEC1(IGRIM) = JL
      GRIM(JL) = .TRUE.
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
    ZVEC2(1:IGRIM) = MAX( 1.00001, MIN(REAL(NGAMINC)-0.00001,           &
                          XRIMINTP1 * LOG( ZVEC1(1:IGRIM) ) + XRIMINTP2 ) )
    IVEC2(1:IGRIM) = INT( ZVEC2(1:IGRIM) )
    ZVEC2(1:IGRIM) = ZVEC2(1:IGRIM) - REAL( IVEC2(1:IGRIM) )

    !
    !        5.1.5  perform the linear interpolation of the normalized
    !               "XBS"-moment of the incomplete gamma function (XGAMINC_RIM2)
    !
    ZVEC1(1:IGRIM) =  XGAMINC_RIM2( IVEC2(1:IGRIM)+1 )* ZVEC2(1:IGRIM)      &
                    - XGAMINC_RIM2( IVEC2(1:IGRIM)   )*(ZVEC2(1:IGRIM) - 1.0)
    ZZW(:) = 0.
    DO JL = 1, IGRIM
      ZZW(IVEC1(JL)) = ZVEC1(JL)
    ENDDO

    !
    !        5.1.6  riming-conversion of the large sized aggregates into graupeln
    !
    !
    WHERE(GRIM(:))
      PRSRIMCG_MR(:) = XSRIMCG * PLBDAS(:)**XEXSRIMCG   & ! RSRIMCG
                               * (1.0 - ZZW(:) )/PRHODREF(:)
      PRSRIMCG_MR(:)=MIN(PRST(:), PRSRIMCG_MR(:))
    END WHERE
  END IF
ENDIF
PB_RS(:) = PB_RS(:) - PRSRIMCG_MR(:)
PB_RG(:) = PB_RG(:) + PRSRIMCG_MR(:)
!
IF (LHOOK) CALL DR_HOOK('ICE4_RSRIMCG_OLD', 1, ZHOOK_HANDLE)
!
END SUBROUTINE ICE4_RSRIMCG_OLD
END MODULE MODE_ICE4_RSRIMCG_OLD
