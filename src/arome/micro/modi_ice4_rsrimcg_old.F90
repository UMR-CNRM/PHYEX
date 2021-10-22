MODULE MODI_ICE4_RSRIMCG_OLD
INTERFACE
SUBROUTINE ICE4_RSRIMCG_OLD(KSIZE, LDSOFT, LDCOMPUTE, &
                           &PRHODREF, &
                           &PLBDAS, &
                           &PT, PRCT, PRST, &
                           &PRSRIMCG_MR, PB_RS, PB_RG)
IMPLICIT NONE
INTEGER, INTENT(IN) :: KSIZE
LOGICAL,                      INTENT(IN)    :: LDSOFT
LOGICAL, DIMENSION(KSIZE),    INTENT(IN)    :: LDCOMPUTE
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRHODREF ! Reference density
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PLBDAS   ! Slope parameter of the aggregate distribution
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PT       ! Temperature
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRCT     ! Cloud water m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(IN)    :: PRST     ! Snow/aggregate m.r. at t
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PRSRIMCG_MR ! Mr change due to cloud droplet riming of the aggregates
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RS
REAL, DIMENSION(KSIZE),       INTENT(INOUT) :: PB_RG
END SUBROUTINE ICE4_RSRIMCG_OLD
END INTERFACE
END MODULE MODI_ICE4_RSRIMCG_OLD
