MODULE MODI_DEEP_CONVECTION
IMPLICIT NONE
INTERFACE
    SUBROUTINE DEEP_CONVECTION( KLON, KLEV, KIDIA, KFDIA, KBDIA, KTDIA,        &
                                PDTCONV, KICE, OREFRESH, ODOWN, OSETTADJ,      &
                                PPABST, PZZ, PDXDY, PTIMEC,                    &
                                PTT, PRVT, PRCT, PRIT, PUT, PVT, PWT,          &
                                KCOUNT, PTTEN, PRVTEN, PRCTEN, PRITEN,         &
                                PPRLTEN, PPRSTEN,                              &
                                KCLTOP, KCLBAS, PPRLFLX, PPRSFLX,              &
                                PUMF, PDMF, PCAPE,                             &
                                OCH1CONV, KCH1, PCH1, PCH1TEN,                 &
                                OUSECHEM, OCH_CONV_SCAV, OCH_CONV_LINOX,       &
                                ODUST, OSALT, PRHODREF, PIC_RATE, PCG_RATE     )
!
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
INTEGER,                    INTENT(IN) :: KIDIA    ! value of the first point in x
INTEGER,                    INTENT(IN) :: KFDIA    ! value of the last point in x
INTEGER,                    INTENT(IN) :: KBDIA    ! vertical  computations start at
!                                                  ! KBDIA that is at least 1
INTEGER,                    INTENT(IN) :: KTDIA    ! vertical computations can be
                                                   ! limited to KLEV + 1 - KTDIA
                                                   ! default=1
REAL,                       INTENT(IN) :: PDTCONV  ! Interval of time between two
                                                   ! calls of the deep convection
                                                   ! scheme
INTEGER,                    INTENT(IN) :: KICE     ! flag for ice ( 1 = yes,
                                                   !                0 = no ice )
LOGICAL,                    INTENT(IN) :: OREFRESH ! refresh or not tendencies
                                                   ! at every call
LOGICAL,                    INTENT(IN) :: ODOWN    ! take or not convective
                                                   ! downdrafts into account
LOGICAL,                    INTENT(IN) :: OSETTADJ ! logical to set convective
                                                   ! adjustment time by user
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PTT      ! grid scale temperature at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRVT     ! grid scale water vapor "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRCT     ! grid scale r_c  "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRIT     ! grid scale r_i "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PUT      ! grid scale horiz. wind u "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PVT      ! grid scale horiz. wind v "
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PWT      ! grid scale vertical
                                                   ! velocity (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PPABST   ! grid scale pressure at t
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PZZ      ! height of model layer (m)
REAL, DIMENSION(KLON),      INTENT(IN) :: PDXDY    ! horizontal grid area (m**2)
REAL, DIMENSION(KLON),      INTENT(IN) :: PTIMEC   ! value of convective adjustment
                                                   ! time if OSETTADJ=.TRUE.
!
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCOUNT ! convective counter (recompute
                                                   ! tendency or keep it)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PTTEN  ! convective temperature
                                                   ! tendency (K/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRVTEN ! convective r_v tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRCTEN ! convective r_c tendency (1/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PRITEN ! convective r_i tendency (1/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRLTEN! liquid surf. precipitation
                                                   ! tendency (m/s)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PPRSTEN! solid surf. precipitation
                                                   ! tendency (m/s)
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLTOP ! cloud top level
INTEGER, DIMENSION(KLON),   INTENT(INOUT):: KCLBAS ! cloud base level
                                                   ! they are given a value of
                                                   ! 0 if no convection
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRLFLX! liquid precip flux (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PPRSFLX! solid  precip flux (m/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PUMF   ! updraft mass flux (kg/s m2)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT):: PDMF   ! downdraft mass flux (kg/s m2)
REAL, DIMENSION(KLON),      INTENT(INOUT):: PCAPE  ! maximum CAPE (J/kg)
!
LOGICAL,                    INTENT(IN) :: OCH1CONV ! include tracer transport
INTEGER,                    INTENT(IN) :: KCH1     ! number of species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(IN) :: PCH1! grid scale chemical species
REAL, DIMENSION(KLON,KLEV,KCH1), INTENT(INOUT):: PCH1TEN! species conv. tendency (1/s)
LOGICAL,                    INTENT(IN) :: OUSECHEM      ! flag for chemistry
LOGICAL,                    INTENT(IN) :: OCH_CONV_SCAV !  & scavenging
LOGICAL,                    INTENT(IN) :: OCH_CONV_LINOX ! & LiNOx
LOGICAL,                    INTENT(IN) :: ODUST         ! flag for dust
LOGICAL,                    INTENT(IN) :: OSALT         ! flag for sea salt
REAL, DIMENSION(KLON,KLEV), INTENT(IN) :: PRHODREF      ! grid scale density
REAL, DIMENSION(KLON), INTENT(INOUT) :: PIC_RATE ! IC lightning frequency
REAL, DIMENSION(KLON), INTENT(INOUT) :: PCG_RATE ! CG lightning frequency
!
END SUBROUTINE DEEP_CONVECTION
!
END INTERFACE
!
END MODULE MODI_DEEP_CONVECTION
