MODULE MODE_SOURCES_NEG_CORRECT
IMPLICIT NONE
CONTAINS
SUBROUTINE SOURCES_NEG_CORRECT(HCLOUD, HBUDNAME, KRR, PTSTEP, PPABST, &
                              &PTHT, PRT, PRTHS, PRRS, PRSVS, PRHODJ)
IMPLICIT NONE
CHARACTER(LEN=*),            INTENT(IN)           :: HCLOUD   ! Kind of cloud parameterization
CHARACTER(LEN=*),            INTENT(IN)           :: HBUDNAME ! Budget name
INTEGER,                     INTENT(IN)           :: KRR      ! Number of moist variables
REAL,                        INTENT(IN)           :: PTSTEP   ! Timestep
REAL, DIMENSION(:, :, :),    INTENT(IN)           :: PPABST   ! Absolute pressure at time t
REAL, DIMENSION(:, :, :),    INTENT(IN)           :: PTHT     ! Theta at time t
REAL, DIMENSION(:, :, :, :), INTENT(IN)           :: PRT      ! Moist variables at time t
REAL, DIMENSION(:, :, :),    INTENT(INOUT)        :: PRTHS    ! Source terms
REAL, DIMENSION(:, :, :, :), INTENT(INOUT)        :: PRRS     ! Source terms
REAL, DIMENSION(:, :, :, :), INTENT(INOUT)        :: PRSVS    ! Source terms
REAL, DIMENSION(:, :, :),    INTENT(IN), OPTIONAL :: PRHODJ   ! Dry density * jacobian
END SUBROUTINE SOURCES_NEG_CORRECT
END MODULE MODE_SOURCES_NEG_CORRECT
