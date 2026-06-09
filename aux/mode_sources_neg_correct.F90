MODULE MODE_SOURCES_NEG_CORRECT
IMPLICIT NONE
CONTAINS
SUBROUTINE SOURCES_NEG_CORRECT_PHY(D, KSV, HCLOUD, HELEC, HBUDNAME, KRR, PTSTEP, PPABST, &
                              &PTHT, PRT, PRTHS, PRRS, PRSVS, PRHODJ)
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
TYPE(DIMPHYEX_t),            INTENT(IN)           :: D
INTEGER,                     INTENT(IN)           :: KSV      ! Number of SV variables
CHARACTER(lEN=*),            INTENT(IN)           :: HELEC    ! Kind of cloud electricity parameterization
CHARACTER(LEN=*),            INTENT(IN)           :: HCLOUD   ! Kind of cloud parameterization
CHARACTER(LEN=*),            INTENT(IN)           :: HBUDNAME ! Budget name
INTEGER,                     INTENT(IN)           :: KRR      ! Number of moist variables
REAL,                        INTENT(IN)           :: PTSTEP   ! Timestep
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)           :: PPABST   ! Absolute pressure at time t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN)           :: PTHT     ! Theta at time t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(IN)           :: PRT      ! Moist variables at time t
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(INOUT)        :: PRTHS    ! Source terms
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KRR), INTENT(INOUT)        :: PRRS     ! Source terms
REAL, DIMENSION(D%NIT,D%NJT,D%NKT,KSV), INTENT(INOUT)        :: PRSVS    ! Source terms
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),    INTENT(IN), OPTIONAL :: PRHODJ   ! Dry density * jacobian
END SUBROUTINE SOURCES_NEG_CORRECT_PHY
END MODULE MODE_SOURCES_NEG_CORRECT
