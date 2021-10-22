INTERFACE
      SUBROUTINE ARO_MNHC(PRSVSIN, PRHODREF, PTSTEP,    &
                          PTHT, PABST,                &
                          PRT, PLAT, PLON,            &
                          PALB_UV, PZS, PZENITH, PZZ, &
                          KYEAR, KMONTH, KDAY, PTIME, &
                          KLON,KLEV,NSV, KRR,         &
                          KTCOUNT, KLUOUT,NDIAG,PPEZDIAG,PRSVS)
!!
!!*** *ARO_MNHC*  monitor of the chemical module of MesoNH-C
!!
!!    PURPOSE
!!    -------
!!       The purpose of this subroutine is to control the chemical module
!!    i.e. to pass the meteorological parameters from MesoNH to its chemical
!!    part and to call the different subroutines (calculation of rate constants,
!!    photolysis rates, stiff solver,..)
!!
!!    METHOD
!!    ------
!!       The calculation  of the chemical terms is performed using a loop
!!    over all spatial dimensions. 
!!
!!       For each single grid point, all necessary meteorological parameters are
!!    passed into the chemical core system (variable TZM). This variable is
!!    then passed on to the subroutines that calculate the reaction and
!!    photolysis rates. Then the chemical solver is called. As the chemistry
!!    part works with different units than MesoNH (MesoNH uses mixing ratio,
!!    the chemisty part uses molec/cm3) some unit conversion is also performed.
!!
!!       Temporal integration is performed over a double timestep 2*PTSTEP
!!    (except in the case of a cold start). If the timestep of MesoNH
!!    is too large for the chemical solver, several smaller steps can
!!    be taken using the NCH_SUBSTEPS parameter.
!!    "SPLIT"  : from PRSVS the scalar variable at t+dt is calculated and
!!               given as input to the solver; the result is rewritten 
!!               into PRSVS; this corresponds to applying first only dynamics
!!               and then only chemistry; this option assures positivity, but
!!               degrades the order of the temporal integration.
!!               In fact, an overhead of a factor two is produced here.
!!               A future solution will be to calculate the dynamics
!!               of the scalar variables not using leapfrog, but forward
!!               temporal integration.
!!     Vectorization Mask : Need ISVECMASK for photolysis (input(IN))
!!     ISVECNMASK = Nb of vector 
!!     IDT1       = Nb of points for the first dimension  (for arome = size of vector)
!!     IDT2       = Nb of points for the second dimension (for arome = 1)
!!     IDT3       = Nb of points in vertical
!!
!!
!!
!!    REFERENCE
!!    ---------
!!    Book 1, 2, 3 of MesoNH-chemistry
!!
!!    AUTHOR
!!    ------
!!    P. Tulet  *CNRM / GMEI* and contributors of MesoNH-C (K. Shure, C. Mari, V. Crassier)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original 10/11/04
!!
!!    EXTERNAL
!!    --------
USE PARKIND1  ,ONLY : JPIM, JPRB
!
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
IMPLICIT NONE
!
!*      0.1    declarations of arguments
INTEGER(KIND=JPIM),   INTENT(IN)   :: KLON     !NPROMA under CPG
INTEGER(KIND=JPIM),   INTENT(IN)   :: KLEV     !Number of vertical levels
INTEGER(KIND=JPIM),   INTENT(IN)   :: NSV     ! Number of passive scalar
INTEGER(KIND=JPIM),   INTENT(IN)   :: KRR      ! Number of moist variables
REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV,NSV),  INTENT(IN) :: PRSVSIN       ! Input source of scalar variable
REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),      INTENT(IN)    :: PRHODREF    ! iteration count
REAL(KIND=JPRB),    INTENT(IN)    :: PTSTEP      ! time step of MesoNH
REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),      INTENT(IN)    :: PTHT, PABST ! theta and pressure at t
REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV,KRR),  INTENT(IN)    :: PRT         ! moist variables at t
REAL(KIND=JPRB), DIMENSION(KLON,1),   INTENT(IN)    :: PLAT, PLON  ! lat, lon of each points
REAL(KIND=JPRB), DIMENSION(KLON,1),   INTENT(IN)    :: PALB_UV, PZS, PZENITH
REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),   INTENT(IN)    :: PZZ
INTEGER(KIND=JPIM), INTENT(IN)    :: KYEAR       ! Current Year
INTEGER(KIND=JPIM), INTENT(IN)    :: KMONTH      ! Current Month
INTEGER(KIND=JPIM), INTENT(IN)    :: KDAY        ! Current Day
REAL(KIND=JPRB),    INTENT(IN)    :: PTIME       ! Current time in second
INTEGER(KIND=JPIM), INTENT(IN)    :: KTCOUNT     ! iteration count
INTEGER(KIND=JPIM), INTENT(IN)    :: KLUOUT      ! unit for output listing count
INTEGER(KIND=JPIM),   INTENT(IN)   :: NDIAG   ! nb of diagnostics
REAL(KIND=JPRB), DIMENSION(KLON,KLEV,NDIAG),INTENT(INOUT)  :: PPEZDIAG
REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV,NSV),  INTENT(OUT) :: PRSVS       ! output source of scalar variable

!
!
END SUBROUTINE ARO_MNHC
END INTERFACE
