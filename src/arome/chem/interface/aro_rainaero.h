INTERFACE
SUBROUTINE ARO_RAINAERO(KLON,KLEV,NSV,KRR, PTSTEP &
       ,PSVTIN     & !I [moments/molec_{air}] Transported moments of dust
       ,PZZ        & !I [m] height of layers
       ,PPABST     & ! Pressure
       ,PTHT       & ! Potential temperature
       ,PRHODREF   & !I [kg/m3] density of air
       ,KTCOUNT    & ! number of time step
       ,PRT        & ! moist field
       ,PEVAP      & ! evaporation profile
       ,KSPLITR    & ! rain sedimentation time splitting
       ,PSVT       & !O [moments/molec_{air}] Transported moments of dust
       )

!
!*** *ARO_RAINAERO*
!
!    PURPOSE
!    -------

!     Interface routine for initialisation of dust optical properties
!     before radiation scheme call  
  
!     AUTHOR.
!     -------
!    P. Tulet

!    MODIFICATIONS
!    -------------
!    Original 10/10/05
!
!    EXTERNAL
!     -------



USE PARKIND1  ,ONLY : JPIM,JPRB


    IMPLICIT NONE
    
    !INPUT
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KLON     !NPROMA under CPG
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KLEV     !Number of vertical levels
    INTEGER(KIND=JPIM),   INTENT(IN)   :: NSV      ! Number of passive scalar
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KRR      ! Number of moist variables
    REAL(KIND=JPRB),      INTENT(IN)   :: PTSTEP   ! Time step

    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV,NSV),INTENT(IN)   :: PSVTIN     !I [moments/molec_{air}] transported moments of dust
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV,KRR),INTENT(IN)      :: PRT   ! Moist variables at time t
   
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PZZ        !I [m] height of layers
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PRHODREF   !I [kg/m3] density of air
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PTHT       !I [K] potentiel temperature
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PPABST     !I [Pa] pressure
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)        :: PEVAP      !I Evaporation
    INTEGER(KIND=JPIM), INTENT(IN)                      :: KTCOUNT    !I  number of time step
    INTEGER(KIND=JPIM), INTENT(IN)                      :: KSPLITR  ! Number of small time step
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV,NSV),INTENT(OUT)   :: PSVT       !O [moments/molec_{air}] transported moments of dust

END SUBROUTINE ARO_RAINAERO
END INTERFACE
