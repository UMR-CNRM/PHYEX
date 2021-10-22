INTERFACE
SUBROUTINE ARO_MNHDUST(KKL,KLON,KLEV,KSV,PTSTEP &
       ,PSVTIN      & !I [moments/molec_{air}] Transported moments of dust
       ,PZZ         & !I [m] height of layers
       ,PDZZ        & !I [m] layers thikness
       ,PPABST      & !I Pressure
       ,PTHT        & !I Potential temperature
       ,PRHODREF    & !I [kg/m3] density of air
       ,KSWB        & !I [nbr] number of shortwave bands
       ,KTCOUNT     & !I number of time step
       ,PSVT        & !O [moments/molec_{air}] Transported moments of dust
       ,PPIZA_WVL   & !O [-] single scattering albedo of dust
       ,PCGA_WVL    & !O [-] assymetry factor for dust layer
       ,PTAUREL_WVL & !O [-] opt.depth/opt.depth(550) for dust
       ,PAER        & !O [-] ext coeff at 550 for dust
       ,NDIAG       & !I [-] nb of diagnostics
       ,PPEZDIAG     & !IO [-] diag Nb/m3,ug/m3,rg(nb;um),rg(m;um),SSA,assym,AOD/550,mode & wvl
        )

!
!*** *ARO_MNHDUST*
!
!    PURPOSE
!    -------

!     Interface routine for initialisation of dust optical properties
!     before radiation scheme call  
  
!     AUTHOR.
!     -------
!      Y. Seity (CNRM/GMAP)
!      10-10-05

!    MODIFICATIONS
!    -------------
!    Original 10/10/05
!
!    EXTERNAL
!     -------

USE PARKIND1  ,ONLY : JPIM,JPRB


    IMPLICIT NONE
    
    !INPUT
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KKL      ! IKL under apl_arome
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KLON     ! NPROMA under CPG
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KLEV     ! Number of vertical levels
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KSV      ! Number of passive scalar
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KSWB     !I  number of shortwave wavelengths
    INTEGER(KIND=JPIM),   INTENT(IN)   :: KTCOUNT  !I  number of time step
    REAL(KIND=JPRB),      INTENT(IN)   :: PTSTEP   ! Time step
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV,KSV),INTENT(IN) :: PSVTIN !I [moments/molec{air}] transported moments of dust
    INTEGER(KIND=JPIM),   INTENT(IN)   :: NDIAG   ! nb of diagnostics
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)    :: PZZ        !I [m] height of layers
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)    :: PDZZ       !I [m] layers thikness
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)    :: PRHODREF   !I [kg/m3] density of air
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)    :: PTHT       !I [K] potentiel temperature
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV),INTENT(IN)    :: PPABST     !I [Pa] pressure

    !OUTPUT
    REAL(KIND=JPRB), DIMENSION(KLON,1,KLEV,KSV),INTENT(OUT) :: PSVT !O [moments/molec{air}] transported moments of dust
    REAL(KIND=JPRB), DIMENSION(KLON,KLEV,KSWB),INTENT(INOUT)   :: PPIZA_WVL   !O [-] SSA of dust layer for all SW wavelengths
    REAL(KIND=JPRB), DIMENSION(KLON,KLEV,KSWB),INTENT(INOUT)   :: PCGA_WVL    !O [-] assymetry for dust layer for all SW wvl
    REAL(KIND=JPRB), DIMENSION(KLON,KLEV,KSWB),INTENT(INOUT)   :: PTAUREL_WVL !O [-] AOD/550 for dust layer for all SW wvl 
    REAL(KIND=JPRB), DIMENSION(KLON,KLEV,NDIAG),INTENT(INOUT)  :: PPEZDIAG
    REAL(KIND=JPRB), DIMENSION(KLON,KLEV),      INTENT(INOUT)  :: PAER

END SUBROUTINE ARO_MNHDUST
END INTERFACE
