!     ######spl
      SUBROUTINE CH_ARES (PAER,PRH, PDENAIR, PPRESSURE, PTEMP, PRC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ##########################################################
!!
!!    PURPOSE
!!    -------
!!
!!    Calculate the aerosol chemical speciation and water content.
!!
!!    REFERENCE
!!    ---------
!!    none
!!
!!    AUTHOR
!!    ------
!!    V. Crassier & K. Suhre (LA)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!!
!*****************************************************************
! Parameters of ARES:
!
!  SO4   : Total sulfate in MICROGRAMS/M**3 as sulfate (IN)
!  HNO3  : Nitric Acid in MICROGRAMS/M**3 as nitric acid (IN)
!  NO3   : Total nitrate in MICROGRAMS/M**3 as nitric acid (IN)
!  NH3   : Total ammonia in MICROGRAMS/M**3 as ammonia (IN)
!  NH4   : Ammonium in MICROGRAMS/M**3 as ammonium (IN)
!  RH    : Fractional relative humidity (IN)
!  TEMP  : Temperature in Kelvin (IN)
!  GNO3  : Gas phase nitric acid in MICROGRAMS/M**3 (OUT)
!  GNH3  : Gas phase ammonia in MICROGRAMS/M**3 (OUT)
!  ASO4  : Aerosol phase sulfate in MICROGRAMS/M**3 (OUT)
!  ANO3  : Aerosol phase nitrate in MICROGRAMS/M**3 (OUT)
!  ANH4  : Aerosol phase ammonium in MICROGRAMS/M**3 (OUT)
!  AH2O  : Aerosol phase water in MICROGRAMS/M**3 (OUT)
!
!***************************************************************
!!
!!   EXTERNAL
!!   -------
!!
!!
!!   IMPLICIT ARGUMENTS
!!   ------------------
!!
IMPLICIT NONE
!...........ARGUMENTS and their descriptions
REAL, DIMENSION(:,:), INTENT(INOUT) :: PAER
REAL, DIMENSION(:),   INTENT(IN)    :: PRH, PDENAIR, PPRESSURE, PTEMP, PRC
!
! PAER(:,1) :: H2SO4    in micrograms / m**3
! PAER(:,2) :: NH3(g)   in micrograms / m**3
! PAER(:,3) :: HNO3(g)  in micrograms / m**3
! PAER(:,4) :: H2O(a)   in micrograms / m**3
! PAER(:,5) :: NO3(a)   in micrograms / m**3
! PAER(:,6) :: NH4(a)   in micrograms / m**3
!
!...........PARAMETERS and their descriptions:

REAL, PARAMETER ::  ZMWNO3 = 62.0049        ! molecular weight for NO3
REAL, PARAMETER ::  ZMWHNO3 = 63.01287      ! molecular weight for HNO3
REAL, PARAMETER ::  ZMWNH3 = 17.03061       ! molecular weight for NH3
REAL, PARAMETER ::  ZMWNH4 = 18.03858       ! molecular weight for NH4

INTEGER :: I

REAL :: SO4,HNO3,NH3,NO3,NH4,TEMP,RH
REAL :: ASO4,ANO3,AH2O,ANH4,GNH3,GNO3

EXTERNAL RPMARES
!-----------------------------------------------------------------------------

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_ARES',0,ZHOOK_HANDLE)
DO I=1,SIZE(PAER,1)

      SO4=PAER(I,1)
      HNO3=PAER(I,3)+PAER(I,5)*ZMWHNO3/ZMWNO3
      NH3=PAER(I,2)+PAER(I,6)*ZMWNH3/ZMWNH4
      NO3=0.
      NH4=0.
      TEMP= PTEMP(I)
      RH  = PRH(I)
       
!   call ARES
      CALL RPMARES ( SO4, HNO3, NO3, NH3, NH4, RH, TEMP,      &
                            ASO4, ANO3, AH2O, ANH4, GNH3, GNO3)

      PAER(I,1)=ASO4
      PAER(I,2)=GNH3
      PAER(I,3)=GNO3
      PAER(I,4)=AH2O
      PAER(I,5)=ANO3
      PAER(I,6)=ANH4

ENDDO

IF (LHOOK) CALL DR_HOOK('CH_ARES',1,ZHOOK_HANDLE)
END SUBROUTINE CH_ARES
