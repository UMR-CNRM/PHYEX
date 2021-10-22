!     ######spl
      SUBROUTINE CH_NNARES (PAER,PRH, PDENAIR, PPRESSURE, PTEMP, PRC)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!!    ##########################################################
!!
!!    PURPOSE
!!    -------
!!
!!    Calculate the aerosol chemical speciation and water content.
!!    THIS IS NEURAL NET VERSION OF THE ORIGINAL ARES CODE,
!!    GENERATED AND TRAINED USING NNFIT
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
!USE MODD_CH_M9
USE MODD_CH_AEROSOL
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

REAL, PARAMETER ::  ZMWH2O = 18.0           ! molecular weight for water      
REAL, PARAMETER ::  ZMWNO3 = 62.0049        ! molecular weight for NO3
REAL, PARAMETER ::  ZMWHNO3 = 63.01287      ! molecular weight for HNO3
REAL, PARAMETER ::  ZMH2SO4 = 98.07354      ! molecular weight for H2SO4
REAL, PARAMETER ::  ZMWNH3 = 17.03061       ! molecular weight for NH3
REAL, PARAMETER ::  ZMWNH4 = 18.03858       ! molecular weight for NH4

!-----------------------------------------------------------------------------
!     .. the independant variables of the NN are:
!     1) mole fraction of total ammonium           X = (NH4 + NH3) / (SO4)
!     2) mole fraction of total nitrate            Y = (NO3 + HNO3) / (SO4 + NO3 +  HNO3)
!     3) fractional relative humidity              R = RH
!     4) temperature                               T = TEMP
!
!     .. the output variables of the NN are:
!     1) mole fraction of water             A = AH2O / (AH20+ASO4+ANO3+ANH4)
!     2) mole fraction of gas phase ammonia B = GNH3 / (GNH3 + ANH4)
!     3) mole fraction of gas phase nitrate C = GNO3 / (GNO3 + ANO3)
!
REAL, DIMENSION(SIZE(PAER,1)) :: TOTAMM,TOTNIT,TOTSUL
!     .. input ..
REAL, DIMENSION(SIZE(PAER,1)) :: X, Y
!
!     .. output ..
REAL, DIMENSION(SIZE(PAER,1)) :: A, B, C
!
!-----------------------------------------------------------------------------
! print *, "PAER(1:5,1) H2SO4 ", PAER(1:5,1)
! print *, "PAER(1:5,2) NH3   ", PAER(1:5,2)
! print *, "PAER(1:5,3) HNO3  ", PAER(1:5,3)

!     convert to NN input variables X, Y, R, T
      REAL(KIND=JPRB) :: ZHOOK_HANDLE
      IF (LHOOK) CALL DR_HOOK('CH_NNARES',0,ZHOOK_HANDLE)
      TOTAMM(:) = PAER(:,2)/ZMWNH3 + PAER(:,6)/ZMWNH4
      TOTNIT(:) = PAER(:,3)/ZMWHNO3 + PAER(:,5)/ZMWNO3
      TOTSUL(:) = PAER(:,1)/ZMH2SO4

      X(:) = TOTAMM(:)/ TOTSUL(:)
      Y(:) = TOTNIT(:) / (TOTNIT(:) + TOTSUL(:))


      
!   call the neural net
      CALL NN ( X, Y, PRH, PTEMP, A, B, C, SIZE(PAER,1))

!     compute return variables in mole concentrations
      PAER(:,1) = TOTSUL(:) 
      PAER(:,2) = B(:) * TOTAMM(:) 
      PAER(:,3) = C(:) * TOTNIT(:) 
      PAER(:,5) = (1-C(:)) * TOTNIT(:) 
      PAER(:,6) = (1.-B(:)) * TOTAMM(:) 
      PAER(:,4) = (1./(1.-A(:))) * (PAER(:,1) + PAER(:,5) + PAER(:,6))

!     convert return variables to mass concentration
      PAER(:,1)=PAER(:,1)*ZMH2SO4
      PAER(:,2)=PAER(:,2)*ZMWNH3
      PAER(:,3)=PAER(:,3)*ZMWHNO3
      PAER(:,4)=PAER(:,4)*ZMWH2O
      PAER(:,5)=PAER(:,5)*ZMWNO3
      PAER(:,6)=PAER(:,6)*ZMWNH4

IF (LHOOK) CALL DR_HOOK('CH_NNARES',1,ZHOOK_HANDLE)
END SUBROUTINE CH_NNARES
