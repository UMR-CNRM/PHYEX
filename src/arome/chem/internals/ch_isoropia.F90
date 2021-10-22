!     ######spl
      SUBROUTINE CH_ISOROPIA (PAER,PRH, PDENAIR, PPRESSURE, PTEMP, PRC)
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
!!    P. Tulet ( Meteo France / GMEI)
!!
!!    MODIFICATIONS
!!    -------------
!!    Original
!!
!*****************************************************************

! ======================== ARGUMENTS / USAGE ===========================
!
!  INPUT:
!  1. [WI] 
!     DOUBLE PRECISION array of length [5].
!     Concentrations, expressed in moles/m3. Depending on the type of
!     problem solved (specified in CNTRL(1)), WI contains either 
!     GAS+AEROSOL or AEROSOL only concentratios.
!     WI(1) - sodium
!     WI(2) - sulfate
!     WI(3) - ammonium
!     WI(4) - nitrate
!     WI(5) - chloride
!
!  2. [RHI] 
!     DOUBLE PRECISION variable.  
!     Ambient relative humidity expressed on a (0,1) scale.
!
!  3. [TEMPI]
!     DOUBLE PRECISION variable. 
!     Ambient temperature expressed in Kelvins. 
!
!  4. [CNTRL]
!     DOUBLE PRECISION array of length [2].
!     Parameters that control the type of problem solved.
!
!     CNTRL(1): Defines the type of problem solved.
!     0 - Forward problem is solved. In this case, array WI contains 
!         GAS and AEROSOL concentrations together.
!     1 - Reverse problem is solved. In this case, array WI contains
!         AEROSOL concentrations only.
!
!     CNTRL(2): Defines the state of the aerosol
!     0 - The aerosol can have both solid+liquid phases (deliquescent)
!     1 - The aerosol is in only liquid state (metastable aerosol)
!
!  OUTPUT:
!  1. [WT] 
!     DOUBLE PRECISION array of length [5].
!     Total concentrations (GAS+AEROSOL) of species, expressed in moles/m3. 
!     If the foreward probelm is solved (CNTRL(1)=0), array WT is 
!     identical to array WI.
!     WT(1) - total sodium
!     WT(2) - total sulfate
!     WT(3) - total ammonium
!     WT(4) - total nitrate
!     WT(5) - total chloride
!
!  2. [GAS]
!     DOUBLE PRECISION array of length [03]. 
!     Gaseous species concentrations, expressed in moles/m3. 
!     GAS(1) - NH3
!     GAS(2) - HNO3
!     GAS(3) - HCl 
!
!  3. [AERLIQ]
!     DOUBLE PRECISION array of length [11]. 
!     Liquid aerosol species concentrations, expressed in moles/m3. 
!     AERLIQ(01) - H+(aq)          
!     AERLIQ(02) - Na+(aq)         
!     AERLIQ(03) - NH4+(aq)
!     AERLIQ(04) - Cl-(aq)         
!     AERLIQ(05) - SO4--(aq)       
!     AERLIQ(06) - HSO4-(aq)       
!     AERLIQ(07) - NO3-(aq)        
!     AERLIQ(08) - H2O             
!     AERLIQ(09) - NH3(aq) (undissociated)
!     AERLIQ(10) - HNCl(aq) (undissociated)
!     AERLIQ(11) - HNO3(aq) (undissociated)
!     AERLIQ(12) - OH-(aq)
!
!  4. [AERSLD]
!     DOUBLE PRECISION array of length [09]. 
!     Solid aerosol species concentrations, expressed in moles/m3. 
!     AERSLD(01) - NaNO3(s)
!     AERSLD(02) - NH4NO3(s)
!     AERSLD(03) - NaCl(s)         
!     AERSLD(04) - NH4Cl(s)
!     AERSLD(05) - Na2SO4(s)       
!     AERSLD(06) - (NH4)2SO4(s)
!     AERSLD(07) - NaHSO4(s)
!     AERSLD(08) - NH4HSO4(s)
!     AERSLD(09) - (NH4)4H(SO4)2(s)
!
!  5. [SCASI]
!     CHARACTER*15 variable.
!     Returns the subcase which the input corresponds to.
!
!  6. [OTHER]
!     DOUBLE PRECISION array of length [6].
!     Returns solution information.
!
!     OTHER(1): Shows if aerosol water exists.
!     0 - Aerosol is WET
!     1 - Aerosol is DRY
!
!     OTHER(2): Aerosol Sulfate ratio, defined as (in moles/m3) :
!               (total ammonia + total Na) / (total sulfate)
!
!     OTHER(3): Sulfate ratio based on aerosol properties that defines 
!               a sulfate poor system:
!               (aerosol ammonia + aerosol Na) / (aerosol sulfate)
!           
!     OTHER(4): Aerosol sodium ratio, defined as (in moles/m3) :
!               (total Na) / (total sulfate)
!      
!     OTHER(5): Ionic strength of the aqueous aerosol (if it exists).
!      
!     OTHER(6): Total number of calls to the activity coefficient 
!               calculation subroutine.
! 
! *** COPYRIGHT 1996-2000, UNIVERSITY OF MIAMI, CARNEGIE MELLON UNIVERSITY
! *** WRITTEN BY ATHANASIOS NENES
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
!!...........ARGUMENTS and their descriptions
!
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
INTEGER, PARAMETER :: NOTHER=6,NCOMP=5,NGASAQ=3,NSLDS=9,NIONS=7

INTEGER :: I
DOUBLE PRECISION,DIMENSION(2) :: CNTRL
DOUBLE PRECISION, DIMENSION(NOTHER) :: OTHER
CHARACTER(LEN=15) :: SCASI
DOUBLE PRECISION, DIMENSION(NCOMP) :: WI, WT
DOUBLE PRECISION, DIMENSION(NGASAQ) :: GAS
DOUBLE PRECISION, DIMENSION(NSLDS) :: AERSLD
DOUBLE PRECISION, DIMENSION(NIONS+NGASAQ+2) :: AERLIQ
DOUBLE PRECISION :: RHI,TEMPI


!-----------------------------------------------------------------------------

! Entry for isoropia

REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('CH_ISOROPIA',0,ZHOOK_HANDLE)
CNTRL(1) = 0
CNTRL(2) = 1


DO I=1,SIZE(PAER,1)

      WI(1)=0.
      WI(2)=PAER(I,1) * 1E-6 / ZMH2SO4
      WI(3)=PAER(I,2) * 1E-6 / ZMWNH3 + PAER(I,6) * 1E-6 / ZMWNH4
      WI(4)=PAER(I,5) * 1E-6 / ZMWNO3 + PAER(I,3) * 1E-6 / ZMWHNO3
      WI(5)=0.
      ! Temperature
      TEMPI=PTEMP(I)
      RHI=MIN(1.,MAX(0.,PRH(I)))

      CALL ISOROPIA (WI,RHI,TEMPI,CNTRL,WT,GAS,AERLIQ,AERSLD,SCASI,OTHER)


      PAER(I,1) = WT(2) * 1E6 * ZMH2SO4
      PAER(I,2) = GAS(1) * 1E6 * ZMWNH3
      PAER(I,3) = GAS(2) * 1E6 * ZMWHNO3
      PAER(I,4) = AERLIQ(8) * 1E6 * ZMWH2O
      PAER(I,5) = (WT(4) - GAS(2) ) * 1E6 * ZMWNO3 
      PAER(I,6) = (WT(3) - GAS(1) ) * 1E6 * ZMWNH4

ENDDO

IF (LHOOK) CALL DR_HOOK('CH_ISOROPIA',1,ZHOOK_HANDLE)
END SUBROUTINE CH_ISOROPIA
