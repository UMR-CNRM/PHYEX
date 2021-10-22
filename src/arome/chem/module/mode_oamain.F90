!     ######spl
MODULE MODE_OAMAIN
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  !Purpose: main subroutine for calling SOA equilibrium routines
  !Similar programs recieved from B. Pun and R. Griffin, 2005
  !Rewritten by Alf Grini, alf.grini@cnrm.meteo.fr, 2005

  USE modd_glo

IMPLICIT NONE

CONTAINS

  SUBROUTINE MPMPO(          &
       TEMPK                   & !I [K] Temperature
       ,RH                     & !I [-] Relative humidity
       ,CB                     & !O [ug/m3] total (aerosol+gas) conc
       ,CPT                    & !I [ug/m3] primary aerosol concentration
       ,GAS                    & !O [ug/m3] gas phase concentrations
       ,AERO_ORG               & !O [ug/m3] organic phase concentrations
       ,AERO_AQ                & !I [ug/m3] aquous phase concentration of acid and ions
       ,PARTORG                & !O [ug/m3] aerosol phase concentration (aq+org)
       ,LWC                    & !I [ug/m3] liquid water content available for partitioning
       ,acHP                   & !I [mol/kg_{water}] H+ concentrations
       ,DELTALWC               & !O [ug/m3] change in LWC
       ,ORGANION               & !O [mole/m3] anion charge
       ,PSOLORG                & !IO [%] soluble SOA mass fraction
       )

    USE mode_soatinit           !Module which calculates parameters only dependent on temperature
    USE mode_firstguess         !Module which calculates first guesses of components
    USE mode_soaeql             !Module which does the equilibrium between aquous phase, organic phase and gas phase
    USE mode_zsrpun             !Module to get water associated with aquous organics
    USE modd_unifacparam       !Module with unifac coefficients
    USE modd_binsolu, only: molalbinAQ

    IMPLICIT NONE

    !INPUTS/OUTPUTS
    REAL, DIMENSION(:), INTENT(IN)       :: TEMPK      !I [K] Temperature
    REAL, DIMENSION(:), INTENT(IN)       :: RH         !I [0-1] Relative humidity
    REAL, DIMENSION(:,:), INTENT(IN)     :: CB         !I [ug/m3] total (g+p) aerosol organic species
    REAL, DIMENSION(:,:), INTENT(IN)     :: CPT        !I [ug/m3]
    REAL, DIMENSION(:,:), INTENT(OUT)    :: GAS        !I [ug/m3] gas of organic species
    REAL, DIMENSION(:,:), INTENT(OUT)    :: AERO_AQ    !I [ug/m3] aqueous species
    REAL, DIMENSION(:,:), INTENT(OUT)    :: AERO_ORG   !I [ug/m3] organic phase species
    REAL, DIMENSION(:,:), INTENT(OUT)    :: PARTORG    !I [ug/m3] aerosol phase organic species
    REAL, DIMENSION(:,:), INTENT(INOUT)  :: PSOLORG    !I soluble SOA mass fraction
    REAL, DIMENSION(:), INTENT(IN)       :: LWC        !I [ug/m3] liquid water content of aerosols
    REAL, DIMENSION(:), INTENT(IN)       :: acHP       !I [mol_{H+}/kg_{water}] concentration of H+ ions
    REAL, DIMENSION(:), INTENT(OUT)      :: DELTALWC   !I [ug/m3] LWC assiciated with organics
    REAL, DIMENSION(:), INTENT(OUT)      :: ORGANION   !I [mole/m3] organic anion concentratios

    !Local variables
    REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2)) :: VP  ![torr] temp. dependent vapor pressures of organic phase
    REAL, DIMENSION(SIZE(AERO_AQ,1),SIZE(AERO_AQ,2))   :: K   ![m3/ug and mole/kg] temp. dependent Henry's law and diss. const
    REAL, DIMENSION(:,:,:), ALLOCATABLE                :: SI_ORG  ![-] temp. dependent unifac coefficient
    REAL, DIMENSION(:,:,:), ALLOCATABLE                :: SI_AQ   ![-] temp. dependent unifac coefficient

    !Small variables
    INTEGER                              :: I         ![idx] counter for main type A components
    INTEGER                              :: J         ![idx] counter for sub type A components
    INTEGER                              :: COMP_IDX  ![idx] index for right type A comp (sub or main)
    INTEGER                              :: COMP_IDX2 ![idx] help counter for idx of main components
    REAL                                 :: XX        ![idx] account for difference in molecular weight

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_OAMAIN:MPMPO',0,ZHOOK_HANDLE)
    ALLOCATE(SI_ORG(SIZE(AERO_ORG,1),NFUNC_ORG,NFUNC_ORG))
    ALLOCATE(SI_AQ(SIZE(AERO_ORG,1),NFUNC_AQ,NFUNC_AQ))

    !Initialize the negative charge associated with type AQ organics (moles/m3)
    ORGANION(:)=0.d0

    !Initialize the liquid water content associated with AQ organics (ug/m3)
    deltaLWC(:)=0.d0

    !Initialize AERO_AQ
    AERO_AQ(:,:)=0.d0

    !Initialize AERO_ORG
    AERO_ORG(:,:)=0.d0

    !Initialize PARTORG / PSOLORG
    PARTORG(:,:)=0.d0
    PSOLORG(:,:)=0.d0

    !Initialize GAS
    GAS(:,:)=0.d0

    !Get termperature corrected vapor pressures 
    CALL VP_GET(           &
         TEMPK             & !I [K] Temperature
         ,VP               & !I [torr] saturation vapor pressures
         )

    !Get temperature corrected henry's law constants 
    CALL AQCONST_GET(            &
         TEMPK                   & !I [K] temperature
         ,K                      & !O [m3/ug and mole/kg] Henry's law and dissociation constants
         ,NKAQ                   & !I [nbr] Number of total and sub-components for one acid
         )

    !Get temperature dependent coefficients needed in UNIFAC parameterization
    CALL SI_GET(                 &
         TEMPK                    & !I [K] temperature
         ,A_ORG                   & !I [units?] term in UNIFAC parameterization
         ,SI_ORG                  & !O [units?] term in UNIFAC parameterization
         ,NFUNC_ORG               & !I [nbr] number of functional group in mixture
         )
    
    !Do the same thing, but now get the coefficients for aquous phase
    CALL SI_GET(                 &
         TEMPK                    & !I [K] temperature
         ,A_AQ                    & !I [units?] term in UNIFAC parameterization
         ,SI_AQ                   & !O [units?] term in UNIFAC parameterization
         ,NFUNC_AQ                & !I [nbr] number of functional group in mixture
         )

    !Get first guesses for aerosols
    CALL GUESS_AERO(  &
         CB             & !I [ug/m3] total (gas + aerosol + ions) of component
         ,AERO_ORG      & !O [ug/m3] liquid aerosol concentrations
         ,AERO_AQ       & !O [ug/m3] solid aerosol concentration
         ,GAS           & !O [ug/m3] gaseous concentration
         ,LWC           & !I [ug/m3] liquid water content already available for partitioning
         ,acHP          & !I [umol/kg_{water}] proton concentration
         ,VP            & !I [torr] saturation vapor pressure of organic precursors
         ,K             & !I [m3/ug and mole/kg] Henry's law and dissociation constants
         )

    !Use the guesses to iterate for the right solution
    CALL SOAEQL(             &
         acHP                & !I [mol/kg_{water}] proton concentration
         ,LWC                & !I [ug/m3] liquid water content
         ,CPT                & !I [ug/m3] primary organic concentration
         ,TEMPK              & !I [K] temperature
         ,GAS                & !I/O [ug/m3] gas phase concentrations of SOA
         ,AERO_AQ            & !I/O [ug/m3] aquous phase concentrations of SOA and ions
         ,AERO_ORG           & !I/O [ug/m3] organic phase concentrations of SOA
         ,CB                 & !I [ug/m3] total SOA concentration (GAS + AQ + ORG)
         ,K                  & !I [m3/ug and mole/kg] Henry's law and dissocation constants at right T
         ,SI_ORG            & !I [?] T-dependent coefficient in UNIFAC parameterization
         ,SI_AQ             & !I [?] T-dependent coefficient in UNIFAC parameterization
         ,VP                & !I [torr] saturation vapor pressure
         )

    !Get the wataer associated with the extra organics using ZSR
    CALL ZSRPUN(               &
         AERO_AQ               & !I [ug/m3] guess for aerosol concentrations
         ,RH                   & !I [0-1] relative humidity
         ,DELTALWC             & !O [ug/m3] liquid water content
         ,MOLALBINAQ           & !I [umol/ug_{water}] molality of binary solutions
         ,NBSP                 & !I [nbr] number of species to take into account
         ,NKAQ                 & !I [nbr] number of sub and main component for main component
         ,MW_SOA               & !I [g/mol] molecular weight of species
         )
    
    COMP_IDX=1
    DO I=1,NBSP

       !Start with summing the organic part
       PARTORG(:,I)=AERO_ORG(:,I)

       !Conserve COMP_IDX of main component
       COMP_IDX2 = COMP_IDX
       
       !Initiate partorg from this component
       PARTORG(:,I)=PARTORG(:,I) + AERO_AQ(:,COMP_IDX)
       PSOLORG(:,I)=  AERO_AQ(:,COMP_IDX)
       COMP_IDX=COMP_IDX+1
       
       !Add the contribution of the ions
       DO J=2,NKAQ(I)
          XX=MW_SOA(COMP_IDX2)  &  !MW of main component
               /MW_SOA(COMP_IDX)   !MW of this component
          
          !Sum PARTORG
          PARTORG(:,I) = PARTORG(:,I) + AERO_AQ(:,COMP_IDX)*XX
          PSOLORG(:,I) = PSOLORG(:,I) + AERO_AQ(:,COMP_IDX)*XX

          !Get anion concentration in mole/m3
          organion(:)=organion(:)  &
               +AERO_AQ(:,COMP_IDX)/MW_SOA(COMP_IDX) & !umole/m3
               *dble(J-1)                            & !number of negative charges
               *1.d-6                                  !==> mole/m3

          !Prepare for next component
          COMP_IDX=COMP_IDX+1
       ENDDO
       PSOLORG(:,I) = PSOLORG(:,I) / PARTORG(:,I)
    ENDDO

    DEALLOCATE(SI_ORG)
    DEALLOCATE(SI_AQ)

  IF (LHOOK) CALL DR_HOOK('MODE_OAMAIN:MPMPO',1,ZHOOK_HANDLE)
  END SUBROUTINE MPMPO

  !*************************************************************

  SUBROUTINE PUN(              &
       TEMPK                   & !I [K] Temperature
       ,RH                     & !I [-] Relative humidity
       ,WORG                   & !I [ug/m3] total (aerosol+gas) conc
       ,GASORG                 & !O [ug/m3] gas phase concentrations
       ,PAOM                   & !I [ug/m3] Primary aerosol organic matter
       ,PARTORG                & !O [ug/m3] aerosol phase concentrations
       ,LWC                    & !I [ug/m3] liquid water content
       ,acHP                   & !I [mol/kg_{water}] H+ concentrations
       ,DELTALWC               & !O [ug/m3] change in LWC
       ,ORGANION               & !O [mole//m3] anion charge
       ,TBOAFLAG               & !I [1/0] Flag for iteration or not
       )

    USE mode_bmain
    USE mode_amain

    IMPLICIT NONE

    !INPUTS/OUTPUTS
    REAL, DIMENSION(:), INTENT(IN)       :: TEMPK      ![K] Temperature
    REAL, DIMENSION(:), INTENT(IN)       :: RH         ![0-1] Relative humidity
    REAL, DIMENSION(:,:), INTENT(IN)     :: WORG       !I[ug/m3] total (g+p) aerosol organic species
    REAL, DIMENSION(:,:), INTENT(OUT)    :: GASORG     !I[ug/m3] gas of organic species
    REAL, DIMENSION(:), INTENT(IN)       :: PAOM       !I[ug/m3] Primary aerosol organic matter
    REAL, DIMENSION(:,:), INTENT(OUT)    :: PARTORG    !I[ug/m3] aerosol phase organic species
    REAL, DIMENSION(:), INTENT(IN)       :: LWC        !I[ug/m3] liquid water content of aerosols
    REAL, DIMENSION(:), INTENT(IN)       :: acHP       !I [mol_{H+}/kg_{water}] concentration of H+ ions
    REAL, DIMENSION(:), INTENT(OUT)      :: DELTALWC   !I [ug/m3] LWC assiciated with organics
    REAL, DIMENSION(:), INTENT(OUT)      :: ORGANION   !I [??] organic anion concentratios
    INTEGER, INTENT(INOUT)               :: TBOAFLAG   !I [flg] do type B calculations (1) or not (0)

    !LOCAL
    REAL, DIMENSION(SIZE(TEMPK),NBSPA)   :: AEROS     ![ug/m3] type A solid species below DRH
    REAL, DIMENSION(SIZE(TEMPK),NAAEROA) :: AERO      ![ug/m3] type A aqueous species
    REAL, DIMENSION(SIZE(TEMPK),NBSPA)   :: GASA      ![ug/m3] type A gas concentrations
    REAL, DIMENSION(SIZE(TEMPK),NBSPA)   :: totA      ![ug/m3] total (gas+aer) type A species
    REAL, DIMENSION(SIZE(TEMPK),NBSPB)   :: AEROB     ![ug/m3] type B concentration
    REAL, DIMENSION(SIZE(TEMPK),NBSPB)   :: GASB      ![ug/m3] type B gas species
    REAL, DIMENSION(SIZE(TEMPK),NBSPB)   :: CB        ![ug/m3] total concentration of type A

    !Small variables
    INTEGER                              :: I         ![idx] counter for main type A components
    INTEGER                              :: J         ![idx] counter for sub type A components
    INTEGER                              :: COMP_IDX  ![idx] index for right type A comp (sub or main)
    INTEGER                              :: COMP_IDX2 ![idx] help counter for idx of main components
    REAL                                 :: XX        ![idx] account for difference in molecular weight

    !***************************************************
    !FROM HERE WE ARE DEALING WITH TYPE A AEROSOL
    !**************************************************
    !Initialize the negative charge associated with type A organics (moles/m3)
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_OAMAIN:PUN',0,ZHOOK_HANDLE)
    ORGANION(:)=0.d0

    !Initialize the liquid water content associated with type A organics (ug/m3)
    deltaLWC(:)=0.d0

    !Initialize AERO (aquous type A species where RH > DRH)
    AERO(:,:)=0.d0
    !Initialize AEROS (solid type A species)
    AEROS(:,:)=0.d0

    !Put total aerosol A content in totA array
    totA(:,:)=worg(:,1:NBSPA)

    !Skipped a lot of checks here to see if aerosol concentrations are smaller or
    !larger than a lot of limits. Since the code is supposed to be called for a 
    !vector it is unlikely that the tests will apply for the whole vector.

    !Call equilibrium for type A aerosols
     CALL AMAIN(      &
          totA        & !I [ug/m3] total (gas+aerosol) organic content
          ,AERO       & !I/O [ug/m3] initial guess for Ai, also contain final output
          ,AEROS      & !I/O [ug/m3] solid aerosol concentration (==0 for the moment in this code)
          ,GASA       & !I/O [ug/m3] gas aerosol concentration
          ,LWC        & !I [ug/m3] input LWC (associated with inorganics), available for partitioning
          ,acHP       & !I [mole_{H+}/kg_{water}] H+ concentration
          ,DELTALWC   & !O [ug/m3] output LWC associated with aerosol
          ,ORGANION   & !O [mol/m3] moles negative charge 
          ,RH         & !I [0-1] relative humidity
          ,TEMPK      & !I [K] temperature
          )

     !Get the ouput values for PARTORG
     COMP_IDX=1
     DO I=1,NBSPA

        !Conserve COMP_IDX of main component
        COMP_IDX2 = COMP_IDX

        !Initiate partorg from this component
        PARTORG(:,I)=AERO(:,COMP_IDX)
        COMP_IDX=COMP_IDX+1

        DO J=2,NKA(I)
           XX=MWA(COMP_IDX2)  &  !MW of main component
                /MWA(COMP_IDX)   !MW of this component

           !Sum PARTORG
           PARTORG(:,I) = PARTORG(:,I) + AERO(:,COMP_IDX)*XX
           COMP_IDX=COMP_IDX+1
        ENDDO
     ENDDO

     !Get the output values for GASORG
     GASORG(:,1:NBSPA)=GASA(:,:)

    !************************END OF TYPE A AEROSOL*************

    !************************THE REST OF THE CALCULATIONS ARE FOR TYPE B AEROSOLS********
    IF(TBOAFLAG==1)THEN  !Check for doing type B aerosol
       
       !Get gas phase concentration of type B aerosol
       !GASB(:,:) = GASORG(:,NBSPA+1:NBSPA+NBSPB)

       !Get aerosol phase concentration of type B aerosol
       !AEROB(:,:) = PARTORG(:,NBSPA+1:NBSPA+NBSPB)

       CB(:,:)=WORG(:,NBSPA+1:NBSPA+NBSPB)

       CALL BMAIN (     &
            CB          &  !I [ug/m3] total concentration
            ,GASB       &  !I/O [ug/m3] concentration of type B gas
            ,AEROB      &  !I/O [ug/m3] concentration of type B aerosol
            ,PAOM       &  !I [ug/m3] Primary aerosol organic matter
            ,TEMPK      &  !I [K] temperature
            )

       !Set values back to the PARTORG table
       PARTORG(:,NBSPA+1:NBSPA+NBSPB) = AEROB(:,:)
       
       !Set values back to the GASORG table
       GASORG(:,NBSPA+1:NBSPA+NBSPB) = GASB(:,:)
       
       !Don't do type B calculations again if done once
       TBOAFLAG=0  

    ENDIF !TBOAFLAG is 1

  IF (LHOOK) CALL DR_HOOK('MODE_OAMAIN:PUN',1,ZHOOK_HANDLE)
  END SUBROUTINE PUN
    
END MODULE MODE_OAMAIN
