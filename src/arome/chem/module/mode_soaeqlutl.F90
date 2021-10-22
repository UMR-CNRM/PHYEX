!     ######spl
module mode_soaeqlutl
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  !Purpose: module to iterate on the SOA concentration until equilibrium
  !Inspired by programs recieved from B. Pun/R. Griffin, 2005
  !Written by Alf Grini, alf.grini@cnrm.meteo.fr, 2005

  use modd_glo
  use modd_unifacparam
  use mode_unifac        !To get activity coefficients

  implicit none
  
contains
  
  !***********************************************
  subroutine XORG_GET(  &
       CPT              & !I [ug/m3] primary organic aerosol
       ,AERO_ORG        & !I [ug/m3] guessed SOA concentrations organic phase
       ,M0_ORG          & !O [ug/m3] total weight of organic aerosol (POA+SOA)
       ,X_ORG           & !O [frc] molar fraction of organics
       ,MW_AVG_ORG      & !O [g/mol] average molecular weight of organic phase
       )
    
    !Purpose: Calculate molar fraction of each component, both SOA and POA

    implicit none
    !Input
    REAL, DIMENSION(:,:), INTENT(IN)      :: CPT         ![ug/m3] concentration of POA
    REAL, DIMENSION(:,:), INTENT(IN)      :: AERO_ORG    ![ug/m3] concentration of SOA organic phase

    !Output
    REAL, DIMENSION(:), INTENT(OUT)       :: M0_ORG      ![ug/m3] total mass of organic aerosol
    REAL, DIMENSION(:,:), INTENT(OUT)     :: X_ORG       ![frc] molar fraction of organic SOA components
    REAL, DIMENSION(:), INTENT(OUT)       :: MW_AVG_ORG  ![g/mol] average molar weight of organic phase

    !Local
    INTEGER                               :: I            ![idx] counter for components
    INTEGER                               :: J            ![idx] counter for ionic species
    INTEGER                               :: COMP_IDX     ![idx] pointer to non-ion component
    REAL, DIMENSION(SIZE(CPT,1))          :: NTOT_ORG     ![umol/m3] total moles of organic phase
    
    !Initialize
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:XORG_GET',0,ZHOOK_HANDLE)
    M0_ORG(:)=0.d0
    NTOT_ORG(:)=0.d0
    MW_AVG_ORG(:)=0.d0

    !Add up total mass and moles of POA species
    DO I=1,NBSPOA
       M0_ORG(:) =  M0_ORG(:)+CPT(:,I)
       NTOT_ORG(:)= NTOT_ORG(:) + CPT(:,I)/MW_POA(I)
    ENDDO
    
    !Add up total mass and moles of SOA species
    COMP_IDX=1
    DO I=1,NBSP
       M0_ORG(:) =  M0_ORG(:) + AERO_ORG(:,I)                       !Total ug/m3
       NTOT_ORG(:)= NTOT_ORG(:) + AERO_ORG(:,I)/MW_SOA(COMP_IDX)    !Total umol/m3
       COMP_IDX=COMP_IDX+1
       !Skip ionic components
       DO J=2,NKAQ(I)
          COMP_IDX=COMP_IDX+1
       ENDDO
    ENDDO

    !Get molecular fraction and average molecular weight of organic phase
    !First calcualte for POA species
    DO I=1,NBSPOA
       X_ORG(:,I) = CPT(:,I)/MW_POA(I) & ![umol/m3]
            /NTOT_ORG(:)                 ![umol/m3]

       MW_AVG_ORG(:) = MW_AVG_ORG(:) + X_ORG(:,I)*MW_POA(I)
    ENDDO

    !Then calculate for SOA species
    COMP_IDX = 1
    DO I=1,NBSP

       !Get molar fraction of SOA component
       X_ORG(:,I+NBSPOA) = AERO_ORG(:,I) &
            /MW_SOA(COMP_IDX)/NTOT_ORG(:)

       !Add up the average molar weight of mixture
       MW_AVG_ORG(:)= MW_AVG_ORG(:)  &
            + X_ORG(:,I+NBSPOA)*MW_SOA(COMP_IDX)

       !Prepare for next compoenents
       COMP_IDX = COMP_IDX +1
       
       !Skip ionic components here
       DO J=2,NKAQ(I)
          COMP_IDX = COMP_IDX +1
       ENDDO !J
       
    ENDDO    !I

  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:XORG_GET',1,ZHOOK_HANDLE)
  end subroutine XORG_GET

  !*********************************************

  subroutine XAQ_GET(      &
       LWC                 & !I [ug/m3] liquid water content
       ,AERO_AQ            & !I [ug/m3] aquous phase aerosol concentrations of acids and ions
       ,X_AQ               & !O [frc] molar fraction of main component in aq phase
       ,TMOL               & !O [umole/m3] total aquous phase aerosol
       ,MW                 & !I [g/mol] molecular weight of species
       ,NK                 & !I [nbr] number of components for main acid
       ,NSPIN              &! I [nbr] number of main species
       )
    
    !Purpose: Get the molar fractions in the aquous phase
    implicit none
    
    !Input:
    REAL, DIMENSION(:), INTENT(IN)      :: LWC         !I [ug/m3] liquid water content
    REAL, DIMENSION(:,:), INTENT(IN)    :: AERO_AQ     !I [ug/m3] aquous phase aerosol concentrations of acids and ions
    INTEGER, DIMENSION(:), INTENT(IN)   :: NK          !I [nbr] number of components for main acid
    INTEGER, INTENT(IN)                 :: NSPIN       !I [nbr] number of main species (not including water!!)
    REAL, DIMENSION(:)                  :: MW          !I [g/mol] molecular weight of species

    !Output
    REAL, DIMENSION(:,:), INTENT(OUT)   :: X_AQ        !O [frc] molar fraction of main component in aquous phase
    REAL, DIMENSION(:), INTENT(OUT)     :: TMOL        !O [umol/m3] total moles of aquous phase

    !Local
    REAL, DIMENSION(SIZE(LWC))          :: TMOLINV     ![m3/umol] inverse total moles of aquous phase
    INTEGER                             :: I           ![idx] counter for aquous phase species
    INTEGER                             :: COMP_IDX    ![idx] pointer to correct aquous phase species
    INTEGER                             :: J           ![idx] counter for ions
    INTEGER                             :: NBR_SUB_AND_MAIN  ![nbr] number of sub and main species

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:XAQ_GET',0,ZHOOK_HANDLE)
    NBR_SUB_AND_MAIN = SUM(NK)

    ! Calculate mole fraction (molecules only)
    TMOL(:) = LWC(:) & ![ug/m3] 
         / MW_WATER ![g/mol] ==> umol/m3
    
    DO I=1,NBR_SUB_AND_MAIN  !All molecules and ions
       TMOL(:) = TMOL(:)   &
            + AERO_AQ(:,I)/MW(I)
    ENDDO
    TMOLINV(:)=1.d0/TMOL(:)
       
    !Get the molar fraction of water as the last components
    X_AQ(:,NSPIN+1) = LWC(:)/MW_WATER &
         *TMOLINV(:)
       
    !Only count X of main components 
    !(i.e. acids and their dissociating products)
    COMP_IDX=1
    DO I=1,NSPIN

       !Start by taking the amount of acid (non-dissociated component, H2A)
       X_AQ(:,I)=AERO_AQ(:,COMP_IDX)/MW(COMP_IDX)*TMOLINV(:)
       COMP_IDX=COMP_IDX + 1
          
       !Add up any concentration of dissociated species (HA- or A2-)
       DO J=2,NK(I)
          X_AQ(:,I)=X_AQ(:,I)+AERO_AQ(:,COMP_IDX)/MW(COMP_IDX)*TMOLINV(:)
          COMP_IDX=COMP_IDX + 1
       ENDDO
    ENDDO
    
  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:XAQ_GET',1,ZHOOK_HANDLE)
  end subroutine XAQ_GET
  
  !*****************************************
  
  !************************************************
  subroutine kpart_get(                    &
       TEMPK                               & !I [K] temperature
       ,GAMMA_ORG                          & !I [-] activity coefficient
       ,MW_AVG_ORG                         & !I [g/mol] average molar weight of organics
       ,KB                                 & !O [m3/ug] partitioning coefficient
       ,VP                                 & !I [torr] vapor pressure of SOA components
       ,NSKIP                              & !I [nbr] number of species not used in the gamma_org array
       )
    
    !Purpose: get the partitioning coefficient
    !Note that the activity coefficients here are for 18 species (POA + SOA)
    !Whereas we calculate the K only for 10 species (SOA)
    implicit none

    !INPUT
    REAL, DIMENSION(:), INTENT(IN)         :: TEMPK       !I [K] temperature
    REAL, DIMENSION(:,:), INTENT(IN)       :: GAMMA_ORG   !I [-] activity coefficients
    REAL, DIMENSION(:)                     :: MW_AVG_ORG  !I [g/mol] molecular weight organics
    REAL, DIMENSION(:,:)                   :: VP          !I [torr] saturation vapor pressure of component
    INTEGER, INTENT(IN)                    :: NSKIP       !I [nbr] number of species not used in the gamma_org array

    !OUTPUT
    REAL, DIMENSION(:,:), INTENT(OUT)      :: KB          !O [m3/ug] partitioning coefficient

    !LOCAL
    INTEGER                                :: I

    !The only goal of this routine is to get the partitioning coeffient:
    !So here is the way to do it:
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:KPART_GET',0,ZHOOK_HANDLE)
    do i=1,NBSP
       KB(:,i)= 760.d0         &!torr/atm
            *R_UNIV            &!(m3*atm)/(mol*K)
            *1.d-6             &!g/ug
            *TEMPK(:)          &!K
            /VP(:,I)           &!torr
            /MW_AVG_ORG(:)     &!g/mol
            /GAMMA_ORG(:,i+NSKIP)      !Note that we use I+NBSPOA here to avoid the POA species
       !==> m3/ug
    enddo

  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:KPART_GET',1,ZHOOK_HANDLE)
  end subroutine kpart_get

  !************************************************************

  subroutine gas_from_henry(    &
       AERO_AQ                  &  !I [ug/m3] aerosol concentration of acid [H2A] and ions
       ,GAMMA_HENRY             &  !I [-] Henry's law activity coeffients
       ,GAS_AQ                  &  !O [ug/m3] gas phase concentrations
       ,K                       &  !I [m3/ug and mole/kg] Henry's law and acid diss. constants
       ,LWC                     &  !I [ug/m3] Liquid water content
       ,NK                      &  !I [nbr] number of sub and main species for component
       ,NSPIN                   &  !I [nbr] number of main species to take into account
       )

    implicit none

    !INPUT
    REAL, DIMENSION(:,:), INTENT(IN)           :: AERO_AQ     !I [ug/m3] aerosol concentration of acid [H2A] and ions
    REAL, DIMENSION(:,:), INTENT(IN)           :: K           !I [m3/ug and mole/kg] Henry's law and acid diss. const
    REAL, DIMENSION(:), INTENT(IN)             :: LWC         !I [ug/m3] liquid water content
    REAL, DIMENSION(:,:), INTENT(IN)           :: GAMMA_HENRY !I [-] Henry's law std. state act. coeff
    INTEGER, DIMENSION(:), INTENT(IN)          :: NK          !I [nbr] number of sub and main species for component
    INTEGER, INTENT(IN)                        :: NSPIN       !I [nbr] number of main species to take into account

    !OUTPUT
    REAL, DIMENSION(:,:), INTENT(OUT)          :: GAS_AQ      !O [ug/m3] gas phase concentration in Henry's law eql with aerosol phase

    !LOCAL
    INTEGER                             :: I           ![idx] counter for aquous phase species
    INTEGER                             :: COMP_IDX    ![idx] pointer to correct aquous phase species
    INTEGER                             :: J           ![idx] counter for ions

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:GAS_FROM_HENRY',0,ZHOOK_HANDLE)
    COMP_IDX=1
    DO I=1,NSPIN

       !Apply Henry's law to get the Gas phase concentrations       
       GAS_AQ(:,I)=       &
            AERO_AQ(:,COMP_IDX)*GAMMA_HENRY(:,COMP_IDX)  &
            /(LWC(:)*K(:,COMP_IDX)) !Henry's law equilibrium

       !Prepare for next component
       COMP_IDX = COMP_IDX +1
       
       !Skip ionic components
       DO J=2,NK(I)
          COMP_IDX = COMP_IDX +1
       ENDDO

    ENDDO
       
  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:GAS_FROM_HENRY',1,ZHOOK_HANDLE)
  end subroutine gas_from_henry

  !***********************************************
  subroutine ACT_COEFF_HENRY_GET(   &
       X_AQ                         & !I [frc] molar fraction of components in aquous phase
       ,TMOL_AQ                     & !I [umol/m3] total number of umol in aquous phase
       ,GAMMA_RAOULT                & !I [-] Raoult's law standard state activity coefficents
       ,SI_AQ                       & !I [?] temperature dependent term for UNIFAC
       ,GAMMA_HENRY                 & !O [-] Henry's law standard state activity coefficients
       )

    !Purpose: transform Raoult's law standard state activity coefficients
    !into Henry' law standard state. This involves calculating activity coefficients
    !for component "i" in infinite dillute solution, and normalize the Raoult's
    !law activity coefficients
    implicit none
    !Input
    REAL, DIMENSION(:,:), INTENT(INOUT)      :: X_AQ           ![frc] fraction of components in aquous phase
    REAL, DIMENSION(:), INTENT(IN)           :: TMOL_AQ        ![umol/m3] total moles of aquous phase mixture
    REAL, DIMENSION(:,:), INTENT(IN)         :: GAMMA_RAOULT   ![-] raoult's law standard state activity coeffients
    REAL, DIMENSION(:,:,:), INTENT(IN)       :: SI_AQ          ![?] temperature dependent term for unifac

    !Output 
    REAL, DIMENSION(:,:), INTENT(OUT)        :: GAMMA_HENRY    ![-] Henry's law standard state activity coeffients

    !Local
    REAL, DIMENSION(SIZE(X_AQ,1))              :: TMPX      ![frc] saved value for mole fractions
    REAL, DIMENSION(SIZE(X_AQ,1),SIZE(X_AQ,2)) :: X2        ![frc] normalized mole fraction when one component is dilluted
    REAL, DIMENSION(SIZE(X_AQ,1),SIZE(X_AQ,2)) :: GAMMA     ![-] activity coeffients when one component is dilluted
    INTEGER                              :: I         ![idx] counter for main components (acids)
    INTEGER                              :: J         ![idx] counter for main components (acids)
    INTEGER                              :: COMP_IDX  ![idx] pointer to correct components
    

    ! calculate solute activity coefficient at infinite dilution 
    ! according to the standard state of Raoults Law 
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:ACT_COEFF_HENRY_GET',0,ZHOOK_HANDLE)
    COMP_IDX = 1
    DO I=1,NBSP
       
       !Save the real value of "X"
       TMPX(:) = X_AQ(:,I)  
       
       !Set component "I" to a value of inifinte dilution
       X_AQ(:,I) =  TMOL_AQ(:) * TINY
       
       !Make sure sum(X)=1 after modifying one value
       DO J=1,NMOL_AQ
          X2(:,J) = X_AQ(:,J)/SUM(X_AQ(:,1:NMOL_AQ))
          !write(6,*)"X2", X2(:,J)
       ENDDO
          
       ! gamma (Raoults law std st) at infinite dilution */
       CALL ACT_COEFF_GET (   &
            NU_AQ             &
            ,X2               &
            ,QG_AQ            &
            ,GAMMA            &
            ,THTAGP_AQ        &
            ,Q_AQ             &
            ,R_AQ             &
            ,L_AQ             &
            ,SI_AQ            &
            ,NMOL_AQ          &  !!note, NAMOL is not equal to NMOL_AQ since NMOL_AQ contains water
            ,NFUNC_AQ         &
            )
       
       !convert gamma for solutes to Henry's Law standard state
       !and set it equal for all main and sub-components "fxm: I am not sure where this equation comes from"
       DO J=1,NKAQ(I)
          GAMMA_HENRY(:,COMP_IDX)=         &  !Value for real "X" and Henry's law std. state
               GAMMA_RAOULT(:,COMP_IDX)    &  !Value for real "X" and Raoult's law std. state
               / GAMMA(:,I)                   !Value for dilute "X" and Raoult's law std. state
          
          COMP_IDX=COMP_IDX+1
       ENDDO
       !write(6,*)"act_coeff_henry", I, gamma_henry(:,COMP_IDX-1)
       
       !Set back X to the value before the "inifinite dilution" value
       X_AQ(:,I) = TMPX(:)
       
    ENDDO !Loop on components
    
  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:ACT_COEFF_HENRY_GET',1,ZHOOK_HANDLE)
  end subroutine ACT_COEFF_HENRY_GET

  !********************************************

  subroutine MB_adjust(        &
       CB                      &  !I [ug/m3] total (aerosol + gas) concentration
       ,GAS_AQ                 &  !I [ug/m3] gas concentration in eql with aquous phase
       ,GAS_ORG                &  !I [ug/m3] gas concentration in eql with org phase
       ,totAER_AQ              &  !I/O [ug/m3] total aerosol concentration (acid+ions) aq phase
       ,AERO_ORG               &  !I/O [ug/m3] total aerosol concentration in org phase
       ,GAS                    &  !O [ug/m3] new guess for gas phase concentration
       )

    !Purpose: Adjust mass balance after that we have calculated aerosol phase 
    !concentrations in equilibrium with each other

    !Author: Alf Grini, alf.grini@cnrm.meteo.fr, 2005
    
    implicit none

    !INPUT
    REAL, INTENT(IN), DIMENSION(:,:)          :: CB         ![ug/m3] total (aq + org + gas) SOA concentration
    REAL, INTENT(IN), DIMENSION(:,:)          :: GAS_AQ     ![ug/m3] gas phase concentration in eql with aq. phase
    REAL, INTENT(IN), DIMENSION(:,:)          :: GAS_ORG    ![ug/m3] gas phase concentration in eql with org. phase

    !INPUT/OUTPUT
    REAL, INTENT(INOUT), DIMENSION(:,:)       :: totAER_AQ  ![ug/m3] total concentration of H2A + HA- + HA(2-)
    REAL, INTENT(INOUT), DIMENSION(:,:)       :: AERO_ORG    ![ug/m3] organic aerosol of organic phase

    !OUTPUT
    REAL, INTENT(OUT), DIMENSION(:,:)         :: GAS        ![ug/m3] new guess for gas concentration

    !LOCAL
    REAL, DIMENSION(SIZE(CB,1),SIZE(CB,2))    :: FRAC_ORG   ![-] adjustment factor for organic aerosol phase 
    REAL, DIMENSION(SIZE(CB,1),SIZE(CB,2))    :: FRAC_AQ    ![-] adjustment factor for aquous aerosol phase
    REAL, DIMENSION(SIZE(CB,1),SIZE(CB,2))    :: XMB        ![ug/m3] sum of component in all phases 
    REAL, DIMENSION(SIZE(CB,1),SIZE(CB,2))    :: GAS_MEAN   ![ug/m3] mean of gas in eql with ORG and AQ phases
    REAL, DIMENSION(SIZE(CB,1),SIZE(CB,2))    :: FRAC_MB    ![-] correction factor for all phases to get mass balance
    INTEGER                                   :: I

    !Get mean gas phase concentration over organic phase and over aquous phase.
    !If one is higher than the other, we know that the aerosol concentration in the 
    !phase with highest gas is too high relative to aerosol conc in the other phase.
    !GAS_MEAN(:,:)=0.5d0*(GAS_AQ(:,:) + GAS_ORG(:,:))
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:MB_ADJUST',0,ZHOOK_HANDLE)
    GAS_MEAN(:,:)=sqrt(GAS_AQ(:,:)*GAS_ORG(:,:))

    IF(LPRINT)THEN
       DO I=1, NBSP
          write(6,*)"GAS_MEAN, GAS_AQ, GAS_ORG", I, GAS_MEAN(:,I), GAS_AQ(:,I), GAS_ORG(:,I)
       ENDDO
    ENDIF
    !Get deviation from mean for organics. 
    !If this is higher than 1 it means the vapor pressure over the organic phase is low
    !compared to the pressure over the aq. phase. (I.e. the aerosol concentration of org. phase is too low)
    FRAC_ORG(:,:)=GAS_MEAN(:,:)/GAS_ORG(:,:)

    !Get deviation from mean for aquous phase
    !If this is higher than 1 it means the vapor pressure over the aq phase is low
    !compared to the pressure over the org. phase. (I.e. the aerosol concentration of aq. phase is too low)
    FRAC_AQ(:,:)=GAS_MEAN(:,:)/GAS_AQ(:,:)

    IF(LPRINT)THEN
       DO I=1,NBSP
          write(6,*)"FRACS ORG/AQ ",I, FRAC_ORG(:,I), FRAC_AQ(:,I)
       ENDDO
    ENDIF

    !Adjust the aq aerosol phase using the information from the gas phase
    totAER_AQ(:,:)=totAER_AQ(:,:)*FRAC_AQ(:,:)

    !Adjust the org aerosol phase using the information from the gas phase
    AERO_ORG(:,:)=AERO_ORG(:,:)*FRAC_ORG(:,:)

    !Get the sum of all mass in the system. This equation will probably deviate from the
    !total mass of the system which is the input 
    XMB(:,:)=GAS_MEAN(:,:) + totAER_AQ(:,:) + AERO_ORG(:,:)
    IF(LPRINT)THEN
       DO I=1,NBSP
          write(6,*)"XMB/CB/Ratio",I, XMB(:,I), CB(:,I), XMB(:,I)/CB(:,I)
       ENDDO
    ENDIF

    !Get the adjustment factor for the mass balance
    FRAC_MB(:,:)=CB(:,:)/XMB(:,:)

    !The new guess for all is that we adjust everything linearly to the mass balance
    GAS(:,:)=GAS_MEAN(:,:)*FRAC_MB(:,:)

    !Total aerosol concentration of aquous phase
    totAER_AQ(:,:)=totAER_AQ(:,:)*FRAC_MB(:,:)

    !Total aerosol concentration of organic phase
    AERO_ORG(:,:)=AERO_ORG(:,:)*FRAC_MB(:,:)

  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:MB_ADJUST',1,ZHOOK_HANDLE)
  end subroutine MB_adjust

  !********************************************************************

  subroutine totAER_to_AERO(       &
       totAER_AQ                   & !I [ug/m3] total aerosol (acid+ions) in aq phase
       ,NOMINATOR                  & !I [-] term which gives acid if we know total (acid+ions)
       ,AERO_AQ                    & !O [ug/m3] acid and aerosol concentrations
       ,acHP                       & !I [mol/kg_{water}] proton concentration
       ,MW                         & !I [g/mol] molecular weight of species
       ,K                          & !I [m3/ug and mole/kg] Henry's law and acid diss. constants
       ,NK                         & !I [nbr] number of sub and main components for a species
       ,NSPIN                      & !I nbr] number of main species to take into account
       )

    !Purpose: Calculate the concentrations of ions and acids given 
    !that we know the gas concentrations of a component

    implicit none

    !INPUT
    REAL, DIMENSION(:,:), INTENT(IN)   :: totAER_AQ   !I [ug/m3] total (H2A+HA-+A(2-)) aq aerosol
    REAL, DIMENSION(:,:), INTENT(IN)   :: NOMINATOR   !I [-] factor which says how much of totAER_AQ which
                                                      !has dissociated
    REAL, DIMENSION(:,:), INTENT(IN)   :: K           !I [m3/ug and mole/kg] Henry and dissociation constants
    REAL, DIMENSION(:), INTENT(IN)     :: MW          !I [g/mol] molecular weight
    REAL, DIMENSION(:), INTENT(IN)     :: acHP        !I [mol/kg_{water}] proton concentration
    INTEGER, DIMENSION(:), INTENT(IN)  :: NK          !I [nbr] number of acid and ions for main component
    INTEGER, INTENT(IN)                :: NSPIN       !I [nbr] number of species to take into account

    !OUTPUT
    REAL, DIMENSION(:,:), INTENT(OUT)  :: AERO_AQ     ![ug/m3] all aq. phase aerosol sub components

    !LOCAL
    INTEGER                            :: COMP_IDX    ![idx] pointer to sub (ions) and main components
    INTEGER                            :: I           ![idx] counter for main components 
    INTEGER                            :: J           ![idx] counter for ions

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:TOTAER_TO_AERO',0,ZHOOK_HANDLE)
    COMP_IDX = 1
    DO I=1,NSPIN
    
       !Get aerosol concentration of first (main) component (H2A)
       !These have the same molar weight, so this is OK
       AERO_AQ(:,COMP_IDX) = totAER_AQ(:,I)/nominator(:,I)
       
       !Prepare for sub components
       COMP_IDX = COMP_IDX +1

       !Get aerosol concentrations of sub-components
       !Simply use equilibrium constants in aquous phase
       DO J=2,NK(I)
          !Note that we skip activity coefficients from this calculations since
          !they are anyway equal for all ions which are "children" of same main comp
          AERO_AQ(:,COMP_IDX) = K(:,COMP_IDX) * AERO_AQ(:,COMP_IDX-1)/acHP(:) &
               * MW(COMP_IDX)/MW(COMP_IDX-1)  !eqn 6/7 in Griffin 2003
          !Prepare for next sub-component
          COMP_IDX = COMP_IDX +1
       ENDDO

    ENDDO

  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:TOTAER_TO_AERO',1,ZHOOK_HANDLE)
  end subroutine totAER_to_AERO

  !**********************************************************************************
  
  subroutine AERO_to_totAER(         &
       AERO_AQ                      & !I [ug/m3] aquous phase concentrations (acid + ions)
       ,totAER_AQ                   & !O [ug/m3] aquous phase concentrations (total = sum of acid+ions)
       ,nominator                   & !I [-] transforming factor between acid and total
       ,NK                          & !I [nbr]  number of sub and mean species
       ,NSPIN                       & !I [nbr] number of main species to take into account
       )
    
    !Purpose: knowing the amount of acid, we can calculate the total amount of aerosol
    implicit none
    REAL, DIMENSION(:,:), INTENT(IN)       :: AERO_AQ   ![ug/m3] 
    REAL, DIMENSION(:,:), INTENT(IN)       :: NOMINATOR ![-] term to transfer between acid and acid+ions
    INTEGER, DIMENSION(:), INTENT(IN)      :: NK        ![nbr] number of sub and mean species
    INTEGER, INTENT(IN)                    :: NSPIN    ![nbr] number of main species to take into account

    REAL, DIMENSION(:,:), INTENT(OUT)      :: totAER_AQ ![ug/m3] 
    INTEGER                                :: I         ![idx] counter for main component
    INTEGER                                :: COMP_IDX  ![idx] pointer to right component including ions
    INTEGER                                :: J         ![idx] counter for ions

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:AERO_TO_TOTAER',0,ZHOOK_HANDLE)
    COMP_IDX=1
    DO I=1,NSPIN
       
       !Get aerosol concentration of first (main) component
       totAER_AQ(:,I) = AERO_AQ(:,COMP_IDX)*nominator(:,I)
       !Prepare for sub components
       COMP_IDX = COMP_IDX +1
       
       !Don't need to care about ion concentrations here because they don't influence
       !anything as long as we know the "nominator" for main component.
       !The "nominator" is only dependent on (H+) concentrations and on dissociation
       !coefficients. Just go ahead and skip the indexes concerning ions
       DO J=2,NK(I)
          COMP_IDX = COMP_IDX +1
       ENDDO

    ENDDO

  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:AERO_TO_TOTAER',1,ZHOOK_HANDLE)
  end subroutine AERO_to_totAER  

  !****************************************************************

  subroutine nominator_get(         &
       acHP                         & !I [mol/kg_{water}] proton concentration
       ,nominator                   & !O [-] term to transfer between acid and total aerosol
       ,K                           & !I [m3/ug and mole/kg] Henry's law and dissociation constants
       ,NK                          & !I [nbr] number of sub and main comp for acid
       ,NSPIN                       & !I [nbr] number of (main) species to take into account
       )

    !Purpose: 
    !In the aquous phase, we have the following equations: 
    !1) K(1) = H2A*gamma(1)/H2A(g)                      (Henry's law)
    !2) [H+]*[HA-]*gamma(2)/[H2A]/gamma(1) = K(2)       (acid dissociation)
    !3) [H+]*[A(2-)*gamma(3)]/[HA-]/gamma(2) = K(3)     (acid dissociation)

    !Remember that gamma(2) is equal to gamma(3) so we get back the original equations for (2) and (3)
    !1) K(1) = H2A*gamma(1)/H2A(g)    (Henry's law)
    !2) [H+]*[HA-]/[H2A] = K(2)       (acid dissociation)
    !3) [H+]*[A(2-)/[HA-] = K(3)      (acid dissociation)
    

    !The mole balance:
    !totAER = H2A + HA- + A(2-)
    !totAER = K(1)*[H2A(g)]/gamma(1) + K(2)[H2A]/[H+] + K(3)*[HA-]/[H+]
    !totAER = K(1)*[H2A(g)]/gamma(1) + K(2)K(1)[H2A(g)]/gamma(1)/[H+] + K(3)K(2)K(1)[H2A(g)]/gamma(1)//([H+])^2
    !totAER = K(1)*[H2A(g)]/gamma(1)* (1 + K(2)/[H+] + K(3)K(2)/([H+])^2)
    !totAER = [H2A]*(1 + K(2)/[H+] + K(3)K(2)/([H+])^2) = [H2A]/"nominator"

    !Remember that this is the mole balance, however, setting this up with mass, it turns
    !out that molecular weights disappears so "nominator is the same"
    !***********************************************************************************************
    implicit none

    REAL, DIMENSION(:), INTENT(IN)      :: acHP           ![mol/kg_{water}] proton concentration
    REAL, DIMENSION(:,:), INTENT(OUT)   :: NOMINATOR      ![-] factor to convert H2A to total aerosol 
    REAL, DIMENSION(:,:), INTENT(IN)    :: K              ![m3/ug and mole/kg] dissociation constants
    INTEGER, DIMENSION(:), INTENT(IN)   :: NK             ![nbr] number of sub and main comp for acid
    INTEGER, INTENT(IN)                 :: NSPIN          ![nbr] number of species to take into account

    REAL, DIMENSION(SIZE(NOMINATOR,1))  :: THISTERM       ![-] term in nominator
    INTEGER                             :: I              ![idx] counter for main component
    INTEGER                             :: J              ![idx] counter for sub components (ions)
    INTEGER                             :: COMP_IDX       ![idx] pointer to correct sub-component
    

    !Initialize
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:NOMINATOR_GET',0,ZHOOK_HANDLE)
    NOMINATOR(:,:)=1.d0
    COMP_IDX=1
    DO I=1, NSPIN

       !Initialize term in nominator calculations
       THISTERM(:)=1.d0
       
       !Prepare for next term
       COMP_IDX = COMP_IDX + 1
       
       !Start loop on ion concentrations
       DO J=2,NK(I)
          THISTERM(:)=THISTERM(:)*K(:,COMP_IDX)/acHP(:)
          NOMINATOR(:,I)=NOMINATOR(:,I) + THISTERM(:)
          COMP_IDX = COMP_IDX + 1
       ENDDO  !Loop on J
        
    ENDDO !Loop on I
    
  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQLUTL:NOMINATOR_GET',1,ZHOOK_HANDLE)
  end subroutine nominator_get
  
end module mode_soaeqlutl
