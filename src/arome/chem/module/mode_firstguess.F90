MODULE MODE_FIRSTGUESS
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  !Purpose: Make a first guess of the aerosol concentration in aquous phase and in organic phase
  !History: Original code recieved fro Robert Griffin, 2005
  !Rewritten by Alf Grini, alf.grini@cnrm.meteo.fr

  USE MODD_GLO, only : HLCRIT

  IMPLICIT NONE

CONTAINS

  SUBROUTINE GUESS_AERO(&
       CB             & !I [ug/m3] total (gas + aerosol + ions) of component
       ,AERO_ORG      & !O [ug/m3] liquid aerosol concentrations
       ,AERO_AQ       & !O [ug/m3] solid aerosol concentration
       ,GAS           & !O [ug/m3] gaseous concentration
       ,LWC           & !I [ug/m3] liquid water content already available for partitioning
       ,acHP          & !I [umol/kg_{water}] proton concentration
       ,VP            & !I [torr] saturation vapor pressure of organic precursors
       ,K             & !I [m3/ug and mole/kg] Henry's law and dissociation constants
       )

    !Purpose: Make initial guess for all concentrations when using Griffin's model
    !**********************************************************************************************
    !Make initial gas phase concentrations: For a set of acids we have the following equilibriums
    !1) H2A(g) <--> H2A
    !2) H2A <--> H+ + HA-
    !3) HA- <--> H+ + A(2-)
    
    !We have the following equations: 
    !1) K(1) = H2A/H2A(g)             (Henry's law)
    !2) [H+]*[HA-]/[H2A] = K(2)       (acid dissociation)
    !3) [H+]*[A(2-)]/[HA-] = K(3)     (acid dissociation)

    !The mole balance:
    !totAER = H2A + HA- + A(2-)
    !totAER = K(1)*[H2A(g)] + K(2)[H2A]/[H+] + K(3)*[HA-]/[H+]
    !totAER = K(1)*[H2A(g)] + K(2)K(1)[H2A(g)]/[H+] + K(3)K(2)K(1)[H2A(g)]/([H+])^2
    !totAER = K(1)*[H2A(g)]* (1 + K(2)/[H+] + K(3)K(2)/([H+])^2)
    !totAER = [H2A]*(1 + K(2)/[H+] + K(3)K(2)/([H+])^2) = [H2A]/"nominator"
    !***********************************************************************************************

    USE MODE_SOAEQLUTL

    IMPLICIT NONE
    
    !INPUT
    REAL, DIMENSION(:,:), INTENT(IN) :: CB      !I [ug/m3] total (gas+aerosol) concentration
    REAL, DIMENSION(:), INTENT(IN)   :: LWC     !I [ug/m3] Liquid water content already ready for partitioning
    REAL, DIMENSION(:), INTENT(IN)   :: acHP    !I [umol/kg_{water}] proton concentration
    REAL, DIMENSION(:,:), INTENT(IN) :: VP      !I [torr] temperature corrected vapor pressure of gas precursors
    REAL, DIMENSION(:,:), INTENT(IN) :: K       !I [m3/ug and mole/kg] Henry's law and dissociation constants

    !OUTPUT
    REAL, DIMENSION(:,:), INTENT(OUT)  :: GAS         !O [ug/m3] gas phase concentrations
    REAL, DIMENSION(:,:), INTENT(OUT)  :: AERO_ORG    !O [ug/m3] liquid aerosol concentrations
    REAL, DIMENSION(:,:), INTENT(OUT)  :: AERO_AQ     !O [ug/m3] solid aerosol concentrations

    !LOCAL, SMALL COUNTERS AND OTHER VARS
    REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2)) :: NOMINATOR ![-] factor to transfer between totAER and acid/ions in aq
    REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2)) :: totAER_AQ ![ug/m3] total SOA in aq. phase (acid+ions)
    INTEGER                      :: I            ![idx] counter for main components (1-NAMOL)
    INTEGER                      :: J            ![idx] counter for sub-componnets or ions
    INTEGER                      :: COMP_IDX     ![idx] index for aerosol components (1-NAAERO)


    !******************************************************************************

    !Start part which will make guesses for the 
    !aerosol concentrations before the iterations start
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_FIRSTGUESS:GUESS_AERO',0,ZHOOK_HANDLE)
    COMP_IDX=1
    DO I=1,NBSP
       
       !If saturation vapor pressure is low, then guess 30% in organic aerosol phase
       !If not, guess 15% in organic aerosol phase
       WHERE(VP(:,I).lt.VPCRIT)
          AERO_ORG(:,I)=0.3d0*CB(:,I)
       ELSEWHERE
          AERO_ORG(:,I)=0.15*CB(:,I)
       ENDWHERE
       
       !If high LWC and high Henry's law constant guess 10% in aerosol phase
       !Else, guess 5% in aerosol phase
       WHERE(LWC(:).gt.LWCCRIT.AND.K(:,COMP_IDX).gt.HLCRIT)
          AERO_AQ(:,COMP_IDX)=0.1*CB(:,I)
       ELSEWHERE
          AERO_AQ(:,COMP_IDX)=0.05*CB(:,I)
       ENDWHERE
       
       !Prepare for next component
       COMP_IDX=COMP_IDX+1

       !Skip ions for correct indexing of K
       DO J=2,NKAQ(I)
          COMP_IDX=COMP_IDX+1
       ENDDO

    ENDDO !I

    !Guess for gas phase has to be whatever is not in aquous or gaseous phase
    !GAS(:,:)=CB(:,:) - AERO_ORG(:,:) - totAER_AQ(:,:)
    !The guess for gas is not really used anywhere since gas is calculated from
    !the aerosol guesses during the iterations.. See figure explaining iteration
    !cycle in the theory paper. Set it to zero anyway just to have a number...
    !fxm: in a future version GAS should not be output from this routine!
    GAS(:,:)=0.d0

    !Get nominator to transfer total aquous to acid/ions
    CALL NOMINATOR_GET(             &
       acHP                         & !I [mol/kg_{water}] proton concentration
       ,nominator                   & !O [-] term to transfer between acid and total aerosol
       ,K                           & !I [m3/ug and mole/kg] Henry's law and dissociation constants
       ,NKAQ                        & !I [nbr] number of main and sub components
       ,NBSP                        & !I [nbr] total number of components
       )

    !Transfer the H2A concentration to total aerosol concentration (acid and ions)
    CALL AERO_TO_totAER(              &
         AERO_AQ                      & !I [ug/m3] aquous phase concentrations (acid + ions)
         ,totAER_AQ                   & !O [ug/m3] aquous phase concentrations (total = sum of acid+ions)
         ,nominator                   & !I [-] transforming factor between acid and total
         ,NKAQ                        & !I [nbr] number of main and sub components
         ,NBSP                        & !I [nbr] total number of components
         )

    !transfer total aerosol to aquous/ions
    CALL totAER_TO_AERO(            &
       totAER_AQ                   & !I [ug/m3] total aerosol (acid+ions) in aq phase
       ,NOMINATOR                  & !I [-] term which gives acid if we know total (acid+ions)
       ,AERO_AQ                    & !O [ug/m3] acid and aerosol concentrations
       ,acHP                       & !I [mole/kg_{water}] proton concentration
       ,MW_SOA                     & !I [g/mol] molecular weight of species
       ,K                          & !I [m3/ug and mole/kg] Henry's law and acid diss. constants
       ,NKAQ                       & !I [nbr] number of sub and main components for a species
       ,NBSP                       & !I nbr] number of main species to take into account
       )

  IF (LHOOK) CALL DR_HOOK('MODE_FIRSTGUESS:GUESS_AERO',1,ZHOOK_HANDLE)
  END SUBROUTINE GUESS_AERO
  !************************************************

END MODULE mode_firstguess
  
     
