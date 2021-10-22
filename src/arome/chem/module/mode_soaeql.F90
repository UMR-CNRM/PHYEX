!     ######spl
module mode_soaeql
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  !Purpose: module to iterate on the SOA concentration until equilibrium
  !Inspired by programs recieved from B. Pun/R. Griffin, 2005
  !Written by Alf Grini, alf.grini@cnrm.meteo.fr, 2005

  use modd_glo
  use modd_unifacparam
  use mode_soaeqlutl
  use mode_unifac
  
  IMPLICIT NONE
  
CONTAINS
  
  subroutine soaeql(       &
       acHP                & !I [mol/kg_{water}] proton concentration
       ,LWC                & !I [ug/m3] liquid water content
       ,CPT                & !I [ug/m3] primary organic concentration
       ,TEMPK              & !I [K] temperature
       ,GAS                & !I/O [ug/m3] gas phase concentrations of SOA
       ,AERO_AQ            & !I/O [ug/m3] aquous phase concentrations of SOA and ions
       ,AERO_ORG           & !I/O [ug/m3] organic phase concentrations of SOA
       ,CB                 & !I [ug/m3] total SOA concentration (GAS + AQ + ORG)
       ,K                  & !I [m3/ug and mole/kg] Henry's law and dissocation constants at right T
       ,SI_ORG             & !I [?] T-dependent coefficient in UNIFAC parameterization
       ,SI_AQ              & !I [?] T-dependent coefficient in UNIFAC parameterization
       ,VP                 & !I [torr] saturation vapor pressures for organic phase SOA
       )

    !Purpose: Iterate the following equations until convergence
    !1) Equilibrium between organic phase and gas phase
    !2) Equilibrium between aquous phase and gas phase
    !3) Dissociation in aquous phase (H2A = H+ + HA- + A--)
    !4) Mass balance for all species
    IMPLICIT NONE

    !Input
    REAL, INTENT(IN), DIMENSION(:,:)        :: CPT               !I [ug/m3] Primary organic aerosols
    REAL, INTENT(IN), DIMENSION(:)          :: TEMPK             !I [K] Temperature
    REAL, INTENT(IN), DIMENSION(:)          :: acHP              !I [mol/kg_{water}] proton concentration
    REAL, INTENT(IN), DIMENSION(:)          :: LWC               !I [ug/m3] liquid water content 
    REAL, INTENT(IN), DIMENSION(:,:,:)      :: SI_AQ             !I [-] term in UNIFAC for aquous component
    REAL, INTENT(IN), DIMENSION(:,:,:)      :: SI_ORG            !I [-] term in UNIFAC for organic phase
    REAL, INTENT(IN), DIMENSION(:,:)        :: VP                !I [torr] saturation vapor pressures for SOA phase
    REAL, INTENT(IN), DIMENSION(:,:)        :: CB                !I [ug/m3] total concentration of (AQ+ORG+GAS)
    REAL, INTENT(IN), DIMENSION(:,:)        :: K                 !I [m3/ug and mole/kg] Henry's law and dissociation constants

    !Input/output
    REAL, INTENT(INOUT), DIMENSION(:,:)     :: AERO_ORG          !I/O [ug/m3] SOA in organic phase
    REAL, INTENT(INOUT), DIMENSION(:,:)     :: AERO_AQ           !I/O [ug/m3] SOA in aquous phase (acids and ions)
    REAL, INTENT(INOUT), DIMENSION(:,:)     :: GAS               !I/O [ug/m3] SOA in gas phase

    !Local
    REAL, DIMENSION(SIZE(TEMPK))                        :: M0_ORG            ![ug/m3] total organic aerosol in grid cell
    REAL, DIMENSION(SIZE(TEMPK))                        :: MW_AVG_ORG        ![g/mol] average molecular weight of organic phase
    REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2))  :: NOMINATOR         ![-] term to transfer between H2A and total AQ aerosol
    REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2))  :: GAS_ORG           ![ug/m3] gas concentration in eql with organic phase
    REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2))  :: GAS_AQ            ![ug/m3] gas concentration in eql with aquous phase
    REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2))  :: KPART             ![m3/ug] partitioning coeffient for organic phase
    REAL, DIMENSION(SIZE(AERO_ORG,1))                   :: TMOL_AQ           ![umol/m3] total moles of aquous phase
    REAL, DIMENSION(SIZE(AERO_ORG,1),SIZE(AERO_ORG,2))  :: totAER_AQ         ![ug/m3] total (acid +ions) concentration of aq. phase
    REAL, DIMENSION(:,:), ALLOCATABLE       :: X_ORG                         ![frc] fraction of component in organic phase POA + SOA
    REAL, DIMENSION(:,:), ALLOCATABLE       :: X_AQ                          ![frc] fraction of main components in aquous phase (SOA + H2O)
    REAL, DIMENSION(:,:), ALLOCATABLE       :: GAMMA_ORG                     ![-] activity coefficients for organic phase (POA and SOA)
    REAL, DIMENSION(:,:), ALLOCATABLE       :: GAMMA_AQ                      ![-] activity coefficients for aquous phase (Raoult and main components)
    REAL, DIMENSION(:,:), ALLOCATABLE       :: GAMMA_AQ_HENRY                ![-] activity coefficients for aquous phase (Henry's law std. state)
    REAL, DIMENSION(:,:), ALLOCATABLE       :: GAMMA_AQ_RAOULT               ![-] activity coefficients for aquous phase (Raoult's law std. state)
    INTEGER                                 :: I                             ![idx] counter for main components (acids)
    INTEGER                                 :: J                             ![idx] counter for sub components (ions)
    INTEGER                                 :: COMP_IDX                      ![idx] pointer to right component
    INTEGER                                 :: ITER_IDX                      ![idx] counter for number of iterations
    INTEGER, PARAMETER                      :: MAX_ITER=6                    !Number of iterations performed

    !Allocate memory
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_SOAEQL:SOAEQL',0,ZHOOK_HANDLE)
    ALLOCATE(X_ORG(SIZE(AERO_ORG,1),NBSP+NBSPOA))
    ALLOCATE(X_AQ(SIZE(AERO_ORG,1),NBSP+1))
    ALLOCATE(GAMMA_ORG(SIZE(AERO_ORG,1),NBSP+NBSPOA))
    ALLOCATE(GAMMA_AQ(SIZE(AERO_ORG,1),NBSP+1))
    ALLOCATE(GAMMA_AQ_RAOULT(SIZE(AERO_ORG,1),NAAERO))
    ALLOCATE(GAMMA_AQ_HENRY(SIZE(AERO_ORG,1),NAAERO))


    !Before starting iterations, get "nominator" which is term to transfer between totAER_AQ and AERO_AQ (acid)
    CALL NOMINATOR_GET(               &
         acHP                         & !I [mol/kg_{water}] proton concentration
         ,nominator                   & !O [-] term to transfer between acid and total aerosol
         ,K                           & !I [m3/ug and mole/kg] Henry's law and acid dissociation constants
         ,NKAQ                        & !I [nbr] number of sub and main comp for acid
         ,NBSP                        & !I [nbr] number of (main) species to take into account
         )

    DO ITER_IDX = 1,MAX_ITER
       IF(LBOX)THEN
          write(6,*)"*******************************************"
          write(6,*)"**************NEW ITERATION #", ITER_IDX,"********"
          write(6,*)"********************************************"
       ENDIF

       !Get molar fraction of species in organic phase
       CALL XORG_GET(     &
            CPT           & !I [ug/m3] concentrations of primary organic aerosols
            ,AERO_ORG     & !I [ug/m3] guessed concentration of aerosols 
            ,M0_ORG       & !O [ug/m3] total organic mass (POA + SOA)
            ,X_ORG        & !O [frc] fractions in organic phase of POA and SOA
            ,MW_AVG_ORG   & !I [g/mol] average molecular weight of organic phase            
            )
       IF(LPRINT)THEN
          DO I=1,NBSPOA
             write(6,*)"CPT", CPT(:,I)
          ENDDO
          DO I=1,NBSP
             write(6,*)"AERO/CB", AERO_ORG(:,I), CB(:,I)
          ENDDO
       ENDIF
       IF(LBOX)THEN
          do I =1,NBSP+NBSPOA
             write(6,*)"X_ORG",I, X_ORG(:,I)
          ENDDO
       ENDIF

       !Get activity coefficients
       CALL ACT_COEFF_GET(  &
            NU_ORG          &
            ,X_ORG          &
            ,QG_ORG         &
            ,GAMMA_ORG      &
            ,THTAGP_ORG     &
            ,Q_ORG          &
            ,R_ORG          &
            ,L_ORG          &
            ,SI_ORG         &
            ,NMOL_ORG       & !POA and SOA species (normally 18 species)
            ,NFUNC_ORG      &
            )

       !Get partitioning coefficients
       CALL KPART_GET(          &
            TEMPK               & !I [K] temperature
            , GAMMA_ORG         & !I [-] Activity coefficient 
            , MW_AVG_ORG        & !I [g/mol] Molecular weight of organic matter
            , KPART             & !O [m3/ug] partitioning coefficient for species
            , VP                & !I [torr] vapor pressures
            , NBSPOA            & !I [nbr] number of activity coefficients to skip in GAMMA_ORG
            )

       IF(LBOX)THEN
          do I=1,NBSP
             write(6,*)"KOM/GAMMA_ORG",I, KPART(:,I),GAMMA_ORG(:,I)
          enddo
       ENDIF

       !Calculate gas phase concentrations from organic aerosol phase
       DO I=1,NBSP
          GAS_ORG(:,I) = AERO_ORG(:,I) &        ![ug/m3]
               /(KPART(:,I)*M0_ORG(:))          ![no unit] ==> [ug/m3]
       ENDDO

       !At this point, we have aerosol and gas phase concentrations which probably
       !are not in accordance with mass balance, we check for this lower down in the
       !subroutine, now take care of aquous phase

       !Get molar fraction of component in aquous phase
       !Only main components, not ions! Ions are summed into main component!
       CALL XAQ_GET(            &
            LWC                 & !I [ug/m3] liquid water content
            ,AERO_AQ            & !I [ug/m3] concentrations of acids and ions in aquous phase
            ,X_AQ               & !O [frc] molar fraction of main component (sum of acid+ions) in aq phase
            ,TMOL_AQ            & !O [umol/m3] total mole of aquous phase
            ,MW_SOA             & !I [g/mol] molecular weight of species
            ,NKAQ               & !I [nbr] number of components for main acid
            ,NBSP               & !I [nbr] number of main species
            )

       IF(LBOX)THEN
          DO I=1,NBSP+1
             write(6,*)"X_AQ" , I, X_AQ(:,I)
          ENDDO
       ENDIF

       IF(LPRINT)THEN
          COMP_IDX=1
          DO I=1,NBSP
             write(6,*)"AERO_AQ",I, AERO_AQ(:,COMP_IDX)
             COMP_IDX=COMP_IDX+1
             DO J=2,NKAQ(I)
                COMP_IDX=COMP_IDX+1
             ENDDO
          ENDDO
          write(6,*)"LWC", LWC(:)
          write(6,*)"HPLUS", acHP(:)
       ENDIF

       !Get raoult's law standard state activity coefficients
       !Here, only treat main components, i.e. (H2A, HA- and A(2-) as one comp)
       CALL ACT_COEFF_GET( &
            NU_AQ          &
            ,X_AQ          &
            ,QG_AQ         &
            ,GAMMA_AQ      &
            ,THTAGP_AQ     &
            ,Q_AQ          &
            ,R_AQ          &
            ,L_AQ          &
            ,SI_AQ         &
            ,NMOL_AQ       &  !Use the 10 main components and water (=11 components)
            ,NFUNC_AQ      &
            )

       IF(LPRINT)THEN
          do I=1,11
             write(6,*)"GAMMA_AQ",I,gamma_aq(:,I)
          enddo
       ENDIF

       !Set activity coeffients for all components, ions and acids
       !to the Raoult's law standard state value for the real X
       COMP_IDX=1
       DO I=1,NBSP
          GAMMA_AQ_RAOULT(:,COMP_IDX)=GAMMA_AQ(:,I)
          COMP_IDX = COMP_IDX +1
          !Set activity coefficients for the sub-components
          !equal to those for the main components
          DO J=2,NKAQ(I)
             GAMMA_AQ_RAOULT(:,COMP_IDX) = GAMMA_AQ(:,I)
             COMP_IDX=COMP_IDX + 1
          ENDDO !loop on ions
       ENDDO    !loop on main comp

       !Get the activity coefficients for Henry's law standard state
       CALL ACT_COEFF_HENRY_GET(         &
            X_AQ                         & !I [frc] molar fraction of components in aquous phase
            ,TMOL_AQ                     & !I [umol/m3] total number of umol in aquous phase
            ,GAMMA_AQ_RAOULT             & !I [-] Raoult's law standard state activity coefficents
            ,SI_AQ                       & !I [?] temperature dependent UNIFAC parameter
            ,GAMMA_AQ_HENRY              & !O [-] Henry's law standard state activity coefficients
            )

       !write things
       IF(LPRINT)THEN
          COMP_IDX=1
          DO I=1,NBSP
             write(6,*)"gamma_henry",I,  GAMMA_AQ_HENRY(:,COMP_IDX)
             COMP_IDX=COMP_IDX+1
             DO J=2,NKAQ(I)
                COMP_IDX=COMP_IDX+1
             ENDDO
          ENDDO
       ENDIF

       !Get the gas phase concentrations
       CALL GAS_FROM_HENRY(          &
            AERO_AQ                  &  !I [ug/m3] aerosol concentration of acid [H2A] and ions
            ,GAMMA_AQ_HENRY          &  !I [-] Henry's law standard state activity coefficient
            ,GAS_AQ                  &  !O [ug/m3] gas phase concentrations
            ,K                       &  !I [m3/ug and mole/kg] Henry's law and acid diss. constants
            ,LWC                     &  !I [ug/m3] Liquid water content
            ,NKAQ                    &  !I [nbr] number of sub and main species for component
            ,NBSP                    &  !I [nbr] number of main species to take into account)
            )

       !Knowing the concentration of H2A and acids (via the "nominator") we can
       !get the concentration of total aerosol by summing them and taking into account the 
       !difference in molecular weight between acids and ions
       CALL AERO_TO_totAER (        &
            AERO_AQ                 & !I [ug/m3] aquous phase concentrations (acid + ions)
            ,totAER_AQ              & !O [ug/m3] aquous phase concentrations (total = sum of acid+ions)
            ,nominator              & !I [-] transforming factor between acid and total
            ,NKAQ                   & !I [nbr]  number of sub and mean species
            ,NBSP                   & !I [nbr] number of main species to take into account
            )

       !Adjust mass balance for all components, both in aq phase
       !and in organic phase
       CALL MB_ADJUST(              &
            CB                      &  !I [ug/m3] total (aerosol + gas) concentration
            ,GAS_AQ                 &  !I [ug/m3] gas concentration in eql with aquous phase
            ,GAS_ORG                &  !I [ug/m3] gas concentration in eql with org phase
            ,totAER_AQ              &  !I/O [ug/m3] total aerosol concentration (acid+ions) aq phase
            ,AERO_ORG               &  !I/O [ug/m3] total aerosol concentration in org phase
            ,GAS                    &  !O [ug/m3] new guess for gas phase concentration
            )

       IF(LBOX)THEN
          DO I=1,NBSP
             write(6,'(a,i5,3e15.4)')"GAS ORG/AQ/NEW", I, GAS_ORG(:,I),GAS_AQ(:,I), GAS(:,I)
          ENDDO
          
          write(6,*)"budget"
          DO I=1,NBSP
             write(6,'(a,i5,3f17.11,e15.4)')"budget", I                     &
                  ,AERO_ORG(:,I)/CB(:,I)*100.&
                  ,totAER_AQ(:,I)/CB(:,I)*100.                             &
                  ,GAS(:,I)/CB(:,I)*100.d0                                 &
                  ,(CB(:,I)-totAER_AQ(:,I)-AERO_ORG(:,I)-GAS(:,I))/CB(:,I)*100.d0
          ENDDO
       ENDIF
       

       !Prepare for new iteration, get concentrations of H2A from knowing the corrected totAER_AQ
       CALL totAER_TO_AERO(         &
            totAER_AQ               & !I [ug/m3] total aerosol (acid+ions) in aq phase
            ,NOMINATOR              & !I [-] term which gives acid if we know total (acid+ions)
            ,AERO_AQ                & !O [ug/m3] acid and aerosol concentrations
            ,acHP                   & !I [mole/kg_{water}] proton concentration
            ,MW_SOA                 & !I [g/mol] molecular weight of species
            ,K                      & !I [m3/ug and mole/kg] Henry's law and acid diss. constants
            ,NKAQ                   & !I [nbr] number of sub and main components for a species
            ,NBSP                   & !I nbr] number of main species to take into account
            )

    ENDDO

    DEALLOCATE(X_ORG)
    DEALLOCATE(X_AQ)
    DEALLOCATE(GAMMA_ORG)
    DEALLOCATE(GAMMA_AQ)
    DEALLOCATE(GAMMA_AQ_RAOULT)
    DEALLOCATE(GAMMA_AQ_HENRY)

  IF (LHOOK) CALL DR_HOOK('MODE_SOAEQL:SOAEQL',1,ZHOOK_HANDLE)
  END SUBROUTINE SOAEQL

end module mode_soaeql
