!     ######spl
module mode_typea
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  use modd_unifacparam
  use mode_unifac
  use mode_soatinit, only: SI_GET
  use modd_glo, only: NBSPA, NAAEROA, NKA

  IMPLICIT NONE

  PUBLIC
  PRIVATE :: PARTFUNC2

  !***************************************************************************
  !Purpose:
  !Type A module: takes input of particle concentrations and returns
  !the values obtained after iterations
  !
  !Preconditions: This subroutine is called from amain, and takes the
  !the qualified guesses for aerosol concentrations as input.
  !Hopefully it finds the solutions for aerosol/gas partitioning
  !
  !Subroutines called:
  !unidriver and partition
  !(unidriver is always called with all 7 molecules that may  be present)
  !
  !Revisions: 1. Developed by Betty Pun, AER, Jan 99, under EPRI funding for
  !              prototype SOA module with 2 condensable compounds (malic acid,
  !              glyoxalic acid) and water
  !
  !           2. Modified November 99 to accept 6 model compounds + water under
  !              CARB funding (for the list of condensables see main)
  !
  !           3. Modified treatment of negative test concentrations so that
  !              program does not exit prematurely.  B. Pun Jan 2000.
  !**************************************************************************

  CONTAINS

  SUBROUTINE TYPEA ( &
       AERO          & !I/O [ug/m3] aerosol phase concentrations
       ,GAS          & !O [ug/m3] gas phase concentrations
       ,totA         & !I [ug/m3] total (aerosol+gas concentration)
       ,TEMPK        & !I [K] temperature
       ,LWC          & !I [ug/m3] Liquid water content of aerosols
       ,acHP         & !I [mol/kg_{water}] proton concentration
       )

    USE mode_soaeqlutl, only: XAQ_GET

    IMPLICIT NONE
    !INPUT
    REAL, DIMENSION(:,:), INTENT(INOUT)       :: AERO    !I/O [ug/m3] aerosol concentration (ions and molecules)
    REAL, DIMENSION(:,:), INTENT(IN)          :: totA    !I [ug/m3] total (gas+aerosol) concentrations
    REAL, DIMENSION(:), INTENT(IN)            :: TEMPK   !I [K]
    REAL, DIMENSION(:), INTENT(IN)            :: LWC     !I [ug/m3] liquid water content of aerosols
    REAL, DIMENSION(:), INTENT(IN)            :: acHP    !I [mol/kg_{water}] proton concentration

    !OUTPUT
    REAL, DIMENSION(:,:), INTENT(OUT)         :: GAS     !I/O [ug/m3] gas phase concentrations

    !Local, automatic arrays
    REAL, DIMENSION(SIZE(TEMPK,1))            :: TMOL      ![umol/m3] total moles in aquous phase
    REAL, DIMENSION(SIZE(TEMPK,1))            :: TMPX      ![-] temporary save variable for mole fraction
    REAL, DIMENSION(SIZE(GAS,1),SIZE(GAS,2))  :: totAER    ![ug/m3] total aerosol
    REAL, DIMENSION(SIZE(GAS,1),NBSPA,3)      :: ACIDFRC   ![frc] fraction of different acids

    !Local, allocatable
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: X        ![-] molar fraction of components
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: X2       ![-] molar fraction of components after one is dilluted
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: GAMMA    ![-] activity coefficients Henry's law std. state all comps
    REAL, ALLOCATABLE, DIMENSION(:,:)   :: GAMMAR   ![-] activity coefficients Raoult's law std. state main comps
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SI_A     ![?] temperature dependent activity coefficient

    !Local small counters etc
    REAL                               :: XX
    INTEGER                            :: COMP_IDX  ![idx] index for component number (1-13 currently)
    INTEGER                            :: I         ![idx] counter for main component (molecules)
    INTEGER                            :: J         ![idx] counter for sub components (ions)
    INTEGER                            :: ITER_IDX  ![idx] counter for number of iterations
    INTEGER, PARAMETER                 :: MAX_ITER=10 ![nbr] max number of iterations

    !Allocate memory
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_TYPEA:TYPEA',0,ZHOOK_HANDLE)
    ALLOCATE(X(SIZE(TEMPK,1),NBSPA+1))      !The extra component is water, could have used NMOL_A here
    ALLOCATE(X2(SIZE(TEMPK,1),NBSPA+1))     !The extra component is water, could have used NMOL_A here
    ALLOCATE(GAMMAR(SIZE(TEMPK,1),NBSPA+1)) !The extra component is water, could have used NMOL_A here
    ALLOCATE(GAMMA(SIZE(TEMPK,1),NAAEROA))
    ALLOCATE(SI_A(SIZE(TEMPK,1),NFUNC_A,NFUNC_A)) !Temperature dependent factor for act. coeff

    IF(LPRINT)THEN
    write(6,*)"                 "
    write(6,*)"*****************"
    write(6,'(a10,13e10.3)')"AERO",AERO
    write(6,'(a10,6e10.3)')"GAS",GAS
    write(6,'(a10,6e10.3)')"totA", totA
    write(6,'(a10,7e14.6)')"X",X
    write(6,*)"*****************"
    ENDIF

    !Get temperature dependent UNIFAC coefficient
    CALL SI_GET(                  &
         TEMPK                    & !I [K] temperature
         ,A_A                     & !I [units?] term in UNIFAC parameterization
         ,SI_A                    & !O [units?] term in UNIFAC parameterization
         ,NFUNC_A                 & !I [nbr] number of functional group in mixture
         )

    !Start the iterations to make sure things converge at this H+/LWC
    DO ITER_IDX = 1, MAX_ITER

       !Get the molar fraction of water
       CALL XAQ_GET(            &
            LWC                 & !I [ug/m3] liquid water content
            ,AERO               & !I [ug/m3] aquous phase aerosol concentrations of acids and ions
            ,X                  & !O [frc] molar fraction of main component in aq phase
            ,TMOL               & !O [umole/m3] total aquous phase aerosol
            ,MWA                & !I [g/mol] molecular weight of species
            ,NKA                & !I [nbr] number of components for main acid
            ,NBSPA              & !I [nbr] number of main species
            )

       !Get first Raoult's law activity coefficients
       CALL ACT_COEFF_GET(            &
            NU_A                      &   !I [nbr] number of functional groups in molecule I
            ,X                        &   !I [frc] mole fraction of molecule I
            ,QG_A                     &   !I [m2(?)] group surface area parameter
            ,GAMMAR                   &   !O [-] activity coefficient
            ,THTAGP_A                 &   !I [-] surface area ratio of groups (j) in pure component (i)
            ,Q_A                      &   !I [m2] surface area of pure component
            ,R_A                      &   !I [m3] total volume of pure component
            ,L_A                      &   !I [?] unifac parameter pure component
            ,SI_A                     &   !I [?] temperature dependent term
            ,NMOL_A                   &   !I [nbr] total number of molecules (including water)
            ,NFUNC_A                  &   !I [nbr] total number of functional groups (e.g. CH2, NO3 ..)
            )

       !Set activity coeffients for all components, ions and acids
       !to the Raoult's law standard state value for the real X
       COMP_IDX=1
       DO I=1,NBSPA
          GAMMA(:,COMP_IDX)=GAMMAR(:,I)
          COMP_IDX = COMP_IDX +1
          !Set activity coefficients for the sub-components
          !equal to those for the main components
          DO J=2,NKA(I)
             GAMMA(:,COMP_IDX) = GAMMAR(:,I)
             COMP_IDX=COMP_IDX + 1
          ENDDO !loop on ions
       ENDDO    !loop on main comp

       ! calculate solute activity coefficient at infinite dilution
       ! according to the standard state of Raoults Law
       COMP_IDX = 1
       DO I=1,NBSPA

          !Save the real value of "X"
          TMPX(:) = X(:,I)

          !Set component "I" to a value of inifinte dilution
          X(:,I) =  TMOL(:) * TINY

          !Make sure sum(X)=1 after modifying one value
          DO J=1,NBSPA+1
             X2(:,J) = X(:,J)/SUM(X(:,1:NBSPA+1))
          ENDDO

          ! gamma (Raoults law std st) at infinite dilution */
          CALL ACT_COEFF_GET(            &
               NU_A                      &   !I [nbr] number of functional groups in molecule I
               ,X2                       &   !I [frc] mole fraction of molecule I
               ,QG_A                     &   !I [m2(?)] group surface area parameter
               ,GAMMAR                   &   !O [-] activity coefficient
               ,THTAGP_A                 &   !I [-] surface area ratio of groups (j) in pure component (i)
               ,Q_A                      &   !I [m2] surface area of pure component
               ,R_A                      &   !I [m3] total volume of pure component
               ,L_A                      &   !I [?] unifac parameter pure component
               ,SI_A                     &   !I [?] temperature dependent term
               ,NMOL_A                   &   !I [nbr] total number of molecules (including water)
               ,NFUNC_A                  &   !I [nbr] total number of functional groups (e.g. CH2, NO3 ..)
               )

          !convert gamma for solutes to Henry's Law standard state
          !and set it equal for all main and sub-components "fxm: I am not sure where this equation comes from"
          DO J=1,NKA(I)
             GAMMA(:,COMP_IDX)=      &  !Value for real "X" and Henry's law std. state
                  GAMMA(:,COMP_IDX)  &  !Value for real "X" and Raoult's law std. state
                  / GAMMAR(:,I)         !Value for dilute "X" and Raoult's law std. state

             COMP_IDX=COMP_IDX+1
          ENDDO

          !Set back X to the value before the "inifinite dilution" value
          X(:,I) = TMPX(:)

       ENDDO !Loop on components

       !Now that we have the gammas, and the mole fractions
       !we can have get the gas and aerosol partitioning in thermodynamic
       !equilibrium with these conditions
       CALL PARTFUNC2(        &
            AERO              & !I/O [ug/m3] aerosol concentration in Henry's law eql
            ,GAS              & !O [ug/m3] gas concentrations
            ,totA             & !I [ug/m3] total aerosol concentration
            ,GAMMA            & !I [-] Henry's law std. state activity coefficients
            ,acHP             & !I [mol/kg_{water}] proton concentration
            ,LWC              & !I [ug/m3] Liquid water content
            )


       IF(LPRINT)THEN
       write(6,*)"*****************"
       write(6,'(a10,7e14.6)')"GAMMAR", GAMMAR
       write(6,'(a10,13e10.3)')"GAMMA", GAMMA
       write(6,'(a10,13e10.3)')"AERO",AERO
       write(6,'(a10,6e10.3)')"GAS",GAS
       write(6,'(a10,6e10.3)')"totA", totA
       write(6,'(a10,7e14.6)')"X",X
       write(6,*)"*****************"
       ENDIF

    ENDDO  !Loop on iterations

    IF(LBOX)THEN
       !Check mass balance after iterations
       COMP_IDX=1
       DO I=1,NBSPA
          totAER(:,I)=0.d0
          DO J=1,NKA(I)
             IF(J.eq.1)THEN
                XX=1.d0
             ELSE
                XX=MWA(COMP_IDX-1)/MWA(COMP_IDX)
             ENDIF
             totAER(:,I) = totAER(:,I) + AERO(:,COMP_IDX)*XX
             COMP_IDX=COMP_IDX+1
          ENDDO
       ENDDO

       ACIDFRC=0.d0
       COMP_IDX=1
       DO I=1,NBSPA
          DO J=1,NKA(I)

             IF(J.eq.1)THEN
                XX=1.d0
             ELSE
                XX=MWA(COMP_IDX-1)/MWA(COMP_IDX)
             ENDIF

             ACIDFRC(:,I,J)=AERO(:,COMP_IDX)*XX/totAER(:,I)
             COMP_IDX=COMP_IDX+1
          ENDDO
       ENDDO

       write(6,*)"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
       write(6,*)"TYPE A TOTALS"
       write(6,*)
       write(6,'(a10,6e10.3)')"GAS",GAS
       write(6,'(a10,6e10.3)')"totAER",totAER
       write(6,'(a10,6e10.3)')"totA",totA
       write(6,'(a10,6e10.3)')"aerfrc",totAER/totA
       write(6,'(a10,6e10.3)')"diff",totA-totAER-GAS
       write(6,*)'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'

       write(6,*)" TYPE A ACID FRACTIONS"
       do I=1,NBSPA
          write(6,'(a10,i5,3f12.5)')"frc : I ",I,acidfrc(:,I,1:3)
       ENDDO
    ENDIF !LBOX

    DEALLOCATE(X)
    DEALLOCATE(X2)
    DEALLOCATE(GAMMAR)
    DEALLOCATE(GAMMA)
    DEALLOCATE(SI_A)


  IF (LHOOK) CALL DR_HOOK('MODE_TYPEA:TYPEA',1,ZHOOK_HANDLE)
  END SUBROUTINE TYPEA


  !****************************************************************************

  SUBROUTINE PARTFUNC2(  &
       AERO              & !I/O [ug/m3] guessed aerosol concentration of main and sub components
       ,GAS              & !O [ug/m3] gas concentrations in eql with aerosol concentrations
       ,totA             & !I [ug/m3] total (gas+aerosol) concentration of main components
       ,GAMMA            & !I [-] Activity coefficients (Henry's law standard state)
       ,acHP             & !I [mol/kg_{water}] proton concentration
       ,LWC              & !I [ug/m3] Liquid water content
       )


    !Purpose: having guessed the aerosol phase concentrations, we want to check if this is
    !consistent with Henry's law equilibrium with gas phase.
    !If the resulting gas phase concentration is too high, we have reduce the aerosol phase concentrations

    !In amain we used K2=[HA-][H+]/[HA] where [X] represented concentrations of X. This was only
    !an approximation to get first guess. Here, we remember that [X] is ACTIVITY of X, and therefore
    ![X]=gamma_X*{X} where {X} is concentration

    !The equations given on top of amain transforms to
    !We have the following equations:
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
    !***********************************************************************************************

    IMPLICIT NONE

    !INPUT
    REAL, DIMENSION(:,:), INTENT(IN)        :: totA      !I [ug/m3] total (gas+aerosol)
    REAL, DIMENSION(:,:), INTENT(IN)        :: GAMMA     !I [-] Henry's law std. state act. coefficients
    REAL, DIMENSION(:), INTENT(IN)          :: acHP      !I [mol/kg_{water}] proton concentration
    REAL, DIMENSION(:), INTENT(IN)          :: LWC       !I [ug/m3] liquid water content

    !INPUT/OUTPUT
    REAL, DIMENSION(:,:), INTENT(INOUT)     :: AERO      !I/O [ug/m3] aerosol concentrations of all (molecs and ions)

    !OUTPUT
    REAL, DIMENSION(:,:), INTENT(OUT)     :: GAS         !I/O [ug/m3] gas phase concentrations

    !Local, automatic arrays
    REAL, DIMENSION(SIZE(AERO,1),SIZE(AERO,2)) :: AERO_NEW     ![ug/m3] aerosol concentration in Henry's law eql
    REAL, DIMENSION(SIZE(GAS,1),SIZE(GAS,2))  :: GAS_NEW1      ![ug/m3] gas in eql with guessed aerosol
    REAL, DIMENSION(SIZE(GAS,1),SIZE(GAS,2))  :: GAS_NEW2      ![ug/m3] gas in eql with guessed aerosol
    REAL, DIMENSION(SIZE(GAS,1),SIZE(GAS,2))  :: totAER_OLD    ![ug/m3] total aerosol (main+sub components), not in MB
    REAL, DIMENSION(SIZE(GAS,1),SIZE(GAS,2))  :: totAER_NEW    ![ug/m3] total aerosol (main+sub components), in MB
    REAL, DIMENSION(SIZE(GAS,1),SIZE(GAS,2))  :: totA_NEW      ![ug/m3] new total which can deviate from mass balance

    REAL, DIMENSION(SIZE(totA,1))           :: THISTERM        ![-] term in solving equilibrium for main component
    REAL, DIMENSION(SIZE(totA,1))           :: NOMINATOR       ![-] term in solving equilibrium for main component

    !Local, small counter etc
    INTEGER                                 :: I          ![idx] counter for main component
    INTEGER                                 :: J          ![idx] counter for sub-components (ions)
    INTEGER                                 :: COMP_IDX   ![idx] pointer to right component in AERO array
    INTEGER                                 :: COMP_IDX2  ![idx] pointer to right index in K-array (eql constants)

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_TYPEA:PARTFUNC2',0,ZHOOK_HANDLE)
    COMP_IDX = 1
    DO I=1,NBSPA

       !***********PART 1**********************************************
       !Guess that the gas is the concentration in Henry's law equilibrium
       !Note that this IS the correct gas, given that our guess for aerosol is OK
       GAS_NEW1(:,I)=MIN(totA(:,I)-TINY, &                            !Total aerosol concentration
            AERO(:,COMP_IDX)*GAMMA(:,COMP_IDX)/(LWC(:)*K_A_298(COMP_IDX))) !Henry's law equilibrium

       !At this point, get first totAER which is input aerosol concentration
       COMP_IDX2=COMP_IDX
       totAER_OLD(:,I)=AERO(:,COMP_IDX)
       DO J=2,NKA(I)
          COMP_IDX2=COMP_IDX2 +1
          totAER_OLD(:,I) = totAER_OLD(:,I) + AERO(:,COMP_IDX2)*MWA(COMP_IDX)/MWA(COMP_IDX2)
       ENDDO

       !totA_NEW is GAS+AEROSOL which might deviate from mass balance
       totA_NEW(:,I)=totAER_OLD(:,I) + GAS_NEW1(:,I)

       !Scale with this one to obtain totAER and GAS_NEW in accordance with mass balance
       totAER_NEW(:,I)=totAER_OLD(:,I)*totA(:,I)/totA_NEW(:,I)

       !Get gas concentration in accordance with mass balance
       GAS_NEW1(:,I)=GAS_NEW1(:,I)*totA(:,I)/totA_NEW(:,I)

       !write(6,'(a,i5, 6f10.5)')"GAS_FRC A0",I,100.*(GAS_NEW1(:,I))/totA(:,I), (GAS_NEW1(:,I)+totAER_NEW(:,I))/totA(:,I)

       !Now, force the aerosol concentrations with this new gas concentrations
       !Get the nominator in the mass balance equation (see mode_amain)
       NOMINATOR(:)=1.d0
       THISTERM(:)=1.d0

       COMP_IDX2 = COMP_IDX       !counter for sub-components for which we sum
       DO J=2,NKA(I)
          COMP_IDX2=COMP_IDX2 + 1
          THISTERM(:)=THISTERM(:)*K_A_298(COMP_IDX2)/acHP(:)*MWA(COMP_IDX2-1)/MWA(COMP_IDX2)
          NOMINATOR(:)=NOMINATOR(:) + THISTERM(:)
       ENDDO

       !Get new aerosol concentration of main component
       !This IS the correct aerosol concentration given that the gas phase was OK
       !However, this aerosol concentration can differ from the orignal one
       AERO_NEW(:,COMP_IDX) = totAER_NEW(:,I)/NOMINATOR(:)
       COMP_IDX = COMP_IDX +1

       !Get aerosol concentrations of sub-components
       !Simply use equilibrium constants in aquous phase
       DO J=2,NKA(I)
          !Note that we skip activity coefficients from this calculations since
          !they are anyway equal for all ions which are "children" of same main comp
          AERO_NEW(:,COMP_IDX) = K_A_298(COMP_IDX) * AERO_NEW(:,COMP_IDX-1)/acHP(:)
          COMP_IDX = COMP_IDX +1
       ENDDO

       !*************END PART 1**********************************************
       !Now we have obtained GAS_NEW1 from eql with original mass balance, and we have
       !Gotten back a new TOTAER_NEW from MB with GAS_NEW1********************
       !**********************************************************************

       !Do an extra calculation here to make sure that the calculations don't oscillate
       !Now, the aerosol concentration is the one in equilibrium with the gas phase (gas_new)
       !Remember that GAS_NEW was obtained from Henry's law with guessed aerosol concentrations.
       !The aerosol concentrations are now reduced if the original guess was too high

       !****************PART 2***********************************************
       !Get another guess for gas, obtained from the new aerosol conc
       !********************************************************************

       !First: go back to main component
       COMP_IDX = COMP_IDX - NKA(I)

       !Re-calculate gas using Henry' law
       GAS_NEW2(:,I)=MIN(totA(:,I)-TINY, &                                !Total aerosol concentration
            AERO_NEW(:,COMP_IDX)*GAMMA(:,COMP_IDX)/(LWC(:)*K_A_298(COMP_IDX)))  !Henry's law equilibrium


       !totA_NEW is GAS+AEROSOL which might deviate from mass balance
       totA_NEW(:,I)=totAER_NEW(:,I) + GAS_NEW2(:,I)

       !Get total aerosol in accordance with mass balance
       totAER_NEW(:,I)=totAER_NEW(:,I)*totA(:,I)/totA_NEW(:,I)

       !Get gas in accordance with mass balance
       GAS_NEW2(:,I)=GAS_NEW2(:,I)*totA(:,I)/totA_NEW(:,I)

       !write(6,'(a,i5,6f10.5)')"GAS_FRC A1",I,100.*(GAS_NEW2(:,I))/totA(:,I),(GAS_NEW2(:,I)+totAER_NEW(:,I))/totA(:,I)

       !************END PART 2 ***************************************

       !***********PART 3*********************************************
       !Here we assume that both guesses for gas are OK
       !so we use the mean of them as the "true" answer
       !and obtain totAER from the mass balance again

       GAS_NEW2(:,I)=0.5*(GAS_NEW2(:,I) + GAS_NEW1(:,I))
       totAER_NEW(:,I)=totA(:,I) - GAS_NEW2(:,I)           !Total aerosol obtained from mass balance

       !*********END PART 3*********************************************

       !**********PART 4 *********************************************
       !Get the aerosol concentrations assuming that totAER_NEW is the true
       !total aerosol concentration
       !****************************************************************

       !Get the nominator in the mass balance equation (see mode_amain)
       NOMINATOR(:)=1.d0
       THISTERM(:)=1.d0
       COMP_IDX2 = COMP_IDX
       DO J=2,NKA(I)
          COMP_IDX2=COMP_IDX2 + 1
          THISTERM(:)=THISTERM(:)*K_A_298(COMP_IDX2)/acHP(:)*MWA(COMP_IDX2-1)/MWA(COMP_IDX2)
          NOMINATOR(:)=NOMINATOR(:) + THISTERM(:)
       ENDDO

       !Get new aerosol concentration of main component
       AERO_NEW(:,COMP_IDX) = totAER_NEW(:,I)/NOMINATOR(:)

       !This is also the output
       AERO(:,COMP_IDX) = AERO_NEW(:,COMP_IDX)

       !Prepare for next component
       COMP_IDX = COMP_IDX +1

       !Get output aerosol concentrations of sub-components
       !Simply use equilibrium constants in aquous phase

       DO J=2,NKA(I)
          !New iterated concentrations of ions.
          !Note that we skip activity coefficients from this calculations since
          !they are anyway equal for all ions which are "children" of same main comp
          AERO_NEW(:,COMP_IDX) = K_A_298(COMP_IDX) * AERO_NEW(:,COMP_IDX-1)/acHP(:)
          !This is also the output
          AERO(:,COMP_IDX) = AERO_NEW(:,COMP_IDX)
          !Prepare for next component
          COMP_IDX = COMP_IDX +1
       ENDDO

       !Check new mass balance
       !COMP_IDX2=COMP_IDX - NK(I)
       !totAER_NEW(:,I)=AERO(:,COMP_IDX2)
       !DO J=2,NK(I)
       !   COMP_IDX2=COMP_IDX2 + 1
       !   totAER_NEW(:,I) = totAER_NEW(:,I) + AERO(:,COMP_IDX2)*MW(COMP_IDX-NK(I))/MW(COMP_IDX2)
       !ENDDO

       !Set output gas concentrations here
       GAS(:,I) = GAS_NEW2(:,I)
       !write(6,'(a,i5,6f10.5)')"GAS_FRC A2",I,100.*(GAS(:,I))/totA(:,I),(GAS(:,I)+totAER_NEW(:,I))/totA(:,I)

       !Set output gas concentration

    ENDDO   !Loop on components
  IF (LHOOK) CALL DR_HOOK('MODE_TYPEA:PARTFUNC2',1,ZHOOK_HANDLE)
  END SUBROUTINE PARTFUNC2
  !********************************************************************

END MODULE mode_typea
