!     ######spl
MODULE MODE_AMAIN
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  USE MODD_BINSOLU, only : &
       molalbinA

  USE MODD_GLO, only: &
       NBSPA     &
       ,MWA      &
       ,VPA_298  &
       ,NKA      &
       ,K_A_298

  IMPLICIT NONE

  PUBLIC
  PRIVATE :: SATURATION

CONTAINS
!***********************************************************************
!Purpose: Starting point of Type A routine.  Performs 2 functions: (1)
!         gas/particle partition of organic compounds based on available
!         water or saturation, (2) calculate water associated with organic
!         solutes, (3) calculate total organic anion concnetrations,
!
!Arguments: 1. *aero: initial guesses of Ai in same units as inputs ug/m3
!              also contain final output from newt
!           2. *aeros: solid concentrations (also used to temporarily store
!                   input PM concentrations in LWC = 0 absorption case with
!              no newton, and PM from absorption before water is added)
!
!Return:    1. deltaLWC is the output for water associated with the organics
!              in microgram/m3 air
!
!Data Needed:
!        Input read in main and passed unsing global variables
!               totA[i] (ug/m3 air), acHP (mole/kg water),
!                LWC (ug/m3 air), RH (between 0 and 1), temperature (K)
!
!    Parameters: NK = no. of eq. relationship for each solute
!                 K = partition parameters H and K
!                     Ai, Gi, W in ug per m3 air units
!                     gamma in mole fraction units, change reference
!                     state to Henry's law
!                     {H+} in moles per kg solvent,
!                MW = Molecular weights in the order defined below
!                     borrowing MW[0] to store MW(water)
!               DRH = deliquescence humidities for molecules
!                     (still need data)
!                VP = vapor pressure in units of mass (ug) per m3 air
!                     pure liquid-gas partitioning:
!                     (1) deliquescent species
!                          (2) input LWC = 0  from inorganic species 5-29-99
!Notes:
!    Molecules (6)
!    1. Propandioic acid (twice dissociating)
!    2. C8 dien-dioic acid with an aldehyde branch (twice dissociating)
!    3. C8 hydroxy-dien-dial (non-dissociating)
!    4. C9 hydroxy-carbonyl acid with one double bond (dissociating)
!    5. C10 hydroxy-carbonyl aldehyde (non-dissociating)
!    6. butandioic acid (twice dissociating)
!
!    Solutes (13)
!    1. Propandioic acid (H2A1)
!    2. Hydrogen Propanoate (HA1-)
!    3. Propanoate (A1=)
!    4. C8 aldehyde-branched dien-dioic acid (H2A2)
!    5. C8 aldehyde-branched hydrogen dienoiate (HA2-)
!    6. C8 aldehyde-branched dienoiate (A2=)
!    7. C8 hydroxy diendial (A3)
!    8. C9 hydroxy carbonyl unsaturated acid (HA4)
!    9. C9 hydroxy carbonyl unsaterated anion (A4-)
!    10. C10 hydroxy carbonyl aldehyde (A5)
!    11. butandioic acid (H2A6)
!    12. hydrogen butanoate (HA6-)
!    13. butanoiate (A6=)
!
!    Cases:
!    1. LWC > 0
!       A) RH > DRH[i] - aqueous phase (with ions)
!       B) RH < DRH[i]
!          (i)  totA > VP - gas phase = VP
!                         - solid phase (no water associated with PM phase i)
!          (ii) totA < VP - gas phase only
!    2. LWC = 0 (option 1, saturation)
!       A) RH > DRH[i]
!          (i)  totA > VP - gas phase = VP
!                    - aq phase molecules only = totA[i] - VP[i]
!          (ii) totA < VP - gas phase only
!       B) RH < DRH[i]
!          (i)  totA > VP - gas phase = VP
!                         - solid phase (no water associated with PM phase i)
!          (ii) totA < VP - gas phase only
!    2. LWC = 0 (option 2, absorption)
!       A) RH > DRH - gas/liquid partition, liquid phase
!                     associated with organic water (molecules, no ions)
!       B) RH < DRH - gas/liquid partition, no water in liquid phase
!
!Revisions: 1. Developed by Betty Pun, AER, Jan 99 under EPRI for prototype
!              Type A module with 2 compounds: malic acid and glyoxalic acid
!              using newt, the globally convergent multi-dimensional
!              Newton's method to solve the non-linear simultaneous equations.
!              NR convension index 1 .. n; most arrays run from 1 to NSP+1
!
!           2. Added code May 99 to deal with LWC input = 0
!              calculate gas-PM partition based on VPsat
!              calculate water associated with organics (H2M RH < DRH)
!
!           3. Under CARB funding, modified code October 99 to
!              perform the partition of 6 compounds.  Removed hard-wired
!              code regarding equations solved
!
!           4. Modified to comply with Models-3 coding standard, Betty Pun,
!              AER, November 99
!
!           5. Included flag to solve for water content based on ZSR with
!              binary solution characteristics stored in file.
!
!           6. Combine Type A and B, Betty Pun Apr 00
!           6a.   bkp 6/00 change criteria to deal with low LWC cases:
!                 old criteria was:
!                 if ((NK[i] >= 2) && (K[jeq]/acHP > Critsol)) {
!           7. Add option to do absorption when LWC = 0 11/00
!           8. For 3-D, add option to solve absorption based on fixed PM
!****************************************************************************/


  SUBROUTINE AMAIN (  &
       totA           & !I [ug/m3] total (gas + aerosol + ions) of component
       ,AERO          & !O [ug/m3] liquid aerosol concentrations
       ,AEROS         & !O [ug/m3] solid aerosol concentration
       ,GAS           & !O [ug/m3] gaseous concentration
       ,LWC           & !I [ug/m3] liquid water content already available for partitioning
       ,acHP          & !I [umol/kg_{water}] proton concentration
       ,deltaLWC      & !O [ug/m3] liquid water content associated with organics
       ,ORGANION      & !O [mole/m3] negative charge associated with organics
       ,RH            & !I [0-1] relative humdity
       ,TEMPK         & !I [K] Temperature
       )

    USE mode_typea
    use mode_zsrpun

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

    IMPLICIT NONE

    !INPUT
    REAL, DIMENSION(:), INTENT(IN)   :: TEMPK   !I [K] temperature
    REAL, DIMENSION(:), INTENT(IN)   :: RH      !I [0-1] relative humidity
    REAL, DIMENSION(:,:), INTENT(IN) :: totA    !I [ug/m3] total (gas+aerosol) concentration
    REAL, DIMENSION(:), INTENT(IN)   :: LWC     !I [ug/m3] Liquid water content already ready for partitioning
    REAL, DIMENSION(:), INTENT(IN)   :: acHP    !I [umol/kg_{water}] proton concentration

    !OUTPUT
    REAL, DIMENSION(:,:), INTENT(OUT)  :: GAS      !O [ug/m3] gas phase concentrations
    REAL, DIMENSION(:,:), INTENT(OUT)  :: AERO     !O [ug/m3] liquid aerosol concentrations
    REAL, DIMENSION(:,:), INTENT(OUT)  :: AEROS    !O [ug/m3] solid aerosol concentrations
    REAL, DIMENSION(:), INTENT(OUT)    :: deltaLWC !O [ug/m3] liquid water content assiciated with organics
    REAL, DIMENSION(:), INTENT(OUT)    :: ORGANION !O [mole/m3] negative charge associated with organics

    !LOCAL, AUTOMATIC ARRAYS
    REAL, DIMENSION(SIZE(TEMPK))       :: THISTERM    !Term in nominator in mass balance eqn (see above)
    REAL, DIMENSION(SIZE(TEMPK))       :: NOMINATOR   !Nominator in mass balance eqn (see above)
    LOGICAL, DIMENSION(SIZE(GAS,1),SIZE(GAS,2)) :: LCRITSOL ![T/F] True if soluble species

    !LOCAL, SMALL COUNTERS AND OTHER VARS
    INTEGER                      :: I            ![idx] counter for main components (1-->NBSPA)
    INTEGER                      :: J            ![idx] counter for sub-componnets or ions
    INTEGER                      :: COMP_IDX     ![idx] index for aerosol components (1-->NAAEROA)
    INTEGER                      :: COMP_IDX2    ![idx] index for aerosol sub-components (1-->NAAEROA)
    REAL, PARAMETER              :: HIGHFRC=0.999999 ![frc] guess for fraction of total in aq. phase if very soluble


    !******************************************************************************

    !Start part which will make guesses for the
    !aerosol concentrations before the iterations start
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_AMAIN:AMAIN',0,ZHOOK_HANDLE)
    COMP_IDX = 1
    DO I=1,NBSPA
       !Check for places where we think that the main acid
       !can dissolve well in the aquous phase, that is:
       !1) [H+] is low
       !2) [LWC] is high
       !3) the main aerosol can dissolve at least once
       LCRITSOL(:,I)=  &
            K_A_298(COMP_IDX)*LWC(:)/acHP(:).GT.CRITSOL.AND.NKA(I).ge.2

       !Prepare for next sub-component
       COMP_IDX=COMP_IDX+1

       !Skip sub-components
       DO J=2,NKA(I)
          COMP_IDX = COMP_IDX +1
       ENDDO
    ENDDO

    IF(LBOX)THEN
       DO I=1,NBSPA
          IF(LCRITSOL(1,I))write(6,*)"comp", I," very soluble"
       ENDDO
    ENDIF

    COMP_IDX = 1
    !Start loop on non-dissociated components
    DO I=1,NBSPA

       !Get the nominator in the mass balance equation above
       !is used to find the concentration of the undissociated species (e.g. H2A)
       NOMINATOR(:)=1.d0
       THISTERM(:)=1.d0

       COMP_IDX2 = COMP_IDX
       DO J=2,NKA(I)
          COMP_IDX2=COMP_IDX2 + 1     !go to first sub-component (ion)
          THISTERM(:)=THISTERM(:)*K_A_298(COMP_IDX2)   &
               /acHP(:)*MWA(COMP_IDX2-1)/MWA(COMP_IDX2)

          NOMINATOR(:)=NOMINATOR(:) + THISTERM(:)
       ENDDO

       !Check if main component is soluble
       WHERE(LCRITSOL(:,I)) !Very soluble species

          !Get concentration of main component from mass balance equation
          AERO(:,COMP_IDX) = HIGHFRC*totA(:,I)  &  !guess: totAER = totA*HIGHFRC ~= totA
               /NOMINATOR(:)

       ELSEWHERE

          !Places where the acid is not so soluble, assume totA ~=gas
          !use Henry's law on non-dissociative acid
          AERO(:,COMP_IDX) = totA(:,I)  &    !Guess: GAS =totA
               * K_A_298(COMP_IDX) * LWC(:) &
               /(1.d0 + K_A_298(COMP_IDX)*LWC(:))

       END WHERE

       IF(LBOX)THEN
          write(6,'(a,2i5,2e15.5)')"guess", I, COMP_IDX, totA(:,I), AERO(:,COMP_IDX)
       ENDIF

       !Prepare for next component
       COMP_IDX = COMP_IDX +1

       !Set the concentrations of the sub-components
       DO J=2,NKA(I)
          AERO(:,COMP_IDX) = K_A_298(COMP_IDX) * AERO(:,COMP_IDX-1)/acHP(:)
          IF(LBOX)THEN
             write(6,'(a,2i5,2e15.5)')"guess", I, COMP_IDX, totA(:,I), AERO(:,COMP_IDX)
          ENDIF
          COMP_IDX = COMP_IDX +1
       ENDDO
       IF(LBOX)write(6,*)"  "

    ENDDO  !Loop on main components

    !Now that we have guesses, call type-A to iterate to solution
    CALL TYPEA(       &
         AERO         & !I/O [ug/m3] aerosol concentrations (main and ion)
         ,GAS         & !O   [ug/m3] gas phase concentrations
         ,totA        & !I   [ug/m3] total (gas+aer) concentrations
         ,TEMPK       & !I [K] temperature
         ,LWC         & !I [ug/m3] liquid water content available for partitioning
         ,acHP        & !I [mole_{H+}/kg_{water}] proton concentration
         )

    !Now we have all components in aquous phase, get the new LWC
    !We need the new (H+ concentration too!)
    CALL ZSRPUN(      &
         AERO         & !I [ug/m3] aquous aerosol concentration
         , RH         & !I [0-1] relative humidity
         ,deltaLWC    & !O [ug/m3] LWC associated with organics
         ,MOLALBINA   & !I [umol/ug_{water}] molality of binary solutions
         ,NBSPA       & !I [nbr] number of species to take into account
         ,NKA         & !I [nbr] number of sub and main component for main component
         ,MWA         & !I [g/mol] molecular weight of species
         )

    !For now, we don't consider solid aerosols, so put it to zero
    AEROS(:,:)=0.d0

    !Get negative charge associated with organics
    ORGANION(:)=0.d0
    COMP_IDX=1
    DO I=1,NBSPA
       !No charge for first component
       COMP_IDX=COMP_IDX+1
       DO J=2,NKA(I)
          ORGANION(:)=ORGANION(:)               &
               + AERO(:,COMP_IDX)/MWA(COMP_IDX)  & !moles of ion component
               *dble(J-1)                         !charge of ion component

          COMP_IDX=COMP_IDX+1
       ENDDO
    ENDDO

    IF(LBOX)write(6,*)"LWC, deltaLWC", LWC, deltaLWC

  IF (LHOOK) CALL DR_HOOK('MODE_AMAIN:AMAIN',1,ZHOOK_HANDLE)
  END SUBROUTINE AMAIN

  !***************************************************************

  SUBROUTINE saturation(&
       totA             &   !I [ug/m3] total (gas+aerosol) concentration
       ,vp              &   !I [ug/m3] saturation vapor pressure
       ,pmconc          &   !O [ug/m3] solid aerosol concentration
       ,gasconc         &   !O [ug/m3] gas phase concentration
       )

    !**************************************************************
    ! purpose: given total amount of a compound,
    !          calculate partition based on saturation
    !          if tot > vp, gas = vp, pm = (tot - vp)
    !          if tot < vp, gas = tot pm = 0
    ! arguments: tot      total amount (microgram/m3 air) of compound
    !            vp       vapor pressure (microgram/m3 air) of compound
    !            *pmconc  pointer to output of particulate-phase concentration
    !            *gasconc pointer to output of gas-phase concentration
    ! history: 1. coded 11/28/00 BKP to replace code in 2 locations:
    !             - LWC > 0 and RH < DRH
    !             - LWC = 0 any RH
    !
    !****************************************************************/

    IMPLICIT NONE
    REAL, DIMENSION(:,:), INTENT(IN)  :: totA       !I [ug/m3] total (gas+aerosol) concentration
    REAL, DIMENSION(:), INTENT(IN)    :: VP         !I [ug/m3] saturation vapor pressure
    REAL, DIMENSION(:,:), INTENT(OUT) :: PMCONC     !O [ug/m3] solid aerosol concentration
    REAL, DIMENSION(:,:), INTENT(OUT) :: GASCONC    !O [ug/m3] gas phase concentration
    INTEGER                           :: I          ![idx] conter for aerosol species

    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_AMAIN:SATURATION',0,ZHOOK_HANDLE)
    DO I=1,NBSPA

       WHERE(totA(:,I) > VP(I))
          PMCONC(:,I) = totA(:,I) - VP(I)
          GASCONC(:,I) = VP(I)
       ELSEWHERE
          GASCONC(:,I)= totA(:,I)
          PMCONC(:,I)=0.d0
       END WHERE

    ENDDO

  IF (LHOOK) CALL DR_HOOK('MODE_AMAIN:SATURATION',1,ZHOOK_HANDLE)
  END SUBROUTINE SATURATION

END MODULE mode_amain
