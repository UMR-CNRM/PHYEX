!     ######spl
module mode_typeb
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK

  use modd_glo
  use modd_unifacparam  ! definition of variables for unifac
  use mode_unifac       ! UNIFAC module
  use mode_soatinit, only: SI_GET

  IMPLICIT NONE

  PUBLIC
  PRIVATE :: KPART_PUN_GET

!***************************************************************************
!Purpose: Partition of Type B condensables,
!        Ai orig used in Unifac --> gamma,
!        gamma --> Ki; Ki --> the distribution of OC
!        Particle phase OC (Ai final) calculated
!        need f = Ai orig - Ai final = 0, f = delta aerosol is output
!
!Notes: 1. defined (hardcoded in glo.h or glodef.h) values used:
!!          MWaom (non-volatile) = 280
!          NBSPAOM = 5,
!          xaom = {0.4, 0.05, 0.15, 0.12, 0.28};
!          default breakdown of AOM
!          compound             mass frac.      mole frac.
!         C24 alkanoic acid         0.50            0.40
!         C18 alkenoic acid         0.05            0.05
!         actonyl syringol          0.11            0.15
!         C20 alkane                0.17            0.12
!         arom. dicarboxylic acid   0.17            0.28
!
!Revision history:  Developed by Betty Pun, AER, Jan 99 Under EPRI
!                   Modified by Betty Pun, AER, Nov 99 Under CARB
!                   1. increase the number of partitioning compound
!                   2. conform to models-3 coding standard
!                   3. allow the selection of equations to solve when
!                      using Newt.
!                   4. Rewritten to f90 and MESONH style: Alf Grini, CNRM
!***************************************************************************/

contains

  subroutine TypeB(    &
        PAOM           & !I [ug/m3] primary aerosol concentrations
       ,AERO           & !I/O [ug/m3] "iteratable" secondary aerosol concentrations
       ,GAS            & !I/O [ug/m3] "iteratable" condenseable gas concentrations
       ,CB             & !I   [ug/m3] total (gas+aerosol) semi volatile species
       ,TEMPK          & !I   [K] temperature
       )

    IMPLICIT NONE

    !INPUT/OUTPUT
    REAL, DIMENSION(:), INTENT(IN)       :: PAOM    ![ug/m3] Non volatile organic absorbing compounds
    REAL, DIMENSION(:,:), INTENT(INOUT)  :: AERO    ![ug/m3] iteratable secondary aerosol concentration
    REAL, DIMENSION(:,:), INTENT(INOUT)  :: GAS     ![ug/m3] iteratable condenseable gas concentration
    REAL, DIMENSION(:,:), INTENT(IN)     :: CB      ![ug/m3] total (gas+aerosol) semi volatile (B)  species
    REAL, DIMENSION(:), INTENT(IN)       :: TEMPK   ![K] temperature

    !LOCAL, WORKING ARRAYS (AUTOMATIC ARRAYS)
    REAL, DIMENSION(SIZE(AERO,1),SIZE(AERO,2)) :: AERO_OLD  ![ug/m3] old secondary aerosol concentrations
    REAL, DIMENSION(SIZE(AERO,1),SIZE(AERO,2)) :: AERO_NEW  ![ug/m3] new secondary aerosol concentrations
    REAL, DIMENSION(SIZE(AERO,1),SIZE(AERO,2)) :: GAS_OLD   ![ug/m3] old condenseable gas concentrations
    REAL, DIMENSION(SIZE(AERO,1),SIZE(AERO,2)) :: GAS_NEW1  ![ug/m3] new guess for gas OK with mass balance
    REAL, DIMENSION(SIZE(AERO,1),SIZE(AERO,2)) :: GAS_NEW2  ![ug/m3] new guess for gas also OK with mass balance
    REAL, DIMENSION(SIZE(AERO,1),SIZE(AERO,2)) :: CB_NEW    ![ug/m3] "mass balance" not equal to CB
    REAL, DIMENSION(SIZE(AERO,1))              :: M0_OLD    ![ug/m3] old value for total aerosol mass
    REAL, DIMENSION(SIZE(AERO,1))              :: M0_NEW    ![ug/m3] new value of total aerosol mass
    REAL, DIMENSION(SIZE(AERO,1),SIZE(AERO,2)) :: KPART     ![m3/ug] partitioning coefficient

    REAL, DIMENSION(SIZE(AERO,1)) :: TMOLAOM                ![umol/m3] total mole primary absorbing organic medium
    REAL, DIMENSION(SIZE(AERO,1)) :: TMOL                   ![umol/m3] total mole of all organic b in aerosol phase
    REAL, DIMENSION(SIZE(AERO,1)) :: TMOLINV                ![m3/umol] 1/TMOL
    REAL, DIMENSION(SIZE(AERO,1)) :: MWOM                   ![g/mol] average molar weight of semi volatile organics

    !LOCAL, WORKING ARRAYS (ALLOCATABLES)
    REAL, ALLOCATABLE, DIMENSION(:,:)  :: X                 ![frc] molar fraction all organics
    REAL, ALLOCATABLE, DIMENSION(:,:)  :: GAMA              ![-] activity coefficient for all species
    REAL, ALLOCATABLE, DIMENSION(:,:)  :: GAMA_B            ![-] activity coefficients for partitioning species
    REAL, ALLOCATABLE, DIMENSION(:,:)  :: ORGP              ![umol/m3] concentration of primary absorbing species
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: SI_B             ![?] temperature dependent coefficient in UNIFAC

    INTEGER                            :: ITER_IDX          !Counter for number of iterations
    INTEGER                            :: I                 !Counter for components
    INTEGER, PARAMETER                 :: MAX_ITER=10       !Max number of iterations

    !Allocate memory for X, GAMA, ORGP
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_TYPEB:TYPEB',0,ZHOOK_HANDLE)
    ALLOCATE(X(SIZE(AERO,1),NBSPAO_PUN+NBSPB))
    ALLOCATE(GAMA(SIZE(AERO,1),NBSPAO_PUN+NBSPB))
    ALLOCATE(GAMA_B(SIZE(AERO,1),NBSPB))
    ALLOCATE(ORGP(SIZE(AERO,1),NBSPAO_PUN))
    ALLOCATE(SI_B(SIZE(AERO,1),NFUNC_B,NFUNC_B))

    !Initialize the "old" arrays with guessed concentrations
    AERO_OLD(:,:)=AERO(:,:)  !initial Guess for aerosol
    GAS_OLD(:,:)=GAS(:,:)    !initial Guess for gas


    IF(LBOX)THEN
       write(6,*)"  "
       write(6,*)"************************START TYPE B VALUES********"
       write(6,'(a10,5e14.5)')"AERO",AERO
       write(6,'(a10,5e14.5)')"GAS",GAS
       write(6,'(a10,5e14.5)')"CB",CB
       write(6,'(a10,e14.5)')"PAOM",PAOM
       write(6,*)"********************************************"
       write(6,*)"  "
    ENDIF

    !Initialize the temperature dependent UNIFAC coefficients outside iterations
    CALL SI_GET(                 &
         TEMPK                   & !I [K] temperature
         ,A_B                    & !I [units?] term in UNIFAC parameterization
         ,SI_B                   & !O [units?] term in UNIFAC parameterization
         ,NFUNC_B                & !I [nbr] number of functional group in mixture
         )

    !Get total moles of primary aerosol (umol/m3)
    DO I=1,NBSPAO_PUN
       TMOLAOM(:) =    &     !umol/m3
            PAOM(:)    &     !ug/m3
            / MWAOM          !g/mol ==> umol/m3
    ENDDO

    !Distribute the PAOM into components (XAOM) which are already fixed in modd_glo
    DO I=1,NBSPAO_PUN
       ORGP(:,I) =       & !umol/m3
            TMOLAOM(:)   & !umol/m3
            * XAOM(I)      !frc  ==> umol/m3
    ENDDO

    !Get the value for M0 (total organic in aerosol phase)
    M0_OLD(:)=PAOM(:)
    DO I=1,NBSPB
       M0_OLD(:)=M0_OLD(:) + AERO_OLD(:,I)  !ug/m3
    ENDDO

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    !From here on, we iterate
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DO ITER_IDX = 1,MAX_ITER

       !Initialize the "total mol" array with only primary aerosols
       TMOL(:) = TMOLAOM(:)

       !Add the guessed aerosol concentrations to the TMOL array
       DO I=1,NBSPAO_PUN
          TMOL(:)= TMOL(:) &
               + AERO_OLD(:,I)/MWB(I)  ![umol/m3] of condenseable organic
       ENDDO

       !Get the inverse of TMOL to multiply with this variable later
       TMOLINV(:) = 1.d0/TMOL(:)

       !Get the mole fraction of each primary constituent
       DO I = 1,NBSPAO_PUN
          X(:,I) = ORGP(:,I) * TMOLINV(:)
       ENDDO

       !Get molar fraction for the partitioning species
       DO I = 1,NBSPB
          X(:,NBSPAO_PUN+I) = AERO_OLD(:,I) / MWB(I) * TMOLINV(:)
       ENDDO

       !Get average molar weight of non partitioning organic species
       !Initalize molecular weight to zero
       MWOM(:)  = 0.d0
       DO I=1, NBSPAO_PUN
          MWOM(:) = MWOM(:) + X(:,i) * MWAOM
       ENDDO

       !Add the partitioning species
       DO I=1, NBSPB
          MWOM(:) = MWOM(:) + X(:,NBSPAO_PUN+I) * MWB(I)
       ENDDO

       !Get activity coefficients given the current X
       CALL ACT_COEFF_GET(&
            NU_B              &   !I [nbr] number of functional groups in molecule I
            ,X                &   !I [frc] mole fraction of molecule I
            ,QG_B             &   !I [m2(?)] group surface area parameter
            ,GAMA             &   !O [-] activity coefficient
            ,THTAGP_B         &   !I [-] surface area ratio of groups (j) in pure component (i)
            ,Q_B              &   !I [m2] surface area of pure component
            ,R_B              &   !I [m3] total volume of pure component
            ,L_B              &   !I [?] unifac parameter pure component
            ,SI_B             &   !I [?] temperature dependent term
            ,NMOL_B           &   !I [nbr] total number of molecules
            ,NFUNC_B          &   !I [nbr] total number of functional groups (e.g. CH2, NO3 ..)
            )

       !We are really only interested in gama for the partitioning species
       !Remember that the first 5 are POA species
       GAMA_B(:,:)=GAMA(:,NBSPAO_PUN+1:NBSPAO_PUN+NBSPB)

       !/* calculate partition coefficient */
       Call KPART_PUN_GET(          &
            TEMPK               & !I [K] temperature
            , GAMA_B            & !I [-] Activity coefficient
            , MWOM              & !I [g/mol] Molecular weight of organic matter
            , KPART             & !O [m3/ug] partitioning coefficient for species
            )

       !Get new aerosol concentrations
       DO I=1,NBSPB
          AERO_NEW(:,I)=        &
               KPART(:,I)       &   !m3/ug
               *GAS_OLD(:,I)    &   !ug/m3
               *M0_OLD(:)           !ug/m3  ==> ug/m3
       ENDDO

       !At this point we have only gas_old (from last iteration) and aero_new
       !We might now be in the situation that AERO_NEW + GAS_OLD .NE.CB
       !correct for it!
       DO I=1,NBSPB
          CB_NEW(:,I) = AERO_NEW(:,I) + GAS_OLD(:,I)  !New mass balance, not equal to CB
          AERO_NEW(:,I)=AERO_NEW(:,I)*CB(:,I)/CB_NEW(:,I)
          GAS_NEW1(:,I)=GAS_OLD(:,I)*CB(:,I)/CB_NEW(:,I)
          !write(6,*)"GASFRCB0", I, GAS_NEW1(:,I)/CB(:,I), (GAS_NEW1(:,I)+AERO_NEW(:,I))/CB(:,I)
       ENDDO

       !Get the new, total aerosol mass
       M0_NEW(:)=PAOM(:)
       DO I=1,NBSPB
          M0_NEW(:) = M0_NEW(:) + AERO_NEW(:,I)
       ENDDO

       !Get the gas phase concentration which would be in equilibrium with what we just found
       DO I=1,NBSPB
          GAS_NEW2(:,I)= AERO_NEW(:,I)/M0_NEW(:)/KPART(:,I)
       ENDDO

       !We might now be in the situation that AERO_NEW + GAS .NE.CB
       !correct for it!
       DO I=1,NBSPB
          CB_NEW(:,I) = AERO_NEW(:,I) + GAS_NEW2(:,I)  !New mass balance, not equal to CB
          AERO_NEW(:,I)=AERO_NEW(:,I)*CB(:,I)/CB_NEW(:,I)
          GAS_NEW2(:,I)=GAS_NEW2(:,I)*CB(:,I)/CB_NEW(:,I)
       ENDDO

       !Use the new value of AERO_NEW to make a new guess for AERO_OLD
       !if the sum of aero+gas is larger than the total, it means that two
       !much is in the aerosol phase, and we guess a concentration which
       !is lower than what we just calculated

       IF(LPRINT)THEN
          write(6,*)"                 "
          write(6,*)"*****************"
          write(6,'(a10,5e14.6)')"KPART",KPART
          write(6,'(a10,5e14.6)')"GAMA_B", GAMA_B
          write(6,'(a10,5e14.6)')"AERO_OLD",AERO_OLD
          write(6,'(a10,5e14.6)')"AERO_NEW",AERO_NEW
          write(6,'(a10,5e14.6)')"GAS_OLD",GAS_OLD
          write(6,'(a10,5e14.6)')"GAS_NEW",GAS_NEW1
          write(6,'(a10,5e14.6)')"GAS_NEW",GAS_NEW2
          write(6,'(a10,5e14.6)')"XPOA",X(:,1:5)
          write(6,'(a10,5e14.6)')"XSOA",X(:,6:10)
          write(6,*)"M0_OLD",M0_OLD, PAOM, M0_OLD-PAOM
       ENDIF

       !New guess for aerosol concentration
       !AERO_OLD(:,:) = AERO_NEW(:,:)

       !Our new guess for gas phase
       !GAS_OLD(:,:)=0.5*(GAS_NEW1(:,:)+GAS_NEW2(:,:))
       GAS_OLD(:,:)=sqrt(GAS_NEW1(:,:)*GAS_NEW2(:,:))

       !New guess for aerosol concentration
       AERO_OLD(:,:)=CB(:,:)-GAS_OLD(:,:)

       !Make a new guess for M0
       M0_OLD(:)=PAOM(:)
       DO I=1,NBSPB
          M0_OLD(:)  =  &
               M0_OLD(:) + AERO_OLD(:,I)
       ENDDO

    ENDDO  !Loop on iterations

    IF(LBOX)write(6,*)"sum X",sum(X)

    !Allocate memory for X, GAMA, ORGP
    DEALLOCATE(X)
    DEALLOCATE(GAMA)
    DEALLOCATE(GAMA_B)
    DEALLOCATE(ORGP)
    DEALLOCATE(SI_B)

    AERO(:,:)=AERO_OLD(:,:)
    GAS(:,:)=GAS_OLD(:,:)

    IF(LBOX)THEN
       write(6,*)"************************END VALUES********"
       write(6,'(a10,5e14.5)')"AERO",AERO
       write(6,'(a10,5e14.5)')"GAS",GAS
       write(6,'(a10,5e14.5)')"CB",CB
       write(6,'(a10,5e14.5)')"DIFF",CB-GAS-AERO
       write(6,'(a10,e14.5)')"PAOM",PAOM
       write(6,*)"********************************************"
       write(6,*)"  "
    ENDIF

  IF (LHOOK) CALL DR_HOOK('MODE_TYPEB:TYPEB',1,ZHOOK_HANDLE)
  END SUBROUTINE TYPEB

  !************************************************
   subroutine kpart_pun_get(                    &
       TEMPK                               & !I [K] temperature
       ,GAMA                               & !I [-] activity coefficient
       ,MWOM                               & !I [g/mol] average molar weight of organics
       ,KB                                 & !O [m3/ug] partitioning coefficient
       )

    implicit none

    !INPUT
    REAL, DIMENSION(:), INTENT(IN)         :: TEMPK !I [K] temperature
    REAL, DIMENSION(:,:), INTENT(IN)       :: GAMA  !I [-] activity coefficients
    REAL, DIMENSION(:)                     :: MWOM  !I [g/mol] molecular weight organics

    !OUTPUT
    REAL, DIMENSION(:,:), INTENT(OUT)      :: KB    !O [m3/ug] partitioning coefficient

    !LOCAL
    INTEGER                                :: I

    !The only goal of this routine is to get the partitioning coeffient:
    !So here is the way to do it:
    REAL(KIND=JPRB) :: ZHOOK_HANDLE
    IF (LHOOK) CALL DR_HOOK('MODE_TYPEB:KPART_PUN_GET',0,ZHOOK_HANDLE)
    do i=1,NBSPB
       KB(:,i)= 760.d0   &!torr/atm
            *R_UNIV      &!(m3*atm)/(mol*K)
            *1.d-6       &!g/ug
            *TEMPK(:)    &!K
            /VPB(i)      &!torr
            /MWOM(:)     &!g/mol
            /GAMA(:,i)      !==> m3/ug
    enddo

  IF (LHOOK) CALL DR_HOOK('MODE_TYPEB:KPART_PUN_GET',1,ZHOOK_HANDLE)
  end subroutine kpart_pun_get
  !************************************************************


END MODULE mode_typeb
