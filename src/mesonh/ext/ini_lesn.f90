!MNH_LIC Copyright 2000-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ####################
      SUBROUTINE  INI_LES_n
!     ####################
!
!
!!****  *INI_LES_n* initializes the LES variables for model _n
!!
!!    PURPOSE
!!    -------
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         07/02/00
!!      Modification     01/02/01 (D.Gazen) add module MODD_NSV for NSV variable
!!                       06/11/02 (V. Masson) add LES budgets
!!                    10/2016 (C.Lac) Add droplet deposition
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!                     02/2019 (C. Lac) Add rain fraction as a LES diagnostic
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 13/09/2019: budget: simplify and modernize date/time management
!  P. Wautelet 12/08/2020: bugfix: use NUNDEF instead of XUNDEF for integer variables
!  P. Wautelet 04/01/2021: bugfix: nles_k was used instead of nspectra_k for a loop index
!  P. Wautelet 30/03/2021: budgets: LES cartesian subdomain limits are defined in the physical domain
!  P. Wautelet 09/07/2021: bugfix: altitude levels are on the correct grid position (mass point)
!  P. Wautelet 22/03/2022: LES averaging periods are more reliable (compute with integers instead of reals)
! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
USE MODE_MSG
USE MODE_MODELN_HANDLER
!
USE MODD_LES
USE MODD_LES_BUDGET
USE MODD_LES_n
!
USE MODD_CONF
USE MODD_PARAMETERS
USE MODD_NESTING
!
USE MODD_LUNIT_n
USE MODD_GRID_n
USE MODD_DYN_n
USE MODD_TIME_n
USE MODD_DIM_n
USE MODD_TURB_n
USE MODD_CONF_n
USE MODD_LBC_n
USE MODD_PARAM_n
USE MODD_DYN
USE MODD_NSV, ONLY: NSV ! update_nsv is done in INI_MODEL
USE MODD_CONDSAMP, ONLY : LCONDSAMP
!
USE MODI_INI_LES_CART_MASKn
USE MODI_COEF_VER_INTERP_LIN
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
!
!
!
!       0.2  declaration of local variables
!
!
!      
INTEGER :: ILUOUT, IRESP
INTEGER :: JI,JJ, JK ! loop counters
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZZ_LES ! LES altitudes 3D array
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZZ_SPEC! " for spectra
!
!
REAL, DIMENSION(:), POINTER :: ZXHAT_ll ! father model coordinates
REAL, DIMENSION(:), POINTER :: ZYHAT_ll !
INTEGER :: IMI
!
!-------------------------------------------------------------------------------
IMI = GET_CURRENT_MODEL_INDEX()
!
ZXHAT_ll => NULL()
ZYHAT_ll => NULL()
!
ILUOUT = TLUOUT%NLU
!
!-------------------------------------------------------------------------------
!
!*      1.   Does LES computations are used?
!            ------------------------------
!
LLES = LLES_MEAN .OR. LLES_RESOLVED  .OR. LLES_SUBGRID .OR. LLES_UPDRAFT &
                 .OR. LLES_DOWNDRAFT .OR. LLES_SPECTRA
!
!
IF (.NOT. LLES) RETURN
!
IF (L1D) THEN
  LLES_RESOLVED    = .FALSE.
  LLES_UPDRAFT     = .FALSE.
  LLES_DOWNDRAFT   = .FALSE.
  LLES_SPECTRA     = .FALSE.
  LLES_NEB_MASK    = .FALSE.
  LLES_CORE_MASK   = .FALSE.
  LLES_CS_MASK     = .FALSE.
  LLES_MY_MASK     = .FALSE.
END IF
!
IF (LLES_RESOLVED )  LLES_MEAN = .TRUE.
IF (LLES_SUBGRID  )  LLES_MEAN = .TRUE.
IF (LLES_UPDRAFT  )  LLES_MEAN = .TRUE.
IF (LLES_DOWNDRAFT)  LLES_MEAN = .TRUE.
IF (LLES_SPECTRA  )  LLES_MEAN = .TRUE.
!
IF (CTURB=='NONE') THEN
  WRITE(ILUOUT,FMT=*) 'LES diagnostics cannot be done without subgrid turbulence.'
  WRITE(ILUOUT,FMT=*) 'You have chosen CTURB="NONE". You must choose a turbulence scheme.'
  call Print_msg( NVERB_FATAL, 'GEN', 'WRITE_LB_n', 'LES diagnostics cannot be done without subgrid turbulence' )
END IF
!-------------------------------------------------------------------------------
!
!*      2.   Number and definition of masks
!            ------------------------------
!
!-------------------------------------------------------------------------------
!
!*      2.1  Cartesian (sub-)domain
!            ----------------------
!
!* updates number of masks
!  -----------------------
!
NLES_MASKS = 1
!
!* For model 1, set default values of cartesian mask, and defines cartesian mask
!  -----------------------------------------------------------------------------
!
IF (IMI==1) THEN
  IF ( LLES_CART_MASK ) THEN
    !Compute LES diagnostics inside a cartesian mask

    !Set default values to physical domain boundaries
    IF ( NLES_IINF == NUNDEF ) NLES_IINF = 1
    IF ( NLES_JINF == NUNDEF ) NLES_JINF = 1
    IF ( NLES_ISUP == NUNDEF ) NLES_ISUP = NIMAX_ll
    IF ( NLES_JSUP == NUNDEF ) NLES_JSUP = NJMAX_ll

    !Check that selected indices are in physical domain
    IF ( NLES_IINF < 1 )         CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_IINF too small (<1)' )
    IF ( NLES_IINF > NIMAX_ll )  CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_IINF too large (>NIMAX)' )
    IF ( NLES_ISUP < 1 )         CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_ISUP too small (<1)' )
    IF ( NLES_ISUP > NIMAX_ll )  CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_ISUP too large (>NIMAX)' )
    IF ( NLES_ISUP < NLES_IINF ) CALL Print_msg( NVERB_ERROR, 'BUD', 'INI_LES_n', 'NLES_ISUP < NLES_IINF' )

    IF ( NLES_JINF < 1 )         CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_JINF too small (<1)' )
    IF ( NLES_JINF > NJMAX_ll )  CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_JINF too large (>NJMAX)' )
    IF ( NLES_JSUP < 1 )         CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_JSUP too small (<1)' )
    IF ( NLES_JSUP > NJMAX_ll )  CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_JSUP too large (>NJMAX)' )
    IF ( NLES_JSUP < NLES_JINF ) CALL Print_msg( NVERB_ERROR, 'BUD', 'INI_LES_n', 'NLES_JSUP < NLES_JINF' )

    !Set LLES_CART_MASK to false if whole domain is selected
    IF (       NLES_IINF == 1        .AND. NLES_JINF == 1        &
         .AND. NLES_ISUP == NIMAX_ll .AND. NLES_ISUP == NJMAX_ll ) THEN
      LLES_CART_MASK = .FALSE.
    END IF
  ELSE
    !Compute LES diagnostics on whole physical domain
    NLES_IINF = 1
    NLES_JINF = 1
    NLES_ISUP = NIMAX_ll
    NLES_JSUP = NJMAX_ll
  END IF
  !
  NLESn_IINF(1)= NLES_IINF
  NLESn_ISUP(1)= NLES_ISUP
  NLESn_JINF(1)= NLES_JINF
  NLESn_JSUP(1)= NLES_JSUP
!
!* For other models, fits cartesian mask on model 1 mask
!  -----------------------------------------------------
!
ELSE
  ZXHAT_ll => XXHAT_ll !Use current (IMI) model XXHAT_ll
  ZYHAT_ll => XYHAT_ll
!
  CALL GOTO_MODEL(NDAD(IMI))
  CALL INI_LES_CART_MASK_n(IMI,ZXHAT_ll,ZYHAT_ll,          &
                         NLESn_IINF(IMI),NLESn_JINF(IMI), &
                         NLESn_ISUP(IMI),NLESn_JSUP(IMI)  )
  CALL GOTO_MODEL(IMI)
END IF
!
!* in non cyclic boundary conditions, limitiation of masks due to u and v grids
!  ----------------------------------------------------------------------------
!
IF ( (.NOT. L1D) .AND. CLBCX(1)/='CYCL') THEN
  NLESn_IINF(IMI) = MAX(NLESn_IINF(IMI),2)
END IF
IF ( (.NOT. L1D) .AND. (.NOT. L2D) .AND. CLBCY(1)/='CYCL') THEN
  NLESn_JINF(IMI) = MAX(NLESn_JINF(IMI),2)
END IF
!
!* X boundary conditions for 2points correlations computations
!  -----------------------------------------------------------
!
IF ( CLBCX(1) == 'CYCL' .AND. NLESn_IINF(IMI) == 1 .AND. NLESn_ISUP(IMI) == NIMAX_ll ) THEN
  CLES_LBCX(:,IMI) = 'CYCL'
ELSE
  CLES_LBCX(:,IMI) = 'OPEN'
END IF
!
!* Y boundary conditions for 2points correlations computations
!  -----------------------------------------------------------
!
IF ( CLBCY(1) == 'CYCL' .AND. NLESn_JINF(IMI) == 1 .AND. NLESn_JSUP(IMI) == NJMAX_ll ) THEN
  CLES_LBCY(:,IMI) = 'CYCL'
ELSE
  CLES_LBCY(:,IMI) = 'OPEN'
END IF
!
!-------------------------------------------------------------------------------
!
!*      2.2  Nebulosity mask
!            ---------------
!
IF (.NOT. LUSERC .AND. .NOT. LUSERI) LLES_NEB_MASK = .FALSE.
!
IF (LLES_NEB_MASK) NLES_MASKS = NLES_MASKS + 2
!
!-------------------------------------------------------------------------------
!
!*      2.3  Cloud core mask
!            ---------------
!
IF (.NOT. LUSERC .AND. .NOT. LUSERI) LLES_CORE_MASK = .FALSE.
!
IF (LLES_CORE_MASK) NLES_MASKS = NLES_MASKS + 2
!
!-------------------------------------------------------------------------------
!
!*      2.4  Conditional sampling mask
!            -------------------------
!
IF (.NOT. LUSERC .AND. .NOT. LCONDSAMP) LLES_CS_MASK = .FALSE.
!
IF (LLES_CS_MASK) NLES_MASKS = NLES_MASKS + 3
!
!-------------------------------------------------------------------------------
!
!*      2.5  User mask
!            ---------
!
IF (LLES_MY_MASK) NLES_MASKS = NLES_MASKS + NLES_MASKS_USER
!
!-------------------------------------------------------------------------------
!
!*      3.   Number of temporal LES samplings
!            --------------------------------
!
!*      3.1  Default value
!            -------------
!
IF (XLES_TEMP_SAMPLING == XUNDEF) THEN
  IF (CTURBDIM=='3DIM') THEN
    XLES_TEMP_SAMPLING =  60.
  ELSE
    XLES_TEMP_SAMPLING = 300.
  END IF
END IF
!
!*      3.2  Number of time steps between two calls
!            --------------------------------------
!
NLES_DTCOUNT = MAX( NINT( XLES_TEMP_SAMPLING / XTSTEP ) , 1)

!
!*      3.3  Redefinition of the LES sampling time coherent with model time-step
!            -------------------------------------------------------------------
!
! Note that this modifies XLES_TEMP_SAMPLING only for father model (model number 1)
! For nested models (for which integration time step is an integer part of  father model)
! the following operation does not change XLES_TEMP_SAMPLING. This way, LEs
! sampling is done at the same instants for all models.
!
XLES_TEMP_SAMPLING = XTSTEP * NLES_DTCOUNT
!
!
!*      3.4  number of temporal calls to LES routines
!            ----------------------------------------
!
!
NLES_TIMES = ( NINT( ( XSEGLEN - DYN_MODEL(1)%XTSTEP ) / XTSTEP ) ) / NLES_DTCOUNT
!
!*      3.5  current LES time counter
!            ------------------------
!
NLES_TCOUNT = 0
!
!*      3.6  dates array for diachro
!            ----------------------
!
allocate( tles_dates( nles_times ) )
allocate( xles_times( nles_times ) )
!
!*      3.7  No data
!            -------
!
IF (NLES_TIMES==0) THEN
  LLES=.FALSE.
  RETURN
END IF
!
!*     3.8  Averaging
!           ---------
IF (     XLES_TEMP_MEAN_END   == XUNDEF &
    .OR. XLES_TEMP_MEAN_START == XUNDEF &
    .OR. XLES_TEMP_MEAN_STEP  == XUNDEF ) THEN
  !No LES temporal averaging
  NLES_MEAN_TIMES = 0
  NLES_MEAN_STEP  = NNEGUNDEF
  NLES_MEAN_START = NNEGUNDEF
  NLES_MEAN_END   = NNEGUNDEF
ELSE
  !LES temporal averaging is enabled
  !Ensure that XLES_TEMP_MEAN_END is not after segment end
  XLES_TEMP_MEAN_END = MIN( XLES_TEMP_MEAN_END, XSEGLEN - DYN_MODEL(1)%XTSTEP )

  NLES_MEAN_START = NINT( XLES_TEMP_MEAN_START / XTSTEP )

  IF ( MODULO( NLES_MEAN_START, NLES_DTCOUNT ) /= 0 ) THEN
    CMNHMSG(1) = 'XLES_TEMP_MEAN_START is not a multiple of XLES_TEMP_SAMPLING'
    CMNHMSG(2) = 'LES averaging periods could be wrong'
    CALL Print_msg( NVERB_WARNING, 'IO', 'INI_LES_n' )
  END IF

  NLES_MEAN_END = NINT( XLES_TEMP_MEAN_END / XTSTEP )

  NLES_MEAN_STEP = NINT( XLES_TEMP_MEAN_STEP / XTSTEP )

  IF ( NLES_MEAN_STEP < NLES_DTCOUNT ) &
    CALL Print_msg( NVERB_ERROR, 'IO', 'INI_LES_n', 'XLES_TEMP_MEAN_STEP < XLES_TEMP_SAMPLING not allowed' )

  IF ( MODULO( NLES_MEAN_STEP, NLES_DTCOUNT ) /= 0 ) THEN
    CMNHMSG(1) = 'XLES_TEMP_MEAN_STEP is not a multiple of XLES_TEMP_SAMPLING'
    CMNHMSG(2) = 'LES averaging periods could be wrong'
    CALL Print_msg( NVERB_WARNING, 'IO', 'INI_LES_n' )
  END IF

  NLES_MEAN_TIMES = ( NLES_MEAN_END - NLES_MEAN_START ) / NLES_MEAN_STEP
  !Add 1 averaging period if the last one is incomplete (for example: start=0., end=10., step=3.)
  IF ( MODULO( NLES_MEAN_END - NLES_MEAN_START, NLES_MEAN_STEP ) > 0 ) NLES_MEAN_TIMES = NLES_MEAN_TIMES + 1
END IF
!-------------------------------------------------------------------------------
!
!*      4.   Number of vertical levels for local diagnostics
!            -----------------------------------------------
!
NLES_K = 0
!
!*      4.1  Case of altitude levels (lowest priority)
!            -----------------------
!
IF (ANY(XLES_ALTITUDES(:)/=XUNDEF)) THEN
  NLES_K = COUNT (XLES_ALTITUDES(:)/=XUNDEF)
  CLES_LEVEL_TYPE='Z'
  !
  ALLOCATE(XCOEFLIN_LES(SIZE(XZZ,1),SIZE(XZZ,2),NLES_K))
  ALLOCATE(NKLIN_LES   (SIZE(XZZ,1),SIZE(XZZ,2),NLES_K))
  !
  ALLOCATE(ZZ_LES      (SIZE(XZZ,1),SIZE(XZZ,2),NLES_K))
  DO JK=1,NLES_K
    DO JJ=1,SIZE(XZZ,2)
      DO JI=1,SIZE(XZZ,1)
        ZZ_LES(JI,JJ,JK) = XLES_ALTITUDES(JK)
      END DO
    END DO
  END DO
  CALL COEF_VER_INTERP_LIN(MZF(XZZ),ZZ_LES,NKLIN_LES,XCOEFLIN_LES)
  !
  DEALLOCATE(ZZ_LES)
END IF
!
!
!*      4.2  Case of model levels (highest priority)
!            --------------------
!
IF (ANY(NLES_LEVELS(:)/=NUNDEF)) THEN
  DO JK = 1, SIZE( NLES_LEVELS )
    IF ( NLES_LEVELS(JK) /= NUNDEF ) THEN
      IF ( NLES_LEVELS(JK) < 1 )     CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_LEVELS too small (<1)' )
      IF ( NLES_LEVELS(JK) > NKMAX ) CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NLES_LEVELS too large (>NKMAX)' )
    END IF
  END DO

  NLES_K = COUNT (NLES_LEVELS(:)/=NUNDEF)
  CLES_LEVEL_TYPE='K'
ELSE
  IF (NLES_K==0) THEN
    NLES_K = MIN(SIZE(NLES_LEVELS),NKMAX)
    CLES_LEVEL_TYPE='K'
    DO JK=1,NLES_K
      NLES_LEVELS(JK) = JK
    END DO
  END IF
END IF
!
!-------------------------------------------------------------------------------
!
!*      5.   Number of vertical levels for non-local diagnostics
!            ---------------------------------------------------
!
NSPECTRA_K = 0
CSPECTRA_LEVEL_TYPE='N'
!
!
!*      5.1  Case of altitude levels (medium priority)
!            -----------------------
!
IF (ANY(XSPECTRA_ALTITUDES(:)/=XUNDEF)) THEN
  NSPECTRA_K = COUNT (XSPECTRA_ALTITUDES(:)/=XUNDEF)
  CSPECTRA_LEVEL_TYPE='Z'
  !
  ALLOCATE(XCOEFLIN_SPEC(SIZE(XZZ,1),SIZE(XZZ,2),NSPECTRA_K))
  ALLOCATE(NKLIN_SPEC   (SIZE(XZZ,1),SIZE(XZZ,2),NSPECTRA_K))
  !
  ALLOCATE(ZZ_SPEC      (SIZE(XZZ,1),SIZE(XZZ,2),NSPECTRA_K))
  DO JK=1,NSPECTRA_K
    DO JJ=1,SIZE(XZZ,2)
      DO JI=1,SIZE(XZZ,1)
        ZZ_SPEC(JI,JJ,JK) = XSPECTRA_ALTITUDES(JK)
      END DO
    END DO
  END DO
  CALL COEF_VER_INTERP_LIN(XZZ,ZZ_SPEC,NKLIN_SPEC,XCOEFLIN_SPEC)
  !
  DEALLOCATE(ZZ_SPEC)
END IF
!
!
!*      5.2  Case of model levels (highest priority)
!            --------------------
!
IF (ANY(NSPECTRA_LEVELS(:)/=NUNDEF)) THEN
  DO JK = 1, SIZE( NSPECTRA_LEVELS )
    IF ( NSPECTRA_LEVELS(JK) /= NUNDEF ) THEN
      IF ( NSPECTRA_LEVELS(JK) < 1 )     CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NSPECTRA_LEVELS too small (<1)' )
      IF ( NSPECTRA_LEVELS(JK) > NKMAX ) CALL Print_msg( NVERB_ERROR, 'GEN', 'INI_LES_n', 'NSPECTRA_LEVELS too large (>NKMAX)' )
    END IF
  END DO

  NSPECTRA_K = COUNT (NSPECTRA_LEVELS(:)/=NUNDEF)
  CSPECTRA_LEVEL_TYPE='K'
END IF
!
!-------------------------------------------------------------------------------
!
!*      6.   Number of horizontal wavelengths for non-local diagnostics
!            ----------------------------------------------------------
!
NSPECTRA_NI = NLESn_ISUP(IMI) - NLESn_IINF(IMI) + 1
NSPECTRA_NJ = NLESn_JSUP(IMI) - NLESn_JINF(IMI) + 1
!
!
!-------------------------------------------------------------------------------
!
!*      7.   Allocations of temporal series of local diagnostics
!            ---------------------------------------------------
!
!*      7.0  Altitude levels
!            ---------------
!
ALLOCATE(XLES_Z  (NLES_K))
!
!*      7.1  Averaging control variables
!            ---------------------------
!
ALLOCATE(NLES_AVG_PTS_ll (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(NLES_UND_PTS_ll (NLES_K,NLES_TIMES,NLES_MASKS))
!
NLES_AVG_PTS_ll = NUNDEF
NLES_UND_PTS_ll = NUNDEF
!
!
!*      7.2  Horizontally mean variables
!            ---------------------------
!
ALLOCATE(XLES_MEAN_U  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_V  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_W  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_P  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_DP  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_TP  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_TR  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_DISS(NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_LM  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_RHO(NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_Th (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_Mf (NLES_K,NLES_TIMES,NLES_MASKS))  
IF (LUSERC ) THEN
  ALLOCATE(XLES_MEAN_Thl(NLES_K,NLES_TIMES,NLES_MASKS))
  ALLOCATE(XLES_MEAN_Rt (NLES_K,NLES_TIMES,NLES_MASKS))
  ALLOCATE(XLES_MEAN_KHt(NLES_K,NLES_TIMES,NLES_MASKS))
  ALLOCATE(XLES_MEAN_KHr(NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Thl(0,0,0))
  ALLOCATE(XLES_MEAN_Rt (0,0,0))
  ALLOCATE(XLES_MEAN_KHt(0,0,0))
  ALLOCATE(XLES_MEAN_KHr(0,0,0))
END IF
IF (LUSERV) THEN
  ALLOCATE(XLES_MEAN_Thv(NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Thv(0,0,0))
END IF
!
IF (LUSERV ) THEN 
  ALLOCATE(XLES_MEAN_Rv (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Rv (0,0,0))
END IF
IF (LUSERV ) THEN 
   ALLOCATE(XLES_MEAN_Rehu (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Rehu (0,0,0))
ENDIF
IF (LUSERV ) THEN 
   ALLOCATE(XLES_MEAN_Qs (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Qs (0,0,0))  
END IF
IF (LUSERC ) THEN
  ALLOCATE(XLES_MEAN_Rc (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Rc (0,0,0))
END IF
IF (LUSERC ) THEN
  ALLOCATE(XLES_MEAN_Cf (NLES_K,NLES_TIMES,NLES_MASKS))
  ALLOCATE(XLES_MEAN_INDCf (NLES_K,NLES_TIMES,NLES_MASKS))
  ALLOCATE(XLES_MEAN_INDCf2 (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Cf (0,0,0))
  ALLOCATE(XLES_MEAN_INDCf (0,0,0))
  ALLOCATE(XLES_MEAN_INDCf2(0,0,0))  
END IF
IF (LUSERR ) THEN
  ALLOCATE(XLES_MEAN_Rr (NLES_K,NLES_TIMES,NLES_MASKS))
  ALLOCATE(XLES_MEAN_RF (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Rr (0,0,0))
  ALLOCATE(XLES_MEAN_RF (0,0,0))
END IF
IF (LUSERI ) THEN
  ALLOCATE(XLES_MEAN_Ri (NLES_K,NLES_TIMES,NLES_MASKS))
  ALLOCATE(XLES_MEAN_If (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Ri (0,0,0))
  ALLOCATE(XLES_MEAN_If (0,0,0))
END IF
IF (LUSERS ) THEN
  ALLOCATE(XLES_MEAN_Rs (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Rs (0,0,0))
END IF
IF (LUSERG ) THEN
  ALLOCATE(XLES_MEAN_Rg (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Rg (0,0,0))
END IF
IF (LUSERH ) THEN
  ALLOCATE(XLES_MEAN_Rh (NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_Rh (0,0,0))
END IF
IF (NSV>0  ) THEN
  ALLOCATE(XLES_MEAN_Sv (NLES_K,NLES_TIMES,NLES_MASKS,NSV))
ELSE
  ALLOCATE(XLES_MEAN_Sv (0,0,0,0))
END IF
ALLOCATE(XLES_MEAN_WIND  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_dUdz  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_dVdz  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_dWdz  (NLES_K,NLES_TIMES,NLES_MASKS))
ALLOCATE(XLES_MEAN_dThldz(NLES_K,NLES_TIMES,NLES_MASKS))
IF (LUSERV) THEN
  ALLOCATE(XLES_MEAN_dRtdz(NLES_K,NLES_TIMES,NLES_MASKS))
ELSE
  ALLOCATE(XLES_MEAN_dRtdz(0,0,0))
END IF
IF (NSV>0) THEN
  ALLOCATE(XLES_MEAN_dSvdz(NLES_K,NLES_TIMES,NLES_MASKS,NSV))
ELSE
  ALLOCATE(XLES_MEAN_dSvdz(0,0,0,0))
END IF
!
IF (LLES_PDF) THEN
!pdf distributions and jpdf distributions
  CALL LES_ALLOCATE('XLES_PDF_TH ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  CALL LES_ALLOCATE('XLES_PDF_W ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  CALL LES_ALLOCATE('XLES_PDF_THV ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  IF (LUSERV) THEN
   CALL LES_ALLOCATE('XLES_PDF_RV ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  ELSE
   CALL LES_ALLOCATE('XLES_PDF_RV ',(/0,0,0,0/))
  END IF
  IF (LUSERC) THEN 
   CALL LES_ALLOCATE('XLES_PDF_RC ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
   CALL LES_ALLOCATE('XLES_PDF_RT ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
   CALL LES_ALLOCATE('XLES_PDF_THL',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  ELSE
   CALL LES_ALLOCATE('XLES_PDF_RC ',(/0,0,0,0/))
   CALL LES_ALLOCATE('XLES_PDF_RT ',(/0,0,0,0/))
   CALL LES_ALLOCATE('XLES_PDF_THL',(/0,0,0,0/))
  ENDIF
  IF (LUSERR) THEN
   CALL LES_ALLOCATE('XLES_PDF_RR ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  ELSE
   CALL LES_ALLOCATE('XLES_PDF_RR ',(/0,0,0,0/))
  ENDIF
  IF (LUSERI) THEN 
   CALL LES_ALLOCATE('XLES_PDF_RI ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  ELSE
   CALL LES_ALLOCATE('XLES_PDF_RI ',(/0,0,0,0/))
  END IF
  IF (LUSERS) THEN
   CALL LES_ALLOCATE('XLES_PDF_RS ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  ELSE
   CALL LES_ALLOCATE('XLES_PDF_RS ',(/0,0,0,0/))
  END IF
  IF (LUSERG) THEN 
   CALL LES_ALLOCATE('XLES_PDF_RG ',(/NLES_K,NLES_TIMES,NLES_MASKS,NPDF/))
  ELSE
   CALL LES_ALLOCATE('XLES_PDF_RG ',(/0,0,0,0/))
  END IF
ENDIF
!  
XLES_MEAN_U  = XUNDEF
XLES_MEAN_V  = XUNDEF
XLES_MEAN_W  = XUNDEF
XLES_MEAN_P  = XUNDEF
XLES_MEAN_DP  = XUNDEF
XLES_MEAN_TP  = XUNDEF
XLES_MEAN_TR  = XUNDEF
XLES_MEAN_DISS= XUNDEF
XLES_MEAN_LM  = XUNDEF
XLES_MEAN_RHO= XUNDEF
XLES_MEAN_Th = XUNDEF
XLES_MEAN_Mf = XUNDEF
IF (LUSERC ) XLES_MEAN_Thl= XUNDEF
IF (LUSERV ) XLES_MEAN_Thv= XUNDEF
IF (LUSERV ) XLES_MEAN_Rv  = XUNDEF
IF (LUSERV ) XLES_MEAN_Rehu  = XUNDEF  
IF (LUSERV ) XLES_MEAN_Qs  = XUNDEF
IF (LUSERC ) XLES_MEAN_KHr = XUNDEF
IF (LUSERC ) XLES_MEAN_KHt = XUNDEF
IF (LUSERC ) XLES_MEAN_Rt  = XUNDEF
IF (LUSERC ) XLES_MEAN_Rc  = XUNDEF
IF (LUSERC ) XLES_MEAN_Cf  = XUNDEF
IF (LUSERC ) XLES_MEAN_RF  = XUNDEF
IF (LUSERC ) XLES_MEAN_INDCf  = XUNDEF
IF (LUSERC ) XLES_MEAN_INDCf2 = XUNDEF
IF (LUSERR ) XLES_MEAN_Rr  = XUNDEF
IF (LUSERI ) XLES_MEAN_Ri  = XUNDEF
IF (LUSERI ) XLES_MEAN_If  = XUNDEF
IF (LUSERS ) XLES_MEAN_Rs  = XUNDEF
IF (LUSERG ) XLES_MEAN_Rg  = XUNDEF
IF (LUSERH ) XLES_MEAN_Rh  = XUNDEF
IF (NSV>0  ) XLES_MEAN_Sv  = XUNDEF
XLES_MEAN_WIND  = XUNDEF
XLES_MEAN_WIND  = XUNDEF
XLES_MEAN_dUdz  = XUNDEF
XLES_MEAN_dVdz  = XUNDEF
XLES_MEAN_dWdz  = XUNDEF
XLES_MEAN_dThldz= XUNDEF
IF (LUSERV) XLES_MEAN_dRtdz = XUNDEF
IF (NSV>0)  XLES_MEAN_dSvdz = XUNDEF
!
IF (LLES_PDF) THEN
 XLES_PDF_TH   = XUNDEF
 XLES_PDF_W    = XUNDEF
 XLES_PDF_THV  = XUNDEF
 IF (LUSERV) THEN
  XLES_PDF_RV   = XUNDEF
 END IF
 IF (LUSERC) THEN
  XLES_PDF_RC    = XUNDEF
  XLES_PDF_RT    = XUNDEF
  XLES_PDF_THL   = XUNDEF
 END IF
 IF (LUSERR) THEN
  XLES_PDF_RR   = XUNDEF
 END IF
 IF (LUSERI) THEN
  XLES_PDF_RI   = XUNDEF
 END IF
 IF (LUSERS) THEN
  XLES_PDF_RS   = XUNDEF
 END IF
 IF (LUSERG) THEN
  XLES_PDF_RG   = XUNDEF
 END IF
END IF
!
!
!
!*      7.3  Resolved quantities
!            -------------------
!
ALLOCATE(XLES_RESOLVED_U2    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'2>
ALLOCATE(XLES_RESOLVED_V2    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'2>
ALLOCATE(XLES_RESOLVED_W2    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'2>
ALLOCATE(XLES_RESOLVED_P2    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <p'2>
ALLOCATE(XLES_RESOLVED_Th2   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Th'2>
IF (LUSERV) THEN
  ALLOCATE(XLES_RESOLVED_ThThv (NLES_K,NLES_TIMES,NLES_MASKS))   ! <Th'Thv'>
ELSE
  ALLOCATE(XLES_RESOLVED_ThThv (0,0,0))
END IF
IF (LUSERC) THEN
  ALLOCATE(XLES_RESOLVED_Thl2  (NLES_K,NLES_TIMES,NLES_MASKS))   ! <Thl'2>
  ALLOCATE(XLES_RESOLVED_ThlThv(NLES_K,NLES_TIMES,NLES_MASKS))   ! <Thl'Thv'>
ELSE
  ALLOCATE(XLES_RESOLVED_Thl2  (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThlThv(0,0,0))
END IF
ALLOCATE(XLES_RESOLVED_Ke    (NLES_K,NLES_TIMES,NLES_MASKS))     ! 0.5 <u'2+v'2+w'2>
ALLOCATE(XLES_RESOLVED_UV    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'v'>
ALLOCATE(XLES_RESOLVED_WU    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'u'>
ALLOCATE(XLES_RESOLVED_WV    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'v'>
ALLOCATE(XLES_RESOLVED_UP    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'p'>
ALLOCATE(XLES_RESOLVED_VP    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'p'>
ALLOCATE(XLES_RESOLVED_WP    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'p'>
ALLOCATE(XLES_RESOLVED_UTh   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'Th'>
ALLOCATE(XLES_RESOLVED_VTh   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'Th'>
ALLOCATE(XLES_RESOLVED_WTh   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Th'>
IF (LUSERC) THEN
  ALLOCATE(XLES_RESOLVED_UThl  (NLES_K,NLES_TIMES,NLES_MASKS))   ! <u'Thl'>
  ALLOCATE(XLES_RESOLVED_VThl  (NLES_K,NLES_TIMES,NLES_MASKS))   ! <v'Thl'>
  ALLOCATE(XLES_RESOLVED_WThl  (NLES_K,NLES_TIMES,NLES_MASKS))   ! <w'Thl'>
ELSE
  ALLOCATE(XLES_RESOLVED_UThl(0,0,0))
  ALLOCATE(XLES_RESOLVED_VThl(0,0,0))
  ALLOCATE(XLES_RESOLVED_WThl(0,0,0))
END IF
IF (LUSERV) THEN
  ALLOCATE(XLES_RESOLVED_UThv  (NLES_K,NLES_TIMES,NLES_MASKS))   ! <u'Thv'>
  ALLOCATE(XLES_RESOLVED_VThv  (NLES_K,NLES_TIMES,NLES_MASKS))   ! <v'Thv'>
  ALLOCATE(XLES_RESOLVED_WThv  (NLES_K,NLES_TIMES,NLES_MASKS))   ! <w'Thv'>
ELSE
  ALLOCATE(XLES_RESOLVED_UThv(0,0,0))
  ALLOCATE(XLES_RESOLVED_VThv(0,0,0))
  ALLOCATE(XLES_RESOLVED_WThv(0,0,0))
END IF
ALLOCATE(XLES_RESOLVED_U3    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'3>
ALLOCATE(XLES_RESOLVED_V3    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'3>
ALLOCATE(XLES_RESOLVED_W3    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'3>
ALLOCATE(XLES_RESOLVED_U4    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'4>
ALLOCATE(XLES_RESOLVED_V4    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'4>
ALLOCATE(XLES_RESOLVED_W4    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'4>
ALLOCATE(XLES_RESOLVED_ThlPz (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thv'dp'/dz>
ALLOCATE(XLES_RESOLVED_WThl2 (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Thl'2>
ALLOCATE(XLES_RESOLVED_W2Thl (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'2Thl'>
ALLOCATE(XLES_RESOLVED_MASSFX(NLES_K,NLES_TIMES,NLES_MASKS))     ! <upward mass flux>
ALLOCATE(XLES_RESOLVED_UKe   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'(u'2+v'2+w'2)>
ALLOCATE(XLES_RESOLVED_VKe   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'(u'2+v'2+w'2)>
ALLOCATE(XLES_RESOLVED_WKe   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'(u'2+v'2+w'2)>

IF (LUSERV ) THEN
  ALLOCATE(XLES_RESOLVED_Rv2   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rv'2>
  ALLOCATE(XLES_RESOLVED_ThRv  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Th'Rv'>
  ALLOCATE(XLES_RESOLVED_ThvRv (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thv'Rv'>
  ALLOCATE(XLES_RESOLVED_URv   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'Rv'>
  ALLOCATE(XLES_RESOLVED_VRv   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'Rv'>
  ALLOCATE(XLES_RESOLVED_WRv   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rv'>
  ALLOCATE(XLES_RESOLVED_WRv2  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rv'2>
  ALLOCATE(XLES_RESOLVED_W2Rv  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'2Rv'>  
  ALLOCATE(XLES_RESOLVED_W2Rt  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'2Rt'> 
  ALLOCATE(XLES_RESOLVED_WRt2  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rt2'>    
  ALLOCATE(XLES_RESOLVED_RvPz  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rv'dp'/dz>
  ALLOCATE(XLES_RESOLVED_WThlRv(NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Thl'Rv'>
  ALLOCATE(XLES_RESOLVED_WThlRt(NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Thl'Rt'>  
ELSE
  ALLOCATE(XLES_RESOLVED_Rv2   (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThRv  (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThvRv (0,0,0))
  ALLOCATE(XLES_RESOLVED_URv   (0,0,0))
  ALLOCATE(XLES_RESOLVED_VRv   (0,0,0))
  ALLOCATE(XLES_RESOLVED_WRv   (0,0,0))
  ALLOCATE(XLES_RESOLVED_WRv2  (0,0,0))
  ALLOCATE(XLES_RESOLVED_W2Rv  (0,0,0))
  ALLOCATE(XLES_RESOLVED_W2Rt  (0,0,0)) 
  ALLOCATE(XLES_RESOLVED_WRt2  (0,0,0))
  ALLOCATE(XLES_RESOLVED_RvPz  (0,0,0))
  ALLOCATE(XLES_RESOLVED_WThlRv(0,0,0))
  ALLOCATE(XLES_RESOLVED_WThlRt(0,0,0))
END IF
IF (LUSERC ) THEN
  ALLOCATE(XLES_RESOLVED_ThlRv (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thl'Rv'>
  !
  ALLOCATE(XLES_RESOLVED_Rc2   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rc'2>
  ALLOCATE(XLES_RESOLVED_ThRc  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Th'Rc'>
  ALLOCATE(XLES_RESOLVED_ThlRc (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thl'Rc'>
  ALLOCATE(XLES_RESOLVED_ThvRc (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thv'Rc'>
  ALLOCATE(XLES_RESOLVED_URc   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'Rc'>
  ALLOCATE(XLES_RESOLVED_VRc   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'Rc'>
  ALLOCATE(XLES_RESOLVED_WRc   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rc'>
  ALLOCATE(XLES_RESOLVED_WRc2  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rc'2>
  ALLOCATE(XLES_RESOLVED_W2Rc  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'2Rc'>
  ALLOCATE(XLES_RESOLVED_RcPz  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rc'dp'/dz>
  ALLOCATE(XLES_RESOLVED_WThlRc(NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Thl'Rc'>
  ALLOCATE(XLES_RESOLVED_WRvRc (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rv'Rc'>
  ALLOCATE(XLES_RESOLVED_WRt   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rt'>
  ALLOCATE(XLES_RESOLVED_Rt2   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rt'2>
  ALLOCATE(XLES_RESOLVED_RtPz  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rv'dp'/dz>
ELSE
  ALLOCATE(XLES_RESOLVED_ThlRv (0,0,0))
  !
  ALLOCATE(XLES_RESOLVED_Rc2   (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThRc  (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThlRc (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThvRc (0,0,0))
  ALLOCATE(XLES_RESOLVED_URc   (0,0,0))
  ALLOCATE(XLES_RESOLVED_VRc   (0,0,0))
  ALLOCATE(XLES_RESOLVED_WRc   (0,0,0))
  ALLOCATE(XLES_RESOLVED_WRc2  (0,0,0))
  ALLOCATE(XLES_RESOLVED_W2Rc  (0,0,0))
  ALLOCATE(XLES_RESOLVED_RcPz  (0,0,0))
  ALLOCATE(XLES_RESOLVED_WThlRc(0,0,0))
  ALLOCATE(XLES_RESOLVED_WRvRc (0,0,0))
  ALLOCATE(XLES_RESOLVED_WRt   (0,0,0))
  ALLOCATE(XLES_RESOLVED_Rt2  (0,0,0))
  ALLOCATE(XLES_RESOLVED_RtPz  (0,0,0))
END IF
IF (LUSERI ) THEN
  ALLOCATE(XLES_RESOLVED_Ri2   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Ri'2>
  ALLOCATE(XLES_RESOLVED_ThRi  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Th'Ri'>
  ALLOCATE(XLES_RESOLVED_ThlRi (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thl'Ri'>
  ALLOCATE(XLES_RESOLVED_ThvRi (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thv'Ri'>
  ALLOCATE(XLES_RESOLVED_URi   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'Ri'>
  ALLOCATE(XLES_RESOLVED_VRi   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'Ri'>
  ALLOCATE(XLES_RESOLVED_WRi   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Ri'>
  ALLOCATE(XLES_RESOLVED_WRi2  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Ri'2>
  ALLOCATE(XLES_RESOLVED_W2Ri  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'2Ri'>
  ALLOCATE(XLES_RESOLVED_RiPz  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Ri'dp'/dz>
  ALLOCATE(XLES_RESOLVED_WThlRi(NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Thl'Ri'>
  ALLOCATE(XLES_RESOLVED_WRvRi (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rv'Ri'>
ELSE
  ALLOCATE(XLES_RESOLVED_Ri2   (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThRi  (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThlRi (0,0,0))
  ALLOCATE(XLES_RESOLVED_ThvRi (0,0,0))
  ALLOCATE(XLES_RESOLVED_URi   (0,0,0))
  ALLOCATE(XLES_RESOLVED_VRi   (0,0,0))
  ALLOCATE(XLES_RESOLVED_WRi   (0,0,0))
  ALLOCATE(XLES_RESOLVED_WRi2  (0,0,0))
  ALLOCATE(XLES_RESOLVED_W2Ri  (0,0,0))
  ALLOCATE(XLES_RESOLVED_RiPz  (0,0,0))
  ALLOCATE(XLES_RESOLVED_WThlRi(0,0,0))
  ALLOCATE(XLES_RESOLVED_WRvRi (0,0,0))
END IF
!
IF (LUSERR) THEN
  ALLOCATE(XLES_RESOLVED_WRr   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rr'>
  ALLOCATE(XLES_INPRR3D        (NLES_K,NLES_TIMES,NLES_MASKS)) !precip flux
  ALLOCATE(XLES_MAX_INPRR3D        (NLES_K,NLES_TIMES,NLES_MASKS)) !precip flux
  ALLOCATE(XLES_EVAP3D        (NLES_K,NLES_TIMES,NLES_MASKS)) ! evap
ELSE
  ALLOCATE(XLES_RESOLVED_WRr   (0,0,0))
  ALLOCATE(XLES_INPRR3D        (0,0,0))
  ALLOCATE(XLES_MAX_INPRR3D        (0,0,0))
  ALLOCATE(XLES_EVAP3D         (0,0,0))
END IF
IF (NSV>0  ) THEN
  ALLOCATE(XLES_RESOLVED_Sv2   (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <Sv'2>
  ALLOCATE(XLES_RESOLVED_ThSv  (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <Th'Sv>
  ALLOCATE(XLES_RESOLVED_USv   (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <u'Sv'>
  ALLOCATE(XLES_RESOLVED_VSv   (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <v'Sv'>
  ALLOCATE(XLES_RESOLVED_WSv   (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <w'Sv'>
  ALLOCATE(XLES_RESOLVED_WSv2  (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <w'Sv'2>
  ALLOCATE(XLES_RESOLVED_W2Sv  (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <w'2Sv'>
  ALLOCATE(XLES_RESOLVED_SvPz  (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <Sv'dp'/dz>
  ALLOCATE(XLES_RESOLVED_WThlSv(NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <w'Thl'Sv'>
  IF (LUSERV) THEN
    ALLOCATE(XLES_RESOLVED_ThvSv (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <Thv'Sv>
    ALLOCATE(XLES_RESOLVED_WRvSv (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <w'Rv'Sv'>
  ELSE
    ALLOCATE(XLES_RESOLVED_ThvSv (0,0,0,0))
    ALLOCATE(XLES_RESOLVED_WRvSv (0,0,0,0))
  END IF
  IF (LUSERC) THEN
    ALLOCATE(XLES_RESOLVED_ThlSv (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <Thl'Sv>
  ELSE
    ALLOCATE(XLES_RESOLVED_ThlSv (0,0,0,0))
  END IF
ELSE
  ALLOCATE(XLES_RESOLVED_Sv2   (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_ThSv  (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_USv   (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_VSv   (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_WSv   (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_WSv2  (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_W2Sv  (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_SvPz  (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_ThvSv (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_ThlSv (0,0,0,0))
  ALLOCATE(XLES_RESOLVED_WThlSv(0,0,0,0))
  ALLOCATE(XLES_RESOLVED_WRvSv (0,0,0,0))
END IF
!
!
XLES_RESOLVED_U2  = XUNDEF
XLES_RESOLVED_V2  = XUNDEF
XLES_RESOLVED_W2  = XUNDEF
XLES_RESOLVED_P2  = XUNDEF
XLES_RESOLVED_Th2 = XUNDEF
IF( LUSERC) THEN
  XLES_RESOLVED_Thl2= XUNDEF
  XLES_RESOLVED_ThlThv= XUNDEF
END IF
IF (LUSERV) THEN
  XLES_RESOLVED_ThThv = XUNDEF
END IF
XLES_RESOLVED_Ke  = XUNDEF
XLES_RESOLVED_UV  = XUNDEF
XLES_RESOLVED_WU  = XUNDEF
XLES_RESOLVED_WV  = XUNDEF
XLES_RESOLVED_UP  = XUNDEF
XLES_RESOLVED_VP  = XUNDEF
XLES_RESOLVED_WP  = XUNDEF
XLES_RESOLVED_UTh = XUNDEF
XLES_RESOLVED_VTh = XUNDEF
XLES_RESOLVED_WTh = XUNDEF
IF (LUSERC) THEN
  XLES_RESOLVED_UThl= XUNDEF
  XLES_RESOLVED_VThl= XUNDEF
  XLES_RESOLVED_WThl= XUNDEF
END IF
IF (LUSERV) THEN
  XLES_RESOLVED_UThv= XUNDEF
  XLES_RESOLVED_VThv= XUNDEF
  XLES_RESOLVED_WThv= XUNDEF
END IF
XLES_RESOLVED_U3  = XUNDEF
XLES_RESOLVED_V3  = XUNDEF
XLES_RESOLVED_W3  = XUNDEF
XLES_RESOLVED_U4  = XUNDEF
XLES_RESOLVED_V4  = XUNDEF
XLES_RESOLVED_W4  = XUNDEF
XLES_RESOLVED_WThl2  = XUNDEF
XLES_RESOLVED_W2Thl  = XUNDEF
XLES_RESOLVED_ThlPz  = XUNDEF
!
XLES_RESOLVED_MASSFX = XUNDEF
XLES_RESOLVED_UKe = XUNDEF
XLES_RESOLVED_VKe = XUNDEF
XLES_RESOLVED_WKe = XUNDEF
IF (LUSERV ) THEN
  XLES_RESOLVED_Rv2  = XUNDEF
  XLES_RESOLVED_ThRv = XUNDEF
  IF (LUSERC) XLES_RESOLVED_ThlRv= XUNDEF
  XLES_RESOLVED_ThvRv= XUNDEF
  XLES_RESOLVED_URv = XUNDEF
  XLES_RESOLVED_VRv = XUNDEF
  XLES_RESOLVED_WRv = XUNDEF
  XLES_RESOLVED_WRv2  = XUNDEF
  XLES_RESOLVED_W2Rv  = XUNDEF
  XLES_RESOLVED_WRt2 = XUNDEF
  XLES_RESOLVED_W2Rt  = XUNDEF
  XLES_RESOLVED_WThlRv= XUNDEF
  XLES_RESOLVED_WThlRt= XUNDEF  
  XLES_RESOLVED_RvPz  = XUNDEF
END IF
IF (LUSERC ) THEN
  XLES_RESOLVED_Rc2  = XUNDEF
  XLES_RESOLVED_ThRc = XUNDEF
  XLES_RESOLVED_ThlRc= XUNDEF
  XLES_RESOLVED_ThvRc= XUNDEF
  XLES_RESOLVED_URc  = XUNDEF
  XLES_RESOLVED_VRc  = XUNDEF
  XLES_RESOLVED_WRc  = XUNDEF
  XLES_RESOLVED_WRc2  = XUNDEF
  XLES_RESOLVED_W2Rc  = XUNDEF
  XLES_RESOLVED_WThlRc= XUNDEF
  XLES_RESOLVED_WRvRc = XUNDEF
  XLES_RESOLVED_RcPz  = XUNDEF
  XLES_RESOLVED_RtPz  = XUNDEF
  XLES_RESOLVED_WRt  = XUNDEF
  XLES_RESOLVED_Rt2  = XUNDEF
END IF
IF (LUSERI ) THEN
  XLES_RESOLVED_Ri2  = XUNDEF
  XLES_RESOLVED_ThRi = XUNDEF
  XLES_RESOLVED_ThlRi= XUNDEF
  XLES_RESOLVED_ThvRi= XUNDEF
  XLES_RESOLVED_URi  = XUNDEF
  XLES_RESOLVED_VRi  = XUNDEF
  XLES_RESOLVED_WRi  = XUNDEF
  XLES_RESOLVED_WRi2  = XUNDEF
  XLES_RESOLVED_W2Ri  = XUNDEF
  XLES_RESOLVED_WThlRi= XUNDEF
  XLES_RESOLVED_WRvRi = XUNDEF
  XLES_RESOLVED_RiPz  = XUNDEF
END IF
!
IF (LUSERR) XLES_RESOLVED_WRr = XUNDEF
IF (LUSERR) XLES_MAX_INPRR3D = XUNDEF
IF (LUSERR) XLES_INPRR3D = XUNDEF
IF (LUSERR) XLES_EVAP3D = XUNDEF
IF (NSV>0  ) THEN
  XLES_RESOLVED_Sv2  = XUNDEF
  XLES_RESOLVED_ThSv = XUNDEF
  IF (LUSERC) XLES_RESOLVED_ThlSv= XUNDEF
  IF (LUSERV) XLES_RESOLVED_ThvSv= XUNDEF
  XLES_RESOLVED_USv = XUNDEF
  XLES_RESOLVED_VSv = XUNDEF
  XLES_RESOLVED_WSv = XUNDEF
  XLES_RESOLVED_WSv2  = XUNDEF
  XLES_RESOLVED_W2Sv  = XUNDEF
  XLES_RESOLVED_WThlSv= XUNDEF
  IF (LUSERV) XLES_RESOLVED_WRvSv = XUNDEF
  XLES_RESOLVED_SvPz  = XUNDEF
END IF
!
!
!*      7.4  interactions of resolved and subgrid quantities
!            -----------------------------------------------
!
ALLOCATE(XLES_RES_U_SBG_Tke         (NLES_K,NLES_TIMES,NLES_MASKS))! <u'Tke>
ALLOCATE(XLES_RES_V_SBG_Tke         (NLES_K,NLES_TIMES,NLES_MASKS))! <v'Tke>
ALLOCATE(XLES_RES_W_SBG_Tke         (NLES_K,NLES_TIMES,NLES_MASKS))! <w'Tke>
!                                                                        ______
ALLOCATE(XLES_RES_W_SBG_WThl        (NLES_K,NLES_TIMES,NLES_MASKS))!  <w'w'Thl'>
!                                                                        _____
ALLOCATE(XLES_RES_W_SBG_Thl2        (NLES_K,NLES_TIMES,NLES_MASKS))!  <w'Thl'2>
!                                                                              _____
ALLOCATE(XLES_RES_ddxa_U_SBG_UaU    (NLES_K,NLES_TIMES,NLES_MASKS))!  <du'/dxa ua'u'>
!                                                                              _____
ALLOCATE(XLES_RES_ddxa_V_SBG_UaV    (NLES_K,NLES_TIMES,NLES_MASKS))!  <dv'/dxa ua'v'>
!                                                                              _____
ALLOCATE(XLES_RES_ddxa_W_SBG_UaW    (NLES_K,NLES_TIMES,NLES_MASKS))!  <dw'/dxa ua'w'>
!                                                                              _______
ALLOCATE(XLES_RES_ddxa_W_SBG_UaThl  (NLES_K,NLES_TIMES,NLES_MASKS))!  <dw'/dxa ua'Thl'>
!                                                                                _____
ALLOCATE(XLES_RES_ddxa_Thl_SBG_UaW  (NLES_K,NLES_TIMES,NLES_MASKS))!  <dThl'/dxa ua'w'>
!                                                                               ___
ALLOCATE(XLES_RES_ddz_Thl_SBG_W2    (NLES_K,NLES_TIMES,NLES_MASKS))!  <dThl'/dz w'2>
!                                                                                _______
ALLOCATE(XLES_RES_ddxa_Thl_SBG_UaThl(NLES_K,NLES_TIMES,NLES_MASKS))!  <dThl'/dxa ua'Thl'>
!
IF (LUSERV) THEN
!                                                                          _____
  ALLOCATE(XLES_RES_W_SBG_WRt         (NLES_K,NLES_TIMES,NLES_MASKS))!  <w'w'Rt'>
!                                                                           ____
  ALLOCATE(XLES_RES_W_SBG_Rt2         (NLES_K,NLES_TIMES,NLES_MASKS))!  <w'Rt'2>
!                                                                          _______
  ALLOCATE(XLES_RES_W_SBG_ThlRt       (NLES_K,NLES_TIMES,NLES_MASKS))!  <w'Thl'Rt'>
!                                                                                ______
  ALLOCATE(XLES_RES_ddxa_W_SBG_UaRt   (NLES_K,NLES_TIMES,NLES_MASKS))!  <dw'/dxa ua'Rt'>
!                                                                                 _____
  ALLOCATE(XLES_RES_ddxa_Rt_SBG_UaW   (NLES_K,NLES_TIMES,NLES_MASKS))!  <dRt'/dxa ua'w'>
!                                                                                ___
  ALLOCATE(XLES_RES_ddz_Rt_SBG_W2     (NLES_K,NLES_TIMES,NLES_MASKS))!  <dRt'/dz w'2>
!                                                                                  ______
  ALLOCATE(XLES_RES_ddxa_Thl_SBG_UaRt (NLES_K,NLES_TIMES,NLES_MASKS))!  <dThl'/dxa ua'Rt'>
!                                                                                 _______
  ALLOCATE(XLES_RES_ddxa_Rt_SBG_UaThl (NLES_K,NLES_TIMES,NLES_MASKS))!  <dRt'/dxa ua'Thl'>
!                                                                                  ______
  ALLOCATE(XLES_RES_ddxa_Rt_SBG_UaRt  (NLES_K,NLES_TIMES,NLES_MASKS)) !  <dRt'/dxa ua'Rt'>
ELSE
  ALLOCATE(XLES_RES_W_SBG_WRt         (0,0,0))
  ALLOCATE(XLES_RES_W_SBG_Rt2         (0,0,0))
  ALLOCATE(XLES_RES_W_SBG_ThlRt       (0,0,0))
  ALLOCATE(XLES_RES_ddxa_W_SBG_UaRt   (0,0,0))
  ALLOCATE(XLES_RES_ddxa_Rt_SBG_UaW   (0,0,0))
  ALLOCATE(XLES_RES_ddz_Rt_SBG_W2     (0,0,0))
  ALLOCATE(XLES_RES_ddxa_Thl_SBG_UaRt (0,0,0))
  ALLOCATE(XLES_RES_ddxa_Rt_SBG_UaThl (0,0,0))
  ALLOCATE(XLES_RES_ddxa_Rt_SBG_UaRt  (0,0,0))
END IF
!
!                                                                                  ______
ALLOCATE(XLES_RES_ddxa_W_SBG_UaSv (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  !  <dw'/dxa ua'Sv'>
!                                                                                   _____
ALLOCATE(XLES_RES_ddxa_Sv_SBG_UaW (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  !  <dSv'/dxa ua'w'>
!                                                                                   ___
ALLOCATE(XLES_RES_ddz_Sv_SBG_W2   (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  !  <dSv'/dz w'2>
!                                                                                  ______
ALLOCATE(XLES_RES_ddxa_Sv_SBG_UaSv(NLES_K,NLES_TIMES,NLES_MASKS,NSV))  !  <dSv'/dxa ua'Sv'>
!                                                                            _____
ALLOCATE(XLES_RES_W_SBG_WSv       (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  !  <w'w'Sv'>
!                                                                            ____
ALLOCATE(XLES_RES_W_SBG_Sv2       (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  !  <w'Sv'2>
!
XLES_RES_U_SBG_Tke= XUNDEF
XLES_RES_V_SBG_Tke= XUNDEF
XLES_RES_W_SBG_Tke= XUNDEF
XLES_RES_W_SBG_WThl         = XUNDEF
XLES_RES_W_SBG_Thl2         = XUNDEF
XLES_RES_ddxa_U_SBG_UaU     = XUNDEF
XLES_RES_ddxa_V_SBG_UaV     = XUNDEF
XLES_RES_ddxa_W_SBG_UaW     = XUNDEF
XLES_RES_ddxa_W_SBG_UaThl   = XUNDEF
XLES_RES_ddxa_Thl_SBG_UaW   = XUNDEF
XLES_RES_ddz_Thl_SBG_W2     = XUNDEF
XLES_RES_ddxa_Thl_SBG_UaThl = XUNDEF
IF (LUSERV) THEN
  XLES_RES_W_SBG_WRt        = XUNDEF
  XLES_RES_W_SBG_Rt2        = XUNDEF
  XLES_RES_W_SBG_ThlRt      = XUNDEF
  XLES_RES_ddxa_W_SBG_UaRt  = XUNDEF
  XLES_RES_ddxa_Rt_SBG_UaW  = XUNDEF
  XLES_RES_ddz_Rt_SBG_W2    = XUNDEF
  XLES_RES_ddxa_Thl_SBG_UaRt= XUNDEF
  XLES_RES_ddxa_Rt_SBG_UaThl= XUNDEF
  XLES_RES_ddxa_Rt_SBG_UaRt = XUNDEF
END IF
IF (NSV>0) THEN
  XLES_RES_ddxa_W_SBG_UaSv = XUNDEF
  XLES_RES_ddxa_Sv_SBG_UaW = XUNDEF
  XLES_RES_ddz_Sv_SBG_W2   = XUNDEF
  XLES_RES_ddxa_Sv_SBG_UaSv= XUNDEF
  XLES_RES_W_SBG_WSv       = XUNDEF
  XLES_RES_W_SBG_Sv2       = XUNDEF
END IF
!
!
!*      7.5  subgrid quantities
!            ------------------
!
ALLOCATE(XLES_SUBGRID_U2    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'2>
ALLOCATE(XLES_SUBGRID_V2    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'2>
ALLOCATE(XLES_SUBGRID_W2    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'2>
ALLOCATE(XLES_SUBGRID_Tke   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <e>
ALLOCATE(XLES_SUBGRID_Thl2  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thl'2>
ALLOCATE(XLES_SUBGRID_UV    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'v'>
ALLOCATE(XLES_SUBGRID_WU    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'u'>
ALLOCATE(XLES_SUBGRID_WV    (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'v'>
ALLOCATE(XLES_SUBGRID_UThl  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'Thl'>
ALLOCATE(XLES_SUBGRID_VThl  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'Thl'>
ALLOCATE(XLES_SUBGRID_WThl  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Thl'>
ALLOCATE(XLES_SUBGRID_WThv     (NLES_K,NLES_TIMES,NLES_MASKS))  ! <w'Thv'>
ALLOCATE(XLES_SUBGRID_ThlThv   (NLES_K,NLES_TIMES,NLES_MASKS))  ! <Thl'Thv'>
ALLOCATE(XLES_SUBGRID_W2Thl    (NLES_K,NLES_TIMES,NLES_MASKS))  ! <w'2Thl>
ALLOCATE(XLES_SUBGRID_WThl2    (NLES_K,NLES_TIMES,NLES_MASKS))  ! <w'Thl'2>
ALLOCATE(XLES_SUBGRID_DISS_Tke (NLES_K,NLES_TIMES,NLES_MASKS))  ! <epsilon>
ALLOCATE(XLES_SUBGRID_DISS_Thl2(NLES_K,NLES_TIMES,NLES_MASKS))  ! <epsilon_Thl2>
ALLOCATE(XLES_SUBGRID_WP       (NLES_K,NLES_TIMES,NLES_MASKS))  ! <w'p'>
ALLOCATE(XLES_SUBGRID_PHI3     (NLES_K,NLES_TIMES,NLES_MASKS))  ! phi3
ALLOCATE(XLES_SUBGRID_LMix     (NLES_K,NLES_TIMES,NLES_MASKS))  ! mixing length
ALLOCATE(XLES_SUBGRID_LDiss    (NLES_K,NLES_TIMES,NLES_MASKS))  ! dissipative length
ALLOCATE(XLES_SUBGRID_Km       (NLES_K,NLES_TIMES,NLES_MASKS))  ! Km
ALLOCATE(XLES_SUBGRID_Kh       (NLES_K,NLES_TIMES,NLES_MASKS))  ! Kh
ALLOCATE(XLES_SUBGRID_ThlPz    (NLES_K,NLES_TIMES,NLES_MASKS))  ! <Thl'dp'/dz>
ALLOCATE(XLES_SUBGRID_UTke  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'Tke>
ALLOCATE(XLES_SUBGRID_VTke  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'Tke>
ALLOCATE(XLES_SUBGRID_WTke  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Tke>
ALLOCATE(XLES_SUBGRID_ddz_WTke  (NLES_K,NLES_TIMES,NLES_MASKS)) ! <dw'Tke/dz>

ALLOCATE(XLES_SUBGRID_THLUP_MF(NLES_K,NLES_TIMES,NLES_MASKS))  ! Thl of the Updraft
ALLOCATE(XLES_SUBGRID_RTUP_MF (NLES_K,NLES_TIMES,NLES_MASKS))  ! Rt of the Updraft
ALLOCATE(XLES_SUBGRID_RVUP_MF (NLES_K,NLES_TIMES,NLES_MASKS))  ! Rv of the Updraft
ALLOCATE(XLES_SUBGRID_RCUP_MF (NLES_K,NLES_TIMES,NLES_MASKS))  ! Rc of the Updraft
ALLOCATE(XLES_SUBGRID_RIUP_MF (NLES_K,NLES_TIMES,NLES_MASKS))  ! Ri of the Updraft
ALLOCATE(XLES_SUBGRID_WUP_MF  (NLES_K,NLES_TIMES,NLES_MASKS))  ! Thl of the Updraft
ALLOCATE(XLES_SUBGRID_MASSFLUX(NLES_K,NLES_TIMES,NLES_MASKS))  ! Mass Flux
ALLOCATE(XLES_SUBGRID_DETR    (NLES_K,NLES_TIMES,NLES_MASKS))  ! Detrainment
ALLOCATE(XLES_SUBGRID_ENTR    (NLES_K,NLES_TIMES,NLES_MASKS))  ! Entrainment
ALLOCATE(XLES_SUBGRID_FRACUP  (NLES_K,NLES_TIMES,NLES_MASKS))  ! Updraft Fraction 
ALLOCATE(XLES_SUBGRID_THVUP_MF(NLES_K,NLES_TIMES,NLES_MASKS))  ! Thv of the Updraft
ALLOCATE(XLES_SUBGRID_WTHLMF  (NLES_K,NLES_TIMES,NLES_MASKS))  ! Flux of thl   
ALLOCATE(XLES_SUBGRID_WRTMF   (NLES_K,NLES_TIMES,NLES_MASKS)) ! Flux of rt
ALLOCATE(XLES_SUBGRID_WTHVMF  (NLES_K,NLES_TIMES,NLES_MASKS)) ! Flux of thv 
ALLOCATE(XLES_SUBGRID_WUMF    (NLES_K,NLES_TIMES,NLES_MASKS))! Flux of u
ALLOCATE(XLES_SUBGRID_WVMF    (NLES_K,NLES_TIMES,NLES_MASKS))! Flux of v

IF (LUSERV ) THEN
  ALLOCATE(XLES_SUBGRID_Rt2   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rt'2>
  ALLOCATE(XLES_SUBGRID_ThlRt (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Thl'Rt'>
  ALLOCATE(XLES_SUBGRID_URt   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'Rt'>
  ALLOCATE(XLES_SUBGRID_VRt   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'Rt'>
  ALLOCATE(XLES_SUBGRID_WRt   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rt'>
  ALLOCATE(XLES_SUBGRID_RtThv (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rt'Thv'>
  ALLOCATE(XLES_SUBGRID_W2Rt  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'2Rt'>
  ALLOCATE(XLES_SUBGRID_WThlRt(NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Thl'Rt'>
  ALLOCATE(XLES_SUBGRID_WRt2  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rt'2>
  ALLOCATE(XLES_SUBGRID_DISS_Rt2 (NLES_K,NLES_TIMES,NLES_MASKS))  ! <epsilon_Rt2>
  ALLOCATE(XLES_SUBGRID_DISS_ThlRt(NLES_K,NLES_TIMES,NLES_MASKS)) ! <epsilon_ThlRt>
  ALLOCATE(XLES_SUBGRID_RtPz  (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rt'dp'/dz>
  ALLOCATE(XLES_SUBGRID_PSI3  (NLES_K,NLES_TIMES,NLES_MASKS))     ! psi3
ELSE
  ALLOCATE(XLES_SUBGRID_Rt2   (0,0,0))
  ALLOCATE(XLES_SUBGRID_ThlRt (0,0,0))
  ALLOCATE(XLES_SUBGRID_URt   (0,0,0))
  ALLOCATE(XLES_SUBGRID_VRt   (0,0,0))
  ALLOCATE(XLES_SUBGRID_WRt   (0,0,0))
  ALLOCATE(XLES_SUBGRID_RtThv (0,0,0))
  ALLOCATE(XLES_SUBGRID_W2Rt  (0,0,0))
  ALLOCATE(XLES_SUBGRID_WThlRt(0,0,0))
  ALLOCATE(XLES_SUBGRID_WRt2  (0,0,0))
  ALLOCATE(XLES_SUBGRID_DISS_Rt2 (0,0,0))
  ALLOCATE(XLES_SUBGRID_DISS_ThlRt(0,0,0))
  ALLOCATE(XLES_SUBGRID_RtPz  (0,0,0))
  ALLOCATE(XLES_SUBGRID_PSI3  (0,0,0))
END IF
IF (LUSERC ) THEN
  ALLOCATE(XLES_SUBGRID_Rc2   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Rc'2>
  ALLOCATE(XLES_SUBGRID_URc   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <u'Rc'>
  ALLOCATE(XLES_SUBGRID_VRc   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <v'Rc'>
  ALLOCATE(XLES_SUBGRID_WRc   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <w'Rc'>
ELSE
  ALLOCATE(XLES_SUBGRID_Rc2   (0,0,0))
  ALLOCATE(XLES_SUBGRID_URc   (0,0,0))
  ALLOCATE(XLES_SUBGRID_VRc   (0,0,0))
  ALLOCATE(XLES_SUBGRID_WRc   (0,0,0))
END IF
IF (LUSERI ) THEN
  ALLOCATE(XLES_SUBGRID_Ri2   (NLES_K,NLES_TIMES,NLES_MASKS))     ! <Ri'2>
ELSE
  ALLOCATE(XLES_SUBGRID_Ri2   (0,0,0))
END IF
IF (NSV>0  ) THEN
  ALLOCATE(XLES_SUBGRID_USv   (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <u'Sv'>
  ALLOCATE(XLES_SUBGRID_VSv   (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <v'Sv'>
  ALLOCATE(XLES_SUBGRID_WSv   (NLES_K,NLES_TIMES,NLES_MASKS,NSV)) ! <w'Sv'>
  ALLOCATE(XLES_SUBGRID_Sv2      (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  ! <Sv'2>
  ALLOCATE(XLES_SUBGRID_SvThv    (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  ! <Sv'Thv'>
  ALLOCATE(XLES_SUBGRID_W2Sv     (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  ! <w'2Sv'>
  ALLOCATE(XLES_SUBGRID_WSv2     (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  ! <w'Sv'2>
  ALLOCATE(XLES_SUBGRID_DISS_Sv2 (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  ! <epsilon_Sv2>
  ALLOCATE(XLES_SUBGRID_SvPz     (NLES_K,NLES_TIMES,NLES_MASKS,NSV))  ! <Sv'dp'/dz>
ELSE
  ALLOCATE(XLES_SUBGRID_USv   (0,0,0,0))
  ALLOCATE(XLES_SUBGRID_VSv   (0,0,0,0))
  ALLOCATE(XLES_SUBGRID_WSv   (0,0,0,0))
  ALLOCATE(XLES_SUBGRID_Sv2     (0,0,0,0))
  ALLOCATE(XLES_SUBGRID_SvThv   (0,0,0,0))
  ALLOCATE(XLES_SUBGRID_W2Sv    (0,0,0,0))
  ALLOCATE(XLES_SUBGRID_WSv2    (0,0,0,0))
  ALLOCATE(XLES_SUBGRID_DISS_Sv2(0,0,0,0))
  ALLOCATE(XLES_SUBGRID_SvPz    (0,0,0,0))
END IF
!
XLES_SUBGRID_U2  = XUNDEF
XLES_SUBGRID_V2  = XUNDEF
XLES_SUBGRID_W2  = XUNDEF
XLES_SUBGRID_Tke = XUNDEF
XLES_SUBGRID_Thl2= XUNDEF
XLES_SUBGRID_UV  = XUNDEF
XLES_SUBGRID_WU  = XUNDEF
XLES_SUBGRID_WV  = XUNDEF
XLES_SUBGRID_UThl= XUNDEF
XLES_SUBGRID_VThl= XUNDEF
XLES_SUBGRID_WThl= XUNDEF
XLES_SUBGRID_WThv= XUNDEF
XLES_SUBGRID_ThlThv= XUNDEF
XLES_SUBGRID_W2Thl= XUNDEF
XLES_SUBGRID_WThl2    = XUNDEF
XLES_SUBGRID_DISS_Tke = XUNDEF
XLES_SUBGRID_DISS_Thl2= XUNDEF
XLES_SUBGRID_WP       = XUNDEF
XLES_SUBGRID_PHI3     = XUNDEF
XLES_SUBGRID_LMix     = XUNDEF
XLES_SUBGRID_LDiss    = XUNDEF
XLES_SUBGRID_Km       = XUNDEF
XLES_SUBGRID_Kh       = XUNDEF
XLES_SUBGRID_ThlPz    = XUNDEF
XLES_SUBGRID_UTke= XUNDEF
XLES_SUBGRID_VTke= XUNDEF
XLES_SUBGRID_WTke= XUNDEF
XLES_SUBGRID_ddz_WTke  = XUNDEF

XLES_SUBGRID_THLUP_MF = XUNDEF
XLES_SUBGRID_RTUP_MF  = XUNDEF
XLES_SUBGRID_RVUP_MF  = XUNDEF
XLES_SUBGRID_RCUP_MF  = XUNDEF
XLES_SUBGRID_RIUP_MF  = XUNDEF
XLES_SUBGRID_WUP_MF   = XUNDEF
XLES_SUBGRID_MASSFLUX = XUNDEF
XLES_SUBGRID_DETR     = XUNDEF
XLES_SUBGRID_ENTR     = XUNDEF
XLES_SUBGRID_FRACUP   = XUNDEF
XLES_SUBGRID_THVUP_MF = XUNDEF
XLES_SUBGRID_WTHLMF   = XUNDEF  
XLES_SUBGRID_WRTMF    = XUNDEF
XLES_SUBGRID_WTHVMF   = XUNDEF
XLES_SUBGRID_WUMF     = XUNDEF
XLES_SUBGRID_WVMF     = XUNDEF

IF (LUSERV ) THEN
  XLES_SUBGRID_Rt2 = XUNDEF
  XLES_SUBGRID_ThlRt= XUNDEF
  XLES_SUBGRID_URt = XUNDEF
  XLES_SUBGRID_VRt = XUNDEF
  XLES_SUBGRID_WRt = XUNDEF
  XLES_SUBGRID_RtThv   = XUNDEF
  XLES_SUBGRID_W2Rt    = XUNDEF
  XLES_SUBGRID_WThlRt  = XUNDEF
  XLES_SUBGRID_WRt2    = XUNDEF
  XLES_SUBGRID_DISS_Rt2= XUNDEF
  XLES_SUBGRID_DISS_ThlRt= XUNDEF
  XLES_SUBGRID_RtPz    = XUNDEF
  XLES_SUBGRID_PSI3     = XUNDEF
END IF
IF (LUSERC ) THEN
  XLES_SUBGRID_Rc2 = XUNDEF
  XLES_SUBGRID_URc = XUNDEF
  XLES_SUBGRID_VRc = XUNDEF
  XLES_SUBGRID_WRc = XUNDEF
END IF
IF (LUSERI ) THEN
  XLES_SUBGRID_Ri2 = XUNDEF
END IF
IF (NSV>0  ) THEN
  XLES_SUBGRID_USv = XUNDEF
  XLES_SUBGRID_VSv = XUNDEF
  XLES_SUBGRID_WSv = XUNDEF
  XLES_SUBGRID_Sv2     = XUNDEF
  XLES_SUBGRID_SvThv   = XUNDEF
  XLES_SUBGRID_W2Sv    = XUNDEF
  XLES_SUBGRID_WSv2    = XUNDEF
  XLES_SUBGRID_DISS_Sv2= XUNDEF
  XLES_SUBGRID_SvPz    = XUNDEF
END IF
!
!
!*      7.6  updraft quantities (only on the cartesian mask)
!            ------------------
!
ALLOCATE(XLES_UPDRAFT       (NLES_K,NLES_TIMES))    ! updraft fraction
ALLOCATE(XLES_UPDRAFT_W     (NLES_K,NLES_TIMES))    ! <w>
ALLOCATE(XLES_UPDRAFT_Th    (NLES_K,NLES_TIMES))    ! <theta>
ALLOCATE(XLES_UPDRAFT_Ke    (NLES_K,NLES_TIMES))    ! <E>
ALLOCATE(XLES_UPDRAFT_WTh   (NLES_K,NLES_TIMES))    ! <w'theta'>
ALLOCATE(XLES_UPDRAFT_Th2   (NLES_K,NLES_TIMES))    ! <th'2>
ALLOCATE(XLES_UPDRAFT_Tke   (NLES_K,NLES_TIMES))    ! <e>

IF (LUSERV) THEN
  ALLOCATE(XLES_UPDRAFT_Thv   (NLES_K,NLES_TIMES))    ! <thetav>
  ALLOCATE(XLES_UPDRAFT_WThv  (NLES_K,NLES_TIMES))    ! <w'thv'>
  ALLOCATE(XLES_UPDRAFT_ThThv (NLES_K,NLES_TIMES))    ! <th'thv'>
ELSE
  ALLOCATE(XLES_UPDRAFT_Thv   (0,0))
  ALLOCATE(XLES_UPDRAFT_WThv  (0,0))
  ALLOCATE(XLES_UPDRAFT_ThThv (0,0))
END IF
!
IF (LUSERC) THEN
  ALLOCATE(XLES_UPDRAFT_Thl   (NLES_K,NLES_TIMES))    ! <thetal>
  ALLOCATE(XLES_UPDRAFT_WThl  (NLES_K,NLES_TIMES))    ! <w'thetal'>
  ALLOCATE(XLES_UPDRAFT_Thl2  (NLES_K,NLES_TIMES))    ! <thl'2>
  ALLOCATE(XLES_UPDRAFT_ThlThv(NLES_K,NLES_TIMES))    ! <thl'thv'>
ELSE
  ALLOCATE(XLES_UPDRAFT_Thl   (0,0))
  ALLOCATE(XLES_UPDRAFT_WThl  (0,0))
  ALLOCATE(XLES_UPDRAFT_Thl2  (0,0))
  ALLOCATE(XLES_UPDRAFT_ThlThv(0,0))
END IF

IF (LUSERV ) THEN
  ALLOCATE(XLES_UPDRAFT_Rv    (NLES_K,NLES_TIMES))    ! <Rv>
  ALLOCATE(XLES_UPDRAFT_WRv   (NLES_K,NLES_TIMES))    ! <w'Rv'>
  ALLOCATE(XLES_UPDRAFT_Rv2   (NLES_K,NLES_TIMES))    ! <Rv'2>
  ALLOCATE(XLES_UPDRAFT_ThRv  (NLES_K,NLES_TIMES))    ! <Th'Rv'>
  ALLOCATE(XLES_UPDRAFT_ThvRv (NLES_K,NLES_TIMES))    ! <Thv'Rv'>
  IF (LUSERC) THEN
    ALLOCATE(XLES_UPDRAFT_ThlRv (NLES_K,NLES_TIMES))    ! <Thl'Rv'>
  ELSE
    ALLOCATE(XLES_UPDRAFT_ThlRv (0,0))
  END IF
ELSE
  ALLOCATE(XLES_UPDRAFT_Rv    (0,0))
  ALLOCATE(XLES_UPDRAFT_WRv   (0,0))
  ALLOCATE(XLES_UPDRAFT_Rv2   (0,0))
  ALLOCATE(XLES_UPDRAFT_ThRv  (0,0))
  ALLOCATE(XLES_UPDRAFT_ThvRv (0,0))
  ALLOCATE(XLES_UPDRAFT_ThlRv (0,0))
END IF
IF (LUSERC ) THEN
  ALLOCATE(XLES_UPDRAFT_Rc    (NLES_K,NLES_TIMES))    ! <Rc>
  ALLOCATE(XLES_UPDRAFT_WRc   (NLES_K,NLES_TIMES))    ! <w'Rc'>
  ALLOCATE(XLES_UPDRAFT_Rc2   (NLES_K,NLES_TIMES))    ! <Rc'2>
  ALLOCATE(XLES_UPDRAFT_ThRc  (NLES_K,NLES_TIMES))    ! <Th'Rc'>
  ALLOCATE(XLES_UPDRAFT_ThvRc (NLES_K,NLES_TIMES))    ! <Thv'Rc'>
  ALLOCATE(XLES_UPDRAFT_ThlRc (NLES_K,NLES_TIMES))    ! <Thl'Rc'>
ELSE
  ALLOCATE(XLES_UPDRAFT_Rc    (0,0))
  ALLOCATE(XLES_UPDRAFT_WRc   (0,0))
  ALLOCATE(XLES_UPDRAFT_Rc2   (0,0))
  ALLOCATE(XLES_UPDRAFT_ThRc  (0,0))
  ALLOCATE(XLES_UPDRAFT_ThvRc (0,0))
  ALLOCATE(XLES_UPDRAFT_ThlRc (0,0))
END IF
IF (LUSERR ) THEN
  ALLOCATE(XLES_UPDRAFT_Rr    (NLES_K,NLES_TIMES))    ! <Rr>
ELSE
  ALLOCATE(XLES_UPDRAFT_Rr    (0,0))
END IF
IF (LUSERI ) THEN
  ALLOCATE(XLES_UPDRAFT_Ri    (NLES_K,NLES_TIMES))    ! <Ri>
  ALLOCATE(XLES_UPDRAFT_WRi   (NLES_K,NLES_TIMES))    ! <w'Ri'>
  ALLOCATE(XLES_UPDRAFT_Ri2   (NLES_K,NLES_TIMES))    ! <Ri'2>
  ALLOCATE(XLES_UPDRAFT_ThRi  (NLES_K,NLES_TIMES))    ! <Th'Ri'>
  ALLOCATE(XLES_UPDRAFT_ThvRi (NLES_K,NLES_TIMES))    ! <Thv'Ri'>
  ALLOCATE(XLES_UPDRAFT_ThlRi (NLES_K,NLES_TIMES))    ! <Thl'Ri'>
ELSE
  ALLOCATE(XLES_UPDRAFT_Ri    (0,0))
  ALLOCATE(XLES_UPDRAFT_WRi   (0,0))
  ALLOCATE(XLES_UPDRAFT_Ri2   (0,0))
  ALLOCATE(XLES_UPDRAFT_ThRi  (0,0))
  ALLOCATE(XLES_UPDRAFT_ThvRi (0,0))
  ALLOCATE(XLES_UPDRAFT_ThlRi (0,0))
END IF
IF (LUSERS ) THEN
  ALLOCATE(XLES_UPDRAFT_Rs    (NLES_K,NLES_TIMES))    ! <Rs>
ELSE
  ALLOCATE(XLES_UPDRAFT_Rs    (0,0))
END IF
IF (LUSERG ) THEN
  ALLOCATE(XLES_UPDRAFT_Rg    (NLES_K,NLES_TIMES))    ! <Rg>
ELSE
  ALLOCATE(XLES_UPDRAFT_Rg    (0,0))
END IF
IF (LUSERH ) THEN
  ALLOCATE(XLES_UPDRAFT_Rh    (NLES_K,NLES_TIMES))    ! <Rh>
ELSE
  ALLOCATE(XLES_UPDRAFT_Rh    (0,0))
END IF
IF (NSV>0  ) THEN
  ALLOCATE(XLES_UPDRAFT_Sv    (NLES_K,NLES_TIMES,NSV))! <Sv>
  ALLOCATE(XLES_UPDRAFT_WSv   (NLES_K,NLES_TIMES,NSV))! <w'Sv'>
  ALLOCATE(XLES_UPDRAFT_Sv2   (NLES_K,NLES_TIMES,NSV))! <Sv'2>
  ALLOCATE(XLES_UPDRAFT_ThSv  (NLES_K,NLES_TIMES,NSV))! <Th'Sv'>
  IF (LUSERV) THEN
    ALLOCATE(XLES_UPDRAFT_ThvSv  (NLES_K,NLES_TIMES,NSV))! <Thv'Sv'>
  ELSE
    ALLOCATE(XLES_UPDRAFT_ThvSv  (0,0,0))
  END IF
  IF (LUSERC) THEN
    ALLOCATE(XLES_UPDRAFT_ThlSv  (NLES_K,NLES_TIMES,NSV))! <Thl'Sv'>
  ELSE
    ALLOCATE(XLES_UPDRAFT_ThlSv  (0,0,0))
  END IF
ELSE
  ALLOCATE(XLES_UPDRAFT_Sv    (0,0,0))
  ALLOCATE(XLES_UPDRAFT_WSv   (0,0,0))
  ALLOCATE(XLES_UPDRAFT_Sv2   (0,0,0))
  ALLOCATE(XLES_UPDRAFT_ThSv  (0,0,0))
  ALLOCATE(XLES_UPDRAFT_ThvSv (0,0,0))
  ALLOCATE(XLES_UPDRAFT_ThlSv (0,0,0))
END IF
!
!
XLES_UPDRAFT     = XUNDEF
XLES_UPDRAFT_W   = XUNDEF
XLES_UPDRAFT_Th  = XUNDEF
XLES_UPDRAFT_Thl = XUNDEF
XLES_UPDRAFT_Tke  = XUNDEF
IF (LUSERV ) XLES_UPDRAFT_Thv = XUNDEF
IF (LUSERC ) XLES_UPDRAFT_Thl = XUNDEF
IF (LUSERV ) XLES_UPDRAFT_Rv  = XUNDEF
IF (LUSERC ) XLES_UPDRAFT_Rc  = XUNDEF
IF (LUSERR ) XLES_UPDRAFT_Rr  = XUNDEF
IF (LUSERI ) XLES_UPDRAFT_Ri  = XUNDEF
IF (LUSERS ) XLES_UPDRAFT_Rs  = XUNDEF
IF (LUSERG ) XLES_UPDRAFT_Rg  = XUNDEF
IF (LUSERH ) XLES_UPDRAFT_Rh  = XUNDEF
IF (NSV>0  ) XLES_UPDRAFT_Sv  = XUNDEF
XLES_UPDRAFT_Ke   = XUNDEF
XLES_UPDRAFT_WTh  = XUNDEF
IF (LUSERV ) XLES_UPDRAFT_WThv = XUNDEF
IF (LUSERC ) XLES_UPDRAFT_WThl = XUNDEF
IF (LUSERV ) XLES_UPDRAFT_WRv  = XUNDEF
IF (LUSERC ) XLES_UPDRAFT_WRc  = XUNDEF
IF (LUSERI ) XLES_UPDRAFT_WRi  = XUNDEF
IF (NSV>0  ) XLES_UPDRAFT_WSv  = XUNDEF
XLES_UPDRAFT_Th2  = XUNDEF
IF (LUSERV ) THEN
  XLES_UPDRAFT_ThThv  = XUNDEF
END IF
IF (LUSERC ) THEN
  XLES_UPDRAFT_Thl2   = XUNDEF
  XLES_UPDRAFT_ThlThv = XUNDEF
END IF
IF (LUSERV ) XLES_UPDRAFT_Rv2  = XUNDEF
IF (LUSERC ) XLES_UPDRAFT_Rc2  = XUNDEF
IF (LUSERI ) XLES_UPDRAFT_Ri2  = XUNDEF
IF (NSV>0  ) XLES_UPDRAFT_Sv2  = XUNDEF
IF (LUSERV ) XLES_UPDRAFT_ThRv = XUNDEF
IF (LUSERC ) XLES_UPDRAFT_ThRc = XUNDEF
IF (LUSERI ) XLES_UPDRAFT_ThRi = XUNDEF
IF (LUSERC ) XLES_UPDRAFT_ThlRv= XUNDEF
IF (LUSERC ) XLES_UPDRAFT_ThlRc= XUNDEF
IF (LUSERI ) XLES_UPDRAFT_ThlRi= XUNDEF
IF (NSV>0  ) XLES_UPDRAFT_ThSv = XUNDEF
IF (LUSERV ) XLES_UPDRAFT_ThvRv= XUNDEF
IF (LUSERC ) XLES_UPDRAFT_ThvRc= XUNDEF
IF (LUSERI ) XLES_UPDRAFT_ThvRi= XUNDEF
IF (NSV>0 .AND. LUSERV) XLES_UPDRAFT_ThvSv = XUNDEF
IF (NSV>0 .AND. LUSERC) XLES_UPDRAFT_ThlSv = XUNDEF
!
!
!*      7.7  downdraft quantities (only on the cartesian mask)
!            --------------------
!
ALLOCATE(XLES_DOWNDRAFT       (NLES_K,NLES_TIMES))    ! updraft fraction
ALLOCATE(XLES_DOWNDRAFT_W     (NLES_K,NLES_TIMES))    ! <w>
ALLOCATE(XLES_DOWNDRAFT_Th    (NLES_K,NLES_TIMES))    ! <theta>
ALLOCATE(XLES_DOWNDRAFT_Ke    (NLES_K,NLES_TIMES))    ! <E>
ALLOCATE(XLES_DOWNDRAFT_WTh   (NLES_K,NLES_TIMES))    ! <w'theta'>
ALLOCATE(XLES_DOWNDRAFT_Th2   (NLES_K,NLES_TIMES))    ! <th'2>
ALLOCATE(XLES_DOWNDRAFT_Tke   (NLES_K,NLES_TIMES))    ! <e>

IF (LUSERV) THEN
  ALLOCATE(XLES_DOWNDRAFT_Thv   (NLES_K,NLES_TIMES))    ! <thetav>
  ALLOCATE(XLES_DOWNDRAFT_WThv  (NLES_K,NLES_TIMES))    ! <w'thv'>
  ALLOCATE(XLES_DOWNDRAFT_ThThv (NLES_K,NLES_TIMES))    ! <th'thv'>
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Thv   (0,0))
  ALLOCATE(XLES_DOWNDRAFT_WThv  (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThThv (0,0))
END IF
!
IF (LUSERC) THEN
  ALLOCATE(XLES_DOWNDRAFT_Thl   (NLES_K,NLES_TIMES))    ! <thetal>
  ALLOCATE(XLES_DOWNDRAFT_WThl  (NLES_K,NLES_TIMES))    ! <w'thetal'>
  ALLOCATE(XLES_DOWNDRAFT_Thl2  (NLES_K,NLES_TIMES))    ! <thl'2>
  ALLOCATE(XLES_DOWNDRAFT_ThlThv(NLES_K,NLES_TIMES))    ! <thl'thv'>
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Thl   (0,0))
  ALLOCATE(XLES_DOWNDRAFT_WThl  (0,0))
  ALLOCATE(XLES_DOWNDRAFT_Thl2  (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThlThv(0,0))
END IF

IF (LUSERV ) THEN
  ALLOCATE(XLES_DOWNDRAFT_Rv    (NLES_K,NLES_TIMES))    ! <Rv>
  ALLOCATE(XLES_DOWNDRAFT_WRv   (NLES_K,NLES_TIMES))    ! <w'Rv'>
  ALLOCATE(XLES_DOWNDRAFT_Rv2   (NLES_K,NLES_TIMES))    ! <Rv'2>
  ALLOCATE(XLES_DOWNDRAFT_ThRv  (NLES_K,NLES_TIMES))    ! <Th'Rv'>
  ALLOCATE(XLES_DOWNDRAFT_ThvRv (NLES_K,NLES_TIMES))    ! <Thv'Rv'>
  IF (LUSERC) THEN
    ALLOCATE(XLES_DOWNDRAFT_ThlRv (NLES_K,NLES_TIMES))    ! <Thl'Rv'>
  ELSE
    ALLOCATE(XLES_DOWNDRAFT_ThlRv (0,0))
  END IF
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Rv    (0,0))
  ALLOCATE(XLES_DOWNDRAFT_WRv   (0,0))
  ALLOCATE(XLES_DOWNDRAFT_Rv2   (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThRv  (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThvRv (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThlRv (0,0))
END IF
IF (LUSERC ) THEN
  ALLOCATE(XLES_DOWNDRAFT_Rc    (NLES_K,NLES_TIMES))    ! <Rc>
  ALLOCATE(XLES_DOWNDRAFT_WRc   (NLES_K,NLES_TIMES))    ! <w'Rc'>
  ALLOCATE(XLES_DOWNDRAFT_Rc2   (NLES_K,NLES_TIMES))    ! <Rc'2>
  ALLOCATE(XLES_DOWNDRAFT_ThRc  (NLES_K,NLES_TIMES))    ! <Th'Rc'>
  ALLOCATE(XLES_DOWNDRAFT_ThvRc (NLES_K,NLES_TIMES))    ! <Thv'Rc'>
  ALLOCATE(XLES_DOWNDRAFT_ThlRc (NLES_K,NLES_TIMES))    ! <Thl'Rc'>
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Rc    (0,0))
  ALLOCATE(XLES_DOWNDRAFT_WRc   (0,0))
  ALLOCATE(XLES_DOWNDRAFT_Rc2   (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThRc  (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThvRc (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThlRc (0,0))
END IF
IF (LUSERR ) THEN
  ALLOCATE(XLES_DOWNDRAFT_Rr    (NLES_K,NLES_TIMES))    ! <Rr>
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Rr    (0,0))
END IF
IF (LUSERI ) THEN
  ALLOCATE(XLES_DOWNDRAFT_Ri    (NLES_K,NLES_TIMES))    ! <Ri>
  ALLOCATE(XLES_DOWNDRAFT_WRi   (NLES_K,NLES_TIMES))    ! <w'Ri'>
  ALLOCATE(XLES_DOWNDRAFT_Ri2   (NLES_K,NLES_TIMES))    ! <Ri'2>
  ALLOCATE(XLES_DOWNDRAFT_ThRi  (NLES_K,NLES_TIMES))    ! <Th'Ri'>
  ALLOCATE(XLES_DOWNDRAFT_ThvRi (NLES_K,NLES_TIMES))    ! <Thv'Ri'>
  ALLOCATE(XLES_DOWNDRAFT_ThlRi (NLES_K,NLES_TIMES))    ! <Thl'Ri'>
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Ri    (0,0))
  ALLOCATE(XLES_DOWNDRAFT_WRi   (0,0))
  ALLOCATE(XLES_DOWNDRAFT_Ri2   (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThRi  (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThvRi (0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThlRi (0,0))
END IF
IF (LUSERS ) THEN
  ALLOCATE(XLES_DOWNDRAFT_Rs    (NLES_K,NLES_TIMES))    ! <Rs>
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Rs    (0,0))
END IF
IF (LUSERG ) THEN
  ALLOCATE(XLES_DOWNDRAFT_Rg    (NLES_K,NLES_TIMES))    ! <Rg>
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Rg    (0,0))
END IF
IF (LUSERH ) THEN
  ALLOCATE(XLES_DOWNDRAFT_Rh    (NLES_K,NLES_TIMES))    ! <Rh>
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Rh    (0,0))
END IF
IF (NSV>0  ) THEN
  ALLOCATE(XLES_DOWNDRAFT_Sv    (NLES_K,NLES_TIMES,NSV))! <Sv>
  ALLOCATE(XLES_DOWNDRAFT_WSv   (NLES_K,NLES_TIMES,NSV))! <w'Sv'>
  ALLOCATE(XLES_DOWNDRAFT_Sv2   (NLES_K,NLES_TIMES,NSV))! <Sv'2>
  ALLOCATE(XLES_DOWNDRAFT_ThSv  (NLES_K,NLES_TIMES,NSV))! <Th'Sv'>
  IF (LUSERV) THEN
    ALLOCATE(XLES_DOWNDRAFT_ThvSv  (NLES_K,NLES_TIMES,NSV))! <Thv'Sv'>
  ELSE
    ALLOCATE(XLES_DOWNDRAFT_ThvSv  (0,0,0))
  END IF
  IF (LUSERC) THEN
    ALLOCATE(XLES_DOWNDRAFT_ThlSv  (NLES_K,NLES_TIMES,NSV))! <Thl'Sv'>
  ELSE
    ALLOCATE(XLES_DOWNDRAFT_ThlSv  (0,0,0))
  END IF
ELSE
  ALLOCATE(XLES_DOWNDRAFT_Sv    (0,0,0))
  ALLOCATE(XLES_DOWNDRAFT_WSv   (0,0,0))
  ALLOCATE(XLES_DOWNDRAFT_Sv2   (0,0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThSv  (0,0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThvSv (0,0,0))
  ALLOCATE(XLES_DOWNDRAFT_ThlSv (0,0,0))
END IF
!
!
XLES_DOWNDRAFT     = XUNDEF
XLES_DOWNDRAFT_W   = XUNDEF
XLES_DOWNDRAFT_Th  = XUNDEF
XLES_DOWNDRAFT_Thl = XUNDEF
XLES_DOWNDRAFT_Tke  = XUNDEF
IF (LUSERV ) XLES_DOWNDRAFT_Thv = XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_Thl = XUNDEF
IF (LUSERV ) XLES_DOWNDRAFT_Rv  = XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_Rc  = XUNDEF
IF (LUSERR ) XLES_DOWNDRAFT_Rr  = XUNDEF
IF (LUSERI ) XLES_DOWNDRAFT_Ri  = XUNDEF
IF (LUSERS ) XLES_DOWNDRAFT_Rs  = XUNDEF
IF (LUSERG ) XLES_DOWNDRAFT_Rg  = XUNDEF
IF (LUSERH ) XLES_DOWNDRAFT_Rh  = XUNDEF
IF (NSV>0  ) XLES_DOWNDRAFT_Sv  = XUNDEF
XLES_DOWNDRAFT_Ke   = XUNDEF
XLES_DOWNDRAFT_WTh  = XUNDEF
IF (LUSERV ) XLES_DOWNDRAFT_WThv = XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_WThl = XUNDEF
IF (LUSERV ) XLES_DOWNDRAFT_WRv  = XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_WRc  = XUNDEF
IF (LUSERI ) XLES_DOWNDRAFT_WRi  = XUNDEF
IF (NSV>0  ) XLES_DOWNDRAFT_WSv  = XUNDEF
XLES_DOWNDRAFT_Th2  = XUNDEF
IF (LUSERV ) THEN
  XLES_DOWNDRAFT_ThThv  = XUNDEF
END IF
IF (LUSERC ) THEN
  XLES_DOWNDRAFT_Thl2   = XUNDEF
  XLES_DOWNDRAFT_ThlThv = XUNDEF
END IF
IF (LUSERV ) XLES_DOWNDRAFT_Rv2  = XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_Rc2  = XUNDEF
IF (LUSERI ) XLES_DOWNDRAFT_Ri2  = XUNDEF
IF (NSV>0  ) XLES_DOWNDRAFT_Sv2  = XUNDEF
IF (LUSERV ) XLES_DOWNDRAFT_ThRv = XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_ThRc = XUNDEF
IF (LUSERI ) XLES_DOWNDRAFT_ThRi = XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_ThlRv= XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_ThlRc= XUNDEF
IF (LUSERI ) XLES_DOWNDRAFT_ThlRi= XUNDEF
IF (NSV>0  ) XLES_DOWNDRAFT_ThSv = XUNDEF
IF (LUSERV ) XLES_DOWNDRAFT_ThvRv= XUNDEF
IF (LUSERC ) XLES_DOWNDRAFT_ThvRc= XUNDEF
IF (LUSERI ) XLES_DOWNDRAFT_ThvRi= XUNDEF
IF (NSV>0 .AND. LUSERV) XLES_DOWNDRAFT_ThvSv = XUNDEF
IF (NSV>0 .AND. LUSERC) XLES_DOWNDRAFT_ThlSv = XUNDEF
!
!*      7.8  production terms
!            ----------------
!
ALLOCATE(XLES_BU_RES_KE   (NLES_K,NLES_TIMES,NLES_TOT))
ALLOCATE(XLES_BU_RES_WThl (NLES_K,NLES_TIMES,NLES_TOT))
ALLOCATE(XLES_BU_RES_Thl2 (NLES_K,NLES_TIMES,NLES_TOT))
ALLOCATE(XLES_BU_SBG_TKE  (NLES_K,NLES_TIMES,NLES_TOT))
XLES_BU_RES_KE    = 0.
XLES_BU_RES_WThl  = 0.
XLES_BU_RES_Thl2  = 0.
XLES_BU_SBG_TKE   = 0.
IF (LUSERV) THEN
  ALLOCATE(XLES_BU_RES_WRt  (NLES_K,NLES_TIMES,NLES_TOT))
  ALLOCATE(XLES_BU_RES_Rt2  (NLES_K,NLES_TIMES,NLES_TOT))
  ALLOCATE(XLES_BU_RES_ThlRt(NLES_K,NLES_TIMES,NLES_TOT))
  XLES_BU_RES_WRt   = 0.
  XLES_BU_RES_Rt2   = 0.
  XLES_BU_RES_ThlRt = 0.
END IF
ALLOCATE(XLES_BU_RES_WSv (NLES_K,NLES_TIMES,NLES_TOT,NSV))
ALLOCATE(XLES_BU_RES_Sv2 (NLES_K,NLES_TIMES,NLES_TOT,NSV))
IF (NSV>0) THEN
  XLES_BU_RES_WSv = 0.
  XLES_BU_RES_Sv2 = 0.
END IF
!
!-------------------------------------------------------------------------------
!
!*      8.   Allocations of the normalization variables temporal series
!            ----------------------------------------------------------
!
ALLOCATE(XLES_UW0       (NLES_TIMES))
ALLOCATE(XLES_VW0       (NLES_TIMES))
ALLOCATE(XLES_USTAR     (NLES_TIMES))
ALLOCATE(XLES_WSTAR     (NLES_TIMES))
ALLOCATE(XLES_Q0        (NLES_TIMES))
ALLOCATE(XLES_E0        (NLES_TIMES))
ALLOCATE(XLES_SV0       (NLES_TIMES,NSV))
ALLOCATE(XLES_BL_HEIGHT (NLES_TIMES))
ALLOCATE(XLES_MO_LENGTH (NLES_TIMES))
ALLOCATE(XLES_ZCB       (NLES_TIMES))
ALLOCATE(XLES_CFtot     (NLES_TIMES))
ALLOCATE(XLES_CF2tot    (NLES_TIMES))
ALLOCATE(XLES_LWP       (NLES_TIMES))
ALLOCATE(XLES_LWPVAR    (NLES_TIMES))
ALLOCATE(XLES_RWP       (NLES_TIMES))
ALLOCATE(XLES_IWP       (NLES_TIMES))
ALLOCATE(XLES_SWP       (NLES_TIMES))
ALLOCATE(XLES_GWP       (NLES_TIMES))
ALLOCATE(XLES_HWP       (NLES_TIMES))
ALLOCATE(XLES_INT_TKE   (NLES_TIMES))
ALLOCATE(XLES_ZMAXCF    (NLES_TIMES))
ALLOCATE(XLES_ZMAXCF2   (NLES_TIMES))  
ALLOCATE(XLES_INPRR     (NLES_TIMES))
ALLOCATE(XLES_INPRC     (NLES_TIMES))
ALLOCATE(XLES_INDEP     (NLES_TIMES))
ALLOCATE(XLES_RAIN_INPRR(NLES_TIMES))
ALLOCATE(XLES_ACPRR     (NLES_TIMES))
ALLOCATE(XLES_PRECFR    (NLES_TIMES))
ALLOCATE(XLES_SWU       (NLES_K,NLES_TIMES))
ALLOCATE(XLES_SWD       (NLES_K,NLES_TIMES))
ALLOCATE(XLES_LWU       (NLES_K,NLES_TIMES))
ALLOCATE(XLES_LWD       (NLES_K,NLES_TIMES))
ALLOCATE(XLES_DTHRADSW  (NLES_K,NLES_TIMES))
ALLOCATE(XLES_DTHRADLW  (NLES_K,NLES_TIMES))
ALLOCATE(XLES_RADEFF    (NLES_K,NLES_TIMES))
!
XLES_UW0       = XUNDEF
XLES_VW0       = XUNDEF
XLES_USTAR     = XUNDEF
XLES_WSTAR     = XUNDEF
XLES_Q0        = XUNDEF
XLES_E0        = XUNDEF
XLES_SV0       = XUNDEF
XLES_BL_HEIGHT = XUNDEF
XLES_MO_LENGTH = XUNDEF
XLES_ZCB       = XUNDEF
XLES_CFtot     = XUNDEF
XLES_CF2tot    = XUNDEF  
XLES_LWP       = XUNDEF
XLES_LWPVAR    = XUNDEF
XLES_RWP       = XUNDEF
XLES_IWP       = XUNDEF
XLES_SWP       = XUNDEF
XLES_GWP       = XUNDEF
XLES_HWP       = XUNDEF
XLES_INT_TKE   = XUNDEF
XLES_ZMAXCF    = XUNDEF
XLES_ZMAXCF2   = XUNDEF  
XLES_PRECFR    = XUNDEF
XLES_ACPRR     = XUNDEF
XLES_INPRR     = XUNDEF
XLES_INPRC     = XUNDEF
XLES_INDEP     = XUNDEF
XLES_RAIN_INPRR = XUNDEF
XLES_SWU        = XUNDEF
XLES_SWD        = XUNDEF
XLES_LWU        = XUNDEF
XLES_LWD        = XUNDEF
XLES_DTHRADSW   = XUNDEF
XLES_DTHRADLW   = XUNDEF
XLES_RADEFF     = XUNDEF
!
!-------------------------------------------------------------------------------
!
!*      9.   Allocations of the normalization variables temporal series
!            ----------------------------------------------------------
!
!       9.1  Two-points correlations in I direction
!            --------------------------------------
!
ALLOCATE(XCORRi_UU    (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between u and u
ALLOCATE(XCORRi_VV    (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between v and v
ALLOCATE(XCORRi_UV    (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between u and v
ALLOCATE(XCORRi_WU    (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between w and u
ALLOCATE(XCORRi_WV    (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between w and v
ALLOCATE(XCORRi_WW    (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between w and w
ALLOCATE(XCORRi_WTh   (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between w and theta
ALLOCATE(XCORRi_ThTh  (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between theta and theta
IF (LUSERC) THEN
  ALLOCATE(XCORRi_WThl  (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))    ! between w and thetal
  ALLOCATE(XCORRi_ThlThl(NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))   ! between thetal and thetal
ELSE
  ALLOCATE(XCORRi_WThl  (0,0,0))
  ALLOCATE(XCORRi_ThlThl(0,0,0))
END IF


IF (LUSERV ) THEN
  ALLOCATE(XCORRi_WRv  (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between w and Rv
  ALLOCATE(XCORRi_ThRv (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between theta and Rv
  IF (LUSERC) THEN
    ALLOCATE(XCORRi_ThlRv(NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between thetal and Rv
  ELSE
    ALLOCATE(XCORRi_ThlRv(0,0,0))
  END IF
  ALLOCATE(XCORRi_RvRv (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between Rv and Rv
ELSE
  ALLOCATE(XCORRi_WRv  (0,0,0))
  ALLOCATE(XCORRi_ThRv (0,0,0))
  ALLOCATE(XCORRi_ThlRv(0,0,0))
  ALLOCATE(XCORRi_RvRv (0,0,0))
END IF

IF (LUSERC ) THEN
  ALLOCATE(XCORRi_WRc  (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between w and Rc
  ALLOCATE(XCORRi_ThRc (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between theta and Rc
  ALLOCATE(XCORRi_ThlRc(NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between thetal and Rc
  ALLOCATE(XCORRi_RcRc (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between Rc and Rc
ELSE
  ALLOCATE(XCORRi_WRc  (0,0,0))
  ALLOCATE(XCORRi_ThRc (0,0,0))
  ALLOCATE(XCORRi_ThlRc(0,0,0))
  ALLOCATE(XCORRi_RcRc (0,0,0))
END IF

IF (LUSERI ) THEN
  ALLOCATE(XCORRi_WRi  (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between w and Ri
  ALLOCATE(XCORRi_ThRi (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between theta and Rc
  ALLOCATE(XCORRi_ThlRi(NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between thetal and Rc
  ALLOCATE(XCORRi_RiRi (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES))     ! between Rc and Rc
ELSE
  ALLOCATE(XCORRi_WRi  (0,0,0))
  ALLOCATE(XCORRi_ThRi (0,0,0))
  ALLOCATE(XCORRi_ThlRi(0,0,0))
  ALLOCATE(XCORRi_RiRi (0,0,0))
END IF

IF (NSV>0  ) THEN
  ALLOCATE(XCORRi_WSv  (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES,NSV)) ! between w and Sv
  ALLOCATE(XCORRi_SvSv (NSPECTRA_NI,NSPECTRA_K,NLES_TIMES,NSV)) ! between Sv and Sv
ELSE
  ALLOCATE(XCORRi_WSv  (0,0,0,0))
  ALLOCATE(XCORRi_SvSv (0,0,0,0))
END IF
!
!
XCORRi_UU  = XUNDEF
XCORRi_VV  = XUNDEF
XCORRi_UV  = XUNDEF
XCORRi_WU  = XUNDEF
XCORRi_WV  = XUNDEF
XCORRi_WW  = XUNDEF
XCORRi_WTh = XUNDEF
IF (LUSERC ) XCORRi_WThl= XUNDEF
IF (LUSERV ) XCORRi_WRv = XUNDEF
IF (LUSERC ) XCORRi_WRc = XUNDEF
IF (LUSERI ) XCORRi_WRi = XUNDEF
IF (NSV>0  ) XCORRi_WSv = XUNDEF
XCORRi_ThTh  = XUNDEF
IF (LUSERC ) XCORRi_ThlThl= XUNDEF
IF (LUSERV ) XCORRi_ThRv = XUNDEF
IF (LUSERC ) XCORRi_ThRc = XUNDEF
IF (LUSERI ) XCORRi_ThRi = XUNDEF
IF (LUSERC ) XCORRi_ThlRv= XUNDEF
IF (LUSERC ) XCORRi_ThlRc= XUNDEF
IF (LUSERI ) XCORRi_ThlRi= XUNDEF
IF (LUSERV ) XCORRi_RvRv = XUNDEF
IF (LUSERC ) XCORRi_RcRc = XUNDEF
IF (LUSERI ) XCORRi_RiRi = XUNDEF
IF (NSV>0  ) XCORRi_SvSv = XUNDEF
!
!
!       9.2  Two-points correlations in J direction
!            --------------------------------------
!
ALLOCATE(XCORRj_UU    (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between u and u
ALLOCATE(XCORRj_VV    (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between v and v
ALLOCATE(XCORRj_UV    (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between u and v
ALLOCATE(XCORRj_WU    (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between w and u
ALLOCATE(XCORRj_WV    (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between w and v
ALLOCATE(XCORRj_WW    (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between w and w
ALLOCATE(XCORRj_WTh   (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between w and theta
ALLOCATE(XCORRj_ThTh  (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between theta and theta
IF (LUSERC) THEN
  ALLOCATE(XCORRj_WThl  (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))    ! between w and thetal
  ALLOCATE(XCORRj_ThlThl(NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))   ! between thetal and thetal
ELSE
  ALLOCATE(XCORRj_WThl  (0,0,0))
  ALLOCATE(XCORRj_ThlThl(0,0,0))
END IF

IF (LUSERV ) THEN
  ALLOCATE(XCORRj_WRv  (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between w and Rv
  ALLOCATE(XCORRj_ThRv (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between theta and Rv
  IF (LUSERC) THEN
    ALLOCATE(XCORRj_ThlRv(NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between thetal and Rv
  ELSE
    ALLOCATE(XCORRj_ThlRv(0,0,0))
  END IF
  ALLOCATE(XCORRj_RvRv (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between Rv and Rv
ELSE
  ALLOCATE(XCORRj_WRv  (0,0,0))
  ALLOCATE(XCORRj_ThRv (0,0,0))
  ALLOCATE(XCORRj_ThlRv(0,0,0))
  ALLOCATE(XCORRj_RvRv (0,0,0))
END IF

IF (LUSERC ) THEN
  ALLOCATE(XCORRj_WRc  (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between w and Rc
  ALLOCATE(XCORRj_ThRc (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between theta and Rc
  ALLOCATE(XCORRj_ThlRc(NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between thetal and Rc
  ALLOCATE(XCORRj_RcRc (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between Rc and Rc
ELSE
  ALLOCATE(XCORRj_WRc  (0,0,0))
  ALLOCATE(XCORRj_ThRc (0,0,0))
  ALLOCATE(XCORRj_ThlRc(0,0,0))
  ALLOCATE(XCORRj_RcRc (0,0,0))
END IF

IF (LUSERI ) THEN
  ALLOCATE(XCORRj_WRi  (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between w and Ri
  ALLOCATE(XCORRj_ThRi (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between theta and Rc
  ALLOCATE(XCORRj_ThlRi(NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between thetal and Rc
  ALLOCATE(XCORRj_RiRi (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES))     ! between Rc and Rc
ELSE
  ALLOCATE(XCORRj_WRi  (0,0,0))
  ALLOCATE(XCORRj_ThRi (0,0,0))
  ALLOCATE(XCORRj_ThlRi(0,0,0))
  ALLOCATE(XCORRj_RiRi (0,0,0))
END IF

IF (NSV>0  ) THEN
  ALLOCATE(XCORRj_WSv  (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES,NSV)) ! between w and Sv
  ALLOCATE(XCORRj_SvSv (NSPECTRA_NJ,NSPECTRA_K,NLES_TIMES,NSV)) ! between Sv and Sv
ELSE
  ALLOCATE(XCORRj_WSv  (0,0,0,0))
  ALLOCATE(XCORRj_SvSv (0,0,0,0))
END IF
!
!
XCORRj_UU  = XUNDEF
XCORRj_VV  = XUNDEF
XCORRj_UV  = XUNDEF
XCORRj_WU  = XUNDEF
XCORRj_WV  = XUNDEF
XCORRj_WW  = XUNDEF
XCORRj_WTh = XUNDEF
IF (LUSERC ) XCORRj_WThl= XUNDEF
IF (LUSERV ) XCORRj_WRv = XUNDEF
IF (LUSERC ) XCORRj_WRc = XUNDEF
IF (LUSERI ) XCORRj_WRi = XUNDEF
IF (NSV>0  ) XCORRj_WSv = XUNDEF
XCORRj_ThTh  = XUNDEF
IF (LUSERC ) XCORRj_ThlThl= XUNDEF
IF (LUSERV ) XCORRj_ThRv = XUNDEF
IF (LUSERC ) XCORRj_ThRc = XUNDEF
IF (LUSERI ) XCORRj_ThRi = XUNDEF
IF (LUSERC ) XCORRj_ThlRv= XUNDEF
IF (LUSERC ) XCORRj_ThlRc= XUNDEF
IF (LUSERI ) XCORRj_ThlRi= XUNDEF
IF (LUSERV ) XCORRj_RvRv = XUNDEF
IF (LUSERC ) XCORRj_RcRc = XUNDEF
IF (LUSERI ) XCORRj_RiRi = XUNDEF
IF (NSV>0  ) XCORRj_SvSv = XUNDEF
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE INI_LES_n
