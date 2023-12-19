!MNH_LIC Copyright 1994-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############
      SUBROUTINE INIT_MNH
!     ###############
!
!!****  *INIT_MNH * - monitor to initialize the variables of the model
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to initialize all the variables
!     used in the  model temporal loop or in the post-processings
!
!!**  METHOD
!!    ------
!!      This initialization  is separated in three parts :
!!         1. A part common to all models where :
!!            - The output-listing file common to all models is opened.
!!            - The physical constants are initialized.
!!            - The other constants for all models are initialized.
!!         2. The treatment of descriptor files model by model :
!!             The DESFM and EXSEG files are read and the EXSEG file is updated
!!         3. The sequential initialization of  nested models :
!!             The initial data fields are read in different files for each
!!             model and  variables which are not in these initial files are
!!             deduced.
!!
!!
!!    EXTERNAL
!!    --------
!!      INI_CST    : to initialize physical constants
!!      INI_CTURB  : to initialize for all models the constants used in the
!!                   turbulence scheme
!!      INI_SEG_n  : to read and update descriptor files  
!!      INI_SIZE   : to initialize the sizes of the different models
!!      INI_MODEL  : to initialize each nested model
!!      INI_PARA_ll: to build the ll data structures
!!      GO_TOMODEL : displace the ll lists to the right nested model
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS : JPMODELMAX
!!
!!      Module MODD_CONF       : NMODEL,NVERB
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (routine INIT_MNH)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    02/06/94
!!      J.Stein     05/01/95  add ini_cturb
!!      J.P. Lafore 18/08/95  Time STEP change
!!      J.P. Lafore 22/07/96  ZTSTEP_ALL introduction for nesting
!!      V. Ducrocq  7/08/98   //
!!      P. Jabouille 7/07/99  split ini_modeln in 2 parts+ cleaning
!!      V. Masson   15/03/99  call to ini_data_cover
!!      P.Jabouille 15/07/99  special initialisation for spawning
!!      J.P Chaboureau 2015   add ini_spectre_n
!!      J.Escobar   2/03/2016 bypass , reset NHALO=1 for SPAWNING
!!  06/2016     (G.Delautier) phasage surfex 8
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_CONF
USE MODD_DYN_n, ONLY: CPRESOPT, NITR ! only for spawning purpose
USE MODD_IO,    ONLY: TFILE_OUTPUTLISTING, TPTR2FILE
USE MODD_LBC_n, ONLY: CLBCX,CLBCY   ! only for spawning purpose
USE MODD_LUNIT
USE MODD_LUNIT_n
USE MODD_MNH_SURFEX_n
USE MODD_NSV,          ONLY: NSV_ASSOCIATE
USE MODD_PARAMETERS
!
use mode_field,            only: Alloc_field_scalars, Fieldlist_goto_model
USE MODE_INI_CST,          ONLY: INI_CST
USE MODE_IO_FILE,          ONLY: IO_File_open
USE MODE_IO_MANAGE_STRUCT, ONLY: IO_File_add2list
USE MODE_ll
USE MODE_MODELN_HANDLER
USE MODE_SPLITTINGZ_ll
!
USE MODI_INI_MODEL_n
USE MODI_INI_SEG_n
USE MODI_INI_SIZE_n
USE MODI_INI_SIZE_SPAWN
USE MODI_INI_SPECTRE_n
USE MODI_READ_ALL_NAMELISTS
USE MODI_RESET_EXSEG
!
IMPLICIT NONE
!
!*       0.1   Local variables
!
INTEGER :: JMI                                        !  Loop index
CHARACTER(LEN=28),DIMENSION(JPMODELMAX)  :: YINIFILEPGD
INTEGER  :: ILUOUT0,IRESP                             ! Logical unit number for
                                                      ! output-listing common
                                                      ! to all models and return
                                                      ! code of file management
REAL, DIMENSION(JPMODELMAX)            :: ZTSTEP_ALL  ! Time STEP of ALL models
INTEGER                                :: IINFO_ll    ! return code of // routines
!
! Dummy pointers needed to correct an ifort Bug
CHARACTER(LEN=4), DIMENSION(:), POINTER :: DPTR_CLBCX,DPTR_CLBCY

!-------------------------------------------------------------------------------
!
!*       1.    INITIALIZATION  COMMON TO ALL MODELS
!              ------------------------------------
!
!*       1.1   initialize // E/S and open  output-listing file
!
!
IF (CPROGRAM/='REAL  ') THEN
  CALL IO_File_add2list(TLUOUT0,'OUTPUT_LISTING0','OUTPUTLISTING','WRITE')
  CALL IO_File_open(TLUOUT0)
  !Set output file for PRINT_MSG
  TFILE_OUTPUTLISTING => TLUOUT0
  ILUOUT0=TLUOUT0%NLU
ELSE
  ILUOUT0=TLUOUT0%NLU
END IF
!
WRITE(UNIT=ILUOUT0,FMT="(50('*'),/,'*',48X,'*',/,                  &
                  &  7('*'),10X, ' MESO-NH MODEL ',10X,8('*'),/,   &
                  & '*',48X,'*',/,                                 &
                  & 7('*'),12X,' CNRM - LA ',12X,8('*'),/,         &
                  & '*',48X,'*',/,  50('*'))")
!
CALL NSV_ASSOCIATE()
!
!
!*      1.2   initialize physical constants
!
CALL INI_CST
!
!
!*       1.3    initialize constants for the turbulence scheme                      
!
!Now done in ini_modeln
!
!
!-------------------------------------------------------------------------------
!
!*       2.    READ AND UPDATE DESCRIPTOR FILES
!              --------------------------------
!
IF (CPROGRAM=='SPAWN ' .OR. CPROGRAM=='DIAG  ' .OR. CPROGRAM=='SPEC  ' .OR. CPROGRAM=='MESONH') THEN
  CALL ALLOC_FIELD_SCALARS()
END IF
!
CALL GOTO_MODEL(1)
CALL INI_SEG_n(1,LUNIT_MODEL(1)%TINIFILE,YINIFILEPGD(1),ZTSTEP_ALL)
!
DO JMI=2,NMODEL
  CALL GOTO_MODEL(JMI)
  CALL INI_SEG_n(JMI,LUNIT_MODEL(JMI)%TINIFILE,YINIFILEPGD(JMI),ZTSTEP_ALL)
END DO
!
IF (CPROGRAM=='SPAWN ') THEN 
  !bypass
  NHALO = 1
END IF
!
IF (CPROGRAM=='DIAG') CALL RESET_EXSEG()
!
!-------------------------------------------------------------------------------
!
!
!*       3.    INITIALIZE EACH MODEL SIZES AND DEPENDENCY
!              ------------------------------------------
!
DO JMI=1,NMODEL
  CALL GOTO_MODEL(JMI)
  CALL INI_SIZE_n(JMI,LUNIT_MODEL(JMI)%TINIFILE,YINIFILEPGD(JMI))
END DO
!
IF (CPROGRAM=='SPAWN ') THEN 
  DPTR_CLBCX=>CLBCX
  DPTR_CLBCY=>CLBCY
  CALL INI_PARAZ_ll(IINFO_ll)
  CALL INI_SIZE_SPAWN(DPTR_CLBCX,DPTR_CLBCY,CPRESOPT,NITR,LUNIT_MODEL(1)%TINIFILE)
END IF
!
!   INITIALIZE data structures of ComLib
!
!JUAN CALL INI_PARA_ll(IINFO_ll)
CALL INI_PARAZ_ll(IINFO_ll)
!
!-------------------------------------------------------------------------------
!
!
!     Allocations of Surfex Types
CALL SURFEX_ALLOC_LIST(NMODEL)
!
DO JMI=1,NMODEL
  YSURF_CUR => YSURF_LIST(JMI)
!
  IF (CPROGRAM=='SPAWN ' .OR. CPROGRAM=='REAL  ') THEN 
    CALL READ_ALL_NAMELISTS(YSURF_CUR,'MESONH','PRE',.FALSE.)
  ELSE
    CALL READ_ALL_NAMELISTS(YSURF_CUR,'MESONH','ALL',.TRUE.)
  ENDIF
ENDDO
!
!
!-------------------------------------------------------------------------------
!
!*       4.    INITIALIZE EACH MODEL
!              ---------------------
!
DO JMI=1,NMODEL
  CALL GO_TOMODEL_ll(JMI,IINFO_ll)
  CALL GOTO_MODEL(JMI)
  IF (CPROGRAM/='SPEC  ') THEN
    CALL INI_MODEL_n(JMI,LUNIT_MODEL(JMI)%TINIFILE)
    !Call necessary to update the TFIELDLIST pointers to the data
    CALL FIELDLIST_GOTO_MODEL(JMI,JMI)
  ELSE
    CALL INI_SPECTRE_n(JMI,LUNIT_MODEL(JMI)%TINIFILE)
  END IF  
END DO
!
!-------------------------------------------------------------------------------
!
!*       5.    WRITE MESSAGE ON OUTPUT-LISTING
!              -------------------------------
!
IF (NVERB >= 5) THEN
  WRITE(UNIT=ILUOUT0,FMT="(50('*'),/,'*',48X,'*',/,                &
                & '*',10X,' INITIALIZATION TERMINATED',10X,'*',/,  &
                & '*',48X,'*',/,50('*'))")
END IF
!
!-------------------------------------------------------------------------------
!
!
END SUBROUTINE INIT_MNH
