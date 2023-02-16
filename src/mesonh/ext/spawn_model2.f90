!MNH_LIC Copyright 1995-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!########################
MODULE MODI_SPAWN_MODEL2
!########################
!
INTERFACE
!
      SUBROUTINE SPAWN_MODEL2 (KRR,KSV_USER,HTURB,HSURF,HCLOUD,    &
                               HCHEM_INPUT_FILE,HSPAFILE,HSPANBR,  &
                               HSONFILE,HINIFILE,HINIFILEPGD,OSPAWN_SURF       )
!
INTEGER,               INTENT(IN)  :: KRR         ! Number of moist variables
INTEGER,               INTENT(IN)  :: KSV_USER    ! Number of Users Scalar Variables
CHARACTER (LEN=4),     INTENT(IN)  :: HTURB       ! Kind of turbulence parameterization
CHARACTER (LEN=4),     INTENT(IN)  :: HSURF       ! Kind of surface parameterization
CHARACTER (LEN=4),     INTENT(IN)  :: HCLOUD      ! Kind of cloud parameterization
                                                  ! model 2 physical domain
CHARACTER (LEN=*),     INTENT(IN) :: HSPAFILE     ! possible name of the output FM-file
CHARACTER (LEN=*),     INTENT(IN) :: HSPANBR      ! NumBeR associated to the SPAwned file
CHARACTER (LEN=*),     INTENT(IN) :: HSONFILE     ! name of the input FM-file SON
CHARACTER (LEN=80),    INTENT(IN) :: HCHEM_INPUT_FILE
CHARACTER (LEN=*),     INTENT(IN) :: HINIFILE     ! Input file
CHARACTER (LEN=*),     INTENT(IN) :: HINIFILEPGD  ! Input pgd file
LOGICAL,               INTENT(IN) :: OSPAWN_SURF  ! flag to spawn surface fields
!
END SUBROUTINE SPAWN_MODEL2
!
END INTERFACE
!
END MODULE MODI_SPAWN_MODEL2
!     ######spl
      SUBROUTINE SPAWN_MODEL2 (KRR,KSV_USER,HTURB,HSURF,HCLOUD,    &
                               HCHEM_INPUT_FILE,HSPAFILE,HSPANBR,  &
                               HSONFILE,HINIFILE,HINIFILEPGD,OSPAWN_SURF       )
!     #######################################################################
!
!!****  *SPAWN_MODEL2 * - subroutine to prepare by horizontal interpolation and
!!                        write an initial FM-file spawned from an other FM-file.
!!
!!    PURPOSE
!!    -------
!!
!!      Initializes by horizontal interpolation, the model 2 in a sub-domain of 
!!    model 1,  possibly overwrites model 2 information by model SON1,
!!    and writes the resulting fields in a FM-file.
!!
!!
!!**  METHOD
!!    ------
!!
!!      In this routine, only the model 2 variables are known through the
!!    MODD_... calls.
!!
!!      The directives to perform the preparation of the initial FM
!!    file are stored in EXSPA.nam file.
!!
!!      The following  SPAWN_MODEL2 routine :
!!
!!             - sets default values of DESFM files
!!             - reads the namelists part of EXSPA file which gives the
!!      directives concerning the spawning to perform
!!             - controls the domain size of model 2 and initializes its 
!!      configuration for parameterizations and LBC
!!             - allocates memory for arrays
!!             - computes the interpolation coefficients needed to spawn model 2 
!!      2 types of interpolations are used:
!!                 1. Clark and Farley (JAS 1984) on 9 points 
!!                 2. Bikhardt on 16 points
!!             - initializes fields
!!             - reads SON1 fields and overwrites on common domain
!!             - writes the DESFM file (variables written have been initialized
!!      by reading the DESFM file concerning the model 1)
!!             - writes the LFIFM file. 
!!
!!       Finally some control prints are performed on the output listing.
!!
!!    EXTERNAL
!!    --------
!!
!!      Module MODE_GRIDPROJ : contains conformal projection routines
!!           SM_GRIDPROJ   : to compute some grid variables, in
!!                           case of conformal projection.
!!      Module MODE_GRIDCART : contains cartesian geometry routines
!!           SM_GRIDCART   : to compute some grid variables, in
!!                           case of cartesian geometry.
!!      SET_REF       : to compute  rhoJ 
!!      TOTAL_DMASS   : to compute the total mass of dry air
!!      ANEL_BALANCE2  : to apply an anelastic correction in the case of changing
!!                      resolution between the two models
!!      IO_File_open : to open a FM-file (DESFM + LFIFM)
!!      WRITE_DESFM   : to write the  DESFM file
!!      WRITE_LFIFM   : to write the  LFIFM file  
!!      IO_File_close : to close a FM-file (DESFM + LFIFM)
!!      INI_BIKHARDT2     : initializes Bikhardt coefficients
!!
!!
!! 
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!
!!      Module MODD_PARAMETERS : contains parameters 
!!      Module MODD_CONF       : contains configuration variables for all models
!!      Module MODD_CTURB :
!!         XTKEMIN : mimimum value for the TKE
!!      Module MODD_GRID       : contains grid variables for all models
!!      Module USE MODD_DYN    : contains configuration for the dynamics
!!      Module MODD_REF        : contains reference state variables for
!!                               all models
!!
!!      Module MODD_DIM2       : contains dimensions 
!!      Module MODD_CONF2      : contains configuration variables 
!!      Module MODD_GRID2      : contains grid variables  
!!      Module MODD_TIME2      : contains time variables and uses MODD_TIME
!!      Module MODD_REF2       : contains reference state variables 
!!      Module MODD_FIELD2     : contains prognostic variables
!!      Module MODD_LSFIELD2   : contains Larger Scale fields
!!      Module MODD_GR_FIELD2  : contains surface fields
!!      Module MODD_DYN2       : contains dynamic control variables for model 2 
!!      Module MODD_LBC2       : contains lbc control variables for model 2
!!      Module MODD_PARAM2     : contains configuration for physical parameterizations
!!
!!    REFERENCE
!!    ---------
!!
!!       PROGRAM SPAWN_MODEL2 (Book2 of the documentation)
!!      
!!
!!    AUTHOR
!!    ------
!!
!!       J.P. Lafore     * METEO-FRANCE *
!!
!!    MODIFICATIONS
!!    -------------
!!
!!      Original     11/01/95 
!!      Modification 27/04/95  (I.Mallet) remove R from the historical variables
!!      Modification 16/04/96  (Lafore) Different resolution ratio case introduction
!!      Modification 24/04/96  (Lafore & Masson) Initialization of LUSERWs
!!      Modification 24/04/96  (Masson) Correction of positivity on Rw and TKE
!!      Modification 25/04/96  (Masson) Copies of internal zs on external points
!!      Modification 02/05/96  (Stein Jabouille) initialize CCONF
!!      Modification 31/05/96  (Lafore) Cumputing time analysis
!!      Modification 10/06/96  (Masson) Call to anel_balance in all cases
!!      Modification 10/06/96  (Masson) Bikhardt and Clark_and_Farley coefficients
!!                                      incorporated in modules
!!      Modification 12/06/96  (Masson) default values of NJMAX and KDYRATIO
!!                                      if 2D version of the model
!!      Modification 13/06/96  (Masson) choice of the name of the spawned file
!!      Modification 30/07/96  (Lafore) MY_NAME and DAD_NAME writing for nesting
!!      Modification 25/09/96  (Masson) grid optionnaly given by a fm file
!!                                      and number of points given relatively
!!                                      to model 1
!!      Modification 10/10/96  (Masson) L1D and L2D verifications
!!      Modification 12/11/96  (Masson) allocations of XSRCM and XSRCT
!!      Modification 19/11/96  (Masson) add deep convection
!!      Modification 26/11/96  (Lafore) spawning configuration writing on the FM-file
!!      Modification 26/11/96  (Lafore) replacing of TOTAL_DMASS by REAL_DMASS
!!      Modification 27/02/97  (Lafore) "surfacic" LS fields
!!      Modification 10/04/97  (Lafore) proper treatment of minima
!!      Modification 09/07/97  (Masson) absolute pressure and directional z0
!!      Modification 10/07/97  (Masson) routines SPAWN_PRESSURE2 and DRY_MASS
!!      Modification 17/07/97  (Masson) vertical interpolations and EPS
!!      Modification 29/07/97  (Masson) split mode_lfifm_pgd
!!      Modification 10/08/97  (Lafore) initialization of LUSERV
!!      Modification 14/09/97  (Masson) use of relative humidity
!!      Modification 08/12/97  (Masson) deallocation of model 1 variables
!!      Modification 24/12/97  (Masson) directional z0 parameters and orographies
!!      Modification 20/07/98  (Stein ) add the LB fields
!!      Modification 15/03/99  (Masson) cover types
!!      Modification 15/07/99  (Jabouille) shift domain initialization in INI_SIZE_SPAWN
!!      Modification 04/01/00  (Masson) removes TSZ0 option
!!      Modification 29/11/02  (Pinty)  add C3R5, ICE2, ICE4
!!      Modification 07/07/05  (D.Barbary) spawn with 2 input files (father+son1)
!!      Modification 20/05/06  Remove EPS, Clark and Farley interpolation
!!                             Replace DRY_MASS by TOTAL_DMASS
!!      Modification 06/12  (M.Tomasini) Interpolation of the advective forcing (ADVFRC)
!!                                       and of the turbulent fluxes (EDDY_FLUX)
!!      Modification 07/13  (Bosseur & Filippi) Adds Forefire
!!                   24/04/2014 (J.escobar) bypass CRAY internal compiler error on IIJ computation
!!      Modification 06/2014   (C.Lac) Initialization of physical param of
!!                                      model2 before the call to ini_nsv
!!      Modification 05/02/2015 (M.Moge) parallelization of SPAWNING
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar   02/05/2016 : test ZZS_MAX in // 
!!      P.Wautelet  08/07/2016 : removed MNH_NCWRIT define
!!      J.Escobar   12/07/2016 : add test on NRIMY & change the one on NRIMX with >=
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!                    10/2016 (C.Lac) Add droplet deposition
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 07/02/2019: force TYPE to a known value for IO_File_add2list
!  S. Bielli      02/2019:  sea salt: significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 22/02/2019: replace Hollerith edit descriptor (deleted from Fortran 95 standard)
!  P. Wautelet 14/03/2019: correct ZWS when variable not present in file
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  P. Wautelet 09/03/2021: move some chemistry initializations to ini_nsv
!  P. Wautelet 24/03/2021: bugfix: allocate XLSRVM, XINPAP and XACPAP to zero size when not needed
!!               03/2021 (JL Redelsperger) Ocean model case 
!-------------------------------------------------------------------------------
!
!*       0.     DECLARATIONS
!               ------------
!
USE MODD_PARAMETERS       ! Declarative modules
USE MODD_CST
USE MODD_CONF
USE MODD_CTURB
USE MODD_GRID 
USE MODD_REF
USE MODD_DYN
USE MODD_NESTING
USE MODD_SPAWN
USE MODD_NSV
USE MODD_PASPOL
!
USE MODD_DIM_n
USE MODD_DYN_n
USE MODD_CONF_n 
USE MODD_LBC_n
USE MODD_GRID_n
USE MODD_TIME_n
USE MODD_REF_n
USE MODD_FIELD_n
USE MODD_LSFIELD_n
USE MODD_DUMMY_GR_FIELD_n
USE MODD_PRECIP_n
USE MODD_ELEC_n
USE MODD_LUNIT_n
USE MODD_PARAM_n
USE MODD_TURB_n
USE MODD_METRICS_n
USE MODD_CH_MNHC_n
USE MODD_PASPOL_n
!$20140515
USE MODD_VAR_ll, ONLY : NPROC
USE MODD_IO, ONLY: TFILEDATA,TFILE_DUMMY,TFILE_SURFEX
use modd_precision, only: MNHREAL_MPI
!
USE MODE_GRIDCART         ! Executive modules
USE MODE_GRIDPROJ
USE MODE_ll
USE MODE_MSG
!
USE MODI_READ_HGRID
USE MODI_SPAWN_GRID2  
USE MODI_SPAWN_FIELD2
USE MODI_SPAWN_SURF
USE MODI_VER_INTERP_FIELD
USE MODI_SPAWN_PRESSURE2
USE MODI_SPAWN_SURF2_RAIN
USE MODI_SET_REF
USE MODI_TOTAL_DMASS
USE MODI_ANEL_BALANCE_n
USE MODI_WRITE_DESFM_n
USE MODI_WRITE_LFIFM_n
USE MODI_METRICS
USE MODI_INI_BIKHARDT_n
USE MODI_DEALLOCATE_MODEL1
USE MODI_BOUNDARIES
USE MODI_INI_NSV
!$20140710
USE MODI_UPDATE_METRICS
!
USE MODE_IO_FIELD_READ,    only: IO_Field_read
USE MODE_IO_FIELD_WRITE,   only: IO_Header_write
USE MODE_IO_FILE,          only: IO_File_close, IO_File_open
USE MODE_IO_MANAGE_STRUCT, only: IO_File_add2list
USE MODE_MODELN_HANDLER
USE MODE_MPPDB
!
USE MODE_THERMO
!
USE MODI_SECOND_MNH
!
! Modules for  EDDY_FLUX
USE MODD_LATZ_EDFLX
USE MODD_DEF_EDDY_FLUX_n           
USE MODD_DEF_EDDYUV_FLUX_n
USE MODD_ADVFRC_n
USE MODD_RELFRC_n
USE MODD_2D_FRC
!
!USE MODE_LB_ll, ONLY : SET_LB_FIELD_ll
USE MODI_GET_SIZEX_LB
USE MODI_GET_SIZEY_LB
!
USE MODD_LIMA_PRECIP_SCAVENGING_n
USE MODD_PARAM_LIMA, ONLY : MDEPOC=>LDEPOC, LSCAV
USE MODD_PARAM_ICE,  ONLY : LDEPOSC
USE MODD_PARAM_C2R2, ONLY : LDEPOC
USE MODD_PASPOL, ONLY : LPASPOL
!
USE MODD_MPIF
USE MODD_VAR_ll
use modd_precision, only: LFIINT
!
IMPLICIT NONE
!
!*       0.1.1  Declarations of global variables not declared in the modules :
!
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZJ ! Jacobian
!
!
!*       0.1.2  Declarations of dummy arguments :
!
INTEGER,               INTENT(IN)  :: KRR         ! Number of moist variables
INTEGER,               INTENT(IN)  :: KSV_USER    ! Number of Users Scalar Variables
CHARACTER (LEN=4),     INTENT(IN)  :: HTURB       ! Kind of turbulence parameterization
CHARACTER (LEN=4),     INTENT(IN)  :: HSURF       ! Kind of surface parameterization
CHARACTER (LEN=4),     INTENT(IN)  :: HCLOUD      ! Kind of cloud parameterization
CHARACTER (LEN=*),     INTENT(IN) :: HSPAFILE     ! possible name of the output FM-file
CHARACTER (LEN=*),     INTENT(IN) :: HSPANBR      ! NumBeR associated to the SPAwned file
CHARACTER (LEN=*),     INTENT(IN) :: HSONFILE     ! name of the input FM-file SON
CHARACTER (LEN=80),    INTENT(IN) :: HCHEM_INPUT_FILE
CHARACTER (LEN=*),     INTENT(IN) :: HINIFILE     ! Input file
CHARACTER (LEN=*),     INTENT(IN) :: HINIFILEPGD  ! Input pgd file
LOGICAL,               INTENT(IN) :: OSPAWN_SURF  ! flag to spawn surface fields
!
!*       0.1.3  Declarations of local variables :
!
!
INTEGER :: ILUOUT   ! Logical unit number for the output listing 
INTEGER(KIND=LFIINT) :: INPRAR ! Number of articles predicted in the LFIFM file
!
!
INTEGER             :: IIU            ! Upper dimension in x direction
INTEGER             :: IJU            ! Upper dimension in y direction
INTEGER             :: IKU            ! Upper dimension in z direction
INTEGER             :: IIB            ! indice I Beginning in x direction
INTEGER             :: IJB            ! indice J Beginning in y direction
INTEGER             :: IKB            ! indice K Beginning in z direction
INTEGER             :: IIE            ! indice I End       in x direction 
INTEGER             :: IJE            ! indice J End       in y direction 
INTEGER             :: IKE            ! indice K End       in z direction 
INTEGER             :: JK             ! Loop index in z direction 
INTEGER             :: JLOOP,JKLOOP   ! Loop indexes 
INTEGER             :: JSV            ! loop index for scalar variables
INTEGER             :: JRR            ! loop index for moist variables
!
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZZS_LS ! large scale interpolated zs
REAL, DIMENSION(:,:),   ALLOCATABLE :: ZZSMT_LS ! large scale interpolated smooth zs
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZZZ_LS ! large scale interpolated z
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZTHVT  ! virtual potential temperature
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZHUT   ! relative humidity
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZSUMRT ! sum of water ratios
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZRHOD  ! dry density
!
REAL    :: ZTIME1,ZTIME2,ZSTART,ZEND,ZTOT,ZALL,ZPERCALL ! for computing time analysis
REAL    ::     ZGRID2,    ZSURF2,    ZFIELD2,     ZVER, &
           ZPRESSURE2,    ZANEL,      ZWRITE,     ZMISC
REAL    :: ZPERCGRID2,ZPERCSURF2,ZPERCFIELD2, ZPERCVER, &
       ZPERCPRESSURE2, ZPERCANEL, ZPERCWRITE,ZPERCMISC
!
INTEGER, DIMENSION(2) :: IIJ
INTEGER               :: IK4000
INTEGER               :: IMI ! Old Model index
!
! Spawning variables for the SON 1 (input one)
INTEGER             :: IIMAXSON,IJMAXSON ! physical dimensions
INTEGER             :: IIUSON,IJUSON     ! upper dimensions
INTEGER             :: IXSIZESON,IYSIZESON ! sizes according to model1 grid
INTEGER             :: IDXRATIOSON,IDYRATIOSON ! x and y-resolution ratios
INTEGER             :: IXORSON,IYORSON   ! horizontal position 
INTEGER             :: IXENDSON,IYENDSON !in x and y directions
! Common indexes for the SON 2 (output one, model2)
INTEGER             :: IIB2           ! indice I Beginning in x direction
INTEGER             :: IJB2           ! indice J Beginning in y direction
INTEGER             :: IIE2           ! indice I End       in x direction
INTEGER             :: IJE2           ! indice J End       in y direction
! Common indexes for the SON 1 (input one)
INTEGER             :: IIB1           ! indice I Beginning in x direction
INTEGER             :: IJB1           ! indice J Beginning in y direction
INTEGER             :: IIE1           ! indice I End       in x direction
INTEGER             :: IJE1           ! indice J End       in y direction
! Logical for no common domain between the 2 sons or no input son
LOGICAL             :: GNOSON = .TRUE.
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORK3D ! working array
CHARACTER(LEN=28)   :: YDAD_SON
!$
INTEGER             :: IINFO_ll
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=>NULL()   ! list of fields to exchange
INTEGER             :: NXOR_TMP, NYOR_TMP, NXEND_TMP, NYEND_TMP 
INTEGER :: IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU     ! dimensions of the
INTEGER :: IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2       ! West-east LB arrays
INTEGER :: IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV     ! dimensions of the
INTEGER :: IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2       ! North-south LB arrays
!
CHARACTER(LEN=4)    :: YLBTYPE
!
INTEGER,DIMENSION(:,:),ALLOCATABLE   :: IJCOUNT 
!
REAL                :: ZZS_MAX, ZZS_MAX_ll
!
TYPE(TFILEDATA),POINTER :: TZFILE      => NULL()
TYPE(TFILEDATA),POINTER :: TZSONFILE   => NULL()
!-------------------------------------------------------------------------------
!
! Save model index and switch to model 2 variables
IMI = GET_CURRENT_MODEL_INDEX()
CALL GOTO_MODEL(2)
CSTORAGE_TYPE='TT'
!
ILUOUT=TLUOUT%NLU
!
!*   1.    INITIALIZATIONS :
!           ---------------
!
!*   1.1   time analysis :
!          -------------
!
ZTIME1 = 0
ZTIME2 = 0
ZSTART = 0
ZEND   = 0 
ZGRID2 = 0
ZSURF2 = 0 
ZFIELD2= 0 
ZANEL  = 0 
ZWRITE = 0 
ZPERCGRID2 = 0
ZPERCSURF2 = 0 
ZPERCFIELD2= 0 
ZPERCANEL  = 0 
ZPERCWRITE = 0 
!
CALL SECOND_MNH(ZSTART)
!
ZTIME1 = ZSTART
!
!*	 1.2   deallocates not used model 1 variables :  
!              --------------------------------------
!
CALL DEALLOCATE_MODEL1(1)
CALL DEALLOCATE_MODEL1(2)
!
!-------------------------------------------------------------------------------
!
!
!*       3.     PROLOGUE:
!               --------
!
!*       3.1    Compute dimensions of model 2 and other indices
!
NIMAX_ll = NXSIZE * NDXRATIO
NJMAX_ll = NYSIZE * NDYRATIO
!
IF (NIMAX_ll==1 .AND. NJMAX_ll==1) THEN
  L1D=.TRUE.
  L2D=.FALSE.
ELSE IF (NJMAX_ll==1) THEN
  L1D=.FALSE.
  L2D=.TRUE.
ELSE
  L1D=.FALSE.
  L2D=.FALSE.
END IF
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
!
NIMAX = IIE-IIB+1
NJMAX = IJE-IJB+1
!$
IKU = SIZE(XTHVREFZ,1)
NKMAX = IKU - 2*JPVEXT           ! initialization of NKMAX (MODD_DIM2)
!
IKB = 1 + JPVEXT
IKE = IKU - JPVEXT
!
!
!*       3.2    Position of model 2 domain relative to model 1 and controls
!
!$20140506 the condition on NXSIZE*NXRATIO ==IIE-IIB+1 only works for monoproc
!$then cancel it
!IF ( (NXSIZE*NDXRATIO) /= (IIE-IIB+1) ) THEN  
!  WRITE(ILUOUT,*) 'SPAWN_MODEL2:  MODEL 2 DOMAIN X-SIZE INCOHERENT WITH THE',  &
!       ' MODEL1 MESH  ',' IIB = ',IIB,' IIE = ', IIE ,'NDXRATIO = ',NDXRATIO
! !callabortstop
!  CALL PRINT_MSG(NVERB_FATAL,'GEN','SPAWN_MODEL2','')
!END IF
!$
!$20140506 the condition on NXSIZE*NXRATIO ==IIE-IIB+1 only works for monoproc
!$then cancel it
!IF ( (NYSIZE*NDYRATIO) /= (IJE-IJB+1) ) THEN  
!  WRITE(ILUOUT,*) 'SPAWN_MODEL2:  MODEL 2 DOMAIN Y-SIZE INCOHERENT WITH THE',  &
!       ' MODEL1 MESH  ',' IJB = ',IJB,' IJE = ', IJE ,'NDYRATIO = ',NDYRATIO
! !callabortstop
!  CALL PRINT_MSG(NVERB_FATAL,'GEN','SPAWN_MODEL2','')
!END IF
!$
!
!*       3.3    Treatement of a SON 1 model (input)
!
IF (LEN_TRIM(HSONFILE) /= 0 ) THEN
!
!        3.3.1  Opening the son input file and reading the grid
! 
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: spawning with a SON input file :',TRIM(HSONFILE)
  CALL IO_File_add2list(TZSONFILE,TRIM(HSONFILE),'MNH','READ',KLFITYPE=2,KLFIVERB=NVERB)
  CALL IO_File_open(TZSONFILE)
  CALL IO_Field_read(TZSONFILE,'DAD_NAME',YDAD_SON)
  CALL IO_Field_read(TZSONFILE,'IMAX',    IIMAXSON)
  CALL IO_Field_read(TZSONFILE,'JMAX',    IJMAXSON)
  CALL IO_Field_read(TZSONFILE,'XOR',     IXORSON)
  CALL IO_Field_read(TZSONFILE,'YOR',     IYORSON)
  CALL IO_Field_read(TZSONFILE,'DXRATIO', IDXRATIOSON)
  CALL IO_Field_read(TZSONFILE,'DYRATIO', IDYRATIOSON)
  !
  IF (ADJUSTL(ADJUSTR(YDAD_SON)).NE.ADJUSTL(ADJUSTR(CMY_NAME(1)))) THEN 
    WRITE(ILUOUT,*) 'SPAWN_MODEL2: DAD of SON file is different from the one of model2'
    WRITE(ILUOUT,*) ' DAD of SON = ',TRIM(YDAD_SON),'  DAD of model2 = ',TRIM(CMY_NAME(1))
    !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SPAWN_MODEL2','DAD of SON file is different from the one of model2')
  END IF
  IF ( IDXRATIOSON /= NDXRATIO ) THEN
    WRITE(ILUOUT,*) 'SPAWN_MODEL2: RATIOX of input SON file is different from the one of model2' ,&
       ' RATIOX SON = ',IDXRATIOSON,' RATIOX model2 = ',NDXRATIO
    !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SPAWN_MODEL2','RATIOX of input SON file is different from the one of model2')
  END IF
  IF ( IDYRATIOSON /= NDYRATIO ) THEN
    WRITE(ILUOUT,*) 'SPAWN_MODEL2: RATIOY of input SON file is different from the one of model2' ,&
       ' RATIOY SON = ',IDYRATIOSON,' RATIOY model2 = ',NDYRATIO
    !callabortstop
    CALL PRINT_MSG(NVERB_FATAL,'GEN','SPAWN_MODEL2','RATIOY of input SON file is different from the one of model2')
  END IF
  !
  IIUSON=IIMAXSON+2*JPHEXT
  IJUSON=IJMAXSON+2*JPHEXT
!
!        3.3.2  Correspondance of indexes between the input SON and model2
! 
  IXSIZESON = IIMAXSON/IDXRATIOSON
  IYSIZESON = IJMAXSON/IDYRATIOSON
  IXENDSON = IXORSON+IXSIZESON
  IYENDSON = IYORSON+IYSIZESON
! Is a common domain between the input SON and the output son (model2)?
  IF( ( MIN(NXEND-1,IXENDSON)-MAX(NXOR,IXORSON) > 0 ) .OR.           &
      ( MIN(NYEND-1,IYENDSON)-MAX(NYOR,IYORSON) > 0 )                ) THEN
    GNOSON=.FALSE.
    ! Common domain for the model2 (output son) indexes
    IIB2 = (MAX(NXOR,IXORSON)-NXOR)*NDXRATIO+1+JPHEXT
    IJB2 = (MAX(NYOR,IYORSON)-NYOR)*NDYRATIO+1+JPHEXT
    IIE2 = (MIN(NXEND-1,IXENDSON)-NXOR)*NDXRATIO+JPHEXT
    IJE2 = (MIN(NYEND-1,IYENDSON)-NYOR)*NDYRATIO+JPHEXT
    ! Common domain for the SON 1 (input one) indexes
    IIB1 = (MAX(NXOR,IXORSON)-IXORSON)*NDXRATIO+1+JPHEXT
    IJB1 = (MAX(NYOR,IYORSON)-IYORSON)*NDYRATIO+1+JPHEXT
    IIE1 = (MIN(NXEND-1,IXENDSON)-IXORSON)*NDXRATIO+JPHEXT
    IJE1 = (MIN(NYEND-1,IYENDSON)-IYORSON)*NDYRATIO+JPHEXT
    ! 
    WRITE(ILUOUT,*) '   common domain in the SON grid (IB,IE=', &
                   1+JPHEXT,'-',IIMAXSON+JPHEXT,' ; JB,JE=',    &
                   1+JPHEXT,'-',IJMAXSON+JPHEXT,'):'
    WRITE(ILUOUT,*) 'I=',IIB1,'->',IIE1,' ; J=',IJB1,'->',IJE1
    WRITE(ILUOUT,*) '   common domain in the model2 grid (IB,IE=',  &
                   1+JPHEXT,'-',NXSIZE*NDXRATIO+JPHEXT,' ; JB,JE=', &
                   1+JPHEXT,'-',NYSIZE*NDYRATIO+JPHEXT,'):'
    WRITE(ILUOUT,*) 'I=',IIB2,'->',IIE2,' ; J=',IJB2,'->',IJE2
  ELSE
    WRITE(ILUOUT,*) 'SPAWN_MODEL2: no common domain between input SON and model2:'
    WRITE(ILUOUT,*) '  the input SON fields are not taken into account, spawned fields are computed from model1'
  END IF
END IF
!
!*       3.4    Initialization of model 2 configuration
! 
NRR = KRR           ! for MODD_CONF2
NSV_USER = KSV_USER
IF (NSV_CHEM>0) THEN
   LUSECHEM=.TRUE.  
   IF (NSV_CHAC>0) THEN
           LUSECHAQ=.TRUE.
   ENDIF
   IF (NSV_CHIC>0) THEN
           LUSECHIC=.TRUE.
   ENDIF 
   CCHEM_INPUT_FILE = HCHEM_INPUT_FILE
END IF
!
CTURB    =  HTURB                 ! for MODD_PARAM2
CRAD     = 'NONE'                 ! radiation will have to be restarted
CSURF    =  HSURF                 ! for surface call
CCLOUD   =  HCLOUD
CDCONV   = 'NONE'                 ! deep convection will have to be restarted
CSCONV   = 'NONE'                 ! shallow convection will have to be restarted
!
! cas LIMA 
!
!IF (HCLOUD=='LIMA') THEN
!  CCLOUD='LIMA'
!  NMOD_CCN=3
!  LSCAV=.FALSE.
!  LAERO_MASS=.FALSE.
!  NMOD_IFN=2
!  NMOD_IMM=1
!  LHHONI=.FALSE.
!ENDIF
!
CALL INI_NSV(2) ! NSV* are set equal for model 2 and model 1. 
                ! NSV is set to the total number of SV for model 2
!
IF (NRR==0) THEN
  LUSERV=.FALSE.        ! as the default is .T.
ELSE
  IDX_RVT = 1
END IF
IF (NRR>1) THEN
  LUSERC=.TRUE.
  IDX_RCT = 2
END IF
IF (NRR>2) THEN
  LUSERR=.TRUE.
  IDX_RRT = 2
END IF
IF (NRR>3) THEN
  LUSERI=.TRUE.
  IDX_RIT = 2
END IF
IF (NRR>4) THEN
  LUSERS=.TRUE.
  IDX_RST = 2
END IF
IF (NRR>5) THEN
  LUSERG=.TRUE.
  IDX_RGT = 2
END IF
IF (NRR>6) THEN
  LUSERH=.TRUE.
  IDX_RHT = 2
END IF
!
!
!
!*       3.5   model 2 configuration in MODD_NESTING to be written
!*                on the FM-file to allow nesting or coupling 
!
CCPLFILE(:) = '    ' 
LSTEADYLS=.TRUE.
!
NDXRATIO_ALL(:) = 0
NDYRATIO_ALL(:) = 0
NDXRATIO_ALL(2) = NDXRATIO
NDYRATIO_ALL(2) = NDYRATIO
NXOR_ALL(2)     = NXOR
NYOR_ALL(2)     = NYOR
NXEND_ALL(2)    = NXEND
NYEND_ALL(2)    = NYEND
!
!*       3.6   size of the RIM area for lbc 
!
NRIMX=MIN(JPRIMMAX,IIU/2-1)
IF ( .NOT. L2D ) THEN
  NRIMY=MIN(JPRIMMAX,IJU/2-1)
ELSE
  NRIMY=0
END IF
IF (NRIMX >= IIU/2-1) THEN      ! Error ! this case is not supported - it should be, but there is a bug
  call Print_msg( NVERB_FATAL, 'GEN', 'SPAWN_MODEL2', 'The size of the LBX zone is too big for the size of the subdomains. '// &
                  'Try with less processes, a smaller LBX size or a bigger grid in X.' )
ENDIF
IF ( ( .NOT. L2D ) .AND. (NRIMY >= IJU/2-1) ) THEN  ! Error ! this case is not supported - it should be, but there is a bug
  call Print_msg( NVERB_FATAL, 'GEN', 'SPAWN_MODEL2', 'The size of the LBY zone is too big for the size of the subdomains. '// &
                  'Try with less processes, a smaller LBY size or a bigger grid in Y.' )
ENDIF
!
LHORELAX_UVWTH=.TRUE.
LHORELAX_RV=LUSERV
LHORELAX_RC=LUSERC
LHORELAX_RR=LUSERR
LHORELAX_RI=LUSERI
LHORELAX_RS=LUSERS
LHORELAX_RG=LUSERG
LHORELAX_RH=LUSERH
!
IF (CTURB/='NONE') LHORELAX_TKE  =.TRUE.
LHORELAX_SV(:)=.FALSE.
DO JSV=1,NSV
  LHORELAX_SV(JSV)=.TRUE.
END DO
IF (NSV_CHEM > 0) LHORELAX_SVCHEM = .TRUE.
IF (NSV_CHIC > 0) LHORELAX_SVCHIC = .TRUE.
IF (NSV_C2R2 > 0) LHORELAX_SVC2R2 = .TRUE.
IF (NSV_C1R3 > 0) LHORELAX_SVC1R3 = .TRUE.
IF (NSV_ELEC > 0) LHORELAX_SVELEC = .TRUE.
IF (NSV_AER  > 0) LHORELAX_SVAER = .TRUE.
IF (NSV_DST  > 0) LHORELAX_SVDST = .TRUE.
IF (NSV_SLT  > 0) LHORELAX_SVSLT = .TRUE.
IF (NSV_PP  > 0) LHORELAX_SVPP   = .TRUE.
#ifdef MNH_FOREFIRE
IF (NSV_FF  > 0) LHORELAX_SVFF   = .TRUE.
#endif
IF (NSV_CS  > 0) LHORELAX_SVCS   = .TRUE.
LHORELAX_SVLG   = .FALSE.
IF (NSV_LIMA > 0) LHORELAX_SVLIMA = .TRUE.
!
!-------------------------------------------------------------------------------
!
!*       4.    ALLOCATE MEMORY FOR ARRAYS :  
!	       -----------------------------
!
!*       4.1  Global variables absent from the modules :
!                  
ALLOCATE(ZJ(IIU,IJU,IKU))                      
!
!*       4.2   Prognostic (and diagnostic) variables (module MODD_FIELD2) :
!
ALLOCATE(XZWS(IIU,IJU)); XZWS(:,:) = XZWS_DEFAULT
ALLOCATE(XLSZWSM(IIU,IJU))
ALLOCATE(XUT(IIU,IJU,IKU))
ALLOCATE(XVT(IIU,IJU,IKU))
ALLOCATE(XWT(IIU,IJU,IKU))
ALLOCATE(XTHT(IIU,IJU,IKU))
IF (CTURB/='NONE') THEN
  ALLOCATE(XTKET(IIU,IJU,IKU))
ELSE
  ALLOCATE(XTKET(0,0,0))
END IF
ALLOCATE(XPABST(IIU,IJU,IKU))
ALLOCATE(XRT(IIU,IJU,IKU,NRR))
ALLOCATE(XSVT(IIU,IJU,IKU,NSV))
!
IF (CTURB /= 'NONE' .AND. NRR>1) THEN
  ALLOCATE(XSRCT(IIU,IJU,IKU))
  ALLOCATE(XSIGS(IIU,IJU,IKU))
ELSE
  ALLOCATE(XSRCT(0,0,0))
  ALLOCATE(XSIGS(0,0,0))
END IF
!
!
!*       4.4   Grid variables (module MODD_GRID2 and MODD_METRICS2):
!
ALLOCATE(XXHAT(IIU),XYHAT(IJU),XZHAT(IKU))
ALLOCATE(XZTOP)
ALLOCATE(XMAP(IIU,IJU))
ALLOCATE(XLAT(IIU,IJU))
ALLOCATE(XLON(IIU,IJU))
ALLOCATE(XDXHAT(IIU),XDYHAT(IJU))
ALLOCATE(XZS(IIU,IJU))
ALLOCATE(XZSMT(IIU,IJU))
ALLOCATE(XZZ(IIU,IJU,IKU))
!
ALLOCATE(XDXX(IIU,IJU,IKU))
ALLOCATE(XDYY(IIU,IJU,IKU))
ALLOCATE(XDZX(IIU,IJU,IKU))
ALLOCATE(XDZY(IIU,IJU,IKU))
ALLOCATE(XDZZ(IIU,IJU,IKU))
!
ALLOCATE(ZZS_LS(IIU,IJU))
ALLOCATE(ZZSMT_LS(IIU,IJU))
ALLOCATE(ZZZ_LS(IIU,IJU,IKU))
!
!*       4.5   Reference state variables (module MODD_REF2):
!
ALLOCATE(XRHODREF(IIU,IJU,IKU),XTHVREF(IIU,IJU,IKU),XRVREF(IIU,IJU,IKU))
ALLOCATE(XRHODJ(IIU,IJU,IKU),XEXNREF(IIU,IJU,IKU))
!
!*       4.6   Larger Scale fields (module MODD_LSFIELD2):
!
                !          LS fields for vertical relaxation and diffusion
ALLOCATE(XLSUM(IIU,IJU,IKU))
ALLOCATE(XLSVM(IIU,IJU,IKU))
ALLOCATE(XLSWM(IIU,IJU,IKU))
ALLOCATE(XLSTHM(IIU,IJU,IKU))
IF ( NRR >= 1) THEN
  ALLOCATE(XLSRVM(IIU,IJU,IKU))
ELSE
  ALLOCATE(XLSRVM(0,0,0))
ENDIF
                !          LB fields for lbc coupling
!
!get the size of the local portion of the LB zone in X and Y direction
CALL GET_SIZEX_LB(NIMAX_ll,NJMAX_ll,NRIMX,               &
                  IISIZEXF,IJSIZEXF,IISIZEXFU,IJSIZEXFU, &
                  IISIZEX4,IJSIZEX4,IISIZEX2,IJSIZEX2)
CALL GET_SIZEY_LB(NIMAX_ll,NJMAX_ll,NRIMY,               &
                  IISIZEYF,IJSIZEYF,IISIZEYFV,IJSIZEYFV, &
                  IISIZEY4,IJSIZEY4,IISIZEY2,IJSIZEY2)
!on fait des choses inutiles avec GET_SIZEX_LB, on pourrait utiliser seulement GET_LOCAL_LB_SIZE_X_ll
!ILOCLBSIZEX = GET_LOCAL_LB_SIZE_X_ll( NRIMX )
!ILOCLBSIZEY = GET_LOCAL_LB_SIZE_Y_ll( NRIMY )
!
  ALLOCATE(XLBXUM(IISIZEXFU,IJU,IKU))
!! ALLOCATE(XLBXUM(2*NRIMX+2*JPHEXT,IJU,IKU))
!
IF ( .NOT. L2D ) THEN
  ALLOCATE(XLBYUM(IIU,IJSIZEYF,IKU))
!!  ALLOCATE(XLBYUM(IIU,2*NRIMY+2*JPHEXT,IKU))
ELSE
  ALLOCATE(XLBYUM(0,0,0))
END IF
!
ALLOCATE(XLBXVM(IISIZEXF,IJU,IKU))
!! ALLOCATE(XLBXVM(2*NRIMX+2*JPHEXT,IJU,IKU))
!
IF ( .NOT. L2D ) THEN
  IF ( NRIMY == 0 ) THEN
    ALLOCATE(XLBYVM(IIU,IJSIZEY4,IKU))
  ELSE
    ALLOCATE(XLBYVM(IIU,IJSIZEYFV,IKU))
!!    ALLOCATE(XLBYVM(IIU,2*NRIMY+2*JPHEXT,IKU))
  END IF
ELSE
  ALLOCATE(XLBYVM(0,0,0))
END IF
!
ALLOCATE(XLBXWM(IISIZEXF,IJU,IKU))
!! ALLOCATE(XLBXWM(2*NRIMX+2*JPHEXT,IJU,IKU))
!
IF ( .NOT. L2D ) THEN
  ALLOCATE(XLBYWM(IIU,IJSIZEYF,IKU))
!!  ALLOCATE(XLBYWM(IIU,2*NRIMY+2*JPHEXT,IKU))
ELSE
  ALLOCATE(XLBYWM(0,0,0))
END IF
!
ALLOCATE(XLBXTHM(IISIZEXF,IJU,IKU))
!!ALLOCATE(XLBXTHM(2*NRIMX+2*JPHEXT,IJU,IKU))
!
IF ( .NOT. L2D )  THEN
  ALLOCATE(XLBYTHM(IIU,IJSIZEYF,IKU))
!!  ALLOCATE(XLBYTHM(IIU,2*NRIMY+2*JPHEXT,IKU))
ELSE
  ALLOCATE(XLBYTHM(0,0,0))
END IF
!
IF (CTURB /= 'NONE') THEN
  ALLOCATE(XLBXTKEM(IISIZEXF,IJU,IKU))
!!  ALLOCATE(XLBXTKEM(2*NRIMX+2*JPHEXT,IJU,IKU))
ELSE
  ALLOCATE(XLBXTKEM(0,0,0))
END IF
!
IF (CTURB /= 'NONE' .AND. (.NOT. L2D)) THEN
  ALLOCATE(XLBYTKEM(IIU,IJSIZEYF,IKU))
!!  ALLOCATE(XLBYTKEM(IIU,2*NRIMY+2*JPHEXT,IKU))
ELSE
  ALLOCATE(XLBYTKEM(0,0,0))
END IF
!
ALLOCATE(XLBXRM(IISIZEXF,IJU,IKU,NRR))
!!ALLOCATE(XLBXRM(2*NRIMX+2*JPHEXT,IJU,IKU,NRR))
!
IF (.NOT. L2D ) THEN
  ALLOCATE(XLBYRM(IIU,IJSIZEYF,IKU,NRR))
!!  ALLOCATE(XLBYRM(IIU,2*NRIMY+2*JPHEXT,IKU,NRR))
ELSE
  ALLOCATE(XLBYRM(0,0,0,0))
END IF
!
ALLOCATE(XLBXSVM(IISIZEXF,IJU,IKU,NSV))
!!ALLOCATE(XLBXSVM(2*NRIMX+2*JPHEXT,IJU,IKU,NSV))
!
IF (.NOT. L2D ) THEN
  ALLOCATE(XLBYSVM(IIU,IJSIZEYF,IKU,NSV))
!!  ALLOCATE(XLBYSVM(IIU,2*NRIMY+2*JPHEXT,IKU,NSV))
ELSE
  ALLOCATE(XLBYSVM(0,0,0,0))
END IF
!
NSIZELBX_ll=2*NRIMX+2*JPHEXT
NSIZELBXU_ll=2*NRIMX+2*JPHEXT
NSIZELBY_ll=2*NRIMY+2*JPHEXT
NSIZELBYV_ll=2*NRIMY+2*JPHEXT
NSIZELBXR_ll=2*NRIMX+2*JPHEXT
NSIZELBXSV_ll=2*NRIMX+2*JPHEXT
NSIZELBXTKE_ll=2*NRIMX+2*JPHEXT
NSIZELBYTKE_ll=2*NRIMY+2*JPHEXT
NSIZELBYR_ll=2*NRIMY+2*JPHEXT
NSIZELBYSV_ll=2*NRIMY+2*JPHEXT
!
!
!        4.8   precipitation variables  ! same allocations than in ini_micron
!
IF (CCLOUD /= 'NONE' .AND. CCLOUD /= 'REVE') THEN
  ALLOCATE(XINPRR(IIU,IJU))
  ALLOCATE(XINPRR3D(IIU,IJU,IKU))
  ALLOCATE(XEVAP3D(IIU,IJU,IKU))
  ALLOCATE(XACPRR(IIU,IJU))
ELSE
  ALLOCATE(XINPRR(0,0))
  ALLOCATE(XINPRR3D(0,0,0))
  ALLOCATE(XEVAP3D(0,0,0))
  ALLOCATE(XACPRR(0,0))
END IF
!
IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C2R2'  &
         .OR. CCLOUD == 'KHKO' .OR. CCLOUD == 'LIMA') THEN
  ALLOCATE(XINPRC(IIU,IJU))
  ALLOCATE(XACPRC(IIU,IJU))
ELSE
  ALLOCATE(XINPRC(0,0))
  ALLOCATE(XACPRC(0,0))
END IF
!
IF (( CCLOUD(1:3) == 'ICE'                                   .AND.LDEPOSC) .OR. &
    ((CCLOUD=='C2R2' .OR. CCLOUD=='KHKO').AND.LDEPOC)  .OR. &
    ( CCLOUD=='LIMA'                                         .AND.MDEPOC))  THEN
  ALLOCATE(XINDEP(IIU,IJU))
  ALLOCATE(XACDEP(IIU,IJU))
ELSE
  ALLOCATE(XINDEP(0,0))
  ALLOCATE(XACDEP(0,0))
END IF
!
IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5'.OR. CCLOUD == 'LIMA') THEN
  ALLOCATE(XINPRS(IIU,IJU))
  ALLOCATE(XACPRS(IIU,IJU))
ELSE
  ALLOCATE(XINPRS(0,0))
  ALLOCATE(XACPRS(0,0))
END IF
!
IF (CCLOUD == 'C3R5' .OR. CCLOUD == 'ICE3' .OR. CCLOUD == 'ICE4'.OR. CCLOUD == 'LIMA' ) THEN
  ALLOCATE(XINPRG(IIU,IJU))
  ALLOCATE(XACPRG(IIU,IJU))
ELSE
  ALLOCATE(XINPRG(0,0))
  ALLOCATE(XACPRG(0,0))
END IF
!
IF (CCLOUD == 'ICE4'.OR. CCLOUD == 'LIMA') THEN
  ALLOCATE(XINPRH(IIU,IJU))
  ALLOCATE(XACPRH(IIU,IJU))
ELSE
  ALLOCATE(XINPRH(0,0))
  ALLOCATE(XACPRH(0,0))
END IF
!
IF ( CCLOUD=='LIMA' .AND. LSCAV ) THEN
  ALLOCATE(XINPAP(IIU,IJU))
  ALLOCATE(XACPAP(IIU,IJU))
  XINPAP(:,:)=0.0
  XACPAP(:,:)=0.0  
ELSE
  ALLOCATE(XINPAP(0,0))
  ALLOCATE(XACPAP(0,0))
END IF
!
!        4.8bis electric variables  
!
IF (CELEC /= 'NONE' ) THEN
  ALLOCATE(XNI_SDRYG(IIU,IJU,IKU))
  ALLOCATE(XNI_IDRYG(IIU,IJU,IKU))
  ALLOCATE(XNI_IAGGS(IIU,IJU,IKU))
  ALLOCATE(XEFIELDU(IIU,IJU,IKU))
  ALLOCATE(XEFIELDV(IIU,IJU,IKU))
  ALLOCATE(XEFIELDW(IIU,IJU,IKU))
  ALLOCATE(XESOURCEFW(IIU,IJU,IKU))
  ALLOCATE(XIND_RATE(IIU,IJU,IKU))
  ALLOCATE(XIONSOURCEFW(IIU,IJU,IKU))
  ALLOCATE(XEW(IIU,IJU,IKU))
  ALLOCATE(XCION_POS_FW(IIU,IJU,IKU))
  ALLOCATE(XCION_NEG_FW(IIU,IJU,IKU))
  ALLOCATE(XMOBIL_POS(IIU,IJU,IKU))
  ALLOCATE(XMOBIL_NEG(IIU,IJU,IKU))
ELSE
  ALLOCATE(XNI_SDRYG(0,0,0))
  ALLOCATE(XNI_IDRYG(0,0,0))
  ALLOCATE(XNI_IAGGS(0,0,0))
  ALLOCATE(XEFIELDU(0,0,0))
  ALLOCATE(XEFIELDV(0,0,0))
  ALLOCATE(XEFIELDW(0,0,0))
  ALLOCATE(XESOURCEFW(0,0,0))
  ALLOCATE(XIND_RATE(0,0,0))
  ALLOCATE(XIONSOURCEFW(0,0,0))
  ALLOCATE(XEW(0,0,0))
  ALLOCATE(XCION_POS_FW(0,0,0))
  ALLOCATE(XCION_NEG_FW(0,0,0))
  ALLOCATE(XMOBIL_POS(0,0,0))
  ALLOCATE(XMOBIL_NEG(0,0,0))
END IF
!
!
!
!        4.9   Passive pollutant variable                                    
!
IF (LPASPOL) THEN
  ALLOCATE( XATC(IIU,IJU,IKU,NSV_PP) )
             ELSE
  ALLOCATE( XATC(0,0,0,0))
END IF
!
!        4.10  Advective forcing variable for 2D (Modif MT)
!
!
IF (L2D_ADV_FRC) THEN
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: L2D_ADV_FRC IS SET TO ',L2D_ADV_FRC,' SO ADVECTIVE FORCING WILL BE SPAWN: NADVFRC=',NADVFRC
  ALLOCATE(TDTADVFRC(NADVFRC))
  ALLOCATE(XDTHFRC(IIU,IJU,IKU,NADVFRC))
  ALLOCATE(XDRVFRC(IIU,IJU,IKU,NADVFRC))
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: ALLOCATION OF ADV FORCING VARIABLES MADE'
ELSE
  ALLOCATE(TDTADVFRC(0))
  ALLOCATE(XDTHFRC(0,0,0,0))
  ALLOCATE(XDRVFRC(0,0,0,0))
END IF
IF (L2D_REL_FRC) THEN
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: L2D_REL_FRC IS SET TO ',L2D_REL_FRC,' SO RELAXATION FORCING WILL BE SPAWN: NRELFRC=',NRELFRC
  ALLOCATE(TDTRELFRC(NRELFRC))
  ALLOCATE(XTHREL(IIU,IJU,IKU,NRELFRC))
  ALLOCATE(XRVREL(IIU,IJU,IKU,NRELFRC))
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: ALLOCATION OF REL FORCING VARIABLES MADE'
ELSE
  ALLOCATE(TDTRELFRC(0))
  ALLOCATE(XTHREL(0,0,0,0))
  ALLOCATE(XRVREL(0,0,0,0))
END IF
!
!        4.11  Turbulent fluxes for 2D (Modif MT)                                    
!
!
IF (LUV_FLX) THEN
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: XUV_FLX1 IS SET TO ',XUV_FLX1,' SO XVU_FLUX WILL BE SPAWN'
  ALLOCATE(XVU_FLUX_M(IIU,IJU,IKU))
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: ALLOCATION OF XVU_FLUX_M  MADE'
ELSE
  ALLOCATE(XVU_FLUX_M(0,0,0))
END IF
!
IF (LTH_FLX) THEN
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: XTH_FLX IS SET TO ',XTH_FLX,' SO XVTH_FLUX and XWTH_FLUX WILL BE SPAWN'
  ALLOCATE(XVTH_FLUX_M(IIU,IJU,IKU))
  ALLOCATE(XWTH_FLUX_M(IIU,IJU,IKU))
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: ALLOCATION OF XVTH_FLUX_M and XWTH_FLUX_M  MADE'
ELSE
  ALLOCATE(XVTH_FLUX_M(0,0,0))
  ALLOCATE(XWTH_FLUX_M(0,0,0))
END IF
!
!-------------------------------------------------------------------------------
!
!*       5.     INITIALIZE ALL THE MODEL VARIABLES
!	        ----------------------------------
!
!*       5.1    Bikhardt interpolation coefficients computation :
!
CALL INI_BIKHARDT_n(NDXRATIO,NDYRATIO,2)
!
CALL SECOND_MNH(ZTIME2)
!
ZMISC = ZTIME2 - ZTIME1
!
!*       5.2    Spatial and Temporal grid (for MODD_GRID2 and MODD_TIME2) :
!
CALL SECOND_MNH(ZTIME1)
!
IF(NPROC.GT.1)THEN
        CALL GO_TOMODEL_ll(2, IINFO_ll)
        CALL GET_FEEDBACK_COORD_ll(NXOR_TMP,NYOR_TMP,NXEND_TMP,NYEND_TMP,IINFO_ll) !phys domain
ELSE
        NXOR_TMP = NXOR
        NYOR_TMP = NYOR
        NXEND_TMP= NXEND
        NYEND_TMP = NYEND
ENDIF
XZS=0.
CALL SPAWN_GRID2 (NXOR,NYOR,NXEND,NYEND,NDXRATIO,NDYRATIO,                    &
                  XLONORI,XLATORI,XXHAT,XYHAT,XZHAT,XZTOP,LSLEVE,XLEN1,XLEN2, &
                  XZS,XZSMT,ZZS_LS,ZZSMT_LS,TDTMOD,TDTCUR                     )
!
CALL MPPDB_CHECK2D(ZZS_LS,"SPAWN_MOD2:ZZS_LS",PRECISION)
CALL MPPDB_CHECK2D(ZZSMT_LS,"SPAWN_MOD2:ZZSMT_LS",PRECISION)
CALL MPPDB_CHECK2D(XZS,"SPAWN_MOD2:XZS",PRECISION)
CALL MPPDB_CHECK2D(XZSMT,"SPAWN_MOD2:XZSMT",PRECISION)
!
CALL SECOND_MNH(ZTIME2)
!
ZGRID2 = ZTIME2 - ZTIME1
!
!*       5.3    Calculation of the grid
!
ZTIME1 = ZTIME2
!
IF (LCARTESIAN) THEN
  CALL SM_GRIDCART(XXHAT,XYHAT,XZHAT,ZZS_LS,LSLEVE,XLEN1,XLEN2,ZZSMT_LS,XDXHAT,XDYHAT,ZZZ_LS,ZJ)
  CALL SM_GRIDCART(XXHAT,XYHAT,XZHAT,XZS   ,LSLEVE,XLEN1,XLEN2,XZSMT   ,XDXHAT,XDYHAT,XZZ   ,ZJ)
ELSE
  CALL SM_GRIDPROJ(XXHAT,XYHAT,XZHAT,ZZS_LS,LSLEVE,XLEN1,XLEN2,ZZSMT_LS,&
                   XLATORI,XLONORI,XMAP,XLAT,XLON,XDXHAT,XDYHAT,ZZZ_LS,ZJ)
  CALL SM_GRIDPROJ(XXHAT,XYHAT,XZHAT,XZS   ,LSLEVE,XLEN1,XLEN2,XZSMT   ,&
                   XLATORI,XLONORI,XMAP,XLAT,XLON,XDXHAT,XDYHAT,XZZ   ,ZJ)
END IF
!
!*       5.4  Compute the metric coefficients
!
CALL ADD3DFIELD_ll( TZFIELDS_ll, XZZ, 'SPAWN_MODEL2::XZZ' )
CALL UPDATE_HALO_ll(TZFIELDS_ll,IINFO_ll)
CALL CLEANLIST_ll(TZFIELDS_ll)
!
CALL METRICS(XMAP,XDXHAT,XDYHAT,XZZ,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
CALL MPPDB_CHECK3D(XDXX,"spawnmod2-beforeupdate_metrics:XDXX",PRECISION)
CALL MPPDB_CHECK3D(XDYY,"spawnmod2-beforeupdate_metrics:XDYY",PRECISION)
CALL MPPDB_CHECK3D(XDZX,"spawnmod2-beforeupdate_metrics:XDZX",PRECISION)
CALL MPPDB_CHECK3D(XDZY,"spawnmod2-beforeupdate_metrics:XDZY",PRECISION)
!
CALL UPDATE_METRICS(CLBCX,CLBCY,XDXX,XDYY,XDZX,XDZY,XDZZ)
!
CALL MPPDB_CHECK3D(XDXX,"spawnmod2-aftrupdate_metrics:XDXX",PRECISION)
CALL MPPDB_CHECK3D(XDYY,"spawnmod2-aftrupdate_metrics:XDYY",PRECISION)
CALL MPPDB_CHECK3D(XDZX,"spawnmod2-aftrupdate_metrics:XDZX",PRECISION)
CALL MPPDB_CHECK3D(XDZY,"spawnmod2-aftrupdate_metrics:XDZY",PRECISION)
!$
!
!*       5.5    3D Reference state variables :
!
CALL SET_REF(0,TFILE_DUMMY,                        &
             XZZ,XZHAT,ZJ,XDXX,XDYY,CLBCX,CLBCY,   &
             XREFMASS,XMASS_O_PHI0,XLINMASS,       &
             XRHODREF,XTHVREF,XRVREF,XEXNREF,XRHODJ)
!
CALL SECOND_MNH(ZTIME2)
!
ZMISC = ZMISC + ZTIME2 - ZTIME1
!
!*       5.6    Prognostic variables and Larger scale fields :
!
ZTIME1 = ZTIME2
!
!* horizontal interpolation
!
ALLOCATE(ZTHVT(IIU,IJU,IKU))
ALLOCATE(ZHUT(IIU,IJU,IKU))
!
MPPDB_CHECK_LB = .TRUE.
IF (GNOSON) THEN
  CALL SPAWN_FIELD2 (NXOR,NYOR,NXEND,NYEND,NDXRATIO,NDYRATIO,CTURB,            &
                 XUT,XVT,XWT,ZTHVT,XRT,ZHUT,XTKET,XSVT,XZWS,XATC,              &
                 XSRCT,XSIGS,                                                  &
                 XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XLSZWSM,                      &
                 XDTHFRC,XDRVFRC,XTHREL,XRVREL,                                &
                 XVU_FLUX_M,XVTH_FLUX_M,XWTH_FLUX_M            )
  CALL MPPDB_CHECK3D(XUT,"SPAWN_M2 after SPAWN_FIELD2:XUT",PRECISION)
ELSE
  CALL MPPDB_CHECK3D(XUT,"SPAWN_M2 before SPAWN_FIELD2:XUT",PRECISION)
  CALL SPAWN_FIELD2 (NXOR,NYOR,NXEND,NYEND,NDXRATIO,NDYRATIO,CTURB,            &
                 XUT,XVT,XWT,ZTHVT,XRT,ZHUT,XTKET,XSVT,XZWS,XATC,              &
                 XSRCT,XSIGS,                                                  &
                 XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM,XLSZWSM,                      &
                 XDTHFRC,XDRVFRC,XTHREL,XRVREL,                                &                 
                 XVU_FLUX_M, XVTH_FLUX_M,XWTH_FLUX_M,                          &
                 TZSONFILE,IIUSON,IJUSON,                                      &
                 IIB2,IJB2,IIE2,IJE2,                                          &
                 IIB1,IJB1,IIE1,IJE1                                           )
  CALL MPPDB_CHECK3D(XUT,"SPAWN_M2 after SPAWN_FIELD2:XUT",PRECISION)
END IF
!
CALL MPPDB_CHECK3D(XUT,"SPAWN_MOD2aftFIELD2:XUT",PRECISION)
CALL MPPDB_CHECK3D(XVT,"SPAWN_MOD2aftFIELD2:XVT",PRECISION)
!$
!* correction of positivity
!
IF (SIZE(XLSRVM,1)>0)      XLSRVM   = MAX(0.,XLSRVM)
IF (SIZE(XRT,1)>0)         XRT      = MAX(0.,XRT)
IF (SIZE(ZHUT,1)>0)        ZHUT     = MIN(MAX(ZHUT,0.),100.)
IF (SIZE(XTKET,1)>0)       XTKET    = MAX(XTKEMIN,XTKET)
!
CALL SECOND_MNH(ZTIME2)
!
ZFIELD2 = ZTIME2 - ZTIME1
!
ZTIME1  = ZTIME2
!
!* vertical interpolation
!
ZZS_MAX = ABS( MAXVAL(XZS(:,:)))
CALL MPI_ALLREDUCE(ZZS_MAX, ZZS_MAX_ll, 1, MNHREAL_MPI, MPI_MAX,  &
                     NMNH_COMM_WORLD,IINFO_ll)
IF ( (ZZS_MAX_ll>0.) .AND. (NDXRATIO/=1 .OR. NDYRATIO/=1) )  THEN
  CALL MPPDB_CHECK3D(XUT,"SPAWN_M2 before VER_INTERP_FIELD:XUT",PRECISION)
  CALL VER_INTERP_FIELD (CTURB,NRR,NSV,ZZZ_LS,XZZ,                             &
               XUT,XVT,XWT,ZTHVT,XRT,ZHUT,XTKET,XSVT,                          &
               XSRCT,XSIGS,                                                    &
               XLSUM,XLSVM,XLSWM,XLSTHM,XLSRVM                                 )
  !
  CALL MPPDB_CHECK3D(XUT,"SPAWN_M2aftVERINTER:XUT",PRECISION)
  CALL MPPDB_CHECK3D(XVT,"SPAWN_M2aftVERINTER:XVT",PRECISION)
  CALL MPPDB_CHECK3D(XWT,"SPAWN_M2aftVERINTER:XWT",PRECISION)
  CALL MPPDB_CHECK3D(ZHUT,"SPAWN_M2aftVERINTER:ZHUT",PRECISION)
  CALL MPPDB_CHECK3D(XTKET,"SPAWN_M2aftVERINTER:XTKET",PRECISION)
  CALL MPPDB_CHECK3D(XSRCT,"SPAWN_M2aftVERINTER:XSRCT",PRECISION)
ENDIF
!
CALL SECOND_MNH(ZTIME2)
!
ZVER = ZTIME2 - ZTIME1
!
!*       5.7    Absolute pressure :
!
ZTIME1 = ZTIME2
!
CALL SPAWN_PRESSURE2(NXOR,NYOR,NXEND,NYEND,NDXRATIO,NDYRATIO,   &
                     ZZZ_LS,XZZ,ZTHVT,XPABST                    )
!
IF (.NOT.GNOSON) THEN
  ALLOCATE(ZWORK3D(IIUSON,IJUSON,IKU))
  CALL IO_Field_read(TZSONFILE,'PABST',ZWORK3D)
  XPABST(IIB2:IIE2,IJB2:IJE2,:) = ZWORK3D(IIB1:IIE1,IJB1:IJE1,:)
  DEALLOCATE(ZWORK3D)
END IF
!
IF (NVERB>=2) THEN
  IK4000 = COUNT(XZHAT(:)<4000.)
  IIJ = MAXLOC(        SUM(ZHUT(IIB:IIE,IJB:IJE,JPVEXT+1:IK4000),3),                  &
                MASK=COUNT(ZHUT(IIB:IIE,IJB:IJE,JPVEXT+1:IKE)                         &
                           >=MAXVAL(ZHUT(IIB:IIE,IJB:IJE,JPVEXT+1:IKE))-0.01,DIM=3 )  &
                      >=1                                                   )           &
        + JPHEXT
  WRITE(ILUOUT,*) ' '
  WRITE(ILUOUT,*) 'humidity     (I=',IIJ(1),';J=',IIJ(2),')'
  DO JK=IKB,IKE
    WRITE(ILUOUT,'(F6.2," %")') ZHUT(IIJ(1),IIJ(2),JK)
  END DO
END IF
!*       5.8    Retrieve model thermodynamical variables :
!
ALLOCATE(ZSUMRT(IIU,IJU,IKU))
ZSUMRT(:,:,:) = 0.
IF (NRR==0) THEN
  XTHT(:,:,:) = ZTHVT(:,:,:)
ELSE
  IF (NDXRATIO/=1 .OR. NDYRATIO/=1) THEN
    XRT(:,:,:,1) = SM_PMR_HU(XPABST(:,:,:),                                 &
                             ZTHVT(:,:,:)*(XPABST(:,:,:)/XP00)**(XRD/XCPD), &
                             ZHUT(:,:,:),XRT(:,:,:,:),KITERMAX=100          )
  END IF
  !
  DO JRR=1,NRR
    ZSUMRT(:,:,:) = ZSUMRT(:,:,:) + XRT(:,:,:,JRR)
  END DO
  XTHT(:,:,:) = ZTHVT(:,:,:)/(1.+XRV/XRD*XRT(:,:,:,1))*(1.+ZSUMRT(:,:,:))
  CALL MPPDB_CHECK3D(XTHT,"SPAWN_MOD2:XTHT",PRECISION)
END IF
!
DEALLOCATE (ZHUT)
!
CALL SECOND_MNH(ZTIME2)
ZPRESSURE2=ZTIME2-ZTIME1
!
!*       5.9   Large Scale field for lbc treatment:
!
!
!*       5.9.1 West-East LB zones
!
!
!JUAN A REVOIR TODO_JPHEXT
! <<<<<<< spawn_model2.f90
  MPPDB_CHECK_LB = .TRUE.
  CALL MPPDB_CHECK3D(XUT,"SPAWN_MOD2 before lbc treatment:XUT",PRECISION)
  CALL MPPDB_CHECK3D(XVT,"SPAWN_MOD2 before lbc treatment:XVT",PRECISION)
  MPPDB_CHECK_LB = .FALSE.
  YLBTYPE = 'LBU'
  CALL SET_LB_FIELD_ll( YLBTYPE, XUT, XLBXUM, XLBYUM, IIB, IJB, IIE, IJE, 1, 0, 0, 0 )
  ! copy XUT(IIB:IIB+NRIMX,:,:) instead of XUT(IIB-1:IIB-1+NRIMX,:,:) in XLBXUM
  CALL SET_LB_FIELD_ll( YLBTYPE, XVT, XLBXVM, XLBYVM, IIB, IJB, IIE, IJE, 0, 0, 1, 0 )
  ! copy XVT(:,IJB:IJB+NRIMY,:) instead of XVT(:,IJB-1:IJB-1+NRIMY,:) in XLBYVM
  CALL SET_LB_FIELD_ll( YLBTYPE, XWT, XLBXWM, XLBYWM, IIB, IJB, IIE, IJE, 0, 0, 0, 0 )
  CALL SET_LB_FIELD_ll( YLBTYPE, XTHT, XLBXTHM, XLBYTHM, IIB, IJB, IIE, IJE, 0, 0, 0, 0 )
  IF (HTURB /= 'NONE') THEN
    CALL SET_LB_FIELD_ll( YLBTYPE, XTKET, XLBXTKEM, XLBYTKEM, IIB, IJB, IIE, IJE, 0, 0, 0, 0 )
  ENDIF
  IF (NRR >= 1) THEN
    DO JRR =1,NRR
      CALL SET_LB_FIELD_ll( YLBTYPE, XRT(:,:,:,JRR), XLBXRM(:,:,:,JRR), XLBYRM(:,:,:,JRR), IIB, IJB, IIE, IJE, 0, 0, 0, 0 )
    END DO
  END IF
  IF (NSV /= 0) THEN
    DO JSV = 1, NSV
      CALL SET_LB_FIELD_ll( YLBTYPE, XSVT(:,:,:,JSV), XLBXSVM(:,:,:,JSV), XLBYSVM(:,:,:,JSV), IIB, IJB, IIE, IJE, 0, 0, 0, 0 )
    END DO
!!$=======
!!$!
!!$XLBXUM(1:NRIMX+JPHEXT,:,:)         = XUT(2:NRIMX+JPHEXT+1,:,:)
!!$XLBXUM(NRIMX+JPHEXT+1:2*NRIMX+2*JPHEXT,:,:) = XUT(IIE+1-NRIMX:IIE+JPHEXT,:,:)
!!$IF( .NOT. L2D ) THEN
!!$  XLBYUM(:,1:NRIMY+JPHEXT,:)         = XUT(:,1:NRIMY+JPHEXT,:)
!!$  XLBYUM(:,NRIMY+JPHEXT+1:2*NRIMY+2*JPHEXT,:) = XUT(:,IJE+1-NRIMY:IJE+JPHEXT,:)
!!$END IF
!!$!
!!$!*       5.9.2  V variable
!!$!
!!$!
!!$XLBXVM(1:NRIMX+JPHEXT,:,:)         = XVT(1:NRIMX+JPHEXT,:,:)
!!$XLBXVM(NRIMX+JPHEXT+1:2*NRIMX+2*JPHEXT,:,:) = XVT(IIE+1-NRIMX:IIE+JPHEXT,:,:)
!!$IF( .NOT. L2D ) THEN
!!$  XLBYVM(:,1:NRIMY+JPHEXT,:)         = XVT(:,2:NRIMY+JPHEXT+1,:)
!!$  XLBYVM(:,NRIMY+JPHEXT+1:2*NRIMY+2*JPHEXT,:) = XVT(:,IJE+1-NRIMY:IJE+JPHEXT,:)
!!$END IF
!!$!
!!$!*       5.9.3  W variable
!!$!
!!$!
!!$XLBXWM(1:NRIMX+JPHEXT,:,:)         = XWT(1:NRIMX+JPHEXT,:,:)
!!$XLBXWM(NRIMX+JPHEXT+1:2*NRIMX+2*JPHEXT,:,:) = XWT(IIE+1-NRIMX:IIE+JPHEXT,:,:)
!!$IF( .NOT. L2D ) THEN
!!$  XLBYWM(:,1:NRIMY+JPHEXT,:)         = XWT(:,1:NRIMY+JPHEXT,:)
!!$  XLBYWM(:,NRIMY+JPHEXT+1:2*NRIMY+2*JPHEXT,:) = XWT(:,IJE+1-NRIMY:IJE+JPHEXT,:)
!!$END IF
!!$!
!!$!*       5.9.4  TH variable
!!$!
!!$!
!!$XLBXTHM(1:NRIMX+JPHEXT,:,:)         = XTHT(1:NRIMX+JPHEXT,:,:)
!!$XLBXTHM(NRIMX+JPHEXT+1:2*NRIMX+2*JPHEXT,:,:) = XTHT(IIE+1-NRIMX:IIE+JPHEXT,:,:)
!!$IF( .NOT. L2D ) THEN
!!$  XLBYTHM(:,1:NRIMY+JPHEXT,:)         = XTHT(:,1:NRIMY+JPHEXT,:)
!!$  XLBYTHM(:,NRIMY+JPHEXT+1:2*NRIMY+2*JPHEXT,:) = XTHT(:,IJE+1-NRIMY:IJE+JPHEXT,:)
!!$END IF
!!$!
!!$!*       5.9.5  TKE variable
!!$!
!!$!
!!$IF (HTURB /= 'NONE') THEN
!!$  XLBXTKEM(1:NRIMX+JPHEXT,:,:)         = XTKET(1:NRIMX+JPHEXT,:,:)
!!$  XLBXTKEM(NRIMX+JPHEXT+1:2*NRIMX+2*JPHEXT,:,:) = XTKET(IIE+1-NRIMX:IIE+JPHEXT,:,:)
!!$  IF( .NOT. L2D ) THEN
!!$    XLBYTKEM(:,1:NRIMY+JPHEXT,:)         = XTKET(:,1:NRIMY+JPHEXT,:)
!!$    XLBYTKEM(:,NRIMY+JPHEXT+1:2*NRIMY+2*JPHEXT,:) = XTKET(:,IJE+1-NRIMY:IJE+JPHEXT,:)
!!$>>>>>>> 1.3.2.4.2.2.2.6.2.3.2.6.2.1
  END IF
!
! <<<<<<< spawn_model2.f90
  CALL MPPDB_CHECKLB(XLBXUM,"SPAWN_MOD2 before SPAWN_SURF2_RAIN",PRECISION,'LBXU',NRIMX)
  CALL MPPDB_CHECKLB(XLBXVM,"SPAWN_MOD2 before SPAWN_SURF2_RAIN:XLBXVM",PRECISION,'LBXU',NRIMX)
  CALL MPPDB_CHECKLB(XLBXWM,"SPAWN_MOD2 before SPAWN_SURF2_RAIN:XLBXWM",PRECISION,'LBXU',NRIMX)
  CALL MPPDB_CHECKLB(XLBYUM,"SPAWN_MOD2 before SPAWN_SURF2_RAIN:XLBYUM",PRECISION,'LBYV',NRIMY)
  CALL MPPDB_CHECKLB(XLBYVM,"SPAWN_MOD2 before SPAWN_SURF2_RAIN:XLBYVM",PRECISION,'LBYV',NRIMY)
  CALL MPPDB_CHECKLB(XLBYWM,"SPAWN_MOD2 before SPAWN_SURF2_RAIN:XLBYWM",PRECISION,'LBYV',NRIMY)
!!$=======
!!$!*       5.9.6  moist variables
!!$!
!!$IF (NRR >= 1) THEN
!!$  DO JRR =1,NRR  
!!$    XLBXRM(1:NRIMX+JPHEXT,:,:,JRR)         = XRT(1:NRIMX+JPHEXT,:,:,JRR)
!!$    XLBXRM(NRIMX+JPHEXT+1:2*NRIMX+2*JPHEXT,:,:,JRR) = XRT(IIE+1-NRIMX:IIE+JPHEXT,:,:,JRR)
!!$    IF( .NOT. L2D ) THEN
!!$      XLBYRM(:,1:NRIMY+JPHEXT,:,JRR)         = XRT(:,1:NRIMY+JPHEXT,:,JRR)
!!$      XLBYRM(:,NRIMY+JPHEXT+1:2*NRIMY+2*JPHEXT,:,JRR) = XRT(:,IJE+1-NRIMY:IJE+JPHEXT,:,JRR)
!!$    END IF
!!$  END DO
!!$END IF
!!$!
!!$!*       5.9.7  scalar variables
!!$!
!!$IF (NSV /= 0) THEN
!!$  DO JSV = 1, NSV
!!$    XLBXSVM(1:NRIMX+JPHEXT,:,:,JSV)         = XSVT(1:NRIMX+JPHEXT,:,:,JSV)
!!$    XLBXSVM(NRIMX+JPHEXT+1:2*NRIMX+2*JPHEXT,:,:,JSV) = XSVT(IIE+1-NRIMX:IIE+JPHEXT,:,:,JSV)
!!$    IF( .NOT. L2D ) THEN
!!$      XLBYSVM(:,1:NRIMY+JPHEXT,:,JSV)         = XSVT(:,1:NRIMY+JPHEXT,:,JSV)
!!$      XLBYSVM(:,NRIMY+JPHEXT+1:2*NRIMY+2*JPHEXT,:,JSV) = XSVT(:,IJE+1-NRIMY:IJE+JPHEXT,:,JSV)
!!$    END IF
!!$  END DO
!!$ENDIF
!!$>>>>>>> 1.3.2.4.2.2.2.6.2.3.2.6.2.1
!
!*       5.10 Surface precipitation computation
!
IF (SIZE(XINPRR) /= 0 ) THEN
  IF (GNOSON) &
    CALL SPAWN_SURF2_RAIN (NXOR,NYOR,NXEND,NYEND,NDXRATIO,NDYRATIO,   &
              XINPRC,XACPRC,XINDEP,XACDEP,XINPRR,XINPRR3D,XEVAP3D,    &
              XACPRR,XINPRS,XACPRS,XINPRG,XACPRG,&
              XINPRH,XACPRH )
  IF (.NOT.GNOSON) &
    CALL SPAWN_SURF2_RAIN (NXOR,NYOR,NXEND,NYEND,NDXRATIO,NDYRATIO,  &
           XINPRC,XACPRC,XINDEP,XACDEP,XINPRR,XINPRR3D,XEVAP3D,      &
           XACPRR,XINPRS,XACPRS,XINPRG,XACPRG,XINPRH,XACPRH,         &
           TZSONFILE,IIUSON,IJUSON,                                  &
           IIB2,IJB2,IIE2,IJE2,                                      &
           IIB1,IJB1,IIE1,IJE1                                       )
ENDIF
!
!*       5.11   Total mass of dry air Md computation :
!
ZTIME1 = ZTIME2
!
ALLOCATE(ZRHOD(IIU,IJU,IKU))
!
IF (LOCEAN) THEN
  ZRHOD(:,:,:)=XRH00OCEAN*(1.-XALPHAOC*(ZTHVT(:,:,:)-XTH00OCEAN)+XBETAOC*(XRT(:,:,:,1)-XSA00OCEAN))
ELSE
  ZRHOD(:,:,:)=XPABST(:,:,:)/(XPABST(:,:,:)/XP00)**(XRD/XCPD) &
              /(XRD*ZTHVT(:,:,:)*(1.+ZSUMRT(:,:,:)))
ENDIF
!$20140709
  CALL MPPDB_CHECK3D(ZRHOD,"SPAWN_MOD2:ZRHOD",PRECISION)
  CALL MPPDB_CHECK3D(XPABST,"SPAWN_MOD2:XPABST",PRECISION)
  CALL MPPDB_CHECK3D(ZSUMRT,"SPAWN_MOD2:ZSUMRT",PRECISION)
!$20140710 until here all ok after UPHALO(XZZ)
!
CALL TOTAL_DMASS(ZJ,ZRHOD,XDRYMASST)
!
DEALLOCATE (ZRHOD)
DEALLOCATE (ZSUMRT,ZTHVT)
!
CALL SECOND_MNH(ZTIME2)
!
ZMISC = ZMISC + ZTIME2 - ZTIME1
!
!*       5.12 Deallocation of model 1 variables : 
!
ZTIME1 = ZTIME2
!
CALL DEALLOCATE_MODEL1(3)
!
CALL SECOND_MNH(ZTIME2)
!
ZMISC = ZMISC + ZTIME2 - ZTIME1
!
!*       5.13 Anelastic correction : 
!
CALL SECOND_MNH(ZTIME1)
!
IF (.NOT. L1D) THEN
  CALL ANEL_BALANCE_n
  CALL BOUNDARIES (                                                 &
            0.,CLBCX,CLBCY,NRR,NSV,1,                               &
            XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,   &
            XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,   &
            XLBXUM,XLBXVM,XLBXWM,XLBXTHM,XLBXTKEM,XLBXRM,XLBXSVM,   &
            XLBYUM,XLBYVM,XLBYWM,XLBYTHM,XLBYTKEM,XLBYRM,XLBYSVM,   &
            XRHODJ,XRHODREF,                                        &
            XUT, XVT, XWT, XTHT, XTKET, XRT, XSVT, XSRCT            )
END IF
!
CALL SECOND_MNH(ZTIME2)
!
ZANEL = ZTIME2 - ZTIME1
!
!
!
!-------------------------------------------------------------------------------
!
!*	 6.    WRITE THE FMFILE
!	       ---------------- 
!
CALL SECOND_MNH(ZTIME1)
!
INPRAR = 22 + 2*(4+NRR+NSV) ! 22 = number of grid variables + reference state
                            ! variables +dimension variables
                            ! 2*(4+NRR+NSV) = number of prognostic variables
                            ! at time t and t-dt
IF ( ( LEN_TRIM(HSPAFILE) /= 0 ) .AND. ( ADJUSTL(HSPAFILE) /= ADJUSTL(CINIFILE) ) ) THEN
  CMY_NAME(2)=HSPAFILE
ELSE
  CMY_NAME(2)=ADJUSTL(ADJUSTR(CINIFILE)//'.spa'//ADJUSTL(HSPANBR))
  IF (.NOT.GNOSON)   &
     CMY_NAME(2)=ADJUSTL(ADJUSTR(CINIFILE)//'.spr'//ADJUSTL(HSPANBR))
END IF
!
CALL IO_File_add2list(TZFILE,CMY_NAME(2),'MNH','WRITE',KLFINPRAR=INPRAR,KLFITYPE=1,KLFIVERB=NVERB)
!
CALL IO_File_open(TZFILE)
!
CALL WRITE_DESFM_n(2,TZFILE)
!
IF (LBAL_ONLY) THEN  ! same relation with its DAD for model2 and for model1
  NDXRATIO_ALL(2) = NDXRATIO_ALL(1)
  NDYRATIO_ALL(2) = NDYRATIO_ALL(1)
  NXOR_ALL(2)     = NXOR_ALL(1)
  NYOR_ALL(2)     = NYOR_ALL(1)
  NXEND_ALL(2)    = NXEND_ALL(1)
  NYEND_ALL(2)    = NYEND_ALL(1)
  CDAD_NAME(2)    = CDAD_NAME(1)
  IF (CDADSPAFILE == '' ) THEN
    IF (NDXRATIO_ALL(1) == 1 .AND. NDYRATIO_ALL(1) == 1 &
      .AND. NXOR_ALL(1) == 1 .AND. NYOR_ALL(1) == 1 ) THEN
      ! for spawning with ratio=1 
      ! if the DAD of model 1 is itself, the DAD of model 2 also.
      CDAD_NAME(2)=CMY_NAME(2)
    ENDIF
  ENDIF
  ! case of model with DAD
  IF (CDADSPAFILE /='') CDAD_NAME(2)=CDADSPAFILE
ELSE
  CDAD_NAME(2)=CMY_NAME(1) ! model 1 becomes the DAD of model 2 (spawned one)
ENDIF
!
CALL IO_Header_write(TZFILE,HDAD_NAME=CDAD_NAME(2))
CALL WRITE_LFIFM_n(TZFILE,CDAD_NAME(2))
!
CALL SECOND_MNH(ZTIME2)
!
ZWRITE = ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*       7.      Surface variables :
!
ZTIME1 = ZTIME2
!
TFILE_SURFEX => TZFILE
CALL SPAWN_SURF(HINIFILE,HINIFILEPGD,TZFILE,OSPAWN_SURF)
NULLIFY(TFILE_SURFEX)
!
CALL SECOND_MNH(ZTIME2)
!
ZSURF2 = ZTIME2 - ZTIME1
!
!-------------------------------------------------------------------------------
!
!*	 8.    CLOSES THE FMFILE
!	       ----------------- 
!
CALL IO_File_close(TZFILE)
IF (ASSOCIATED(TZSONFILE)) THEN
  CALL IO_File_close(TZSONFILE)
END IF
!
!-------------------------------------------------------------------------------
!
!*       9.    PRINTS ON OUTPUT-LISTING
!              ------------------------
!
WRITE(ILUOUT,FMT=9900) XZHAT(1)
!
DO JLOOP = 2,IKU
 WRITE(ILUOUT,FMT=9901) JLOOP,XZHAT(JLOOP),XZHAT(JLOOP)-XZHAT(JLOOP-1)
END DO
!
IF (NVERB >= 5) THEN
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: LUSERV,LUSERC=',LUSERV,LUSERC 
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: LUSERR,LUSERI,LUSERS=',LUSERR,LUSERI,LUSERS
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: LUSERG,LUSERH,NSV=',LUSERG,LUSERH,NSV
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: NRR=',NRR
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: NVERB=',NVERB
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: XLON0,XLAT0,XBETA=',XLON0,XLAT0,XBETA
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: LCARTESIAN=',LCARTESIAN
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: LOCEAN,LCOUPLES=',LOCEAN,LCOUPLES
  IF(LCARTESIAN) THEN
    WRITE(ILUOUT,*) 'SPAWN_MODEL2: No map projection used.'
  ELSE
    WRITE(ILUOUT,*) 'SPAWN_MODEL2: XRPK,XLONORI,XLATORI=',XRPK,XLONORI,XLATORI
    IF (ABS(XRPK) == 1.) THEN
      WRITE(ILUOUT,*) 'SPAWN_MODEL2: Polar stereo used.'
    ELSE IF (XRPK == 0.) THEN
      WRITE(ILUOUT,*) 'SPAWN_MODEL2: Mercator used.'
    ELSE
      WRITE(ILUOUT,*) 'SPAWN_MODEL2: Lambert used, cone factor=',XRPK 
    END IF
  END IF
END IF
!
IF (NVERB >= 10) THEN
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: IIB, IJB, IKB=',IIB,IJB,IKB
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: IIU, IJU, IKU=',IIU,IJU,IKU
END IF
!
IF(NVERB >= 10) THEN                               !Value control
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: Some XZS values:'
  WRITE(ILUOUT,*) XZS(1,IJU),XZS((IIU-1)/2,IJU),XZS(IIU,IJU)  
  WRITE(ILUOUT,*) XZS(1,(IJU-1)/2),XZS((IIU-1)/2,(IJU-1)/2),XZS(IIU,(IJU-1)/2)  
  WRITE(ILUOUT,*) XZS(1,1)      ,XZS((IIU-1)/2,1)      ,XZS(IIU,1)  
END IF
!
IF(NVERB >= 10) THEN       !Value control
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: Some XUT values:'
  WRITE(ILUOUT,*) ' (1,IJU/2,JK) (IIU/2,1,JK) (IIU/2,IJU/2,JK) &
                   &(IIU/2,IJU,JK) (IIU,IJU/2,JK)'
  DO JKLOOP=1,IKU
    WRITE(ILUOUT,*) 'JK = ',JKLOOP
    WRITE(ILUOUT,*) XUT(1,IJU/2,JKLOOP),XUT(IIU/2,1,JKLOOP),       &
                    XUT(IIU/2,IJU/2,JKLOOP),XUT(IIU/2,IJU,JKLOOP), &
                    XUT(IIU,IJU/2,JKLOOP)
  END DO
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: Some XVT values:'
  WRITE(ILUOUT,*) ' (1,IJU/2,JK) (IIU/2,1,JK) (IIU/2,IJU/2,JK) &
                   &(IIU/2,IJU,JK) (IIU,IJU/2,JK)'
  DO JKLOOP=1,IKU
    WRITE(ILUOUT,*) 'JK = ',JKLOOP
    WRITE(ILUOUT,*) XVT(1,IJU/2,JKLOOP),XVT(IIU/2,1,JKLOOP),       &
                    XVT(IIU/2,IJU/2,JKLOOP),XVT(IIU/2,IJU,JKLOOP), &
                    XVT(IIU,IJU/2,JKLOOP)
  END DO
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: Some XWT values:'
  WRITE(ILUOUT,*) ' (1,IJU/2,JK) (IIU/2,1,JK) (IIU/2,IJU/2,JK) &
                   &(IIU/2,IJU,JK) (IIU,IJU/2,JK)'
  DO JKLOOP=1,IKU
    WRITE(ILUOUT,*) 'JK = ',JKLOOP
    WRITE(ILUOUT,*) XWT(1,IJU/2,JKLOOP),XWT(IIU/2,1,JKLOOP),       &
                    XWT(IIU/2,IJU/2,JKLOOP),XWT(IIU/2,IJU,JKLOOP), &
                    XWT(IIU,IJU/2,JKLOOP)
  END DO
  WRITE(ILUOUT,*) 'SPAWN_MODEL2: Some XTHT values:'
  WRITE(ILUOUT,*) ' (1,IJU/2,JK) (IIU/2,1,JK) (IIU/2,IJU/2,JK) &
                   &(IIU/2,IJU,JK) (IIU,IJU/2,JK)'
  DO JKLOOP=1,IKU
    WRITE(ILUOUT,*) 'JK = ',JKLOOP
    WRITE(ILUOUT,*) XTHT(1,IJU/2,JKLOOP),XTHT(IIU/2,1,JKLOOP),       &
                    XTHT(IIU/2,IJU/2,JKLOOP),XTHT(IIU/2,IJU,JKLOOP), &
                    XTHT(IIU,IJU/2,JKLOOP)
  END DO
  IF(NRR >= 1) THEN
    WRITE(ILUOUT,*) 'SPAWN_MODEL2: Some XRT values:'
    WRITE(ILUOUT,*) ' (1,IJU/2,JK) (IIU/2,1,JK) (IIU/2,IJU/2,JK) &
                     &(IIU/2,IJU,JK) (IIU,IJU/2,JK)'
    DO JKLOOP=1,IKU
      WRITE(ILUOUT,*) 'JK = ',JKLOOP
      WRITE(ILUOUT,*) XRT(1,IJU/2,JKLOOP,1),XRT(IIU/2,1,JKLOOP,1),       &
                      XRT(IIU/2,IJU/2,JKLOOP,1),XRT(IIU/2,IJU,JKLOOP,1), &
                      XRT(IIU,IJU/2,JKLOOP,1)
    END DO
  END IF
  !
  IF (LUV_FLX) THEN
    WRITE(ILUOUT,*)'SPAWN_MODEL2: Some EDDY_FLUX values XVU_FLUX(IIU/2,2,:)=',XVU_FLUX_M(IIU/2,2,:)
  END IF
  !
  IF (LTH_FLX) THEN
    WRITE(ILUOUT,*)'SPAWN_MODEL2: Some EDDY_FLUX values XVTH_FLUX(IIU/2,2,:)=',XVTH_FLUX_M(IIU/2,2,:)
    WRITE(ILUOUT,*)'SPAWN_MODEL2: Some EDDY_FLUX values XWTH_FLUX(IIU/2,2,:)=',XWTH_FLUX_M(IIU/2,2,:)
  END IF
  !
END IF
!
WRITE(ILUOUT,*) 'SPAWN_MODEL2: SPAWN_MODEL2 ENDS CORRECTLY.'
!
CALL SECOND_MNH (ZEND)
!
ZTOT        = ZEND - ZSTART          ! for computing time analysis
!
ZALL = ZGRID2 + ZSURF2 + ZMISC + ZFIELD2 + ZVER + ZPRESSURE2 + ZANEL + ZWRITE
!
ZPERCALL  = 100.*ZALL/ZTOT
!
ZPERCGRID2  = 100.*ZGRID2/ZTOT
ZPERCSURF2  = 100.*ZSURF2/ZTOT
ZPERCMISC   = 100.*ZMISC/ZTOT
ZPERCFIELD2 = 100.*ZFIELD2/ZTOT
ZPERCVER    = 100.*ZVER/ZTOT
ZPERCPRESSURE2 = 100.*ZPRESSURE2/ZTOT
ZPERCANEL   = 100.*ZANEL/ZTOT
ZPERCWRITE  = 100.*ZWRITE/ZTOT
!
WRITE(ILUOUT,*)
WRITE(ILUOUT,*) ' ------------------------------------------------------------ '
WRITE(ILUOUT,*) '|                                                            |'
WRITE(ILUOUT,*) '|           COMPUTING TIME ANALYSIS in SPAWN_MODEL2          |'
WRITE(ILUOUT,*) '|                                                            |'
WRITE(ILUOUT,*) '|------------------------------------------------------------|'
WRITE(ILUOUT,*) '|                     |                   |                  |'
WRITE(ILUOUT,*) '|    ROUTINE NAME     |     CPU-TIME      |   PERCENTAGE %   |'
WRITE(ILUOUT,*) '|                     |                   |                  |'
WRITE(ILUOUT,*) '|---------------------|-------------------|------------------|'
WRITE(ILUOUT,*) '|                     |                   |                  |'
WRITE(UNIT=ILUOUT,FMT=1) ZGRID2 ,ZPERCGRID2
WRITE(UNIT=ILUOUT,FMT=3) ZFIELD2,ZPERCFIELD2
WRITE(UNIT=ILUOUT,FMT=8) ZVER,ZPERCVER
WRITE(UNIT=ILUOUT,FMT=7) ZPRESSURE2,ZPERCPRESSURE2
WRITE(UNIT=ILUOUT,FMT=2) ZSURF2 ,ZPERCSURF2
WRITE(UNIT=ILUOUT,FMT=4) ZANEL  ,ZPERCANEL
WRITE(UNIT=ILUOUT,FMT=5) ZWRITE ,ZPERCWRITE
WRITE(UNIT=ILUOUT,FMT=9) ZMISC  ,ZPERCMISC
WRITE(UNIT=ILUOUT,FMT=6) ZTOT   ,ZPERCALL
WRITE(ILUOUT,*) ' ------------------------------------------------------------ '
!
!                  FORMATS
!                  -------
!
1  FORMAT(' |    SPAWN_GRID2      |     ',F8.3,'      |     ',F8.3,'     |')
3  FORMAT(' |    SPAWN_FIELD2     |     ',F8.3,'      |     ',F8.3,'     |')
8  FORMAT(' |  VER_INTERP_FIELD   |     ',F8.3,'      |     ',F8.3,'     |')
7  FORMAT(' |  SPAWN_PRESSURE2    |     ',F8.3,'      |     ',F8.3,'     |')
2  FORMAT(' |    SPAWN_SURF2      |     ',F8.3,'      |     ',F8.3,'     |')
4  FORMAT(' |   ANEL_BALANCE2     |     ',F8.3,'      |     ',F8.3,'     |')
5  FORMAT(' |      WRITE          |     ',F8.3,'      |     ',F8.3,'     |')
9  FORMAT(' |   MISCELLANEOUS     |     ',F8.3,'      |     ',F8.3,'     |')
6  FORMAT(' |    SPAWN_MODEL2     |     ',F8.3,'      |     ',F8.3,'     |')
!
!
CALL IO_File_close(TLUOUT)
!
9900  FORMAT(' K = 001    ZHAT = ',E14.7)
9901  FORMAT(' K = ',I3.3,'    ZHAT = ',E14.7,'    DZ = ' ,E14.7)
!
!-------------------------------------------------------------------------------
!
!
! Switch back to model index of calling routine
CALL GOTO_MODEL(IMI)
!
END SUBROUTINE SPAWN_MODEL2
