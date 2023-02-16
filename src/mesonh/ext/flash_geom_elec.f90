!MNH_LIC Copyright 2010-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     #############################
      MODULE MODI_FLASH_GEOM_ELEC_n
!     #############################
!
INTERFACE
    SUBROUTINE FLASH_GEOM_ELEC_n (KTCOUNT, KMI, KRR, PTSTEP, OEXIT,                      &
                                  PRHODJ, PRHODREF, PRT, PCIT, PRSVS, PRS, PTHT, PPABST, &
                                  PEFIELDU, PEFIELDV, PEFIELDW, PZZ, PSVS_LINOX,         &
                                  TPFILE_FGEOM_DIAG, TPFILE_FGEOM_COORD, TPFILE_LMA,     &
                                  PTOWN, PSEA                                            )
!
USE MODD_IO, ONLY: TFILEDATA
!
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
INTEGER,                  INTENT(IN)    :: KMI      ! current model index
INTEGER,                  INTENT(IN)    :: KRR      ! number of moist variables
REAL,                     INTENT(IN)    :: PTSTEP   ! Double time step except for
                                                    ! cold start
LOGICAL,                  INTENT(IN)    :: OEXIT    ! switch for the end of the temporal loop
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ   ! Dry density * Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT      ! Moist variables at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT     ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS    ! Scalar variables source term
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDU ! x-component of the electric field
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDV ! y-component of the electric field
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDW ! z-component of the electric field
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS      ! Moist variables vol. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta (K) at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST   ! Absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ      ! height
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSVS_LINOX ! NOx source term
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE_FGEOM_DIAG
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE_FGEOM_COORD
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE_LMA
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN ! town fraction
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA  ! Land-sea mask
!
END SUBROUTINE FLASH_GEOM_ELEC_n
END INTERFACE
END MODULE MODI_FLASH_GEOM_ELEC_n
!
!
!   ######################################################################################
    SUBROUTINE FLASH_GEOM_ELEC_n (KTCOUNT, KMI, KRR, PTSTEP, OEXIT,                      &
                                  PRHODJ, PRHODREF, PRT, PCIT, PRSVS, PRS, PTHT, PPABST, &
                                  PEFIELDU, PEFIELDV, PEFIELDW, PZZ, PSVS_LINOX,         &
                                  TPFILE_FGEOM_DIAG, TPFILE_FGEOM_COORD, TPFILE_LMA,     &
                                  PTOWN, PSEA                                            )
!   ######################################################################################
!
!!****  * -
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to compute the lightning flash path,
!!    and to neutralize the electric charge along the lightning channel.
!!
!!
!!    METHOD
!!    ------
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
!!      C. Barthe   * LACy * 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original : Jan. 2010
!!      Modifications:
!!      M. Chong  * LA *  Juin 2010 : add small ions
!!      J-P Pinty * LA *  Feb. 2013 : add LMA storage
!!      J-P Pinty * LA *  Nov. 2013 : add flash map storage
!!      M. Chong  * LA *  Juin 2010 : add LiNOx
!!      C. Barthe * LACy * Jan. 2015 : convert trig. pt into lat,lon in ascii file
!!      J.Escobar : 18/12/2015 : Correction of bug in bound in // for NHALO <>1 
!!      J.Escobar : 28/03/2018 : Correction of multiple // bug & compiler indepedent mnh_random_number
!!      J.Escobar : 20/06/2018 : Correction of computation of global index I8VECT
!!      J.Escobar : 10/12/2018 : // Correction , mpi_bcast CG & CG_POS parameter 
!!                               & initialize INBLIGHT on all proc for filling/saving AREA* arrays
!  P. Wautelet 10/01/2019: use NEWUNIT argument of OPEN
!  P. Wautelet 22/01/2019: use standard FLUSH statement instead of non standard intrinsics!!
!  P. Wautelet 22/02/2019: use MOD intrinsics with same kind for all arguments (to respect Fortran standard)
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 19/04/2019: use modd_precision kinds
!  P. Wautelet 26/04/2019: use modd_precision parameters for datatypes of MPI communications
!  P. Wautelet 26/04/2019: replace non-standard FLOAT function by REAL function
!  P. Wautelet 20/05/2019: add name argument to ADDnFIELD_ll + new ADD4DFIELD_ll subroutine
!  P. Wautelet 18/09/2019: correct support of 64bit integers (MNH_INT=8)
!-------------------------------------------------------------------------------
!
!*      0.      DECLARATIONS
!               ------------
!
USE MODD_ARGSLIST_ll,    ONLY: LIST_ll
USE MODD_CONF,           ONLY: CEXP, LCARTESIAN
USE MODD_CST,            ONLY: XAVOGADRO, XMD
USE MODD_DYN_n,          ONLY: XDXHATM, XDYHATM, NSTOP
USE MODD_ELEC_DESCR
USE MODD_ELEC_FLASH
USE MODD_ELEC_PARAM,     ONLY: XFQLIGHTR, XEXQLIGHTR, &
                               XFQLIGHTI, XEXQLIGHTI, &
                               XFQLIGHTS, XEXQLIGHTS, &
                               XFQLIGHTG, XEXQLIGHTG, &
                               XFQLIGHTH, XEXQLIGHTH, &
                               XFQLIGHTC
USE MODD_GRID,           ONLY: XLATORI,XLONORI
USE MODD_GRID_n,         ONLY: XXHAT, XYHAT, XZHAT
USE MODD_IO,             ONLY: TFILEDATA
USE MODD_LMA_SIMULATOR
USE MODD_METRICS_n,      ONLY: XDXX, XDYY, XDZZ ! in linox_production
USE MODD_NSV,            ONLY: NSV_ELECBEG, NSV_ELECEND, NSV_ELEC
USE MODD_PARAMETERS,     ONLY: JPHEXT, JPVEXT
use MODD_PRECISION,      only: MNHINT_MPI, MNHLOG_MPI, MNHREAL_MPI
USE MODD_RAIN_ICE_DESCR, ONLY: XLBR, XLBEXR, XLBS, XLBEXS, &
                               XLBG, XLBEXG, XLBH, XLBEXH, &
                               XRTMIN
USE MODD_SUB_ELEC_n
USE MODD_TIME_n
USE MODD_VAR_ll,         ONLY: NPROC,NMNH_COMM_WORLD
!
USE MODE_ELEC_ll
USE MODE_GRIDPROJ
USE MODE_ll
USE MODE_MPPDB
#ifdef MNH_PGI
USE MODE_PACK_PGI
#endif
!
USE MODI_ION_ATTACH_ELEC
USE MODI_SHUMAN
USE MODI_TO_ELEC_FIELD_n
!
IMPLICIT NONE
!
!
!       0.1     Declaration of arguments
!
INTEGER,                  INTENT(IN)    :: KTCOUNT  ! Temporal loop counter
INTEGER,                  INTENT(IN)    :: KMI      ! current model index
INTEGER,                  INTENT(IN)    :: KRR      ! number of moist variables
REAL,                     INTENT(IN)    :: PTSTEP   ! Double time step except for
                                                    ! cold start
LOGICAL,                  INTENT(IN)    :: OEXIT    ! switch for the end of the temporal loop
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODREF ! Reference dry air density
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PRHODJ   ! Dry density * Jacobian
REAL, DIMENSION(:,:,:,:), INTENT(IN)    :: PRT      ! Moist variables at time t
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PCIT     ! Pristine ice n.c. at t
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRSVS    ! Scalar variables source term
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDU ! x-component of the electric field
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDV ! y-component of the electric field
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PEFIELDW ! z-component of the electric field
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRS      ! Moist variables vol. source
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PTHT     ! Theta (K) at time t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PPABST   ! Absolute pressure at t
REAL, DIMENSION(:,:,:),   INTENT(IN)    :: PZZ      ! height
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PSVS_LINOX ! NOx source term
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE_FGEOM_DIAG
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE_FGEOM_COORD
TYPE(TFILEDATA),          INTENT(IN)    :: TPFILE_LMA
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PTOWN ! town fraction
REAL, DIMENSION(:,:), OPTIONAL, INTENT(IN) :: PSEA  ! Land-sea mask
!
!
!       0.2     Declaration of local variables
!
INTEGER :: IIB, IIE  ! index values of the first and last inner mass points along x
INTEGER :: IJB, IJE  ! index values of the first and last inner mass points along y
INTEGER :: IKB, IKE  ! index values of the first and last inner mass points along z
INTEGER :: II, IJ, IK, IL, IM, IPOINT  ! loop indexes
INTEGER :: IX, IY, IZ
INTEGER :: IXOR, IYOR  ! origin of the extended subdomain
INTEGER :: INB_CELL    ! Number of detected electrified cells
INTEGER :: IPROC_CELL  ! Proc with the center of the cell
INTEGER :: IICOORD, IJCOORD, IKCOORD ! local indexes of the cell center / max electric field
INTEGER :: IPROC       ! my proc number
INTEGER :: IINFO_ll    ! return code of parallel routine
INTEGER :: COUNT_BEF   ! nb of pts in zcell before testing neighbour pts
INTEGER :: COUNT_AFT   ! nb of pts in zcell after testing neighbour pts
INTEGER :: INBFTS_MAX  ! Max number of flashes per time step / cell 
INTEGER :: IIBL_LOC    ! local i index of the ongoing bi-leader segment
INTEGER :: IJBL_LOC    ! local j index of the ongoing bi-leader segment
INTEGER :: IKBL        ! k index of the ongoing bi-leader segment
INTEGER :: II_TRIG_LOC  ! local i index of the triggering point
INTEGER :: IJ_TRIG_LOC  ! local j index of the triggering point
INTEGER :: II_TRIG_GLOB ! global i index of the potential triggering pt
INTEGER :: IJ_TRIG_GLOB ! global j index of the potential triggering pt
INTEGER :: IK_TRIG      ! k index of the triggering point
INTEGER :: ISIGN_LEADER ! sign of the leader
INTEGER :: IPROC_AUX    ! proc number for max_ll and min_ll
INTEGER :: IIND_MAX   ! max nb of indexes between the trig. pt and the possible branches
INTEGER :: IIND_MIN   ! min nb of indexes between the trig. pt and the possible branches
INTEGER :: IDELTA_IND   ! number of indexes between iind_max and iind_min
INTEGER :: IPT_DIST     ! nb of possible pts for branching on each proc
INTEGER :: IPT_DIST_GLOB  ! global nb of possible pts for branching
INTEGER :: IFOUND         ! if =1, then the random selection is successful 
INTEGER :: ICHOICE_LOCX   ! local i indice for random choice
INTEGER :: ICHOICE_LOCY   ! local j indice for random choice
INTEGER :: ICHOICE_Z      !       k indice for random choice
INTEGER :: INB_PROP       ! nb of pts where the flash can propagate
INTEGER :: INB_NEUT       ! nb of pts to neutralize
INTEGER :: INB_NEUT_OK    ! nb of effective flash neutralization
INTEGER :: ISTOP
INTEGER :: IERR         ! error status
INTEGER :: IWORK
INTEGER :: ICHOICE
INTEGER :: IIMIN, IIMAX, IJMIN, IJMAX, IKMIN, IKMAX
INTEGER :: IPOS_LEADER, INEG_LEADER
INTEGER :: INBLIGHT
INTEGER, DIMENSION(:), ALLOCATABLE, SAVE :: ITYPE   ! flash type (IC, CGN or CGP)
INTEGER, DIMENSION(:), ALLOCATABLE :: INBSEG_LEADER ! number of segments in the leader
INTEGER, DIMENSION(:), ALLOCATABLE :: ISIGNE_EZ     ! sign of the vertical electric field 
                                                    ! component at the trig. pt
INTEGER, DIMENSION(:), ALLOCATABLE :: IPROC_TRIG    ! proc that contains the triggering point
INTEGER, DIMENSION(:), ALLOCATABLE :: INBSEG        ! Number of segments per flash
INTEGER, DIMENSION(:), ALLOCATABLE :: INBSEG_ALL    ! Number of segments, all processes
INTEGER, DIMENSION(NPROC)          :: INBSEG_PROC   ! ------------------ per process
INTEGER, DIMENSION(:), ALLOCATABLE :: INB_FLASH     ! Number of flashes per time step / cell
INTEGER, DIMENSION(:), ALLOCATABLE :: INB_FL_REAL   ! Effective Number of flashes per timestep/cell
INTEGER, DIMENSION(:), ALLOCATABLE :: IHIST_LOC     ! local nb of possible branches at [r,r+dr]
INTEGER, DIMENSION(:), ALLOCATABLE :: IHIST_GLOB    ! global nb of possible branches at [r,r+dr]
                                                    ! at [r,r+dr] on each proc
INTEGER, DIMENSION(:), ALLOCATABLE :: IMAX_BRANCH   ! max nb of branches at [r,r+dr]
                                                    ! proportional to the percentage of 
                                                    ! available pts / proc at this distance
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISEG_LOC    ! Local indexes of the flash segments
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ICELL_LOC   ! local indexes + proc of the cell 'center'
INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: IMASKQ_DIST ! contains the distance/indice 
                                                      ! from the triggering pt
!
LOGICAL :: GPOSITIVE    ! if T, positive charge regions where the negative part 
                        ! of the leader propagates
LOGICAL :: GEND_DOMAIN  ! no more points with E > E_threshold
LOGICAL :: GEND_CELL    ! if T, end of the cell
LOGICAL :: GCG          ! if true, the flash is a CG
LOGICAL :: GCG_POS      ! if true, the flash is a +CG
LOGICAL :: GNEUTRALIZATION
LOGICAL :: GNEW_FLASH_GLOB
LOGICAL, DIMENSION(:), ALLOCATABLE :: GNEW_FLASH
LOGICAL, DIMENSION(:,:,:),   ALLOCATABLE :: GATTACH  ! if T, ion recombination and
                                                     ! attachment
LOGICAL, DIMENSION(:,:,:),   ALLOCATABLE :: GPOSS    ! if T, new cell possible at this pt
LOGICAL, DIMENSION(:,:,:,:), ALLOCATABLE :: GPROP    ! if T, propagation possible at this pt
!
REAL :: ZE_TRIG_THRES ! Triggering Electric field threshold corrected for
                      !  pressure   
REAL :: ZMAXE         ! Max electric field module (V/m)
REAL :: ZEMOD_BL      ! E module at the tip of the last segment of the leader (V/m)
REAL :: ZMEAN_GRID    ! mean grid size
REAL :: ZMAX_DIST     ! max distance between the triggering pt and the possible branches
REAL :: ZMIN_DIST     ! min distance between the triggering pt and the possible branches
REAL :: ZRANDOM       ! random number
REAL :: ZQNET         ! net charge carried by the flash (C/kg)
REAL :: ZCLOUDLIM     ! cloud limit
REAL :: ZSIGMIN       ! min efficient cross section
REAL :: ZLAT, ZLON    ! lat,lon coordinates of the triggering points if not lcartesian
!
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZQMT   ! mass charge density (C/kg)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZCELL  ! define the electrified cells
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZSIGMA ! efficient cross section of hydrometeors
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZDQDT  ! charge to neutralize at each pt (C/kg)
REAL, DIMENSION(:,:,:,:), ALLOCATABLE :: ZFLASH ! = 1 if the flash leader reaches this pt
                                                ! = 2 if the flash branch is concerned
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLBDAR   ! Lambda for rain
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLBDAS   ! Lambda for snow
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLBDAG   ! Lambda for graupel
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLBDAH   ! Lambda for hail
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZQMTOT   ! total mass charge density (C/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZCLOUD   ! total mixing ratio (kg/kg)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZEMODULE ! Electric field module (V/m)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZDIST    ! distance between the trig. pt and the cell pts (m)
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZSIGLOB  ! sum of the cross sections
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZQFLASH  ! total charge in excess of xqexcess (C/kg)
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOORD_TRIG ! Global coordinates of triggering point
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOORD_SEG ! Global coordinates of segments
REAL, DIMENSION(:), ALLOCATABLE :: ZEM_TRIG     ! Electric field module at the triggering pt
REAL, DIMENSION(:), ALLOCATABLE :: ZNEUT_POS    ! Positive charge neutralized at each segment
REAL, DIMENSION(:), ALLOCATABLE :: ZNEUT_NEG    ! Negative charge neutralized at each segment
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ISEG_GLOB  ! Global indexes of LMA segments
INTEGER, DIMENSION(:,:), ALLOCATABLE :: ILMA_SEG_ALL   ! Global indexes of LMA segments
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLMA_QMT ! Particle charge at neutralization point
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLMA_PRT ! Particle mixing ratio at neutralization point
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLMA_NEUT_POS
REAL, DIMENSION(:,:), ALLOCATABLE :: ZLMA_NEUT_NEG
REAL, DIMENSION(:,:), ALLOCATABLE :: ZCOORD_SEG_ALL
REAL, DIMENSION(:), ALLOCATABLE :: ZEMAX        ! Max electric field in each cell
REAL, DIMENSION(:), ALLOCATABLE :: ZHIST_PERCENT ! percentage of possible branches at [r,r+dr] on each proc
REAL, DIMENSION(:), ALLOCATABLE :: ZMAX_BRANCH  ! max nb of branches at [r,r+dr]
REAL, DIMENSION(:), ALLOCATABLE :: ZVECT
!
! Storage for nflash_write flashes before writing output files (denoted xSxxx)
INTEGER, SAVE :: ISAVE_STATUS ! 0: print and save
                              ! 1: save only
                              ! 2: print only
!
TYPE(LIST_ll), POINTER :: TZFIELDS_ll=> NULL()   ! list of fields to exchange
!
! Storage for the localization of the flashes
LOGICAL :: GFIRSTFLASH
INTEGER,DIMENSION(SIZE(PRT,1),SIZE(PRT,2)) :: IMAP2D
!
!  Storage for the NOx production terms
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZLNOX
REAL    :: ZLGHTLENGTH, ZCOEF
INTEGER :: IFLASH_COUNT, IFLASH_COUNT_GLOB  ! Total number of flashes within the timestep
!
REAL,DIMENSION(SIZE(PRT,1),SIZE(PRT,2)) :: ZCELL_NEW
!
INTEGER :: ILJ
INTEGER :: NIMAX_ll, NJMAX_ll,IIU_ll,IJU_ll    ! dimensions of global domain
!
!-------------------------------------------------------------------------------
!
!*      1.      INITIALIZATION
!               --------------
CALL MYPROC_ELEC_ll(IPROC)
!
!*      1.1     subdomains indexes
!
! beginning and end indexes of the physical subdomain
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB = 1 + JPVEXT
IKE = SIZE(PRT,3) - JPVEXT
!
! global indexes of the local subdomains origin
CALL GET_GLOBALDIMS_ll (NIMAX_ll,NJMAX_ll)
CALL GET_OR_ll('B',IXOR,IYOR)
IIU_ll = NIMAX_ll + 2*JPHEXT
IJU_ll = NJMAX_ll + 2*JPHEXT
!
!
!*      1.2     allocations and initializations
!
!
! from the litterature, the max number of flash per minute is ~ 1000
! this value is used here as the max number of flash per minute per cell
INBFTS_MAX = ANINT(1000 * PTSTEP / 60)
!
IF (GEFIRSTCALL) THEN
  GEFIRSTCALL = .FALSE.
  ALLOCATE (ZXMASS(SIZE(XXHAT)))
  ALLOCATE (ZYMASS(SIZE(XYHAT)))
  ALLOCATE (ZZMASS(SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3)))
  ALLOCATE (ZPRES_COEF(SIZE(PZZ,1), SIZE(PZZ,2), SIZE(PZZ,3)))
  IF(LLMA) THEN
    ALLOCATE (ZLMA_LAT(NFLASH_WRITE, NBRANCH_MAX))
    ALLOCATE (ZLMA_LON(NFLASH_WRITE, NBRANCH_MAX))
    ALLOCATE (ZSLMA_NEUT_POS(NFLASH_WRITE, NBRANCH_MAX))
    ALLOCATE (ZSLMA_NEUT_NEG(NFLASH_WRITE, NBRANCH_MAX))
    ALLOCATE (ISLMA_SEG_GLOB(NFLASH_WRITE, NBRANCH_MAX, 3))
    ALLOCATE (ZSLMA_QMT(NFLASH_WRITE, NBRANCH_MAX, SIZE(PRSVS,4)))
    ALLOCATE (ZSLMA_PRT(NFLASH_WRITE, NBRANCH_MAX, SIZE(PRSVS,4)))
    ISLMA_SEG_GLOB(:,:,:) = 0
  END IF
  ALLOCATE (ZSCOORD_SEG(NFLASH_WRITE, NBRANCH_MAX, 3))  ! NFLASH_WRITE nb of flash to be stored
                                            ! before writing in files
                                            ! NBRANCH_MAX=5000 default
  ALLOCATE (ISFLASH_NUMBER(0:NFLASH_WRITE))
  ALLOCATE (ISNB_FLASH(NFLASH_WRITE))
  ALLOCATE (ISCELL_NUMBER(NFLASH_WRITE))
  ALLOCATE (ISNBSEG(NFLASH_WRITE))
  ALLOCATE (ISTCOUNT_NUMBER(NFLASH_WRITE))
  ALLOCATE (ISTYPE(NFLASH_WRITE))
  ALLOCATE (ZSEM_TRIG(NFLASH_WRITE))
  ALLOCATE (ZSNEUT_POS(NFLASH_WRITE))
  ALLOCATE (ZSNEUT_NEG(NFLASH_WRITE))
!
  ZXMASS(IIB:IIE) = 0.5 * (XXHAT(IIB:IIE) + XXHAT(IIB+1:IIE+1))
  ZYMASS(IJB:IJE) = 0.5 * (XYHAT(IJB:IJE) + XYHAT(IJB+1:IJE+1))
  ZZMASS = MZF(PZZ)
  ZPRES_COEF = EXP(ZZMASS/8400.)
  ZSCOORD_SEG(:,:,:) = 0.0
  ISAVE_STATUS = 1
  ISFLASH_NUMBER(:) = 0
END IF
!
ALLOCATE (ZQMT(SIZE(PRSVS,1),SIZE(PRSVS,2),SIZE(PRSVS,3),SIZE(PRSVS,4)))
ALLOCATE (ZQMTOT(SIZE(PRSVS,1),SIZE(PRSVS,2),SIZE(PRSVS,3)))
ALLOCATE (ZCLOUD(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
ALLOCATE (GPOSS(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
ALLOCATE (ZEMODULE(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
ALLOCATE (ZCELL(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),NMAX_CELL))

!
ZQMT(:,:,:,:) = 0.
ZQMTOT(:,:,:) = 0.
ZCLOUD(:,:,:) = 0.
GPOSS(:,:,:) = .FALSE.
GPOSS(IIB:IIE,IJB:IJE,IKB:IKE) = .TRUE.
ZEMODULE(:,:,:) = 0.
ZCELL(:,:,:,:) = 0.
!
!
!*      1.3     point discharge (Corona)
!
PRSVS(:,:,:,1) = XECHARGE * PRSVS(:,:,:,1)             ! C /(m3 s)
PRSVS(:,:,:,NSV_ELEC) = -1. * XECHARGE * PRSVS(:,:,:,NSV_ELEC)  ! C /(m3 s)
!
CALL PT_DISCHARGE
!
!
!*      1.4     total charge density and mixing ratio
!
DO II = 1, NSV_ELEC
! transform the source term (C/s) into the updated charge density (C/kg)
  ZQMT(:,:,:,II) = PRSVS(:,:,:,II) * PTSTEP / PRHODJ(:,:,:)
!
! total mass charge density (C/kg)
  ZQMTOT(:,:,:)  = ZQMTOT(:,:,:) + PRSVS(:,:,:,II) * PTSTEP / PRHODJ(:,:,:)
END DO
!
! total mixing ratio (g/kg)
DO II = 2, KRR
  ZCLOUD(:,:,:) = ZCLOUD(:,:,:) + PRT(:,:,:,II)
END DO
!
!
!*      1.5     constants
!
ZCLOUDLIM = 1.E-5 
ZSIGMIN   = 1.E-12
!
!
!-------------------------------------------------------------------------------
!
!*      2.      FIND AND COUNT THE ELECTRIFIED CELLS
!               ------------------------------------
!
ALLOCATE (ZEMAX(NMAX_CELL))
ALLOCATE (ICELL_LOC(4,NMAX_CELL))
!
ZEMAX(:) = 0.
ICELL_LOC(:,:) = 0
!
WHERE (ZCLOUD(IIB:IIE,IJB:IJE,IKB:IKE) .LE. ZCLOUDLIM)
  GPOSS(IIB:IIE,IJB:IJE,IKB:IKE) = .FALSE.
END WHERE
!
!
!*      2.1     find the maximum electric field
!
GEND_DOMAIN = .FALSE.
GEND_CELL = .FALSE.
INB_CELL = 0
ZE_TRIG_THRES = XETRIG * (1. - XEBALANCE)
!
CALL MPPDB_CHECK3DM("flash:: PRHODJ,PRT",PRECISION,&
     PRHODJ,PRT(:,:,:,1),PRT(:,:,:,2),PRT(:,:,:,3),PRT(:,:,:,4),&
     PRT(:,:,:,5),PRT(:,:,:,6))
CALL MPPDB_CHECK3DM("flash:: ZQMT",PRECISION,&
     ZQMT(:,:,:,1),ZQMT(:,:,:,2),ZQMT(:,:,:,3),ZQMT(:,:,:,4),&
     ZQMT(:,:,:,5),ZQMT(:,:,:,6),ZQMT(:,:,:,7))

CALL TO_ELEC_FIELD_n (PRT, ZQMT, PRHODJ, KTCOUNT, KRR, &
                      PEFIELDU, PEFIELDV, PEFIELDW)
CALL MPPDB_CHECK3DM("flash:: PEFIELDU, PEFIELDV, PEFIELDW",PRECISION,&
                     PEFIELDU, PEFIELDV, PEFIELDW)
!
! electric field module including pressure effect
ZEMODULE(IIB:IIE,IJB:IJE,IKB:IKE) = ZPRES_COEF(IIB:IIE,IJB:IJE,IKB:IKE)*    &
                                    (PEFIELDU(IIB:IIE,IJB:IJE,IKB:IKE)**2 + &
                                     PEFIELDV(IIB:IIE,IJB:IJE,IKB:IKE)**2 + &
                                     PEFIELDW(IIB:IIE,IJB:IJE,IKB:IKE)**2)**0.5 
!
DO WHILE (.NOT. GEND_DOMAIN .AND. INB_CELL .LT. NMAX_CELL)  
!
! find the maximum electric field on each proc
  IF (COUNT(GPOSS(IIB:IIE,IJB:IJE,IKB:IKE)) .GT. 0) THEN
    ZMAXE = MAXVAL(ZEMODULE(IIB:IIE,IJB:IJE,IKB:IKE), MASK=GPOSS(IIB:IIE,IJB:IJE,IKB:IKE))
  ELSE
    ZMAXE = 0.
  END IF
!
! find the max electric field on the whole domain + the proc that contains this value
  CALL MAX_ELEC_ll (ZMAXE, IPROC_CELL)
!
  IF (ZMAXE .GT. ZE_TRIG_THRES) THEN
    INB_CELL = INB_CELL + 1  ! one cell is detected
    ZEMAX(INB_CELL) = ZMAXE
! local coordinates of the maximum electric field
    ICELL_LOC(1:3,INB_CELL) = MAXLOC(ZEMODULE, MASK=GPOSS )
    IICOORD = ICELL_LOC(1,INB_CELL)
    IJCOORD = ICELL_LOC(2,INB_CELL)
    ICELL_LOC(1,INB_CELL) = IICOORD + IXOR -1
    ICELL_LOC(2,INB_CELL) = IJCOORD + IYOR -1
    IKCOORD = ICELL_LOC(3,INB_CELL) 
    ICELL_LOC(4,INB_CELL) = IPROC_CELL
! 
! Broadcast the center of the cell to all procs
    CALL MPI_BCAST (ICELL_LOC(:,INB_CELL), 4, MNHINT_MPI, IPROC_CELL, &
                    NMNH_COMM_WORLD, IERR)
!
!
!*      2.2     horizontal extension of the cell 
!
    DO IK = IKB, IKE 
      IF (IPROC_CELL .EQ. IPROC) THEN
        IF (GPOSS(IICOORD,IJCOORD,IK)) THEN
          ZCELL(IICOORD,IJCOORD,IK,INB_CELL) = 1.
          GPOSS(IICOORD,IJCOORD,IK) = .FALSE.
        END IF
      END IF
!
!*      2.2.1   do the neighbour points have q_tot > q_thresh?
!
      GEND_CELL = .FALSE.
      DO WHILE (.NOT. GEND_CELL)
!
        CALL ADD2DFIELD_ll  ( TZFIELDS_ll, ZCELL(:,:,IK,INB_CELL), 'FLASH_GEOM_ELEC_n::ZCELL(:,:,IK,INB_CELL)' )
        CALL UPDATE_HALO_ll ( TZFIELDS_ll, IINFO_ll )
        CALL CLEANLIST_ll   ( TZFIELDS_ll )
!
        COUNT_BEF = COUNT(ZCELL(IIB:IIE,IJB:IJE,IK,INB_CELL) .EQ. 1.)
        CALL SUM_ELEC_ll (COUNT_BEF)
!
        ZCELL_NEW = ZCELL(:,:,IK,INB_CELL)
        DO II = IIB, IIE
          DO IJ = IJB, IJE
            IF ((ZCELL(II,IJ,IK,INB_CELL) .EQ. 0.) .AND.  &
                (GPOSS(II,IJ,IK)) .AND.                   &
                (ZCLOUD(II,IJ,IK) .GT. 1.E-5) .AND.       &
                ((ABS(ZQMT(II,IJ,IK,2)) * PRHODREF(II,IJ,IK) .GT. XQEXCES).OR. &
                 (ABS(ZQMT(II,IJ,IK,3)) * PRHODREF(II,IJ,IK) .GT. XQEXCES).OR. &
                 (ABS(ZQMT(II,IJ,IK,4)) * PRHODREF(II,IJ,IK) .GT. XQEXCES).OR. &
                 (ABS(ZQMT(II,IJ,IK,5)) * PRHODREF(II,IJ,IK) .GT. XQEXCES).OR. &
                 (ABS(ZQMT(II,IJ,IK,6)) * PRHODREF(II,IJ,IK) .GT. XQEXCES)) )THEN
!
              IF ((ZCELL(II-1,IJ,  IK,INB_CELL) .EQ. 1.) .OR. &
                  (ZCELL(II+1,IJ,  IK,INB_CELL) .EQ. 1.) .OR. &
                  (ZCELL(II,  IJ-1,IK,INB_CELL) .EQ. 1.) .OR. &
                  (ZCELL(II,  IJ+1,IK,INB_CELL) .EQ. 1.) .OR. &
                  (ZCELL(II-1,IJ-1,IK,INB_CELL) .EQ. 1.) .OR. &
                  (ZCELL(II-1,IJ+1,IK,INB_CELL) .EQ. 1.) .OR. &
                  (ZCELL(II+1,IJ+1,IK,INB_CELL) .EQ. 1.) .OR. &
                  (ZCELL(II+1,IJ-1,IK,INB_CELL) .EQ. 1.)) THEN
                GPOSS(II,IJ,IK) = .FALSE.
                ZCELL_NEW(II,IJ) = 1.
              END IF
            END IF
          END DO
        END DO
        ZCELL(:,:,IK,INB_CELL) = ZCELL_NEW  
!
        COUNT_AFT = COUNT(ZCELL(IIB:IIE,IJB:IJE,IK,INB_CELL) .EQ. 1.)
        CALL SUM_ELEC_ll(COUNT_AFT)
!
        IF (COUNT_BEF .EQ. COUNT_AFT) THEN
          GEND_CELL = .TRUE.  ! no more point in the cell at this level
        ELSE
          GEND_CELL = .FALSE.
        END IF
      END DO  ! end loop gend_cell
    END DO  ! end loop ik
!
! avoid cell detection in the colums where a previous cell is already present 
    DO II = IIB, IIE
      DO IJ = IJB, IJE
        DO IK = IKB, IKE
          IF (ZCELL(II,IJ,IK,INB_CELL) .EQ. 1.) GPOSS(II,IJ,:) = .FALSE.
        END DO
      END DO
    END DO
  ELSE  
    GEND_DOMAIN = .TRUE.    ! no more points with E > E_threshold
  END IF  ! max E
END DO  ! end loop gend_domain
!
DEALLOCATE (GPOSS)
DEALLOCATE (ZEMAX)
!
!
!*      2.3     if at least 1 cell, allocate arrays
!
IF (INB_CELL .GE. 1) THEN
!
! mean mesh size
  ZMEAN_GRID = (XDXHATM**2 + XDYHATM**2 +                            &
               (SUM(XZHAT(2:SIZE(PRT,3)) - XZHAT(1:SIZE(PRT,3)-1)) / &
               (SIZE(PRT,3)-1.))**2)**0.5
! chaque proc calcule son propre zmean_grid
! mais cette valeur peut etre differente sur chaque proc (ex: relief)
! laisse tel quel pour le moment
!
  ALLOCATE (ISEG_LOC(3*SIZE(PRT,3), INB_CELL)) ! 3 coord indices of the leader
  ALLOCATE (ZCOORD_TRIG(3, INB_CELL))
  ALLOCATE (ZCOORD_SEG(NBRANCH_MAX*3, INB_CELL))
                                         ! NBRANCH_MAX=5000 default
                                         ! 3= 3 coord index
  ALLOCATE (ZCOORD_SEG_ALL(NBRANCH_MAX*3, INB_CELL))
  ALLOCATE (ISEG_GLOB(NBRANCH_MAX*3, INB_CELL))
  ISEG_GLOB(:,:) = 0
!
  IF(LLMA) THEN
    ALLOCATE (ILMA_SEG_ALL (NBRANCH_MAX*3, INB_CELL))
    ALLOCATE (ZLMA_QMT(NBRANCH_MAX*NSV_ELEC, INB_CELL))  ! charge des part.
                                                  ! a neutraliser
    ALLOCATE (ZLMA_PRT(NBRANCH_MAX*NSV_ELEC, INB_CELL))  ! mixing ratio
    ALLOCATE (ZLMA_NEUT_POS(NBRANCH_MAX, INB_CELL))
    ALLOCATE (ZLMA_NEUT_NEG(NBRANCH_MAX, INB_CELL))
    ZLMA_QMT(:,:) = 0.
    ZLMA_PRT(:,:) = 0.
    ZLMA_NEUT_POS(:,:) = 0.
    ZLMA_NEUT_NEG(:,:) = 0.
  END IF
!
  IF (LLNOX_EXPLICIT) THEN
    ALLOCATE (ZLNOX(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
    ZLNOX(:,:,:) = 0.
  END IF
!
  ALLOCATE (ZEM_TRIG(INB_CELL))
  ALLOCATE (INB_FLASH(INB_CELL))
  ALLOCATE (INB_FL_REAL(INB_CELL))
  ALLOCATE (INBSEG(INB_CELL)) 
  ALLOCATE (INBSEG_ALL(INB_CELL))
  ALLOCATE (ITYPE(INB_CELL)) 
  ALLOCATE (INBSEG_LEADER(INB_CELL))
  ALLOCATE (ZDQDT(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),SIZE(PRT,4)+1))
  ALLOCATE (ZSIGMA(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),SIZE(PRT,4)-1))
  ALLOCATE (ZLBDAR(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
  ALLOCATE (ZLBDAS(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
  ALLOCATE (ZLBDAG(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
  IF (KRR == 7) ALLOCATE (ZLBDAH(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
  ALLOCATE (ZSIGLOB(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
  ALLOCATE (ZFLASH(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),INB_CELL))
  ALLOCATE (ZDIST(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
  ALLOCATE (ZQFLASH(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
  ALLOCATE (GATTACH(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
!
  ISEG_LOC(:,:) = 0
  ZCOORD_TRIG(:,:) = 0.
  ZCOORD_SEG(:,:) = 0.
  ZDQDT(:,:,:,:) = 0.
  ZSIGMA(:,:,:,:) = 0.
  ZLBDAR(:,:,:) = 0.
  ZLBDAS(:,:,:) = 0.
  ZLBDAG(:,:,:) = 0. 
  ZSIGLOB(:,:,:) = 0.
  ZFLASH(:,:,:,:) = 0.
  ZDIST(:,:,:) = 0.
  ZQFLASH(:,:,:) = 0.
  ZEM_TRIG(:) = 0.
  INB_FLASH(:) = 0
  INB_FL_REAL(:) = 0
  INBSEG(:) = 0
  INBSEG_ALL(:) = 0 
  INBSEG_PROC(:) = 0 
  INBSEG_LEADER(:) = 0
  ITYPE(:) = 1  ! default = IC
!
!
!-------------------------------------------------------------------------------
!
!*      3.      COMPUTE THE EFFICIENT CROSS SECTIONS OF HYDROMETEORS
!               ----------------------------------------------------
!
!*      3.1     for cloud droplets
!
  WHERE (PRT(:,:,:,2) > ZCLOUDLIM)
    ZSIGMA(:,:,:,1) = XFQLIGHTC * PRHODREF(:,:,:) * PRT(:,:,:,2)
  ENDWHERE
!
!
!*      3.2     for raindrops
!
  WHERE (PRT(:,:,:,3) > 0.0)
    ZLBDAR(:,:,:) = XLBR * (PRHODREF(:,:,:) * &
                            MAX(PRT(:,:,:,3),XRTMIN(3)))**XLBEXR
  END WHERE
!
  WHERE (PRT(:,:,:,3) > ZCLOUDLIM .AND. ZLBDAR(:,:,:) < XLBDAR_MAXE .AND. &
                                        ZLBDAR(:,:,:) > 0.)
    ZSIGMA(:,:,:,2) = XFQLIGHTR * ZLBDAR(:,:,:)**XEXQLIGHTR
  END WHERE
!
!
!*      3.3     for ice crystals
!
  WHERE (PRT(:,:,:,4) > ZCLOUDLIM .AND. PCIT(:,:,:) > 1.E4)
    ZSIGMA(:,:,:,3) = XFQLIGHTI * PCIT(:,:,:)**(1.-XEXQLIGHTI) * &
                     ((PRHODREF(:,:,:) * PRT(:,:,:,4))**XEXQLIGHTI)
  ENDWHERE
!
!
!*      3.4     for snow
!
  WHERE (PRT(:,:,:,5) > 0.0)
    ZLBDAS(:,:,:) = MIN(XLBDAS_MAXE,               &
                        XLBS * (PRHODREF(:,:,:) *  &
                        MAX(PRT(:,:,:,5),XRTMIN(5)))**XLBEXS)
  END WHERE
!
  WHERE (PRT(:,:,:,5) > ZCLOUDLIM .AND. ZLBDAS(:,:,:) < XLBDAS_MAXE .AND. &
                                        ZLBDAS(:,:,:) > 0.)
    ZSIGMA(:,:,:,4) = XFQLIGHTS * ZLBDAS(:,:,:)**XEXQLIGHTS
  ENDWHERE
!
!
!*      3.5     for graupel
!
  WHERE (PRT(:,:,:,6) > 0.0)
    ZLBDAG(:,:,:) = XLBG * (PRHODREF(:,:,:) * MAX(PRT(:,:,:,6),XRTMIN(6)))**XLBEXG
  END WHERE
!
  WHERE (PRT(:,:,:,6) > ZCLOUDLIM .AND. ZLBDAG(:,:,:) < XLBDAG_MAXE .AND. &
                                        ZLBDAG(:,:,:) > 0.)
    ZSIGMA(:,:,:,5) = XFQLIGHTG * ZLBDAG(:,:,:)**XEXQLIGHTG
  ENDWHERE
!
!
!*      3.6     for hail
!
  IF (KRR == 7) THEN
    WHERE (PRT(:,:,:,7) > 0.0)
      ZLBDAH(:,:,:) = XLBH * (PRHODREF(:,:,:) * &
                      MAX(PRT(:,:,:,7), XRTMIN(7)))**XLBEXH
    END WHERE
!
    WHERE (PRT(:,:,:,7) > ZCLOUDLIM .AND. ZLBDAH(:,:,:) < XLBDAH_MAXE .AND. &
                                          ZLBDAH(:,:,:) > 0.)
      ZSIGMA(:,:,:,6) = XFQLIGHTH * ZLBDAH(:,:,:)**XEXQLIGHTH
    ENDWHERE
  END IF
!
!
!*      3.7     sum of the efficient cross sections
!
  ZSIGLOB(:,:,:) = ZSIGMA(:,:,:,1) + ZSIGMA(:,:,:,2) + ZSIGMA(:,:,:,3) + &
                   ZSIGMA(:,:,:,4) + ZSIGMA(:,:,:,5)
!
  IF (KRR == 7) ZSIGLOB(:,:,:) = ZSIGLOB(:,:,:) + ZSIGMA(:,:,:,6)
!
IF (KRR == 7) THEN
   CALL MPPDB_CHECK3DM("flash:: ZLBDAR,ZLBDAS,ZLBDAG,ZLBDAH",PRECISION,&
        ZLBDAR,ZLBDAS,ZLBDAG,ZLBDAH,&
        ZSIGMA(:,:,:,1),ZSIGMA(:,:,:,2),ZSIGMA(:,:,:,3),ZSIGMA(:,:,:,4),&
        ZSIGMA(:,:,:,5),ZSIGMA(:,:,:,6))
ELSE
   CALL MPPDB_CHECK3DM("flash:: ZLBDAR,ZLBDAS,ZLBDAG",PRECISION,&
        ZLBDAR,ZLBDAS,ZLBDAG,&
        ZSIGMA(:,:,:,1),ZSIGMA(:,:,:,2),ZSIGMA(:,:,:,3),ZSIGMA(:,:,:,4),&
        ZSIGMA(:,:,:,5))
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      4.      FIND THE TRIGGERING POINT IN EACH CELL
!               --------------------------------------
!
  ALLOCATE (IPROC_TRIG(INB_CELL))
  ALLOCATE (ISIGNE_EZ(INB_CELL))
  ALLOCATE (GNEW_FLASH(INB_CELL))
  ALLOCATE (ZNEUT_POS(INB_CELL))
  ALLOCATE (ZNEUT_NEG(INB_CELL))
!
  IPROC_TRIG(:) = 0
  ISIGNE_EZ(:) = 0
  GNEW_FLASH(:) = .FALSE.
  ZNEUT_POS(:) = 0.
  ZNEUT_NEG(:) = 0.
!
  CALL TRIG_POINT
!
!
!-------------------------------------------------------------------------------
!
!*      4.      FLASH TRIGGERING
!               ----------------
!
  IFLASH_COUNT = 0
  IFLASH_COUNT_GLOB = 0
!
  DO WHILE (GNEW_FLASH_GLOB)
!
    GATTACH(:,:,:) = .FALSE.
!
    DO IL = 1, INB_CELL
      IF (GNEW_FLASH(IL)) THEN
        ZFLASH(:,:,:,IL) = 0.
! update lightning informations
        INB_FLASH(IL) = INB_FLASH(IL) + 1   ! nb of flashes / cell / time step
        INB_FL_REAL(IL) = INB_FL_REAL(IL) + 1   ! nb of flashes / cell / time step
        INBSEG(IL) = 0        ! nb of segments / flash
        ITYPE(IL) = 1
!
        IF (IPROC .EQ. IPROC_TRIG(IL)) THEN 
           ZEMOD_BL = ZEM_TRIG(IL)
           IIBL_LOC = ISEG_LOC(1,IL)    
           IJBL_LOC = ISEG_LOC(2,IL)
           IKBL     = ISEG_LOC(3,IL) 
!
           INBSEG(IL) = 1        ! nb of segments / flash
           ZFLASH(IIBL_LOC,IJBL_LOC,IKBL,IL) = 1.
        ENDIF
!
        GCG = .FALSE.
        GCG_POS = .FALSE.

        CALL MPPDB_CHECK3DM("flash:: 4. ZFLASH(IL)",PRECISION,&
             ZFLASH(:,:,:,IL))
!
!
!-------------------------------------------------------------------------------
!
!*      5.      PROPAGATE THE BIDIRECTIONAL LEADER
!               ----------------------------------
!
! it is assumed that the leader propagates only along the vertical
!
!*      5.1     positive segments
!
! the positive leader propagates parallel to the electric field
        ISIGN_LEADER = 1
        CALL ONE_LEADER
        IPOS_LEADER = INBSEG(IL) -1
!
!
!*      5.2     negative segments
!
! the negative leader propagates anti-parallel to the electric field
        ZEMOD_BL = ZEM_TRIG(IL)
        IKBL     = ISEG_LOC(3,IL)
        ISIGN_LEADER = -1
        CALL ONE_LEADER
!
        INBSEG_LEADER(IL) = INBSEG(IL)
        INEG_LEADER = INBSEG_LEADER(IL) - IPOS_LEADER - 1
!
! Eliminate this flash if only positive or negative leader exists        
        IF (IPROC .EQ. IPROC_TRIG(IL)) THEN 
          IF (IPOS_LEADER .EQ. 0 .OR. INEG_LEADER .EQ. 0) THEN
            ZFLASH(IIBL_LOC,IJBL_LOC,IKB:IKE,IL) = 0.
            INB_FL_REAL(IL) = INB_FL_REAL(IL) - 1
            GNEW_FLASH(IL) = .FALSE.
          ELSE    ! return to actual Triggering electrical field
            IIBL_LOC = ISEG_LOC(1,IL)
            IJBL_LOC = ISEG_LOC(2,IL)
            IKBL     = ISEG_LOC(3,IL)
            ZEM_TRIG(IL) = ZEM_TRIG(IL)/ZPRES_COEF(IIBL_LOC,IJBL_LOC,IKBL)
          ENDIF
        ENDIF

       CALL MPPDB_CHECK3DM("flash:: 5. ZFLASH(IL)",PRECISION,&
             ZFLASH(:,:,:,IL))
!
        CALL MPI_BCAST (GNEW_FLASH(IL),1,   MNHLOG_MPI, IPROC_TRIG(IL), &
                        NMNH_COMM_WORLD, IERR)
        CALL MPI_BCAST (ZEM_TRIG(IL), 1,    MNHREAL_MPI, IPROC_TRIG(IL), &
                        NMNH_COMM_WORLD, IERR)
        CALL MPI_BCAST (INB_FL_REAL(IL), 1, MNHINT_MPI, IPROC_TRIG(IL), &
                        NMNH_COMM_WORLD, IERR)
      END IF
    END DO  ! end loop il
!
!
!-------------------------------------------------------------------------------
!
!*      6.      POSITIVE AND NEGATIVE REGIONS WHERE THE FLASH CAN PROPAGATE
!               -----------------------------------------------------------
!
! Note: this is done to avoid branching in a third charge region:
! the branches 'stay' in the 2 charge regions where the bileader started to propagate
!
!*      6.1     positive charge region associated to the negative leader
!
    ALLOCATE (GPROP(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),INB_CELL))
    GPROP(:,:,:,:) = .FALSE.
!
    GPOSITIVE = .TRUE.
    CALL CHARGE_POCKET
!
!
!*      6.2     negative charge region associated to the positive leader
!
    GPOSITIVE = .FALSE.
    CALL CHARGE_POCKET
!
! => a point can be added to the flash only if gprop = true
!
!
!-------------------------------------------------------------------------------
!
!*      7.      NUMBER OF POINTS TO REDISTRIBUTE AT DISTANCE D
!               ----------------------------------------------
!
!*      7.1     distance between the triggering point and each point of the mask
!*              global coordinates: only points possibly contributing to branches
!
    INB_NEUT_OK = 0
!
    DO IL = 1, INB_CELL
      IF (GNEW_FLASH(IL)) THEN
        INB_PROP = COUNT(GPROP(IIB:IIE,IJB:IJE,IKB:IKE,IL))
        CALL SUM_ELEC_ll(INB_PROP)
!
        IF (INB_PROP .GT. 0) THEN
          ZDIST(:,:,:) = 0.
          DO II = IIB, IIE
            DO IJ = IJB, IJE
              DO IK = IKB, IKE
                IF (GPROP(II,IJ,IK,IL)) THEN
                  ZDIST(II,IJ,IK) = ((ZXMASS(II) - ZCOORD_TRIG(1,IL))**2 + &
                                     (ZYMASS(IJ) - ZCOORD_TRIG(2,IL))**2 + &
                                     (ZZMASS(II,IJ,IK) - ZCOORD_TRIG(3,IL))**2)**0.5
                END IF
              END DO
            END DO
          END DO
!
!
!*      7.3     compute the min and max distance from the triggering point - global
!
          ZMIN_DIST = 0.0
          ZMAX_DIST = MAX_ll(ZDIST,IPROC_AUX)
!
! transform the min and max distances into min and max increments 
          IIND_MIN = 1
          IIND_MAX = MAX(1, INT((ZMAX_DIST-ZMIN_DIST)/ZMEAN_GRID +1.))
          IDELTA_IND = IIND_MAX + 1
!
          ALLOCATE (IHIST_LOC(IDELTA_IND))
          ALLOCATE (ZHIST_PERCENT(IDELTA_IND))
          ALLOCATE (IHIST_GLOB(IDELTA_IND))
          ALLOCATE (ZMAX_BRANCH(IDELTA_IND))
          ALLOCATE (IMAX_BRANCH(IDELTA_IND))
          ALLOCATE (IMASKQ_DIST(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)))
!
          IHIST_LOC(:) = 0
          ZHIST_PERCENT(:) = 0.
          IHIST_GLOB(:) = 0
          ZMAX_BRANCH(:) = 0.
          IMAX_BRANCH(:) = 0
          IMASKQ_DIST(:,:,:) = 0
!
!
!*      7.4     histogram: number of points between r and r+dr
!*              for each proc
!
! build an array with the possible points: IMASKQ_DIST contains the distance 
! rank of points contributing to branches, excluding the leader points
!
          DO II = IIB, IIE 
            DO IJ = IJB, IJE
              DO IK = IKB, IKE
                IF (ZDIST(II,IJ,IK) .NE. 0.) THEN
                  IM = INT( (ZDIST(II,IJ,IK)-ZMIN_DIST)/ZMEAN_GRID + 1.)
                  IHIST_LOC(IM) = IHIST_LOC(IM) + 1
                  IMASKQ_DIST(II,IJ,IK) = IM
                ENDIF
              END DO
            END DO
          END DO
!
!
!*      7.5     global histogram
!
          IHIST_GLOB(:) = IHIST_LOC(:) 
          CALL SUM_ELEC_ll(IHIST_GLOB)
!
!
!*      7.6     normalization
!
          ZHIST_PERCENT(:) = 0.
          ZMAX_BRANCH(:) = 0.
          IMAX_BRANCH(:) = 0
!
          DO IM = 1, IDELTA_IND
            IF (IHIST_GLOB(IM) .GT. 0) THEN
              ZHIST_PERCENT(IM) = REAL(IHIST_LOC(IM)) / REAL(IHIST_GLOB(IM))
            END IF
!
!
!-------------------------------------------------------------------------------
!
!*      8.      BRANCHES
!               --------
!
!*      8.1     max number of branches at distance d from the triggering point
!
            ZMAX_BRANCH(IM) = (XDFRAC_L / ZMEAN_GRID) * &
                              REAL(IIND_MIN+IM-1)**(XDFRAC_ECLAIR - 1.)
            ZMAX_BRANCH(IM) = ANINT(ZMAX_BRANCH(IM))
! all procs know the max total number of branches at distance d
! => the max number of branches / proc is proportional to the percentage of 
! available points / proc at this distance
!
            IMAX_BRANCH(IM) = INT(ANINT(ZMAX_BRANCH(IM)))
          END DO
!
          DEALLOCATE (IHIST_LOC)
          DEALLOCATE (ZHIST_PERCENT)
          DEALLOCATE (IHIST_GLOB)
          DEALLOCATE (ZMAX_BRANCH)
!
!
!*      8.3     distribute the branches
!
!
          CALL BRANCH_GEOM(IKB, IKE)
!
          DEALLOCATE (IMAX_BRANCH)
          DEALLOCATE (IMASKQ_DIST)
        END IF   ! end if count(gprop)
!
!
!-------------------------------------------------------------------------------
!
!*      9.      NEUTRALIZATION
!               --------------
        CALL MPPDB_CHECK3DM("flash:: 9. ZQMTOT",PRECISION,ZQMTOT)
        CALL MPPDB_CHECK3DM("flash:: 9. ZFLASH",PRECISION,ZFLASH(:,:,:,IL))
!
!*      9.1     charge carried by the lightning flash
!
        ZQFLASH(:,:,:) = 0.
        WHERE (ZFLASH(IIB:IIE,IJB:IJE,IKB:IKE,IL) .GT. 0. .AND.          &
               ABS(ZQMTOT(IIB:IIE,IJB:IJE,IKB:IKE) *                     &
                   PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE)) .GT. XQNEUT .AND. &
               ZSIGLOB(IIB:IIE,IJB:IJE,IKB:IKE) .GE. ZSIGMIN)
          ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) = -1. *               &             
                         (ABS(ZQMTOT(IIB:IIE,IJB:IJE,IKB:IKE)) / &
                              ZQMTOT(IIB:IIE,IJB:IJE,IKB:IKE)) * &
                         (ABS(ZQMTOT(IIB:IIE,IJB:IJE,IKB:IKE)) - &
                         (XQNEUT / PRHODREF(IIB:IIE,IJB:IJE,IKB:IKE)))
          GATTACH(IIB:IIE,IJB:IJE,IKB:IKE) = .TRUE.

        END WHERE
!
! net charge carried by the flash (for charge conservation / IC)
        ZQNET = SUM3D_ll(ZQFLASH*PRHODJ, IINFO_ll)
!
!
!*      9.2     number of points to neutralize
!
        INB_NEUT = COUNT(ZSIGLOB(IIB:IIE,IJB:IJE,IKB:IKE) .GE. ZSIGMIN .AND. &
                         ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) .NE. 0.)
        CALL SUM_ELEC_ll(INB_NEUT)

!
!
!*      9.3     ensure total charge conservation for IC
!
        IF (INB_NEUT .GE. 3) THEN
          GNEUTRALIZATION = .TRUE.
        ELSE
          GNEUTRALIZATION = .FALSE.
          GNEW_FLASH(IL) = .FALSE. 
          INB_FL_REAL(IL) = INB_FL_REAL(IL) - 1
        END IF
!
        IF (GNEUTRALIZATION .AND. (.NOT. GCG) .AND. ZQNET .NE. 0.) THEN
          ZQNET = ZQNET / REAL(INB_NEUT) 
          WHERE (ZSIGLOB(IIB:IIE,IJB:IJE,IKB:IKE) .GE. ZSIGMIN .AND. &
                 ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) .NE. 0.)
            ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) = ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) - &
                                       ZQNET / PRHODJ(IIB:IIE,IJB:IJE,IKB:IKE) 
          ENDWHERE
        END IF
!
!
!*      9.4     charge neutralization 
!
        CALL MPPDB_CHECK3DM("flash:: 9.4 ZQFLASH,ZSIGLOB",PRECISION,&
             ZQFLASH,ZSIGLOB)

        ZDQDT(:,:,:,:) = 0.
!  
        IF (GNEUTRALIZATION) THEN
          IF (ITYPE(IL) .EQ. 1.) THEN         
            WHERE (ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) < 0.)
                       !  increase negative ion charge
              ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_ELEC) =         &
                     ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_ELEC) +  &
                     ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE)
            ENDWHERE
!
            WHERE (ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) > 0.)
                     ! Increase positive ion charge
              ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,1) =         &
                     ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,1) +  &
                     ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE)
            ENDWHERE
!
!
!*      9.4.2   cloud-to-ground flashes
!
          ELSE   
!
! Neutralization of the charge on positive CG flashes
            IF (ITYPE(IL) .EQ. 3) THEN   
              DO II = 1, NSV_ELEC
                WHERE (ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) > 0.)
                  ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,II) =    &
                     ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,II) - &
                     ZQMT(IIB:IIE,IJB:IJE,IKB:IKE,II)
                END WHERE
              ENDDO
!
              WHERE (ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) > 0.) 
                ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE)=0.
              END WHERE
!
              WHERE (ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) < 0.)
! Increase negative ion charge
                ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_ELEC) =         &
                       ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,NSV_ELEC) +  &
                       ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE)
              ENDWHERE
            ELSE
!
! Neutralization of the charge on negative CG flashes
!
              DO II = 1, NSV_ELEC
                WHERE (ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) < 0.)
                  ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,II) =      &
                       ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,II) - &
                       ZQMT(IIB:IIE,IJB:IJE,IKB:IKE,II)
                END WHERE
              ENDDO
!
              WHERE (ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) < 0.)
                ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE)=0.
              END WHERE
!
              WHERE (ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE) > 0.)
                        ! Increase positive ion charge
                ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,1) =         &
                       ZDQDT(IIB:IIE,IJB:IJE,IKB:IKE,1) +  &
                       ZQFLASH(IIB:IIE,IJB:IJE,IKB:IKE)
              ENDWHERE
            END IF        ! GCG_POS
          END IF          ! NOT(GCG)
!
! Counting the total number of points neutralized in the cell
          IF (IPROC .EQ. IPROC_TRIG(IL)) THEN
             INB_NEUT_OK = INB_NEUT_OK + INB_NEUT
          END IF
!
          CALL MPI_BCAST (INB_NEUT_OK,1, MNHINT_MPI, IPROC_TRIG(IL), &
                    NMNH_COMM_WORLD, IERR)
!
!*      9.5     Gather lightning information from all processes
!*              Save the particule charge and total pos/neg charge neutralization points.
!*                   the coordinates of all flash branch points
!
          CALL MPI_ALLGATHER(INBSEG(IL), 1, MNHINT_MPI, &
                             INBSEG_PROC,  1, MNHINT_MPI, NMNH_COMM_WORLD, IERR)

          INBSEG_ALL(IL) = INBSEG(IL)
          CALL SUM_ELEC_ll(INBSEG_ALL(IL))

          CALL GATHER_ALL_BRANCH
!
!*      9.6     update the source term
!
          CALL MPPDB_CHECK3DM("flash:: 9.6 PRSVS",PRECISION,&
               PRSVS(:,:,:,1),PRSVS(:,:,:,2),PRSVS(:,:,:,3),PRSVS(:,:,:,4),&
               PRSVS(:,:,:,5),PRSVS(:,:,:,6),PRSVS(:,:,:,7))
          CALL MPPDB_CHECK3DM("flash:: 9.6 ZDQDT",PRECISION,&
               ZDQDT(:,:,:,1),ZDQDT(:,:,:,2),ZDQDT(:,:,:,3),ZDQDT(:,:,:,4),&
               ZDQDT(:,:,:,5),ZDQDT(:,:,:,6),ZDQDT(:,:,:,7))

          DO II = IIB, IIE
            DO IJ = IJB, IJE
              DO IK = IKB, IKE
                DO IM = 1, NSV_ELEC
                  IF (ZDQDT(II,IJ,IK,IM) .NE. 0.) THEN
                    PRSVS(II,IJ,IK,IM) = PRSVS(II,IJ,IK,IM) + &
                                         ZDQDT(II,IJ,IK,IM) * &
                                         PRHODJ(II,IJ,IK) / PTSTEP
                  END IF
!
!
!*      9.7     update the positive and negative charge neutralized
!
                  IF (ZDQDT(II,IJ,IK,IM) .LT. 0.) THEN
                    ZNEUT_NEG(IL) = ZNEUT_NEG(IL) + ZDQDT(II,IJ,IK,IM) * &
                                                    PRHODJ(II,IJ,IK) 
                  ELSE IF (ZDQDT(II,IJ,IK,IM) .GT. 0.) THEN
                    ZNEUT_POS(IL) = ZNEUT_POS(IL) + ZDQDT(II,IJ,IK,IM) * &
                                                    PRHODJ(II,IJ,IK) 
                  END IF
                END DO
              END DO
            END DO
          END DO
!
          CALL SUM_ELEC_ll(ZNEUT_POS(IL))
          CALL SUM_ELEC_ll(ZNEUT_NEG(IL))
!
!
!*      9.8     compute the NOx production
!
!!   The lightning length is first computed. The number of NOx molecules per
!! meter of lightning flash is taken from Wang et al. (1998). It is a linear
!! function of the pressure. No distinction is made between ICs and CGs.

          IF (LLNOX_EXPLICIT) THEN
            IFLASH_COUNT_GLOB = IFLASH_COUNT_GLOB + 1
            IF (INBSEG(IL) .NE. 0) THEN
              DO II = 0, INBSEG(IL)-1
                IM = 3 * II
                IX = ISEG_GLOB(IM+1,IL) - IXOR + 1
                IY = ISEG_GLOB(IM+2,IL) - IYOR + 1
                IZ = ISEG_GLOB(IM+3,IL)
                ZLGHTLENGTH = (XDXX(IX,IY,IZ) * XDYY(IX,IY,IZ) * &
                               XDZZ(IX,IY,IZ))**(1./3.)
                ZLNOX(IX, IY, IZ) = ZLNOX(IX, IY, IZ) +                     &
                                   (XWANG_A + XWANG_B * PPABST(IX,IY,IZ)) * &
                                    ZLGHTLENGTH
              ENDDO
              IFLASH_COUNT = IFLASH_COUNT + 1
            END IF
          END IF
        END IF            !  GNEUTRALIZATION
      END IF    ! end if gnew_flash
    END DO    ! end loop il
!
    DEALLOCATE (GPROP)
!
!
!----------------------------------------------------------------------------
!
!*      10.     PRINT OR SAVE (before print) LIGHTNING INFORMATIONS
!               ---------------------------------------------------
!
! Synchronizing all processes
!   CALL MPI_BARRIER(NMNH_COMM_WORLD, IERR)   ! A ACTIVER SI PB.
!
    INBLIGHT = COUNT(GNEW_FLASH(1:INB_CELL))
    IF (IPROC .EQ. 0) THEN
      IF (INBLIGHT .NE. 0) THEN
        IF ((NNBLIGHT+INBLIGHT) .LE. NFLASH_WRITE) THEN       ! SAVE
          ISAVE_STATUS = 1
          DO IL = 1, INB_CELL
            IF (GNEW_FLASH(IL)) THEN
              NNBLIGHT = NNBLIGHT + 1
              ISFLASH_NUMBER(NNBLIGHT) = ISFLASH_NUMBER(NNBLIGHT-1) + 1
              ISNB_FLASH(NNBLIGHT) = INB_FL_REAL(IL)
              ISNBSEG(NNBLIGHT) = INBSEG_ALL(IL)
              ISCELL_NUMBER(NNBLIGHT) = IL
              ISTCOUNT_NUMBER(NNBLIGHT) = KTCOUNT
              ISTYPE(NNBLIGHT) = ITYPE(IL)
              ZSEM_TRIG(NNBLIGHT) = ZEM_TRIG(IL) / 1000.
              ZSNEUT_POS(NNBLIGHT) = ZNEUT_POS(IL) 
              ZSNEUT_NEG(NNBLIGHT) = ZNEUT_NEG(IL)
!
              DO II = 1, INBSEG_ALL(IL)
                IM = 3 * (II - 1)
                ZSCOORD_SEG(NNBLIGHT,II,1:3) = ZCOORD_SEG_ALL(IM+1:IM+3,IL)
              ENDDO
!
              IF(LLMA) THEN
                DO II = 1, INBSEG_ALL(IL)
                  IM = 3 * (II - 1)
                  ISLMA_SEG_GLOB(NNBLIGHT,II,1:3) = ILMA_SEG_ALL(IM+1:IM+3,IL)
                  IM = NSV_ELEC * (II - 1)
                  ZSLMA_QMT(NNBLIGHT,II,2:6) = ZLMA_QMT(IM+2:IM+6,IL)
                  ZSLMA_PRT(NNBLIGHT,II,2:6) = ZLMA_PRT(IM+2:IM+6,IL)
                  ZSLMA_NEUT_POS(NNBLIGHT,II) = ZLMA_NEUT_POS(II,IL)
                  ZSLMA_NEUT_NEG(NNBLIGHT,II) = ZLMA_NEUT_NEG(II,IL)
                END DO
              END IF   ! llma
            END IF   ! gnew_flash
          END DO   ! end loop il
!
          IF (NNBLIGHT .EQ. NFLASH_WRITE) ISAVE_STATUS = 0
!
        ELSE    ! Print in output files
          ISAVE_STATUS = 2
        END IF
!      
        IF (ISAVE_STATUS .EQ. 0 .OR. ISAVE_STATUS .EQ. 2) THEN
          CALL WRITE_OUT_ASCII
          IF(LLMA) THEN
            CALL WRITE_OUT_LMA
          END IF
          ISFLASH_NUMBER(0) = ISFLASH_NUMBER(NNBLIGHT)
        END IF
!
        IF (ISAVE_STATUS .EQ. 2) THEN   ! Save flashes of the temporal loop
          NNBLIGHT = 0
          DO IL = 1, INB_CELL
            IF (GNEW_FLASH(IL)) THEN
              NNBLIGHT = NNBLIGHT + 1
              ISFLASH_NUMBER(NNBLIGHT) = ISFLASH_NUMBER(NNBLIGHT-1) + 1
              ISNB_FLASH(NNBLIGHT) = INB_FL_REAL(IL)
              ISNBSEG(NNBLIGHT) = INBSEG_ALL(IL)
              ISCELL_NUMBER(NNBLIGHT) = IL
              ISTCOUNT_NUMBER(NNBLIGHT) = KTCOUNT
              ISTYPE(NNBLIGHT) = ITYPE(IL)
              ZSEM_TRIG(NNBLIGHT) = ZEM_TRIG(IL) / 1000.
              ZSNEUT_POS(NNBLIGHT) = ZNEUT_POS(IL) 
              ZSNEUT_NEG(NNBLIGHT) = ZNEUT_NEG(IL)
!
              DO II = 1, INBSEG_ALL(IL)
                IM = 3 * (II - 1)
                ZSCOORD_SEG(NNBLIGHT, II, 1:3) = ZCOORD_SEG_ALL(IM+1:IM+3, IL)
              ENDDO
!
              IF(LLMA) THEN
                DO II = 1, INBSEG_ALL(IL)
                  IM = 3 * (II - 1)
                  ISLMA_SEG_GLOB(NNBLIGHT,II,1:3) = ILMA_SEG_ALL(IM+1:IM+3,IL)
                  IM = NSV_ELEC*(II-1)
                  ZSLMA_QMT(NNBLIGHT,II,2:6) = ZLMA_QMT(IM+2:IM+6,IL)
                  ZSLMA_PRT(NNBLIGHT,II,2:6) = ZLMA_PRT(IM+2:IM+6,IL)
                  ZSLMA_NEUT_POS(NNBLIGHT,II) = ZLMA_NEUT_POS(II,IL)
                  ZSLMA_NEUT_NEG(NNBLIGHT,II) = ZLMA_NEUT_NEG(II,IL)
                END DO
              END IF
            END IF
          ENDDO
        END IF
!
        IF (ISAVE_STATUS .EQ. 0) THEN
          NNBLIGHT = 0
        END IF
      END IF   ! INBLIGHT
    END IF   ! IPROC
!
! Save flash location statistics in all processes
    IF (INBLIGHT .NE. 0) THEN
      DO IL = 1, INB_CELL
        IF (GNEW_FLASH(IL)) THEN
          IMAP2D(:,:)   = 0
          DO IK = IKB, IKE
            IMAP2D(:,:) = IMAP2D(:,:) + ZFLASH(:,:,IK,IL)
          END DO
!
! Detect Trig/Impact X,Y location
          IX = 0
          IY = 0
          GFIRSTFLASH = .FALSE.
          DO II = IIB, IIE
            DO IJ = IJB, IJE
              DO IK = IKB, IKE
                IF (GFIRSTFLASH) EXIT
                IF (ZFLASH(II,IJ,IK,IL)==1.) THEN
                  IX = II
                  IY = IJ
                  GFIRSTFLASH = .TRUE.
                END IF
              END DO
            END DO
          END DO
!
! Store
          IF (ITYPE(IL)==1) THEN ! IC
            IF (IX*IY/=0) NMAP_TRIG_IC(IX,IY) = NMAP_TRIG_IC(IX,IY) + 1
            NMAP_2DAREA_IC(:,:) = NMAP_2DAREA_IC(:,:) + MIN(1,IMAP2D(:,:))
            NMAP_3DIC(:,:,:)    = NMAP_3DIC(:,:,:) + ZFLASH(:,:,:,IL)
          ELSE ! CGN & CGP
            IF (IX*IY/=0) NMAP_IMPACT_CG(IX,IY) = NMAP_IMPACT_CG(IX,IY) + 1
            NMAP_2DAREA_CG(:,:) = NMAP_2DAREA_CG(:,:) + MIN(1,IMAP2D(:,:))
            NMAP_3DCG(:,:,:)    = NMAP_3DCG(:,:,:) + ZFLASH(:,:,:,IL)
          END IF
        END IF
      ENDDO
    END IF   ! INBLIGHT
!
!------------------------------------------------------------------------------
!
!*    11.      ATTACHMENT AFTER CHARGE NEUTRALIZATION
!              --------------------------------------
!
!*    11.1     ion attachment
!
    IF (INB_NEUT_OK .NE. 0) THEN

       CALL MPPDB_CHECK3DM("flash:: PRSVS",PRECISION,&
            PRSVS(:,:,:,1),PRSVS(:,:,:,2),PRSVS(:,:,:,3),PRSVS(:,:,:,4),&
            PRSVS(:,:,:,5),PRSVS(:,:,:,6),PRSVS(:,:,:,7))

      PRSVS(:,:,:,1) = PRSVS(:,:,:,1) / XECHARGE
      PRSVS(:,:,:,NSV_ELEC) = - PRSVS(:,:,:,NSV_ELEC) / XECHARGE
!
      CALL ION_ATTACH_ELEC(KTCOUNT, KRR, PTSTEP, PRHODREF,                   &
                           PRHODJ, PRSVS, PRS, PTHT, PCIT, PPABST, PEFIELDU, &
                           PEFIELDV, PEFIELDW, GATTACH, PTOWN, PSEA          )
!
      PRSVS(:,:,:,1) = PRSVS(:,:,:,1) * XECHARGE
      PRSVS(:,:,:,NSV_ELEC) = - PRSVS(:,:,:,NSV_ELEC) * XECHARGE

       CALL MPPDB_CHECK3DM("flash:: after ION PRSVS",PRECISION,&
            PRSVS(:,:,:,1),PRSVS(:,:,:,2),PRSVS(:,:,:,3),PRSVS(:,:,:,4),&
            PRSVS(:,:,:,5),PRSVS(:,:,:,6),PRSVS(:,:,:,7))
    ENDIF
!
!
!*    11.2    update the charge density to check if another flash can be triggered
!
    ZQMTOT(:,:,:) = 0.  
    DO II = 1, NSV_ELEC
! transform the source term (C/s) into the updated charge density (C/kg)
      ZQMT(:,:,:,II) = PRSVS(:,:,:,II) * PTSTEP / PRHODJ(:,:,:)
!
! total charge density (C/kg)
      ZQMTOT(:,:,:)  = ZQMTOT(:,:,:) + PRSVS(:,:,:,II) * PTSTEP / PRHODJ(:,:,:)
    END DO
!
!
!-------------------------------------------------------------------------------
!
!*      12.     CHECK IF ANOTHER FLASH CAN BE TRIGGERED
!               ---------------------------------------
!

    IF ((MAXVAL(INB_FLASH(:))+1) < INBFTS_MAX) THEN
      IF (INB_NEUT_OK .NE. 0) THEN
         CALL MPPDB_CHECK3DM("flash:: PRHODJ,PRT",PRECISION,&
              PRHODJ,PRT(:,:,:,1),PRT(:,:,:,2),PRT(:,:,:,3),PRT(:,:,:,4),&
              PRT(:,:,:,5),PRT(:,:,:,6))
         CALL MPPDB_CHECK3DM("flash:: ZQMT",PRECISION,&
              ZQMT(:,:,:,1),ZQMT(:,:,:,2),ZQMT(:,:,:,3),ZQMT(:,:,:,4),&
              ZQMT(:,:,:,5),ZQMT(:,:,:,6),ZQMT(:,:,:,7))
        CALL TO_ELEC_FIELD_n (PRT, ZQMT, PRHODJ, KTCOUNT, KRR, &
                              PEFIELDU, PEFIELDV, PEFIELDW)
        CALL MPPDB_CHECK3DM("flash:: PEFIELDU, PEFIELDV, PEFIELDW",PRECISION,&
             PEFIELDU, PEFIELDV, PEFIELDW)
! electric field module including pressure effect
        ZEMODULE(IIB:IIE,IJB:IJE,IKB:IKE) = ZPRES_COEF(IIB:IIE,IJB:IJE,IKB:IKE)*    &
                                            (PEFIELDU(IIB:IIE,IJB:IJE,IKB:IKE)**2 + &
                                             PEFIELDV(IIB:IIE,IJB:IJE,IKB:IKE)**2 + &
                                             PEFIELDW(IIB:IIE,IJB:IJE,IKB:IKE)**2)**0.5
      ENDIF
!
      ISEG_LOC(:,:) = 0
      ZCOORD_TRIG(:,:) = 0.
      ZCOORD_SEG(:,:) = 0.
      IPROC_TRIG(:) = 0
      ISIGNE_EZ(:) = 0
!
      CALL TRIG_POINT
    ELSE
      GNEW_FLASH_GLOB = .FALSE.
    END IF
!
    ZNEUT_POS(:) = 0.
    ZNEUT_NEG(:) = 0.
!
    IF (LLMA) THEN
      ZLMA_NEUT_POS(:,:) = 0.
      ZLMA_NEUT_NEG(:,:) = 0.
    END IF
  END DO   ! end loop do while
!
!
!-------------------------------------------------------------------------------
!
!*      13.     COMPUTE THE NOX SOURCE TERM
!               ---------------------------
!
  IF (LLNOX_EXPLICIT) THEN
    IF (IFLASH_COUNT_GLOB .NE. 0) THEN
      ZCOEF = XMD / XAVOGADRO
      XLNOX_ECLAIR = 0.
      IF (IFLASH_COUNT .NE. 0) THEN
        XLNOX_ECLAIR = SUM(ZLNOX(:,:,:))
        PSVS_LINOX(:,:,:) = PSVS_LINOX(:,:,:) + ZLNOX(:,:,:) * ZCOEF ! PRHODJ is
                                                                     ! implicit
      END IF
      CALL SUM_ELEC_ll (XLNOX_ECLAIR)
      XLNOX_ECLAIR = XLNOX_ECLAIR / (XAVOGADRO * REAL(IFLASH_COUNT_GLOB))
    END IF
    DEALLOCATE (ZLNOX)
  END IF
!
  DEALLOCATE (ZNEUT_POS)
  DEALLOCATE (ZNEUT_NEG)
  DEALLOCATE (ZSIGMA)
  DEALLOCATE (ZLBDAR)
  DEALLOCATE (ZLBDAS)
  DEALLOCATE (ZLBDAG)
  IF (KRR == 7) DEALLOCATE (ZLBDAH)
  DEALLOCATE (ZSIGLOB)
  DEALLOCATE (ZDQDT)
  DEALLOCATE (ZDIST)
  DEALLOCATE (ZFLASH)
  DEALLOCATE (ZQFLASH)
  DEALLOCATE (IPROC_TRIG)
  DEALLOCATE (ISIGNE_EZ)
  DEALLOCATE (GNEW_FLASH)
  DEALLOCATE (INBSEG)
  DEALLOCATE (INBSEG_ALL)
  DEALLOCATE (INBSEG_LEADER)
  DEALLOCATE (INB_FLASH)
  DEALLOCATE (INB_FL_REAL)
  DEALLOCATE (ZEM_TRIG)
  DEALLOCATE (ITYPE)
  DEALLOCATE (ISEG_LOC)
  DEALLOCATE (ZCOORD_TRIG)
  DEALLOCATE (ZCOORD_SEG)
  DEALLOCATE (ZCOORD_SEG_ALL)
  DEALLOCATE (ISEG_GLOB)
  DEALLOCATE (GATTACH)
  IF(LLMA) THEN
    DEALLOCATE (ILMA_SEG_ALL)
    DEALLOCATE (ZLMA_QMT)
    DEALLOCATE (ZLMA_PRT)
    DEALLOCATE (ZLMA_NEUT_POS)
    DEALLOCATE (ZLMA_NEUT_NEG)
  END IF
END IF   ! (inb_cell .ge. 1)
!
!
!-------------------------------------------------------------------------------
!
!*      13.     PRINT LIGHTNING INFORMATIONS FOR THE LAST TIMESTEP
!               OR LMA_TIME_SAVE IS REACHED IF LLMA OPTION IS USED
!               --------------------------------------------------
!
IF (LLMA) THEN
  IF( IPROC .EQ. 0 .AND. TDTCUR%xtime >= TDTLMA%xtime - PTSTEP ) THEN
    CALL WRITE_OUT_ASCII
    CALL WRITE_OUT_LMA
    ISFLASH_NUMBER(0) = ISFLASH_NUMBER(NNBLIGHT)
    NNBLIGHT = 0
  END IF
END IF
!
IF (NNBLIGHT .NE. 0 .AND. ((IPROC .EQ. 0 .AND. OEXIT) .OR. &
                           (KTCOUNT == NSTOP .AND. KMI==1))) THEN
  CALL WRITE_OUT_ASCII
  IF(LLMA) CALL WRITE_OUT_LMA
END IF
!
!
!-------------------------------------------------------------------------------
!
!*      14.     DEALLOCATE
!               ----------
!
DEALLOCATE (ICELL_LOC)
DEALLOCATE (ZQMT)
DEALLOCATE (ZQMTOT)
DEALLOCATE (ZCLOUD)
DEALLOCATE (ZCELL)
DEALLOCATE (ZEMODULE)
!
!
!-------------------------------------------------------------------------------
!
!*      14.     BACK TO INPUT UNITS (per kg and per (m3 s)) FOR IONS
!               ----------------------------------------------------
!
PRSVS(:,:,:,1) = PRSVS(:,:,:,1) / XECHARGE          ! 1 /(m3 s)
PRSVS(:,:,:,NSV_ELEC) = -PRSVS(:,:,:,NSV_ELEC) / XECHARGE    ! 1 /(m3 s)
!
!
!-------------------------------------------------------------------------------
!
CONTAINS
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE TRIG_POINT ()
!
! Goal : find randomly a triggering point where E > E_trig
!
!*      0.      DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
!*      0.1     declaration of dummy arguments
!
!*      0.2     declaration of local variables
!
LOGICAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3),INB_CELL) :: &
           GTRIG  ! mask for the triggering pts
INTEGER :: INB_TRIG  ! Nb of pts where triggering is possible
INTEGER :: IWEST_GLOB_TRIG   ! western  global limit of possible triggering
INTEGER :: IEAST_GLOB_TRIG   ! eastern  global limit of possible triggering
INTEGER :: ISOUTH_GLOB_TRIG  ! southern global limit of possible triggering
INTEGER :: INORTH_GLOB_TRIG  ! northern global limit of possible triggering
INTEGER :: IUP_TRIG          ! upper limit of possible triggering
INTEGER :: IDOWN_TRIG        ! down limit of possible triggering
!
!
!*      1.      INITIALIZATIONS 
!               -----------
!
GTRIG(:,:,:,:) = .FALSE.
GNEW_FLASH(:) = .FALSE.
GNEW_FLASH_GLOB = .FALSE.
!
!
!*      2.      FIND THE POSSIBLE TRIGGERING POINTS 
!               -----------------------------------
!
DO IL = 1, INB_CELL
  WHERE (ZEMODULE(IIB:IIE,IJB:IJE,IKB:IKE) > ZE_TRIG_THRES .AND. &
         ZCELL(IIB:IIE,IJB:IJE,IKB:IKE,IL) .GT. 0.)
    GTRIG(IIB:IIE,IJB:IJE,IKB:IKE,IL) = .TRUE.
  ENDWHERE
END DO
!
!
!*      3.     CHOICE OF THE TRIGGERING POINT 
!              ------------------------------
!
!*      3.1    number and coordinates of the possible triggering points
!
INB_TRIG = 0
DO IL = 1, INB_CELL
  INB_TRIG = COUNT(GTRIG(IIB:IIE,IJB:IJE,IKB:IKE,IL))
  CALL SUM_ELEC_ll(INB_TRIG)
!
!
!*      3.2    random choice of the triggering point
!
 IF (INB_TRIG .GT. 0) THEN
    IFOUND = 0
!
! find the global limits where GTRIG = T
    CALL EXTREMA_ELEC_ll(GTRIG(:,:,:,IL), IWEST_GLOB_TRIG,  IEAST_GLOB_TRIG,  &
                                          ISOUTH_GLOB_TRIG, INORTH_GLOB_TRIG, &
                                          IDOWN_TRIG, IUP_TRIG)
!
    DO WHILE (IFOUND .NE. 1)
!
! random choice of the 3 global ind.
      CALL MNH_RANDOM_NUMBER(ZRANDOM)
      II_TRIG_GLOB = IWEST_GLOB_TRIG + &
                      INT(ANINT(ZRANDOM * (IEAST_GLOB_TRIG - IWEST_GLOB_TRIG)))
      CALL MNH_RANDOM_NUMBER(ZRANDOM)
      IJ_TRIG_GLOB = ISOUTH_GLOB_TRIG + &
                      INT(ANINT(ZRANDOM * (INORTH_GLOB_TRIG - ISOUTH_GLOB_TRIG)))
      CALL MNH_RANDOM_NUMBER(ZRANDOM)
      IK_TRIG = IDOWN_TRIG + INT(ANINT(ZRANDOM * (IUP_TRIG - IDOWN_TRIG)))
!
! global ind. --> local ind. of the potential triggering pt
      II_TRIG_LOC = II_TRIG_GLOB - IXOR + 1
      IJ_TRIG_LOC = IJ_TRIG_GLOB - IYOR + 1
!
! test if the randomly chosen pt meets all conditions for triggering
      IF ((II_TRIG_LOC .LE. IIE) .AND. (II_TRIG_LOC .GE. IIB) .AND. &
          (IJ_TRIG_LOC .LE. IJE) .AND. (IJ_TRIG_LOC .GE. IJB) .AND. &
          (IK_TRIG     .LE. IKE) .AND. (IK_TRIG     .GE. IKB)) THEN
        IF (GTRIG(II_TRIG_LOC,IJ_TRIG_LOC,IK_TRIG,IL)) THEN
          IFOUND = 1
!
! update the local coordinates of the flash segments
          ISEG_LOC(1,IL) = II_TRIG_LOC
          ISEG_LOC(2,IL) = IJ_TRIG_LOC
          ISEG_LOC(3,IL) = IK_TRIG
!
          ISEG_GLOB(1,IL) = II_TRIG_GLOB
          ISEG_GLOB(2,IL) = IJ_TRIG_GLOB
          ISEG_GLOB(3,IL) = IK_TRIG
!
          ZCOORD_TRIG(1,IL) = ZXMASS(II_TRIG_LOC)
          ZCOORD_TRIG(2,IL) = ZYMASS(IJ_TRIG_LOC)
          ZCOORD_TRIG(3,IL) = ZZMASS(II_TRIG_LOC, IJ_TRIG_LOC, IK_TRIG)
!
          ZCOORD_SEG(1:3,IL) = ZCOORD_TRIG(1:3,IL)
!
! electric field module at the triggering point
          ZEM_TRIG(IL) = ZEMODULE(II_TRIG_LOC,IJ_TRIG_LOC,IK_TRIG)
!
! sign of Ez at the triggering point
          ISIGNE_EZ(IL) = 0
          IF (PEFIELDW(II_TRIG_LOC,IJ_TRIG_LOC,IK_TRIG) .GT. 0.) THEN
            ISIGNE_EZ(IL) = 1
          ELSE IF (PEFIELDW(II_TRIG_LOC,IJ_TRIG_LOC,IK_TRIG) .LT. 0.) THEN
            ISIGNE_EZ(IL) = -1
          END IF
        END IF
      END IF
!
! broadcast IFOUND and find the proc where IFOUND = 1
      CALL MAX_ELEC_ll (IFOUND, IPROC_TRIG(IL))
!
    END DO
!
!
!
!*      4.      BROADCAST USEFULL PARAMETERS
!               ----------------------------
!
    CALL MPI_BCAST (ZEM_TRIG(IL), 1, &
                    MNHREAL_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)
    CALL MPI_BCAST (ISEG_LOC(:,IL), 3*SIZE(PRT,3), &
                    MNHINT_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)
    CALL MPI_BCAST (ZCOORD_TRIG(:,IL), 3, &
                    MNHREAL_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)
    CALL MPI_BCAST (ISIGNE_EZ(IL), 1, &
                    MNHINT_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)
!
!
!*      5.      CHECK IF THE FLASH CAN DEVELOP
!               ------------------------------
!
    IF (INB_FLASH(IL) < INBFTS_MAX) THEN
      IF (IPROC.EQ.IPROC_TRIG(IL)) THEN  
        ZCELL(II_TRIG_LOC,IJ_TRIG_LOC,IK_TRIG,IL) = 0.
      END IF
!
      GNEW_FLASH(IL) = .TRUE.
      GNEW_FLASH_GLOB = .TRUE.
      CALL MPI_BCAST (GNEW_FLASH(IL),1, MNHLOG_MPI, IPROC_TRIG(IL), &
                      NMNH_COMM_WORLD, IERR)
      CALL MPI_BCAST (GNEW_FLASH_GLOB,1, MNHLOG_MPI, IPROC_TRIG(IL), &
                      NMNH_COMM_WORLD, IERR)
    END IF
  END IF 
END DO
!
!
END SUBROUTINE TRIG_POINT
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE ONE_LEADER ()
!
!! Purpose: propagates the bidirectional leader along the vertical
!
!*      0.      DECLARATIONS
!               ------------
!
IMPLICIT NONE
!
INTEGER :: IKSTEP, IIDECAL         
!
!*      1.      BUILD THE POSITIVE/NEGATIVE LEADER
!               ----------------------------------
CALL MPPDB_CHECK3DM("flash:: one_leader ZFLASH",PRECISION,ZFLASH(:,:,:,IL))
!
IKSTEP = ISIGN_LEADER * ISIGNE_EZ(IL)
    ! the positive leader propagates parallel to the electric field
    ! while the negative leader propagates anti// to the electric field
ISTOP = 0
!
!
IF (IPROC .EQ. IPROC_TRIG(IL)) THEN

  DO WHILE (ZEMOD_BL > XEPROP .AND. IKBL > IKB .AND. &
            IKBL < IKE .AND. ISTOP .EQ. 0 .AND.      &
            INBSEG(IL) .LE. (NLEADER_MAX-1))  
!
! local coordinates of the new segment
    IIBL_LOC = ISEG_LOC(1,IL)
    IJBL_LOC = ISEG_LOC(2,IL)
    IKBL     = IKBL + IKSTEP
    IIDECAL = INBSEG(IL) * 3
!
    ISEG_LOC(IIDECAL+1,IL) = IIBL_LOC
    ISEG_LOC(IIDECAL+2,IL) = IJBL_LOC
    ISEG_LOC(IIDECAL+3,IL) = IKBL
!
    ISEG_GLOB(IIDECAL+1,IL) = IIBL_LOC + IXOR - 1
    ISEG_GLOB(IIDECAL+2,IL) = IJBL_LOC + IYOR - 1
    ISEG_GLOB(IIDECAL+3,IL) = IKBL
!
    ZCOORD_SEG(IIDECAL+1,IL) = ZXMASS(IIBL_LOC)
    ZCOORD_SEG(IIDECAL+2,IL) = ZYMASS(IJBL_LOC)
    ZCOORD_SEG(IIDECAL+3,IL) = ZZMASS(IIBL_LOC, IJBL_LOC, IKBL)
!
    INBSEG(IL) = INBSEG(IL) + 1
!
!
!*      1.3     test if Ez keeps the same sign
!
    IF (PEFIELDW(IIBL_LOC,IJBL_LOC,IKBL) .EQ. 0. .OR. &
        INT(ABS(PEFIELDW(IIBL_LOC,IJBL_LOC,IKBL)) / &
                PEFIELDW(IIBL_LOC,IJBL_LOC,IKBL)) /= ISIGNE_EZ(IL) .OR. &
                ZCELL(IIBL_LOC,IJBL_LOC,IKBL,IL) .EQ. 0.) THEN
      ISTOP = 1
! then this segment is not part of the leader
      INBSEG(IL) = INBSEG(IL) - 1
    END IF
!
!
!*      1.4     sign of the induced charge
!
    IF (ISTOP .EQ. 0) THEN
      ZFLASH(IIBL_LOC,IJBL_LOC,IKBL,IL) = 1.
      ZCELL(IIBL_LOC,IJBL_LOC,IKBL,IL) = 0.   
!
!
!*      1.6     electric field module at the tip of the leader
!
      ZEMOD_BL = ZEMODULE(IIBL_LOC,IJBL_LOC,IKBL)
!
!
!*      1.7     test if the domain boundaries are reached
!
      IF ((IIBL_LOC < IIB .AND. LWEST_ll())  .OR. &
          (IIBL_LOC > IIE .AND. LEAST_ll())  .OR. &
          (IJBL_LOC < IJB .AND. LSOUTH_ll()) .OR. &
          (IJBL_LOC > IJE .AND. LNORTH_ll())) THEN
        PRINT*,'DOMAIN BOUNDARIES REACHED BY THE LIGHTNING ' 
        ISTOP = 1
      ENDIF
!
      IF (IKBL .LE. IKB) THEN
        PRINT*,'THE LIGHTNING FLASH HAS REACHED THE GROUND ' 
        ISTOP = 1
        GCG = .TRUE.
        NNB_CG = NNB_CG + 1
        IF (ISIGN_LEADER > 0) THEN
          GCG_POS = .TRUE.
          ITYPE(IL) = 3 ! CGP
          NNB_CG_POS = NNB_CG_POS + 1 
        ELSE
          ITYPE(IL) = 2 ! CGN
        END IF
      ENDIF
!
      IF (IKBL .GE. IKE) THEN
        PRINT*,'THE LIGHTNING FLASH HAS REACHED THE TOP OF THE DOMAIN ' 
        ISTOP = 1
      ENDIF
!
!
!*      2.      TEST IF THE FLASH IS A CG
!               -------------------------
!
      IF (.NOT. GCG) THEN
        IF ( (ZZMASS(IIBL_LOC,IJBL_LOC,IKBL)-PZZ(IIBL_LOC,IJBL_LOC,IKB)) <=   &
             XALT_CG .AND. INBSEG(IL) .GT. 1  .AND. IKSTEP .LT. 0) THEN
!
!
!*      2.1    the channel is prolongated to the ground if 
!*             one segment reaches the altitude XALT_CG
!
          DO WHILE (IKBL > IKB)
            IKBL = IKBL - 1
!
! local coordinates of the new segment
            IIDECAL = INBSEG(IL) * 3
!
            ISEG_LOC(IIDECAL+1,IL) = IIBL_LOC
            ISEG_LOC(IIDECAL+2,IL) = IJBL_LOC
            ISEG_LOC(IIDECAL+3,IL) = IKBL
!
            ISEG_GLOB(IIDECAL+1:IIDECAL+2,IL) = ISEG_GLOB(IIDECAL-2:IIDECAL-1,IL)
            ISEG_GLOB(IIDECAL+3,IL) = IKBL
!
            ZCOORD_SEG(IIDECAL+1:IIDECAL+2,IL) = ZCOORD_SEG(IIDECAL-2:IIDECAL-1,IL)
            ZCOORD_SEG(IIDECAL+3,IL) = ZZMASS(IIBL_LOC, IJBL_LOC, IKBL)
!
!  Increment number of segments
            INBSEG(IL) = INBSEG(IL) + 1 ! Nb of segments
            ZFLASH(IIBL_LOC,IJBL_LOC,IKBL,IL) = 1.
            ZCELL(IIBL_LOC,IJBL_LOC,IKBL,IL) = 0.   
          END DO
!
!
!*      2.2    update the number of CG flashes
!
          GCG = .TRUE.
          NNB_CG = NNB_CG + 1 
          ISTOP = 1
!
          IF (ISIGN_LEADER > 0) THEN
            GCG_POS = .TRUE.
            NNB_CG_POS = NNB_CG_POS + 1 
            ITYPE(IL) = 3
          ELSE
            ITYPE(IL) = 2
          END IF
        END IF
      END IF
    END IF     ! end if ISTOP=0
  END DO   ! end loop leader
END IF  ! only iproc_trig was working
!
!
!*      3.     BROADCAST THE INFORMATIONS TO ALL PROCS
!              ---------------------------------------
!
CALL MPI_BCAST (ISEG_LOC(:,IL), 3*SIZE(PRT,3), &  
                MNHINT_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)
CALL MPI_BCAST (ITYPE(IL), 1, &
                MNHINT_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)

CALL MPI_BCAST (GCG, 1, &
                MNHLOG_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)
CALL MPI_BCAST (GCG_POS, 1, &
                MNHLOG_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)
CALL MPI_BCAST (NNB_CG, 1, &
                MNHINT_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)
CALL MPI_BCAST (NNB_CG_POS, 1, &
                MNHINT_MPI, IPROC_TRIG(IL), NMNH_COMM_WORLD, IERR)

!
CALL MPPDB_CHECK3DM("flash:: one_leader end ZFLASH",PRECISION,ZFLASH(:,:,:,IL))
!
END SUBROUTINE ONE_LEADER
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE CHARGE_POCKET
!
!!
!! Purpose: limit flash propagation into the positive and negative charge layers
!!          located immediatly above and below the triggering point
!!
!*       0.      DECLARATIONS
!                ------------
!
IMPLICIT NONE
!
REAL, DIMENSION(SIZE(PRT,1),SIZE(PRT,2),SIZE(PRT,3)) :: ZSIGN_AREA,ZSIGN_AREA_NEW

REAL, DIMENSION(INB_CELL) :: ZSIGN  ! sign of the charge immediatly below/above the triggering pt 
!
INTEGER, DIMENSION(INB_CELL) :: IEND  ! if 1, no more neighbour pts meeting the conditions
INTEGER, DIMENSION(INB_CELL) :: COUNT_BEF2
INTEGER, DIMENSION(INB_CELL) :: COUNT_AFT2
INTEGER :: IPROC_END
INTEGER :: IEND_GLOB
INTEGER :: IIDECAL, IKMIN, IKMAX
REAL :: ZFACT
!
!
!*       1.      SEARCH THE POINTS BELONGING TO THE LAYERS 
!                -----------------------------------------
!
ZFACT = -1.
IF(GPOSITIVE) ZFACT = 1.

ZSIGN_AREA(:,:,:) = 0.
ZSIGN(:) = 0.
IEND(:) = 0
IEND_GLOB = 0
!
!
DO IL = 1, INB_CELL
  IF (.NOT. GNEW_FLASH(IL)) THEN
    IEND(IL) = 1
    IEND_GLOB = IEND_GLOB + IEND(IL)
  END IF
  IF (GNEW_FLASH(IL) .AND. IPROC .EQ. IPROC_TRIG(IL)) THEN
    DO II = 1, INBSEG(IL)
      IIDECAL = 3 * (II - 1)
      IIBL_LOC = ISEG_LOC(IIDECAL+1,IL)
      IJBL_LOC = ISEG_LOC(IIDECAL+2,IL)
      IKBL     = ISEG_LOC(IIDECAL+3,IL)
!
      IF (ZQMTOT(IIBL_LOC,IJBL_LOC,IKBL) .GT. 0. .AND. GPOSITIVE) THEN
        ZSIGN_AREA(IIBL_LOC,IJBL_LOC,IKBL) = 1. * REAL(IL)
        ZSIGN(IL) = ZSIGN_AREA(IIBL_LOC,IJBL_LOC,IKBL)
      ELSE IF (ZQMTOT(IIBL_LOC,IJBL_LOC,IKBL) .LT. 0. .AND. .NOT.GPOSITIVE) THEN
        ZSIGN_AREA(IIBL_LOC,IJBL_LOC,IKBL) = -1. * REAL(IL)
        ZSIGN(IL) = ZSIGN_AREA(IIBL_LOC,IJBL_LOC,IKBL)
      END IF
    END DO
  END IF
!
  CALL MPI_BCAST (ZSIGN(IL), 1, MNHREAL_MPI, IPROC_TRIG(IL), &
                  NMNH_COMM_WORLD, IERR)
END DO
!
DO WHILE (IEND_GLOB .NE. INB_CELL)
  DO IL = 1, INB_CELL
    CALL ADD3DFIELD_ll  ( TZFIELDS_ll, ZSIGN_AREA, 'FLASH_GEOM_ELEC_n::ZSIGN_AREA' )
    CALL UPDATE_HALO_ll ( TZFIELDS_ll, IINFO_ll)
    CALL CLEANLIST_ll   ( TZFIELDS_ll)
!
    IF (GNEW_FLASH(IL) .AND. (IEND(IL) .NE. 1)) THEN
      COUNT_BEF2(IL) = COUNT(ZSIGN_AREA(IIB:IIE,IJB:IJE,IKB:IKE) .EQ. ZSIGN(IL))
      CALL SUM_ELEC_ll (COUNT_BEF2(IL))
!
      IF (ISIGNE_EZ(IL).EQ.1) THEN
        IF (GPOSITIVE) THEN
          IKMIN = IKB
          IKMAX = ISEG_LOC(3, IL)
        ELSE
          IKMIN = ISEG_LOC(3, IL)
          IKMAX = IKE
        ENDIF
      ENDIF
!
      IF (ISIGNE_EZ(IL).EQ.-1) THEN
        IF (GPOSITIVE) THEN
          IKMIN = ISEG_LOC(3, IL)
          IKMAX = IKE
        ELSE
          IKMIN = IKB
          IKMAX = ISEG_LOC(3, IL)
        ENDIF
      ENDIF
!
      ZSIGN_AREA_NEW(:,:,IKMIN:IKMAX) = ZSIGN_AREA (:,:,IKMIN:IKMAX)
      DO II = IIB, IIE
        DO IJ = IJB, IJE
          DO IK = IKMIN, IKMAX
            IF ((ZSIGN_AREA(II,  IJ,  IK)   .EQ. 0.) .AND.   &
                (ZCELL(II,IJ,IK,IL) .EQ. 1.) .AND.           &
                (.NOT. GPROP(II,IJ,IK,IL)) .AND.             &
                (ZQMTOT(II,IJ,IK)*ZFACT .GT. 0.) .AND.       &
                (ABS(ZQMTOT(II,IJ,IK) *                      &
                     PRHODREF(II,IJ,IK)) .GT. XQNEUT)) THEN
!
              IF ((ZSIGN_AREA(II-1,IJ,  IK)   .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ,  IK)   .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II,  IJ-1,IK)   .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II,  IJ+1,IK)   .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II-1,IJ-1,IK)   .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II-1,IJ+1,IK)   .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ+1,IK)   .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ-1,IK)   .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II,  IJ,  IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II-1,IJ,  IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ,  IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II,  IJ-1,IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II,  IJ+1,IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II-1,IJ-1,IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II-1,IJ+1,IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ+1,IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ-1,IK+1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II,  IJ,  IK-1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II-1,IJ,  IK-1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ,  IK-1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II,  IJ-1,IK-1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II,  IJ+1,IK-1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II-1,IJ-1,IK-1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II-1,IJ+1,IK-1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ+1,IK-1) .EQ. ZSIGN(IL)) .OR. &
                  (ZSIGN_AREA(II+1,IJ-1,IK-1) .EQ. ZSIGN(IL))) THEN
                ZSIGN_AREA_NEW(II,IJ,IK) = ZSIGN(IL)
                GPROP(II,IJ,IK,IL) = .TRUE.
              END IF
            END IF
          END DO
        END DO
      END DO
      ZSIGN_AREA (:,:,IKMIN:IKMAX) = ZSIGN_AREA_NEW(:,:,IKMIN:IKMAX)
!
      COUNT_AFT2(IL) = COUNT(ZSIGN_AREA(IIB:IIE,IJB:IJE,IKB:IKE) .EQ. ZSIGN(IL))
      CALL SUM_ELEC_ll(COUNT_AFT2(IL))
!
      IF (COUNT_BEF2(IL) .EQ. COUNT_AFT2(IL)) THEN
        IEND(IL) = 1
      ELSE
        IEND(IL) = 0
      END IF
! broadcast IEND and find the proc where IEND = 1
      CALL MAX_ELEC_ll (IEND(IL), IPROC_END)
      IEND_GLOB = IEND_GLOB + IEND(IL)
    END IF
  END DO
END DO  ! end do while
!
END SUBROUTINE CHARGE_POCKET
!
!-------------------------------------------------------------------------------
!
  SUBROUTINE BRANCH_GEOM (IKMIN, IKMAX)
!
! Goal : find randomly flash branch points
!
!*      0.      DECLARATIONS
!               ------------
!
use modd_precision, only: MNHINT64, MNHINT64_MPI

IMPLICIT NONE
!
!*      0.1     declaration of dummy arguments
!
INTEGER, INTENT(IN) :: IKMIN, IKMAX
!
!*      0.2     declaration of local variables
!
INTEGER :: IIDECALB
INTEGER :: IPLOOP  ! loop index for the proc number
INTEGER :: IMIN, IMAX
INTEGER :: IAUX
INTEGER :: INB_SEG_BEF  ! nb of segments before branching
INTEGER :: INB_SEG_AFT  ! nb of segments after branching
INTEGER :: INB_SEG_TO_BRANCH ! = NBRANCH_MAX-INB_SEG_BEF
LOGICAL :: GRANDOM           ! T = the gridpoints are chosen randomly
INTEGER, DIMENSION(NPROC) :: INBPT_PROC
REAL, DIMENSION(:), ALLOCATABLE :: ZAUX
!
INTEGER                                           :: JI,JJ,JK,JIL , ICHOICE,IPOINT
INTEGER,                DIMENSION(NPROC+1)        :: IDISPL
INTEGER(kind=MNHINT64), DIMENSION(:), ALLOCATABLE :: I8VECT , I8VECT_LL
INTEGER,                DIMENSION(:), ALLOCATABLE :: IRANK  , IRANK_LL , IORDER_LL
!
!
!
!*      1.      ON EACH PROC, COUNT THE NUMBER OF POINTS AT DISTANCE D 
!*              THAT CAN RECEIVE A BRANCH
!               ------------------------------------------------------
CALL MPPDB_CHECK3DM("flash:: branch ZFLASH,IMASKQ_DIST",PRECISION,&
     ZFLASH(:,:,:,IL),IMASKQ_DIST*1.0)
!
IM = 1
ISTOP = 0
INB_SEG_BEF = COUNT(ZFLASH(IIB:IIE,IJB:IJE,IKB:IKE,IL) .NE. 0.)
CALL SUM_ELEC_ll(INB_SEG_BEF)
!
INB_SEG_TO_BRANCH = NBRANCH_MAX - INB_SEG_BEF
!
DO WHILE (IM .LE. IDELTA_IND .AND. ISTOP .NE. 1)
! number of points that can receive a branch in each proc
  IPT_DIST = COUNT(IMASKQ_DIST(IIB:IIE,IJB:IJE,IKB:IKE) .EQ. IM)
! global number of points that can receive a branch
  IPT_DIST_GLOB = IPT_DIST
  CALL SUM_ELEC_ll (IPT_DIST_GLOB)
!
  IF (IPT_DIST_GLOB .LE. INB_SEG_TO_BRANCH) THEN
    IF (IPT_DIST_GLOB .LE. IMAX_BRANCH(IM)) THEN
      GRANDOM = .FALSE.
    ELSE
      GRANDOM = .TRUE.
    END IF
  ELSE
    GRANDOM = .TRUE.
  END IF
!
!
!*      2.      DISTRIBUTE THE BRANCHES
!               -----------------------
!
  IF (IPT_DIST_GLOB .GT. 0 .AND. INB_SEG_TO_BRANCH .NE. 0) THEN
    IF (.NOT. GRANDOM) THEN
      INB_SEG_TO_BRANCH = INB_SEG_TO_BRANCH - IPT_DIST_GLOB
!
!*      2.1     all points are selected
!
      IF(IPT_DIST .GT. 0) THEN 
        WHERE (IMASKQ_DIST(IIB:IIE,IJB:IJE,IKB:IKE) .EQ. IM)
          ZFLASH(IIB:IIE,IJB:IJE,IKB:IKE,IL) = 2.
          ZCELL(IIB:IIE,IJB:IJE,IKB:IKE,IL) = 0.
        END WHERE
      END IF
    ELSE
!
!*      2.2      the gridpoints are chosen randomly
!
      IF (IMAX_BRANCH(IM) .GT. 0) THEN
        INBPT_PROC(:) = 0
        CALL MPI_ALLGATHER(IPT_DIST, 1, MNHINT_MPI, &
                   INBPT_PROC, 1, MNHINT_MPI, NMNH_COMM_WORLD, IERR)
!
        IDISPL(1) = 0
        DO JI=2, NPROC+1
           IDISPL(JI) = IDISPL(JI-1)+INBPT_PROC(JI-1)
        ENDDO
!
        ALLOCATE (I8VECT(IPT_DIST))
        ALLOCATE (IRANK(IPT_DIST))
        IF (IPT_DIST .GT. 0) THEN 
           JIL=0
           DO JK=IKB,IKE
              DO JJ=IJB,IJE
                 DO JI=IIB,IIE
                    IF (IMASKQ_DIST(JI,JJ,JK) .EQ. IM) THEN
                       JIL = JIL + 1
                       I8VECT(JIL) = IJU_ll*IIU_ll*(JK-1) + IIU_ll*(JJ-1 +IYOR-1) + (JI +IXOR-1)
                       !print*,"IN  => I8VECT(JIL    )=",I8VECT(JIL),JI,JJ,JK,JIL
                    END IF
                 END DO
              END DO
           END DO
           !
           IRANK(:)  = IPROC           
        END IF
!
        ALLOCATE(I8VECT_LL(IPT_DIST_GLOB))
        ALLOCATE(IRANK_LL(IPT_DIST_GLOB))
        ALLOCATE(IORDER_LL(IPT_DIST_GLOB))
        CALL MPI_ALLGATHERV(I8VECT,IPT_DIST, MNHINT64_MPI,I8VECT_LL , &
                        INBPT_PROC, IDISPL, MNHINT64_MPI, NMNH_COMM_WORLD, IERR)
        CALL MPI_ALLGATHERV(IRANK,IPT_DIST, MNHINT_MPI,IRANK_LL , &
                        INBPT_PROC, IDISPL, MNHINT_MPI, NMNH_COMM_WORLD, IERR)
        CALL N8QUICK_SORT(I8VECT_LL, IORDER_LL)
!
        DO IPOINT = 1, MIN(IMAX_BRANCH(IM), INB_SEG_TO_BRANCH)
           IFOUND = 0
           DO WHILE (IFOUND .NE. 1)
              ! randomly chose points in zvect
              CALL MNH_RANDOM_NUMBER(ZRANDOM)
              ICHOICE = INT(ANINT(ZRANDOM * IPT_DIST_GLOB)) 
              IF (ICHOICE .EQ. 0) ICHOICE = 1
              IF (I8VECT_LL(ICHOICE) .NE. 0 ) THEN
                 IFOUND = 1               
                 ! The points is in this processors  , get is coord and set it
                 IF (IRANK_LL(IORDER_LL(ICHOICE)) .EQ. IPROC) THEN
                    JK = 1 +     (I8VECT_LL(ICHOICE)-1) / ( IJU_ll*IIU_ll ) 
                    JJ = 1 + (   (I8VECT_LL(ICHOICE)-1) - IJU_ll*IIU_ll*(JK-1) ) / IIU_ll  - IYOR +1
                    JI = 1 + MOD((I8VECT_LL(ICHOICE)-1)                          , int(IIU_ll,kind(I8VECT_LL(1)))) - IXOR +1
                    !print*,"OUT => I8VECT_LL(ICHOICE)=",I8VECT_ll(ICHOICE),JI,JJ,JK,ICHOICE
                    ZFLASH(JI,JJ,JK,IL) = 2.
                 END IF
                 I8VECT_LL(ICHOICE) = 0
              ENDIF
           END DO
        END DO
!
        INB_SEG_TO_BRANCH = INB_SEG_TO_BRANCH - MIN(IMAX_BRANCH(IM), INB_SEG_TO_BRANCH)
!
        DEALLOCATE(I8VECT,I8VECT_LL,IRANK,IRANK_LL,IORDER_LL)
        CALL MPPDB_CHECK3DM("flash:: branch IPT_DIST ZFLASH",PRECISION,&
             ZFLASH(:,:,:,IL))
      END IF
    END IF                   !IPT_DIST .LE. IMAX_BRANCH(IM)
  ELSE
! if no pt available at r, then no branching possible at r+dr !
      ISTOP = 1
  END IF  ! end if ipt_dist > 0
!
! next distance
  CALL MPPDB_CHECK3DM("flash:: branch IM+1 ZFLASH",PRECISION,ZFLASH(:,:,:,IL))
  IM = IM + 1
END DO   ! end loop / do while / radius IM
!
INB_SEG_AFT = COUNT (ZFLASH(IIB:IIE,IJB:IJE,IKB:IKE,IL) .NE. 0.)
CALL SUM_ELEC_ll(INB_SEG_AFT)
!
IF (INB_SEG_AFT .GT. INB_SEG_BEF) THEN
  DO II = IIB, IIE
    DO IJ = IJB, IJE
      DO IK = IKB, IKE
        IF (ZFLASH(II,IJ,IK,IL) .EQ. 2.) THEN 
          IIDECALB = INBSEG(IL) * 3
!
          ISEG_GLOB(IIDECALB+1,IL) = II + IXOR - 1
          ISEG_GLOB(IIDECALB+2,IL) = IJ + IYOR - 1
          ISEG_GLOB(IIDECALB+3,IL) = IK
!
          ZCOORD_SEG(IIDECALB+1,IL) = ZXMASS(II)
          ZCOORD_SEG(IIDECALB+2,IL) = ZYMASS(IJ)
          ZCOORD_SEG(IIDECALB+3,IL) = ZZMASS(II,IJ,IK)
          INBSEG(IL) = INBSEG(IL) + 1
        END IF
      END DO
    END DO
  END DO
END IF
!
CALL MPPDB_CHECK3DM("flash:: end branch ZFLASH",PRECISION,ZFLASH(:,:,:,IL))
!
END SUBROUTINE BRANCH_GEOM
!
!--------------------------------------------------------------------------------
!
  SUBROUTINE GATHER_ALL_BRANCH
!
!!
!! Purpose:
!!
!
!*      0.     DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
INTEGER :: INSEGPROC, INSEGCELL  ! number of segments in the process,
                                 ! and number of segments in the cell
INTEGER :: ISAVEDECAL
INTEGER :: INSEGTRIG, IPROCTRIG
REAL, DIMENSION(:), ALLOCATABLE :: ZLMAQMT, ZLMAPRT, ZLMAPOS, ZLMANEG
REAL, DIMENSION(:), ALLOCATABLE :: ZSEND, ZRECV
INTEGER, DIMENSION(:), ALLOCATABLE :: ISEND, IRECV
INTEGER, DIMENSION(NPROC) :: IDECAL, IDECAL3, IDECALN
INTEGER, DIMENSION(NPROC) :: INBSEG_PROC_X3, INBSEG_PROC_XNSV
!
!
IPROCTRIG = IPROC_TRIG(IL)
INSEGCELL = INBSEG_ALL(IL)
INSEGPROC = INBSEG_PROC(IPROC+1)
INSEGTRIG = INBSEG_PROC(IPROCTRIG+1)
!
IDECAL(1) = INSEGTRIG
DO IK = 2, NPROC
  IDECAL(IK) = IDECAL(IK-1) + INBSEG_PROC(IK-1)
END DO
!
IF(IPROCTRIG .EQ. 0) ISAVEDECAL = IDECAL(IPROCTRIG+1)
!
IDECAL(IPROCTRIG+1) = 0
DO IK = IPROCTRIG+2, NPROC
  IF(IPROCTRIG .EQ. 0) THEN
    IDECAL(IK) = IDECAL(IK) - ISAVEDECAL
  ELSE
    IDECAL(IK) = IDECAL(IK) - IDECAL(1)
  END IF
END DO
!
IDECAL3(:) = 3 * IDECAL(:)
!
!
!*      1.     BRANCH COORDINATES
!
ALLOCATE (ZRECV(INSEGCELL*3))
ALLOCATE (ZSEND(INSEGPROC*3))
!
IF (INSEGPROC .NE. 0) THEN
  ZSEND(1:3*INSEGPROC) = ZCOORD_SEG(1:3*INSEGPROC,IL)
END IF
!
IF (IPROC .EQ. 0) THEN
  INBSEG_PROC_X3(:) = 3 * INBSEG_PROC(:)
END IF
!
CALL MPI_GATHERV (ZSEND, 3*INSEGPROC, MNHREAL_MPI, ZRECV, INBSEG_PROC_X3, &
                  IDECAL3, MNHREAL_MPI, 0, NMNH_COMM_WORLD, IERR)
!
IF (IPROC .EQ. 0) THEN
  ZCOORD_SEG_ALL(1:3*INSEGCELL,IL) = ZRECV(1:3*INSEGCELL)
END IF
!
DEALLOCATE (ZRECV)
DEALLOCATE (ZSEND)
!
!
!*      2.     FOR LMA-LIKE RESULTS: Charge, mixing ratio,
!*                        neutralized positive/negative charge
!*                        and grid index
!
IF (LLMA) THEN
  ALLOCATE (ISEND(3*INSEGPROC))
  ALLOCATE (ZLMAQMT(INSEGPROC*NSV_ELEC))
  ALLOCATE (ZLMAPRT(INSEGPROC*NSV_ELEC))
  ALLOCATE (ZLMAPOS(INSEGPROC))
  ALLOCATE (ZLMANEG(INSEGPROC))
!
  ISEND  (:) = 0
  ZLMAPOS(:) = 0.
  ZLMANEG(:) = 0.
  ZLMAQMT(:) = 0.
  ZLMAPRT(:) = 0.
!
  IF (INSEGPROC .NE. 0) THEN
    DO II = 1, INSEGPROC
      IM = 3 * (II - 1)
      IX = ISEG_GLOB(IM+1,IL) - IXOR + 1
      IY = ISEG_GLOB(IM+2,IL) - IYOR + 1
      IZ = ISEG_GLOB(IM+3,IL)
!
      IM = NSV_ELEC * (II - 1)
      IF (IX .LE. IIE .AND. IX .GE. IIB .AND. &
          IY .LE. IJE .AND. IY .GE. IJB) THEN
        ZLMAQMT(IM+2:IM+6) = ZQMT(IX,IY,IZ,2:6)
        ZLMAPRT(IM+2:IM+6) =  PRT(IX,IY,IZ,2:6)
        DO IJ = 1, NSV_ELEC
          IF (ZDQDT(IX,IY,IZ,IJ) .GT. 0.) THEN
            ZLMAPOS(II) = ZLMAPOS(II) + &
                          ZDQDT(IX,IY,IZ,IJ) * PRHODJ(IX,IY,IZ)
          ELSE IF (ZDQDT(IX,IY,IZ,IJ) .LT. 0.) THEN
            ZLMANEG(II) = ZLMANEG(II) + &
                          ZDQDT(IX,IY,IZ,IJ) * PRHODJ(IX,IY,IZ)
          END IF
        END DO
      END IF
    END DO
!
    ISEND(1:3*INSEGPROC) = ISEG_GLOB(1:3*INSEGPROC, IL)
  END IF
!
! Grid Indexes
!
  ALLOCATE (IRECV(3*INSEGCELL))
!
  CALL MPI_GATHERV (ISEND, 3*INSEGPROC, MNHINT_MPI, IRECV, INBSEG_PROC_X3, &
                    IDECAL3, MNHINT_MPI, 0, NMNH_COMM_WORLD, IERR)
!
  IF (IPROC .EQ. 0) THEN
    ILMA_SEG_ALL(1:3*INSEGCELL,IL) = IRECV(1:3*INSEGCELL)
  END IF
!
  DEALLOCATE (IRECV)
  DEALLOCATE (ISEND)
!
! Neutralized charge at grid points
!
  ALLOCATE (ZRECV(INSEGCELL))
!
  CALL MPI_GATHERV (ZLMAPOS, INSEGPROC, MNHREAL_MPI, ZRECV, INBSEG_PROC,  &
                    IDECAL, MNHREAL_MPI, 0, NMNH_COMM_WORLD, IERR)
!
  IF (IPROC .EQ. 0) THEN
    ZLMA_NEUT_POS(1:INSEGCELL,IL) = ZRECV(1:INSEGCELL)
  END IF
!
  CALL MPI_GATHERV (ZLMANEG, INSEGPROC, MNHREAL_MPI, ZRECV, INBSEG_PROC,  &
                    IDECAL, MNHREAL_MPI, 0, NMNH_COMM_WORLD, IERR)
!
  IF (IPROC .EQ. 0) THEN
    ZLMA_NEUT_NEG(1:INSEGCELL,IL) = ZRECV(1:INSEGCELL)
  END IF
!
  DEALLOCATE (ZLMAPOS)
  DEALLOCATE (ZLMANEG)
  DEALLOCATE (ZRECV)
!
! Charge and mixing ratios at neutralized points
!
  ALLOCATE (ZRECV(NSV_ELEC*INSEGCELL))
!
  IDECALN(:) = IDECAL(:) * NSV_ELEC
!
  IF (IPROC .EQ. 0) THEN
    INBSEG_PROC_XNSV(:) = NSV_ELEC * INBSEG_PROC(:)
  END IF
!
  CALL MPI_GATHERV (ZLMAQMT, NSV_ELEC*INSEGPROC, MNHREAL_MPI, ZRECV, &
                    INBSEG_PROC_XNSV,                                &
                    IDECALN, MNHREAL_MPI, 0, NMNH_COMM_WORLD, IERR   )
!
  IF (IPROC .EQ. 0) THEN
    ZLMA_QMT(1:NSV_ELEC*INSEGCELL,IL) = ZRECV(1:NSV_ELEC*INSEGCELL)
  END IF
!
  CALL MPI_GATHERV (ZLMAPRT, NSV_ELEC*INSEGPROC, MNHREAL_MPI, ZRECV, &
                    INBSEG_PROC_XNSV,                                &
                    IDECALN, MNHREAL_MPI, 0, NMNH_COMM_WORLD, IERR   )
!
  IF (IPROC .EQ. 0) THEN
    ZLMA_PRT(1:NSV_ELEC*INSEGCELL,IL) = ZRECV(1:NSV_ELEC*INSEGCELL)
  END IF
!
  DEALLOCATE (ZLMAQMT)
  DEALLOCATE (ZLMAPRT)
  DEALLOCATE (ZRECV)
!
END IF
!
END SUBROUTINE GATHER_ALL_BRANCH
!
!--------------------------------------------------------------------------------
!
  SUBROUTINE PT_DISCHARGE
!
!!
!! Purpose:
!!
!
!*      0.     DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
!
WHERE (ABS(PEFIELDW(:,:,IKB)) > XECORONA .AND. PEFIELDW(:,:,IKB) > 0.)
  PRSVS(:,:,IKB,1) = PRSVS(:,:,IKB,1) +                                       &
                     XFCORONA * PEFIELDW(:,:,IKB) * (ABS(PEFIELDW(:,:,IKB)) - &
                     XECORONA)**2 / (PZZ(:,:,IKB+1) - PZZ(:,:,IKB))
ENDWHERE
!
WHERE (ABS(PEFIELDW(:,:,IKB)) > XECORONA .AND. PEFIELDW(:,:,IKB) < 0.)
  PRSVS(:,:,IKB,NSV_ELEC) = PRSVS(:,:,IKB,NSV_ELEC) +                         &
                     XFCORONA * PEFIELDW(:,:,IKB) * (ABS(PEFIELDW(:,:,IKB)) - &
                     XECORONA)**2 / (PZZ(:,:,IKB+1) - PZZ(:,:,IKB))
ENDWHERE
!
END SUBROUTINE PT_DISCHARGE
!
!----------------------------------------------------------------------------------
!
  SUBROUTINE WRITE_OUT_ASCII
!
!!
!! Purpose:
!!
!
!*      0.     DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
INTEGER :: I1, I2
INTEGER :: ILU ! unit number for IO
!
!
!*      1.     FLASH PARAMETERS
!              ----------------
!
ILU = TPFILE_FGEOM_DIAG%NLU
!
! Ecriture ascii dans CEXP//'_fgeom_diag.asc" defini dans RESOLVED_ELEC
!
IF (LCARTESIAN) THEN
  DO I1 = 1, NNBLIGHT
    WRITE (UNIT=ILU,FMT='(I8,F9.1,I4,I6,I4,I6,F9.3,F12.3,F12.3,F9.3,F8.2,F9.2,f9.4)') &
          ISFLASH_NUMBER(I1),              &
          ISTCOUNT_NUMBER(I1) * PTSTEP,    &
          ISCELL_NUMBER(I1),               &
          ISNB_FLASH(I1),                  &
          ISTYPE(I1),                      &
          ISNBSEG(I1),                     &
          ZSEM_TRIG(I1),                   &
          ZSCOORD_SEG(I1,1,1)*1.E-3,       &
          ZSCOORD_SEG(I1,1,2)*1.E-3,       &
          ZSCOORD_SEG(I1,1,3)*1.E-3,       &
          ZSNEUT_POS(I1),                  &
          ZSNEUT_NEG(I1), ZSNEUT_POS(I1)+ZSNEUT_NEG(I1)
  END DO
ELSE
  DO I1 = 1, NNBLIGHT
! compute latitude and longitude of the triggering point
    CALL SM_LATLON(XLATORI,XLONORI,ZSCOORD_SEG(I1,1,1),&
                                   ZSCOORD_SEG(I1,1,2),&
                                   ZLAT,ZLON)
!
    WRITE (UNIT=ILU,FMT='(I8,F9.1,I4,I6,I4,I6,F9.3,F12.3,F12.3,F9.3,F8.2,F9.2,f9.4)') &
          ISFLASH_NUMBER(I1),              &
          ISTCOUNT_NUMBER(I1) * PTSTEP,    &
          ISCELL_NUMBER(I1),               &
          ISNB_FLASH(I1),                  &
          ISTYPE(I1),                      &
          ISNBSEG(I1),                     &
          ZSEM_TRIG(I1),                   &
          ZLAT,                            &
          ZLON,                            &
          ZSCOORD_SEG(I1,1,3)*1.E-3,       &
          ZSNEUT_POS(I1),                  &
          ZSNEUT_NEG(I1), ZSNEUT_POS(I1)+ZSNEUT_NEG(I1)
  END DO
END IF
!
FLUSH(UNIT=ILU)
!
!
!*      2.     FLASH SEGMENT COORDINATES
!              -------------------------
!
IF (LSAVE_COORD) THEN
!
! Ecriture ascii dans CEXP//'_fgeom_coord.asc" defini dans RESOLVED_ELEC
!
  ILU = TPFILE_FGEOM_COORD%NLU
!
  DO I1 = 1, NNBLIGHT
    DO I2 = 1, ISNBSEG(I1)
      WRITE (ILU, FMT='(I4,F9.1,I4,F12.3,F12.3,F12.3)') &
                 ISFLASH_NUMBER(I1),           & 
                 ISTCOUNT_NUMBER(I1) * PTSTEP, & 
                 ISTYPE(I1),                   &
                 ZSCOORD_SEG(I1,I2,1)*1.E-3,   &
                 ZSCOORD_SEG(I1,I2,2)*1.E-3,   &
                 ZSCOORD_SEG(I1,I2,3)*1.E-3
    END DO
  END DO
!
  FLUSH(UNIT=ILU)
END IF
!
END SUBROUTINE WRITE_OUT_ASCII
!
!-------------------------------------------------------------------------------
!
SUBROUTINE WRITE_OUT_LMA
!
!!
!! Purpose:
!!
!
!*      0.     DECLARATIONS
!              ------------
!
IMPLICIT NONE
!
INTEGER :: I1, I2
INTEGER :: ILU ! unit number for IO
!
!
!*      1.     LMA SIMULATOR
!              -------------
!
CALL SM_LATLON(XLATORI,XLONORI,ZSCOORD_SEG(:,:,1),ZSCOORD_SEG(:,:,2), &
                               ZLMA_LAT(:,:),ZLMA_LON(:,:))
!
ILU = TPFILE_LMA%NLU
!
DO I1 = 1, NNBLIGHT
  DO I2 = 1, ISNBSEG(I1)
    WRITE (UNIT=ILU,FMT='(I6,F12.1,I6,2(F15.6),3(F15.3),3(I6),12(E15.4))') &
               ISFLASH_NUMBER(I1),           &
               ISTCOUNT_NUMBER(I1) * PTSTEP, &
               ISTYPE(I1),                   &
               ZLMA_LAT(I1,I2),              &
               ZLMA_LON(I1,I2),              &
               ZSCOORD_SEG(I1,I2,1)*1.E-3,   &
               ZSCOORD_SEG(I1,I2,2)*1.E-3,   &
               ZSCOORD_SEG(I1,I2,3)*1.E-3,   &
               ISLMA_SEG_GLOB(I1,I2,1),      &
               ISLMA_SEG_GLOB(I1,I2,2),      &
               ISLMA_SEG_GLOB(I1,I2,3),      &
               ZSLMA_PRT(I1,I2,2),           &
               ZSLMA_PRT(I1,I2,3),           &
               ZSLMA_PRT(I1,I2,4),           &
               ZSLMA_PRT(I1,I2,5),           &
               ZSLMA_PRT(I1,I2,6),           &
               ZSLMA_QMT(I1,I2,2),           &
               ZSLMA_QMT(I1,I2,3),           &
               ZSLMA_QMT(I1,I2,4),           &
               ZSLMA_QMT(I1,I2,5),           &
               ZSLMA_QMT(I1,I2,6),           &
               ZSLMA_NEUT_POS(I1,I2),        &
               ZSLMA_NEUT_NEG(I1,I2)
  END DO
END DO
!
FLUSH(UNIT=ILU)
!
END SUBROUTINE WRITE_OUT_LMA
!
!-------------------------------------------------------------------------------
!
RECURSIVE SUBROUTINE N8QUICK_SORT(PLIST, KORDER)

! Quick sort routine from:
! Brainerd, W.S., Goldberg, C.H. & Adams, J.C. (1990) "Programmer's Guide to
! Fortran 90", McGraw-Hill  ISBN 0-07-000248-7, pages 149-150.
! Modified by Alan Miller to include an associated integer array which gives
! the positions of the elements in the original order.
!
use modd_precision, only: MNHINT64

IMPLICIT NONE
!
INTEGER(kind=MNHINT64), DIMENSION (:), INTENT(INOUT)  :: PLIST
INTEGER, DIMENSION (:), INTENT(OUT)  :: KORDER
!
! Local variable
INTEGER :: JI

DO JI = 1, SIZE(PLIST)
  KORDER(JI) = JI
END DO

CALL N8QUICK_SORT_1(1, SIZE(PLIST), PLIST, KORDER)

END SUBROUTINE N8QUICK_SORT
!
!-------------------------------------------------------------------------------
!
RECURSIVE SUBROUTINE N8QUICK_SORT_1(KLEFT_END, KRIGHT_END, PLIST1, KORDER1)

use modd_precision, only: MNHINT64

implicit none

INTEGER,                               INTENT(IN)    :: KLEFT_END, KRIGHT_END
INTEGER(kind=MNHINT64), DIMENSION (:), INTENT(INOUT) :: PLIST1
INTEGER,                DIMENSION (:), INTENT(INOUT) :: KORDER1
!     Local variables
INTEGER, PARAMETER     :: IMAX_SIMPLE_SORT_SIZE = 6

INTEGER                :: JI, JJ, ITEMP
INTEGER(kind=MNHINT64) :: ZREF, ZTEMP

IF (KRIGHT_END < KLEFT_END + IMAX_SIMPLE_SORT_SIZE) THEN
  ! Use interchange sort for small PLISTs
  CALL N8INTERCHANGE_SORT(KLEFT_END, KRIGHT_END, PLIST1, KORDER1)
  !
ELSE
  !
  ! Use partition ("quick") sort
  ! valeur au centre du tableau
  ZREF = PLIST1((KLEFT_END + KRIGHT_END)/2)
  JI = KLEFT_END - 1
  JJ = KRIGHT_END + 1

  DO
    ! Scan PLIST from left end until element >= ZREF is found
    DO
      JI = JI + 1
      IF (PLIST1(JI) >= ZREF) EXIT
    END DO
    ! Scan PLIST from right end until element <= ZREF is found
    DO
      JJ = JJ - 1
      IF (PLIST1(JJ) <= ZREF) EXIT
    END DO


    IF (JI < JJ) THEN
      ! Swap two out-of-order elements
      ZTEMP = PLIST1(JI)
      PLIST1(JI) = PLIST1(JJ)
      PLIST1(JJ) = ZTEMP
      ITEMP = KORDER1(JI)
      KORDER1(JI) = KORDER1(JJ)
      KORDER1(JJ) = ITEMP
    ELSE IF (JI == JJ) THEN
      JI = JI + 1
      EXIT
    ELSE
      EXIT
    END IF
  END DO

  IF ( KLEFT_END < JJ )         CALL N8QUICK_SORT_1( KLEFT_END, JJ,         PLIST1, KORDER1 )
  IF ( JI        < KRIGHT_END ) CALL N8QUICK_SORT_1( JI,        KRIGHT_END, PLIST1, KORDER1 )
END IF

END SUBROUTINE N8QUICK_SORT_1
!
!-------------------------------------------------------------------------------
!
SUBROUTINE N8INTERCHANGE_SORT(KLEFT_END, KRIGHT_END, PLIST2, KORDER2)

use modd_precision, only: MNHINT64

implicit none

INTEGER,                              INTENT(IN)    :: KLEFT_END, KRIGHT_END
INTEGER(kind=MNHINT64), DIMENSION(:), INTENT(INOUT) :: PLIST2
INTEGER,                DIMENSION(:), INTENT(INOUT) :: KORDER2
!     Local variables
INTEGER                :: JI, JJ, ITEMP
INTEGER(kind=MNHINT64) :: ZTEMP

! boucle sur tous les points
DO JI = KLEFT_END, KRIGHT_END - 1
  !
  ! boucle sur les points suivants le point JI
  DO JJ = JI+1, KRIGHT_END
    !
    ! si la distance de JI au point est plus grande que celle de JJ
    IF (PLIST2(JI) > PLIST2(JJ)) THEN
      ! distance de JI au point (la plus grande)
      ZTEMP = PLIST2(JI)
      ! le point JJ est déplacé à l'indice JI dans le tableau 
      PLIST2(JI) = PLIST2(JJ)
      ! le point JI est déplacé à l'indice JJ dans le tableau
      PLIST2(JJ) = ZTEMP
      ! indice du point JI dans le tableau
      ITEMP = KORDER2(JI)
      ! l'indice du point JJ est mis à la place JI
      KORDER2(JI) = KORDER2(JJ)
      ! l'indice du point JI est mis à la place JJ
      KORDER2(JJ) = ITEMP
    END IF
    !
  END DO
  !
END DO

END SUBROUTINE N8INTERCHANGE_SORT
!-------------------------------------------------------------------------------
  SUBROUTINE MNH_RANDOM_NUMBER(ZRANDOM)

    use modd_precision, only: MNHINT32

    REAL                         :: ZRANDOM
    INTEGER(kind=MNHINT32), SAVE :: NSEED_MNH = 26032012_MNHINT32

    ZRANDOM = real( r8_uniform_01( NSEED_MNH ), kind(ZRANDOM) )

  END SUBROUTINE MNH_RANDOM_NUMBER

!------------------------------------------------------------------------------------------

  FUNCTION r8_uniform_01 ( seed )

    !*****************************************************************************80
    !
    !! R8_UNIFORM_01 returns a unit pseudorandom R8.
    !
    !  Discussion:
    !
    !    An R8 is a real ( kind = 8 ) value.
    !
    !    For now, the input quantity SEED is an integer variable.
    !
    !    This routine implements the recursion
    !
    !      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
    !      r8_uniform_01 = seed / ( 2^31 - 1 )
    !
    !    The integer arithmetic never requires more than 32 bits,
    !    including a sign bit.
    !
    !    If the initial seed is 12345, then the first three computations are
    !
    !      Input     Output      R8_UNIFORM_01
    !      SEED      SEED
    !
    !         12345   207482415  0.096616
    !     207482415  1790989824  0.833995
    !    1790989824  2035175616  0.947702
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.     
    !    Souce here : https://people.sc.fsu.edu/~jburkardt/f_src/uniform/uniform.f90
    !
    !  Modified:
    !
    !    31 May 2007
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Reference:
    !
    !    Paul Bratley, Bennett Fox, Linus Schrage,
    !    A Guide to Simulation,
    !    Second Edition,
    !    Springer, 1987,
    !    ISBN: 0387964673,
    !    LC: QA76.9.C65.B73.
    !
    !    Bennett Fox,
    !    Algorithm 647:
    !    Implementation and Relative Efficiency of Quasirandom
    !    Sequence Generators,
    !    ACM Transactions on Mathematical Software,
    !    Volume 12, Number 4, December 1986, pages 362-376.
    !
    !    Pierre L'Ecuyer,
    !    Random Number Generation,
    !    in Handbook of Simulation,
    !    edited by Jerry Banks,
    !    Wiley, 1998,
    !    ISBN: 0471134031,
    !    LC: T57.62.H37.
    !
    !    Peter Lewis, Allen Goodman, James Miller,
    !    A Pseudo-Random Number Generator for the System/360,
    !    IBM Systems Journal,
    !    Volume 8, Number 2, 1969, pages 136-143.
    !
    !  Parameters:
    !
    !    Input/output, integer ( kind = MNHINT32 ) SEED, the "seed" value, which should
    !    NOT be 0. On output, SEED has been updated.
    !
    !    Output, real ( kind = MNHREAL64 ) R8_UNIFORM_01, a new pseudorandom variate,
    !    strictly between 0 and 1.
    !
    use modd_precision, only: MNHINT32, MNHREAL64

    implicit none

    integer(kind = MNHINT32), intent(inout) :: seed
    real(kind=MNHREAL64)                    :: r8_uniform_01

    integer(kind = MNHINT32), parameter :: i4_huge = 2147483647_MNHINT32

    integer(kind = MNHINT32) :: k

    if ( seed == 0_MNHINT32 ) THEN
      call Print_msg( NVERB_FATAL, 'GEN', 'r8_uniform_01', 'seed dummy argument must be different of 0' )
    end if

    k = seed / 127773_MNHINT32

    seed = 16807_MNHINT32 * ( seed - k * 127773_MNHINT32 ) - k * 2836_MNHINT32

    if ( seed < 0_MNHINT32 ) then
       seed = seed + i4_huge
    end if

    r8_uniform_01 = real(seed) * 4.656612875d-10

    return
  end function r8_uniform_01
!
END SUBROUTINE FLASH_GEOM_ELEC_n
!
!-------------------------------------------------------------------------------
