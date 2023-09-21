!MNH_LIC Copyright 1998-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_INI_LB
!     ######################
!
INTERFACE 
!
SUBROUTINE INI_LB(TPINIFILE,OLSOURCE,KSV,                                   &
     KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,                     &
     KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                         &
     KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,                 &
     HGETTKEM,HGETRVM,HGETRCM,HGETRRM,HGETRIM,HGETRSM,                      &
     HGETRGM,HGETRHM,HGETSVM,                                               &
     PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,                  &
     PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,                  &
     PLBXUMM,PLBXVMM,PLBXWMM,PLBXTHMM,PLBXTKEMM,PLBXRMM,PLBXSVMM,           &
     PLBYUMM,PLBYVMM,PLBYWMM,PLBYTHMM,PLBYTKEMM,PLBYRMM,PLBYSVMM,           &
     PLENG )
!
USE MODD_IO, ONLY: TFILEDATA
!
TYPE(TFILEDATA),       INTENT(IN)   :: TPINIFILE ! Initial file
LOGICAL,               INTENT(IN)   :: OLSOURCE  ! switch for the source term
! Larger Scale fields (source if OLSOURCE=T,  fields at time t-dt if OLSOURCE=F) :
INTEGER,               INTENT(IN)   :: KSV       ! number of passive variables
! sizes of the West-east total LB area
INTEGER, INTENT(IN) :: KSIZELBX_ll,KSIZELBXU_ll      ! for T,V,W and u 
INTEGER, INTENT(IN) :: KSIZELBXTKE_ll                ! for TKE 
INTEGER, INTENT(IN) :: KSIZELBXR_ll,KSIZELBXSV_ll    ! for Rx and SV    
! sizes of the North-south total LB area
INTEGER, INTENT(IN) :: KSIZELBY_ll,KSIZELBYV_ll      ! for T,U,W  and v
INTEGER, INTENT(IN) :: KSIZELBYTKE_ll                ! for TKE
INTEGER, INTENT(IN) :: KSIZELBYR_ll,KSIZELBYSV_ll    ! for Rx and SV
! Get indicators
CHARACTER (LEN=*),         INTENT(IN)  :: HGETTKEM,                          &
                                          HGETRVM,HGETRCM,HGETRRM,           &
                                          HGETRIM,HGETRSM,HGETRGM,HGETRHM
CHARACTER (LEN=*), DIMENSION(:),INTENT(IN)  :: HGETSVM
! LB  fields (source if OLSOURCE=T,  fields at time t-dt if OLSOURCE=F) :
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXUM,PLBXVM,PLBXWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYUM,PLBYVM,PLBYWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTKEM          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTKEM
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBXRM  ,PLBXSVM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBYRM  ,PLBYSVM  ! in x and y-dir.
! LB arrays at time t-dt (if OLSOURCE=T) : 
REAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL  :: PLBXUMM,PLBXVMM,PLBXWMM ! Wind
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBXTHMM              ! Mass
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBYUMM,PLBYVMM,PLBYWMM ! Wind
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBYTHMM              ! Mass
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBXTKEMM           ! TKE
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBYTKEMM
REAL, DIMENSION(:,:,:,:),INTENT(IN), OPTIONAL  :: PLBXRMM  ,PLBXSVMM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),INTENT(IN), OPTIONAL  :: PLBYRMM  ,PLBYSVMM  ! in x and y-dir.
REAL,                  INTENT(IN),   OPTIONAL :: PLENG    ! Interpolation length
!
END SUBROUTINE INI_LB
!
END INTERFACE
!
END MODULE MODI_INI_LB
!     ############################################################
SUBROUTINE INI_LB(TPINIFILE,OLSOURCE,KSV,                          &
     KSIZELBX_ll,KSIZELBXU_ll,KSIZELBY_ll,KSIZELBYV_ll,            &
     KSIZELBXTKE_ll,KSIZELBYTKE_ll,                                &
     KSIZELBXR_ll,KSIZELBYR_ll,KSIZELBXSV_ll,KSIZELBYSV_ll,        &
     HGETTKEM,HGETRVM,HGETRCM,HGETRRM,HGETRIM,HGETRSM,             &
     HGETRGM,HGETRHM,HGETSVM,                                      &
     PLBXUM,PLBXVM,PLBXWM,PLBXTHM,PLBXTKEM,PLBXRM,PLBXSVM,         &
     PLBYUM,PLBYVM,PLBYWM,PLBYTHM,PLBYTKEM,PLBYRM,PLBYSVM,         &
     PLBXUMM,PLBXVMM,PLBXWMM,PLBXTHMM,PLBXTKEMM,PLBXRMM,PLBXSVMM,  &
     PLBYUMM,PLBYVMM,PLBYWMM,PLBYTHMM,PLBYTKEMM,PLBYRMM,PLBYSVMM,  &
     PLENG )
!     ############################################################
!
!!****  *INI_LB* - routine to initialize  LB fields
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to read the LB fields and to distribute
! on subdomain which have a non-nul intersection with the LB areas.
!       In case of OLSOURCE=T, it initializes the LB sources instead of the
!   LB fields at time t-dt
!
!!**  METHOD
!!    ------
!!    The LB fields are read in file and distributed by FMREAD_LB
!!
!!    In case of OLSOURCE=T (INI_LB called by INI_CPL or LS_COUPLING), the LB sources
!!   are computed
!!     
!!
!!    EXTERNAL
!!    --------
!!      FMREAD    : to read data in LFIFM file
!!      FMREAD_LB : to read LB data in LFIFM file
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------ 
!!      Module MODD_CONF   : NVERB
!!
!!    REFERENCE
!!    ---------
!!      Book2 of the documentation (routine INI_LB)
!!      
!!
!!    AUTHOR
!!    ------
!!  	V. Ducrocq       * Meteo France *
!!      D. Gazen         L.A. 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        22/09/98    FMREAD_LB handle LBs fields
!!      J. Stein        18/09/99    problem with the dry case
!!      D. Gazen        22/01/01    treat NSV_* with floating indices
!!      F Gheusi        29/10/03    bug in LB sources for NSV
!!      J.-P. Pinty     06/05/04    treat NSV_* for C1R3 and ELEC
!!                      20/05/06    Remove KEPS
!!      C.Lac           20/03/08    Add passive pollutants
!!      M.Leriche       16/07/10    Add ice phase chemical species
!!      Pialat/tulet    15/02/12    Add ForeFire scalars 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      M.Leriche       09/02/16    Treat gas and aq. chemicals separately
!!      J.Escobar : 27/04/2016 : bug , test only on ANY(HGETSVM({{1:KSV}})=='READ'
!!      J.-P. Pinty     09/02/16    Add LIMA that is LBC for CCN and IFN
!!      M.Leriche       09/02/16    Treat gas and aq. chemicals separately
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 13/02/2019: initialize PLBXSVM and PLBYSVM in all cases
!  S. Bielli      02/2019: Sea salt: significant sea wave height influences salt emission; 5 salt modes
!  P. Wautelet 04/02/2022: use TSVLIST to manage metadata of scalar variables
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_TURB_n,          ONLY: XTKEMIN
USE MODD_CONF,            ONLY: LCPL_AROME
use modd_field,           only: NMNHDIM_UNKNOWN, tfieldmetadata, TYPELOG, TYPEREAL
USE MODD_IO,              ONLY: TFILEDATA
USE MODD_NSV,             ONLY: NSV, NSV_CS, NSV_CSBEG, NSV_CSEND, NSV_LIMA_BEG, NSV_LIMA_END,      &
#ifdef MNH_FOREFIRE
                                NSV_FF, NSV_FFBEG, NSV_FFEND,                                       &
#endif
                                NSV_LIMA_CCN_FREE, NSV_LIMA_IFN_FREE, NSV_PP, NSV_PPBEG, NSV_PPEND, &
                                NSV_SNWBEG, NSV_SNWEND, NSV_USER, TSVLIST
USE MODD_PARAMETERS,      ONLY: JPHEXT, JPSVNAMELGTMAX, NLONGNAMELGTMAX, NMNHNAMELGTMAX
USE MODD_PARAM_LIMA,      ONLY: NMOD_CCN, NMOD_IFN
!
USE MODE_IO_FIELD_READ,   only: IO_Field_read, IO_Field_read_lb
USE MODE_MSG
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
!
TYPE(TFILEDATA),       INTENT(IN)   :: TPINIFILE ! Initial file
LOGICAL,               INTENT(IN)   :: OLSOURCE  ! switch for the source term
! Larger Scale fields (source if OLSOURCE=T,  fields at time t-dt if OLSOURCE=F) :
INTEGER,               INTENT(IN)   :: KSV       ! number of passive variables
! sizes of the West-east total LB area
INTEGER, INTENT(IN) :: KSIZELBX_ll,KSIZELBXU_ll      ! for T,V,W and u 
INTEGER, INTENT(IN) :: KSIZELBXTKE_ll                !  for TKE 
INTEGER, INTENT(IN) :: KSIZELBXR_ll,KSIZELBXSV_ll    ! for Rx and SV    
! sizes of the North-south total LB area
INTEGER, INTENT(IN) :: KSIZELBY_ll,KSIZELBYV_ll      ! for T,U,W  and v
INTEGER, INTENT(IN) :: KSIZELBYTKE_ll                ! for TKE 
INTEGER, INTENT(IN) :: KSIZELBYR_ll,KSIZELBYSV_ll    ! for Rx and SV
! Get indicators
CHARACTER (LEN=*),         INTENT(IN)  :: HGETTKEM,                          &
                                          HGETRVM,HGETRCM,HGETRRM,           &
                                          HGETRIM,HGETRSM,HGETRGM,HGETRHM
CHARACTER (LEN=*), DIMENSION(:),INTENT(IN)  :: HGETSVM
! LB  fields (source if OLSOURCE=T,  fields at time t-dt if OLSOURCE=F) :
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXUM,PLBXVM,PLBXWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYUM,PLBYVM,PLBYWM ! Wind
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTHM              ! Mass
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBXTKEM          ! TKE
REAL, DIMENSION(:,:,:),          INTENT(OUT) :: PLBYTKEM          ! 
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBXRM  ,PLBXSVM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),        INTENT(OUT) :: PLBYRM  ,PLBYSVM  ! in x and y-dir.
! LB arrays at time t-dt (if OLSOURCE=T) : 
REAL, DIMENSION(:,:,:), INTENT(IN), OPTIONAL  :: PLBXUMM,PLBXVMM,PLBXWMM ! Wind
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBXTHMM              ! Mass
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBYUMM,PLBYVMM,PLBYWMM ! Wind
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBYTHMM              ! Mass
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBXTKEMM           ! TKE
REAL, DIMENSION(:,:,:),INTENT(IN), OPTIONAL  :: PLBYTKEMM
REAL, DIMENSION(:,:,:,:),INTENT(IN), OPTIONAL  :: PLBXRMM  ,PLBXSVMM  ! Moisture and SV
REAL, DIMENSION(:,:,:,:),INTENT(IN), OPTIONAL  :: PLBYRMM  ,PLBYSVMM  ! in x and y-dir.
REAL,                  INTENT(IN),   OPTIONAL :: PLENG    ! Interpolation length
!
!
!*       0.2   declarations of local variables
!
INTEGER             :: ILBSIZEX,ILBSIZEY   ! depth  of the LB area in the RIM direction 
                                           !  written in FM file
INTEGER       :: IL3DX,IL3DY ! Size of the LB arrays in FM file
                             ! in  the RIM direction
INTEGER       :: IL3DXU,IL3DYV  ! Size of the LB arrays in FM file
                                ! in  the RIM direction for the normal wind
INTEGER             :: IRIMX,IRIMY ! Total size of the LB area (for the RIM direction)
INTEGER             :: IRIMXU,IRIMYV ! Total size of the LB area (for the RIM direction)
                                     ! for the normal wind (spatial gradient needed)

INTEGER             :: JSV,JRR                    ! Loop index for MOIST AND 
                                                  !  additional scalar variables 
INTEGER             :: IRR                        !  counter for moist variables
INTEGER             :: IRESP
LOGICAL :: GHORELAX_UVWTH  ! switch for the horizontal relaxation for U,V,W,TH in the FM file
LOGICAL :: GHORELAX_TKE    ! switch for the horizontal relaxation for tke in the FM file
LOGICAL :: GHORELAX_R, GHORELAX_SV ! switch for the horizontal relaxation 
                                   ! for moist and scalar variables
LOGICAL :: GIS551 ! True if file was written with MNH 5.5.1
LOGICAL :: GOLDFILEFORMAT
CHARACTER (LEN= LEN(HGETRVM)), DIMENSION (7) :: YGETRXM ! Arrays with  the get indicators 
                                                        !  for the moist variables
CHARACTER (LEN=1), DIMENSION (7) :: YC    ! array with the prefix of the moist variables
CHARACTER(LEN=NMNHNAMELGTMAX)  :: YMNHNAME_BASE
CHARACTER(LEN=NLONGNAMELGTMAX) :: YLONGNAME_BASE
TYPE(TFIELDMETADATA) :: TZFIELD
!-------------------------------------------------------------------------------
!
!
!*       0.    READ CPL_AROME to know which LB_fileds there are to read
!              --------------------
IF ((TPINIFILE%NMNHVERSION(1)==4 .AND. TPINIFILE%NMNHVERSION(2)>8) .OR. TPINIFILE%NMNHVERSION(1)>4) THEN
  CALL IO_Field_read(TPINIFILE,'CPL_AROME',LCPL_AROME)
ELSE
  LCPL_AROME=.FALSE.
ENDIF
!
!
!*       1.    SOME INITIALIZATIONS
!              --------------------
!
!If TPINIFILE file was written with a MesoNH version < 5.6, some variables had different names or were not available
GOLDFILEFORMAT = (        TPINIFILE%NMNHVERSION(1) < 5                                       &
                   .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) < 6 ) )
GIS551 = TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) == 5 .AND. TPINIFILE%NMNHVERSION(3) == 1
!
!
!-------------------------------------------------------------------------------
!
!*       2.    READ 2D "surfacic" LB fields 
!              ----------------------------
!
!*       2.1   read the number of available points for the horizontal relaxation
! for basic variables 
CALL IO_Field_read(TPINIFILE,'RIMX',ILBSIZEX)
CALL IO_Field_read(TPINIFILE,'RIMY',ILBSIZEY)
!
!*        2.2 Basic variables
! 
CALL IO_Field_read(TPINIFILE,'HORELAX_UVWTH',GHORELAX_UVWTH)
                                !
IF (GHORELAX_UVWTH) THEN 
  IRIMX =(KSIZELBX_ll-2*JPHEXT)/2   
  IRIMXU=(KSIZELBXU_ll-2*JPHEXT)/2  
  IRIMY =(KSIZELBY_ll-2*JPHEXT)/2
  IRIMYV=(KSIZELBYV_ll-2*JPHEXT)/2
  IL3DX=2*ILBSIZEX+2*JPHEXT
  IL3DXU=IL3DX
  IL3DY=2*ILBSIZEY+2*JPHEXT
  IL3DYV=IL3DY
ELSE
  IRIMX=0
  IRIMXU=1
  IRIMY=0
  IRIMYV=1
  IL3DX=2*JPHEXT ! 2
  IL3DY=2*JPHEXT ! 2
  IL3DXU=2 + 2*JPHEXT ! 4 
  IL3DYV=2 + 2*JPHEXT ! 4 
ENDIF
!
IF ( KSIZELBXU_ll /= 0 ) CALL IO_Field_read_lb( TPINIFILE, 'LBXUM',  IL3DXU, IRIMXU, PLBXUM  )
IF ( KSIZELBX_ll  /= 0 ) CALL IO_Field_read_lb( TPINIFILE, 'LBXVM',  IL3DX,  IRIMX,  PLBXVM  )
IF ( KSIZELBX_ll  /= 0 ) CALL IO_Field_read_lb( TPINIFILE, 'LBXWM',  IL3DX,  IRIMX,  PLBXWM  )
IF ( KSIZELBY_ll  /= 0 ) CALL IO_Field_read_lb( TPINIFILE, 'LBYUM',  IL3DY,  IRIMY,  PLBYUM  )
IF ( KSIZELBYV_ll /= 0 ) CALL IO_Field_read_lb( TPINIFILE, 'LBYVM',  IL3DYV, IRIMYV, PLBYVM  )
IF ( KSIZELBY_ll  /= 0 ) CALL IO_Field_read_lb( TPINIFILE, 'LBYWM',  IL3DY,  IRIMY,  PLBYWM  )
IF ( KSIZELBX_ll  /= 0 ) CALL IO_Field_read_lb( TPINIFILE, 'LBXTHM', IL3DX,  IRIMX,  PLBXTHM )
IF ( KSIZELBY_ll  /= 0 ) CALL IO_Field_read_lb( TPINIFILE, 'LBYTHM', IL3DY,  IRIMY,  PLBYTHM )
!
!*        2.3  LB-TKE
!
SELECT CASE(HGETTKEM)                 
CASE('READ') 
  IF (.NOT. LCPL_AROME .AND. OLSOURCE) THEN
    IF (PRESENT(PLBXTKEMM).AND.PRESENT(PLBYTKEMM)) THEN
      CALL PRINT_MSG( NVERB_INFO, 'IO', 'INI_LB', 'LBXTKES and LBYTKE are initialized to PLBXTKEMM and PLBYTKEMM' )
      PLBXTKEM(:,:,:) = PLBXTKEMM(:,:,:)    
      PLBYTKEM(:,:,:) = PLBYTKEMM(:,:,:)
    ELSE
      CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize LBXTKES and LBYTKES')
    ENDIF
  ELSE
    CALL IO_Field_read(TPINIFILE,'HORELAX_TKE',GHORELAX_TKE)
    IF (GHORELAX_TKE) THEN 
      IRIMX=(KSIZELBXTKE_ll-2*JPHEXT)/2   
      IRIMY=(KSIZELBYTKE_ll-2*JPHEXT)/2
      IL3DX=2*ILBSIZEX+2*JPHEXT
      IL3DY=2*ILBSIZEY+2*JPHEXT
    ELSE
      IRIMX=0
      IRIMY=0
      IL3DX=2*JPHEXT ! 2
      IL3DY=2*JPHEXT ! 2
    ENDIF
!
    IF (KSIZELBXTKE_ll /= 0) THEN
      CALL IO_Field_read_lb(TPINIFILE,'LBXTKEM',IL3DX,IRIMX,PLBXTKEM)
    END IF
!
    IF (KSIZELBYTKE_ll /= 0) THEN  
      CALL IO_Field_read_lb(TPINIFILE,'LBYTKEM',IL3DY,IRIMY,PLBYTKEM)
    END IF
  ENDIF
CASE('INIT')
  IF (SIZE(PLBXTKEM,1) /= 0) PLBXTKEM(:,:,:) = XTKEMIN
  IF (SIZE(PLBYTKEM,1) /= 0) PLBYTKEM(:,:,:) = XTKEMIN
END SELECT
!
! 
!*        2.5 LB-Rx
!
IF(KSIZELBXR_ll  > 0 ) THEN
  TZFIELD = TFIELDMETADATA(                                      &
    CMNHNAME   = 'HORELAX_R',                                    &
    CSTDNAME   = '',                                             &
    CLONGNAME  = 'HORELAX_R',                                    &
    CUNITS     = '',                                             &
    CDIR       = '--',                                           &
    CCOMMENT   = 'Switch to activate the HOrizontal RELAXation', &
    CLBTYPE    = 'NONE',                                         &
    NGRID      = 1,                                              &
    NTYPE      = TYPELOG,                                        &
    NDIMS      = 0,                                              &
    LTIMEDEP   = .FALSE.                                         )
  !
  CALL IO_Field_read(TPINIFILE,TZFIELD,GHORELAX_R)
  !
  YGETRXM(:)=(/HGETRVM,HGETRCM,HGETRRM,HGETRIM,HGETRSM,HGETRGM,HGETRHM/)
  YC(:)=(/"V","C","R","I","S","G","H"/)
  IF (GHORELAX_R) THEN 
    IRIMX=(KSIZELBXR_ll-2*JPHEXT)/2  
    IRIMY= (KSIZELBYR_ll-2*JPHEXT)/2  
    IL3DX=2*ILBSIZEX+2*JPHEXT
    IL3DY=2*ILBSIZEY+2*JPHEXT
  ELSE
    IRIMX=0
    IRIMY=0
    IL3DX=2*JPHEXT ! 2
    IL3DY=2*JPHEXT ! 2
  END IF
  !
  TZFIELD = TFIELDMETADATA( &
    CUNITS     = 'kg kg-1', &
    CDIR       = '',        &
    NGRID      = 1,         &
    NTYPE      = TYPEREAL,  &
    NDIMS      = 3,         &
    LTIMEDEP   = .TRUE.     )
  !
  IRR=0
  JRR=1
  SELECT CASE(YGETRXM(1))                 
    CASE('READ') 
      IRR=IRR+1 
      IF ( KSIZELBXR_ll  /= 0 ) THEN
        TZFIELD%CMNHNAME   = 'LBXR'//YC(JRR)//'M'
        TZFIELD%CLONGNAME  = 'LBXR'//YC(JRR)//'M'
        TZFIELD%CLBTYPE    = 'LBX'
        TZFIELD%CCOMMENT   = '2_Y_Z_LBXR'//YC(JRR)//'M'
        CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXRM(:,:,:,IRR))
      END IF
      !
      IF ( KSIZELBYR_ll /= 0 ) THEN
        TZFIELD%CMNHNAME   = 'LBYR'//YC(JRR)//'M'
        TZFIELD%CLONGNAME  = 'LBYR'//YC(JRR)//'M'
        TZFIELD%CLBTYPE    = 'LBY'
        TZFIELD%CCOMMENT   = '2_Y_Z_LBYR'//YC(JRR)//'M'
        CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYRM(:,:,:,IRR))
      END IF
    CASE('INIT')
      IRR=IRR+1 
      IF ( SIZE(PLBXRM,1) /= 0 )  PLBXRM(:,:,:,IRR) = 0.
      IF ( SIZE(PLBYRM,1) /= 0 )  PLBYRM(:,:,:,IRR) = 0.
  END SELECT
    !
    !
  DO JRR=2,7
    SELECT CASE(YGETRXM(JRR))                 
    CASE('READ') 
      IRR=IRR+1 
      IF ( KSIZELBXR_ll  /= 0 ) THEN
        IF (.NOT. LCPL_AROME .AND. OLSOURCE) THEN
            IF (PRESENT(PLBXRMM)) THEN
              PLBXRM(:,:,:,IRR)=PLBXRMM(:,:,:,IRR)
              CALL PRINT_MSG( NVERB_INFO, 'IO', 'INI_LB', 'PLBXRM is initialized to PLBXRMM for LBXR'//YC(JRR)//'M' )
            ELSE
              CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize PLBXRM for LBXR'//YC(JRR)//'M')
            ENDIF
        ELSE
          TZFIELD%CMNHNAME   = 'LBXR'//YC(JRR)//'M'
          TZFIELD%CLONGNAME  = 'LBXR'//YC(JRR)//'M'
          TZFIELD%CLBTYPE    = 'LBX'
          TZFIELD%CCOMMENT   = '2_Y_Z_LBXR'//YC(JRR)//'M'
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXRM(:,:,:,IRR))
        ENDIF
      END IF
      !
      IF ( KSIZELBYR_ll /= 0 ) THEN
        IF (.NOT. LCPL_AROME .AND. OLSOURCE) THEN
            IF (PRESENT(PLBYRMM)) THEN
              PLBYRM(:,:,:,IRR)=PLBYRMM(:,:,:,IRR)
              CALL PRINT_MSG( NVERB_INFO, 'IO', 'INI_LB', 'PLBYRM is initialized to PLBYRMM for LBYR'//YC(JRR)//'M' )
            ELSE
              CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize PLBYRM for LBYR'//YC(JRR)//'M')
            ENDIF
         ELSE
           TZFIELD%CMNHNAME   = 'LBYR'//YC(JRR)//'M'
           TZFIELD%CLONGNAME  = 'LBYR'//YC(JRR)//'M'
           TZFIELD%CLBTYPE    = 'LBY'
           TZFIELD%CCOMMENT   = '2_Y_Z_LBYR'//YC(JRR)//'M'
           CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYRM(:,:,:,IRR))
         ENDIF
       END IF
    CASE('INIT')
      IRR=IRR+1 
      IF ( SIZE(PLBXRM,1) /= 0 )  PLBXRM(:,:,:,IRR) = 0.
      IF ( SIZE(PLBYRM,1) /= 0 )  PLBYRM(:,:,:,IRR) = 0.
    END SELECT
  END DO
END IF
!
!*        2.6    LB-Scalar Variables
!
IF (KSV > 0) THEN
  IF (ANY(HGETSVM(1:KSV)=='READ')) THEN
    TZFIELD = TFIELDMETADATA(    &
      CMNHNAME   = 'HORELAX_SV', &
      CSTDNAME   = '',           &
      CLONGNAME  = 'HORELAX_SV', &
      CUNITS     = '',           &
      CDIR       = '--',         &
      CCOMMENT   = '',           &
      CLBTYPE    = 'NONE',       &
      NGRID      = 0,            &
      NTYPE      = TYPELOG,      &
      NDIMS      = 0,            &
      LTIMEDEP   = .FALSE.       )
    CALL IO_Field_read( TPINIFILE, TZFIELD, GHORELAX_SV )

    IF ( GHORELAX_SV ) THEN
      IRIMX=(KSIZELBXSV_ll-2*JPHEXT)/2
      IRIMY=(KSIZELBYSV_ll-2*JPHEXT)/2
      IL3DX=2*ILBSIZEX+2*JPHEXT
      IL3DY=2*ILBSIZEY+2*JPHEXT
    ELSE
      IRIMX=0
      IRIMY=0
      IL3DX=2*JPHEXT
      IL3DY=2*JPHEXT
    END IF
  END IF
END IF

! Scalar variables
DO JSV = 1, NSV
  SELECT CASE( HGETSVM(JSV) )
    CASE ( 'READ' )
      TZFIELD = TSVLIST(JSV)
      TZFIELD%CDIR = ''
      TZFIELD%NDIMLIST(:) = NMNHDIM_UNKNOWN
      YMNHNAME_BASE  = TRIM( TZFIELD%CMNHNAME  )
      YLONGNAME_BASE = TRIM( TZFIELD%CLONGNAME )

      IF ( KSIZELBXSV_ll /= 0 ) THEN
        TZFIELD%CMNHNAME  = 'LBX_' // TRIM( YMNHNAME_BASE  )
        TZFIELD%CLONGNAME = 'LBX_' // TRIM( YLONGNAME_BASE )

        !Some variables were written with an other name in MesoNH < 5.6
        IF ( GOLDFILEFORMAT ) THEN
          IF ( JSV >= 1 .AND. JSV <= NSV_USER ) THEN
            WRITE( TZFIELD%CMNHNAME, '( A6, I3.3 )' ) 'LBXSVM',JSV
            TZFIELD%CSTDNAME  = ''
            TZFIELD%CLONGNAME = TRIM( TZFIELD%CMNHNAME )
          ELSE IF ( JSV >= NSV_LIMA_BEG .AND. JSV <= NSV_LIMA_END ) THEN
            ! Name was corrected in MNH 5.5.1
            IF ( .NOT. GIS551 ) CALL OLD_CMNHNAME_GENERATE_INTERN( TZFIELD%CMNHNAME, TZFIELD%CLONGNAME )
            TZFIELD%CSTDNAME  = ''
          ELSE IF ( JSV >= NSV_PPBEG .AND. JSV <= NSV_PPEND ) THEN
            TZFIELD%CMNHNAME  = 'LBX_PP'
            TZFIELD%CSTDNAME  = ''
            TZFIELD%CLONGNAME = 'LBX_PP'
            IF ( JSV == NSV_PPBEG .AND. NSV_PP > 1 ) THEN
              CMNHMSG(1) = 'reading older file (<5.6) for LBX_PP scalar variables'
              CMNHMSG(2) = 'they are bugged: there should be several LBX_PP variables'
              CMNHMSG(3) = 'but they were all written with the same name ''LBX_PP'''
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB' )
            END IF
#ifdef MNH_FOREFIRE
          ELSE IF ( JSV >= NSV_FFBEG .AND. JSV <= NSV_FFEND ) THEN
            TZFIELD%CMNHNAME  = 'LBX_FF'
            TZFIELD%CSTDNAME  = ''
            TZFIELD%CLONGNAME = 'LBX_FF'
            IF ( JSV == NSV_FFBEG .AND. NSV_FF > 1 ) THEN
              CMNHMSG(1) = 'reading older file (<5.6) for LBX_FF scalar variables'
              CMNHMSG(2) = 'they are bugged: there should be several LBX_FF variables'
              CMNHMSG(3) = 'but they were all written with the same name ''LBX_FF'''
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB' )
            END IF
#endif
          ELSE IF ( JSV >= NSV_CSBEG .AND. JSV <= NSV_CSEND ) THEN
            TZFIELD%CMNHNAME  = 'LBX_CS'
            TZFIELD%CSTDNAME  = ''
            TZFIELD%CLONGNAME = 'LBX_CS'
            IF ( JSV == NSV_CSBEG .AND. NSV_CS > 1 ) THEN
              CMNHMSG(1) = 'reading older file (<5.6) for LBX_CS scalar variables'
              CMNHMSG(2) = 'they are bugged: there should be several LBX_CS variables'
              CMNHMSG(3) = 'but they were all written with the same name ''LBX_CS'''
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB' )
            END IF
          END IF
        END IF

        WRITE( TZFIELD%CCOMMENT, '( A6, A6, I3.3 )' ) '2_Y_Z_', 'LBXSVM', JSV
        TZFIELD%CLBTYPE     = 'LBX'

        CALL IO_Field_read_lb( TPINIFILE, TZFIELD, IL3DX, IRIMX, PLBXSVM(:,:,:,JSV), IRESP )

        IF ( IRESP /= 0 ) THEN
          IF ( PRESENT( PLBXSVMM ) ) THEN
            PLBXSVM(:,:,:,JSV) = PLBXSVMM(:,:,:,JSV)
            CALL PRINT_MSG( NVERB_INFO, 'IO', 'INI_LB', 'PLBXSVM is initialized to PLBXSVMM for ' // TRIM( YMNHNAME_BASE ) )
          ELSE
            IF ( JSV >= NSV_LIMA_BEG .AND. JSV <= NSV_LIMA_END ) THEN
               PLBXSVM(:,:,:,JSV) = 0.
               CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB', 'PLBXSVM is initialized to 0 for ' // TRIM( YMNHNAME_BASE ) )
            ELSE IF ( ( JSV >= NSV_PPBEG  .AND. JSV <= NSV_PPEND ) .OR. &
#ifdef MNH_FOREFIRE
                      ( JSV >= NSV_FFBEG  .AND. JSV <= NSV_FFEND ) .OR. &
#endif
                      ( JSV >= NSV_CSBEG  .AND. JSV <= NSV_CSEND ) .OR. &
                      ( JSV >= NSV_SNWBEG .AND. JSV <= NSV_SNWEND .AND. GOLDFILEFORMAT ) ) THEN !Snow was not written in <5.6
              PLBXSVM(:,:,:,JSV) = 0.
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB', 'PLBXSVM is initialized to 0 for ' // TRIM( YMNHNAME_BASE ) )
            ELSE
              CALL PRINT_MSG( NVERB_FATAL, 'IO', 'INI_LB', 'problem to initialize PLBXSVM for ' // TRIM( YMNHNAME_BASE ) )
            END IF
          END IF
        END IF
      END IF

      IF ( KSIZELBYSV_ll /= 0 ) THEN
        TZFIELD%CMNHNAME  = 'LBY_' // TRIM( YMNHNAME_BASE  )
        TZFIELD%CLONGNAME = 'LBY_' // TRIM( YLONGNAME_BASE )

        !Some variables were written with an other name in MesoNH < 5.6
        IF ( GOLDFILEFORMAT ) THEN
          IF ( JSV >= 1 .AND. JSV <= NSV_USER ) THEN
            WRITE( TZFIELD%CMNHNAME, '( A6, I3.3 )' ) 'LBYSVM',JSV
            TZFIELD%CSTDNAME  = ''
            TZFIELD%CLONGNAME = TRIM( TZFIELD%CMNHNAME )
          ELSE IF ( JSV >= NSV_LIMA_BEG .AND. JSV <= NSV_LIMA_END ) THEN
            ! Name was corrected in MNH 5.5.1
            IF ( .NOT. GIS551 ) CALL OLD_CMNHNAME_GENERATE_INTERN( TZFIELD%CMNHNAME, TZFIELD%CLONGNAME )
            TZFIELD%CSTDNAME  = ''
          ELSE IF ( JSV >= NSV_PPBEG .AND. JSV <= NSV_PPEND ) THEN
            TZFIELD%CMNHNAME  = 'LBY_PP'
            TZFIELD%CSTDNAME  = ''
            TZFIELD%CLONGNAME = 'LBY_PP'
            IF ( JSV == NSV_PPBEG .AND. NSV_PP > 1 ) THEN
              CMNHMSG(1) = 'reading older file (<5.6) for LBY_PP scalar variables'
              CMNHMSG(2) = 'they are bugged: there should be several LBY_PP variables'
              CMNHMSG(3) = 'but they were all written with the same name ''LBY_PP'''
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB' )
            END IF
#ifdef MNH_FOREFIRE
          ELSE IF ( JSV >= NSV_FFBEG .AND. JSV <= NSV_FFEND ) THEN
            TZFIELD%CMNHNAME  = 'LBY_FF'
            TZFIELD%CSTDNAME  = ''
            TZFIELD%CLONGNAME = 'LBY_FF'
            IF ( JSV == NSV_FFBEG .AND. NSV_FF > 1 ) THEN
              CMNHMSG(1) = 'reading older file (<5.6) for LBY_FF scalar variables'
              CMNHMSG(2) = 'they are bugged: there should be several LBY_FF variables'
              CMNHMSG(3) = 'but they were all written with the same name ''LBY_FF'''
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB' )
            END IF
#endif
          ELSE IF ( JSV >= NSV_CSBEG .AND. JSV <= NSV_CSEND ) THEN
            TZFIELD%CMNHNAME  = 'LBY_CS'
            TZFIELD%CSTDNAME  = ''
            TZFIELD%CLONGNAME = 'LBY_CS'
            IF ( JSV == NSV_CSBEG .AND. NSV_CS > 1 ) THEN
              CMNHMSG(1) = 'reading older file (<5.6) for LBY_CS scalar variables'
              CMNHMSG(2) = 'they are bugged: there should be several LBY_CS variables'
              CMNHMSG(3) = 'but they were all written with the same name ''LBY_CS'''
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB' )
            END IF
          END IF
        END IF
        WRITE( TZFIELD%CCOMMENT, '( A6, A6, I3.3 )' ) 'X_2_Z_', 'LBYSVM', JSV
        TZFIELD%CLBTYPE     = 'LBY'

        CALL IO_Field_read_lb( TPINIFILE, TZFIELD, IL3DY, IRIMY, PLBYSVM(:,:,:,JSV), IRESP )

        IF ( IRESP /= 0 ) THEN
          IF ( PRESENT( PLBYSVMM ) ) THEN
            PLBYSVM(:,:,:,JSV) = PLBYSVMM(:,:,:,JSV)
            CALL PRINT_MSG( NVERB_INFO, 'IO', 'INI_LB', 'PLBYSVM is initialized to PLBYSVMM for ' // TRIM( YMNHNAME_BASE ) )
          ELSE
            IF ( JSV >= NSV_LIMA_BEG .AND. JSV <= NSV_LIMA_END ) THEN
              PLBYSVM(:,:,:,JSV) = 0.
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB', 'PLBYSVM is initialized to 0 for ' // TRIM( YMNHNAME_BASE ) )
            ELSE IF ( ( JSV >= NSV_PPBEG  .AND. JSV <= NSV_PPEND ) .OR. &
#ifdef MNH_FOREFIRE
                      ( JSV >= NSV_FFBEG  .AND. JSV <= NSV_FFEND ) .OR. &
#endif
                      ( JSV >= NSV_CSBEG  .AND. JSV <= NSV_CSEND ) .OR. &
                      ( JSV >= NSV_SNWBEG .AND. JSV <= NSV_SNWEND .AND. GOLDFILEFORMAT ) ) THEN !Snow was not written in <5.6
              PLBYSVM(:,:,:,JSV) = 0.
              CALL PRINT_MSG( NVERB_WARNING, 'IO', 'INI_LB', 'PLBYSVM is initialized to 0 for ' // TRIM( YMNHNAME_BASE ) )
            ELSE
              CALL PRINT_MSG( NVERB_FATAL, 'IO', 'INI_LB', 'problem to initialize PLBYSVM for ' // TRIM( YMNHNAME_BASE ) )
            END IF
          END IF
        END IF
      END IF

    CASE( 'INIT' )
      IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
      IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
  END SELECT
END DO
!-------------------------------------------------------------------------------
!
!*       3.    COMPUTE THE LB SOURCES
!              -----------------------
!
! IN case of initialization of LB source terms (OLSOURCE=T) :
! xxxM are LB source terms 
! xxxMM are LB fields at time t -dt 
IF (OLSOURCE) THEN 
  IF (PRESENT(PLBXUMM).AND.PRESENT(PLBYUMM)) THEN
    PLBXUM(:,:,:) = (PLBXUM(:,:,:) - PLBXUMM(:,:,:))   / PLENG
    PLBYUM(:,:,:) = (PLBYUM(:,:,:) - PLBYUMM(:,:,:))   / PLENG
  ENDIF
  IF (PRESENT(PLBXVMM).AND.PRESENT(PLBYVMM)) THEN
    PLBXVM(:,:,:) = (PLBXVM(:,:,:) - PLBXVMM(:,:,:))   / PLENG
    PLBYVM(:,:,:) = (PLBYVM(:,:,:) - PLBYVMM(:,:,:))   / PLENG
  ENDIF
  IF (PRESENT(PLBXWMM).AND.PRESENT(PLBYWMM)) THEN 
    PLBXWM(:,:,:) = (PLBXWM(:,:,:) - PLBXWMM(:,:,:))   / PLENG
    PLBYWM(:,:,:) = (PLBYWM(:,:,:) - PLBYWMM(:,:,:))   / PLENG
  ENDIF
   IF (PRESENT(PLBXTHMM).AND.PRESENT(PLBYTHMM)) THEN 
    PLBXTHM(:,:,:) = (PLBXTHM(:,:,:) - PLBXTHMM(:,:,:))   / PLENG
    PLBYTHM(:,:,:) = (PLBYTHM(:,:,:) - PLBYTHMM(:,:,:))   / PLENG
  ENDIF
  IF (HGETTKEM =='READ') THEN
    IF (PRESENT(PLBXTKEMM).AND.PRESENT(PLBYTKEMM)) THEN 
      PLBXTKEM(:,:,:) = (PLBXTKEM(:,:,:) - PLBXTKEMM(:,:,:))   / PLENG
      PLBYTKEM(:,:,:) = (PLBYTKEM(:,:,:) - PLBYTKEMM(:,:,:))   / PLENG
    ENDIF
  ENDIF
  IF (HGETTKEM =='INIT') THEN
      PLBXTKEM(:,:,:) = 0.
      PLBYTKEM(:,:,:) = 0.
  ENDIF
! LB moist variables 
  IRR=0
  IF (PRESENT(PLBXRMM).AND.PRESENT(PLBYRMM))   THEN      
    DO JRR=1,7
      IF (YGETRXM(JRR) == 'READ') THEN      
        IRR=IRR+1  
        PLBXRM(:,:,:,IRR) = (PLBXRM(:,:,:,IRR) - PLBXRMM(:,:,:,IRR))   / PLENG
        PLBYRM(:,:,:,IRR) = (PLBYRM(:,:,:,IRR) - PLBYRMM(:,:,:,IRR))   / PLENG  
      ENDIF
    END DO
  ENDIF
! LB-scalar variables
  DO JSV=1,KSV
    IF (HGETSVM(JSV) == 'READ') THEN   
      PLBXSVM(:,:,:,JSV) = (PLBXSVM(:,:,:,JSV) - PLBXSVMM(:,:,:,JSV))   / PLENG
      PLBYSVM(:,:,:,JSV) = (PLBYSVM(:,:,:,JSV) - PLBYSVMM(:,:,:,JSV))   / PLENG 
    ENDIF
  END DO
!
ENDIF
!
CONTAINS

  SUBROUTINE OLD_CMNHNAME_GENERATE_INTERN( YMNHNAME, YLONGNAME )

    CHARACTER(LEN=*), INTENT(INOUT) :: YMNHNAME
    CHARACTER(LEN=*), INTENT(INOUT) :: YLONGNAME

    INTEGER :: IPOS
    INTEGER :: JI

    !Try to generate CMNHNAME with old format
    !In the old format, an indice of 2 numbers was written after the name but without trimming it
    IPOS = SCAN( YMNHNAME, '0123456789' )

    !Unmodified part YMNHNAME(1:IPOS-1) = YMNHNAME(1:IPOS-1)

    !Move number part at the new end
    IF ( 4+JPSVNAMELGTMAX+2 > LEN( YMNHNAME ) ) &
      CALL PRINT_MSG(NVERB_FATAL,'GEN','OLD_CMNHNAME_GENERATE_INTERN','CMNHNAME too small')
    YMNHNAME(4+JPSVNAMELGTMAX+1 : 4+JPSVNAMELGTMAX+2) = YMNHNAME(IPOS : IPOS+1)
    DO JI = IPOS, 4+JPSVNAMELGTMAX
      YMNHNAME(JI:JI) = ' '
    END DO

    YLONGNAME = TRIM( YMNHNAME )

  END SUBROUTINE OLD_CMNHNAME_GENERATE_INTERN

END SUBROUTINE INI_LB
