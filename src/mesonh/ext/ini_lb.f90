!MNH_LIC Copyright 1998-2021 CNRS, Meteo-France and Universite Paul Sabatier
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
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 13/02/2019: initialize PLBXSVM and PLBYSVM in all cases
!  P. Wautelet 14/02/2019: move UPCASE function to tools.f90
!!      Bielli S. 02/2019  Sea salt : significant sea wave height influences salt emission; 5 salt modes
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
USE MODD_CH_AEROSOL
USE MODD_CH_M9_n,         ONLY: CNAMES, CICNAMES
USE MODD_CTURB
USE MODD_CONF
USE MODD_DUST
USE MODD_ELEC_DESCR,      ONLY: CELECNAMES
use modd_field,           only: tfielddata, TYPELOG, TYPEREAL
USE MODD_ICE_C1R3_DESCR,  ONLY: C1R3NAMES
USE MODD_IO,              ONLY: TFILEDATA
USE MODD_LG,              ONLY: CLGNAMES
USE MODD_LUNIT_n,         ONLY: TLUOUT
USE MODD_NSV
USE MODD_PARAMETERS,      ONLY: JPHEXT,NMNHNAMELGTMAX
USE MODD_PARAM_LIMA
USE MODD_PARAM_LIMA_COLD, ONLY: CLIMA_COLD_NAMES
USE MODD_PARAM_LIMA_WARM, ONLY: CLIMA_WARM_NAMES
USE MODD_PARAM_n
USE MODD_RAIN_C2R2_DESCR, ONLY: C2R2NAMES
USE MODD_SALT
!
USE MODE_IO_FIELD_READ,   only: IO_Field_read, IO_Field_read_lb
USE MODE_MSG
USE MODE_TOOLS, ONLY: UPCASE
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
INTEGER             :: ILUOUT   !  Logical unit number associated with TLUOUT
LOGICAL :: GHORELAX_UVWTH  ! switch for the horizontal relaxation for U,V,W,TH in the FM file 
LOGICAL :: GHORELAX_TKE    ! switch for the horizontal relaxation for tke in the FM file
LOGICAL :: GHORELAX_R, GHORELAX_SV ! switch for the horizontal relaxation 
                                   ! for moist and scalar variables
CHARACTER (LEN= LEN(HGETRVM)), DIMENSION (7) :: YGETRXM ! Arrays with  the get indicators 
                                                        !  for the moist variables
CHARACTER (LEN=1), DIMENSION (7) :: YC    ! array with the prefix of the moist variables
CHARACTER(LEN=2)  :: INDICE ! to index CCN and IFN fields of LIMA scheme
TYPE(TFIELDDATA)  :: TZFIELD
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
ILUOUT = TLUOUT%NLU
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
IF (KSIZELBXU_ll/= 0) THEN
  CALL IO_Field_read_lb(TPINIFILE,'LBXUM',IL3DXU,IRIMXU,PLBXUM)
END IF

IF ( KSIZELBX_ll /= 0) THEN
  CALL IO_Field_read_lb(TPINIFILE,'LBXVM',IL3DX,IRIMX,PLBXVM)
ENDIF

IF ( KSIZELBX_ll /= 0) THEN
  CALL IO_Field_read_lb(TPINIFILE,'LBXWM',IL3DX,IRIMX,PLBXWM)
END IF

IF ( KSIZELBY_ll /= 0) THEN
  CALL IO_Field_read_lb(TPINIFILE,'LBYUM',IL3DY,IRIMY,PLBYUM)
END IF

IF ( KSIZELBYV_ll  /= 0) THEN
  CALL IO_Field_read_lb(TPINIFILE,'LBYVM',IL3DYV,IRIMYV,PLBYVM)
END IF

IF (KSIZELBY_ll /= 0) THEN
  CALL IO_Field_read_lb(TPINIFILE,'LBYWM',IL3DY,IRIMY,PLBYWM)
END IF

IF (KSIZELBX_ll /= 0) THEN
  CALL IO_Field_read_lb(TPINIFILE,'LBXTHM',IL3DX,IRIMX,PLBXTHM)
END IF

IF ( KSIZELBY_ll /= 0) THEN
  CALL IO_Field_read_lb(TPINIFILE,'LBYTHM',IL3DY,IRIMY,PLBYTHM)
END IF
!
!*        2.3  LB-TKE
!
SELECT CASE(HGETTKEM)                 
CASE('READ') 
  IF (.NOT. LCPL_AROME .AND. OLSOURCE) THEN
    IF (PRESENT(PLBXTKEMM).AND.PRESENT(PLBYTKEMM)) THEN
      WRITE ( ILUOUT,*) 'LBXTKES AND LBYTKES WILL BE INITIALIZED TO 0'
      PLBXTKEM(:,:,:) = PLBXTKEMM(:,:,:)    
      PLBYTKEM(:,:,:) = PLBYTKEMM(:,:,:)
    ELSE
!callabortstop
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
  TZFIELD%CMNHNAME   = 'HORELAX_R'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'HORELAX_R'
  TZFIELD%CUNITS     = ''
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = 'Switch to activate the HOrizontal RELAXation'
  TZFIELD%CLBTYPE    = 'NONE'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPELOG
  TZFIELD%NDIMS      = 0
  TZFIELD%LTIMEDEP   = .FALSE.
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
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
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
              WRITE(ILUOUT,*) 'PLBXRS  will be initialized to 0 for LBXR'//YC(JRR)//'M'
            ELSE
              !callabortstop
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
              WRITE(ILUOUT,*) 'PLBYRS  will be initialized to 0 for LBYR'//YC(JRR)//'M'
            ELSE
              !callabortstop
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
PLBXSVM(:,:,:,:) = 0.
PLBYSVM(:,:,:,:) = 0.
!
IF (KSV > 0) THEN
  IF (ANY(HGETSVM(1:KSV)=='READ')) THEN
    TZFIELD%CMNHNAME   = 'HORELAX_SV'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'HORELAX_SV'
    TZFIELD%CUNITS     = ''
    TZFIELD%CDIR       = '--'
    TZFIELD%CCOMMENT   = ''
    TZFIELD%CLBTYPE    = 'NONE'
    TZFIELD%NGRID      = 0
    TZFIELD%NTYPE      = TYPELOG
    TZFIELD%NDIMS      = 0
    TZFIELD%LTIMEDEP   = .FALSE.
    CALL IO_Field_read(TPINIFILE,TZFIELD,GHORELAX_SV)
    IF ( GHORELAX_SV ) THEN
      IRIMX=(KSIZELBXSV_ll-2*JPHEXT)/2   
      IRIMY=(KSIZELBYSV_ll-2*JPHEXT)/2
      IL3DX=2*ILBSIZEX+2*JPHEXT
      IL3DY=2*ILBSIZEY+2*JPHEXT
    ELSE
      IRIMX=0
      IRIMY=0
      IL3DX=2*JPHEXT !2
      IL3DY=2*JPHEXT !2
    END IF
  END IF
END IF
! User scalar variables
IF (NSV_USER>0) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = 1, NSV_USER
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          WRITE(TZFIELD%CMNHNAME,'(A6,I3.3)')'LBXSVM',JSV
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3,A8)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'PLXYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          WRITE(TZFIELD%CMNHNAME,'(A6,I3.3)')'LBYSVM',JSV
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3,A8)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! C2R2 scalar variables
IF (NSV_C2R2END>=NSV_C2R2BEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'm-3'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_C2R2BEG, NSV_C2R2END
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'C2R2 PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize C2R2 PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'C2R2 PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize C2R2 PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! C1R3 scalar variables
IF (NSV_C1R3END>=NSV_C1R3BEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'm-3'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_C1R3BEG, NSV_C1R3END
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'C1R3 PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize C1R3 PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'C1R3 PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize C1R3 PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
!
! LIMA: CCN and IFN scalar variables
!
IF (CCLOUD=='LIMA' ) THEN
  IF (NSV_LIMA_CCN_FREE+NMOD_CCN-1 >= NSV_LIMA_CCN_FREE) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'kg-1'
    TZFIELD%CDIR       = ''
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_LIMA_CCN_FREE,NSV_LIMA_CCN_FREE+NMOD_CCN-1
      SELECT CASE(HGETSVM(JSV))
        CASE ('READ')
          WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_CCN_FREE + 1)
          IF ( KSIZELBXSV_ll /= 0 ) THEN
            IF (        TPINIFILE%NMNHVERSION(1) < 5                                       &
                 .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) < 5 ) &
                 .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) == 5  &
                        .AND. TPINIFILE%NMNHVERSION(3) < 1 ) ) THEN
              TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CLIMA_WARM_NAMES(3)))//INDICE
            ELSE
              TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CLIMA_WARM_NAMES(3)))//INDICE
            END IF
            TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
            TZFIELD%CLBTYPE    = 'LBX'
            WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
            CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
            IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
              IF (IRESP/=0) THEN
                IF (PRESENT(PLBXSVMM)) THEN
                  PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                  WRITE(ILUOUT,*) 'CCN PLBXSVM   will be initialized to 0'
                ELSE
!callabortstop
                  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize CCN PLBXSVM')
                ENDIF
              END IF
            END IF
          END IF
          !
          IF (KSIZELBYSV_ll  /= 0 ) THEN
            IF (        TPINIFILE%NMNHVERSION(1) < 5                                       &
                 .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) < 5 ) &
                 .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) == 5  &
                        .AND. TPINIFILE%NMNHVERSION(3) < 1 ) ) THEN
              TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CLIMA_WARM_NAMES(3)))//INDICE
            ELSE
              TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CLIMA_WARM_NAMES(3)))//INDICE
            END IF
            TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
            TZFIELD%CLBTYPE    = 'LBY'
            WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
            CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
            IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
              IF (IRESP/=0) THEN
                IF (PRESENT(PLBYSVMM)) THEN
                  PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                  WRITE(ILUOUT,*) 'CCN PLBYSVM   will be initialized to 0'
                ELSE
!callabortstop
                  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize CCN PLBYSVM')
                ENDIF
              END IF
            END IF
          END IF
        CASE('INIT')
          IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
          IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
      END SELECT
    END DO
  END IF
  !
  IF (NSV_LIMA_IFN_FREE+NMOD_IFN-1 >= NSV_LIMA_IFN_FREE) THEN
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CUNITS     = 'kg-1'
    TZFIELD%CDIR       = ''
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    !
    DO JSV = NSV_LIMA_IFN_FREE,NSV_LIMA_IFN_FREE+NMOD_IFN-1
      SELECT CASE(HGETSVM(JSV))
        CASE ('READ')
          WRITE(INDICE,'(I2.2)')(JSV - NSV_LIMA_IFN_FREE + 1)
          IF ( KSIZELBXSV_ll /= 0 ) THEN
            IF (        TPINIFILE%NMNHVERSION(1) < 5                                       &
                 .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) < 5 ) &
                 .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) == 5  &
                        .AND. TPINIFILE%NMNHVERSION(3) < 1 ) ) THEN
              TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CLIMA_COLD_NAMES(5)))//INDICE
            ELSE
              TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CLIMA_COLD_NAMES(5)))//INDICE
            END IF
            TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
            TZFIELD%CLBTYPE    = 'LBX'
            WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
            CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
            IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
              IF (IRESP/=0) THEN
                IF (PRESENT(PLBXSVMM)) THEN
                  PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                  WRITE(ILUOUT,*) 'IFN PLBXSVM   will be initialized to 0'
                ELSE
!callabortstop
                  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize IFN')
                ENDIF
              END IF
            END IF
          END IF
          !
          IF (KSIZELBYSV_ll  /= 0 ) THEN
            IF (        TPINIFILE%NMNHVERSION(1) < 5                                       &
                 .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) < 5 ) &
                 .OR. ( TPINIFILE%NMNHVERSION(1) == 5 .AND. TPINIFILE%NMNHVERSION(2) == 5  &
                        .AND. TPINIFILE%NMNHVERSION(3) < 1 ) ) THEN
              TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CLIMA_COLD_NAMES(5)))//INDICE
            ELSE
              TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CLIMA_COLD_NAMES(5)))//INDICE
            END IF
            TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
            TZFIELD%CLBTYPE    = 'LBY'
            WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
            CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
            IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
              IF (IRESP/=0) THEN
                IF (PRESENT(PLBYSVMM)) THEN
                  PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                  WRITE(ILUOUT,*) 'IFN PLBYSVM   will be initialized to 0'
                ELSE
!callabortstop
                  CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize IFN')
                ENDIF
              END IF
            END IF
          END IF
        CASE('INIT')
          IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
          IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
      END SELECT
    END DO
  END IF
ENDIF
! ELEC scalar variables
IF (NSV_ELECEND>=NSV_ELECBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_ELECBEG, NSV_ELECEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'ELEC PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize ELEC PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'ELEC PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize ELEC PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Chemical gas phase scalar variables
IF (NSV_CHGSEND>=NSV_CHGSBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'kg-1'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_CHGSBEG, NSV_CHGSEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CNAMES(JSV-NSV_CHGSBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Chemical PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize gas phase chemical PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CNAMES(JSV-NSV_CHGSBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Chemical PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize gas phase chemical PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Chemical aqueous phase scalar variables
IF (NSV_CHACEND>=NSV_CHACBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppp'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_CHACBEG, NSV_CHACEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CNAMES(JSV-NSV_CHACBEG+NSV_CHGS+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Chemical PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize aqueous phase chemical PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CNAMES(JSV-NSV_CHACBEG+NSV_CHGS+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Chemical PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize aqueous phase chemical PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Chemical ice phase scalar variables
IF (NSV_CHICEND>=NSV_CHICBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppp'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_CHICBEG, NSV_CHICEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CICNAMES(JSV-NSV_CHICBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Ice phase chemical PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize ice phase chemical PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CICNAMES(JSV-NSV_CHICBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Ice phase chemical PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize ice phase chemical PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Orilam aerosol scalar variables
IF (NSV_AEREND>=NSV_AERBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppp'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_AERBEG, NSV_AEREND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CAERONAMES(JSV-NSV_AERBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Aerosol PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize aerosol PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CAERONAMES(JSV-NSV_AERBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Aerosol PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize aerosol PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Orilam aerosols moist scalar variables
IF (NSV_AERDEPEND>=NSV_AERDEPBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppp'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_AERDEPBEG, NSV_AERDEPEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CDEAERNAMES(JSV-NSV_AERDEPBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Aerosol PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize aerosol PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(CDEAERNAMES(JSV-NSV_AERDEPBEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Aerosol PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize aerosol PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Dust scalar variables
IF (NSV_DSTEND>=NSV_DSTBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppp'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_DSTBEG, NSV_DSTEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CDUSTNAMES(JSV-NSV_DSTBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Dust PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize dust PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CDUSTNAMES(JSV-NSV_DSTBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Dust PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize dust PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
!
IF (NSV_DSTDEPEND>=NSV_DSTDEPBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppp'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_DSTDEPBEG, NSV_DSTDEPEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CDEDSTNAMES(JSV-NSV_DSTDEPBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Dust Desposition PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize dust PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(UPCASE(CDEDSTNAMES(JSV-NSV_DSTDEPBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Dust Depoistion  PLBYSVM   will be initialized to 0'
              ELSE
                WRITE(ILUOUT,*) 'Pb to initialize dust PLBYSVM '
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize dust PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Sea salt scalar variables
IF (NSV_SLTEND>=NSV_SLTBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppp'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_SLTBEG, NSV_SLTEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(UPCASE(CSALTNAMES(JSV-NSV_SLTBEG+1)))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Sea Salt PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize sea salt PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(CSALTNAMES(JSV-NSV_SLTBEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Sea Salt PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize sea salt PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Passive pollutant variables
IF (NSV_PPEND>=NSV_PPBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_PPBEG, NSV_PPEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_PP'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Passive pollutant PLBXSVM   will be initialized to 0'
              ELSE
                PLBXSVM(:,:,:,JSV)=0.
                WRITE(ILUOUT,*) 'Passive pollutant PLBXSVM   will be initialized to 0'
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_PP'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Passive pollutant PLBYSVM   will be initialized to 0'
              ELSE
                PLBYSVM(:,:,:,JSV)=0.
                WRITE(ILUOUT,*) 'Passive pollutant PLBYSVM   will be initialized to 0'
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
#ifdef MNH_FOREFIRE
! ForeFire scalar variables
IF (NSV_FFEND>=NSV_FFBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_FFBEG, NSV_FFEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_FF'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          WRITE(ILUOUT,*) 'ForeFire LBX_FF ', IRESP
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'ForeFire pollutant PLBXSVM   will be initialized to 0'
              ELSE
                PLBXSVM(:,:,:,JSV)=0.
                WRITE(ILUOUT,*) 'ForeFire pollutant PLBXSVM   will be initialized to 0'
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_FF'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'ForeFire scalar variable PLBYSVM will be initialized to 0'
              ELSE
                PLBYSVM(:,:,:,JSV)=0.
                WRITE(ILUOUT,*) 'ForeFire scalar variable PLBYSVM will be initialized to 0'
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
#endif
! Conditional sampling variables
IF (NSV_CSEND>=NSV_CSBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_CSBEG, NSV_CSEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_CS'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Conditional sampling LBXSVM   will be initialized to 0'
              ELSE
                PLBXSVM(:,:,:,JSV)=0.
                WRITE(ILUOUT,*) 'Conditional sampling PLBXSVM   will be initialized to 0'
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_CS'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Conditional sampling PLBYSVM   will be initialized to 0'
              ELSE
                PLBYSVM(:,:,:,JSV)=0.
                WRITE(ILUOUT,*) 'Conditional sampling PLBYSVM   will be initialized to 0'
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Linox scalar variables
IF (NSV_LNOXEND>=NSV_LNOXBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'kg kg-1'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_LNOXBEG, NSV_LNOXEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_LINOX'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Linox PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize linox PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_LINOX'
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'Linox PLBYSVM   will be initialized to 0'
              ELSE
!calla bortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize linox PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
! Lagrangian variables
IF (NSV_LGEND>=NSV_LGBEG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = ''
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_LGBEG, NSV_LGEND
    SELECT CASE(HGETSVM(JSV))
      CASE ('READ')
        IF ( KSIZELBXSV_ll /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBX_'//TRIM(CLGNAMES(JSV-NSV_LGBEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBX'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'2_Y_Z_','LBXSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DX,IRIMX,PLBXSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBXSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBXSVMM)) THEN
                PLBXSVM(:,:,:,JSV)=PLBXSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'lagrangian PLBXSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize lagrangian PLBXSVM')
              ENDIF
            END IF
          END IF
        END IF
        !
        IF (KSIZELBYSV_ll  /= 0 ) THEN
          TZFIELD%CMNHNAME   = 'LBY_'//TRIM(CLGNAMES(JSV-NSV_LGBEG+1))
          TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
          TZFIELD%CLBTYPE    = 'LBY'
          WRITE(TZFIELD%CCOMMENT,'(A6,A6,I3.3)')'X_2_Z_','LBYSVM',JSV
          CALL IO_Field_read_lb(TPINIFILE,TZFIELD,IL3DY,IRIMY,PLBYSVM(:,:,:,JSV),IRESP)
          IF ( SIZE(PLBYSVM,1) /= 0 ) THEN
            IF (IRESP/=0) THEN
              IF (PRESENT(PLBYSVMM)) THEN
                PLBYSVM(:,:,:,JSV)=PLBYSVMM(:,:,:,JSV)
                WRITE(ILUOUT,*) 'lagrangian PLBYSVM   will be initialized to 0'
              ELSE
!callabortstop
                CALL PRINT_MSG(NVERB_FATAL,'GEN','INI_LB','problem to initialize lagrangian PLBYSVM')
              ENDIF
            END IF
          END IF
        END IF
      !
      CASE('INIT')
        IF ( SIZE(PLBXSVM,1) /= 0 ) PLBXSVM(:,:,:,JSV) = 0.
        IF ( SIZE(PLBYSVM,1) /= 0 ) PLBYSVM(:,:,:,JSV) = 0.
    END SELECT
  END DO
END IF
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
END SUBROUTINE INI_LB
