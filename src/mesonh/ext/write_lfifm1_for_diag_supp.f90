!MNH_LIC Copyright 2000-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################################
      MODULE MODI_WRITE_LFIFM1_FOR_DIAG_SUPP
!     ######################################
INTERFACE
!
   SUBROUTINE WRITE_LFIFM1_FOR_DIAG_SUPP(TPFILE)
!
USE MODD_IO, ONLY: TFILEDATA
!
!*       0.1   Declarations of arguments
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! Output file
!
END SUBROUTINE WRITE_LFIFM1_FOR_DIAG_SUPP
!
END INTERFACE
!
END MODULE MODI_WRITE_LFIFM1_FOR_DIAG_SUPP
!
!     ##############################################
      SUBROUTINE WRITE_LFIFM1_FOR_DIAG_SUPP(TPFILE)
!     ##############################################
!
!!****  *WRITE_LFIFM1_FOR_DIAG_SUPP* - write records in the diag file
!!
!!    PURPOSE
!!    -------
!        The purpose of this routine is to write in the file
!     of name YFMFILE//'.lfi' with the FM routines.  
!
!!**  METHOD
!!    ------
!!      The data are written in the LFIFM file :
!!        - diagnostics from the convection
!!        - diagnostics from the radiatif transfer code
!!
!!      The localization on the model grid is also indicated :
!!        IGRID = 1 for mass grid point
!!        IGRID = 2 for U grid point
!!        IGRID = 3 for V grid point
!!        IGRID = 4 for w grid point
!!        IGRID = 0 for meaningless case
!!
!!    EXTERNAL
!!    --------
!!      FMWRIT : FM-routine to write a record
!!
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!  	J. Stein   *Meteo France* 
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    13/09/00
!!      N. Asencio  15/09/00 computation of temperature and height of clouds is moved
!!                           here and deleted in WRITE_LFIFM1_FOR_DIAG routine
!!      I. Mallet   02/11/00 add the call to RADTR_SATEL
!!      J.-P. Chaboureau 11/12/03 add call the CALL_RTTOV (table NRTTOVINFO to
!!              choose the platform, the satellite, the sensor for all channels 
!!              (see the table in rttov science and validation report) and the
!!              type of calculations in the namelist: 0 = tb, 1 = tb + jacobian,
!!              2 = tb + adjoint, 3 = tb + jacobian + adjoint)
!!      V. Masson   01/2004  removes surface (externalization)
!!      October 2009 (G. Tanguy) add ILENCH=LEN(YCOMMENT) after
!!                                              change of YCOMMENT
!!      October 2011 (C.Lac) FF10MAX  : interpolation of 10m wind
!!        between 2 Meso-NH levels if 10m is above the first atmospheric level
!!      2015 : D.Ricard add UM10/VM10 for LCARTESIAN=T cases
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      P.Tulet : Diag for salt and orilam
!!      J.-P. Chaboureau 07/03/2016 fix the dimensions of local arrays
!!      P.Wautelet : 11/07/2016 : removed MNH_NCWRIT define
!!      J.-P. Chaboureau 31/10/2016 add the call to RTTOV11
!!      F. Brosse 10/2016 add chemical production destruction terms outputs
!!      M.Leriche 01/07/2017 Add DIAG chimical surface fluxes
!!      J.-P. Chaboureau 01/2018 add altitude interpolation
!!      J.-P. Chaboureau 01/2018 add coarse graining
!!  Philippe Wautelet: 05/2016-04/2018: new data structures and calls for I/O
!!      J.-P. Chaboureau 07/2018 bug fix on XEMIS when calling CALL_RTTOVxx
!!      J.-P. Chaboureau 09/04/2021 add the call to RTTOV13
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODE_ll
USE MODD_CST
use modd_field,           only: tfielddata, tfieldlist, TYPEINT, TYPEREAL
USE MODD_IO, ONLY: TFILEDATA
USE MODD_PARAMETERS
USE MODD_CONF_n
USE MODD_CONF
USE MODD_DEEP_CONVECTION_n
USE MODD_DIM_n
USE MODD_FIELD_n
USE MODD_GRID_n
USE MODD_LUNIT_n
USE MODD_PARAM_n
USE MODD_PARAM_KAFR_n
USE MODD_PARAM_RAD_n
USE MODD_RADIATIONS_n
USE MODD_TIME_n
USE MODD_TURB_n
USE MODD_REF_n, ONLY: XRHODREF
USE MODD_DIAG_FLAG
USE MODD_NSV, ONLY : NSV,NSV_USER,NSV_C2R2BEG,NSV_C2R2END,             &
                     NSV_C1R3BEG, NSV_C1R3END,NSV_ELECBEG,NSV_ELECEND, &
                     NSV_CHEMBEG, NSV_CHEMEND,NSV_LGBEG,  NSV_LGEND
USE MODD_CH_M9_n,         ONLY: CNAMES
USE MODD_RAIN_C2R2_DESCR, ONLY: C2R2NAMES
USE MODD_ICE_C1R3_DESCR,  ONLY: C1R3NAMES
USE MODD_ELEC_DESCR,      ONLY: CELECNAMES
USE MODD_LG,              ONLY: CLGNAMES
USE MODD_DUST,            ONLY: LDUST
USE MODD_SALT,            ONLY: LSALT
USE MODD_CH_AEROSOL,      ONLY: LORILAM
USE MODD_CH_MNHC_n
USE MODD_CH_BUDGET_n
USE MODD_CH_PRODLOSSTOT_n
USE MODD_CH_FLX_n,          ONLY: XCHFLX
USE MODD_RAD_TRANSF
USE MODD_DIAG_IN_RUN, ONLY: XCURRENT_ZON10M,XCURRENT_MER10M,           &
                            XCURRENT_SFCO2,XCURRENT_SWD, XCURRENT_LWD, &
                            XCURRENT_SWU, XCURRENT_LWU
!
USE MODD_DYN_n
USE MODD_CURVCOR_n
USE MODD_METRICS_n
USE MODD_DIAG_BLANK
USE MODI_PINTER
USE MODI_ZINTER
USE MODI_GRADIENT_M
USE MODI_GRADIENT_W
USE MODI_GRADIENT_U
USE MODI_GRADIENT_V
USE MODI_GRADIENT_UV
!
USE MODI_SHUMAN
USE MODE_NEIGHBORAVG
#ifdef MNH_RTTOV_8
USE MODI_CALL_RTTOV8
#endif
#ifdef MNH_RTTOV_11
USE MODI_CALL_RTTOV11
#endif
#ifdef MNH_RTTOV_13
USE MODI_CALL_RTTOV13
#endif
USE MODI_RADTR_SATEL
USE MODI_UV_TO_ZONAL_AND_MERID
!
use mode_field,          only: Find_field_id_from_mnhname
USE MODE_IO_FIELD_WRITE, only: IO_Field_write
!
USE MODI_GET_SURF_UNDEF
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments
!
TYPE(TFILEDATA),   INTENT(IN) :: TPFILE ! Output file
!
!*       0.2   Declarations of local variables
!
INTEGER           :: IIU,IJU,IKU,IIB,IJB,IKB,IIE,IJE,IKE ! Arrays bounds
INTEGER           :: IKRAD  
! 
INTEGER           :: JI,JJ,JK,JSV   ! loop index
! 
! variables for Diagnostic variables related to deep convection
REAL,DIMENSION(:,:), ALLOCATABLE              :: ZWORK21,ZWORK22
!
! variables for computation of temperature and height of clouds
REAL :: ZCLMR ! value of mixing ratio tendency  for detection of cloud top
LOGICAL, DIMENSION(:,:), ALLOCATABLE          :: GMASK2 
INTEGER, DIMENSION(:,:), ALLOCATABLE          :: IWORK1, IWORK2
INTEGER, DIMENSION(:,:), ALLOCATABLE          :: ICL_HE_ST
REAL,    DIMENSION(:,:,:), ALLOCATABLE        :: ZWORK31,ZTEMP
!
! variables needed for the transfer radiatif diagnostic code
INTEGER :: ITOTGEO
INTEGER, DIMENSION (JPGEOST) :: INDGEO
CHARACTER(LEN=8), DIMENSION (JPGEOST) :: YNAM_SAT
REAL, DIMENSION(:,:), ALLOCATABLE :: ZIRBT, ZWVBT
REAL  :: ZUNDEF ! undefined value in SURFEX
!
! variables needed for 10m wind                                 
INTEGER :: ILEVEL
!
INTEGER :: IPRES, ITH
CHARACTER(LEN=4) :: YCAR4
CHARACTER(LEN=4), DIMENSION(SIZE(XISOPR)) :: YPRES
CHARACTER(LEN=4), DIMENSION(SIZE(XISOTH)) :: YTH
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWORK32,ZWORK33,ZWORK34,ZWRES,ZPRES,ZWTH
REAL, DIMENSION(:), ALLOCATABLE :: ZTH
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZPOVO
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZVOX,ZVOY,ZVOZ
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3))  :: ZCORIOZ
TYPE(TFIELDDATA)              :: TZFIELD
TYPE(TFIELDDATA),DIMENSION(2) :: TZFIELD2
!
! variables needed for altitude interpolation                                 
INTEGER :: IAL
REAL :: ZFILLVAL
REAL, DIMENSION(:), ALLOCATABLE :: ZAL
REAL, DIMENSION(:,:,:), ALLOCATABLE :: ZWAL
!
! variables needed for coarse graining
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ZUT_PRM,ZVT_PRM,ZWT_PRM
REAL,DIMENSION(SIZE(XTHT,1),SIZE(XTHT,2),SIZE(XTHT,3)) :: ZUU_AVG,ZVV_AVG,ZWW_AVG
INTEGER :: IDX, IID, IRESP
CHARACTER(LEN=3) :: YDX
!-------------------------------------------------------------------------------
!
!*       0.     ARRAYS BOUNDS INITIALIZATION
!
IIU=SIZE(XTHT,1)
IJU=SIZE(XTHT,2)
IKU=SIZE(XTHT,3)
CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
IKB=1+JPVEXT
IKE=IKU-JPVEXT
!
ALLOCATE(ZWORK21(IIU,IJU))
ALLOCATE(ZWORK31(IIU,IJU,IKU))
ALLOCATE(ZTEMP(IIU,IJU,IKU))
ZTEMP(:,:,:)=XTHT(:,:,:)*(XPABST(:,:,:)/ XP00) **(XRD/XCPD)
!
!-------------------------------------------------------------------------------
!
!*       1.     DIAGNOSTIC RELATED TO CONVECTION
!               -------------------------------- 
!
!* Diagnostic variables related to deep convection
!
IF (NCONV_KF >= 0) THEN
!
  CALL IO_Field_write(TPFILE,'CAPE',XCAPE)
!
  ! top height (km) of convective clouds
  ZWORK21(:,:)= 0.
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IF (NCLTOPCONV(JI,JJ)/=0) ZWORK21(JI,JJ)= XZZ(JI,JJ,NCLTOPCONV(JI,JJ))/1.E3
    END DO
  END DO
  TZFIELD%CMNHNAME   = 'CLTOPCONV'
  TZFIELD%CSTDNAME   = 'convective_cloud_top_altitude'
  TZFIELD%CLONGNAME  = 'CLTOPCONV'
  TZFIELD%CUNITS     = 'km'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Top of Convective Cloud'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
!
  ! base height (km) of convective clouds
  ZWORK21(:,:)= 0.
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IF (NCLBASCONV(JI,JJ)/=0) ZWORK21(JI,JJ)= XZZ(JI,JJ,NCLBASCONV(JI,JJ))/1.E3
    END DO
  END DO
  TZFIELD%CMNHNAME   = 'CLBASCONV'
  TZFIELD%CSTDNAME   = 'convective_cloud_base_altitude'
  TZFIELD%CLONGNAME  = 'CLBASCONV'
  TZFIELD%CUNITS     = 'km'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Base of Convective Cloud'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
!
END IF
IF (NCONV_KF >= 1) THEN
!
  CALL IO_Field_write(TPFILE,'DTHCONV',XDTHCONV)
  CALL IO_Field_write(TPFILE,'DRVCONV',XDRVCONV)
  CALL IO_Field_write(TPFILE,'DRCCONV',XDRCCONV)
  CALL IO_Field_write(TPFILE,'DRICONV',XDRICONV)
!  
  IF ( LCHTRANS .AND. NSV > 0 ) THEN
    ! User scalar variables
    IF (NSV_USER>0) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 's-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = 1, NSV_USER
        WRITE(TZFIELD%CMNHNAME,'(A7,I3.3)')'DSVCONV',JSV
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        WRITE(TZFIELD%CCOMMENT,'(A6,A2,I3.3,A20)')'X_Y_Z_','SV',JSV,' CONVective tendency'
        CALL IO_Field_write(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
    END IF
    ! microphysical C2R2 scheme scalar variables
    IF (NSV_C2R2END>=NSV_C2R2BEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 's-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_C2R2BEG, NSV_C2R2END
        TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(C2R2NAMES(JSV-NSV_C2R2BEG+1))//' CONVective tendency'
        CALL IO_Field_write(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
    END IF
    ! microphysical C3R5 scheme additional scalar variables
    IF (NSV_C1R3END>=NSV_C1R3BEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 's-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_C1R3BEG,NSV_C1R3END
        TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(C1R3NAMES(JSV-NSV_C1R3BEG+1))//' CONVective tendency'
        CALL IO_Field_write(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
    END IF
    ! electrical scalar variables
    IF (NSV_ELECEND>=NSV_ELECBEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 's-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_ELECBEG,NSV_ELECEND
        TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(CELECNAMES(JSV-NSV_ELECBEG+1))//' CONVective tendency'
        CALL IO_Field_write(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
    END IF
    ! chemical scalar variables
    IF (NSV_CHEMEND>=NSV_CHEMBEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 's-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_CHEMBEG, NSV_CHEMEND
        TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CNAMES(JSV-NSV_CHEMBEG+1))
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(CNAMES(JSV-NSV_CHEMBEG+1))//' CONVective tendency'
        CALL IO_Field_write(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
    END IF
    ! lagrangian variables
    IF (NSV_LGEND>=NSV_LGBEG) THEN
      TZFIELD%CSTDNAME   = ''
      TZFIELD%CUNITS     = 's-1'
      TZFIELD%CDIR       = 'XY'
      TZFIELD%NGRID      = 1
      TZFIELD%NTYPE      = TYPEREAL
      TZFIELD%NDIMS      = 3
      TZFIELD%LTIMEDEP   = .TRUE.
      !
      DO JSV = NSV_LGBEG,NSV_LGEND
        TZFIELD%CMNHNAME   = 'DSVCONV_'//TRIM(CLGNAMES(JSV-NSV_LGBEG+1))
        TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
        TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(CLGNAMES(JSV-NSV_LGBEG+1))//' CONVective tendency'
        CALL IO_Field_write(TPFILE,TZFIELD,XDSVCONV(:,:,:,JSV))
      END DO
    END IF
  END IF
!
END IF
IF (NCONV_KF >= 2) THEN
  CALL IO_Field_write(TPFILE,'PRLFLXCONV',XPRLFLXCONV)
  CALL IO_Field_write(TPFILE,'PRSFLXCONV',XPRSFLXCONV)
  CALL IO_Field_write(TPFILE,'UMFCONV',   XUMFCONV)
  CALL IO_Field_write(TPFILE,'DMFCONV',   XDMFCONV)
END IF
!-------------------------------------------------------------------------------
!
!* Height and temperature of clouds top
!
IF (LCLD_COV .AND. LUSERC) THEN
  ALLOCATE(IWORK1(IIU,IJU),IWORK2(IIU,IJU))
  ALLOCATE(ICL_HE_ST(IIU,IJU))
  ALLOCATE(GMASK2(IIU,IJU))
  ALLOCATE(ZWORK22(IIU,IJU))
!
! Explicit clouds
!
  ICL_HE_ST(:,:)=IKB  !initialization
  IWORK1(:,:)=IKB     ! with the
  IWORK2(:,:)=IKB     ! ground values
  ZCLMR=1.E-4         ! detection of clouds for cloud mixing ratio > .1g/kg
!
  GMASK2(:,:)=.TRUE.
  ZWORK31(:,:,:)= MZM( XRT(:,:,:,2) ) ! cloud mixing ratio at zz levels
  DO JK=IKE,IKB,-1
    WHERE ( (GMASK2(:,:)).AND.(ZWORK31(:,:,JK)>ZCLMR) )
      GMASK2(:,:)=.FALSE.
      IWORK1(:,:)=JK
    END WHERE
  END DO
!
  IF (LUSERI) THEN
    GMASK2(:,:)=.TRUE.
    ZWORK31(:,:,:)= MZM( XRT(:,:,:,4) ) ! cloud mixing ratio at zz levels
    DO JK=IKE,IKB,-1
      WHERE ( (GMASK2(:,:)).AND.(ZWORK31(:,:,JK)>ZCLMR) )
        GMASK2(:,:)=.FALSE.
        IWORK2(:,:)=JK
      END WHERE
    END DO
  END IF
!
  ZWORK21(:,:)=0.
  DO JJ=IJB,IJE
   DO JI=IIB,IIE
     ICL_HE_ST(JI,JJ)=MAX(IWORK1(JI,JJ),IWORK2(JI,JJ) )
     ZWORK21(JI,JJ)  =XZZ(JI,JJ,ICL_HE_ST(JI,JJ)) ! height (m) of explicit clouds
   END DO
  END DO 
!
  WHERE ( ZWORK21(:,:)==XZZ(:,:,IKB) ) ZWORK21=0. ! set the height to 
                                                  ! 0 if there is no cloud
  ZWORK21(:,:)=ZWORK21(:,:)/1.E3            ! height (km) of explicit clouds
!
  TZFIELD%CMNHNAME   = 'HECL'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'HECL'
  TZFIELD%CUNITS     = 'km'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Height of Explicit CLoud top'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
!
!  Higher top of the different species of clouds
!
  IWORK1(:,:)=IKB  ! initialization with the ground values
  ZWORK31(:,:,:)=MZM(ZTEMP(:,:,:)) ! temperature (K) at zz levels
  IF(CRAD/='NONE')  ZWORK31(:,:,IKB)=XTSRAD(:,:)
  ZWORK21(:,:)=0.
  ZWORK22(:,:)=0.
  DO JJ=IJB,IJE
    DO JI=IIB,IIE
      IWORK1(JI,JJ)=ICL_HE_ST(JI,JJ)
      IF (NCONV_KF >=0) &
      IWORK1(JI,JJ)= MAX(ICL_HE_ST(JI,JJ),NCLTOPCONV(JI,JJ))
      ZWORK21(JI,JJ)= XZZ(JI,JJ,IWORK1(JI,JJ))         ! max. cloud height (m)
      ZWORK22(JI,JJ)= ZWORK31(JI,JJ,IWORK1(JI,JJ))-XTT ! cloud temperature (C)
    END DO
  END DO 
!
  IF (NCONV_KF <0) THEN
    PRINT*,'YOU DO NOT ASK FOR CONVECTIVE DIAGNOSTICS (NCONV_KF<0), SO'
    PRINT*,'  HC not written in FM-file (equal to HEC)'
  ELSE
    WHERE ( ZWORK21(:,:)==XZZ(:,:,IKB) ) ZWORK21(:,:)=0. ! set the height to 
                                                         ! 0 if there is no cloud
    ZWORK21(:,:)=ZWORK21(:,:)/1.E3                 ! max. cloud height (km)
!
    TZFIELD%CMNHNAME   = 'HCL'
    TZFIELD%CSTDNAME   = 'cloud_top_altitude'
    TZFIELD%CLONGNAME  = 'HCL'
    TZFIELD%CUNITS     = 'km'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Height of CLoud top'
    TZFIELD%NGRID      = 4
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  ENDIF
!
  TZFIELD%CMNHNAME   = 'TCL'
  TZFIELD%CSTDNAME   = 'air_temperature_at_cloud_top'
  TZFIELD%CLONGNAME  = 'TCL'
  TZFIELD%CUNITS     = 'celsius'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Height of CLoud top'
  TZFIELD%NGRID      = 4
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK22)
!
  CALL IO_Field_write(TPFILE,'CLDFR',XCLDFR)
  CALL IO_Field_write(TPFILE,'ICEFR',XICEFR)
!
!  Visibility                                    
!
  ZWORK31(:,:,:)= 1.E4                ! 10 km for clear sky
  WHERE (XRT(:,:,:,2) > 0.)
    ZWORK31(:,:,:)=3.9E3/(144.7*(XRHODREF(:,:,:)*1.E3*XRT(:,:,:,2)/(1.+XRT(:,:,:,2)))**0.88)
  END WHERE
!
  TZFIELD%CMNHNAME   = 'VISI_HOR'
  TZFIELD%CSTDNAME   = 'visibility_in_air'
  TZFIELD%CLONGNAME  = 'VISI_HOR'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Z_VISI_HOR'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!
  DEALLOCATE(IWORK1,IWORK2,ICL_HE_ST,GMASK2,ZWORK22)
END IF
!
!-------------------------------------------------------------------------------
!
!*       2.    DIAGNOSTIC RELATED TO RADIATIONS
!              --------------------------------
!
IF (NRAD_3D >= 0) THEN
  IF (CRAD /= 'NONE') THEN
    CALL IO_Field_write(TPFILE,'DTHRAD',      XDTHRAD)
    CALL IO_Field_write(TPFILE,'FLALWD',      XFLALWD)
    CALL IO_Field_write(TPFILE,'DIRFLASWD',   XDIRFLASWD)
    CALL IO_Field_write(TPFILE,'SCAFLASWD',   XSCAFLASWD)
    CALL IO_Field_write(TPFILE,'DIRSRFSWD',   XDIRSRFSWD)
    CALL IO_Field_write(TPFILE,'CLEARCOL_TM1',NCLEARCOL_TM1)
    CALL IO_Field_write(TPFILE,'ZENITH',      XZENITH)
    CALL IO_Field_write(TPFILE,'AZIM',        XAZIM)
    CALL IO_Field_write(TPFILE,'DIR_ALB',     XDIR_ALB)
    CALL IO_Field_write(TPFILE,'SCA_ALB',     XSCA_ALB)
    !
    CALL PRINT_MSG(NVERB_INFO,'IO','WRITE_LFIFM1_FOR_DIAG_SUPP','EMIS: writing only first band')
    CALL FIND_FIELD_ID_FROM_MNHNAME('EMIS',IID,IRESP)
    TZFIELD = TFIELDLIST(IID)
    TZFIELD%NDIMS = 2
    CALL IO_Field_write(TPFILE,TZFIELD,XEMIS(:,:,1))
    !
    CALL IO_Field_write(TPFILE,'TSRAD',       XTSRAD)
  ELSE
    PRINT*,'YOU WANT DIAGNOSTICS RELATED TO RADIATION'
    PRINT*,' BUT NO RADIATIVE SCHEME WAS ACTIVATED IN THE MODEL'
  END IF
END IF
IF (NRAD_3D >= 1) THEN
  IF (LDUST) THEN
!Dust optical depth between two vertical levels
    ZWORK31(:,:,:)=0.
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      ZWORK31(:,:,JK)= XAER(:,:,IKRAD,3)
    END DO
    TZFIELD%CMNHNAME   = 'DSTAOD3D'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'DSTAOD3D'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_DuST Aerosol Optical Depth'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!Dust optical depth
    ZWORK21(:,:)=0.0
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZWORK21(JI,JJ)=ZWORK21(JI,JJ)+XAER(JI,JJ,IKRAD,3)
        ENDDO
      ENDDO
    ENDDO
    TZFIELD%CMNHNAME   = 'DSTAOD2D'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'DSTAOD2D'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_DuST Aerosol Optical Depth'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
!Dust extinction (optical depth per km)
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      ZWORK31(:,:,JK)= XAER(:,:,IKRAD,3)/(XZZ(:,:,JK+1)-XZZ(:,:,JK))*1.D3
    ENDDO
    TZFIELD%CMNHNAME   = 'DSTEXT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'DSTEXT'
    TZFIELD%CUNITS     = 'km-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_DuST EXTinction'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END IF
  IF (LSALT) THEN
!Salt optical depth between two vertical levels
    ZWORK31(:,:,:)=0.
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      ZWORK31(:,:,JK)= XAER(:,:,IKRAD,2)
    END DO
    TZFIELD%CMNHNAME   = 'SLTAOD3D'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'SLTAOD3D'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Salt Aerosol Optical Depth'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!Salt optical depth
    ZWORK21(:,:)=0.0
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZWORK21(JI,JJ)=ZWORK21(JI,JJ)+XAER(JI,JJ,IKRAD,2)
        ENDDO
      ENDDO
    ENDDO
    TZFIELD%CMNHNAME   = 'SLTAOD2D'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'SLTAOD2D'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Salt Aerosol Optical Depth'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
!Salt extinction (optical depth per km)
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      ZWORK31(:,:,JK)= XAER(:,:,IKRAD,2)/(XZZ(:,:,JK+1)-XZZ(:,:,JK))*1.D3
    ENDDO
    TZFIELD%CMNHNAME   = 'SLTEXT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'SLTEXT'
    TZFIELD%CUNITS     = 'km-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Salt EXTinction'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END IF
  IF (LORILAM) THEN
!Orilam anthropogenic optical depth between two vertical levels
    ZWORK31(:,:,:)=0.
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      ZWORK31(:,:,JK)= XAER(:,:,IKRAD,4)
    END DO
    TZFIELD%CMNHNAME   = 'AERAOD3D'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'AERAOD3D'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Anthropogenic Aerosol Optical Depth'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!Orilam anthropogenic optical depth
    ZWORK21(:,:)=0.0
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      DO JJ=IJB,IJE
        DO JI=IIB,IIE
          ZWORK21(JI,JJ)=ZWORK21(JI,JJ)+XAER(JI,JJ,IKRAD,4)
        ENDDO
      ENDDO
    ENDDO
    TZFIELD%CMNHNAME   = 'AERAOD2D'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'AERAOD2D'
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Anthropogenic Aerosol Optical Depth'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
!Orilam anthropogenic extinction (optical depth per km)
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      ZWORK31(:,:,JK)= XAER(:,:,IKRAD,4)/(XZZ(:,:,JK+1)-XZZ(:,:,JK))*1.D3
    ENDDO
    TZFIELD%CMNHNAME   = 'AEREXT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'AEREXT'
    TZFIELD%CUNITS     = 'km-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_Z_Anthropogenic EXTinction'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
  END IF
END IF
!
!-------------------------------------------------------------------------------
! Net surface gaseous fluxes
!print*,'LCHEMDIAG, NSV_CHEMBEG, NSV_CHEMEND=',&
!LCHEMDIAG, NSV_CHEMBEG, NSV_CHEMEND

IF (LCHEMDIAG) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppb m s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = NSV_CHEMBEG, NSV_CHEMEND
    TZFIELD%CMNHNAME   = 'FLX_'//TRIM(CNAMES(JSV-NSV_CHEMBEG+1))
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    WRITE(TZFIELD%CCOMMENT,'(A6,A,A)')'X_Y_Z_',TRIM(CNAMES(JSV-NSV_CHEMBEG+1)),' Net chemical flux'
    CALL IO_Field_write(TPFILE,TZFIELD,XCHFLX(:,:,JSV-NSV_CHEMBEG+1) * 1E9)
  END DO
END IF
!-------------------------------------------------------------------------------
!
!* Brightness temperatures from the radiatif transfer code (Morcrette, 1991)
!
IF (LEN_TRIM(CRAD_SAT) /= 0 .AND. NRR /=0) THEN
  ALLOCATE (ZIRBT(IIU,IJU),ZWVBT(IIU,IJU))
  ITOTGEO=0
  IF (INDEX(CRAD_SAT,'GOES-E')   /= 0) THEN
    ITOTGEO= ITOTGEO+1
    INDGEO(ITOTGEO) = 1
    YNAM_SAT(ITOTGEO) = 'GOES-E'
  END IF
  IF (INDEX(CRAD_SAT,'GOES-W')   /= 0) THEN
    ITOTGEO= ITOTGEO+1
    INDGEO(ITOTGEO) = 2
    YNAM_SAT(ITOTGEO) = 'GOES-W'
  END IF
  IF (INDEX(CRAD_SAT,'GMS')      /= 0) THEN
    ITOTGEO= ITOTGEO+1
    INDGEO(ITOTGEO) = 3
    YNAM_SAT(ITOTGEO) = 'GMS'
  END IF
  IF (INDEX(CRAD_SAT,'INDSAT')   /= 0) THEN
    ITOTGEO= ITOTGEO+1
    INDGEO(ITOTGEO) = 4
    YNAM_SAT(ITOTGEO) = 'INDSAT'
  END IF
  IF (INDEX(CRAD_SAT,'METEOSAT') /= 0) THEN
    ITOTGEO= ITOTGEO+1
    INDGEO(ITOTGEO) = 5
    YNAM_SAT(ITOTGEO) = 'METEOSAT'
  END IF
  PRINT*,'YOU ASK FOR BRIGHTNESS TEMPERATURES FOR ',ITOTGEO,' SATELLITE(S)'
  IF (NRR==1) THEN
    PRINT*,' THERE IS ONLY VAPOR WATER IN YOUR ATMOSPHERE'
    PRINT*,' IRBT WILL NOT TAKE INTO ACCOUNT CLOUDS.'
  END IF
  !
  DO JI=1,ITOTGEO
    ZIRBT(:,:) = XUNDEF
    ZWVBT(:,:) = XUNDEF
    CALL RADTR_SATEL( TDTCUR%nyear, TDTCUR%nmonth, TDTCUR%nday, TDTCUR%xtime, &
                      NDLON, NFLEV, NSTATM, NRAD_COLNBR, XEMIS(:,:,1),        &
                      XCCO2, XTSRAD, XSTATM, XTHT, XRT, XPABST, XZZ,          &
                      XSIGS, XMFCONV, MAX(XCLDFR,XICEFR), LUSERI, LSIGMAS,    &
                      LSUBG_COND, LRAD_SUBG_COND, ZIRBT, ZWVBT,               &
                      INDGEO(JI), VSIGQSAT                                    )
    !
    TZFIELD%CMNHNAME   = TRIM(YNAM_SAT(JI))//'_IRBT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = TRIM(YNAM_SAT(JI))//' Infra-Red Brightness Temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZIRBT)
    !
    TZFIELD%CMNHNAME   = TRIM(YNAM_SAT(JI))//'_WVBT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = TRIM(YNAM_SAT(JI))//' Water-Vapor Brightness Temperature'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWVBT)
  END DO
  DEALLOCATE(ZIRBT,ZWVBT)
END IF
!
!-------------------------------------------------------------------------------
!
!* Brightness temperatures from the Radiatif Transfer for Tiros Operational
! Vertical Sounder (RTTOV) code
!
IF (NRTTOVINFO(1,1) /= NUNDEF) THEN
! PRINT*,'YOU ASK FOR BRIGHTNESS TEMPERATURE COMPUTED BY THE RTTOV CODE'
#if defined(MNH_RTTOV_8)
  CALL CALL_RTTOV8(NDLON, NFLEV, NSTATM, XEMIS(:,:,1), XTSRAD, XSTATM, XTHT, XRT, &
                  XPABST, XZZ, XMFCONV, MAX(XCLDFR,XICEFR), XUT(:,:,IKB), XVT(:,:,IKB),   &
                  LUSERI, NRTTOVINFO, TPFILE                                  )
#elif defined(MNH_RTTOV_11)
  CALL CALL_RTTOV11(NDLON, NFLEV, XEMIS(:,:,1), XTSRAD, XTHT, XRT,            &
                  XPABST, XZZ, XMFCONV, MAX(XCLDFR,XICEFR), XUT(:,:,IKB), XVT(:,:,IKB),   &
                  LUSERI, NRTTOVINFO, TPFILE                                  )
#elif defined(MNH_RTTOV_13)
  CALL CALL_RTTOV13(NDLON, NFLEV, XEMIS(:,:,1), XTSRAD, XTHT, XRT,            &
                  XPABST, XZZ, XMFCONV, MAX(XCLDFR,XICEFR), XUT(:,:,IKB), XVT(:,:,IKB),   &
                  LUSERI, NRTTOVINFO, TPFILE                                  )
#else
PRINT *, "RTTOV LIBRARY NOT AVAILABLE = ###CALL_RTTOV####"
#endif
END IF
!
!-------------------------------------------------------------------------------
!
!*       3.    DIAGNOSTIC RELATED TO SURFACE
!              -----------------------------
!
IF (CSURF=='EXTE') THEN
!! Since SURFEX7 (masdev49) XCURRENT_ZON10M and XCURRENT_MER10M
!! are equal to XUNDEF of SURFEX if the first atmospheric level
!! is under 10m
  CALL GET_SURF_UNDEF(ZUNDEF)
!
  ILEVEL=IKB 
  !While there are XUNDEF values and we aren't at model's top
  DO WHILE(ANY(XCURRENT_ZON10M(IIB:IIE,IJB:IJE)==ZUNDEF) .AND. (ILEVEL/=IKE-1) )

    !Where interpolation is needed and possible
    !(10m is between ILEVEL and ILEVEL+1 or 10m is below the bottom level)
    WHERE(XCURRENT_ZON10M(IIB:IIE,IJB:IJE)==ZUNDEF .AND. &
                   ( XZHAT(ILEVEL+1) + XZHAT(ILEVEL+2)) /2. >10.)

      !Interpolation between ILEVEL and ILEVEL+1
      XCURRENT_ZON10M(IIB:IIE,IJB:IJE)=XUT(IIB:IIE,IJB:IJE,ILEVEL) + &
            (XUT(IIB:IIE,IJB:IJE,ILEVEL+1)-XUT(IIB:IIE,IJB:IJE,ILEVEL)) * &
            ( 10.- (XZHAT(ILEVEL)+XZHAT(ILEVEL+1))/2. ) / &
           ( (XZHAT(ILEVEL+2)-XZHAT(ILEVEL)) /2.)
      XCURRENT_MER10M(IIB:IIE,IJB:IJE)=XVT(IIB:IIE,IJB:IJE,ILEVEL) + &
            (XVT(IIB:IIE,IJB:IJE,ILEVEL+1)-XVT(IIB:IIE,IJB:IJE,ILEVEL)) * &
            (10.- (XZHAT(ILEVEL)+XZHAT(ILEVEL+1))/2. ) / &                                    
           ( (XZHAT(ILEVEL+2)-XZHAT(ILEVEL)) /2.)
    END WHERE
    ILEVEL=ILEVEL+1 !level just higher
  END DO
  !
  ! in this case (argument KGRID=0), input winds are ZONal and MERidian 
  !          and, output ones are in MesoNH grid   
  IF (.NOT. LCARTESIAN) THEN
    TZFIELD2(1)%CMNHNAME   = 'UM10'
    TZFIELD2(1)%CSTDNAME   = ''
    TZFIELD2(1)%CLONGNAME  = 'UM10'
    TZFIELD2(1)%CUNITS     = 'm s-1'
    TZFIELD2(1)%CDIR       = 'XY'
    TZFIELD2(1)%CCOMMENT   = 'Zonal wind at 10m'
    TZFIELD2(1)%NGRID      = 1
    TZFIELD2(1)%NTYPE      = TYPEREAL
    TZFIELD2(1)%NDIMS      = 2
    TZFIELD2(1)%LTIMEDEP   = .TRUE.
    !
    TZFIELD2(2)%CMNHNAME   = 'VM10'
    TZFIELD2(2)%CSTDNAME   = ''
    TZFIELD2(2)%CLONGNAME  = 'VM10'
    TZFIELD2(2)%CUNITS     = 'm s-1'
    TZFIELD2(2)%CDIR       = 'XY'
    TZFIELD2(2)%CCOMMENT   = 'Meridian wind at 10m'
    TZFIELD2(2)%NGRID      = 1
    TZFIELD2(2)%NTYPE      = TYPEREAL
    TZFIELD2(2)%NDIMS      = 2
    TZFIELD2(2)%LTIMEDEP   = .TRUE.
    !
    CALL UV_TO_ZONAL_AND_MERID(XCURRENT_ZON10M,XCURRENT_MER10M,KGRID=0,TPFILE=TPFILE,TZFIELDS=TZFIELD2)
  ELSE
    TZFIELD%CMNHNAME   = 'UM10'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'UM10'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'Zonal wind at 10m'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XCURRENT_ZON10M)
    !
    TZFIELD%CMNHNAME   = 'VM10'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'VM10'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'Meridian wind at 10m'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XCURRENT_MER10M)
  ENDIF
  !
  IF (SIZE(XTKET)>0) THEN
    ZWORK21(:,:) = SQRT(XCURRENT_ZON10M(:,:)**2+XCURRENT_MER10M(:,:)**2)
    ZWORK21(:,:) = ZWORK21(:,:) + 4. * SQRT(XTKET(:,:,IKB))
    TZFIELD%CMNHNAME   = 'FF10MAX'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'FF10MAX'
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_FF10MAX'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWORK21)
  END IF
  !
  IF(ANY(XCURRENT_SFCO2/=XUNDEF))THEN
    TZFIELD%CMNHNAME   = 'SFCO2'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'SFCO2'
    TZFIELD%CUNITS     = 'mg m-2 s-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'CO2 Surface flux'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XCURRENT_SFCO2)
  END IF
  !
  IF(ANY(XCURRENT_SWD/=XUNDEF))THEN
    TZFIELD%CMNHNAME   = 'SWD'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'SWD'
    TZFIELD%CUNITS     = 'W m-2'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'incoming ShortWave at the surface'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XCURRENT_SWD)
  END IF
  !
  IF(ANY(XCURRENT_SWU/=XUNDEF))THEN
    TZFIELD%CMNHNAME   = 'SWU'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'SWU'
    TZFIELD%CUNITS     = 'W m-2'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'outcoming ShortWave at the surface'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XCURRENT_SWU)
  END IF
!
  IF(ANY(XCURRENT_LWD/=XUNDEF))THEN
    TZFIELD%CMNHNAME   = 'LWD'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'LWD'
    TZFIELD%CUNITS     = 'W m-2'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'incoming LongWave at the surface'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XCURRENT_LWD)
  END IF
!
  IF(ANY(XCURRENT_LWU/=XUNDEF))THEN
    TZFIELD%CMNHNAME   = 'LWU'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'LWU'
    TZFIELD%CUNITS     = 'W m-2'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'outcoming LongWave at the surface'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 2
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,XCURRENT_LWU)
  END IF
END IF

! MODIF FP NOV 2012
!-------------------------------------------------------------------------------
!
!*       4.     DIAGNOSTIC ON PRESSURE LEVELS
!               -----------------------------
!
IF (LISOPR .AND. XISOPR(1)/=0.) THEN
!
!
ALLOCATE(ZWORK32(IIU,IJU,IKU))
ALLOCATE(ZWORK33(IIU,IJU,IKU))
ALLOCATE(ZWORK34(IIU,IJU,IKU))
!
! *************************************************
! Determine the pressure level where to interpolate
! *************************************************
  IPRES=0
  DO JI=1,SIZE(XISOPR)
    IF (XISOPR(JI)<=10..OR.XISOPR(JI)>1000.) EXIT
    IPRES=IPRES+1
    WRITE(YCAR4,'(I4)') INT(XISOPR(JI))
    YPRES(IPRES)=ADJUSTL(YCAR4)
  END DO

  ALLOCATE(ZWRES(IIU,IJU,IPRES))
  ZWRES(:,:,:)=XUNDEF
  ALLOCATE(ZPRES(IIU,IJU,IPRES))
  IPRES=0
  DO JI=1,SIZE(XISOPR)
    IF (XISOPR(JI)<=10..OR.XISOPR(JI)>1000.) EXIT
    IPRES=IPRES+1
    ZPRES(:,:,IPRES)=XISOPR(JI)*100.
  END DO
  PRINT *,'PRESSURE LEVELS WHERE TO INTERPOLATE=',ZPRES(1,1,:)
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  !
!
!*       Standard Variables
!
! *********************
! Potential Temperature
! *********************
  CALL PINTER(XTHT, XPABST, XZZ, ZTEMP, ZWRES, ZPRES, &
         IIU, IJU, IKU, IKB, IPRES, 'LOG', 'RHU.')
  DO JK=1,IPRES
    TZFIELD%CMNHNAME   = 'THT'//TRIM(YPRES(JK))//'HPA'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'K'
    TZFIELD%CCOMMENT   = 'X_Y_potential temperature '//TRIM(YPRES(JK))//' hPa'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWRES(:,:,JK))
  END DO
! *********************
! Wind
! *********************
  ZWORK31(:,:,:) = MXF(XUT(:,:,:))
  CALL PINTER(ZWORK31, XPABST, XZZ, ZTEMP, ZWRES, ZPRES, &
         IIU, IJU, IKU, IKB, IPRES, 'LOG', 'RHU.')
  DO JK=1,IPRES
    TZFIELD%CMNHNAME   = 'UT'//TRIM(YPRES(JK))//'HPA'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CCOMMENT   = 'X_Y_U component of wind '//TRIM(YPRES(JK))//' hPa'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWRES(:,:,JK))
  END DO
  !
  ZWORK31(:,:,:) = MYF(XVT(:,:,:))
  CALL PINTER(ZWORK31, XPABST, XZZ, ZTEMP, ZWRES, ZPRES, &
          IIU, IJU, IKU, IKB, IPRES, 'LOG', 'RHU.')
  DO JK=1,IPRES
    TZFIELD%CMNHNAME   = 'VT'//TRIM(YPRES(JK))//'HPA'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CCOMMENT   = 'X_Y_V component of wind '//TRIM(YPRES(JK))//' hPa'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWRES(:,:,JK))
  END DO
! *********************
! Water Vapour Mixing Ratio
! *********************
  CALL PINTER(XRT(:,:,:,1), XPABST, XZZ, ZTEMP, ZWRES, ZPRES, &
         IIU, IJU, IKU, IKB, IPRES, 'LOG', 'RHU.')
  DO JK=1,IPRES
    TZFIELD%CMNHNAME   = 'MRV'//TRIM(YPRES(JK))//'HPA'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'g kg-1'
    TZFIELD%CCOMMENT   = 'X_Y_Vapor Mixing Ratio '//TRIM(YPRES(JK))//' hPa'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWRES(:,:,JK)*1.E3)
  END DO
! *********************
! Geopotential in meters
! *********************
  ZWORK31(:,:,:) = MZF(XZZ(:,:,:))
  CALL PINTER(ZWORK31, XPABST, XZZ, ZTEMP, ZWRES, ZPRES, &
           IIU, IJU, IKU, IKB, IPRES, 'LOG', 'RHU.')
  DO JK=1,IPRES
    TZFIELD%CMNHNAME   = 'ALT'//TRIM(YPRES(JK))//'HPA'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm'
    TZFIELD%CCOMMENT   = 'X_Y_ALTitude '//TRIM(YPRES(JK))//' hPa'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWRES(:,:,JK))
  END DO
!
  DEALLOCATE(ZWRES,ZPRES,ZWORK32,ZWORK33,ZWORK34)
END IF
!
!-------------------------------------------------------------------------------
!
!*       5.     DIAGNOSTIC ON POTENTIEL TEMPERATURE LEVELS
!               -----------------------------
!
IF (LISOTH .AND.XISOTH(1)/=0.) THEN
!
!
ALLOCATE(ZWORK32(IIU,IJU,IKU))
ALLOCATE(ZWORK33(IIU,IJU,IKU))
ALLOCATE(ZWORK34(IIU,IJU,IKU))
!
! *************************************************
! Determine the potentiel temperature level where to interpolate
! *************************************************
  ITH=0
  DO JI=1,SIZE(XISOTH)
    IF (XISOTH(JI)<=100..OR.XISOTH(JI)>1000.) EXIT
    ITH=ITH+1
    WRITE(YCAR4,'(I4)') INT(XISOTH(JI))
    YTH(ITH)=ADJUSTL(YCAR4)
  END DO

  ALLOCATE(ZWTH(IIU,IJU,ITH))
  ZWTH(:,:,:)=XUNDEF
  ALLOCATE(ZTH(ITH))
  ZTH(:) = XISOTH(1:ITH)

  PRINT *,'POTENTIAL TEMPERATURE LEVELS WHERE TO INTERPOLATE=',ZTH(:)
  !
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 2
  TZFIELD%LTIMEDEP   = .TRUE.
  !
!
!*       Standard Variables
!
! *********************
! Pressure
! *********************
  CALL ZINTER(XPABST, XTHT, ZWTH, ZTH, IIU, IJU, IKU, IKB, ITH, XUNDEF)
  DO JK=1,ITH
    TZFIELD%CMNHNAME   = 'PABST'//TRIM(YTH(JK))//'K'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'Pa'
    TZFIELD%CCOMMENT   = 'X_Y_pressure '//TRIM(YTH(JK))//' K'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWTH(:,:,JK))
  END DO
! *********************
! Potential Vorticity
! *********************
  ZCORIOZ(:,:,:)=SPREAD( XCORIOZ(:,:),DIM=3,NCOPIES=IKU )
  ZVOX(:,:,:)=GY_W_VW(XWT,XDYY,XDZZ,XDZY)-GZ_V_VW(XVT,XDZZ)
  ZVOX(:,:,2)=ZVOX(:,:,3)
  ZVOY(:,:,:)=GZ_U_UW(XUT,XDZZ)-GX_W_UW(XWT,XDXX,XDZZ,XDZX)
  ZVOY(:,:,2)=ZVOY(:,:,3)
  ZVOZ(:,:,:)=GX_V_UV(XVT,XDXX,XDZZ,XDZX)-GY_U_UV(XUT,XDYY,XDZZ,XDZY)
  ZVOZ(:,:,2)=ZVOZ(:,:,3)
  ZVOZ(:,:,1)=ZVOZ(:,:,3)
  ZWORK31(:,:,:)=GX_M_M(XTHT,XDXX,XDZZ,XDZX)
  ZWORK32(:,:,:)=GY_M_M(XTHT,XDYY,XDZZ,XDZY)
  ZWORK33(:,:,:)=GZ_M_M(XTHT,XDZZ)
  ZPOVO(:,:,:)= ZWORK31(:,:,:)*MZF(MYF(ZVOX(:,:,:)))     &
  + ZWORK32(:,:,:)*MZF(MXF(ZVOY(:,:,:)))     &
   + ZWORK33(:,:,:)*(MYF(MXF(ZVOZ(:,:,:))) + ZCORIOZ(:,:,:))
  ZPOVO(:,:,:)= ZPOVO(:,:,:)*1E6/XRHODREF(:,:,:)
  ZPOVO(:,:,1)  =-1.E+11
  ZPOVO(:,:,IKU)=-1.E+11
  CALL ZINTER(ZPOVO, XTHT, ZWTH, ZTH, IIU, IJU, IKU, IKB, ITH, XUNDEF)
  DO JK=1,ITH
    TZFIELD%CMNHNAME   = 'POVOT'//TRIM(YTH(JK))//'K'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'PVU'
    TZFIELD%CCOMMENT   = 'X_Y_POtential VOrticity '//TRIM(YTH(JK))//' K'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWTH(:,:,JK))
  END DO
! *********************
! Wind
! *********************
  ZWORK31(:,:,:) = MXF(XUT(:,:,:))
  CALL ZINTER(ZWORK31, XTHT, ZWTH, ZTH, IIU, IJU, IKU, IKB, ITH, XUNDEF)
  DO JK=1,ITH
    TZFIELD%CMNHNAME   = 'UT'//TRIM(YTH(JK))//'K'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CCOMMENT   = 'X_Y_U component of wind '//TRIM(YTH(JK))//' K'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWTH(:,:,JK))
  END DO
  !
  ZWORK31(:,:,:) = MYF(XVT(:,:,:))
  CALL ZINTER(ZWORK31, XTHT, ZWTH, ZTH, IIU, IJU, IKU, IKB, ITH, XUNDEF)
  DO JK=1,ITH
    TZFIELD%CMNHNAME   = 'VT'//TRIM(YTH(JK))//'K'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CUNITS     = 'm s-1'
    TZFIELD%CCOMMENT   = 'X_Y_V component of wind '//TRIM(YTH(JK))//' K'
    CALL IO_Field_write(TPFILE,TZFIELD,ZWTH(:,:,JK))
  END DO
!
  DEALLOCATE(ZWTH,ZTH,ZWORK32,ZWORK33,ZWORK34)
END IF
!-------------------------------------------------------------------------------
!
!*       6.     DIAGNOSTIC ON ALTITUDE LEVELS
!               -----------------------------
!
IF (LISOAL .AND.XISOAL(1)/=0.) THEN
!
!
  ZFILLVAL = -99999.
  ALLOCATE(ZWORK32(IIU,IJU,IKU))
  ALLOCATE(ZWORK33(IIU,IJU,IKU))
!
! *************************************************
! Determine the altitude level where to interpolate
! *************************************************
  IAL=0
  DO JI=1,SIZE(XISOAL)
    IF (XISOAL(JI)<0.) EXIT
    IAL=IAL+1
  END DO
  ALLOCATE(ZWAL(IIU,IJU,IAL))
  ZWAL(:,:,:)=XUNDEF
  ALLOCATE(ZAL(IAL))
  ZAL(:) = XISOAL(1:IAL)
  PRINT *,'ALTITUDE LEVELS WHERE TO INTERPOLATE=',ZAL(:)
! *********************
! Altitude
! *********************
  TZFIELD%CMNHNAME   = 'ALT_ALT'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ALT_ALT'
  TZFIELD%CUNITS     = 'm'
  TZFIELD%CDIR       = '--'
  TZFIELD%CCOMMENT   = 'Z_alt ALT'
  TZFIELD%NGRID      = 0
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 1
  TZFIELD%LTIMEDEP   = .FALSE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZAL)
!
!*       Standard Variables
!
! *********************
! Cloud
! *********************
  ZWORK31(:,:,:) = 0.
  IF (SIZE(XRT,4) >= 2) ZWORK31(:,:,:) = XRT(:,:,:,2) ! Rc
  IF (SIZE(XRT,4) >= 4) ZWORK31(:,:,:) = ZWORK31(:,:,:) + XRT(:,:,:,4) !Ri
  ZWORK31(:,:,:) = ZWORK31(:,:,:)*1.E3
  CALL ZINTER(ZWORK31, XZZ, ZWAL, ZAL, IIU, IJU, IKU, IKB, IAL, XUNDEF)
  WHERE(ZWAL.EQ.XUNDEF) ZWAL=ZFILLVAL
  TZFIELD%CMNHNAME   = 'ALT_CLOUD'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ALT_CLOUD'
  TZFIELD%CUNITS     = 'g kg-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_cloud ALT'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWAL)
! *********************
! Precipitation
! *********************
  ZWORK31(:,:,:) = 0.
  IF (SIZE(XRT,4) >= 3) ZWORK31(:,:,:) = XRT(:,:,:,3) ! Rr
  IF (SIZE(XRT,4) >= 5) ZWORK31(:,:,:) = ZWORK31(:,:,:) + XRT(:,:,:,5) !Rsnow
  IF (SIZE(XRT,4) >= 6) ZWORK31(:,:,:) = ZWORK31(:,:,:) + XRT(:,:,:,6) !Rgraupel
  IF (SIZE(XRT,4) >= 7) ZWORK31(:,:,:) = ZWORK31(:,:,:) + XRT(:,:,:,7) !Rhail
  ZWORK31(:,:,:) = ZWORK31(:,:,:)*1.E3
  CALL ZINTER(ZWORK31, XZZ, ZWAL, ZAL, IIU, IJU, IKU, IKB, IAL, XUNDEF)
  WHERE(ZWAL.EQ.XUNDEF) ZWAL=ZFILLVAL
  TZFIELD%CMNHNAME   = 'ALT_PRECIP'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ALT_PRECIP'
  TZFIELD%CUNITS     = 'g kg-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_precipitation ALT'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWAL)
! *********************
! Potential temperature
! *********************
  CALL ZINTER(XTHT, XZZ, ZWAL, ZAL, IIU, IJU, IKU, IKB, IAL, XUNDEF)
  WHERE(ZWAL.EQ.XUNDEF) ZWAL=ZFILLVAL
  TZFIELD%CMNHNAME   = 'ALT_THETA'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ALT_THETA'
  TZFIELD%CUNITS     = 'K'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_potential temperature ALT'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWAL)
! *********************
! Pressure
! *********************
  CALL ZINTER(XPABST, XZZ, ZWAL, ZAL, IIU, IJU, IKU, IKB, IAL, XUNDEF)
  WHERE(ZWAL.EQ.XUNDEF) ZWAL=ZFILLVAL
  TZFIELD%CMNHNAME   = 'ALT_PRESSURE'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ALT_PRESSURE'
  TZFIELD%CUNITS     = 'Pa'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_pressure ALT'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWAL)
! *********************
! Potential Vorticity
! *********************
  ZCORIOZ(:,:,:)=SPREAD( XCORIOZ(:,:),DIM=3,NCOPIES=IKU )
  ZVOX(:,:,:)=GY_W_VW(XWT,XDYY,XDZZ,XDZY)-GZ_V_VW(XVT,XDZZ)
  ZVOX(:,:,2)=ZVOX(:,:,3)
  ZVOY(:,:,:)=GZ_U_UW(XUT,XDZZ)-GX_W_UW(XWT,XDXX,XDZZ,XDZX)
  ZVOY(:,:,2)=ZVOY(:,:,3)
  ZVOZ(:,:,:)=GX_V_UV(XVT,XDXX,XDZZ,XDZX)-GY_U_UV(XUT,XDYY,XDZZ,XDZY)
  ZVOZ(:,:,2)=ZVOZ(:,:,3)
  ZVOZ(:,:,1)=ZVOZ(:,:,3)
  ZWORK31(:,:,:)=GX_M_M(XTHT,XDXX,XDZZ,XDZX)
  ZWORK32(:,:,:)=GY_M_M(XTHT,XDYY,XDZZ,XDZY)
  ZWORK33(:,:,:)=GZ_M_M(XTHT,XDZZ)
  ZPOVO(:,:,:)= ZWORK31(:,:,:)*MZF(MYF(ZVOX(:,:,:)))     &
  + ZWORK32(:,:,:)*MZF(MXF(ZVOY(:,:,:)))     &
   + ZWORK33(:,:,:)*(MYF(MXF(ZVOZ(:,:,:))) + ZCORIOZ(:,:,:))
  ZPOVO(:,:,:)= ZPOVO(:,:,:)*1E6/XRHODREF(:,:,:)
  ZPOVO(:,:,1)  =-1.E+11
  ZPOVO(:,:,IKU)=-1.E+11
  CALL ZINTER(ZPOVO, XZZ, ZWAL, ZAL, IIU, IJU, IKU, IKB, IAL, XUNDEF)
  WHERE(ZWAL.EQ.XUNDEF) ZWAL=ZFILLVAL
  TZFIELD%CMNHNAME   = 'ALT_PV'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ALT_PV'
  TZFIELD%CUNITS     = 'PVU'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_Potential Vorticity ALT'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWAL)
! *********************
! Wind
! *********************
  ZWORK31(:,:,:) = MXF(XUT(:,:,:))
  CALL ZINTER(ZWORK31, XZZ, ZWAL, ZAL, IIU, IJU, IKU, IKB, IAL, XUNDEF)
  WHERE(ZWAL.EQ.XUNDEF) ZWAL=ZFILLVAL
  TZFIELD%CMNHNAME   = 'ALT_U'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ALT_U'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_U component of wind ALT'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWAL)
  !
  ZWORK31(:,:,:) = MYF(XVT(:,:,:))
  CALL ZINTER(ZWORK31, XZZ, ZWAL, ZAL, IIU, IJU, IKU, IKB, IAL, XUNDEF)
  WHERE(ZWAL.EQ.XUNDEF) ZWAL=ZFILLVAL
  TZFIELD%CMNHNAME   = 'ALT_V'
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'ALT_V'
  TZFIELD%CUNITS     = 'm s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'X_Y_V component of wind ALT'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWAL)
! *********************
! Dust extinction (optical depth per km)
! *********************
  IF (NRAD_3D >= 1.AND.LDUST) THEN
    DO JK=IKB,IKE
      IKRAD = JK - JPVEXT
      ZWORK31(:,:,JK)= XAER(:,:,IKRAD,3)/(XZZ(:,:,JK+1)-XZZ(:,:,JK))*1.D3
    ENDDO
    CALL ZINTER(ZWORK31, XZZ, ZWAL, ZAL, IIU, IJU, IKU, IKB, IAL, XUNDEF)
    WHERE(ZWAL.EQ.XUNDEF) ZWAL=ZFILLVAL
    TZFIELD%CMNHNAME   = 'ALT_DSTEXT'
    TZFIELD%CSTDNAME   = ''
    TZFIELD%CLONGNAME  = 'ALT_DSTEXT'
    TZFIELD%CUNITS     = 'km-1'
    TZFIELD%CDIR       = 'XY'
    TZFIELD%CCOMMENT   = 'X_Y_DuST EXTinction ALT'
    TZFIELD%NGRID      = 1
    TZFIELD%NTYPE      = TYPEREAL
    TZFIELD%NDIMS      = 3
    TZFIELD%LTIMEDEP   = .TRUE.
    CALL IO_Field_write(TPFILE,TZFIELD,ZWAL)
  END IF
!
! *********************
  DEALLOCATE(ZWAL,ZAL,ZWORK32,ZWORK33)
END IF
!
!-------------------------------------------------------------------------------
!
!*       7.     COARSE GRAINING DIAGNOSTIC
!               --------------------------
!
IF (LCOARSE) THEN
  IDX = NDXCOARSE
!-------------------------------
! AVERAGE OF TKE BY BLOCK OF IDX POINTS
  CALL BLOCKAVG(XUT,IDX,IDX,ZWORK31)
  ZUT_PRM=XUT-ZWORK31
  CALL BLOCKAVG(XVT,IDX,IDX,ZWORK31)
  ZVT_PRM=XVT-ZWORK31
  CALL BLOCKAVG(XWT,IDX,IDX,ZWORK31)
  ZWT_PRM=XWT-ZWORK31
!
  ZWORK31=MXF(ZUT_PRM*ZUT_PRM)
  CALL BLOCKAVG(ZWORK31,IDX,IDX,ZUU_AVG)
  ZWORK31=MYF(ZVT_PRM*ZVT_PRM)
  CALL BLOCKAVG(ZWORK31,IDX,IDX,ZVV_AVG)
  ZWORK31=MZF(ZWT_PRM*ZWT_PRM)
  CALL BLOCKAVG(ZWORK31,IDX,IDX,ZWW_AVG)
  CALL BLOCKAVG(XTKET,IDX,IDX,ZWORK31)
  ZWORK31=0.5*( ZUU_AVG+ZVV_AVG+ZWW_AVG ) + ZWORK31
  WRITE (YDX,FMT='(I3.3)') IDX
  TZFIELD%CMNHNAME   = 'TKEBAVG'//YDX
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'TKEBAVG'//YDX
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'TKE_BLOCKAVG'//YDX
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
!---------------------------------
! MOVING AVERAGE OF TKE OVER IDX+1 POINTS
  IDX = IDX/2
  CALL MOVINGAVG(XUT,IDX,IDX,ZWORK31)
  ZUT_PRM=XUT-ZWORK31
  CALL MOVINGAVG(XVT,IDX,IDX,ZWORK31)
  ZVT_PRM=XVT-ZWORK31
  CALL MOVINGAVG(XWT,IDX,IDX,ZWORK31)
  ZWT_PRM=XWT-ZWORK31
!
  ZWORK31=MXF(ZUT_PRM*ZUT_PRM)
  CALL MOVINGAVG(ZWORK31,IDX,IDX,ZUU_AVG)
  ZWORK31=MYF(ZVT_PRM*ZVT_PRM)
  CALL MOVINGAVG(ZWORK31,IDX,IDX,ZVV_AVG)
  ZWORK31=MZF(ZWT_PRM*ZWT_PRM)
  CALL MOVINGAVG(ZWORK31,IDX,IDX,ZWW_AVG)
  CALL MOVINGAVG(XTKET,IDX,IDX,ZWORK31)
  ZWORK31=0.5*( ZUU_AVG+ZVV_AVG+ZWW_AVG ) + ZWORK31
  WRITE (YDX,FMT='(I3.3)') 2*IDX+1
  TZFIELD%CMNHNAME   = 'TKEMAVG'//YDX
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CLONGNAME  = 'TKEMAVG'//YDX
  TZFIELD%CUNITS     = 'm2 s-2'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%CCOMMENT   = 'TKE_MOVINGAVG'//YDX
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  CALL IO_Field_write(TPFILE,TZFIELD,ZWORK31)
END IF
!
!-------------------------------------------------------------------------------
!
!*       8.     DIAGNOSTIC RELATED TO CHEMISTRY
!               -------------------------------
!
IF (NEQ_BUDGET>0) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  !
  TZFIELD%CUNITS     = 'ppp s-1'
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 4
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = 1, NEQ_BUDGET
    TZFIELD%CMNHNAME   = TRIM(CNAMES_BUDGET(JSV))//'_BUDGET'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(CNAMES_BUDGET(JSV))//'_BUDGET'
    CALL IO_Field_write(TPFILE,TZFIELD,XTCHEM(JSV)%XB_REAC(:,:,:,:))
  END DO
  !
  TZFIELD%CUNITS     = ''
  TZFIELD%NTYPE      = TYPEINT
  TZFIELD%NDIMS      = 1
  !
  DO JSV=1, NEQ_BUDGET
    TZFIELD%CMNHNAME   = TRIM(CNAMES_BUDGET(JSV))//'_CHREACLIST'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CCOMMENT   = TRIM(CNAMES_BUDGET(JSV))//'_REACTION_LIST'
    CALL IO_Field_write(TPFILE,TZFIELD,XTCHEM(JSV)%NB_REAC(:))
  END DO
END IF
!
!
! chemical prod/loss terms
IF (NEQ_PLT>0) THEN
  TZFIELD%CSTDNAME   = ''
  TZFIELD%CUNITS     = 'ppp s-1'
  TZFIELD%CDIR       = 'XY'
  TZFIELD%NGRID      = 1
  TZFIELD%NTYPE      = TYPEREAL
  TZFIELD%NDIMS      = 3
  TZFIELD%LTIMEDEP   = .TRUE.
  !
  DO JSV = 1, NEQ_PLT
    TZFIELD%CMNHNAME   = TRIM(CNAMES_PRODLOSST(JSV))//'_PROD'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(CNAMES_PRODLOSST(JSV))//'_PROD'
    CALL IO_Field_write(TPFILE,TZFIELD,XPROD(:,:,:,JSV))
    !
    TZFIELD%CMNHNAME   = TRIM(CNAMES_PRODLOSST(JSV))//'_LOSS'
    TZFIELD%CLONGNAME  = TRIM(TZFIELD%CMNHNAME)
    TZFIELD%CCOMMENT   = 'X_Y_Z_'//TRIM(CNAMES_PRODLOSST(JSV))//'_LOSS'
    CALL IO_Field_write(TPFILE,TZFIELD,XLOSS(:,:,:,JSV))
  END DO
END IF
!
!
DEALLOCATE(ZWORK21,ZWORK31,ZTEMP)
!
END SUBROUTINE WRITE_LFIFM1_FOR_DIAG_SUPP 
