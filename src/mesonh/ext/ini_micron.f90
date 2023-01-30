!MNH_LIC Copyright 2002-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!      ########################
       MODULE MODI_INI_MICRO_n 
!      ########################
!
INTERFACE
      SUBROUTINE INI_MICRO_n  ( TPINIFILE,KLUOUT )
!
USE MODD_IO, ONLY: TFILEDATA
!
TYPE(TFILEDATA), INTENT(IN) :: TPINIFILE ! Initial file
INTEGER,         INTENT(IN) :: KLUOUT    ! Logical unit number for prints
!
END SUBROUTINE INI_MICRO_n 
!
END INTERFACE
!
END MODULE MODI_INI_MICRO_n 
!     ############################################
      SUBROUTINE INI_MICRO_n  ( TPINIFILE,KLUOUT )
!     ############################################
!
!
!!****  *INI_MICRO_n* allocates and fills MODD_PRECIP_n variables 
!!                    and initialize parameter for microphysical scheme
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
!!      P. Jabouille
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         27/11/02
!!      O.Geoffroy (03/2006) : Add KHKO scheme
!!      Modification    01/2016  (JP Pinty) Add LIMA
!!      C.LAc          10/2016   Add budget for droplet deposition
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet    01/2019: bugfix: add missing allocations
!  C. Lac         02/2020: add missing allocation of INPRC and ACPRC with deposition
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!  P. Wautelet 04/06/2020: bugfix: correct bounds of passed arrays
! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
!
USE MODD_CONF, ONLY : CCONF,CPROGRAM       
USE MODD_IO, ONLY : TFILEDATA
USE MODD_GET_n, ONLY : CGETRCT,CGETRRT, CGETRST, CGETRGT, CGETRHT, CGETCLOUD
USE MODD_DIM_n, ONLY : NIMAX_ll, NJMAX_ll
USE MODD_PARAMETERS, ONLY : JPVEXT, JPHEXT
USE MODD_PARAM_n, ONLY : CCLOUD
USE MODD_PRECIP_n, ONLY : XINPRR, XACPRR, XINPRS, XACPRS, XINPRG, XACPRG, &
                          XINPRH, XACPRH, XINPRC, XACPRC, XINPRR3D, XEVAP3D,&
                          XINDEP,XACDEP
USE MODD_FIELD_n, ONLY : XRT, XSVT, XTHT, XPABST, XTHM, XRCM
USE MODD_GRID_n, ONLY : XZZ
USE MODD_METRICS_n, ONLY : XDXX,XDYY,XDZZ,XDZX,XDZY
USE MODD_REF_n, ONLY : XRHODREF
USE MODD_DYN_n, ONLY : XTSTEP
USE MODD_CLOUDPAR_n, ONLY : NSPLITR, NSPLITG
USE MODD_PARAM_n, ONLY : CELEC
USE MODD_PARAM_ICE,  ONLY : LSEDIC, LDEPOSC
USE MODD_PARAM_C2R2, ONLY : LSEDC, LACTIT, LDEPOC
USE MODD_BLOWSNOW
USE MODD_BLOWSNOW_n
!
USE MODI_READ_PRECIP_FIELD
USE MODI_INI_CLOUD
USE MODI_INI_RAIN_ICE
USE MODI_INI_RAIN_C2R2
USE MODI_INI_ICE_C1R3
USE MODI_CLEAN_CONC_RAIN_C2R2
USE MODI_SET_CONC_RAIN_C2R2
USE MODI_CLEAN_CONC_ICE_C1R3
USE MODI_SET_CONC_ICE_C1R3
!
USE MODE_ll
USE MODE_MODELN_HANDLER
USE MODE_BLOWSNOW_SEDIM_LKT
USE MODE_SET_CONC_LIMA
!
USE MODD_NSV,        ONLY : NSV,NSV_CHEM,NSV_C2R2BEG,NSV_C2R2END, &
                            NSV_C1R3BEG,NSV_C1R3END,              &
                            NSV_LIMA_BEG, NSV_LIMA_END
USE MODD_PARAM_LIMA, ONLY : LSCAV, MSEDC=>LSEDC, MACTIT=>LACTIT, MDEPOC=>LDEPOC
USE MODD_LIMA_PRECIP_SCAVENGING_n
!
USE MODI_INIT_AEROSOL_CONCENTRATION
USE MODI_INI_LIMA
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(TFILEDATA), INTENT(IN) :: TPINIFILE ! Initial file
INTEGER,         INTENT(IN) :: KLUOUT    ! Logical unit number for prints
!
!       0.2  declaration of local variables
!
!
!
INTEGER             :: IIU     ! Upper dimension in x direction (local)
INTEGER             :: IJU     ! Upper dimension in y direction (local)
INTEGER             :: IKU     ! Upper dimension in z direction
INTEGER             :: JK      ! loop vertical index
INTEGER             :: IINFO_ll! Return code of //routines
INTEGER             :: IKB,IKE
!
REAL, DIMENSION(:,:,:), ALLOCATABLE  :: ZDZ    ! mesh size
REAL :: ZDZMIN
INTEGER :: IMI
!
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE
!
!
CALL GET_DIM_EXT_ll('B',IIU,IJU)
IKU=SIZE(XZZ,3)
IMI = GET_CURRENT_MODEL_INDEX()
!
!
!*       2.    ALLOCATE  Module MODD_PRECIP_n
!              ------------------------------
!
IF (CCLOUD /= 'NONE' .AND. CCLOUD /= 'REVE') THEN
  ALLOCATE(XINPRR(IIU,IJU))
  ALLOCATE(XINPRR3D(IIU,IJU,IKU))
  ALLOCATE(XEVAP3D(IIU,IJU,IKU))
  ALLOCATE(XACPRR(IIU,IJU))
  XINPRR(:,:)=0.0
  XACPRR(:,:)=0.0
  XINPRR3D(:,:,:)=0.0
  XEVAP3D(:,:,:)=0.0
ELSE
  ALLOCATE(XINPRR(0,0))
  ALLOCATE(XINPRR3D(0,0,0))
  ALLOCATE(XEVAP3D(0,0,0))
  ALLOCATE(XACPRR(0,0))
END IF
!
IF (( CCLOUD(1:3) == 'ICE'                                   .AND.(LSEDIC .OR. LDEPOSC)) .OR. &
    ((CCLOUD=='C2R2' .OR. CCLOUD=='C3R5' .OR. CCLOUD=='KHKO').AND.(LSEDC .OR. LDEPOC))  .OR. &
    ( CCLOUD=='LIMA'                                         .AND.(MSEDC .OR. MDEPOC)))  THEN
  ALLOCATE(XINPRC(IIU,IJU))
  ALLOCATE(XACPRC(IIU,IJU))
  XINPRC(:,:)=0.0
  XACPRC(:,:)=0.0
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
  XINDEP(:,:)=0.0
  XACDEP(:,:)=0.0
ELSE
  ALLOCATE(XINDEP(0,0))
  ALLOCATE(XACDEP(0,0))
END IF
!
IF (CCLOUD(1:3) == 'ICE' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'LIMA') THEN
  ALLOCATE(XINPRS(IIU,IJU))
  ALLOCATE(XACPRS(IIU,IJU))
  XINPRS(:,:)=0.0
  XACPRS(:,:)=0.0
ELSE
  ALLOCATE(XINPRS(0,0))
  ALLOCATE(XACPRS(0,0))
 END IF
!
IF (CCLOUD == 'C3R5' .OR. CCLOUD(1:3) == 'ICE'.OR. CCLOUD == 'LIMA') THEN
  ALLOCATE(XINPRG(IIU,IJU))
  ALLOCATE(XACPRG(IIU,IJU))
  XINPRG(:,:)=0.0
  XACPRG(:,:)=0.0
ELSE
  ALLOCATE(XINPRG(0,0))
  ALLOCATE(XACPRG(0,0))
END IF
!
IF (CCLOUD =='ICE4' .OR. CCLOUD == 'LIMA') THEN
  ALLOCATE(XINPRH(IIU,IJU))
  ALLOCATE(XACPRH(IIU,IJU))
  XINPRH(:,:)=0.0
  XACPRH(:,:)=0.0
ELSE
  ALLOCATE(XINPRH(0,0))
  ALLOCATE(XACPRH(0,0))
END IF
!
IF(LBLOWSNOW) THEN
  ALLOCATE(XSNWSUBL3D(IIU,IJU,IKU))
  XSNWSUBL3D(:,:,:) = 0.0
  IF(CSNOWSEDIM=='TABC') THEN
!Read in look up tables of snow particles properties
!No arguments, all look up tables are defined in module
!mode_snowdrift_sedim_lkt
    CALL BLOWSNOW_SEDIM_LKT_SET
  END IF
ELSE
  ALLOCATE(XSNWSUBL3D(0,0,0))
END IF
!
!*       2b.    ALLOCATION for Radiative cooling 
!              ------------------------------
IF (LACTIT .OR. MACTIT) THEN
  ALLOCATE( XTHM(IIU,IJU,IKU) )
  ALLOCATE( XRCM(IIU,IJU,IKU) )
  XTHM = XTHT
  XRCM(:,:,:) = XRT(:,:,:,2)
ELSE
  ALLOCATE( XTHM(0,0,0) )
  ALLOCATE( XRCM(0,0,0) )
END IF
!
!*       2.bis ALLOCATE  Module MODD_PRECIP_SCAVENGING_n
!              ------------------------------
!
IF ( (CCLOUD=='LIMA') .AND. LSCAV ) THEN
  ALLOCATE(XINPAP(IIU,IJU))
  ALLOCATE(XACPAP(IIU,IJU))
  XINPAP(:,:)=0.0
  XACPAP(:,:)=0.0
ELSE
  ALLOCATE(XINPAP(0,0))
  ALLOCATE(XACPAP(0,0))
END IF
!
IF(SIZE(XINPRR) == 0) RETURN
!
!*       3.    INITIALIZE MODD_PRECIP_n variables
!              ----------------------------------
!
CALL READ_PRECIP_FIELD(TPINIFILE,CPROGRAM,CCONF,                      &
                  CGETRCT,CGETRRT,CGETRST,CGETRGT,CGETRHT,            &
                  XINPRC,XACPRC,XINDEP,XACDEP,XINPRR,XINPRR3D,XEVAP3D,&
                  XACPRR,XINPRS,XACPRS,XINPRG,XACPRG, XINPRH,XACPRH )
!
!
!*       4.    INITIALIZE THE PARAMETERS FOR THE MICROPHYSICS
!              ----------------------------------------------
!
!
!*       4.1    Compute the minimun vertical mesh size
!
ALLOCATE(ZDZ(IIU,IJU,IKU))
ZDZ=0.
IKB = 1 + JPVEXT
IKE = SIZE(XZZ,3)- JPVEXT
DO JK = IKB,IKE
  ZDZ(:,:,JK) = XZZ(:,:,JK+1) - XZZ(:,:,JK)
END DO
ZDZMIN = MIN_ll (ZDZ,IINFO_ll,1,1,IKB,NIMAX_ll+2*JPHEXT,NJMAX_ll+2*JPHEXT,IKE )
DEALLOCATE(ZDZ)
!
IF (CCLOUD(1:3) == 'KES') THEN
  CALL INI_CLOUD(XTSTEP,ZDZMIN,NSPLITR)                  ! Warm cloud only
ELSE IF (CCLOUD(1:3) == 'ICE'  ) THEN
  CALL INI_RAIN_ICE(KLUOUT,XTSTEP,ZDZMIN,NSPLITR,CCLOUD) ! Mixed phase cloud
                                                         ! including hail
ELSE IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO') THEN
  CALL INI_RAIN_C2R2(XTSTEP,ZDZMIN,NSPLITR,CCLOUD)       ! 1/2 spectral warm cloud
  IF (CCLOUD == 'C3R5') THEN
    CALL INI_ICE_C1R3(XTSTEP,ZDZMIN,NSPLITG)       ! 1/2 spectral cold cloud
  END IF
ELSE IF (CCLOUD == 'LIMA') THEN
  IF (CGETCLOUD /= 'READ') CALL INIT_AEROSOL_CONCENTRATION( XRHODREF, XSVT(:, :, :, :), XZZ(:, :, :) )
  CALL INI_LIMA(XTSTEP,ZDZMIN,NSPLITR, NSPLITG)   ! 1/2 spectral warm cloud
END IF
!
IF (CCLOUD == 'C2R2' .OR. CCLOUD == 'C3R5' .OR. CCLOUD == 'KHKO') THEN
  IF (CGETCLOUD=='READ') THEN
    CALL CLEAN_CONC_RAIN_C2R2 (XRT,XSVT(:,:,:,NSV_C2R2BEG:NSV_C2R2END))
  ELSE IF (CGETCLOUD=='INI1'.OR.CGETCLOUD=='INI2') THEN
    CALL SET_CONC_RAIN_C2R2 (CGETCLOUD,XRHODREF,&
         &XRT,XSVT(:,:,:,NSV_C2R2BEG:NSV_C2R2END))
  ENDIF
  IF (CCLOUD == 'C3R5' ) THEN
    IF (CGETCLOUD=='READ') THEN
      CALL CLEAN_CONC_ICE_C1R3 (XRT,XSVT(:,:,:,NSV_C2R2BEG:NSV_C1R3END))
    ELSE
      CALL SET_CONC_ICE_C1R3 (XRHODREF,XRT,XSVT(:,:,:,NSV_C2R2BEG:NSV_C1R3END))
    ENDIF
  ENDIF
ENDIF
!
IF (CCLOUD == 'LIMA') THEN
  IF (CGETCLOUD/='READ') THEN
    CALL SET_CONC_LIMA(IMI,CGETCLOUD,XRHODREF,XRT,XSVT(:,:,:,NSV_LIMA_BEG:NSV_LIMA_END))
  END IF
END IF
!
!
!*       5.    INITIALIZE ATMOSPHERIC ELECTRICITY
!              ----------------------------------
!
!
!IF (CELEC /= 'NONE') THEN
!  CALL INI_ELEC(IMI,TPINIFILE,XTSTEP,ZDZMIN,NSPLITR, &
!                XDXX,XDYY,XDZZ,XDZX,XDZY            )
!END IF
!
!
END SUBROUTINE INI_MICRO_n
