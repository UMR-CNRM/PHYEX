!MNH_LIC Copyright 1999-2023 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!########################
MODULE MODE_INI_ONE_WAY_n
!########################

use mode_msg

implicit none

private

public :: INI_ONE_WAY_n, Compute_ini_LB

contains

!     ####################################################################
SUBROUTINE INI_ONE_WAY_n(KDAD,KMI,                                       &
                    PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,     &
                    PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,     &
                    KDXRATIO,KDYRATIO,                                   &
                    HLBCX,HLBCY,KRIMX,KRIMY,                             &
                    KKLIN_LBXU,PCOEFLIN_LBXU,KKLIN_LBYU,PCOEFLIN_LBYU,   &
                    KKLIN_LBXV,PCOEFLIN_LBXV,KKLIN_LBYV,PCOEFLIN_LBYV,   &
                    KKLIN_LBXW,PCOEFLIN_LBXW,KKLIN_LBYW,PCOEFLIN_LBYW,   &
                    KKLIN_LBXM,PCOEFLIN_LBXM,KKLIN_LBYM,PCOEFLIN_LBYM,   &
                    HCLOUD,OUSECHAQ,OUSECHIC,                            &
                    PLBXUM,PLBYUM,PLBXVM,PLBYVM,PLBXWM,PLBYWM,           &
                    PLBXTHM,PLBYTHM,                                     &
                    PLBXTKEM,PLBYTKEM,                                   &
                    PLBXRM,PLBYRM,PLBXSVM,PLBYSVM                        )
!     ####################################################################
!
!!****  *INI_ONE_WAY$n* - INItializing a nested model Large Scale sources
!!
!!    PURPOSE
!!    -------
!!      The purpose of INI_ONE_WAY$n is to initialize Large Scale sources
!!    of all the prognostic variables of the current nested model when the
!!    current time step is in phase with its outer (DAD) model $n.
!
!
!!**  METHOD
!!    ------
!!      The basic task consists in interpolating fields from outer model $n
!!    to present inner model, using horizontal Bikhardt interpolation.
!!
!!    EXTERNAL
!!    --------
!!
!!        Function  VER_INTERP_LIN : performs the vertical interpolation
!!
!!        Subroutine  BIKHARDT : performs the horizontal interpolation
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: JPHEXT,JPVEXT
!!
!!      Module MODD_CST: XRD,XRV,XCPD,XP00,XTH00
!!
!!      Module MODD_CONF: CEQNSYS
!!
!!      Module MODD_FIELD$n : XUM,XVM,XWM,XRM,XTHM
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!    J. P. Lafore and J. Stein  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     22/09/99
!!    J.-P. Pinty    29/11/02 modify the LB*SVS for the C3R5 scheme
!!                            and add ICE2, ICE4, CELEC
!!    Modification   03/2006   (O.Geoffroy) add KHKO schem
!!    Modification   05/2006   Remove KEPS
!!    M. Leriche     11/2009  modify the LB*SVS for the aqueous phase chemistry
!!                   07/2010  idem for ice phase chemical species
!!    Bosseur & Filippi 07/2013 Adds Forefire
!!    J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!!      J.Escobar : 18/12/2015 : Correction of bug in bound in // for NHALO <>1 
!!      B.VIE   2016 : LIMA
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 14/02/2019: remove CLUOUT/CLUOUT0 and associated variables
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!  P. Wautelet 03/05/2019: restructuration of one_wayn and ini_one_wayn
!  P. Wautelet 04/06/2020: correct call to Set_conc_lima + initialize ZCONCM
!  J-L Redelsperger 06/2021: add Ocean coupling
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
!
USE MODD_CH_MNHC_n,      only: LUSECHAQ, LUSECHIC
USE MODD_CST,            only: XTH00
USE MODD_FIELD_n,        only: XRT, XSVT, XUT, XVT, XWT, XTHT, XTKET
USE MODD_NSV,            only: NSV_A, NSV_C1R3BEG_A, NSV_C1R3_A, NSV_C2R2BEG_A, NSV_C2R2_A, NSV_CHEMBEG_A, NSV_CHEMEND_A,     &
                               NSV_CHEM_A, NSV_CHICBEG_A, NSV_CHIC_A, NSV_DSTBEG_A, NSV_DSTDEPBEG_A, NSV_DSTDEP_A, NSV_DST_A, &
                               NSV_ELECBEG_A, NSV_ELEC_A, NSV_LGBEG_A, NSV_LG_A, NSV_LIMA_A, NSV_LIMA_BEG_A,                  &
                               NSV_LNOXBEG_A, NSV_LNOX_A, NSV_PPBEG_A, NSV_PP_A,                                              &

#ifdef MNH_FOREFIRE
                               NSV_FFBEG_A, NSV_FF_A, &
#endif                                              
                               NSV_SLTBEG_A, NSV_SLTDEPBEG_A, NSV_SLTDEP_A, NSV_SLT_A, NSV_USER_A

USE MODD_PARAM_n,        only: CCLOUD
USE MODD_REF,            ONLY: LCOUPLES
USE MODD_REF_n,          only: XRHODJ, XRHODREF
!
use mode_bikhardt
use mode_ll,             only: GET_CHILD_DIM_ll, GET_DIM_EXT_ll, GO_TOMODEL_ll,         &
                               LS_FORCING_ll, LWEST_ll, LEAST_ll, LNORTH_ll, LSOUTH_ll, &
                               SET_LSFIELD_1WAY_ll, UNSET_LSFIELD_1WAY_ll
USE MODE_MODELN_HANDLER, only: GOTO_MODEL
USE MODE_SET_CONC_LIMA
!
USE MODI_SET_CHEMAQ_1WAY
USE MODI_SET_CONC_ICE_C1R3
USE MODI_SET_CONC_RAIN_C2R2
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
!
INTEGER,                          INTENT(IN)  :: KDAD     !  Number of the DAD model
INTEGER,                          INTENT(IN)  :: KMI      ! model number
! interpolation coefficients
REAL,    DIMENSION(:),            INTENT(IN)  :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL,    DIMENSION(:),            INTENT(IN)  :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL,    DIMENSION(:),            INTENT(IN)  :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL,    DIMENSION(:),            INTENT(IN)  :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
!
INTEGER,                          INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,                          INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN)  :: HLBCX      ! type of lateral
CHARACTER (LEN=4), DIMENSION (2), INTENT(IN)  :: HLBCY      ! boundary conditions
INTEGER,          INTENT(IN)    :: KRIMX,KRIMY              ! size of the RIM area
!  coefficients for the vertical interpolation of the LB fields
INTEGER, DIMENSION(:,:,:),        INTENT(IN)  :: KKLIN_LBXU,KKLIN_LBYU
REAL,    DIMENSION(:,:,:),        INTENT(IN)  :: PCOEFLIN_LBXU,PCOEFLIN_LBYU
INTEGER, DIMENSION(:,:,:),        INTENT(IN)  :: KKLIN_LBXV,KKLIN_LBYV
REAL,    DIMENSION(:,:,:),        INTENT(IN)  :: PCOEFLIN_LBXV,PCOEFLIN_LBYV
INTEGER, DIMENSION(:,:,:),        INTENT(IN)  :: KKLIN_LBXW,KKLIN_LBYW
REAL,    DIMENSION(:,:,:),        INTENT(IN)  :: PCOEFLIN_LBXW,PCOEFLIN_LBYW
INTEGER, DIMENSION(:,:,:),        INTENT(IN)  :: KKLIN_LBXM,KKLIN_LBYM
REAL,    DIMENSION(:,:,:),        INTENT(IN)  :: PCOEFLIN_LBXM,PCOEFLIN_LBYM
CHARACTER (LEN=4),                INTENT(IN)  :: HCLOUD        ! Indicator of the cloud scheme
LOGICAL,                          INTENT(IN)  :: OUSECHAQ      ! logical for aqueous phase
LOGICAL,                          INTENT(IN)  :: OUSECHIC      ! logical for ice phase chemistry
!
REAL,    DIMENSION(:,:,:),        INTENT(OUT) :: PLBXUM,PLBXVM,PLBXWM ! Large Scale fields at t-dt
REAL,    DIMENSION(:,:,:),        INTENT(OUT) :: PLBYUM,PLBYVM,PLBYWM
REAL,    DIMENSION(:,:,:),        INTENT(OUT) :: PLBXTHM ,PLBYTHM  ! Large Scale fields at t-dt
REAL,    DIMENSION(:,:,:),        INTENT(OUT) :: PLBXTKEM,PLBYTKEM ! Theta, TKE
REAL,    DIMENSION(:,:,:,:),      INTENT(OUT) :: PLBXRM  ,PLBYRM   ! Moisture and SV
REAL,    DIMENSION(:,:,:,:),      INTENT(OUT) :: PLBXSVM ,PLBYSVM  ! in x and y-dir.
!
!
!*       0.2   declarations of local variables
!
INTEGER                                  :: IIB,IIE,IJB,IJE,IIU,IJU
REAL,   DIMENSION(:,:,:),    ALLOCATABLE :: ZWORK
!
INTEGER                                  :: IRR,ISV_USER !  Number of moist and user scalar variables
INTEGER                                  :: JRR,JSV      !  Loop index
INTEGER                                  :: IGRID
!
! reduced array for the interpolation coefficients
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZCOEFLIN_LBXM_RED,ZCOEFLIN_LBYM_RED
INTEGER, DIMENSION(:,:,:),   ALLOCATABLE :: IKLIN_LBXM_RED,IKLIN_LBYM_RED
!
! Variables used for LS communications
INTEGER                                  :: IINFO_ll, IDIMX, IDIMY
REAL,    DIMENSION(:,:,:),   ALLOCATABLE :: ZTUM, ZTVM, ZTWM, ZTTHM, ZTTKEM
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZTRM,ZTSVM
!
CHARACTER(LEN=4)                         :: ZINIT_TYPE ! type of C2R2 initialisation
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZCONCM  ! C2R2 concentrations
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZCHEMM  ! chemical concentrations
REAL,    DIMENSION(:,:,:,:), ALLOCATABLE :: ZCHEMMI  ! chemical ice phase concentrations
!-------------------------------------------------------------------------------
!
IF (LCOUPLES) THEN
   PLBXUM=0.
   PLBXVM=0.
   PLBXWM=0.
   PLBXTHM=0.
   PLBYTHM=0.
   PLBXTKEM=0.
   PLBYTKEM =0.
   PLBXRM =0.
   PLBYRM=0.
   PLBXSVM =0.
   PLBYSVM=0. 
RETURN
ENDIF
!*      0.   INITIALISATION
! 
CALL GOTO_MODEL(KDAD)
!
!!$CALL GET_INDICE_ll (IIB,IJB,IIE,IJE)
CALL GET_DIM_EXT_ll ('B',IIU,IJU)
IIB=1
IIE=IIU
IJB=1
IJE=IJU

ALLOCATE(ZWORK(IIB:IIE,IJB:IJE,SIZE(PLBXTHM,3)))  ! can be smaller than child extended subdomain
! LS_FORCING routine can not correctly manage extra halo zone
! LB will be filled only with one layer halo zone for the moment
!
!
!
IF (LWEST_ll() .AND. LEAST_ll()) THEN
  ALLOCATE (ZCOEFLIN_LBXM_RED(2,SIZE(PLBXTHM,2),SIZE(PLBXTHM,3)))
  ALLOCATE (  IKLIN_LBXM_RED(2,SIZE(PLBXTHM,2),SIZE(PLBXTHM,3)))
ELSE
  ALLOCATE (ZCOEFLIN_LBXM_RED(1,SIZE(PLBXTHM,2),SIZE(PLBXTHM,3)))
  ALLOCATE (  IKLIN_LBXM_RED(1,SIZE(PLBXTHM,2),SIZE(PLBXTHM,3)))
ENDIF
!
IF (LSOUTH_ll() .AND. LNORTH_ll()) THEN
  ALLOCATE (ZCOEFLIN_LBYM_RED(SIZE(PLBYTHM,1),2,SIZE(PLBYTHM,3)))
  ALLOCATE (  IKLIN_LBYM_RED(SIZE(PLBYTHM,1),2,SIZE(PLBYTHM,3)))
ELSE
  ALLOCATE (ZCOEFLIN_LBYM_RED(SIZE(PLBYTHM,1),1,SIZE(PLBYTHM,3)))
  ALLOCATE (  IKLIN_LBYM_RED(SIZE(PLBYTHM,1),1,SIZE(PLBYTHM,3)))
ENDIF
!
!
IRR=MIN(SIZE(XRT,4),SIZE(PLBXRM,4))
ISV_USER=MIN(NSV_USER_A(KDAD),NSV_USER_A(KMI))
!
IF(LWEST_ll()) THEN
  ZCOEFLIN_LBXM_RED(1,:,:)=PCOEFLIN_LBXM(1,:,:)
  IKLIN_LBXM_RED(1,:,:)=KKLIN_LBXM(1,:,:)
ENDIF
IF(LEAST_ll()) THEN
  ZCOEFLIN_LBXM_RED(SIZE(ZCOEFLIN_LBXM_RED,1),:,:) = &
             PCOEFLIN_LBXM(SIZE(PCOEFLIN_LBXM,1),:,:)

  IKLIN_LBXM_RED(SIZE(IKLIN_LBXM_RED,1),:,:) = &
             KKLIN_LBXM(SIZE(IKLIN_LBXM_RED,1),:,:)
ENDIF
IF ( SIZE(PLBYTHM,2) /= 0 ) THEN
  IF(LSOUTH_ll()) THEN
    ZCOEFLIN_LBYM_RED(:,1,:)=PCOEFLIN_LBYM(:,1,:)
    IKLIN_LBYM_RED(:,1,:)=KKLIN_LBYM(:,1,:)
  ENDIF
  IF(LNORTH_ll()) THEN
    ZCOEFLIN_LBYM_RED(:,SIZE(ZCOEFLIN_LBYM_RED,2),:) = &
               PCOEFLIN_LBYM(:,SIZE(PCOEFLIN_LBYM,2),:)
    IKLIN_LBYM_RED(:,SIZE(IKLIN_LBYM_RED,2),:) = &
               KKLIN_LBYM(:,SIZE(IKLIN_LBYM_RED,2),:)
  ENDIF
END IF
!
!*      1 GATHER LS FIELD FOR THE CHILD MODEL KMI
!
!       1.1  Must be on the father model to call get_child_dim
!
CALL GO_TOMODEL_ll(KDAD, IINFO_ll)
CALL GET_CHILD_DIM_ll(KMI, IDIMX, IDIMY, IINFO_ll)
!
!       1.2  Allocate array which will receive coarse grid points
!
ALLOCATE(ZTUM(IDIMX,IDIMY,SIZE(XUT,3)))
ALLOCATE(ZTVM(IDIMX,IDIMY,SIZE(XVT,3)))
ALLOCATE(ZTWM(IDIMX,IDIMY,SIZE(XWT,3)))
ALLOCATE(ZTTHM(IDIMX,IDIMY,SIZE(XTHT,3)))
IF (SIZE(XTKET) /= 0) ALLOCATE(ZTTKEM(IDIMX,IDIMY,SIZE(XTKET,3)))
IF (IRR /= 0) ALLOCATE(ZTRM(IDIMX,IDIMY,SIZE(XRT,3),IRR))
IF (NSV_A(KMI)/= 0) ALLOCATE(ZTSVM(IDIMX,IDIMY,SIZE(XRT,3),NSV_A(KMI)))
!
!         1.3  Specify the ls "source" fields and receiver fields
!
CALL SET_LSFIELD_1WAY_ll(XUT,ZTUM,KMI)
CALL SET_LSFIELD_1WAY_ll(XVT,ZTVM,KMI)
CALL SET_LSFIELD_1WAY_ll(XWT,ZTWM,KMI)
CALL SET_LSFIELD_1WAY_ll(XTHT,ZTTHM,KMI)
IF (ALLOCATED(ZTTKEM)) CALL SET_LSFIELD_1WAY_ll(XTKET,ZTTKEM,KMI)
!
DO JRR=1,IRR
  CALL SET_LSFIELD_1WAY_ll(XRT(:,:,:,JRR),ZTRM(:,:,:,JRR),KMI)
ENDDO
!
! USERs scalar variables
!
IF (ALLOCATED(ZTSVM)) ZTSVM=0.
DO JSV=1,ISV_USER
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV),ZTSVM(:,:,:,JSV),KMI)
ENDDO
!  Checking if it is necessary to compute the Nc and Nr
!  concentrations to use the C2R2 microphysical scheme
!  (FATHER does not use C2R2(or KHKO) and CHILD uses C2R2(or KHKO))
IF ( HCLOUD=="C2R2" .OR. HCLOUD=="KHKO" ) THEN    
 IF ( CCLOUD/="NONE" .AND. CCLOUD/="C2R2" .AND. CCLOUD/="KHKO" ) THEN
  ZINIT_TYPE="NONE"
  ALLOCATE(ZCONCM(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),3))
  IF (CCLOUD == "REVE") THEN
    ZINIT_TYPE = "INI1"
  ELSE IF (CCLOUD == "KESS" ) THEN
    ZINIT_TYPE = "INI2"
  END IF
  CALL SET_CONC_RAIN_C2R2(ZINIT_TYPE,XRHODREF,XRT,ZCONCM)
  DO JSV=1,3
    CALL SET_LSFIELD_1WAY_ll(ZCONCM(:,:,:,JSV),&
         &ZTSVM(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
  ENDDO
 ELSE
  DO JSV=1,NSV_C2R2_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
  END DO
 ENDIF
ENDIF
!
!  Checking also if it is necessary to compute the Ni
!  concentrations to use the C3R5 microphysical scheme
!  (FATHER does not use C3R5 and CHILD uses C3R5)
!
IF (HCLOUD=="C3R5"  ) THEN
 IF (CCLOUD(1:3)=="ICE") THEN
  ZINIT_TYPE="NONE"
  ALLOCATE(ZCONCM(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),5))
  IF (CCLOUD == "REVE") THEN
      ZINIT_TYPE = "INI1"
  ELSE IF (CCLOUD == "KESS" ) THEN
    ZINIT_TYPE = "INI2"
  END IF
  CALL SET_CONC_RAIN_C2R2(ZINIT_TYPE,XRHODREF,XRT,ZCONCM)
  DO JSV=1,3
    CALL SET_LSFIELD_1WAY_ll(ZCONCM(:,:,:,JSV),&
         &ZTSVM(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
  ENDDO
  ZINIT_TYPE="INI3"
  CALL SET_CONC_ICE_C1R3 (XRHODREF,XRT,ZCONCM)
  DO JSV=4,5
    CALL SET_LSFIELD_1WAY_ll(ZCONCM(:,:,:,JSV), &
         ZTSVM(:,:,:,JSV-4+NSV_C1R3BEG_A(KMI)),KMI)
  ENDDO
 ELSE
  DO JSV=1,NSV_C2R2_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C2R2BEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_C2R2BEG_A(KMI)),KMI)
  END DO
  DO JSV=1,NSV_C1R3_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_C1R3BEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_C1R3BEG_A(KMI)),KMI)
  END DO
 ENDIF
ENDIF
!
!  Checking if it is necessary to compute the Nc, Nr, Ni
!  concentrations to use the LIMA microphysical scheme
!  (FATHER does not use LIMA and CHILD uses LIMA)
!
IF (HCLOUD=="LIMA"  ) THEN
   IF (CCLOUD/="LIMA") THEN
      ALLOCATE(ZCONCM(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),NSV_LIMA_A(KMI)))
      ZCONCM(:, :, :, :) = 0.
      IF (CCLOUD == "REVE") THEN
         ZINIT_TYPE = "INI1"
      ELSE
         ZINIT_TYPE = "NONE"
      END IF
      CALL SET_CONC_LIMA (KMI,ZINIT_TYPE,XRHODREF,XRT,ZCONCM)
      DO JSV=1,NSV_LIMA_A(KMI)
         CALL SET_LSFIELD_1WAY_ll(ZCONCM(:,:,:,JSV),&
              &ZTSVM(:,:,:,JSV-1+NSV_LIMA_BEG_A(KMI)),KMI)
      ENDDO   
   ELSE
      IF (NSV_LIMA_A(KMI)/=NSV_LIMA_A(KDAD)) call Print_msg(NVERB_FATAL,'GEN','INI_ONE_WAY_n','NSV_LIMA_A(KMI)/=NSV_LIMA_A(KDAD)')
      DO JSV=1,NSV_LIMA_A(KMI)
         CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_LIMA_BEG_A(KDAD)),&
              &ZTSVM(:,:,:,JSV-1+NSV_LIMA_BEG_A(KMI)),KMI)
      END DO
   END IF
ENDIF
!
! electrical variables
!
DO JSV=1,MIN(NSV_ELEC_A(KMI),NSV_ELEC_A(KDAD))
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_ELECBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_ELECBEG_A(KMI)),KMI)
END DO
!
! chemical Scalar variables
!  Checking if it is necessary to compute the Caq
!  concentrations to use the aqueous phase chemistry
!  (FATHER does not use aqueous phase chemistry and CHILD uses it)
!
IF (OUSECHAQ) THEN
  IF (.NOT.(LUSECHAQ)) THEN
    ALLOCATE(ZCHEMM(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),&
                   NSV_CHEM_A(KMI)))
    CALL SET_CHEMAQ_1WAY(XRHODREF,&
         XSVT(:,:,:,NSV_CHEMBEG_A(KDAD):NSV_CHEMEND_A(KDAD)),ZCHEMM)
    DO JSV=1,NSV_CHEM_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(ZCHEMM(:,:,:,JSV),&
         &ZTSVM(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
    ENDDO
  ELSE
    DO JSV=1,NSV_CHEM_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KDAD)),&
            &ZTSVM(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
    END DO
  ENDIF
ELSE
  DO JSV=1,NSV_CHEM_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHEMBEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_CHEMBEG_A(KMI)),KMI)
  END DO
ENDIF
!  Checking if it is necessary to compute the Cic
!  concentrations to use the ice phase chemistry
!  (FATHER does not use ice phase chemistry and CHILD uses it)
!
IF (OUSECHIC) THEN
  IF (.NOT.(LUSECHIC)) THEN
    ALLOCATE(ZCHEMMI(SIZE(XRHODJ,1),SIZE(XRHODJ,2),SIZE(XRHODJ,3),&
                   NSV_CHIC_A(KMI)))
    ZCHEMMI(:,:,:,:) = 0.
    DO JSV=1,NSV_CHIC_A(KMI)
      CALL SET_LSFIELD_1WAY_ll(ZCHEMMI(:,:,:,JSV),&
       &ZTSVM(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
    ENDDO
  ELSE
    DO JSV=1,NSV_CHIC_A(KMI)
       CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
    END DO
  ENDIF
ELSE
  DO JSV=1,NSV_CHIC_A(KMI)
    CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_CHICBEG_A(KDAD)),&
         &ZTSVM(:,:,:,JSV-1+NSV_CHICBEG_A(KMI)),KMI)
  END DO
ENDIF
!
!
! lagrangian variables
DO JSV=1,NSV_LG_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_LGBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_LGBEG_A(KMI)),KMI)
END DO
!
! NOX                     
DO JSV=1,NSV_LNOX_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_LNOXBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_LNOXBEG_A(KMI)),KMI)
END DO
!
! Dust Scalar variables
DO JSV=1,NSV_DST_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_DSTBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_DSTBEG_A(KMI)),KMI)
END DO
!
! Moist Dust Scalar variables
DO JSV=1,NSV_DSTDEP_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_DSTDEPBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_DSTDEPBEG_A(KMI)),KMI)
END DO

! Sea Salt Scalar variables
DO JSV=1,NSV_SLT_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_SLTBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_SLTBEG_A(KMI)),KMI)
END DO
!
! Moist Sea Salt Scalar variables
DO JSV=1,NSV_SLTDEP_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_SLTDEPBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_SLTDEPBEG_A(KMI)),KMI)
END DO
!
!
! Passive pollutant      
DO JSV=1,NSV_PP_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_PPBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_PPBEG_A(KMI)),KMI)
END DO
#ifdef MNH_FOREFIRE
! ForeFire variables      
DO JSV=1,NSV_FF_A(KMI)
  CALL SET_LSFIELD_1WAY_ll(XSVT(:,:,:,JSV-1+NSV_FFBEG_A(KDAD)),&
       &ZTSVM(:,:,:,JSV-1+NSV_FFBEG_A(KMI)),KMI)
END DO
#endif
!        1.4  Communication
!
CALL LS_FORCING_ll(KMI,IINFO_ll)
!
!        1.5  Back to the (current) child model
!
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
!
CALL UNSET_LSFIELD_1WAY_ll()
IF (ALLOCATED(ZCONCM)) DEALLOCATE(ZCONCM)
IF (ALLOCATED(ZCHEMM)) DEALLOCATE(ZCHEMM)
IF (ALLOCATED(ZCHEMMI)) DEALLOCATE(ZCHEMMI)
!
!
!-------------------------------------------------------------------------------
!
!*      1.   U FIELD TREATMENT
!            -----------------
!
IGRID = 2
CALL Compute_ini_LB( PLBXUM, PLBYUM, ZTUM, ZWORK,                                  &
                     PBMX1, PBMX2, PBMX3, PBMX4, PBMY1, PBMY2, PBMY3, PBMY4,       &
                     PBFX1, PBFX2, PBFX3, PBFX4, PBFY1, PBFY2, PBFY3, PBFY4,       &
                     IIB, IIE, IJB, IJE, IGRID,                                    &
                     IDIMX, IDIMY, KDXRATIO, KDYRATIO, HLBCX, HLBCY, KRIMX, KRIMY, &
                     KKLIN_LBXU, KKLIN_LBYU,                                       &
                     PCOEFLIN_LBXU, PCOEFLIN_LBYU )
DEALLOCATE(ZTUM)
!
!-------------------------------------------------------------------------------
!
!*      2.   V FIELD TREATMENT
!            -----------------
!
IGRID = 3
CALL Compute_ini_LB( PLBXVM, PLBYVM, ZTVM, ZWORK,                                  &
                     PBMX1, PBMX2, PBMX3, PBMX4, PBMY1, PBMY2, PBMY3, PBMY4,       &
                     PBFX1, PBFX2, PBFX3, PBFX4, PBFY1, PBFY2, PBFY3, PBFY4,       &
                     IIB, IIE, IJB, IJE, IGRID,                                    &
                     IDIMX, IDIMY, KDXRATIO, KDYRATIO, HLBCX, HLBCY, KRIMX, KRIMY, &
                     KKLIN_LBXV, KKLIN_LBYV,                                       &
                     PCOEFLIN_LBXV, PCOEFLIN_LBYV )
DEALLOCATE(ZTVM)
!
!-------------------------------------------------------------------------------
!
!*      3.   W FIELD TREATMENT
!            -----------------
!
IGRID = 4
CALL Compute_ini_LB( PLBXWM, PLBYWM, ZTWM, ZWORK,                                  &
                     PBMX1, PBMX2, PBMX3, PBMX4, PBMY1, PBMY2, PBMY3, PBMY4,       &
                     PBFX1, PBFX2, PBFX3, PBFX4, PBFY1, PBFY2, PBFY3, PBFY4,       &
                     IIB, IIE, IJB, IJE, IGRID,                                    &
                     IDIMX, IDIMY, KDXRATIO, KDYRATIO, HLBCX, HLBCY, KRIMX, KRIMY, &
                     KKLIN_LBXW, KKLIN_LBYW,                                       &
                     PCOEFLIN_LBXW, PCOEFLIN_LBYW )
DEALLOCATE(ZTWM)
!
!-------------------------------------------------------------------------------
!
!*      4.   COMPUTE LARGE SCALE SOURCES FOR POTENTIAL TEMPERATURE
!            -----------------------------------------------------
!
IGRID = 1
CALL Compute_ini_LB( PLBXTHM, PLBYTHM, ZTTHM, ZWORK,                                            &
                     PBMX1, PBMX2, PBMX3, PBMX4, PBMY1, PBMY2, PBMY3, PBMY4,                    &
                     PBFX1, PBFX2, PBFX3, PBFX4, PBFY1, PBFY2, PBFY3, PBFY4,                    &
                     IIB, IIE, IJB, IJE, IGRID,                                                 &
                     IDIMX, IDIMY, KDXRATIO, KDYRATIO, HLBCX, HLBCY, KRIMX, KRIMY,              &
                     KKLIN_LBXM, KKLIN_LBYM,                                                    &
                     PCOEFLIN_LBXM, PCOEFLIN_LBYM,                                              &
                     PTH00 = XTH00,                                                             &
                     KKLIN_LBX_RED    = IKLIN_LBXM_RED,    KKLIN_LBY_RED    = IKLIN_LBYM_RED,   &
                     PCOEFLIN_LBX_RED = ZCOEFLIN_LBXM_RED, PCOEFLIN_LBY_RED = ZCOEFLIN_LBYM_RED )
!
DEALLOCATE(ZTTHM)
!
!
!-------------------------------------------------------------------------------
!
!*      5.   COMPUTE LARGE SCALE SOURCES FOR TURBULENT KINETIC ENERGY
!            --------------------------------------------------------
!
!
IF (SIZE(XTKET,3) == 0 .OR. SIZE(PLBXTKEM,3) == 0) THEN
  PLBXTKEM(:,:,:) = 0.                      ! turbulence not activated
  PLBYTKEM(:,:,:) = 0.
ELSE
  IGRID = 1
  CALL Compute_ini_LB( PLBXTKEM, PLBYTKEM, ZTTKEM, ZWORK,                                         &
                       PBMX1, PBMX2, PBMX3, PBMX4, PBMY1, PBMY2, PBMY3, PBMY4,                    &
                       PBFX1, PBFX2, PBFX3, PBFX4, PBFY1, PBFY2, PBFY3, PBFY4,                    &
                       IIB, IIE, IJB, IJE, IGRID,                                                 &
                       IDIMX, IDIMY, KDXRATIO, KDYRATIO, HLBCX, HLBCY, KRIMX, KRIMY,              &
                       KKLIN_LBXM, KKLIN_LBYM,                                                    &
                       PCOEFLIN_LBXM, PCOEFLIN_LBYM,                                              &
                       KKLIN_LBX_RED    = IKLIN_LBXM_RED,    KKLIN_LBY_RED    = IKLIN_LBYM_RED,   &
                       PCOEFLIN_LBX_RED = ZCOEFLIN_LBXM_RED, PCOEFLIN_LBY_RED = ZCOEFLIN_LBYM_RED )
  DEALLOCATE(ZTTKEM)
ENDIF
!
!-------------------------------------------------------------------------------
!
!*      6.   COMPUTE LARGE SCALE SOURCES FOR MOIST VARIABLES
!            -----------------------------------------------
!
!
IF (IRR == 0 ) THEN
  PLBXRM(:,:,:,:) = 0.                      ! water cycle not activated
  PLBYRM(:,:,:,:) = 0.
ELSE
  IGRID = 1
  DO JRR = 1,IRR
    CALL Compute_ini_LB( PLBXRM(:,:,:,JRR), PLBYRM(:,:,:,JRR), ZTRM(:,:,:,JRR), ZWORK,              &
                         PBMX1, PBMX2, PBMX3, PBMX4, PBMY1, PBMY2, PBMY3, PBMY4,                    &
                         PBFX1, PBFX2, PBFX3, PBFX4, PBFY1, PBFY2, PBFY3, PBFY4,                    &
                         IIB, IIE, IJB, IJE, IGRID,                                                 &
                         IDIMX, IDIMY, KDXRATIO, KDYRATIO, HLBCX, HLBCY, KRIMX, KRIMY,              &
                         KKLIN_LBXM, KKLIN_LBYM,                                                    &
                         PCOEFLIN_LBXM, PCOEFLIN_LBYM,                                              &
                         KKLIN_LBX_RED    = IKLIN_LBXM_RED,    KKLIN_LBY_RED    = IKLIN_LBYM_RED,   &
                         PCOEFLIN_LBX_RED = ZCOEFLIN_LBXM_RED, PCOEFLIN_LBY_RED = ZCOEFLIN_LBYM_RED )
  END DO
  DEALLOCATE(ZTRM)
!
  IF ( SIZE(PLBXRM,1) /= 0 ) PLBXRM(:,:,:,IRR+1:SIZE(PLBXRM,4)) = 0.
  IF ( SIZE(PLBYRM,1) /= 0 ) PLBYRM(:,:,:,IRR+1:SIZE(PLBYRM,4)) = 0.
!
END IF
!
!-------------------------------------------------------------------------------
!
!*      7.   COMPUTE LARGE SCALE SOURCES FOR SCALAR VARIABLES
!            ------------------------------------------------
!
!
IF (NSV_A(KMI) > 0) THEN
  IGRID = 1
  DO JSV = 1,NSV_A(KMI)
    CALL Compute_ini_LB( PLBXSVM(:,:,:,JSV),PLBYSVM(:,:,:,JSV),ZTSVM(:,:,:,JSV), ZWORK,             &
                         PBMX1, PBMX2, PBMX3, PBMX4, PBMY1, PBMY2, PBMY3, PBMY4,                    &
                         PBFX1, PBFX2, PBFX3, PBFX4, PBFY1, PBFY2, PBFY3, PBFY4,                    &
                         IIB, IIE, IJB, IJE, IGRID,                                                 &
                         IDIMX, IDIMY, KDXRATIO, KDYRATIO, HLBCX, HLBCY, KRIMX, KRIMY,              &
                         KKLIN_LBXM, KKLIN_LBYM,                                                    &
                         PCOEFLIN_LBXM, PCOEFLIN_LBYM,                                              &
                         KKLIN_LBX_RED    = IKLIN_LBXM_RED,    KKLIN_LBY_RED    = IKLIN_LBYM_RED,   &
                         PCOEFLIN_LBX_RED = ZCOEFLIN_LBXM_RED, PCOEFLIN_LBY_RED = ZCOEFLIN_LBYM_RED )
  END DO
  DEALLOCATE(ZTSVM)
ELSE
  PLBXSVM(:,:,:,:) = 0.
  PLBYSVM(:,:,:,:) = 0.
END IF
!
DEALLOCATE(ZWORK)
DEALLOCATE(ZCOEFLIN_LBXM_RED,ZCOEFLIN_LBYM_RED,IKLIN_LBXM_RED,IKLIN_LBYM_RED)
!
CALL GOTO_MODEL(KMI)
!------------------------------------------------------------------------------
!
END SUBROUTINE INI_ONE_WAY_n



!#################################################################################
SUBROUTINE Compute_ini_LB(PLBX,PLBY,PTFIELD,PWORK,                               &
                          PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4,       &
                          PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4,       &
                          KIB,KIE,KJB,KJE, KGRID,                                &
                          KDIMX,KDIMY,KDXRATIO,KDYRATIO,HLBCX,HLBCY,KRIMX,KRIMY, &
                          KKLIN_LBX,KKLIN_LBY,                                   &
                          PCOEFLIN_LBX,PCOEFLIN_LBY,                             &
                          PTH00,                                                 &
                          KKLIN_LBX_RED,KKLIN_LBY_RED,                           &
                          PCOEFLIN_LBX_RED,PCOEFLIN_LBY_RED )
!#################################################################################
!
use modd_parameters, only: jphext

use mode_bikhardt
use mode_ll,         only: LWEST_ll, LEAST_ll, LNORTH_ll, LSOUTH_ll

use modi_ver_interp_lin

IMPLICIT NONE
!
!*       0.1   declarations of arguments
!
REAL,    DIMENSION(:,:,:),           INTENT(OUT) :: PLBX,PLBY ! source term
REAL,    DIMENSION(:,:,:),           INTENT(IN)  :: PTFIELD   ! ls forcing array
REAL,    DIMENSION(:,:,:),           INTENT(OUT) :: PWORK
! interpolation coefficients
REAL,    DIMENSION(:),               INTENT(IN)  :: PBMX1,PBMX2,PBMX3,PBMX4 ! Mass points in X-direc.
REAL,    DIMENSION(:),               INTENT(IN)  :: PBMY1,PBMY2,PBMY3,PBMY4 ! Mass points in Y-direc.
REAL,    DIMENSION(:),               INTENT(IN)  :: PBFX1,PBFX2,PBFX3,PBFX4 ! Flux points in X-direc.
REAL,    DIMENSION(:),               INTENT(IN)  :: PBFY1,PBFY2,PBFY3,PBFY4 ! Flux points in Y-direc.
INTEGER,                             INTENT(IN)  :: KIB,KIE,KJB,KJE
INTEGER,                             INTENT(IN)  :: KGRID      ! code of grid point
INTEGER,                             INTENT(IN)  :: KDIMX, KDIMY
INTEGER,                             INTENT(IN)  :: KDXRATIO   !  x and y-direction resolution RATIO
INTEGER,                             INTENT(IN)  :: KDYRATIO   ! between inner model and outer model
CHARACTER (LEN=4), DIMENSION (2),    INTENT(IN)  :: HLBCX   ! type of lateral
CHARACTER (LEN=4), DIMENSION (2),    INTENT(IN)  :: HLBCY   ! boundary conditions
INTEGER,                             INTENT(IN)  :: KRIMX,KRIMY ! size of the RIM area
INTEGER, DIMENSION(:,:,:),           INTENT(IN)  :: KKLIN_LBX,KKLIN_LBY
REAL,    DIMENSION(:,:,:),           INTENT(IN)  :: PCOEFLIN_LBX,PCOEFLIN_LBY
REAL,                      OPTIONAL, INTENT(IN)  :: PTH00 ! reference temperature
INTEGER, DIMENSION(:,:,:), optional, INTENT(IN)  :: KKLIN_LBX_RED,KKLIN_LBY_RED
REAL,    DIMENSION(:,:,:), optional, INTENT(in)  :: PCOEFLIN_LBX_RED,PCOEFLIN_LBY_RED
!
!*       0.2   declarations of local variables
!
INTEGER :: ILBX, ILBY, ILBX2, ILBY2
INTEGER :: IW, IE, IN, IS
!

if ( kgrid<1 .or. kgrid>4 ) call Print_msg( NVERB_FATAL, 'GEN', 'Compute_LB', 'invalid kgrid dummy argument' )

IF(PRESENT(PTH00)) THEN
  PLBX(:,:,:) = PTH00 ! to avoid undefined computation
  PLBY(:,:,:) = PTH00
ELSE
  PLBX=0.
  PLBY=0.
ENDIF
!
!*    Horizontal Bikhardt interpolation
!
!
CALL BIKHARDT (PBMX1,PBMX2,PBMX3,PBMX4,PBMY1,PBMY2,PBMY3,PBMY4, &
               PBFX1,PBFX2,PBFX3,PBFX4,PBFY1,PBFY2,PBFY3,PBFY4, &
               2,2,KDIMX-1,KDIMY-1,KDXRATIO,KDYRATIO,KGRID,     &
               HLBCX,HLBCY,PTFIELD,PWORK)
!
ILBX2=SIZE(PLBX,1)
IF(LWEST_ll( ).AND.LEAST_ll( )) THEN
  ILBX=ILBX2/2
ELSE
  ILBX=ILBX2
ENDIF
!
IF(LWEST_ll( ) .AND. ILBX/=0) THEN
  iw = kib
  ie = kib -1 + ilbx
  if ( kgrid == 2 ) then
    iw = iw + 1
    ie = ie + 1
  end if
  PLBX(1:ILBX,KJB:KJE,:) = PWORK(iw:ie,KJB:KJE,: )
ENDIF
!
IF(LEAST_ll( ) .AND. ILBX/=0) THEN
  iw = kie + 1 - ilbx
  ie = kie
  PLBX(ILBX2-ILBX+1:ILBX2,KJB:KJE,:) = PWORK(iw:ie,KJB:KJE,:)
ENDIF
!
ILBY2=SIZE(PLBY,2)
IF(LSOUTH_ll( ).AND.LNORTH_ll( )) THEN
  ILBY=ILBY2/2
ELSE
  ILBY=ILBY2
ENDIF
!
IF(LSOUTH_ll( ) .AND. ILBY/=0) THEN
  is = kjb
  in = kjb - 1 + ilby
  if ( kgrid == 3 ) then
    is = is + 1
    in = in + 1
  end if
  PLBY(KIB:KIE,1:ILBY,:) = PWORK(KIB:KIE,is:in,:)
ENDIF
!
IF(LNORTH_ll( ) .AND. ILBY/=0) THEN
  is = kje + 1 - ilby
  in = kje
  PLBY(KIB:KIE,ILBY2-ILBY+1:ILBY2,:) = PWORK(KIB:KIE,is:in,:)
ENDIF
!
!*    Vertical interpolation
!
IF ( SIZE(PLBX,1) /= 0 ) THEN
  if ( present( kklin_lbx_red ) .and. ilbx /= krimx+jphext ) then
    PLBX(:,:,:) = VER_INTERP_LIN( PLBX(:,:,:), KKLIN_LBX_RED(:,:,:), PCOEFLIN_LBX_RED(:,:,:) )
  else
    PLBX(:,:,:) = VER_INTERP_LIN( PLBX(:,:,:), KKLIN_LBX(:,:,:),     PCOEFLIN_LBX(:,:,:) )
  end if
END IF
!
IF ( SIZE(PLBY,1) /= 0 ) THEN
  if ( present( kklin_lby_red ) .and. ilby /= krimy+jphext ) then
    PLBY(:,:,:) = VER_INTERP_LIN( PLBY(:,:,:), KKLIN_LBY_RED(:,:,:), PCOEFLIN_LBY_RED(:,:,:) )
  else
    PLBY(:,:,:) = VER_INTERP_LIN( PLBY(:,:,:), KKLIN_LBY(:,:,:),     PCOEFLIN_LBY(:,:,:) )
  end if
END IF
!
END SUBROUTINE Compute_ini_LB

END MODULE MODE_INI_ONE_WAY_n
