!MNH_LIC Copyright 1997-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###################
      MODULE MODI_TWO_WAY_n
!     ###################
!
INTERFACE 
!
      SUBROUTINE TWO_WAY_n (KRR,KSV,PRHODJ,KMI,PTSTEP,                             &
                            PUM ,PVM, PWM, PTHM, PRM, PSVM,                        &
                            PRUS,PRVS,PRWS,PRTHS,PRRS,PRSVS,                       &
                            PINPRC,PINPRR,PINPRS,PINPRG,PINPRH,PPRCONV,PPRSCONV,   &
                            PDIRFLASWD,PSCAFLASWD,PDIRSRFSWD,OMASKkids             )
! 
INTEGER,                  INTENT(IN)  :: KRR     ! Number of moist variables
INTEGER,                  INTENT(IN)  :: KSV     ! Number of Scalar Variables
INTEGER,                  INTENT(IN)  :: KMI     ! Model index     
!
REAL,                     INTENT(IN)  :: PTSTEP  ! Timestep duration
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PRHODJ         ! (Rho) dry * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PUM, PVM, PWM  ! Variables at t-dt 
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHM
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PRM, PSVM
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS         ! Source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS, PRSVS              !  terms
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRC,PINPRR,PINPRS,PINPRG,PINPRH, &
                                           PPRCONV,PPRSCONV     !  precipitating variables
LOGICAL, DIMENSION(:,:), INTENT(INOUT)  :: OMASKkids ! true where kids exist
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PDIRFLASWD,PSCAFLASWD   ! Long wave radiation
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PDIRSRFSWD              ! Long wave radiation
!
END SUBROUTINE TWO_WAY_n
!
END INTERFACE
!
END MODULE MODI_TWO_WAY_n
!     #######################################################################
      SUBROUTINE TWO_WAY_n (KRR,KSV,PRHODJ,KMI,PTSTEP,                             &
                            PUM ,PVM, PWM, PTHM, PRM, PSVM,                        &
                            PRUS,PRVS,PRWS,PRTHS,PRRS,PRSVS,                       &
                            PINPRC,PINPRR,PINPRS,PINPRG,PINPRH,PPRCONV,PPRSCONV,   &
                            PDIRFLASWD,PSCAFLASWD,PDIRSRFSWD,OMASKkids             )
!     #######################################################################
!
!!****  *TWO_WAY_n* - Relaxation of all fields toward the average value obtained 
!!****                by the nested model $n for TWO_WAY interactive gridnesting
!!
!!    PURPOSE
!!    -------
!!      The purpose of TWO_WAY_n is:
!!          - first to average the fine scale fields of the inner model $n to
!!            the coarse mesh scale of the present outer model DAD($n).
!!          - second to apply the relaxation toward these average fields over the
!!            intersecting domain
!
!
!!**  METHOD
!!    ------
!!      Use a simple top hat horizontal average applied in the inner domain 
!!    except in a halo inner band of IHALO width (default value 0).
!!      The relaxation equation writes:
!!                                                       ___          t-1    
!!                                              |        \   rhodj * a   |
!!               d (RHODJ * A)                  | t-1    /__             |
!!               -------------- = -K  * RHODJ * |A   - ----------------- |
!!                     dt           2W          |          ___           |
!!                                              |          \   rhodj     |
!!                                              |          /__           |
!!
!!      In this routine $n denotes the nested model (with all variables X...,N...).
!!      KMI is the number of father model (all variables P..., K...)
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    MODULE MODD_CONF_n :   all
!!
!!    MODULE MODD_NESTING:   NDT_2_WAY
!!
!!    REFERENCE
!!    ---------
!!    
!!
!!    AUTHOR
!!    ------
!!    J. P. Lafore  *Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     12/11/97 
!!                   20/01/98  remove the TKE and EPS change
!!      P. Jabouille 03/04/00  parallelisation
!!      N. Asencio   18/07/05  Add the surface parameters : precipitating
!!                             hydrometeors, the Short and Long Wave 
!!                              + MASKkids array
!!                   20/05/06  Remove EPS
!!      M. Leriche   16/07/10  Add ice phase chemical species
!!      V.Masson, C.Lac 08/10  Corrections in relaxation
!!      J. Escobar   27/06/2011 correction for gridnesting with different SHAPE
!!      Bosseur & Filippi 07/2013 Adds Forefire
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1
!!      Modification    01/2016  (JP Pinty) Add LIMA
!  P. Wautelet 05/2016-04/2018: new data structures and calls for I/O
!  P. Wautelet 29/03/2019: bugfix: use correct sizes for 3rd dimension in allocation and loops when CRAD/='NONE'
!------------------------------------------------------------------------------
!
!*      0.   DECLARATIONS
!            ------------
USE MODE_ll
USE MODE_MODELN_HANDLER
!
USE MODD_PARAMETERS       ! Declarative modules
USE MODD_NESTING
USE MODD_CONF
USE MODD_NSV
USE MODD_PARAM_ICE_n,      ONLY : LSEDIC
USE MODD_PARAM_C2R2,       ONLY : LSEDC
USE MODD_PARAM_LIMA,       ONLY : NSEDC => LSEDC
!
USE MODD_FIELD_n          ! modules relative to the inner (fine scale) model $n
USE MODD_PRECIP_n          , ONLY : XINPRC,XINPRR,XINPRS,XINPRG,XINPRH
USE MODD_RADIATIONS_n      ,ONLY:XDIRFLASWD,XSCAFLASWD,XDIRSRFSWD
USE MODD_DEEP_CONVECTION_n ,ONLY : XPRCONV,XPRSCONV
USE MODD_REF_n
USE MODD_CONF_n
USE MODD_PARAM_n
USE MODI_SHUMAN
!
!
IMPLICIT NONE
!
!*       0.1   declarations of arguments 
! 
INTEGER,                  INTENT(IN)  :: KRR     ! Number of moist variables
INTEGER,                  INTENT(IN)  :: KSV     ! Number of SV (father model)
INTEGER,                  INTENT(IN)  :: KMI     ! Model index     
!
REAL,                     INTENT(IN)  :: PTSTEP  ! Timestep duration
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PRHODJ         ! (Rho) dry * Jacobian
!
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PUM, PVM, PWM  ! Variables at t-dt 
REAL, DIMENSION(:,:,:),   INTENT(IN)  :: PTHM
REAL, DIMENSION(:,:,:,:), INTENT(IN)  :: PRM, PSVM
!
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRUS, PRVS, PRWS         ! Source
REAL, DIMENSION(:,:,:),   INTENT(INOUT) :: PRTHS
REAL, DIMENSION(:,:,:,:), INTENT(INOUT) :: PRRS, PRSVS              !  terms
REAL, DIMENSION(:,:), INTENT(INOUT)     :: PINPRC,PINPRR,PINPRS,PINPRG,PINPRH &
                                          ,PPRCONV,PPRSCONV     !  precipitating variables
LOGICAL, DIMENSION(:,:), INTENT(INOUT)  :: OMASKkids ! true where kids exist
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PDIRFLASWD,PSCAFLASWD   ! Long wave radiation
REAL, DIMENSION(:,:,:), INTENT(INOUT)   :: PDIRSRFSWD              ! Long wave radiation
!
!*       0.2   declarations of local variables
!
!
INTEGER  :: IIB,IJB,IIE,IJE
INTEGER  :: IKU,IKB
INTEGER  :: II1,II2,IJ1,IJ2,II1U,IJ1V,IWEST,ISOUTH,IDIST
INTEGER  :: IXOR,IXEND !  horizontal position (i,j) of the ORigin and END
INTEGER  :: IYOR,IYEND ! of the inner model $n domain, relative to outer model subdomain
INTEGER  :: IXORU,IYORV ! particular case dure to C grid
INTEGER  :: IDXRATIO,IDYRATIO  !  x and y-direction resolution RATIO
INTEGER  :: IXOR_ll,IYOR_ll ! origin's coordinates of extended subdomain
INTEGER  :: IXDIM,IYDIM ! size of the extended dad subdomain
!
INTEGER  :: JX,JY,JVAR ! loop index
INTEGER  :: IRR,ISV_USER          ! number of moist and scalar var commun to both models
!
REAL     :: ZK2W             ! Relaxation value
REAL, DIMENSION(SIZE(PTHM,1),SIZE(PTHM,2),SIZE(PTHM,3)) :: ZAVE_RHODJ
!
! intermediate arrays for model communication
REAL, DIMENSION(:, :, :), ALLOCATABLE :: ZTUM, ZTVM, ZTWM, ZTTHM
REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: ZTRM, ZTSVM
REAL, DIMENSION(:, :, :), ALLOCATABLE :: ZUM, ZVM, ZWM, ZTHM
REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: ZRM, ZSVM
REAL, DIMENSION(:, :, :), ALLOCATABLE :: ZTRHODJ, ZTRHODJU, ZTRHODJV
REAL, DIMENSION(:, :, :), ALLOCATABLE :: ZRHODJ, ZRHODJU, ZRHODJV
REAL, DIMENSION(:, :), ALLOCATABLE ::ZTINPRC,ZTINPRR,ZTINPRS,ZTINPRG,ZTINPRH,&
                                     ZTPRCONV,ZTPRSCONV
REAL, DIMENSION(:, :,:), ALLOCATABLE :: ZTDIRFLASWD,ZTSCAFLASWD
REAL, DIMENSION(:, :,:), ALLOCATABLE :: ZTDIRSRFSWD
REAL, DIMENSION(:, :), ALLOCATABLE ::ZINPRC,ZINPRR,ZINPRS,ZINPRG,ZINPRH,&
                                     ZPRCONV,ZPRSCONV
REAL, DIMENSION(:, :,:), ALLOCATABLE :: ZDIRFLASWD,ZSCAFLASWD
REAL, DIMENSION(:, :,:), ALLOCATABLE :: ZDIRSRFSWD
!
INTEGER :: IINFO_ll, IDIMX, IDIMY ! size of intermediate arrays
INTEGER :: IHALO  ! band size where relaxation is not performed
LOGICAL :: LINTER ! flag for intersection or not with the child domain    
INTEGER :: IMI    ! Current model index KMI==NDAD(IMI)
!
INTEGER  :: IIBC,IJBC,IIEC,IJEC
!
!-------------------------------------------------------------------------------
!
!*       1.    PROLOGUE:
!
IMI = GET_CURRENT_MODEL_INDEX()
!
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
CALL GET_INDICE_ll(IIB,IJB,IIE,IJE)
!
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
CALL GET_CHILD_DIM_ll(IMI, IDIMX, IDIMY, IINFO_ll)
!
!  here we need to go back to SON domain for boundaries test
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
!
IKU = SIZE(PTHM,3)
IKB = JPVEXT+1
!
IDXRATIO = NDXRATIO_ALL(IMI)
IDYRATIO = NDYRATIO_ALL(IMI)
!
IRR = MIN(KRR,NRR)
ISV_USER = MIN(NSV_USER_A(KMI),NSV_USER_A(IMI))
!
!             1.1 Allocate array of horizontal average fields
!
ALLOCATE(ZTUM(IDIMX, IDIMY, SIZE(PUM, 3)))
ALLOCATE(ZTVM(IDIMX, IDIMY, SIZE(PUM, 3)))
ALLOCATE(ZTWM(IDIMX, IDIMY, SIZE(PUM, 3)))
ALLOCATE(ZTTHM(IDIMX, IDIMY, SIZE(PUM, 3)))
IF (IRR /= 0) THEN
  ALLOCATE(ZTRM(IDIMX, IDIMY, SIZE(PUM, 3),IRR))
              ELSE
  ALLOCATE(ZTRM(0,0,0,0))
ENDIF
IF (KSV /= 0) THEN
  ALLOCATE(ZTSVM(IDIMX, IDIMY, SIZE(PUM, 3),KSV))
ELSE
  ALLOCATE(ZTSVM(0,0,0,0))
ENDIF
!
IF (LUSERC .AND. ( (LSEDIC .AND. CCLOUD(1:3) == 'ICE') .OR.                    &
                   (LSEDC  .AND. (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO')) .OR.&
                   (NSEDC  .AND. CCLOUD == 'LIMA')                       )) THEN
  ALLOCATE(ZTINPRC(IDIMX, IDIMY))
ELSE
  ALLOCATE(ZTINPRC(0,0))
ENDIF
IF (LUSERR) THEN
  ALLOCATE(ZTINPRR(IDIMX, IDIMY))
ELSE
  ALLOCATE(ZTINPRR(0,0))
ENDIF
IF (LUSERS) THEN
  ALLOCATE(ZTINPRS(IDIMX, IDIMY))
ELSE
  ALLOCATE(ZTINPRS(0,0))
ENDIF
IF (LUSERG) THEN
  ALLOCATE(ZTINPRG(IDIMX, IDIMY))
ELSE
  ALLOCATE(ZTINPRG(0,0))
ENDIF
IF (LUSERH) THEN
  ALLOCATE(ZTINPRH(IDIMX, IDIMY))
ELSE
  ALLOCATE(ZTINPRH(0,0))
ENDIF
IF (CDCONV /= 'NONE') THEN
  ALLOCATE(ZTPRCONV (IDIMX, IDIMY))
  ALLOCATE(ZTPRSCONV(IDIMX, IDIMY))
                      ELSE
  ALLOCATE(ZTPRCONV (0,0))
  ALLOCATE(ZTPRSCONV(0,0))
END IF
IF (CRAD /= 'NONE') THEN
  ALLOCATE(ZTDIRFLASWD(IDIMX, IDIMY, SIZE(PDIRFLASWD,3)))
  ALLOCATE(ZTSCAFLASWD(IDIMX, IDIMY, SIZE(PSCAFLASWD,3)))
  ALLOCATE(ZTDIRSRFSWD(IDIMX, IDIMY, SIZE(PDIRSRFSWD,3)))
ELSE
  ALLOCATE(ZTDIRFLASWD(0,0,0))
  ALLOCATE(ZTSCAFLASWD(0,0,0))
  ALLOCATE(ZTDIRSRFSWD(0,0,0))
ENDIF
!
ALLOCATE(ZTRHODJ (IDIMX, IDIMY, SIZE(PUM, 3)))
ALLOCATE(ZTRHODJU(IDIMX, IDIMY, SIZE(PUM, 3)))
ALLOCATE(ZTRHODJV(IDIMX, IDIMY, SIZE(PUM, 3)))
!
!
ZK2W = 1. / (PTSTEP * NDT_2_WAY(NDAD(IMI)))
!
!-------------------------------------------------------------------------------
!
!*      2.    AVERAGE OF SCALAR VARIABLES
!             ---------------------------
!
IIBC=JPHEXT+2
IIEC=IDIMX-JPHEXT-1
IJBC=JPHEXT+2
IJEC=IDIMY-JPHEXT-1
!
!*      2.1  summation of rhodj
!
ZTRHODJ(:,:,:) = 0.
DO JX=1,IDXRATIO
  DO JY=1,IDYRATIO
    II1 = IIB+JX-1
    II2 = IIE+JX-IDXRATIO
    IJ1 = IJB+JY-1
    IJ2 = IJE+JY-IDYRATIO
    ZTRHODJ(IIBC:IIEC,IJBC:IJEC,:) = ZTRHODJ(IIBC:IIEC,IJBC:IJEC,:) &
                               +XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)
  END DO
END DO
!
!*      2.2  temperature
!
ZTTHM(:,:,:) = 0.
DO JX=1,IDXRATIO
  DO JY=1,IDYRATIO
    II1 = IIB+JX-1
    II2 = IIE+JX-IDXRATIO
    IJ1 = IJB+JY-1
    IJ2 = IJE+JY-IDYRATIO
    ZTTHM(IIBC:IIEC,IJBC:IJEC,:) = ZTTHM(IIBC:IIEC,IJBC:IJEC,:) &
                           +XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:) &
                             *XTHT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)
!
  END DO
END DO
!
!
!*      2.5  moist variables
!
DO JVAR=1,IRR
  ZTRM(:,:,:,JVAR) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTRM(IIBC:IIEC,IJBC:IJEC,:,JVAR) =  ZTRM(IIBC:IIEC,IJBC:IJEC,:,JVAR)  &
                             +XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)    &
                                *XRT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR)
    END DO
  END DO 
END DO
!
!*      2.6  scalar variables SV
!
! User scalar variables
IF (KSV /= 0) THEN
 DO JVAR=1,ISV_USER
  ZTSVM(:,:,:,JVAR) = 0.
   DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR) = ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR)  &
                             +XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)    &
                               *XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR)
    END DO
   END DO
 END DO
! C2R2 scalar variables 
IF (NSV_C2R2_A(IMI) > 0) THEN
  ! nested model uses C2R2 microphysical scheme
  DO JVAR=1,NSV_C2R2_A(KMI)
    ZTSVM(:,:,:,JVAR-1+NSV_C2R2BEG_A(KMI)) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_C2R2BEG_A(KMI)) = &
             &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_C2R2BEG_A(KMI))+&
             &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
             &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_C2R2BEG_A(IMI))
      END DO
    END DO
  END DO
END IF
! C1R3 scalar variables 
IF (NSV_C1R3_A(IMI) > 0) THEN
  ! nested model uses C1R3 microphysical scheme
  DO JVAR=1,NSV_C1R3_A(KMI)
    ZTSVM(:,:,:,JVAR-1+NSV_C1R3BEG_A(KMI)) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_C1R3BEG_A(KMI)) = &
             &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_C1R3BEG_A(KMI))+&
             &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
             &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_C1R3BEG_A(IMI))
      END DO
    END DO
  END DO
END IF
! LIMA scalar variables 
IF (NSV_LIMA_A(IMI) > 0) THEN
  ! nested model uses LIMA microphysical scheme
  DO JVAR=1,NSV_LIMA_A(KMI)
    ZTSVM(:,:,:,JVAR-1+NSV_LIMA_BEG_A(KMI)) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_LIMA_BEG_A(KMI)) = &
             &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_LIMA_BEG_A(KMI))+&
             &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
             &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_LIMA_BEG_A(IMI))
      END DO
    END DO
  END DO
END IF
! Electrical scalar variables 
IF (NSV_ELEC_A(IMI) > 0) THEN
  ! nested model uses electrical scheme
  DO JVAR=1,NSV_ELEC_A(KMI)
    ZTSVM(:,:,:,JVAR-1+NSV_ELECBEG_A(KMI)) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_ELECBEG_A(KMI)) = &
             &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_ELECBEG_A(KMI))+&
             &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
             &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_ELECBEG_A(IMI))
      END DO
    END DO
  END DO
END IF
! Chemical scalar variables
DO JVAR=1,NSV_CHEM_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_CHEMBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_CHEMBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_CHEMBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_CHEMBEG_A(IMI))
    END DO
  END DO
END DO
! Ice phase chemical scalar variables
IF (NSV_CHIC_A(IMI) > 0) THEN
  ! nested model uses aqueous chemistry and ice3/4 scheme
  DO JVAR=1,NSV_CHIC_A(KMI)
    ZTSVM(:,:,:,JVAR-1+NSV_CHICBEG_A(KMI)) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_CHICBEG_A(KMI)) = &
             &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_CHICBEG_A(KMI))+&
             &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
             &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_CHICBEG_A(IMI))
      END DO
    END DO
  END DO
END IF
! NOX  variables
DO JVAR=1,NSV_LNOX_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_LNOXBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_LNOXBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_LNOXBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_LNOXBEG_A(IMI))
    END DO
  END DO
END DO
! Orilam  scalar variables
DO JVAR=1,NSV_AER_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_AERBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_AERBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_AERBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_AERBEG_A(IMI))
    END DO
  END DO
END DO
DO JVAR=1,NSV_AERDEP_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_AERDEPBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_AERDEPBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_AERDEPBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_AERDEPBEG_A(IMI))
    END DO
  END DO
END DO
! Dust scalar variables
DO JVAR=1,NSV_DST_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_DSTBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_DSTBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_DSTBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_DSTBEG_A(IMI))
    END DO
  END DO
END DO
DO JVAR=1,NSV_DSTDEP_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_DSTDEPBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_DSTDEPBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_DSTDEPBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_DSTDEPBEG_A(IMI))
    END DO
  END DO
END DO
! Salt scalar variables
DO JVAR=1,NSV_SLT_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_SLTBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_SLTBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_SLTBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_SLTBEG_A(IMI))
    END DO
  END DO
END DO
DO JVAR=1,NSV_SLTDEP_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_SLTDEPBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_SLTDEPBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_SLTDEPBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_SLTDEPBEG_A(IMI))
    END DO
  END DO
END DO
! lagrangian variables
DO JVAR=1,NSV_LG_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_LGBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_LGBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_LGBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_LGBEG_A(IMI))
    END DO
  END DO
END DO
END IF
! Passive scalar variables
IF (NSV_PP_A(IMI) > 0) THEN
DO JVAR=1,NSV_PP_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_PPBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_PPBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_PPBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_PPBEG_A(IMI))
    END DO
  END DO 
END DO
END IF
#ifdef MNH_FOREFIRE
! ForeFire variables
IF (NSV_FF_A(IMI) > 0) THEN
DO JVAR=1,NSV_FF_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_FFBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_FFBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_FFBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_FFBEG_A(IMI))
    END DO
  END DO
END DO
END IF
#endif
! Conditional sampling variables
IF (NSV_CS_A(IMI) > 0) THEN
DO JVAR=1,NSV_CS_A(KMI)
  ZTSVM(:,:,:,JVAR-1+NSV_CSBEG_A(KMI)) = 0.
  DO JX=1,IDXRATIO
    DO JY=1,IDYRATIO
      II1 = IIB+JX-1
      II2 = IIE+JX-IDXRATIO
      IJ1 = IJB+JY-1
      IJ2 = IJE+JY-IDYRATIO
      ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_CSBEG_A(KMI)) = &
           &ZTSVM(IIBC:IIEC,IJBC:IJEC,:,JVAR-1+NSV_CSBEG_A(KMI))+&
           &XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)*&
           &XSVT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:,JVAR-1+NSV_CSBEG_A(IMI))
    END DO
  END DO
END DO
END IF
! Precipitating variables
  IF (LUSERR) THEN
    ZTINPRR(:,:) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTINPRR(IIBC:IIEC,IJBC:IJEC) = ZTINPRR(IIBC:IIEC,IJBC:IJEC) &
                              +XINPRR(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO)
      END DO
    END DO  
    ZTINPRR(IIBC:IIEC,IJBC:IJEC)=ZTINPRR(IIBC:IIEC,IJBC:IJEC)/(IDXRATIO*IDYRATIO)
  END IF
!
  IF (LUSERC .AND. ((LSEDIC .AND. CCLOUD(1:3) == 'ICE')      .OR.              &
                    (LSEDC .AND. (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO')) .OR.&
                    (NSEDC .AND. CCLOUD == 'LIMA')                       )) THEN
    ZTINPRC(:,:) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTINPRC(IIBC:IIEC,IJBC:IJEC) = ZTINPRC(IIBC:IIEC,IJBC:IJEC) &
                              +XINPRC(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO)
      END DO
    END DO  
    ZTINPRC(IIBC:IIEC,IJBC:IJEC)=ZTINPRC(IIBC:IIEC,IJBC:IJEC)/(IDXRATIO*IDYRATIO)
  END IF
!
  IF (LUSERS) THEN
    ZTINPRS(:,:) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTINPRS(IIBC:IIEC,IJBC:IJEC) = ZTINPRS(IIBC:IIEC,IJBC:IJEC) &
                              +XINPRS(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO)
      END DO
    END DO
    ZTINPRS(IIBC:IIEC,IJBC:IJEC) = ZTINPRS(IIBC:IIEC,IJBC:IJEC)/(IDXRATIO*IDYRATIO)
  END IF
!
  IF (LUSERG) THEN
    ZTINPRG(:,:) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTINPRG(IIBC:IIEC,IJBC:IJEC) = ZTINPRG(IIBC:IIEC,IJBC:IJEC) &
                              +XINPRG(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO)
      END DO
    END DO  
    ZTINPRG(IIBC:IIEC,IJBC:IJEC) =ZTINPRG(IIBC:IIEC,IJBC:IJEC)/(IDXRATIO*IDYRATIO)
  END IF
!
  IF (LUSERH) THEN
    ZTINPRH(:,:) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTINPRH(IIBC:IIEC,IJBC:IJEC) = ZTINPRH(IIBC:IIEC,IJBC:IJEC) &
                              +XINPRH(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO)
      END DO
    END DO  
    ZTINPRH(IIBC:IIEC,IJBC:IJEC) =ZTINPRH(IIBC:IIEC,IJBC:IJEC)/(IDXRATIO*IDYRATIO)
  END IF
!
  IF (CDCONV /= 'NONE') THEN
    ZTPRCONV(:,:) = 0.
    ZTPRSCONV(:,:) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTPRCONV(IIBC:IIEC,IJBC:IJEC) = ZTPRCONV(IIBC:IIEC,IJBC:IJEC) &
                              +XPRCONV(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO)
        ZTPRSCONV(IIBC:IIEC,IJBC:IJEC) = ZTPRSCONV(IIBC:IIEC,IJBC:IJEC) &
                              +XPRSCONV(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO)
      END DO
    END DO  
    ZTPRCONV(IIBC:IIEC,IJBC:IJEC) = ZTPRCONV(IIBC:IIEC,IJBC:IJEC)/(IDXRATIO*IDYRATIO)
    ZTPRSCONV(IIBC:IIEC,IJBC:IJEC) = ZTPRSCONV(IIBC:IIEC,IJBC:IJEC)/(IDXRATIO*IDYRATIO)
  END IF                           
! Short Wave and Long Wave variables
  IF (CRAD /= 'NONE') THEN
    ZTDIRFLASWD(:,:,:) = 0.
    ZTSCAFLASWD(:,:,:) = 0.
    ZTDIRSRFSWD(:,:,:) = 0.
    DO JX=1,IDXRATIO
      DO JY=1,IDYRATIO
        II1 = IIB+JX-1
        II2 = IIE+JX-IDXRATIO
        IJ1 = IJB+JY-1
        IJ2 = IJE+JY-IDYRATIO
        ZTDIRFLASWD(IIBC:IIEC,IJBC:IJEC,:) = ZTDIRFLASWD(IIBC:IIEC,IJBC:IJEC,:)&
                           +XDIRFLASWD(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)
        ZTSCAFLASWD(IIBC:IIEC,IJBC:IJEC,:) = ZTSCAFLASWD(IIBC:IIEC,IJBC:IJEC,:)&
                           +XSCAFLASWD(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)
        ZTDIRSRFSWD(IIBC:IIEC,IJBC:IJEC,:) = ZTDIRSRFSWD(IIBC:IIEC,IJBC:IJEC,:)&
                           +XDIRSRFSWD(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,:)
      END DO
    END DO  
    ZTDIRFLASWD(IIBC:IIEC,IJBC:IJEC,:) = ZTDIRFLASWD(IIBC:IIEC,IJBC:IJEC,:)/(IDXRATIO*IDYRATIO)
    ZTSCAFLASWD(IIBC:IIEC,IJBC:IJEC,:) = ZTSCAFLASWD(IIBC:IIEC,IJBC:IJEC,:)/(IDXRATIO*IDYRATIO)
    ZTDIRSRFSWD(IIBC:IIEC,IJBC:IJEC,:) = ZTDIRSRFSWD(IIBC:IIEC,IJBC:IJEC,:)/(IDXRATIO*IDYRATIO) 
  END IF
!
!-------------------------------------------------------------------------------
!
!*      3.    AVERAGE OF WIND VARIABLES
!             -------------------------
!
!*      3.1   vertical wind W
!
ZTWM(:,:,:) = 0.
DO JX=1,IDXRATIO
  DO JY=1,IDYRATIO
    II1 = IIB+JX-1
    II2 = IIE+JX-IDXRATIO
    IJ1 = IJB+JY-1
    IJ2 = IJE+JY-IDYRATIO
    ZTWM(IIBC:IIEC,IJBC:IJEC,IKB) = ZTWM(IIBC:IIEC,IJBC:IJEC,IKB)           &
                          +2.*XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,IKB) &
                                *XWT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,IKB)
!
    ZTWM(IIBC:IIEC,IJBC:IJEC,IKB+1:IKU) = ZTWM(IIBC:IIEC,IJBC:IJEC,IKB+1:IKU)    &
                            +(XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,IKB+1:IKU  ) &
                            + XRHODJ(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,IKB  :IKU-1))&
                                *XWT(II1:II2:IDXRATIO,IJ1:IJ2:IDYRATIO,IKB+1:IKU)
  END DO
END DO
!
!*      3.2   horizontal wind U
!
ZTRHODJU(:,:,:) = 0.
!
IF(LWEST_ll()) THEN
  II1U = IIB+IDXRATIO  !C grid
  IWEST=JPHEXT+3
ELSE
  II1U = IIB
  IWEST=JPHEXT+2
ENDIF
!
II2 = IIE+1-IDXRATIO
!
DO JY=1,IDYRATIO
  IJ1 = IJB+JY-1
  IJ2 = IJE+JY-IDYRATIO
  ZTRHODJU(IWEST:IIEC,IJBC:IJEC,:) = ZTRHODJU(IWEST:IIEC,IJBC:IJEC,:) &
                              +XRHODJ(II1U  :II2  :IDXRATIO,IJ1:IJ2:IDYRATIO,:) &
                              +XRHODJ(II1U-1:II2-1:IDXRATIO,IJ1:IJ2:IDYRATIO,:)
END DO
!
!
ZTUM(:,:,:) = 0.
DO JY=1,IDYRATIO
  IJ1 = IJB+JY-1
  IJ2 = IJE+JY-IDYRATIO
  ZTUM(IWEST:IIEC,IJBC:IJEC,:) = ZTUM(IWEST:IIEC,IJBC:IJEC,:) &
                              +(XRHODJ(II1U  :II2  :IDXRATIO,IJ1:IJ2:IDYRATIO,:) &
                              +XRHODJ(II1U-1:II2-1:IDXRATIO,IJ1:IJ2:IDYRATIO,:)) &
                                 *XUT(II1U  :II2  :IDXRATIO,IJ1:IJ2:IDYRATIO,:)
END DO
!
!
!*      3.3   horizontal wind V
!
ZTRHODJV(:,:,:) = 0.
!
IF(LSOUTH_ll() .AND. .NOT. L2D) THEN
  IJ1V = IJB+IDYRATIO  !C grid
  ISOUTH=JPHEXT+3
ELSE
  IJ1V = IJB
  ISOUTH=JPHEXT+2
ENDIF
!
IJ2 = IJE+1-IDYRATIO
!
DO JX=1,IDXRATIO
  II1 = IIB+JX-1
  II2 = IIE+JX-IDXRATIO
  ZTRHODJV(IIBC:IIEC,ISOUTH:IJEC,:) = ZTRHODJV(IIBC:IIEC,ISOUTH:IJEC,:) &
                              +XRHODJ(II1:II2:IDXRATIO,IJ1V  :IJ2  :IDYRATIO,:) &
                              +XRHODJ(II1:II2:IDXRATIO,IJ1V-1:IJ2-1:IDYRATIO,:)
END DO
!
!
ZTVM(:,:,:) = 0.
DO JX=1,IDXRATIO
  II1 = IIB+JX-1
  II2 = IIE+JX-IDXRATIO
  ZTVM(IIBC:IIEC,ISOUTH:IJEC,:) = ZTVM(IIBC:IIEC,ISOUTH:IJEC,:) &
                              +(XRHODJ(II1:II2:IDXRATIO,IJ1V  :IJ2  :IDYRATIO,:) &
                             + XRHODJ(II1:II2:IDXRATIO,IJ1V-1:IJ2-1:IDYRATIO,:)) &
                                 *XVT(II1:II2:IDXRATIO,IJ1V  :IJ2  :IDYRATIO,:)
END DO
!
!
!*        4.   EXCHANGE OF DATA
!              ----------------
!
!
CALL GO_TOMODEL_ll(IMI, IINFO_ll)
CALL GET_FEEDBACK_COORD_ll(IXOR,IYOR,IXEND,IYEND,IINFO_ll) ! physical domain's origine
!
!
IF (IINFO_ll == 0) THEN
  LINTER=.TRUE.
ELSE
  LINTER=.FALSE.
ENDIF
!
! Allocate array which will receive average child fields
!
IF (LINTER) THEN 
  ALLOCATE(ZUM(IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3)))
  ALLOCATE(ZVM(IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3)))
  ALLOCATE(ZWM(IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3)))
  ALLOCATE(ZTHM(IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3)))
  ALLOCATE(ZRHODJ (IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3)))
  ALLOCATE(ZRHODJU(IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3)))
  ALLOCATE(ZRHODJV(IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3)))
  IF (IRR /= 0) THEN
    ALLOCATE(ZRM(IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3),IRR))
  END IF
  IF (KSV /= 0) THEN
    ALLOCATE(ZSVM(IXOR:IXEND,IYOR:IYEND, SIZE(PUM, 3),KSV))
  ENDIF
  IF (LUSERR) THEN
    ALLOCATE(ZINPRR(IXOR:IXEND,IYOR:IYEND))
  ELSE
    ALLOCATE(ZINPRR(0,0))
  END IF
  IF (LUSERC .AND. ((LSEDIC .AND. CCLOUD(1:3) == 'ICE') .OR.                   &
                    (LSEDC .AND. (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO')) .OR.&
                    (NSEDC .AND. CCLOUD == 'LIMA')                       )) THEN
    ALLOCATE(ZINPRC(IXOR:IXEND,IYOR:IYEND))
  ELSE
    ALLOCATE(ZINPRC(0,0))
  END IF
  IF (LUSERS) THEN
    ALLOCATE(ZINPRS(IXOR:IXEND,IYOR:IYEND))
  ELSE
    ALLOCATE(ZINPRS(0,0))
  END IF
  IF (LUSERG) THEN
    ALLOCATE(ZINPRG(IXOR:IXEND,IYOR:IYEND))
  ELSE
    ALLOCATE(ZINPRG(0,0))
  END IF
  IF (LUSERH) THEN
    ALLOCATE(ZINPRH(IXOR:IXEND,IYOR:IYEND))
  ELSE
    ALLOCATE(ZINPRH(0,0))
  END IF
  IF (CDCONV /= 'NONE') THEN
    ALLOCATE(ZPRCONV(IXOR:IXEND,IYOR:IYEND))
    ALLOCATE(ZPRSCONV(IXOR:IXEND,IYOR:IYEND))
  ELSE
    ALLOCATE(ZPRCONV(0,0))
    ALLOCATE(ZPRSCONV(0,0))
  END IF
  IF (CRAD /= 'NONE') THEN
    ALLOCATE(ZDIRFLASWD(IXOR:IXEND,IYOR:IYEND, SIZE(PDIRFLASWD, 3)))  
    ALLOCATE(ZSCAFLASWD(IXOR:IXEND,IYOR:IYEND, SIZE(PSCAFLASWD, 3))) 
    ALLOCATE(ZDIRSRFSWD(IXOR:IXEND,IYOR:IYEND, SIZE(PDIRSRFSWD, 3)))  
  ELSE
    !3rd dimension size can also be allocated with a zero size
    ALLOCATE( ZDIRFLASWD(0, 0, SIZE( PDIRFLASWD, 3 )) )
    ALLOCATE( ZSCAFLASWD(0, 0, SIZE( PSCAFLASWD, 3 )) )
    ALLOCATE( ZDIRSRFSWD(0, 0, SIZE( PDIRSRFSWD, 3 )) )
  ENDIF
ELSE
  ALLOCATE(ZUM(0,0,0))
  ALLOCATE(ZVM(0,0,0))
  ALLOCATE(ZWM(0,0,0))
  ALLOCATE(ZTHM(0,0,0))
   IF (IRR /= 0) ALLOCATE(ZRM(0,0,0,IRR))
   IF (KSV /= 0) ALLOCATE(ZSVM(0,0,0,KSV))
  ALLOCATE(ZRHODJ (0,0,0))
  ALLOCATE(ZRHODJU(0,0,0))
  ALLOCATE(ZRHODJV(0,0,0))
  ALLOCATE(ZINPRC(0,0))
  ALLOCATE(ZINPRR(0,0))
  ALLOCATE(ZINPRS(0,0))
  ALLOCATE(ZINPRG(0,0))
  ALLOCATE(ZINPRH(0,0))
  ALLOCATE(ZPRCONV(0,0))
  ALLOCATE(ZPRSCONV(0,0))
  !3rd dimension of ZDIRFLASWD, ZSCAFLASWD and ZDIRSRFSWD is allocated with a not necessarily zero size
  !because it needs to be to this size for the SET_LSFIELD_2WAY_ll loops if CRAD/='NONE'
  ALLOCATE( ZDIRFLASWD(0, 0, SIZE( PDIRFLASWD, 3 )) )
  ALLOCATE( ZSCAFLASWD(0, 0, SIZE( PSCAFLASWD, 3 )) )
  ALLOCATE( ZDIRSRFSWD(0, 0, SIZE( PDIRSRFSWD, 3 )) )
ENDIF
!
! Initialize the list for the forcing
!
CALL SET_LSFIELD_2WAY_ll(ZUM, ZTUM)
CALL SET_LSFIELD_2WAY_ll(ZVM, ZTVM)
CALL SET_LSFIELD_2WAY_ll(ZWM, ZTWM)
CALL SET_LSFIELD_2WAY_ll(ZTHM, ZTTHM)
DO JVAR=1,IRR
  CALL SET_LSFIELD_2WAY_ll(ZRM(:,:,:,JVAR), ZTRM(:,:,:,JVAR))
ENDDO
DO JVAR=1,KSV
  CALL SET_LSFIELD_2WAY_ll(ZSVM(:,:,:,JVAR), ZTSVM(:,:,:,JVAR))
ENDDO
IF (LUSERR) THEN
  CALL SET_LSFIELD_2WAY_ll(ZINPRR , ZTINPRR) 
END IF
!
IF (LUSERC .AND. ((LSEDIC .AND. CCLOUD(1:3) == 'ICE') .OR.                   &
                  (LSEDC .AND. (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO')) .OR.&
                  (NSEDC .AND. CCLOUD == 'LIMA')                       )) THEN
  CALL SET_LSFIELD_2WAY_ll(ZINPRC , ZTINPRC) 
END IF
IF (LUSERS) THEN
  CALL SET_LSFIELD_2WAY_ll(ZINPRS , ZTINPRS) 
END IF
IF (LUSERG) THEN
  CALL SET_LSFIELD_2WAY_ll(ZINPRG , ZTINPRG) 
END IF
IF (LUSERH) THEN
  CALL SET_LSFIELD_2WAY_ll(ZINPRH , ZTINPRH) 
END IF
IF (CDCONV /= 'NONE') THEN
  CALL SET_LSFIELD_2WAY_ll(ZPRCONV , ZTPRCONV) 
  CALL SET_LSFIELD_2WAY_ll(ZPRSCONV , ZTPRSCONV) 
END IF
IF (CRAD /= 'NONE') THEN
  DO JVAR = 1, SIZE( PDIRFLASWD, 3 )
    CALL SET_LSFIELD_2WAY_ll(ZDIRFLASWD(:,:,JVAR) , ZTDIRFLASWD(:,:,JVAR))
  END DO
  DO JVAR = 1, SIZE( PSCAFLASWD, 3 )
    CALL SET_LSFIELD_2WAY_ll(ZSCAFLASWD(:,:,JVAR) , ZTSCAFLASWD(:,:,JVAR))
  END DO
  DO JVAR = 1, SIZE( PDIRSRFSWD, 3 )
    CALL SET_LSFIELD_2WAY_ll(ZDIRSRFSWD(:,:,JVAR) , ZTDIRSRFSWD(:,:,JVAR))
  END DO
END IF
CALL SET_LSFIELD_2WAY_ll(ZRHODJ, ZTRHODJ)
CALL SET_LSFIELD_2WAY_ll(ZRHODJU, ZTRHODJU)
CALL SET_LSFIELD_2WAY_ll(ZRHODJV, ZTRHODJV)
!
CALL LS_FEEDBACK_ll(IINFO_ll)
CALL GO_TOMODEL_ll(KMI, IINFO_ll)
CALL UNSET_LSFIELD_2WAY_ll(IMI)
!
DEALLOCATE(ZTUM,ZTVM,ZTWM,ZTTHM,ZTRHODJ,ZTRHODJU,ZTRHODJV)
DEALLOCATE(ZTRM,ZTSVM)
DEALLOCATE(ZTINPRC,ZTINPRR,ZTINPRS,ZTINPRG,ZTINPRH,ZTPRCONV,ZTPRSCONV)
DEALLOCATE(ZTDIRFLASWD,ZTSCAFLASWD,ZTDIRSRFSWD)
!
IF (.NOT. LINTER) THEN  ! no computation for the dad subdomain
  DEALLOCATE(ZUM,ZVM,ZWM,ZTHM,ZRHODJ,ZRHODJU,ZRHODJV)
  IF (IRR /= 0) DEALLOCATE(ZRM) 
  IF (KSV /= 0) DEALLOCATE(ZSVM)
  DEALLOCATE(ZINPRC,ZINPRR,ZINPRS,ZINPRG,ZINPRH,ZPRCONV,ZPRSCONV) 
  DEALLOCATE(ZDIRFLASWD,ZSCAFLASWD,ZDIRSRFSWD)
RETURN
ENDIF
!
!
!          5.  RELAXATION
!             -----------
!    5.1   Compute the bounds of relaxation area
!
IHALO=2
!!$IF (JPHEXT/=1) STOP ! boundaries are hard coded supposing JPHEXT=1
!
CALL GET_OR_ll('B',IXOR_ll,IYOR_ll)
CALL GET_DIM_EXT_ll('B',IXDIM,IYDIM)
!
IF(LWEST_ll()) THEN
  IDIST=IXOR_ll+1-(NXOR_ALL(IMI)+1) ! comparison of first physical
                                    ! points of subdomain and current processor
ELSE
  IDIST=IXOR_ll+NHALO-(NXOR_ALL(IMI)+1)! comparison of first physical
                                       ! points of subdomain and current processor
ENDIF
!
IF(IDIST<=0) THEN       !  west side of the child domain
  IXOR=IXOR+IHALO
ENDIF
!
IF(IDIST>=1 .AND. IDIST<=IHALO-1) THEN
  IXOR=IXOR+IHALO-IDIST
ENDIF
!
! C grid for v component
IF(IDIST >=IHALO+1)             IXORU=IXOR   ! interior child domain
IF(IDIST>=1 .AND. IDIST<=IHALO) IXORU=IXOR+1 ! partial overlapping of the relaxation area
IF(IDIST<=0)                    IXORU=IXOR+1
!
IF(LEAST_ll()) THEN
  IDIST=(NXEND_ALL(IMI)-1)-(IXOR_ll-1+IXDIM-1)  ! comparison of last physical
                                                ! points of subdomain and current processor
ELSE
  IDIST=(NXEND_ALL(IMI)-1)-(IXOR_ll-1+IXDIM-NHALO)! comparison of last physical
                                                  ! points of subdomain and current processor
ENDIF
!
IF(IDIST<=0) IXEND=IXEND-IHALO ! east side of the child domain
IF(IDIST>=1 .AND. IDIST<=IHALO-1) IXEND=IXEND-IHALO+IDIST
!
!
IF(.NOT.L2D) THEN
  IF(LSOUTH_ll()) THEN
    IDIST=IYOR_ll+1-(NYOR_ALL(IMI)+1)! comparison of first physical
                                     ! points of subdomain and current processor
  ELSE
   IDIST=IYOR_ll+NHALO-(NYOR_ALL(IMI)+1)! comparison of first physical
                                        ! points of subdomain and current processor
  ENDIF
!
  IF(IDIST<=0) THEN       ! south side of the child domain
    IYOR=IYOR+IHALO
  ENDIF
!
  IF(IDIST>=1 .AND. IDIST<=IHALO-1) THEN
    IYOR=IYOR+IHALO-IDIST
  ENDIF
!
! C grid for v component
  IF(IDIST >=IHALO+1)             IYORV=IYOR   ! interior child domain
  IF(IDIST>=1 .AND. IDIST<=IHALO) IYORV=IYOR+1 ! partial overlapping of the relaxation area
  IF(IDIST<=0)                    IYORV=IYOR+1
!
!
!
  IF(LNORTH_ll()) THEN
    IDIST=(NYEND_ALL(IMI)-1)-(IYOR_ll-1+IYDIM-1)! comparison of last physical
                                                ! points of subdomain and current processor
  ELSE
    IDIST=(NYEND_ALL(IMI)-1)-(IYOR_ll-1+IYDIM-NHALO)! comparison of last physical
                                                    ! points of subdomain and current processor
  ENDIF
  IF(IDIST<=0) IYEND=IYEND-IHALO ! north side of the child domain
  IF(IDIST>=1 .AND. IDIST<=IHALO-1) IYEND=IYEND-IHALO+IDIST
!
ELSE
  IYORV=IYOR+1   ! no parallelized
ENDIF

!  at this point, IXOR:IXEND,IYOR:IYEND define the 2way area outside
!  the relaxation area
  IF (LUSERR) THEN
    PINPRR(IXOR:IXEND,IYOR:IYEND)=ZINPRR(IXOR:IXEND,IYOR:IYEND)
  ENDIF
  IF (LUSERC .AND. ((LSEDIC .AND. CCLOUD(1:3) == 'ICE') .OR.                   &
                    (LSEDC .AND. (CCLOUD == 'C2R2' .OR. CCLOUD == 'KHKO')) .OR.&
                    (NSEDC .AND. CCLOUD == 'LIMA')                       )) THEN
    PINPRC(IXOR:IXEND,IYOR:IYEND)=ZINPRC(IXOR:IXEND,IYOR:IYEND)
  ENDIF
  IF (LUSERS) THEN
    PINPRS(IXOR:IXEND,IYOR:IYEND)=ZINPRS(IXOR:IXEND,IYOR:IYEND)
  ENDIF
  IF (LUSERG) THEN
    PINPRG(IXOR:IXEND,IYOR:IYEND)=ZINPRG(IXOR:IXEND,IYOR:IYEND)
  ENDIF
  IF (LUSERH) THEN
    PINPRH(IXOR:IXEND,IYOR:IYEND)=ZINPRH(IXOR:IXEND,IYOR:IYEND)
  ENDIF
  IF (CDCONV /= 'NONE') THEN
    PPRCONV(IXOR:IXEND,IYOR:IYEND)=ZPRCONV(IXOR:IXEND,IYOR:IYEND)
    PPRSCONV(IXOR:IXEND,IYOR:IYEND)=ZPRSCONV(IXOR:IXEND,IYOR:IYEND)
  END IF
  IF (CRAD /= 'NONE') THEN
    PDIRFLASWD(IXOR:IXEND,IYOR:IYEND,:)=ZDIRFLASWD(IXOR:IXEND,IYOR:IYEND,:)
    PSCAFLASWD(IXOR:IXEND,IYOR:IYEND,:)=ZSCAFLASWD(IXOR:IXEND,IYOR:IYEND,:)
    PDIRSRFSWD(IXOR:IXEND,IYOR:IYEND,:)=ZDIRSRFSWD(IXOR:IXEND,IYOR:IYEND,:)
  ENDIF
  DEALLOCATE(ZINPRC,ZINPRR,ZINPRS,ZINPRG,ZINPRH,ZPRCONV,ZPRSCONV)
  DEALLOCATE(ZDIRFLASWD,ZSCAFLASWD,ZDIRSRFSWD)
!
!* initialize the OMASKkids array
!
OMASKkids(IXOR:IXEND,IYOR:IYEND)=.TRUE.
!
!
!    5.2   relaxation computation
!
PRTHS(IXOR:IXEND,IYOR:IYEND,:) = PRTHS(IXOR:IXEND,IYOR:IYEND,:)        &
     - ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * ( PTHM(IXOR:IXEND,IYOR:IYEND,:) &
                 -ZTHM(IXOR:IXEND,IYOR:IYEND,:)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
!
DO JVAR=1,IRR
  PRRS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRRS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PRM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZRM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
!
! User scalar variables
DO JVAR=1,ISV_USER
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! C2R2 scalar variables
DO JVAR=NSV_C2R2BEG_A(KMI),NSV_C2R2END_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! C1R3 scalar variables
DO JVAR=NSV_C1R3BEG_A(KMI),NSV_C1R3END_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! LIMA scalar variables
DO JVAR=NSV_LIMA_BEG_A(KMI),NSV_LIMA_END_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! Electrical scalar variables
DO JVAR=NSV_ELECBEG_A(KMI),NSV_ELECEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! Chemical scalar variables
DO JVAR=NSV_CHEMBEG_A(KMI),NSV_CHEMEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! Ice phase chemical scalar variables
DO JVAR=NSV_CHICBEG_A(KMI),NSV_CHICEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! NOX variables
DO JVAR=NSV_LNOXBEG_A(KMI),NSV_LNOXEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! Orilam scalar variables
DO JVAR=NSV_AERBEG_A(KMI),NSV_AEREND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
DO JVAR=NSV_AERDEPBEG_A(KMI),NSV_AERDEPEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! Dust scalar variables
DO JVAR=NSV_DSTBEG_A(KMI),NSV_DSTEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
DO JVAR=NSV_DSTDEPBEG_A(KMI),NSV_DSTDEPEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! Salt scalar variables
DO JVAR=NSV_SLTBEG_A(KMI),NSV_SLTEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
DO JVAR=NSV_SLTDEPBEG_A(KMI),NSV_SLTDEPEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! Lagrangian scalar variables
DO JVAR=NSV_LGBEG_A(KMI),NSV_LGEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
! Passive pollutant variables
DO JVAR=NSV_PPBEG_A(KMI),NSV_PPEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
#ifdef MNH_FOREFIRE

! ForeFire variables
DO JVAR=NSV_FFBEG_A(KMI),NSV_FFEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
#endif
! Conditional sampling variables
DO JVAR=NSV_CSBEG_A(KMI),NSV_CSEND_A(KMI)
  PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) = PRSVS(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
     -  ZK2W * PRHODJ(IXOR:IXEND,IYOR:IYEND,:) * (PSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR) &
                 -ZSVM(IXOR:IXEND,IYOR:IYEND,:,JVAR)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
ENDDO
!
ZRHODJ(IXOR:IXEND,IYOR:IYEND,IKB) = 2.*ZRHODJ(IXOR:IXEND,IYOR:IYEND,IKB)
ZRHODJ(IXOR:IXEND,IYOR:IYEND,IKB+1:IKU) = ZRHODJ(IXOR:IXEND,IYOR:IYEND,IKB+1:IKU)   &
                                         +ZRHODJ(IXOR:IXEND,IYOR:IYEND,IKB:IKU-1)
!
ZAVE_RHODJ=MZM(PRHODJ)
PRWS(IXOR:IXEND,IYOR:IYEND,:) = PRWS(IXOR:IXEND,IYOR:IYEND,:)        &
     - ZK2W * ZAVE_RHODJ(IXOR:IXEND,IYOR:IYEND,:) * ( PWM(IXOR:IXEND,IYOR:IYEND,:) &
                 -ZWM(IXOR:IXEND,IYOR:IYEND,:)/ZRHODJ(IXOR:IXEND,IYOR:IYEND,:) )
!
ZAVE_RHODJ=MXM(PRHODJ)
PRUS(IXORU:IXEND,IYOR:IYEND,:) = PRUS(IXORU:IXEND,IYOR:IYEND,:)        &
     - ZK2W * ZAVE_RHODJ(IXORU:IXEND,IYOR:IYEND,:) * ( PUM(IXORU:IXEND,IYOR:IYEND,:) &
                 -ZUM(IXORU:IXEND,IYOR:IYEND,:)/ZRHODJU(IXORU:IXEND,IYOR:IYEND,:) )
!
ZAVE_RHODJ=MYM(PRHODJ)
PRVS(IXOR:IXEND,IYORV:IYEND,:) = PRVS(IXOR:IXEND,IYORV:IYEND,:)        &
     - ZK2W * ZAVE_RHODJ(IXOR:IXEND,IYORV:IYEND,:) * ( PVM(IXOR:IXEND,IYORV:IYEND,:) &
                 -ZVM(IXOR:IXEND,IYORV:IYEND,:)/ZRHODJV(IXOR:IXEND,IYORV:IYEND,:) )
!
DEALLOCATE(ZUM,ZVM,ZWM,ZTHM,ZRHODJ,ZRHODJU,ZRHODJV)
IF (IRR /= 0) DEALLOCATE(ZRM) 
IF (KSV /= 0) DEALLOCATE(ZSVM)
!------------------------------------------------------------------------------
!
END SUBROUTINE TWO_WAY_n
