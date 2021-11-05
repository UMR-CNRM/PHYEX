!MNH_LIC Copyright 2009-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######spl
     MODULE MODI_COMPUTE_MF_CLOUD
!    ############################
!
INTERFACE
!     #################################################################
      SUBROUTINE COMPUTE_MF_CLOUD(KKA,KKB,KKE,KKU,KKL,KRR,KRRL,KRRI,HMF_CLOUD,&
                                  PFRAC_ICE,                                &
                                  PRC_UP,PRI_UP,PEMF,                       &
                                  PTHL_UP, PRT_UP, PFRAC_UP,                &
                                  PTHV_UP, PFRAC_ICE_UP, PRSAT_UP,          &
                                  PEXNM, PTHLM, PRTM, PTHM, PTHVM, PRM,     &
                                  PDZZ, PZZ, KKLCL,                         &
                                  PPABSM, PRHODREF,                         &
                                  PRC_MF, PRI_MF, PCF_MF, PSIGMF, PDEPTH    )
!     #################################################################
!!
!               
!*               1.1  Declaration of Arguments
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   ::  KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
INTEGER,                INTENT(IN)   ::  KRRL         ! number of liquid water var.
INTEGER,                INTENT(IN)   ::  KRRI         ! number of ice water var.
CHARACTER (LEN=4),      INTENT(IN)   ::  HMF_CLOUD    ! Type of statistical cloud scheme
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_ICE    ! liquid/ice fraction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRC_UP,PRI_UP,PEMF ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHL_UP, PRT_UP  
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_UP
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHV_UP           ! updraft thetaV
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_ICE_UP      ! liquid/solid fraction in updraft
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRSAT_UP          ! Rsat in updraft
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PEXNM             ! exner function
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHLM, PRTM ! cons. var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHM, PTHVM ! theta and thetaV
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRM         ! water var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDZZ, PZZ
INTEGER, DIMENSION(:),  INTENT(IN)   ::  KKLCL       ! index of updraft condensation level
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PPABSM, PRHODREF ! environement
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PRC_MF, PRI_MF   ! cloud content and
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PCF_MF           ! cloud fraction for MF scheme
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PSIGMF ! SQRT(variance) for statistical cloud scheme
REAL, DIMENSION(:),     INTENT(IN)   ::  PDEPTH ! Deepness of cloud

END SUBROUTINE COMPUTE_MF_CLOUD

END INTERFACE
!
END MODULE MODI_COMPUTE_MF_CLOUD
!     ######spl
      SUBROUTINE COMPUTE_MF_CLOUD(KKA,KKB,KKE,KKU,KKL,KRR,KRRL,KRRI,HMF_CLOUD,          &
                                  PFRAC_ICE,                                &
                                  PRC_UP,PRI_UP,PEMF,                 &
                                  PTHL_UP, PRT_UP, PFRAC_UP,                &
                                  PTHV_UP, PFRAC_ICE_UP, PRSAT_UP,          &
                                  PEXNM, PTHLM, PRTM, PTHM, PTHVM, PRM,     &
                                  PDZZ, PZZ, KKLCL,                         &
                                  PPABSM, PRHODREF,                         &
                                  PRC_MF, PRI_MF, PCF_MF, PSIGMF, PDEPTH    )
!     #################################################################
!!
!!****  *COMPUTE_MF_CLOUD* - 
!!       compute diagnostic subgrid cumulus cloud caracteristics
!!
!!    PURPOSE
!!    -------
!!****  The purpose of this routine is to compute the cloud fraction and 
!!      the mean cloud content associated with clouds described by the 
!!      mass flux scheme
!!
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!     REFERENCE
!!     ---------
!!
!!
!!     AUTHOR
!!     ------
!!
!!    MODIFICATIONS
!!    -------------
!!      Original
!!      S. Riette Dec 2010 BIGA case
!!      S. Riette Aug 2011 code is split into subroutines
!!      S. Riette Jan 2012: support for both order of vertical levels
!  P. Wautelet 10/04/2019: replace ABORT and STOP calls by Print_msg
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
!
use mode_msg
!
USE MODI_COMPUTE_MF_CLOUD_BIGAUS
USE MODI_COMPUTE_MF_CLOUD_DIRECT
USE MODI_COMPUTE_MF_CLOUD_STAT
!

IMPLICIT NONE

!*                    1.1  Declaration of Arguments
!
!
!
INTEGER,                INTENT(IN)   :: KKA          ! near ground array index
INTEGER,                INTENT(IN)   :: KKB          ! near ground physical index
INTEGER,                INTENT(IN)   :: KKE          ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKU          ! uppest atmosphere array index
INTEGER,                INTENT(IN)   ::  KKL          ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER,                INTENT(IN)   ::  KRR          ! number of moist var.
INTEGER,                INTENT(IN)   ::  KRRL         ! number of liquid water var.
INTEGER,                INTENT(IN)   ::  KRRI         ! number of ice water var.
CHARACTER (LEN=4),      INTENT(IN)   ::  HMF_CLOUD    ! Type of statistical cloud scheme
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_ICE    ! liquid/ice fraction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRC_UP,PRI_UP,PEMF! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHL_UP, PRT_UP   ! rc,w,Mass Flux,Thetal,rt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_UP          ! Updraft Fraction
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHV_UP           ! updraft thetaV
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PFRAC_ICE_UP      ! liquid/solid fraction in updraft
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PRSAT_UP          ! Rsat in updraft
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PEXNM             ! exner function
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHLM, PRTM       ! cons. var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PTHM, PTHVM       ! theta and thetaV
REAL, DIMENSION(:,:,:), INTENT(IN)   ::  PRM               ! water var. at t-dt
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PDZZ, PZZ
INTEGER, DIMENSION(:),  INTENT(IN)   ::  KKLCL             ! index of updraft condensation level
REAL, DIMENSION(:,:),   INTENT(IN)   ::  PPABSM, PRHODREF  ! environement
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PRC_MF, PRI_MF    ! cloud content (INPUT=environment, OUTPUT=conv. cloud)
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PCF_MF            ! and cloud fraction for MF scheme
REAL, DIMENSION(:,:),   INTENT(OUT)  ::  PSIGMF            ! SQRT(variance) for statistical cloud scheme
REAL, DIMENSION(:),     INTENT(IN)   ::  PDEPTH            ! Deepness of cloud

!
!                       1.2  Declaration of local variables
!
!------------------------------------------------------------------------

!                     1. INITIALISATION
!
!
!                     2.1 Internal domain

PRC_MF = 0.
PRI_MF = 0.
PCF_MF = 0.
PSIGMF = 0.

IF (HMF_CLOUD == 'DIRE') THEN
  !Direct cloud scheme
  CALL COMPUTE_MF_CLOUD_DIRECT(KKE, KKL, &
                              &KKLCL(:), PFRAC_UP(:,:), PRC_UP(:,:), PRI_UP(:,:),&
                              &PRC_MF(:,:), PRI_MF(:,:), PCF_MF(:,:))
  !
ELSEIF (HMF_CLOUD == 'STAT') THEN
  !Statistical scheme using the PDF proposed by Bougeault (81, 82) and
  !Bechtold et al (95).
  CALL COMPUTE_MF_CLOUD_STAT(KKA, KKB, KKE, KKU, KKL, KRR, KRRL, KRRI,&
                            &PFRAC_ICE,&
                            &PTHLM, PRTM, PPABSM, PRM,&
                            &PDZZ, PTHM, PEXNM,&
                            &PEMF, PTHL_UP, PRT_UP,&
                            &PSIGMF)
ELSEIF (HMF_CLOUD == 'BIGA') THEN
  !Statistical scheme using the bi-gaussian PDF proposed by E. Perraud.
  CALL COMPUTE_MF_CLOUD_BIGAUS(KKA, KKB, KKE, KKU, KKL,&
                              &PEMF, PDEPTH,&
                              &PRT_UP, PTHV_UP, PFRAC_ICE_UP, PRSAT_UP,&
                              &PRTM, PTHM, PTHVM,&
                              &PDZZ, PZZ, PRHODREF,&
                              &PRC_MF, PRI_MF, PCF_MF)
  !
ELSEIF  (HMF_CLOUD == 'NONE') THEN
  ! No CONVECTIVE CLOUD SCHEME
  ! Nothing to do: PRC_MF, PRI_MF, PCF_MF, PSIGMF are already filled with zero
ELSE
  call Print_msg(NVERB_FATAL,'GEN','COMPUTE_MF_CLOUD','Shallow convection cloud scheme not valid: HMF_CLOUD='//TRIM(HMF_CLOUD))
ENDIF

END SUBROUTINE COMPUTE_MF_CLOUD
