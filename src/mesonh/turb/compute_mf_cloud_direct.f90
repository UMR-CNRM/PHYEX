!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!     ######spl
     MODULE MODI_COMPUTE_MF_CLOUD_DIRECT
!    ###################################
!
INTERFACE
!     #################################################################
      SUBROUTINE COMPUTE_MF_CLOUD_DIRECT(KKE, KKL, &
                                        &KKLCL, PFRAC_UP, PRC_UP, PRI_UP,&
                                        &PRC_MF, PRI_MF, PCF_MF)
!     #################################################################
!!
!
!*               1.1  Declaration of Arguments
!
INTEGER,                INTENT(IN)   :: KKE            ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKL            ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER, DIMENSION(:),  INTENT(IN)   :: KKLCL          ! index of updraft condensation level
REAL, DIMENSION(:,:),   INTENT(IN)   :: PFRAC_UP       ! Updraft Fraction
REAL, DIMENSION(:,:),   INTENT(IN)   :: PRC_UP,PRI_UP  ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PRC_MF, PRI_MF ! cloud content (INPUT=environment, OUTPUT=conv. cloud)
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PCF_MF         ! and cloud fraction for MF scheme

END SUBROUTINE COMPUTE_MF_CLOUD_DIRECT

END INTERFACE
!
END MODULE MODI_COMPUTE_MF_CLOUD_DIRECT
!     ######spl
      SUBROUTINE COMPUTE_MF_CLOUD_DIRECT(KKE, KKL, &
                                        &KKLCL, PFRAC_UP, PRC_UP, PRI_UP,&
                                        &PRC_MF, PRI_MF, PCF_MF)
!     #################################################################
!!
!!****  *COMPUTE_MF_CLOUD_DIRECT* -
!!       compute diagnostic subgrid cumulus cloud caracteristics with a direct scheme
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
!!      Updraft variables are used directly to diagnose subgrid clouds
!!      This scheme may be activated only if the selected updraft model
!!      gives the updraft fraction as an output
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
!!     S. Riette moving of code previously in compute_mf_cloud code
!!
!!    MODIFICATIONS
!!    -------------
!!      Original 25 Aug 2011
!!      S. Riette Jan 2012: support for both order of vertical levels
!! --------------------------------------------------------------------------
!
!*      0. DECLARATIONS
!          ------------
USE MODD_PARAM_MFSHALL_n, ONLY : XKCF_MF
!
IMPLICIT NONE
!
!*                    0.1  Declaration of Arguments
!
INTEGER,                INTENT(IN)   :: KKE            ! uppest atmosphere physical index
INTEGER,                INTENT(IN)   :: KKL            ! +1 if grid goes from ground to atmosphere top, -1 otherwise
INTEGER, DIMENSION(:),  INTENT(IN)   :: KKLCL          ! index of updraft condensation level
REAL, DIMENSION(:,:),   INTENT(IN)   :: PFRAC_UP       ! Updraft Fraction
REAL, DIMENSION(:,:),   INTENT(IN)   :: PRC_UP,PRI_UP  ! updraft characteritics
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PRC_MF, PRI_MF ! cloud content
REAL, DIMENSION(:,:),   INTENT(OUT)  :: PCF_MF         ! and cloud fraction for MF scheme
!
!*                    0.1  Declaration of local variables
!
INTEGER  :: JI,JK
!
!*                    0.2 Initialisation
!
!
!*      1. COMPUTATION OF SUBGRID CLOUD
!          ----------------------------

!
! Warning: updraft variables are on flux levels
! and PRC_MF, PRI_MF and PCF_MF are on mass levels
PRC_MF(:,:)=0.
PRI_MF(:,:)=0.
PCF_MF(:,:)=0.

DO JI=1,SIZE(PCF_MF,1)
  DO JK=KKLCL(JI),KKE-KKL,KKL
    PCF_MF(JI,JK ) = MAX( 0., MIN(1.,XKCF_MF *0.5* (       &
                &    PFRAC_UP(JI,JK) +  PFRAC_UP(JI,JK+KKL) ) ))
    PRC_MF(JI,JK)  = 0.5* XKCF_MF * ( PFRAC_UP(JI,JK)*PRC_UP(JI,JK)  &
                         + PFRAC_UP(JI,JK+KKL)*PRC_UP(JI,JK+KKL) )
    PRI_MF(JI,JK)  = 0.5* XKCF_MF * ( PFRAC_UP(JI,JK)*PRI_UP(JI,JK)  &
                         + PFRAC_UP(JI,JK+KKL)*PRI_UP(JI,JK+KKL) )
  END DO
END DO


END SUBROUTINE COMPUTE_MF_CLOUD_DIRECT
