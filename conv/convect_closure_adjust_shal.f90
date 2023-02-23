!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!--------------- special set of characters for RCS information
!-----------------------------------------------------------------
! $Source$ $Revision$
! MASDEV4_7 conv 2006/05/18 13:07:25
!-----------------------------------------------------------------
!     #################
      MODULE MODI_CONVECT_CLOSURE_ADJUST_SHAL
!     #################
!
INTERFACE 
!
       SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, PADJ,                      &
                                             PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR  )
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
REAL, DIMENSION(KLON),      INTENT(IN) :: PADJ     ! mass adjustment factor
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUMF ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUER ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUDR ! initial value of  "
!
END SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL
!
END INTERFACE
!
END MODULE MODI_CONVECT_CLOSURE_ADJUST_SHAL
!    ################################################################################
     SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL( KLON, KLEV, PADJ,                      &
                                             PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR  )
!    ################################################################################
!
!!**** Uses closure adjustment factor to adjust mass flux and to modify
!!     precipitation efficiency  when necessary. The computations are
!!     similar to routine CONVECT_PRECIP_ADJUST.
!!
!!
!!    PURPOSE
!!    -------
!!      The purpose of this routine is to adjust the mass flux using the
!!      factor PADJ computed in CONVECT_CLOSURE
!!
!!
!!**  METHOD
!!    ------
!!      Computations are done at every model level starting from bottom.
!!      The use of masks allows to optimise the inner loops (horizontal loops).
!!      
!!
!!    EXTERNAL
!!    --------
!!     Module MODD_CONVPAREXT
!!          JCVEXB, JCVEXT     ! extra levels on the vertical boundaries
!!     
!!    None
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    None
!!
!!    REFERENCE
!!    ---------
!!
!!      Book1,2 of documentation ( routine CONVECT_CLOSURE_ADJUST)
!!
!!    AUTHOR
!!    ------
!!      P. BECHTOLD       * Laboratoire d'Aerologie *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    26/03/96 
!!   Last modified  15/11/96
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONVPAREXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
!
INTEGER,                    INTENT(IN) :: KLON     ! horizontal dimension
INTEGER,                    INTENT(IN) :: KLEV     ! vertical dimension
REAL, DIMENSION(KLON),      INTENT(IN) :: PADJ     ! mass adjustment factor
!
!
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUMF  ! updraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUMF ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUER  ! updraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUER ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PUDR  ! updraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZUDR ! initial value of  "
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB, IKE                 ! vert. loop bounds
INTEGER :: JK                       ! vertical loop index
!
!
!-------------------------------------------------------------------------------
!
!*       0.3   Compute loop bounds
!              -------------------
!
IKB  = 1 + JCVEXB 
IKE  = KLEV - JCVEXT
!
!
!*       1.     Adjust mass flux by the factor PADJ to converge to
!               specified degree of stabilization
!               ----------------------------------------------------
!
     DO JK = IKB + 1, IKE
	  PUMF(:,JK)  = PZUMF(:,JK)   * PADJ(:)
          PUER(:,JK)  = PZUER(:,JK)   * PADJ(:)
          PUDR(:,JK)  = PZUDR(:,JK)   * PADJ(:)
     END DO
!
END SUBROUTINE CONVECT_CLOSURE_ADJUST_SHAL
