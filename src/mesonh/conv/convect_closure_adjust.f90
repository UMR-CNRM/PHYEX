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
      MODULE MODI_CONVECT_CLOSURE_ADJUST
!     #################
!
INTERFACE
!
       SUBROUTINE CONVECT_CLOSURE_ADJUST( KLON, KLEV, PADJ,                      &
                                        PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR, &
                                        PDMF, PZDMF, PDER, PZDER, PDDR, PZDDR, &
                                        PPRMELT, PZPRMELT, PDTEVR, PZDTEVR,    &
                                        PTPR, PZTPR,                           &
                                        PPRLFLX, PZPRLFL, PPRSFLX, PZPRSFL     )
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
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF  ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDMF ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER  ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDER ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR  ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDDR ! initial value of  "
REAL, DIMENSION(KLON),   INTENT(INOUT):: PTPR     ! total precipitation (kg/s)
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZTPR    ! initial value of "
REAL, DIMENSION(KLON),   INTENT(INOUT):: PDTEVR   ! donwndraft evapor. (kg/s)
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZDTEVR  ! initial value of " 
REAL, DIMENSION(KLON),   INTENT(INOUT):: PPRMELT  ! melting of precipitation
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZPRMELT ! initial value of " 
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRLFLX! liquid precip flux
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRLFL! initial value "
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRSFLX! solid  precip flux
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRSFL! initial value "
!
END SUBROUTINE CONVECT_CLOSURE_ADJUST
!
END INTERFACE
!
END MODULE MODI_CONVECT_CLOSURE_ADJUST
!    ###########################################################################
     SUBROUTINE CONVECT_CLOSURE_ADJUST( KLON, KLEV, PADJ,                      &
                                        PUMF, PZUMF, PUER, PZUER, PUDR, PZUDR, &
                                        PDMF, PZDMF, PDER, PZDER, PDDR, PZDDR, &
                                        PPRMELT, PZPRMELT, PDTEVR, PZDTEVR,    &
                                        PTPR, PZTPR,                           &
                                        PPRLFLX, PZPRLFL, PPRSFLX, PZPRSFL     )
!    ###########################################################################
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
!!   Last modified  04/10/97
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
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDMF  ! downdraft mass flux (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDMF ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDER  ! downdraft entrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDER ! initial value of  "
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PDDR  ! downdraft detrainment (kg/s)
REAL, DIMENSION(KLON,KLEV), INTENT(INOUT) :: PZDDR ! initial value of  "
REAL, DIMENSION(KLON),   INTENT(INOUT):: PTPR     ! total precipitation (kg/s)
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZTPR    ! initial value of "
REAL, DIMENSION(KLON),   INTENT(INOUT):: PDTEVR   ! donwndraft evapor. (kg/s)
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZDTEVR  ! initial value of " 
REAL, DIMENSION(KLON),   INTENT(INOUT):: PPRMELT  ! melting of precipitation
REAL, DIMENSION(KLON),   INTENT(INOUT):: PZPRMELT ! initial value of " 
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRLFLX! liquid precip flux
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRLFL! initial value "
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PPRSFLX! solid  precip flux
REAL, DIMENSION(KLON,KLEV),INTENT(INOUT)  :: PZPRSFL! initial value "
!
!
!*       0.2   Declarations of local variables :
!
INTEGER :: IKB, IKE                 !  vert. loop bounds
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
          PPRMELT(:)  = PZPRMELT(:)   * PADJ(:)
          PDTEVR(:)   = PZDTEVR(:)    * PADJ(:)
          PTPR(:)     = PZTPR(:)      * PADJ(:)
!
     DO JK = IKB + 1, IKE
	  PUMF(:,JK)  = PZUMF(:,JK)   * PADJ(:)
          PUER(:,JK)  = PZUER(:,JK)   * PADJ(:)
          PUDR(:,JK)  = PZUDR(:,JK)   * PADJ(:)
          PDMF(:,JK)  = PZDMF(:,JK)   * PADJ(:)
          PDER(:,JK)  = PZDER(:,JK)   * PADJ(:)
          PDDR(:,JK)  = PZDDR(:,JK)   * PADJ(:)
          PPRLFLX(:,JK) = PZPRLFL(:,JK) * PADJ(:)
          PPRSFLX(:,JK) = PZPRSFL(:,JK) * PADJ(:)
     END DO
!
END SUBROUTINE CONVECT_CLOSURE_ADJUST
