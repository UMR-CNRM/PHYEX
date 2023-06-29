!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_COMPUTE_CLOUD_FRACTIONS
  IMPLICIT NONE
CONTAINS
!################################################################
  SUBROUTINE LIMA_COMPUTE_CLOUD_FRACTIONS (D,                     &
                                           PCCT, PRCT,            &
                                           PCRT, PRRT,            &
                                           PCIT, PRIT,            &
                                           PCST, PRST,            &
                                           PCGT, PRGT,            &
                                           PCHT, PRHT,            &
                                           PCLDFR, PICEFR, PPRCFR )
!################################################################
!
!!
!!    PURPOSE
!!    -------
!!      Compute cloud, ice and precipitating fractions
!!
!!    AUTHOR
!!    ------
!!      B.    Vié        * CNRM *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original             06/03/2019 
!!
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAM_LIMA,      ONLY : XCTMIN, XRTMIN, &
                                 NMOM_C, NMOM_R, NMOM_I, NMOM_S, NMOM_G, NMOM_H
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(DIMPHYEX_t),      INTENT(IN)    :: D
!
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PCCT          !
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRCT          !
!
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PCRT          !
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRRT          !
!
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PCIT          !
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRIT          !
!
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PCST          !
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRST          !
!
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PCGT          !
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRGT          !
!
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PCHT          !
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRHT          !
!
REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PCLDFR        ! 
REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PICEFR        ! 
REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PPRCFR        ! 
!
!*       0.2   Declarations of local variables :
!
!
!-------------------------------------------------------------------------------
!
!              CLOUD FRACTIONS
!              ---------------
!
! Liquid cloud fraction is kept from input data, except where PCLDFR=0 and rc>0
WHERE(PCLDFR(:,:,:)<1.E-10 .AND. PRCT(:,:,:)>XRTMIN(2) .AND. (NMOM_C.EQ.1 .OR. PCCT(:,:,:)>XCTMIN(2))) PCLDFR(:,:,:)=1.
!
! Ice cloud fraction is currently 0 or 1
PICEFR(:,:,:)=0.
WHERE(PICEFR(:,:,:)<1.E-10 .AND. PRIT(:,:,:)>XRTMIN(4) .AND. (NMOM_I.EQ.1 .OR. PCIT(:,:,:)>XCTMIN(4))) PICEFR(:,:,:)=1.
!
! Precipitation fraction
!!$PPRCFR(:,:,:) = MAX(PCLDFR(:,:,:),PICEFR(:,:,:))
!!$DO JI = D%NIB,D%NIE
!!$   DO JJ = D%NJB, D%NJE
!!$      DO JK=D%NKE-D%NKL, D%NKB, -D%NKL
!!$         IF ( (PRRT(JI,JJ,JK).GT.XRTMIN(3) .AND. PCRT(JI,JJ,JK).GT.XCTMIN(3)) .OR. &
!!$               PRST(JI,JJ,JK).GT.XRTMIN(5)                                    .OR. &
!!$               PRGT(JI,JJ,JK).GT.XRTMIN(6)                                    .OR. &
!!$               PRHT(JI,JJ,JK).GT.XRTMIN(7)                                         ) THEN
!!$            PPRCFR(JI,JJ,JK)=MAX(PPRCFR(JI,JJ,JK),PPRCFR(JI,JJ,JK+D%NKL))
!!$            IF (PPRCFR(JI,JJ,JK)==0) THEN
!!$               PPRCFR(JI,JJ,JK)=1.
!!$            END IF
!!$         ELSE
!!$            !PPRCFR(JI,JJ,JK)=0.
!!$         END IF
!!$      END DO
!!$   END DO
!!$END DO
!!$
!!$PPRCFR(:,:,:) = MAX(PCLDFR(:,:,:),PICEFR(:,:,:))
!!$DO JI = D%NIB,D%NIE
!!$   DO JJ = D%NJB, D%NJE
!!$      DO JK=D%NKE-D%NKL, D%NKB, -D%NKL
!!$         IF ( (PRRT(JI,JJ,JK).GT.0. .AND. PCRT(JI,JJ,JK).GT.0.) .OR. &
!!$               PRST(JI,JJ,JK).GT.0.                             .OR. &
!!$               PRGT(JI,JJ,JK).GT.0.                             .OR. &
!!$               PRHT(JI,JJ,JK).GT.0.                                  ) THEN
!!$            PPRCFR(JI,JJ,JK)=MAX(PPRCFR(JI,JJ,JK),PPRCFR(JI,JJ,JK+D%NKL))
!!$            IF (PPRCFR(JI,JJ,JK)==0) THEN
!!$               PPRCFR(JI,JJ,JK)=1.
!!$            END IF
!!$         ELSE
!!$            !PPRCFR(JI,JJ,JK)=0.
!!$         END IF
!!$      END DO
!!$   END DO
!!$END DO
!!$
!!$PPRCFR(:,:,:) = 0.
!!$WHERE ( (PRRT(:,:,:).GT.XRTMIN(3) .AND. PCRT(:,:,:).GT.XCTMIN(3)) .OR. &
!!$         PRST(:,:,:).GT.XRTMIN(5)                                 .OR. &
!!$         PRGT(:,:,:).GT.XRTMIN(6)                                 .OR. &
!!$         PRHT(:,:,:).GT.XRTMIN(7)                                      )  PPRCFR(:,:,:) = 1.
!!$
PPRCFR(:,:,:) = 0.
WHERE ( (PRRT(:,:,:).GT.0. .AND. (NMOM_R.EQ.1 .OR. PCRT(:,:,:).GT.0.) ) .OR. &
        (PRST(:,:,:).GT.0. .AND. (NMOM_S.EQ.1 .OR. PCST(:,:,:).GT.0.) ) .OR. &
        (PRGT(:,:,:).GT.0. .AND. (NMOM_G.EQ.1 .OR. PCGT(:,:,:).GT.0.) ) .OR. &
        (PRHT(:,:,:).GT.0. .AND. (NMOM_H.EQ.1 .OR. PCHT(:,:,:).GT.0.) ) )  PPRCFR(:,:,:) = 1.
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_COMPUTE_CLOUD_FRACTIONS
END MODULE MODE_LIMA_COMPUTE_CLOUD_FRACTIONS
