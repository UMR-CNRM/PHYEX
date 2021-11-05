!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
!#######################################
MODULE MODI_LIMA_COMPUTE_CLOUD_FRACTIONS
!#######################################
  INTERFACE
     SUBROUTINE LIMA_COMPUTE_CLOUD_FRACTIONS (KIB, KIE, KJB, KJE, KKB, KKE, KKL, &
                                              PCCT, PRCT,                        &
                                              PCRT, PRRT,                        &
                                              PCIT, PRIT,                        &
                                              PRST, PRGT, PRHT,                  &
                                              PCLDFR, PICEFR, PPRCFR             )
       INTEGER,               INTENT(IN)    :: KIB           !  
       INTEGER,               INTENT(IN)    :: KIE           !  
       INTEGER,               INTENT(IN)    :: KJB           !  
       INTEGER,               INTENT(IN)    :: KJE           !  
       INTEGER,               INTENT(IN)    :: KKB           !  
       INTEGER,               INTENT(IN)    :: KKE           !  
       INTEGER,               INTENT(IN)    :: KKL           !  
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
       REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRST          !
       REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRGT          !
       REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRHT          !
       !
       REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PCLDFR        ! 
       REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PICEFR        ! 
       REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PPRCFR        ! 
       !
     END SUBROUTINE LIMA_COMPUTE_CLOUD_FRACTIONS
  END INTERFACE
END MODULE MODI_LIMA_COMPUTE_CLOUD_FRACTIONS
!
!
!################################################################
SUBROUTINE LIMA_COMPUTE_CLOUD_FRACTIONS (KIB, KIE, KJB, KJE, KKB, KKE, KKL, &
                                         PCCT, PRCT,                        &
                                         PCRT, PRRT,                        &
                                         PCIT, PRIT,                        &
                                         PRST, PRGT, PRHT,                  &
                                         PCLDFR, PICEFR, PPRCFR             )
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
USE MODD_PARAM_LIMA,      ONLY : XCTMIN, XRTMIN
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
INTEGER,               INTENT(IN)    :: KIB           !  
INTEGER,               INTENT(IN)    :: KIE           !  
INTEGER,               INTENT(IN)    :: KJB           !  
INTEGER,               INTENT(IN)    :: KJE           !  
INTEGER,               INTENT(IN)    :: KKB           !  
INTEGER,               INTENT(IN)    :: KKE           !  
INTEGER,               INTENT(IN)    :: KKL           !  
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
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRST          !
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRGT          !
REAL, DIMENSION(:,:,:),INTENT(IN)    :: PRHT          !
!
REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PCLDFR        ! 
REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PICEFR        ! 
REAL, DIMENSION(:,:,:),INTENT(INOUT) :: PPRCFR        ! 
!
!*       0.2   Declarations of local variables :
!
INTEGER :: JI, JJ, JK
!
!-------------------------------------------------------------------------------
!
!              CLOUD FRACTIONS
!              ---------------
!
! Liquid cloud fraction is kept from input data, except where PCLDFR=0 and rc>0
WHERE(PCLDFR(:,:,:)<1.E-10 .AND. PRCT(:,:,:)>XRTMIN(2) .AND. PCCT(:,:,:)>XCTMIN(2)) PCLDFR(:,:,:)=1.
!
! Ice cloud fraction is currently 0 or 1
PICEFR(:,:,:)=0.
WHERE(PICEFR(:,:,:)<1.E-10 .AND. PRIT(:,:,:)>XRTMIN(4) .AND. PCIT(:,:,:)>XCTMIN(4)) PICEFR(:,:,:)=1.
!
! Precipitation fraction
!!$PPRCFR(:,:,:) = MAX(PCLDFR(:,:,:),PICEFR(:,:,:))
!!$DO JI = KIB,KIE
!!$   DO JJ = KJB, KJE
!!$      DO JK=KKE-KKL, KKB, -KKL
!!$         IF ( (PRRT(JI,JJ,JK).GT.XRTMIN(3) .AND. PCRT(JI,JJ,JK).GT.XCTMIN(3)) .OR. &
!!$               PRST(JI,JJ,JK).GT.XRTMIN(5)                                    .OR. &
!!$               PRGT(JI,JJ,JK).GT.XRTMIN(6)                                    .OR. &
!!$               PRHT(JI,JJ,JK).GT.XRTMIN(7)                                         ) THEN
!!$            PPRCFR(JI,JJ,JK)=MAX(PPRCFR(JI,JJ,JK),PPRCFR(JI,JJ,JK+KKL))
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
!!$DO JI = KIB,KIE
!!$   DO JJ = KJB, KJE
!!$      DO JK=KKE-KKL, KKB, -KKL
!!$         IF ( (PRRT(JI,JJ,JK).GT.0. .AND. PCRT(JI,JJ,JK).GT.0.) .OR. &
!!$               PRST(JI,JJ,JK).GT.0.                             .OR. &
!!$               PRGT(JI,JJ,JK).GT.0.                             .OR. &
!!$               PRHT(JI,JJ,JK).GT.0.                                  ) THEN
!!$            PPRCFR(JI,JJ,JK)=MAX(PPRCFR(JI,JJ,JK),PPRCFR(JI,JJ,JK+KKL))
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
WHERE ( (PRRT(:,:,:).GT.0. .AND. PCRT(:,:,:).GT.0.) .OR. &
         PRST(:,:,:).GT.0.                          .OR. &
         PRGT(:,:,:).GT.0.                          .OR. &
         PRHT(:,:,:).GT.0.                               )  PPRCFR(:,:,:) = 1.
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE LIMA_COMPUTE_CLOUD_FRACTIONS
