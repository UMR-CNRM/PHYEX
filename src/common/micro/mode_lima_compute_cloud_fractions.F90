!MNH_LIC Copyright 2019-2021 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-------------------------------------------------------------------------------
MODULE MODE_LIMA_COMPUTE_CLOUD_FRACTIONS
  IMPLICIT NONE
CONTAINS
!################################################################
  SUBROUTINE LIMA_COMPUTE_CLOUD_FRACTIONS (LIMAP, D,                     &
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
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_T
USE MODD_PARAM_LIMA, ONLY:PARAM_LIMA_T
USE YOMHOOK, ONLY:LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!
!*       0.1   Declarations of dummy arguments :
!
TYPE(PARAM_LIMA_T),INTENT(IN)::LIMAP
TYPE(DIMPHYEX_T),      INTENT(IN)    :: D
!
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PCCT          !
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PRCT          !
!
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PCRT          !
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PRRT          !
!
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PCIT          !
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PRIT          !
!
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PCST          !
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PRST          !
!
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PCGT          !
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PRGT          !
!
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PCHT          !
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(IN)    :: PRHT          !
!
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(INOUT) :: PCLDFR        ! 
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(INOUT) :: PICEFR        ! 
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(INOUT) :: PPRCFR        ! 
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
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
IF (LHOOK) CALL DR_HOOK('LIMA_COMPUTE_CLOUD_FRACTIONS', 0, ZHOOK_HANDLE)
WHERE(PCLDFR(:,:)<1.E-10 .AND. PRCT(:,:)>LIMAP%XRTMIN(2) .AND. (LIMAP%NMOM_C.EQ.1 .OR. PCCT(:,:)>LIMAP%XCTMIN(2))) PCLDFR(:,:)=1.
!
! Ice cloud fraction is currently 0 or 1
PICEFR(:,:)=0.
WHERE(PICEFR(:,:)<1.E-10 .AND. PRIT(:,:)>LIMAP%XRTMIN(4) .AND. (LIMAP%NMOM_I.EQ.1 .OR. PCIT(:,:)>LIMAP%XCTMIN(4))) PICEFR(:,:)=1.
!
! Precipitation fraction
!!$PPRCFR(:,:) = MAX(PCLDFR(:,:),PICEFR(:,:))
!!$DO JI = D%NIJB,D%NIJE
!!$      DO JK=D%NKE-D%NKL, D%NKB, -D%NKL
!!$         IF ( (PRRT(JI,JK).GT.LIMAP%XRTMIN(3) .AND. PCRT(JI,JK).GT.LIMAP%XCTMIN(3)) .OR. &
!!$               PRST(JI,JK).GT.LIMAP%XRTMIN(5)                                    .OR. &
!!$               PRGT(JI,JK).GT.LIMAP%XRTMIN(6)                                    .OR. &
!!$               PRHT(JI,JK).GT.LIMAP%XRTMIN(7)                                         ) THEN
!!$            PPRCFR(JI,JK)=MAX(PPRCFR(JI,JK),PPRCFR(JI,JK+D%NKL))
!!$            IF (PPRCFR(JI,JK)==0) THEN
!!$               PPRCFR(JI,JK)=1.
!!$            END IF
!!$         ELSE
!!$            !PPRCFR(JI,JK)=0.
!!$         END IF
!!$      END DO
!!$END DO
!!$
!!$PPRCFR(:,:) = MAX(PCLDFR(:,:),PICEFR(:,:))
!!$DO JI = D%NIJB,D%NIJE
!!$      DO JK=D%NKE-D%NKL, D%NKB, -D%NKL
!!$         IF ( (PRRT(JI,JK).GT.0. .AND. PCRT(JI,JK).GT.0.) .OR. &
!!$               PRST(JI,JK).GT.0.                             .OR. &
!!$               PRGT(JI,JK).GT.0.                             .OR. &
!!$               PRHT(JI,JK).GT.0.                                  ) THEN
!!$            PPRCFR(JI,JK)=MAX(PPRCFR(JI,JK),PPRCFR(JI,JK+D%NKL))
!!$            IF (PPRCFR(JI,JK)==0) THEN
!!$               PPRCFR(JI,JK)=1.
!!$            END IF
!!$         ELSE
!!$            !PPRCFR(JI,JK)=0.
!!$         END IF
!!$      END DO
!!$END DO
!!$
!!$PPRCFR(:,:) = 0.
!!$WHERE ( (PRRT(:,:).GT.LIMAP%XRTMIN(3) .AND. PCRT(:,:).GT.LIMAP%XCTMIN(3)) .OR. &
!!$         PRST(:,:).GT.LIMAP%XRTMIN(5)                                 .OR. &
!!$         PRGT(:,:).GT.LIMAP%XRTMIN(6)                                 .OR. &
!!$         PRHT(:,:).GT.LIMAP%XRTMIN(7)                                      )  PPRCFR(:,:) = 1.
!!$
PPRCFR(:,:) = 0.
WHERE ( (PRRT(:,:).GT.0. .AND. (LIMAP%NMOM_R.EQ.1 .OR. PCRT(:,:).GT.0.) ) .OR. &
        (PRST(:,:).GT.0. .AND. (LIMAP%NMOM_S.EQ.1 .OR. PCST(:,:).GT.0.) ) .OR. &
        (PRGT(:,:).GT.0. .AND. (LIMAP%NMOM_G.EQ.1 .OR. PCGT(:,:).GT.0.) ) .OR. &
        (PRHT(:,:).GT.0. .AND. (LIMAP%NMOM_H.EQ.1 .OR. PCHT(:,:).GT.0.) ) )  PPRCFR(:,:) = 1.
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('LIMA_COMPUTE_CLOUD_FRACTIONS', 1, ZHOOK_HANDLE)
END SUBROUTINE LIMA_COMPUTE_CLOUD_FRACTIONS
END MODULE MODE_LIMA_COMPUTE_CLOUD_FRACTIONS
