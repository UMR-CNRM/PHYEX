!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_BL_DEPTH_DIAG
!
INTERFACE BL_DEPTH_DIAG  
      MODULE PROCEDURE BL_DEPTH_DIAG_3D
      MODULE PROCEDURE BL_DEPTH_DIAG_1D
END INTERFACE
!
CONTAINS
!
FUNCTION BL_DEPTH_DIAG_3D(KKB,KKE,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
!
!!****  *SBL_DEPTH* - computes SBL depth
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!    SBL is defined as the layer where momentum flux is equal to XSBL_FRAC of its surface value
!!    
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!
!!    REFERENCE
!!    ---------
!!
!!    AUTHOR
!!    ------
!!      V. Masson * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         nov. 2005
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
!*      0.1  declarations of arguments
!
IMPLICIT NONE
!
INTEGER,                INTENT(IN)           :: KKB          ! bottom point
INTEGER,                INTENT(IN)           :: KKE          ! top point
REAL, DIMENSION(:,:),   INTENT(IN)           :: PSURF        ! surface flux
REAL, DIMENSION(:,:),   INTENT(IN)           :: PZS          ! orography
REAL, DIMENSION(:,:,:), INTENT(IN)           :: PFLUX        ! flux
REAL, DIMENSION(:,:,:), INTENT(IN)           :: PZZ          ! altitude of flux points
REAL,                   INTENT(IN)           :: PFTOP_O_FSURF! Flux at BL top / Surface flux
REAL, DIMENSION(SIZE(PSURF,1),SIZE(PSURF,2)) :: BL_DEPTH_DIAG_3D
!
!
!       0.2  declaration of local variables
!
INTEGER :: JI,JJ,JK ! loop counters
INTEGER :: IKL      ! +1 : MesoNH levels -1: Arome
REAL    :: ZFLX     ! flux at top of BL
!
!----------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_3D',0,ZHOOK_HANDLE)
IF (KKB < KKE) THEN
  IKL=1
ELSE
  IKL=-1
ENDIF

BL_DEPTH_DIAG_3D(:,:) = 0.
!

DO JJ=1,SIZE(PSURF,2)
  DO JI=1,SIZE(PSURF,1)
    IF (PSURF(JI,JJ)==0.) CYCLE
    DO JK=KKB,KKE,IKL
      IF (PZZ(JI,JJ,JK-IKL)<=PZS(JI,JJ)) CYCLE
      ZFLX = PSURF(JI,JJ) * PFTOP_O_FSURF
      IF ( (PFLUX(JI,JJ,JK)-ZFLX)*(PFLUX(JI,JJ,JK-IKL)-ZFLX) <= 0. ) THEN
        BL_DEPTH_DIAG_3D(JI,JJ) = (PZZ  (JI,JJ,JK-IKL) - PZS(JI,JJ))     &
                         + (PZZ  (JI,JJ,JK) - PZZ  (JI,JJ,JK-IKL))    &
                         * (ZFLX            - PFLUX(JI,JJ,JK-IKL)  )  &
                         / (PFLUX(JI,JJ,JK) - PFLUX(JI,JJ,JK-IKL)   )
        EXIT
      END IF
    END DO
  END DO
END DO
!
BL_DEPTH_DIAG_3D(:,:) = BL_DEPTH_DIAG_3D(:,:) / (1. - PFTOP_O_FSURF)
!
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_3D',1,ZHOOK_HANDLE)
END FUNCTION BL_DEPTH_DIAG_3D
!
FUNCTION BL_DEPTH_DIAG_1D(KKB,KKE,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!
IMPLICIT NONE
!
INTEGER,                INTENT(IN)           :: KKB          ! bottom point
INTEGER,                INTENT(IN)           :: KKE          ! top point
REAL,                   INTENT(IN)           :: PSURF        ! surface flux
REAL,                   INTENT(IN)           :: PZS          ! orography
REAL, DIMENSION(:),     INTENT(IN)           :: PFLUX        ! flux
REAL, DIMENSION(:),     INTENT(IN)           :: PZZ          ! altitude of flux points
REAL,                   INTENT(IN)           :: PFTOP_O_FSURF! Flux at BL top / Surface flux
REAL                                         :: BL_DEPTH_DIAG_1D
!
REAL, DIMENSION(1,1)             :: ZSURF
REAL, DIMENSION(1,1)             :: ZZS
REAL, DIMENSION(1,1,SIZE(PFLUX)) :: ZFLUX
REAL, DIMENSION(1,1,SIZE(PZZ))   :: ZZZ
REAL, DIMENSION(1,1)             :: ZBL_DEPTH_DIAG
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_1D',0,ZHOOK_HANDLE)
ZSURF        = PSURF
ZZS          = PZS
ZFLUX(1,1,:) = PFLUX(:)
ZZZ  (1,1,:) = PZZ  (:)
!
ZBL_DEPTH_DIAG = BL_DEPTH_DIAG_3D(KKB,KKE,ZSURF,ZZS,ZFLUX,ZZZ,PFTOP_O_FSURF)
!
BL_DEPTH_DIAG_1D = ZBL_DEPTH_DIAG(1,1)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_1D',1,ZHOOK_HANDLE)
END FUNCTION BL_DEPTH_DIAG_1D
END MODULE MODE_BL_DEPTH_DIAG
