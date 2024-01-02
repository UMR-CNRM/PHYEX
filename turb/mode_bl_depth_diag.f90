!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_BL_DEPTH_DIAG
IMPLICIT NONE
!
INTERFACE BL_DEPTH_DIAG  
      MODULE PROCEDURE BL_DEPTH_DIAG_3D
      MODULE PROCEDURE BL_DEPTH_DIAG_1D
END INTERFACE
!
CONTAINS
!
SUBROUTINE BL_DEPTH_DIAG_3D(D,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF,BL_DEPTH_DIAG3D)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),              INTENT(IN)           :: D
REAL, DIMENSION(D%NIJT),       INTENT(IN)           :: PSURF        ! surface flux
REAL, DIMENSION(D%NIJT),       INTENT(IN)           :: PZS          ! orography
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)           :: PFLUX        ! flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)           :: PZZ          ! altitude of flux points
REAL,                          INTENT(IN)           :: PFTOP_O_FSURF! Flux at BL top / Surface flux
REAL, DIMENSION(D%NIJT),       INTENT(OUT)          :: BL_DEPTH_DIAG3D
!
!
!       0.2  declaration of local variables
!
INTEGER :: JIJ,JK ! loop counters
INTEGER :: IKB,IKE,IIJB,IIJE,IKL
REAL    :: ZFLX     ! flux at top of BL
!
!----------------------------------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_3D',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IKL=D%NKL
IIJE=D%NIJE
IIJB=D%NIJB 
!
BL_DEPTH_DIAG3D(:) = 0.
!

DO JIJ=IIJB,IIJE
    IF (PSURF(JIJ)/=0.) THEN
    DO JK=IKB,IKE,IKL
      IF (PZZ(JIJ,JK-IKL)>PZS(JIJ)) THEN
        ZFLX = PSURF(JIJ) * PFTOP_O_FSURF
        IF ( (PFLUX(JIJ,JK)-ZFLX)*(PFLUX(JIJ,JK-IKL)-ZFLX) <= 0. ) THEN
          BL_DEPTH_DIAG3D(JIJ) = (PZZ  (JIJ,JK-IKL) - PZS(JIJ))     &
                         + (PZZ  (JIJ,JK) - PZZ  (JIJ,JK-IKL))    &
                         * (ZFLX          - PFLUX(JIJ,JK-IKL)  )  &
                         / (PFLUX(JIJ,JK) - PFLUX(JIJ,JK-IKL)   )
        END IF
      END IF
    END DO
  END IF
END DO
!
!$mnh_expand_array(JIJ=IIJB:IIJE)
BL_DEPTH_DIAG3D(:) = BL_DEPTH_DIAG3D(:) / (1. - PFTOP_O_FSURF)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_3D',1,ZHOOK_HANDLE)
END SUBROUTINE BL_DEPTH_DIAG_3D
!
SUBROUTINE BL_DEPTH_DIAG_1D(D,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF,BL_DEPTH_DIAG1D)
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),       INTENT(IN)           :: D
REAL,                   INTENT(IN)           :: PSURF        ! surface flux
REAL,                   INTENT(IN)           :: PZS          ! orography
REAL, DIMENSION(D%NKT), INTENT(IN)           :: PFLUX        ! flux
REAL, DIMENSION(D%NKT), INTENT(IN)           :: PZZ          ! altitude of flux points
REAL,                   INTENT(IN)           :: PFTOP_O_FSURF! Flux at BL top / Surface flux
REAL,                   INTENT(OUT)          :: BL_DEPTH_DIAG1D
!
REAL, DIMENSION(1,1)             :: ZSURF
REAL, DIMENSION(1,1)             :: ZZS
REAL, DIMENSION(1,1,D%NKT)       :: ZFLUX
REAL, DIMENSION(1,1,D%NKT)       :: ZZZ
REAL, DIMENSION(1,1)             :: ZBL_DEPTH_DIAG
!
INTEGER :: IKT
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_1D',0,ZHOOK_HANDLE)
IKT=D%NKT
ZSURF        = PSURF
ZZS          = PZS
ZFLUX(1,1,:) = PFLUX(:)
ZZZ  (1,1,:) = PZZ  (:)
!
CALL BL_DEPTH_DIAG_3D(D,ZSURF,ZZS,ZFLUX,ZZZ,PFTOP_O_FSURF,ZBL_DEPTH_DIAG)
!
BL_DEPTH_DIAG1D = ZBL_DEPTH_DIAG(1,1)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_1D',1,ZHOOK_HANDLE)
END SUBROUTINE BL_DEPTH_DIAG_1D
END MODULE MODE_BL_DEPTH_DIAG
