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
SUBROUTINE BL_DEPTH_DIAG_3D(D,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF,BL_DEPTH_DIAG3D)
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
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),                   INTENT(IN)           :: D
REAL, DIMENSION(D%NIT,D%NJT),       INTENT(IN)           :: PSURF        ! surface flux
REAL, DIMENSION(D%NIT,D%NJT),       INTENT(IN)           :: PZS          ! orography
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)           :: PFLUX        ! flux
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)           :: PZZ          ! altitude of flux points
REAL,                               INTENT(IN)           :: PFTOP_O_FSURF! Flux at BL top / Surface flux
REAL, DIMENSION(D%NIT,D%NJT),       INTENT(OUT)          :: BL_DEPTH_DIAG3D
!
!
!       0.2  declaration of local variables
!
INTEGER :: JI,JJ,JK ! loop counters
INTEGER :: IKB,IKE,IIB,IIE,IJB,IJE   ! index value for the Beginning
REAL    :: ZFLX     ! flux at top of BL
!
!----------------------------------------------------------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_3D',0,ZHOOK_HANDLE)
IKB=D%NKTB
IKE=D%NKTE
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
BL_DEPTH_DIAG3D(:,:) = 0.
!

DO JJ=1,IJE
  DO JI=1,IIE
    IF (PSURF(JI,JJ)==0.) CYCLE
    DO JK=IKB,IKE,D%NKL
      IF (PZZ(JI,JJ,JK-D%NKL)<=PZS(JI,JJ)) CYCLE
      ZFLX = PSURF(JI,JJ) * PFTOP_O_FSURF
      IF ( (PFLUX(JI,JJ,JK)-ZFLX)*(PFLUX(JI,JJ,JK-D%NKL)-ZFLX) <= 0. ) THEN
        BL_DEPTH_DIAG3D(JI,JJ) = (PZZ  (JI,JJ,JK-D%NKL) - PZS(JI,JJ))     &
                         + (PZZ  (JI,JJ,JK) - PZZ  (JI,JJ,JK-D%NKL))    &
                         * (ZFLX            - PFLUX(JI,JJ,JK-D%NKL)  )  &
                         / (PFLUX(JI,JJ,JK) - PFLUX(JI,JJ,JK-D%NKL)   )
        EXIT
      END IF
    END DO
  END DO
END DO
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
BL_DEPTH_DIAG3D(IIB:IIE,IJB:IJE) = BL_DEPTH_DIAG3D(IIB:IIE,IJB:IJE) / (1. - PFTOP_O_FSURF)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_3D',1,ZHOOK_HANDLE)
END SUBROUTINE BL_DEPTH_DIAG_3D
!
SUBROUTINE BL_DEPTH_DIAG_1D(D,PSURF,PZS,PFLUX,PZZ,PFTOP_O_FSURF,BL_DEPTH_DIAG1D)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('BL_DEPTH_DIAG_1D',0,ZHOOK_HANDLE)
ZSURF        = PSURF
ZZS          = PZS
ZFLUX(1,1,1:D%NKT) = PFLUX(1:D%NKT)
ZZZ  (1,1,1:D%NKT) = PZZ  (1:D%NKT)
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
