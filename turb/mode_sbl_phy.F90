!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ###############
      MODULE MODE_SBL_PHY
!     ###############
!
!!****  *MODE_SBL * - contains Surface Boundary Layer characteristics functions
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!
!!
!!    EXTERNAL
!!    --------
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!
!!    REFERENCE
!!    ---------
!!
!!    Businger et al 1971,   Wyngaard and Cote 1974
!!
!!
!!    AUTHOR
!!    ------
!!      V. Masson       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original        13/10/99
!!      V. Masson       06/11/02 optimization and add Businger fonction for TKE
!!      V. Masson       01/01/03 use PAULSON_PSIM function
!-----------------------------------------------------------------------------
USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!
IMPLICIT NONE
!-------------------------------------------------------------------------------
CONTAINS
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIM(D,PZ_O_LMO,BUSINGERPHIM)
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),                   INTENT(IN)  :: D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PZ_O_LMO
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: BUSINGERPHIM
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ,JK,IIJB,IIJE,IKT
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    IF ( PZ_O_LMO(JIJ, JK) < 0. ) THEN
      BUSINGERPHIM(JIJ, JK) = (1.-15.*PZ_O_LMO(JIJ, JK))**(-0.25)
    ELSE
      BUSINGERPHIM(JIJ, JK) = 1. + 4.7 * PZ_O_LMO(JIJ, JK)
    END IF
  END DO
END DO
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIM
!
!-------------------------------------------------------------------------------
!
SUBROUTINE BUSINGER_PHIH(D,PZ_O_LMO,BUSINGERPHIH)
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),                   INTENT(IN)  :: D
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PZ_O_LMO
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: BUSINGERPHIH
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ,JK,IIJB,IIJE,IKT
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    IF ( PZ_O_LMO(JIJ, JK) < 0. ) THEN
      BUSINGERPHIH(JIJ, JK) = 0.74 * (1.-9.*PZ_O_LMO(JIJ, JK))**(-0.5)
    ELSE
      BUSINGERPHIH(JIJ, JK) = 0.74 + 4.7 * PZ_O_LMO(JIJ, JK)
    END IF
  END DO
END DO
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIH
!
!-------------------------------------------------------------------------------
SUBROUTINE BUSINGER_PHIE(D,CSTURB,PZ_O_LMO,BUSINGERPHIE)
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_CTURB, ONLY: CSTURB_t
!
IMPLICIT NONE
!
TYPE(DIMPHYEX_t),                   INTENT(IN)   :: D
TYPE(CSTURB_t),                     INTENT(IN)   :: CSTURB
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)   :: PZ_O_LMO
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT)  :: BUSINGERPHIE
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
INTEGER :: JIJ,JK,IIJB,IIJE,IKT
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIE',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    IF ( PZ_O_LMO(JIJ, JK) < 0. ) THEN
      BUSINGERPHIE(JIJ, JK)=(1.+(-PZ_O_LMO(JIJ, JK))**(2./3.)/CSTURB%XALPSBL)&
                                * (1.-15.*PZ_O_LMO(JIJ, JK))**(0.5)
    ELSE
      BUSINGERPHIE(JIJ, JK) = 1./(1. + 4.7 * PZ_O_LMO(JIJ, JK))**2
    END IF
  END DO
END DO
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIE',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIE
!
SUBROUTINE LMO(D,CST,PUSTAR,PTHETA,PRV,PSFTH,PSFRV,PLMO)
  USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
  USE MODD_CST, ONLY: CST_t
  USE MODD_PARAMETERS, ONLY: XUNDEF
  !
  TYPE(DIMPHYEX_t),        INTENT(IN)  :: D
  TYPE(CST_t),             INTENT(IN)  :: CST
  REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PUSTAR
  REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PTHETA
  REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PRV
  REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PSFTH
  REAL, DIMENSION(D%NIJT), INTENT(IN)  :: PSFRV
  REAL, DIMENSION(D%NIJT),INTENT(OUT)   :: PLMO
!
  REAL, DIMENSION(D%NIJT)   :: ZTHETAV, ZQ0
  REAL                            :: ZEPS
  INTEGER :: IIJB,IIJE, JIJ,IKT
!
  REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO',0,ZHOOK_HANDLE)
!
  IIJE=D%NIJE
  IIJB=D%NIJB
  IKT=D%NKT
  ZEPS=(CST%XRV-CST%XRD)/CST%XRD
!
  DO JIJ=IIJB, IIJE
    ZTHETAV(JIJ) = PTHETA(JIJ) * ( 1. +ZEPS * PRV(JIJ))
    ZQ0(JIJ) = PSFTH(JIJ) + ZTHETAV(JIJ) * ZEPS * PSFRV(JIJ)
  !
    PLMO(JIJ) = XUNDEF
  END DO
  DO JIJ=IIJB, IIJE
    IF ( ZQ0(JIJ)/=0. ) THEN
      PLMO(JIJ) = - MAX(PUSTAR(JIJ),1.E-6)**3                &
                    / ( CST%XKARMAN * CST%XG / ZTHETAV(JIJ) *ZQ0(JIJ) )
    END IF
  END DO
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO',1,ZHOOK_HANDLE)
END SUBROUTINE LMO
!
END MODULE MODE_SBL_PHY
