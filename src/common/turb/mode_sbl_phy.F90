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
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PZ_O_LMO
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: BUSINGERPHIM
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JI,JJ,JK
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM',0,ZHOOK_HANDLE)
!$mnh_expand_where(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
WHERE ( PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) < 0. )
  BUSINGERPHIM(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) = (1.-15.*PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT))**(-0.25)
ELSEWHERE
  BUSINGERPHIM(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) = 1. + 4.7 * PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT)
END WHERE
!$mnh_end_expand_where(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PZ_O_LMO
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: BUSINGERPHIH
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JI,JJ,JK
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH',0,ZHOOK_HANDLE)
!$mnh_expand_where(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
WHERE ( PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) < 0. )
  BUSINGERPHIH(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) = 0.74 * (1.-9.*PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT))**(-0.5)
ELSEWHERE
  BUSINGERPHIH(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) = 0.74 + 4.7 * PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT)
END WHERE
!$mnh_end_expand_where(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)   :: PZ_O_LMO
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT)  :: BUSINGERPHIE
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JI,JJ,JK
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIE',0,ZHOOK_HANDLE)
!$mnh_expand_where(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
WHERE ( PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) < 0. )
  BUSINGERPHIE(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT)=(1.+(-PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT))**(2./3.)/CSTURB%XALPSBL)&
                            * (1.-15.*PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT))**(0.5)
ELSEWHERE
  BUSINGERPHIE(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) = 1./(1. + 4.7 * PZ_O_LMO(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT))**2
END WHERE
!$mnh_end_expand_where(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIE',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIE
END MODULE MODE_SBL_PHY
