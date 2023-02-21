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
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PZ_O_LMO
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: BUSINGERPHIM
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JIJ,JK,IIJB,IIJE,IKT
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIM',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE ( PZ_O_LMO(IIJB:IIJE,1:IKT) < 0. )
  BUSINGERPHIM(IIJB:IIJE,1:IKT) = (1.-15.*PZ_O_LMO(IIJB:IIJE,1:IKT))**(-0.25)
ELSEWHERE
  BUSINGERPHIM(IIJB:IIJE,1:IKT) = 1. + 4.7 * PZ_O_LMO(IIJB:IIJE,1:IKT)
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JIJ,JK,IIJB,IIJE,IKT
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIH',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE ( PZ_O_LMO(IIJB:IIJE,1:IKT) < 0. )
  BUSINGERPHIH(IIJB:IIJE,1:IKT) = 0.74 * (1.-9.*PZ_O_LMO(IIJB:IIJE,1:IKT))**(-0.5)
ELSEWHERE
  BUSINGERPHIH(IIJB:IIJE,1:IKT) = 0.74 + 4.7 * PZ_O_LMO(IIJB:IIJE,1:IKT)
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
INTEGER :: JIJ,JK,IIJB,IIJE,IKT
!
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIE',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE ( PZ_O_LMO(IIJB:IIJE,1:IKT) < 0. )
  BUSINGERPHIE(IIJB:IIJE,1:IKT)=(1.+(-PZ_O_LMO(IIJB:IIJE,1:IKT))**(2./3.)/CSTURB%XALPSBL)&
                            * (1.-15.*PZ_O_LMO(IIJB:IIJE,1:IKT))**(0.5)
ELSEWHERE
  BUSINGERPHIE(IIJB:IIJE,1:IKT) = 1./(1. + 4.7 * PZ_O_LMO(IIJB:IIJE,1:IKT))**2
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
IF (LHOOK) CALL DR_HOOK('MODE_SBL:BUSINGER_PHIE',1,ZHOOK_HANDLE)
END SUBROUTINE BUSINGER_PHIE
!
SUBROUTINE LMO(D,CST,PUSTAR,PTHETA,PRV,PSFTH,PSFRV,PLMO)
  USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
  USE MODD_CST, ONLY: CST_t
  USE MODD_PARAMETERS, ONLY: JPVEXT_TURB,XUNDEF
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
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO',0,ZHOOK_HANDLE)
!
  IIJE=D%NIJE
  IIJB=D%NIJB
  IKT=D%NKT
  ZEPS=(CST%XRV-CST%XRD)/CST%XRD
!
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  ZTHETAV(IIJB:IIJE) = PTHETA(IIJB:IIJE) * ( 1. +ZEPS * PRV(IIJB:IIJE))
  ZQ0(IIJB:IIJE) = PSFTH(IIJB:IIJE) + ZTHETAV(IIJB:IIJE) * ZEPS * PSFRV(IIJB:IIJE)
!
  PLMO(IIJB:IIJE) = XUNDEF
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  WHERE ( ZQ0(IIJB:IIJE)/=0. )
    PLMO(IIJB:IIJE) = - MAX(PUSTAR(IIJB:IIJE),1.E-6)**3                &
                  / ( CST%XKARMAN * CST%XG / ZTHETAV(IIJB:IIJE) *ZQ0(IIJB:IIJE) )
  END WHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
  IF (LHOOK) CALL DR_HOOK('MODE_SBL:LMO',1,ZHOOK_HANDLE)
END SUBROUTINE LMO
!
END MODULE MODE_SBL_PHY
