!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_TM06
IMPLICIT NONE
CONTAINS
SUBROUTINE TM06(D,CST,PTHVREF,PBL_DEPTH,PZZ,PSFTH,PMWTH,PMTH2)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #################################################################
!
!
!!****  *TM06* - computes the Third Order Moments
!!
!!    PURPOSE
!!    -------
!
!!**  METHOD
!!    ------
!!
!!     TOMs are deduced from convective normalized TOMs according to Tomas and
!!     Masson 2006
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
!!      V. MAsson and S. Tomas  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         sept. 2005
!!
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CST,        ONLY: CST_t
USE MODD_DIMPHYEX,   ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: XUNDEF
!
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN) :: D
TYPE(CST_t),            INTENT(IN) :: CST
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PTHVREF    ! reference potential temperature
REAL, DIMENSION(D%NIJT),   INTENT(IN) :: PBL_DEPTH ! boundary layer height
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN) :: PZZ        ! altitude of flux levels
REAL, DIMENSION(D%NIJT),   INTENT(IN) :: PSFTH      ! surface heat flux
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT):: PMWTH      ! w'2th'
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT):: PMTH2      ! w'th'2
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
REAL, DIMENSION(D%NIJT,D%NKT):: ZZ_O_H ! normalized height z/h (where h=BL height)
REAL, DIMENSION(D%NIJT)            :: ZWSTAR ! normalized convective velocity w*
REAL, DIMENSION(D%NIJT)            :: ZTSTAR ! normalized temperature velocity w*
!
INTEGER                                             :: JK,JIJ     ! loop counter
INTEGER                                             :: IIJE,IIJB
INTEGER                                             :: IKTB,IKTE,IKB,IKE,IKT,IKU ! vertical levels
!----------------------------------------------------------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TM06',0,ZHOOK_HANDLE)
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IKT=D%NKT
IKE=D%NKE
IIJE=D%NIJE
IIJB=D%NIJB
IKU=D%NKU
!
!
!* w* and T*
!
DO JIJ=IIJB, IIJE
  IF (PSFTH(JIJ)>0.) THEN
    ZWSTAR(JIJ) = ((CST%XG/PTHVREF(JIJ, IKB))*PSFTH(JIJ)*PBL_DEPTH(JIJ))**(1./3.)
    ZTSTAR(JIJ) = PSFTH(JIJ) / ZWSTAR(JIJ)
  ELSE
    ZWSTAR(JIJ) = 0.
    ZTSTAR(JIJ) = 0.
  END IF
END DO
!
!
!* normalized height
!
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZZ_O_H(JIJ, JK) = XUNDEF
  END DO
END DO
DO JK=1,IKT
  DO JIJ=IIJB, IIJE
    IF (PBL_DEPTH(JIJ)/=XUNDEF) THEN
      ZZ_O_H(JIJ, JK) = (PZZ(JIJ, JK)-PZZ(JIJ, IKB)) / PBL_DEPTH(JIJ)
    END IF
  END DO
END DO
!
!* w'th'2
!
PMTH2(IIJB:IIJE,1:IKT) = 0.
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    IF (ZZ_O_H(JIJ, JK) < 0.95 .AND. ZZ_O_H(JIJ, JK)/=XUNDEF) THEN
      PMTH2(JIJ, JK) = 4.*(MAX(ZZ_O_H(JIJ, JK),0.))**0.4*(ZZ_O_H(JIJ, JK)-0.95)**2
    END IF
  END DO
END DO
DO JK=IKTB+1,IKTE-1
  DO JIJ=IIJB, IIJE
    PMTH2(JIJ, JK) = PMTH2(JIJ, JK) * ZTSTAR(JIJ)**2*ZWSTAR(JIJ)
  END DO
END DO
DO JIJ=IIJB, IIJE
  PMTH2(JIJ, IKE)=PMTH2(JIJ, IKE) * ZTSTAR(JIJ)**2*ZWSTAR(JIJ)
  PMTH2(JIJ, IKU)=PMTH2(JIJ, IKU) * ZTSTAR(JIJ)**2*ZWSTAR(JIJ)
END DO
!
!
!* w'2th'
!
PMWTH(IIJB:IIJE,1:IKT) = 0.
DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    IF (ZZ_O_H(JIJ, JK) <0.9 .AND. ZZ_O_H(JIJ, JK)/=XUNDEF) THEN
      PMWTH(JIJ, JK) = MAX(-7.9*(ABS(ZZ_O_H(JIJ, JK)-0.35))**2.9 &
                               * (ABS(ZZ_O_H(JIJ, JK)-1.))**0.58 + 0.37, 0.)
    END IF
  END DO
END DO

DO JK=IKTB+1,IKTE-1
  DO JIJ=IIJB, IIJE
    PMWTH(JIJ, JK) = PMWTH(JIJ, JK) * ZWSTAR(JIJ)**2*ZTSTAR(JIJ)
  END DO
END DO
DO JIJ=IIJB, IIJE
  PMWTH(JIJ, IKE) = PMWTH(JIJ, IKE) * ZWSTAR(JIJ)**2*ZTSTAR(JIJ)
  PMWTH(JIJ, IKU) = PMWTH(JIJ, IKU) * ZWSTAR(JIJ)**2*ZTSTAR(JIJ)
END DO
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM06',1,ZHOOK_HANDLE)
END SUBROUTINE TM06
END MODULE MODE_TM06
