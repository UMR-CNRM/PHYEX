!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_TM06
IMPLICIT NONE
CONTAINS
SUBROUTINE TM06(D,CST,PTHVREF,PBL_DEPTH,PZZ,PSFTH,PMWTH,PMTH2)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
USE MODD_PARAMETERS, ONLY: XUNDEF,JPVEXT_TURB
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
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
!$mnh_expand_where(JIJ=IIJB:IIJE)
WHERE(PSFTH(IIJB:IIJE)>0.)
  ZWSTAR(IIJB:IIJE) = ((CST%XG/PTHVREF(IIJB:IIJE,IKB))*PSFTH(IIJB:IIJE)*PBL_DEPTH(IIJB:IIJE))**(1./3.)
  ZTSTAR(IIJB:IIJE) = PSFTH(IIJB:IIJE) / ZWSTAR(IIJB:IIJE)
ELSEWHERE
  ZWSTAR(IIJB:IIJE) = 0.
  ZTSTAR(IIJB:IIJE) = 0.
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE)
!
!
!* normalized height
!
!$mnh_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
ZZ_O_H(IIJB:IIJE,1:IKT) = XUNDEF
!$mnh_end_expand_array(JIJ=IIJB:IIJE,JK=1:IKT)
DO JK=1,IKT
  !$mnh_expand_where(JIJ=IIJB:IIJE)
  WHERE (PBL_DEPTH(IIJB:IIJE)/=XUNDEF)
    ZZ_O_H(IIJB:IIJE,JK) = (PZZ(IIJB:IIJE,JK)-PZZ(IIJB:IIJE,IKB)) / PBL_DEPTH(IIJB:IIJE)
  END WHERE
  !$mnh_end_expand_where(JIJ=IIJB:IIJE)
END DO
!
!* w'th'2
!
PMTH2(IIJB:IIJE,1:IKT) = 0.
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE(ZZ_O_H(IIJB:IIJE,1:IKT) < 0.95 .AND. ZZ_O_H(IIJB:IIJE,1:IKT)/=XUNDEF)
  PMTH2(IIJB:IIJE,1:IKT) = 4.*(MAX(ZZ_O_H(IIJB:IIJE,1:IKT),0.))**0.4*(ZZ_O_H(IIJB:IIJE,1:IKT)-0.95)**2
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
DO JK=IKTB+1,IKTE-1
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  PMTH2(IIJB:IIJE,JK) = PMTH2(IIJB:IIJE,JK) * ZTSTAR(IIJB:IIJE)**2*ZWSTAR(IIJB:IIJE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END DO
!$mnh_expand_array(JIJ=IIJB:IIJE)
PMTH2(IIJB:IIJE,IKE)=PMTH2(IIJB:IIJE,IKE) * ZTSTAR(IIJB:IIJE)**2*ZWSTAR(IIJB:IIJE)
PMTH2(IIJB:IIJE,IKU)=PMTH2(IIJB:IIJE,IKU) * ZTSTAR(IIJB:IIJE)**2*ZWSTAR(IIJB:IIJE)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!
!* w'2th'
!
PMWTH(IIJB:IIJE,1:IKT) = 0.
!$mnh_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)
WHERE(ZZ_O_H(IIJB:IIJE,1:IKT) <0.9 .AND. ZZ_O_H(IIJB:IIJE,1:IKT)/=XUNDEF)
  PMWTH(IIJB:IIJE,1:IKT) = MAX(-7.9*(ABS(ZZ_O_H(IIJB:IIJE,1:IKT)-0.35))**2.9 &
                           * (ABS(ZZ_O_H(IIJB:IIJE,1:IKT)-1.))**0.58 + 0.37, 0.)
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE,JK=1:IKT)

DO JK=IKTB+1,IKTE-1
  !$mnh_expand_array(JIJ=IIJB:IIJE)
  PMWTH(IIJB:IIJE,JK) = PMWTH(IIJB:IIJE,JK) * ZWSTAR(IIJB:IIJE)**2*ZTSTAR(IIJB:IIJE)
  !$mnh_end_expand_array(JIJ=IIJB:IIJE)
END DO
!$mnh_expand_array(JIJ=IIJB:IIJE)
PMWTH(IIJB:IIJE,IKE) = PMWTH(IIJB:IIJE,IKE) * ZWSTAR(IIJB:IIJE)**2*ZTSTAR(IIJB:IIJE)
PMWTH(IIJB:IIJE,IKU) = PMWTH(IIJB:IIJE,IKU) * ZWSTAR(IIJB:IIJE)**2*ZTSTAR(IIJB:IIJE)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM06',1,ZHOOK_HANDLE)
END SUBROUTINE TM06
END MODULE MODE_TM06
