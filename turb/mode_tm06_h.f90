!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1
MODULE MODE_TM06_H
IMPLICIT NONE
CONTAINS
SUBROUTINE TM06_H(D,PTSTEP,PZZ,PFLXZ,PBL_DEPTH)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #################################################################
!
!
!!****  *TM06_H* - computes the Third Order Moments
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
USE MODD_PARAMETERS, ONLY : XUNDEF
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
REAL,                   INTENT(IN)    :: PTSTEP    ! Double time step
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PZZ       ! altitude of flux levels
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PFLXZ     ! heat flux
REAL, DIMENSION(D%NIJT),   INTENT(INOUT) :: PBL_DEPTH ! boundary layer height
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER                                  :: JK,JIJ     ! loop counter
INTEGER :: IKB,IKTB,IKTE,IIJB,IIJE
REAL, DIMENSION(D%NIJT) :: ZFLXZMIN ! minimum of temperature flux 
REAL, DIMENSION(D%NIJT) :: ZBL_DEPTH! BL depth at previous time-step
REAL                                     :: ZGROWTH  ! maximum BL growth rate
!----------------------------------------------------------------------------
!
!* mixed boundary layer cannot grow more rapidly than 1800m/h
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TM06_H',0,ZHOOK_HANDLE)
ZGROWTH = 2.0 ! (m/s)
!
!----------------------------------------------------------------------------
!
IKTB=D%NKTB          
IKTE=D%NKTE
IKB=D%NKB
IIJE=D%NIJE
IIJB=D%NIJB

!
ZBL_DEPTH(:) = PBL_DEPTH(:)
!
!$mnh_expand_where(JIJ=IIJB:IIJE)
WHERE(ZBL_DEPTH(:)==XUNDEF)
  ZBL_DEPTH(:)=0.
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE)
!
!$mnh_expand_array(JIJ=IIJB:IIJE)
PBL_DEPTH(:) = XUNDEF
ZFLXZMIN(:) = PFLXZ(:,IKB)
!$mnh_end_expand_array(JIJ=IIJB:IIJE)
!
DO JK=IKTB,IKTE
!$mnh_expand_where(JIJ=IIJB:IIJE)
  WHERE(PFLXZ(:,IKB)>0. .AND. PFLXZ(:,JK)<ZFLXZMIN(:))
    PBL_DEPTH(:) = PZZ(:,JK) - PZZ(:,IKB)
    ZFLXZMIN(:) = PFLXZ(:,JK)
  END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE)
END DO
!
!$mnh_expand_where(JIJ=IIJB:IIJE)
WHERE(PBL_DEPTH(:)/=XUNDEF) 
  PBL_DEPTH(:)=MIN(PBL_DEPTH(:),ZBL_DEPTH(:)+ZGROWTH*PTSTEP)
END WHERE
!$mnh_end_expand_where(JIJ=IIJB:IIJE)
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM06_H',1,ZHOOK_HANDLE)
END SUBROUTINE TM06_H
END MODULE MODE_TM06_H
