!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1
MODULE MODE_TM06_H
IMPLICIT NONE
CONTAINS
SUBROUTINE TM06_H(KKB,KKTB,KKTE,PTSTEP,PZZ,PFLXZ,PBL_DEPTH)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
!
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
INTEGER,                INTENT(IN)    :: KKB       ! index of 1st physical level
                                                   ! close to ground 
INTEGER,                INTENT(IN)    :: KKTB      ! first physical level in k
INTEGER,                INTENT(IN)    :: KKTE      ! last physical level in k
REAL,                   INTENT(IN)    :: PTSTEP    ! Double time step
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ       ! altitude of flux levels
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFLXZ     ! heat flux
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PBL_DEPTH ! boundary layer height
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER                                  :: JK     ! loop counter
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZFLXZMIN ! minimum of temperature flux 
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZBL_DEPTH! BL depth at previous time-step
REAL                                     :: ZGROWTH  ! maximum BL growth rate
!----------------------------------------------------------------------------
!
!* mixed boundary layer cannot grow more rapidly than 1800m/h
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('TM06_H',0,ZHOOK_HANDLE)
ZGROWTH = 2.0 ! (m/s)
!
!----------------------------------------------------------------------------
!
ZBL_DEPTH(:,:) = PBL_DEPTH(:,:)
WHERE(ZBL_DEPTH(:,:)==XUNDEF) ZBL_DEPTH(:,:)=0.
!
PBL_DEPTH(:,:) = XUNDEF
ZFLXZMIN (:,:) = PFLXZ(:,:,KKB)
!
DO JK=KKTB,KKTE
  WHERE (PFLXZ(:,:,KKB)>0. .AND. PFLXZ(:,:,JK)<ZFLXZMIN(:,:))
    PBL_DEPTH(:,:) = PZZ  (:,:,JK) - PZZ(:,:,KKB)
    ZFLXZMIN (:,:) = PFLXZ(:,:,JK)
  END WHERE
END DO
!
WHERE(PBL_DEPTH(:,:)/=XUNDEF) PBL_DEPTH(:,:)=MIN(PBL_DEPTH(:,:),ZBL_DEPTH(:,:)+ZGROWTH*PTSTEP)
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM06_H',1,ZHOOK_HANDLE)
END SUBROUTINE TM06_H
END MODULE MODE_TM06_H
