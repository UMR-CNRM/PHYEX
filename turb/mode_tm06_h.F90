!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1
MODULE MODE_TM06_H
IMPLICIT NONE
CONTAINS
SUBROUTINE TM06_H(D,PTSTEP,PZZ,PFLXZ,PBL_DEPTH)
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
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
ZBL_DEPTH(IIJB:IIJE) = PBL_DEPTH(IIJB:IIJE)
!
DO JIJ=IIJB,IIJE 
  IF(ZBL_DEPTH(JIJ)==XUNDEF)THEN
    ZBL_DEPTH(JIJ)=0.
  ENDIF
ENDDO
!
DO JIJ=IIJB,IIJE 
  PBL_DEPTH(JIJ) = XUNDEF
  ZFLXZMIN(JIJ) = PFLXZ(JIJ,IKB)
ENDDO
!
DO JK=IKTB,IKTE
  DO JIJ=IIJB,IIJE 
    IF(PFLXZ(JIJ,IKB)>0. .AND. PFLXZ(JIJ,JK)<ZFLXZMIN(JIJ))THEN
      PBL_DEPTH(JIJ) = PZZ(JIJ,JK) - PZZ(JIJ,IKB)
      ZFLXZMIN(JIJ) = PFLXZ(JIJ,JK)
    ENDIF
  ENDDO
END DO
!
DO JIJ=IIJB,IIJE 
  IF(PBL_DEPTH(JIJ)/=XUNDEF) THEN
    PBL_DEPTH(JIJ)=MIN(PBL_DEPTH(JIJ),ZBL_DEPTH(JIJ)+ZGROWTH*PTSTEP)
  ENDIF
ENDDO
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('TM06_H',1,ZHOOK_HANDLE)
END SUBROUTINE TM06_H
END MODULE MODE_TM06_H
