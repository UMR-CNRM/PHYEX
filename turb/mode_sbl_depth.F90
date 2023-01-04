!MNH_LIC Copyright 1994-2022 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
MODULE MODE_SBL_DEPTH
IMPLICIT NONE
CONTAINS
!     ######spl
      SUBROUTINE SBL_DEPTH(D,CSTURB,PZZ,PFLXU,PFLXV,PWTHV,PLMO,PSBL_DEPTH)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #################################################################
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
!!      V. Masson  * Meteo-France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original         nov. 2005
!!      26/02/2020   T.Nagel Correction of SBL depth computation in neutral stratification
!! --------------------------------------------------------------------------
!       
!*      0. DECLARATIONS
!          ------------
!
USE MODD_CTURB, ONLY: CSTURB_t
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY : XUNDEF
!
USE MODE_BL_DEPTH_DIAG, ONLY : BL_DEPTH_DIAG
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
TYPE(DIMPHYEX_t),       INTENT(IN)    :: D
TYPE(CSTURB_t),         INTENT(IN)    :: CSTURB
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PZZ       ! altitude of flux levels
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PFLXU     ! u'w'
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PFLXV     ! v'w'
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)    :: PWTHV     ! buoyancy flux
REAL, DIMENSION(D%NIJT),   INTENT(IN)    :: PLMO      ! Monin-Obukhov length
REAL, DIMENSION(D%NIJT),   INTENT(INOUT) :: PSBL_DEPTH! boundary layer height
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER                                  :: JLOOP,JIJ,JK    ! loop counter
INTEGER :: IKB,IKE,IIJB,IIJE,IKT   ! index value for the Beginning
REAL, DIMENSION(D%NIJT) :: ZQ0      ! surface buoyancy flux
REAL, DIMENSION(D%NIJT) :: ZWU      ! surface friction u'w'
REAL, DIMENSION(D%NIJT) :: ZWV      ! surface friction v'w'
REAL, DIMENSION(D%NIJT) :: ZUSTAR2  ! surface friction
REAL, DIMENSION(D%NIJT) :: ZSBL_DYN ! SBL wih dynamical criteria
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWIND
                                         ! intermediate wind for SBL calculation
REAL, DIMENSION(D%NIJT) :: ZSBL_THER! SBL wih thermal   criteria
REAL, DIMENSION(D%NIJT) :: ZA       ! ponderation coefficient
!----------------------------------------------------------------------------
!
!* initialisations
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('SBL_DEPTH',0,ZHOOK_HANDLE)
!
IKB=D%NKTB
IKE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
!
DO JIJ=IIJB,IIJE 
  ZWU(JIJ) = PFLXU(JIJ,IKB)
  ZWV(JIJ) = PFLXV(JIJ,IKB)
  ZQ0(JIJ) = PWTHV(JIJ,IKB)
!
  ZUSTAR2(JIJ) = SQRT(ZWU(JIJ)**2+ZWV(JIJ)**2)
!
ENDDO
!----------------------------------------------------------------------------
!
!* BL and SBL diagnosed with friction criteria
!
DO JK=1,IKT 
  DO JIJ=IIJB,IIJE 
    ZWIND(JIJ,JK)=SQRT(PFLXU(JIJ,JK)**2+PFLXV(JIJ,JK)**2)
  ENDDO
ENDDO
CALL BL_DEPTH_DIAG(D,ZUSTAR2,PZZ(:,IKB),ZWIND,PZZ,CSTURB%XFTOP_O_FSURF,ZSBL_DYN)
DO JIJ=IIJB,IIJE 
  ZSBL_DYN(JIJ) = CSTURB%XSBL_O_BL * ZSBL_DYN(JIJ)
ENDDO
!
!----------------------------------------------------------------------------
!
!* BL and SBL diagnosed with buoyancy flux criteria
!
CALL BL_DEPTH_DIAG(D,ZQ0,PZZ(:,IKB),PWTHV,PZZ,CSTURB%XFTOP_O_FSURF,ZSBL_THER)
DO JIJ=IIJB,IIJE 
  ZSBL_THER(JIJ)= CSTURB%XSBL_O_BL * ZSBL_THER(JIJ)
ENDDO
!
!----------------------------------------------------------------------------
!
!* SBL depth
!
PSBL_DEPTH(:) = 0.
DO JIJ=IIJB,IIJE 
  IF (ZSBL_THER(JIJ)> 0. .AND. ZSBL_DYN(JIJ)> 0.) THEN
    PSBL_DEPTH = MIN(ZSBL_THER(JIJ),ZSBL_DYN(JIJ))
  ENDIF
ENDDO
!
DO JIJ=IIJB,IIJE 
  IF (ZSBL_THER(JIJ)> 0. .AND. ZSBL_DYN(JIJ)==0.) THEN
    PSBL_DEPTH(JIJ) = ZSBL_THER(JIJ)
  ENDIF
ENDDO
!
DO JIJ=IIJB,IIJE 
  IF (ZSBL_THER(JIJ)==0. .AND. ZSBL_DYN(JIJ)> 0.) THEN
    PSBL_DEPTH(JIJ) = ZSBL_DYN(JIJ)
  ENDIF
ENDDO
!
DO JLOOP=1,5
  DO JIJ=IIJB,IIJE 
    IF (PLMO(JIJ)/=XUNDEF .AND. ABS(PLMO(JIJ))>=0.01 )THEN
      ZA(JIJ) = TANH(2.*PSBL_DEPTH(JIJ)/PLMO(JIJ))**2
      PSBL_DEPTH(JIJ) = 0.2 * PSBL_DEPTH(JIJ) + 0.8 * ((1.-ZA(JIJ)) &
      * ZSBL_DYN(JIJ) + ZA(JIJ) * ZSBL_THER(JIJ) )
    ENDIF
  ENDDO
END DO
DO JIJ=IIJB,IIJE 
  IF (ABS(PLMO(JIJ))<=0.01 ) THEN
    PSBL_DEPTH(JIJ) = ZSBL_THER(JIJ)
  ENDIF
ENDDO
DO JIJ=IIJB,IIJE 
  IF (PLMO(JIJ)==XUNDEF)THEN
    PSBL_DEPTH(JIJ) = ZSBL_DYN(JIJ)
  ENDIF
ENDDO
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SBL_DEPTH',1,ZHOOK_HANDLE)
END SUBROUTINE SBL_DEPTH
END MODULE MODE_SBL_DEPTH
