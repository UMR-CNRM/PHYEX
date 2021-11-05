!MNH_LIC Copyright 1994-2014 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt  
!MNH_LIC for details. version 1.
!    ################ 
     MODULE MODI_SBL_DEPTH  
!    ################ 
!
INTERFACE
!
      SUBROUTINE SBL_DEPTH(KKB,KKE,PZZ,PFLXU,PFLXV,PWTHV,PLMO,PSBL_DEPTH)
!
INTEGER,                INTENT(IN)    :: KKB       ! first physical level
INTEGER,                INTENT(IN)    :: KKE       ! upper physical level
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ       ! altitude of flux levels
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFLXU     ! u'w'
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFLXV     ! v'w'
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PWTHV     ! buoyancy flux
REAL, DIMENSION(:,:),   INTENT(IN)    :: PLMO      ! Monin-Obukhov length
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSBL_DEPTH! boundary layer height
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE SBL_DEPTH
!
END INTERFACE
!
END MODULE MODI_SBL_DEPTH
!
!     #################################################################
      SUBROUTINE SBL_DEPTH(KKB,KKE,PZZ,PFLXU,PFLXV,PWTHV,PLMO,PSBL_DEPTH)
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
USE MODD_PARAMETERS, ONLY : XUNDEF
USE MODD_CTURB,      ONLY : XFTOP_O_FSURF, XSBL_O_BL
!
USE MODI_BL_DEPTH_DIAG
!
IMPLICIT NONE
!
!
!*      0.1  declarations of arguments
!
INTEGER,                INTENT(IN)    :: KKB       ! first physical level
INTEGER,                INTENT(IN)    :: KKE       ! upper physical level
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PZZ       ! altitude of flux levels
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFLXU     ! u'w'
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PFLXV     ! v'w'
REAL, DIMENSION(:,:,:), INTENT(IN)    :: PWTHV     ! buoyancy flux
REAL, DIMENSION(:,:),   INTENT(IN)    :: PLMO      ! Monin-Obukhov length
REAL, DIMENSION(:,:),   INTENT(INOUT) :: PSBL_DEPTH! boundary layer height
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER                                  :: JLOOP    ! loop counter
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZQ0      ! surface buoyancy flux
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZWU      ! surface friction u'w'
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZWV      ! surface friction v'w'
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZUSTAR2  ! surface friction
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZSBL_DYN ! SBL wih dynamical criteria
REAL, DIMENSION(SIZE(PFLXU,1),SIZE(PFLXU,2),SIZE(PFLXU,3)) :: ZWIND
                                         ! intermediate wind for SBL calculation
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZSBL_THER! SBL wih thermal   criteria
REAL, DIMENSION(SIZE(PZZ,1),SIZE(PZZ,2)) :: ZA       ! ponderation coefficient
!----------------------------------------------------------------------------
!
!* initialisations
!
!
ZWU (:,:) = PFLXU(:,:,KKB)
ZWV (:,:) = PFLXV(:,:,KKB)
ZQ0 (:,:) = PWTHV(:,:,KKB)
!
ZUSTAR2(:,:) = SQRT(ZWU**2+ZWV**2)
!
!----------------------------------------------------------------------------
!
!* BL and SBL diagnosed with friction criteria
!
ZWIND=SQRT(PFLXU**2+PFLXV**2)
ZSBL_DYN = XSBL_O_BL * BL_DEPTH_DIAG(KKB,KKE,ZUSTAR2,PZZ(:,:,KKB),ZWIND,PZZ,XFTOP_O_FSURF) 
!
!----------------------------------------------------------------------------
!
!* BL and SBL diagnosed with buoyancy flux criteria
!
ZSBL_THER= XSBL_O_BL * BL_DEPTH_DIAG(KKB,KKE,ZQ0,PZZ(:,:,KKB),PWTHV,PZZ,XFTOP_O_FSURF)
!
!----------------------------------------------------------------------------
!
!* SBL depth
!
PSBL_DEPTH = 0.
WHERE (ZSBL_THER> 0. .AND. ZSBL_DYN> 0.) PSBL_DEPTH = MIN(ZSBL_THER(:,:),ZSBL_DYN(:,:))
WHERE (ZSBL_THER> 0. .AND. ZSBL_DYN==0.) PSBL_DEPTH = ZSBL_THER(:,:)
WHERE (ZSBL_THER==0. .AND. ZSBL_DYN> 0.) PSBL_DEPTH = ZSBL_DYN(:,:)
!
DO JLOOP=1,5
  WHERE (PLMO(:,:)/=XUNDEF .AND. ABS(PLMO(:,:))>=0.01 )
    ZA = TANH(2.*PSBL_DEPTH/PLMO)**2
    PSBL_DEPTH = 0.2 * PSBL_DEPTH + 0.8 * ((1.-ZA) * ZSBL_DYN + ZA * ZSBL_THER )
  END WHERE
END DO
WHERE (ABS(PLMO(:,:))<=0.01 )     PSBL_DEPTH = ZSBL_THER
WHERE (PLMO(:,:)==XUNDEF) PSBL_DEPTH = ZSBL_DYN
!
!----------------------------------------------------------------------------
END SUBROUTINE SBL_DEPTH
