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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PZZ       ! altitude of flux levels
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PFLXU     ! u'w'
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PFLXV     ! v'w'
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)    :: PWTHV     ! buoyancy flux
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(IN)    :: PLMO      ! Monin-Obukhov length
REAL, DIMENSION(D%NIT,D%NJT),   INTENT(INOUT) :: PSBL_DEPTH! boundary layer height
!
!-------------------------------------------------------------------------------
!
!       0.2  declaration of local variables
!
!
INTEGER                                  :: JLOOP,JI,JJ,JK    ! loop counter
INTEGER :: IKB,IKE,IIB,IIE,IJB,IJE   ! index value for the Beginning
REAL, DIMENSION(D%NIT,D%NJT) :: ZQ0      ! surface buoyancy flux
REAL, DIMENSION(D%NIT,D%NJT) :: ZWU      ! surface friction u'w'
REAL, DIMENSION(D%NIT,D%NJT) :: ZWV      ! surface friction v'w'
REAL, DIMENSION(D%NIT,D%NJT) :: ZUSTAR2  ! surface friction
REAL, DIMENSION(D%NIT,D%NJT) :: ZSBL_DYN ! SBL wih dynamical criteria
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWIND
                                         ! intermediate wind for SBL calculation
REAL, DIMENSION(D%NIT,D%NJT) :: ZSBL_THER! SBL wih thermal   criteria
REAL, DIMENSION(D%NIT,D%NJT) :: ZA       ! ponderation coefficient
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
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZWU(IIB:IIE,IJB:IJE) = PFLXU(IIB:IIE,IJB:IJE,IKB)
ZWV(IIB:IIE,IJB:IJE) = PFLXV(IIB:IIE,IJB:IJE,IKB)
ZQ0(IIB:IIE,IJB:IJE) = PWTHV(IIB:IIE,IJB:IJE,IKB)
!
ZUSTAR2(IIB:IIE,IJB:IJE) = SQRT(ZWU(IIB:IIE,IJB:IJE)**2+ZWV(IIB:IIE,IJB:IJE)**2)
!
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!----------------------------------------------------------------------------
!
!* BL and SBL diagnosed with friction criteria
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
ZWIND(IIB:IIE,IJB:IJE,:)=SQRT(PFLXU(IIB:IIE,IJB:IJE,:)**2+PFLXV(IIB:IIE,IJB:IJE,:)**2)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
CALL BL_DEPTH_DIAG(D,ZUSTAR2,PZZ(:,:,IKB),ZWIND,PZZ,CSTURB%XFTOP_O_FSURF,ZSBL_DYN)
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZSBL_DYN(IIB:IIE,IJB:IJE) = CSTURB%XSBL_O_BL * ZSBL_DYN(IIB:IIE,IJB:IJE)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!----------------------------------------------------------------------------
!
!* BL and SBL diagnosed with buoyancy flux criteria
!
CALL BL_DEPTH_DIAG(D,ZQ0,PZZ(:,:,IKB),PWTHV,PZZ,CSTURB%XFTOP_O_FSURF,ZSBL_THER)
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
ZSBL_THER(IIB:IIE,IJB:IJE)= CSTURB%XSBL_O_BL * ZSBL_THER(IIB:IIE,IJB:IJE)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!----------------------------------------------------------------------------
!
!* SBL depth
!
PSBL_DEPTH(:,:) = 0.
!$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
WHERE (ZSBL_THER(IIB:IIE,IJB:IJE)> 0. .AND. ZSBL_DYN(IIB:IIE,IJB:IJE)> 0.) 
  PSBL_DEPTH = MIN(ZSBL_THER(IIB:IIE,IJB:IJE),ZSBL_DYN(IIB:IIE,IJB:IJE))
END WHERE
!$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
!
!$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
WHERE (ZSBL_THER(IIB:IIE,IJB:IJE)> 0. .AND. ZSBL_DYN(IIB:IIE,IJB:IJE)==0.) 
  PSBL_DEPTH(IIB:IIE,IJB:IJE) = ZSBL_THER(IIB:IIE,IJB:IJE)
END WHERE
!$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
!
!$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
WHERE (ZSBL_THER(IIB:IIE,IJB:IJE)==0. .AND. ZSBL_DYN(IIB:IIE,IJB:IJE)> 0.) 
  PSBL_DEPTH(IIB:IIE,IJB:IJE) = ZSBL_DYN(IIB:IIE,IJB:IJE)
END WHERE
!$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
!
DO JLOOP=1,5
  !$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
  WHERE (PLMO(IIB:IIE,IJB:IJE)/=XUNDEF .AND. ABS(PLMO(IIB:IIE,IJB:IJE))>=0.01 )
    ZA(IIB:IIE,IJB:IJE) = TANH(2.*PSBL_DEPTH(IIB:IIE,IJB:IJE)/PLMO(IIB:IIE,IJB:IJE))**2
    PSBL_DEPTH(IIB:IIE,IJB:IJE) = 0.2 * PSBL_DEPTH(IIB:IIE,IJB:IJE) + 0.8 * ((1.-ZA(IIB:IIE,IJB:IJE)) &
                                * ZSBL_DYN(IIB:IIE,IJB:IJE) + ZA(IIB:IIE,IJB:IJE) * ZSBL_THER(IIB:IIE,IJB:IJE) )
  END WHERE
  !$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
END DO
!$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
WHERE (ABS(PLMO(IIB:IIE,IJB:IJE))<=0.01 ) 
  PSBL_DEPTH(IIB:IIE,IJB:IJE) = ZSBL_THER(IIB:IIE,IJB:IJE)
END WHERE
!$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
!$mnh_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
WHERE (PLMO(IIB:IIE,IJB:IJE)==XUNDEF)
  PSBL_DEPTH(IIB:IIE,IJB:IJE) = ZSBL_DYN(IIB:IIE,IJB:IJE)
END WHERE
!$mnh_end_expand_where(JI=IIB:IIE,JJ=IJB:IJE)
!
!----------------------------------------------------------------------------
IF (LHOOK) CALL DR_HOOK('SBL_DEPTH',1,ZHOOK_HANDLE)
END SUBROUTINE SBL_DEPTH
END MODULE MODE_SBL_DEPTH
