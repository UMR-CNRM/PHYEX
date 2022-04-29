MODULE SHUMAN_PHY
IMPLICIT NONE
CONTAINS
!     ###############################
      SUBROUTINE MYF_PHY(D,PA,PMYF)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *MYF* -  Shuman operator : mean operator in y direction for a
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the y direction (J index) for a field PA localized at a y-flux
!     point (v point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------
!!        The result PMYF(i,:,:) is defined by 0.5*(PA(:,j,:)+PA(:,j+1,:))
!!        At j=size(PA,2), PMYF(:,j,:) are replaced by the values of PMYF,
!!    which are the right values in the y-cyclic case
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PMYF   ! result at flux localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! upper bound in y direction of PA
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MYF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MYF',0,ZHOOK_HANDLE)
IJU = SIZE(PA,2)

!POUR AROME

PMYF=PA

!
!DO JJ=1,IJU-1
!  PMYF(:,JJ,:) = 0.5*( PA(:,JJ,:)+PA(:,JJ+1,:) )
!END DO
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MYF',1,ZHOOK_HANDLE)
END SUBROUTINE MYF_PHY
!
!
!     ###############################
      SUBROUTINE MYM_PHY(D,PA,PMYM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *MYM* -  Shuman operator : mean operator in y direction for a
!!                                 mass variable
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the y direction (J index) for a field PA localized at a mass
!     point. The result is localized at a y-flux point (v point).
!
!!**  METHOD
!!    ------
!!        The result PMYM(:,j,:) is defined by 0.5*(PA(:,j,:)+PA(:,j-1,:))
!!    At j=1, PMYM(:,j,:) are replaced by the values of PMYM,
!!    which are the right values in the y-cyclic case.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PMYM   ! result at flux localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MYM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MYM',0,ZHOOK_HANDLE)
IJU=SIZE(PA,2)

!POUR AROME

PMYM=PA

!
!DO JJ=2,IJU
!  PMYM(:,JJ,:) = 0.5*( PA(:,JJ,:)+PA(:,JJ-1,:) )
!END DO
!
!PMYM(:,1,:)    = PMYM(:,IJU-2*JPHEXT+1,:)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MYM',1,ZHOOK_HANDLE)
END SUBROUTINE MYM_PHY
!     ###############################
      SUBROUTINE MZM_PHY(D,PA,PMZM)
!     ###############################
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *MZM* -  Shuman operator : mean operator in z direction for a
!!                                 mass variable
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the z direction (K index) for a field PA localized at a mass
!     point. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------
!!        The result PMZM(:,:,k) is defined by 0.5*(PA(:,:,k)+PA(:,:,k-1))
!!        At k=1, PMZM(:,:,1) is defined by -999.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PMZM   ! result at flux localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK             ! Loop index in z direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MZM',0,ZHOOK_HANDLE)
DO JK=2,SIZE(PA,3)-1
  PMZM(:,:,JK) = 0.5*( PA(:,:,JK)+PA(:,:,JK-D%NKL) )
END DO
PMZM(:,:,D%NKA)    = -999.
PMZM(:,:,D%NKU) = 0.5*( PA(:,:,D%NKU)+PA(:,:,D%NKU-D%NKL) )
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MZM',1,ZHOOK_HANDLE)
END SUBROUTINE MZM_PHY
!     ###############################
      SUBROUTINE DZM_PHY(D,PA,PDZM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *DZM* -  Shuman operator : finite difference operator in z direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference
!     along the z direction (K index) for a field PA localized at a mass
!     point. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------
!!        The result PDZM(:,j,:) is defined by (PA(:,:,k)-PA(:,:,k-1))
!!        At k=1, PDZM(:,:,k) is defined by -999.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDZM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK            ! Loop index in z direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DZM',0,ZHOOK_HANDLE)
DO JK=2,SIZE(PA,3)-1
  PDZM(:,:,JK)          = PA(:,:,JK) -  PA(:,:,JK-D%NKL)
END DO
PDZM(:,:,D%NKA)    =  -999.
PDZM(:,:,D%NKU)    = PA(:,:,D%NKU) -  PA(:,:,D%NKU-D%NKL)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DZM',1,ZHOOK_HANDLE)
END SUBROUTINE DZM_PHY

!     ###############################
      SUBROUTINE MXM_PHY(D,PA,PMXM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *MXM* -  Shuman operator : mean operator in x direction for a
!!                                 mass variable
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the x direction (I index) for a field PA localized at a mass
!     point. The result is localized at a x-flux point (u point).
!
!!**  METHOD
!!    ------
!!        The result PMXM(i,:,:) is defined by 0.5*(PA(i,:,:)+PA(i-1,:,:))
!!    At i=1, PMXM(1,:,:) are replaced by the values of PMXM,
!!    which are the right values in the x-cyclic case.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PMXM   ! result at flux localization
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! Size of the array in the x direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MXM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MXM',0,ZHOOK_HANDLE)
IIU = SIZE(PA,1)
!
!POUR AROME

PMXM=PA

!
!DO JI=2,IIU
!  PMXM(JI,:,:) = 0.5*( PA(JI,:,:)+PA(JI-1,:,:) )
!END DO
!
!PMXM(1,:,:)    = PMXM(IIU-2*JPHEXT+1,:,:)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MXM',1,ZHOOK_HANDLE)
END SUBROUTINE MXM_PHY
!     ###############################
      SUBROUTINE MXF_PHY(D,PA,PMXF)      
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *MXF* -  Shuman operator : mean operator in x direction for a
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the x direction (I index) for a field PA localized at a x-flux
!     point (u point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------
!!        The result PMXF(i,:,:) is defined by 0.5*(PA(i,:,:)+PA(i+1,:,:))
!!        At i=size(PA,1), PMXF(i,:,:) are replaced by the values of PMXF,
!!    which are the right values in the x-cyclic case
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PMXF   ! result at flux localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MXF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MXF',0,ZHOOK_HANDLE)
IIU = SIZE(PA,1)
!
!POUR AROME

PMXF=PA

!
!DO JI=1,IIU-1
!  PMXF(JI,:,:) = 0.5*( PA(JI,:,:)+PA(JI+1,:,:) )
!END DO
!
!PMXF(IIU,:,:)    = PMXF(2*JPHEXT,:,:)
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MXF',1,ZHOOK_HANDLE)
END SUBROUTINE MXF_PHY
!     ###############################
      SUBROUTINE MZF_PHY(D,PA,PMZF)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *MZF* -  Shuman operator : mean operator in z direction for a
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the z direction (K index) for a field PA localized at a z-flux
!     point (w point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------
!!        The result PMZF(:,:,k) is defined by 0.5*(PA(:,:,k)+PA(:,:,k+1))
!!        At k=size(PA,3), PMZF(:,:,k) is defined by -999.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PMZF   ! result at mass localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK             ! Loop index in z direction
INTEGER :: IKT          ! upper bound in z direction of PA
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MZF',0,ZHOOK_HANDLE)
IKT = SIZE(PA,3)
DO JK=2,IKT-1
  PMZF(:,:,JK) = 0.5*( PA(:,:,JK)+PA(:,:,JK+D%NKL) )
END DO
PMZF(:,:,D%NKU) = -999.
PMZF(:,:,D%NKA) = 0.5*( PA(:,:,D%NKA)+PA(:,:,D%NKA+D%NKL) )
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MZF',1,ZHOOK_HANDLE)
END SUBROUTINE MZF_PHY
!     ###############################
      SUBROUTINE DZF_PHY(D,PA,PDZF)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *DZF* -  Shuman operator : finite difference operator in z direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference
!     along the z direction (K index) for a field PA localized at a z-flux
!     point (w point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------
!!        The result PDZF(:,:,k) is defined by (PA(:,:,k+1)-PA(:,:,k))
!!        At k=size(PA,3), PDZF(:,:,k) is defined by -999.
!!
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PDZF   ! result at mass localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK           ! Loop index in z direction
INTEGER :: IKT          ! upper bound in z direction of PA
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DZF',0,ZHOOK_HANDLE)
IKT = SIZE(PA,3)
DO JK=2,IKT-1
  PDZF(:,:,JK)          = PA(:,:,JK+D%NKL) -  PA(:,:,JK)
END DO
PDZF(:,:,D%NKA)    = PA(:,:,D%NKA+D%NKL) -  PA(:,:,D%NKA)
PDZF(:,:,D%NKU)    = -999.
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DZF',1,ZHOOK_HANDLE)
END SUBROUTINE DZF_PHY
END MODULE SHUMAN_PHY
