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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMYF   ! result at flux localization 
!
! 1.    DEFINITION OF MYF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MYF',0,ZHOOK_HANDLE)

!POUR AROME
!
PMYF=PA
!
IF (LHOOK) CALL DR_HOOK('MYF',1,ZHOOK_HANDLE)
END SUBROUTINE MYF_PHY
!
!     ###############################
      SUBROUTINE MYF2D_PHY(D,PA,PMYF)
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
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PMYF   ! result at flux localization 
!
! 1.    DEFINITION OF MYF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MYF',0,ZHOOK_HANDLE)

!POUR AROME
!
PMYF=PA
!
IF (LHOOK) CALL DR_HOOK('MYF',1,ZHOOK_HANDLE)
END SUBROUTINE MYF2D_PHY
!
!     ###############################
      SUBROUTINE MYM2D_PHY(D,PA,PMYM)
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
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PMYM   ! result at flux localization 
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MYM',0,ZHOOK_HANDLE)

!POUR AROME

PMYM=PA

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MYM',1,ZHOOK_HANDLE)
END SUBROUTINE MYM2D_PHY
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMYM   ! result at flux localization 
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
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PMZM   ! result at flux localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK,JIJ,IIJB,IIJE             ! Loop index
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MZM',0,ZHOOK_HANDLE)
IIJB = D%NIJB
IIJE = D%NIJE
DO JK=2,D%NKT-1
  DO JIJ=IIJB,IIJE 
    PMZM(JIJ,JK) = 0.5*( PA(JIJ,JK)+PA(JIJ,JK-D%NKL) )
  ENDDO
END DO
DO JIJ=IIJB,IIJE 
  PMZM(JIJ,D%NKA)    = -999.
  PMZM(JIJ,D%NKU) = 0.5*( PA(JIJ,D%NKU)+PA(JIJ,D%NKU-D%NKL) )
ENDDO
!
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
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PDZM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK,JIJ,IIJB,IIJE             ! Loop index
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DZM',0,ZHOOK_HANDLE)
IIJB = D%NIJB
IIJE = D%NIJE
DO JK=2,D%NKT-1
  DO JIJ=IIJB,IIJE 
    PDZM(JIJ,JK)          = PA(JIJ,JK) -  PA(JIJ,JK-D%NKL)
  ENDDO
END DO
DO JIJ=IIJB,IIJE 
  PDZM(JIJ,D%NKA)    =  -999.
  PDZM(JIJ,D%NKU)    = PA(JIJ,D%NKU) -  PA(JIJ,D%NKU-D%NKL)
ENDDO
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMXM   ! result at flux localization
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
!
!     ###############################
      SUBROUTINE MXM2D_PHY(D,PA,PMXM)
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
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PMXM   ! result at flux localization
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MXM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MXM',0,ZHOOK_HANDLE)

!POUR AROME

PMXM=PA

IF (LHOOK) CALL DR_HOOK('MXM',1,ZHOOK_HANDLE)
END SUBROUTINE MXM2D_PHY
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMXF   ! result at flux localization 
!
!*       1.    DEFINITION OF MXF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MXF',0,ZHOOK_HANDLE)
!POUR AROME
!
PMXF=PA
!
IF (LHOOK) CALL DR_HOOK('MXF',1,ZHOOK_HANDLE)
END SUBROUTINE MXF_PHY
!     ###############################
      SUBROUTINE MXF2D_PHY(D,PA,PMXF)      
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
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PMXF   ! result at flux localization 
!
!*       1.    DEFINITION OF MXF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MXF',0,ZHOOK_HANDLE)
!POUR AROME
!
PMXF=PA
!
IF (LHOOK) CALL DR_HOOK('MXF',1,ZHOOK_HANDLE)
END SUBROUTINE MXF2D_PHY
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
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PMZF   ! result at mass localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK,JIJ,IIJB,IIJE             ! Loop index
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MZF',0,ZHOOK_HANDLE)
IIJB = D%NIJB
IIJE = D%NIJE
DO JK=2,D%NKT-1
  DO JIJ=IIJB,IIJE 
    PMZF(JIJ,JK) = 0.5*( PA(JIJ,JK)+PA(JIJ,JK+D%NKL) )
  ENDDO
END DO
DO JIJ=IIJB,IIJE 
  PMZF(JIJ,D%NKU) = -999.
  PMZF(JIJ,D%NKA) = 0.5*( PA(JIJ,D%NKA)+PA(JIJ,D%NKA+D%NKL) )
ENDDO
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
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PDZF   ! result at mass localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK,JIJ,IIJB,IIJE             ! Loop index
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DZF',0,ZHOOK_HANDLE)
IIJB = D%NIJB
IIJE = D%NIJE
DO JK=2,D%NKT-1
  DO JIJ=IIJB,IIJE 
    PDZF(JIJ,JK)          = PA(JIJ,JK+D%NKL) -  PA(JIJ,JK)
  ENDDO
END DO
DO JIJ=IIJB,IIJE 
  PDZF(JIJ,D%NKA)    = PA(JIJ,D%NKA+D%NKL) -  PA(JIJ,D%NKA)
  PDZF(JIJ,D%NKU)    = -999.
ENDDO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DZF',1,ZHOOK_HANDLE)
END SUBROUTINE DZF_PHY
!
!     ###############################
      SUBROUTINE DYM_PHY(D,PA,PDYM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *DYM* -  Shuman operator : finite difference operator in y direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference
!     along the y direction (J index) for a field PA localized at a mass
!     point. The result is localized at a y-flux point (v point).
!
!!**  METHOD
!!    ------
!!        The result PDYM(:,j,:) is defined by (PA(:,j,:)-PA(:,j-1,:))
!!    At j=1, PDYM(:,1,:) are replaced by the values of PDYM,
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
!!      Original    05/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDYM     ! result at flux
                                                            ! side
!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DYM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DYM',0,ZHOOK_HANDLE)
IJU=SIZE(PA,2)
!
DO JJ=2,IJU
  PDYM(:,JJ,:)          = PA(:,JJ,:) -  PA(:,JJ-1,:)
END DO
!
PDYM(:,1,:)    =  PDYM(:,IJU-2*JPHEXT+1,:)
CALL ABORT ! AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DYM',1,ZHOOK_HANDLE)
END SUBROUTINE DYM_PHY
!
!     ###############################
      SUBROUTINE DXM_PHY(D,PA,PDXM)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *DXM* -  Shuman operator : finite difference operator in x direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference
!     along the x direction (I index) for a field PA localized at a mass
!     point. The result is localized at a x-flux point (u point).
!
!!**  METHOD
!!    ------
!!        The result PDXM(i,:,:) is defined by (PA(i,:,:)-PA(i-1,:,:))
!!    At i=1, PDXM(1,:,:) are replaced by the values of PDXM,
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
!!      Original    05/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDXM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! Size of the array in the x direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DXM
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DXM',0,ZHOOK_HANDLE)
IIU = SIZE(PA,1)
!
DO JI=2,IIU
  PDXM(JI,:,:)          = PA(JI,:,:) -  PA(JI-1,:,:)
END DO
!
PDXM(1,:,:)    =  PDXM(IIU-2*JPHEXT+1,:,:)
!
CALL ABORT ! AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DXM',1,ZHOOK_HANDLE)
END SUBROUTINE DXM_PHY
!     ###############################
      SUBROUTINE DXF_PHY(D,PA,PDXF)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *DXF* -  Shuman operator : finite difference operator in x direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference
!     along the x direction (I index) for a field PA localized at a x-flux
!     point (u point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------
!!        The result PDXF(i,:,:) is defined by (PA(i+1,:,:)-PA(i,:,:))
!!        At i=size(PA,1), PDXF(i,:,:) are replaced by the values of PDXF,
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
!!      Original    05/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at flux
                                                          !  side
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDXF   ! result at mass
                                                          ! localization
!
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DXF',0,ZHOOK_HANDLE)
!
CALL ABORT ! AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DXF',1,ZHOOK_HANDLE)
END SUBROUTINE DXF_PHY
!
!     ###############################
      SUBROUTINE DYF_PHY(D,PA,PDYF)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ###############################
!
!!****  *DYF* -  Shuman operator : finite difference operator in y direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference
!     along the y direction (J index) for a field PA localized at a y-flux
!     point (v point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------
!!        The result PDYF(:,j,:) is defined by (PA(:,j+1,:)-PA(:,j,:))
!!        At j=size(PA,2), PDYF(:,j,:) are replaced by the values of PDYM,
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
!!      Original    05/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at flux
                                                          !  side
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDYF   ! result at mass
                                                          ! localization
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DYF
!              ------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DYF',0,ZHOOK_HANDLE)
!
CALL ABORT ! AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DYF',1,ZHOOK_HANDLE)
END SUBROUTINE DYF_PHY
!
END MODULE SHUMAN_PHY
