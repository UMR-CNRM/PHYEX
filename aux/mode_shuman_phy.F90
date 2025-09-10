MODULE MODE_SHUMAN_PHY
IMPLICIT NONE
CONTAINS
!     ###############################
      SUBROUTINE MYF_PHY(D,PA,PMYF)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
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
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
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
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MYM',0,ZHOOK_HANDLE)

!POUR AROME

PMYM=PA

!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MYM',1,ZHOOK_HANDLE)
END SUBROUTINE MYM2D_PHY
!     ###############################
      SUBROUTINE MYM_PHY(D,PA,PMYM)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
INTEGER :: IJU            ! Size of the array in the y direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MYM
!              ------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
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
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
INTEGER :: JK,JIJ,IIJB,IIJE,IKT             ! Loop index
INTEGER :: IKL,IKA,IKU
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZM
!              ------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MZM',0,ZHOOK_HANDLE)
IIJB = D%NIJB
IIJE = D%NIJE
IKT=D%NKT
IKL=D%NKL
IKA=D%NKA
IKU=D%NKU
DO JK=2,IKT-1
  DO JIJ=IIJB, IIJE
    PMZM(JIJ, JK) = 0.5*( PA(JIJ, JK)+PA(JIJ, JK-IKL) )
  END DO
END DO
DO JIJ=IIJB, IIJE
  PMZM(JIJ, IKA)    = -999.
  PMZM(JIJ, IKU) = 0.5*( PA(JIJ, IKU)+PA(JIJ, IKU-IKL) )
END DO
!
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MZM',1,ZHOOK_HANDLE)
END SUBROUTINE MZM_PHY
!     ###############################
      SUBROUTINE DZM_PHY(D,PA,PDZM)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
INTEGER :: JK,JIJ,IIJB,IIJE,IKT             ! Loop index
INTEGER :: IKL, IKA, IKU
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZM
!              ------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DZM',0,ZHOOK_HANDLE)
IIJB = D%NIJB
IIJE = D%NIJE
IKT=D%NKT
IKL=D%NKL
IKA=D%NKA
IKU=D%NKU
DO JK=2,IKT-1
  DO JIJ=IIJB, IIJE
    PDZM(JIJ, JK)          = PA(JIJ, JK) -  PA(JIJ, JK-IKL)
  END DO
END DO
DO JIJ=IIJB, IIJE
  PDZM(JIJ, IKA)    =  -999.
  PDZM(JIJ, IKU)    = PA(JIJ, IKU) -  PA(JIJ, IKU-IKL)
END DO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DZM',1,ZHOOK_HANDLE)
END SUBROUTINE DZM_PHY

!     ###############################
      SUBROUTINE MXM_PHY(D,PA,PMXM)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
INTEGER :: IIU            ! Size of the array in the x direction
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MXM
!              ------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
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
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MXM',0,ZHOOK_HANDLE)

!POUR AROME

PMXM=PA

IF (LHOOK) CALL DR_HOOK('MXM',1,ZHOOK_HANDLE)
END SUBROUTINE MXM2D_PHY
!

!     ###############################
      SUBROUTINE MXF_PHY(D,PA,PMXF)      
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MXF',0,ZHOOK_HANDLE)
!POUR AROME
!
PMXF=PA
!
IF (LHOOK) CALL DR_HOOK('MXF',1,ZHOOK_HANDLE)
END SUBROUTINE MXF_PHY
!     ###############################
      SUBROUTINE MXF2D_PHY(D,PA,PMXF)      
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MXF',0,ZHOOK_HANDLE)
!POUR AROME
!
PMXF=PA
!
IF (LHOOK) CALL DR_HOOK('MXF',1,ZHOOK_HANDLE)
END SUBROUTINE MXF2D_PHY
!     ###############################
      SUBROUTINE MZF_PHY(D,PA,PMZF)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
INTEGER :: JK,JIJ,IIJB,IIJE,IKT             ! Loop index
INTEGER :: IKL, IKA, IKU
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZF
!              ------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('MZF',0,ZHOOK_HANDLE)
IIJB = D%NIJB
IIJE = D%NIJE
IKT=D%NKT
IKL=D%NKL
IKA=D%NKA
IKU=D%NKU
DO JK=2,IKT-1
  DO JIJ=IIJB, IIJE
    PMZF(JIJ, JK) = 0.5*( PA(JIJ, JK)+PA(JIJ, JK+IKL) )
  END DO
END DO
DO JIJ=IIJB, IIJE
  PMZF(JIJ, IKU) = -999.
  PMZF(JIJ, IKA) = 0.5*( PA(JIJ, IKA)+PA(JIJ, IKA+IKL) )
END DO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('MZF',1,ZHOOK_HANDLE)
END SUBROUTINE MZF_PHY
!
!     ###############################
      SUBROUTINE DXF2D_PHY(D,PA,PDXF)
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
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1    
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL

!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PDXF   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!             
INTEGER :: JJ,IJU
!             
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DXF
!              ------------------
!
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'DXF2D_PHY', 'AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE')

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DXF2D_PHY
!     ###############################
      SUBROUTINE DZF_PHY(D,PA,PDZF)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
INTEGER :: JK,JIJ,IIJB,IIJE,IKT             ! Loop index
INTEGER :: IKL, IKA, IKU
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZF
!              ------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DZF',0,ZHOOK_HANDLE)
IIJB = D%NIJB
IIJE = D%NIJE
IKT=D%NKT
IKL=D%NKL
IKA=D%NKA
IKU=D%NKU
DO JK=2,IKT-1
  DO JIJ=IIJB, IIJE
    PDZF(JIJ, JK)          = PA(JIJ, JK+IKL) -  PA(JIJ, JK)
  END DO
END DO
DO JIJ=IIJB, IIJE
  PDZF(JIJ, IKA)    = PA(JIJ, IKA+IKL) -  PA(JIJ, IKA)
  PDZF(JIJ, IKU)    = -999.
END DO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DZF',1,ZHOOK_HANDLE)
END SUBROUTINE DZF_PHY
      SUBROUTINE DXM2D_PHY(D,PA,PDXM)
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
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PDXM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!             
INTEGER :: JJ,IJU
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DXM
!              ------------------
!
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'DXM2D_PHY', 'AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE')
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DXM2D_PHY
!
!     ###############################
      SUBROUTINE DYM_PHY(D,PA,PDYM)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DYM',0,ZHOOK_HANDLE)
IJU=SIZE(PA,2)
!
DO JJ=2,IJU
  PDYM(:,JJ,:)          = PA(:,JJ,:) -  PA(:,JJ-1,:)
END DO
!
PDYM(:,1,:)    =  PDYM(:,IJU-2*JPHEXT+1,:)
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'DYM_PHY', 'AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE')
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DYM',1,ZHOOK_HANDLE)
END SUBROUTINE DYM_PHY
!
!     ###############################
      SUBROUTINE DYM2D_PHY(D,PA,PDYM)
!     ###############################
!
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PDYM   ! result at flux

CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'DYM2D_PHY', 'AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE')

END SUBROUTINE DYM2D_PHY

!     ###############################
      SUBROUTINE DXM_PHY(D,PA,PDXM)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DXM',0,ZHOOK_HANDLE)
IIU = SIZE(PA,1)
!
DO JI=2,IIU
  PDXM(JI,:,:)          = PA(JI,:,:) -  PA(JI-1,:,:)
END DO
!
PDXM(1,:,:)    =  PDXM(IIU-2*JPHEXT+1,:,:)
!
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'DXM_PHY', 'AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE')
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DXM',1,ZHOOK_HANDLE)
END SUBROUTINE DXM_PHY
!     ###############################
      SUBROUTINE DXF_PHY(D,PA,PDXF)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DXF',0,ZHOOK_HANDLE)
!
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'DXF_PHY', 'AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE')
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DXF',1,ZHOOK_HANDLE)
END SUBROUTINE DXF_PHY
!
!     ###############################
      SUBROUTINE DYF_PHY(D,PA,PDYF)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('DYF',0,ZHOOK_HANDLE)
!
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'DYF_PHY', 'AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE')
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('DYF',1,ZHOOK_HANDLE)
END SUBROUTINE DYF_PHY
!
!     ###############################
      SUBROUTINE DYF2D_PHY(D,PA,PDYF)
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
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODE_MSG, ONLY: PRINT_MSG, NVERB_FATAL
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PDYF   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ ,JI           ! Loop index in y direction
INTEGER :: IJU           ! upper bound in y direction of PA 
!
!          
INTEGER :: IIU
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DYF
!              ------------------
!
CALL PRINT_MSG(NVERB_FATAL, 'GEN', 'DYF2D_PHY', 'AROME SHOULD NOT CALLED HORIZONTAL FINITE DIFFERENCE')
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DYF2D_PHY
!
END MODULE MODE_SHUMAN_PHY
