MODULE MODE_GRADIENT_W_PHY
IMPLICIT NONE
CONTAINS
      SUBROUTINE GX_W_UW_PHY(D,OFLAT,PA,PDXX,PDZZ,PDZX,PGX_W_UW)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #########################################################
!
!!****  *GX_W_UW* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian X
!!                          direction for a variable placed at the 
!!                          V point and the result is placed at
!!                          the UW vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the X cartesian direction for a field PA placed at the 
!     W point. The result is placed at the UW vorticity point.
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MXM,MZM,MZF     : Shuman functions (mean operators)
!!      DXM,DZM         : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (GRAD_CAR operators)
!!      A Turbulence scheme for the Meso-NH model (Chapter 6)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart        *INM and Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/07/94
!!                  18/10/00 (V.Masson) add LFLAT switch
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODE_SHUMAN_PHY,    ONLY: MZF_PHY, MXM_PHY, DXM_PHY, MZM_PHY, DZM_PHY
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
LOGICAL, INTENT(IN) :: OFLAT
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),INTENT(OUT) :: PGX_W_UW ! result UW point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5
INTEGER :: IIB,IJB,IIE,IJE,IKT
INTEGER :: JI,JJ,JK
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GX_W_UW
!              ---------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GX_W_UW',0,ZHOOK_HANDLE)
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
IKT=D%NKT
CALL MZM_PHY(D,PDXX,ZWORK1)
CALL DXM_PHY(D,PA,ZWORK2)
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
ZWORK3(IIB:IIE,IJB:IJE,:) = ZWORK2(IIB:IIE,IJB:IJE,:) / ZWORK1(IIB:IIE,IJB:IJE,:)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
!
IF (.NOT. OFLAT) THEN
  CALL MZF_PHY(D,PA,ZWORK2)
  CALL MXM_PHY(D,ZWORK2,ZWORK4)
  CALL DZM_PHY(D,ZWORK4,ZWORK5)
  !
  CALL MXM_PHY(D,PDZZ,ZWORK2)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  PGX_W_UW(IIB:IIE,IJB:IJE,:)= ZWORK3(IIB:IIE,IJB:IJE,:)    &
                 - ZWORK5(IIB:IIE,IJB:IJE,:)*PDZX(IIB:IIE,IJB:IJE,:)  &
                 / (ZWORK1(IIB:IIE,IJB:IJE,:)*ZWORK2(IIB:IIE,IJB:IJE,:))
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
ELSE
  PGX_W_UW = ZWORK3
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GX_W_UW',1,ZHOOK_HANDLE)
END SUBROUTINE GX_W_UW_PHY
!
      SUBROUTINE GY_W_VW_PHY(D,OFLAT,PA,PDYY,PDZZ,PDZY,PGY_W_VW)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #########################################################
!
!!****  *GY_W_VW* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Y
!!                          direction for a variable placed at the 
!!                          W point and the result is placed at
!!                          the VW vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Y cartesian direction for a field PA placed at the 
!     W point. The result is placed at the VW vorticity point.
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MYM,MZM,MZF     : Shuman functions (mean operators)
!!      DYM,DZM         : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (GRAD_CAR operators)
!!      A Turbulence scheme for the Meso-NH model (Chapter 6)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart        *INM and Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    20/07/94
!!                  18/10/00 (V.Masson) add LFLAT switch
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODE_SHUMAN_PHY,    ONLY: MZF_PHY, MYM_PHY, DYM_PHY, MZM_PHY, DZM_PHY
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDYY    ! metric coefficient dxx
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZY    ! metric coefficient dzx
LOGICAL, INTENT(IN) :: OFLAT
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),INTENT(OUT) :: PGY_W_VW ! result UW point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5
INTEGER :: IIB,IJB,IIE,IJE,IKT
INTEGER :: JI,JJ,JK
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_W_VW
!              ---------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GY_W_VW',0,ZHOOK_HANDLE)
!IF (.NOT. LFLAT) THEN
!  PGY_W_VW(:,:,:)= DYM(PA(:,:,:))/(MZM(PDYY(:,:,:), KKA, KKU, KL))    &
!                 -DZM(MYM(MZF(PA(:,:,:), KKA, KKU, KL)), KKA, KKU, KL)*PDZY(:,:,:)  &
!                  /( MZM(PDYY(:,:,:), KKA, KKU, KL)*MYM(PDZZ(:,:,:)) )
!ELSE
!  PGY_W_VW(:,:,:)= DYM(PA(:,:,:))/(MZM(PDYY(:,:,:), KKA, KKU, KL))
!END IF
!
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
IKT=D%NKT
CALL MZM_PHY(D,PDYY,ZWORK1)
CALL DYM_PHY(D,PA,ZWORK2)
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
ZWORK3(IIB:IIE,IJB:IJE,:) = ZWORK2(IIB:IIE,IJB:IJE,:) / ZWORK1(IIB:IIE,IJB:IJE,:)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
!
IF (.NOT. OFLAT) THEN
  CALL MZF_PHY(D,PA,ZWORK2)
  CALL MYM_PHY(D,ZWORK2,ZWORK4)
  CALL DZM_PHY(D,ZWORK4,ZWORK5)
  !
  CALL MYM_PHY(D,PDZZ,ZWORK2)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  PGY_W_VW(IIB:IIE,IJB:IJE,:)= ZWORK3(IIB:IIE,IJB:IJE,:)    &
                 - ZWORK5(IIB:IIE,IJB:IJE,:)*PDZY(IIB:IIE,IJB:IJE,:)  &
                 / (ZWORK1(IIB:IIE,IJB:IJE,:)*ZWORK2(IIB:IIE,IJB:IJE,:))
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
ELSE
  PGY_W_VW = ZWORK3
END IF

!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GY_W_VW',1,ZHOOK_HANDLE)
END SUBROUTINE GY_W_VW_PHY
!
      SUBROUTINE GZ_W_M_PHY(D,PA,PDZZ,PGZ_W_M)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #######################################################
!
!!****  *GZ_W_M* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Z
!!                          direction for a variable placed at the 
!!                          W point and the result is placed at
!!                          the mass point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Z cartesian direction for a field PA placed at the 
!     W point. The result is placed at the mass point.
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MZF     : Shuman functions (mean operators)
!!      DZF     : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (GRAD_CAR operators)
!!      A Turbulence scheme for the Meso-NH model (Chapter 6)
!!
!!    AUTHOR
!!    ------
!!      Joan Cuxart        *INM and Meteo-France*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    19/07/94
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODE_SHUMAN_PHY,    ONLY: MZF_PHY, DZF_PHY
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) , INTENT(OUT):: PGZ_W_M ! result mass point
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZWORK1, ZWORK2
INTEGER :: IIB,IJB,IIE,IJE,IKT
INTEGER :: JI,JJ,JK
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_W_M
!              --------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GZ_W_M',0,ZHOOK_HANDLE)
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
IKT=D%NKT
CALL DZF_PHY(D,PA,ZWORK1)
CALL MZF_PHY(D,PDZZ,ZWORK2)
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
PGZ_W_M(IIB:IIE,IJB:IJE,:)= ZWORK1(IIB:IIE,IJB:IJE,:)/ZWORK2(IIB:IIE,IJB:IJE,:)    
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GZ_W_M',1,ZHOOK_HANDLE)
END SUBROUTINE GZ_W_M_PHY
!
END MODULE MODE_GRADIENT_W_PHY
