MODULE MODE_GRADIENT_U_PHY
IMPLICIT NONE
CONTAINS
!     #######################################################
      SUBROUTINE GZ_U_UW_PHY(D,PA,PDZZ,PGZ_U_UW)
!     #######################################################
!
!!****  *GZ_U_UW - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Z
!!                          direction for a variable placed at the 
!!                          U point and the result is placed at
!!                          the UW vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Z cartesian direction for a field PA placed at the 
!     U point. The result is placed at the UW vorticity point.
!
!                   dzm(PA) 
!      PGZ_U_UW =   ------  
!                    ____x
!                    d*zz   
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MXM     : Shuman functions (mean operators)
!!      DZM     : Shuman functions (finite difference operators)
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
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODE_SHUMAN_PHY, ONLY : DZM_PHY, MXM_PHY
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PGZ_U_UW ! result UW point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: PA_WORK, PDZZ_WORK
!
INTEGER :: JI,JJ,JK, IIB, IIE, IJB, IJE,IKT

!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_U_UW
!              ---------------------
!
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
IKT=D%NKT
!
CALL DZM_PHY(D,PA,PA_WORK)
CALL MXM_PHY(D,PDZZ,PDZZ_WORK)
!
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
PGZ_U_UW(IIB:IIE,IJB:IJE,1:IKT)= PA_WORK(IIB:IIE,IJB:IJE,1:IKT) &
                                               / PDZZ_WORK(IIB:IIE,IJB:IJE,1:IKT)
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
!
!----------------------------------------------------------------------------
!
END SUBROUTINE GZ_U_UW_PHY
!
      SUBROUTINE GX_U_M_PHY(D,OFLAT,PA,PDXX,PDZZ,PDZX,PGX_U_M)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #######################################################
!
!!****  *GX_U_M* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian X
!!                          direction for a variable placed at the 
!!                          U point and the result is placed at
!!                          the mass point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the X cartesian direction for a field PA placed at the 
!     U point. The result is placed at the mass point.
!
!
!                       (          ______________z )
!                       (          (___________x ) )
!                    1  (          (d*zx dzm(PA) ) ) 
!      PGX_U_M =   ---- (dxf(PA) - (------------)) )
!                  ___x (          (             ) )
!                  d*xx (          (      d*zz   ) )     
!
!       
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MXF,MZF         : Shuman functions (mean operators)
!!      DXF,DZF         : Shuman functions (finite difference operators)
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
!!                  18/10/00 (V.Masson) add LFLAT switch
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODE_SHUMAN_PHY, ONLY : DZM_PHY, DXF_PHY, MXF_PHY, MZF_PHY
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
!
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: PGX_U_M ! result mass point
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT)  :: ZWORK1, ZWORK2, ZWORK3, ZWORK4
INTEGER :: IIB,IJB,IIE,IJE,IKT
INTEGER :: JI,JJ,JK
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GX_U_M
!              --------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GX_U_M',0,ZHOOK_HANDLE)
!
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
IKT=D%NKT
!
CALL DXF_PHY(D,PA,ZWORK1)
CALL MXF_PHY(D,PDXX,ZWORK2)

IF (.NOT. OFLAT) THEN
  CALL DZM_PHY(D,PA,ZWORK3)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  ZWORK3(IIB:IIE,IJB:IJE,1:IKT) = ZWORK3(IIB:IIE,IJB:IJE,1:IKT) * PDZX(IIB:IIE,IJB:IJE,1:IKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  CALL MXF_PHY(D,ZWORK3,ZWORK4)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  ZWORK4(IIB:IIE,IJB:IJE,1:IKT) = ZWORK4(IIB:IIE,IJB:IJE,1:IKT) / PDZZ(IIB:IIE,IJB:IJE,1:IKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  CALL MZF_PHY(D,ZWORK4,ZWORK3)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  PGX_U_M(IIB:IIE,IJB:IJE,1:IKT) = ( ZWORK1(IIB:IIE,IJB:IJE,1:IKT) - ZWORK3(IIB:IIE,IJB:IJE,1:IKT)) &
                                     / ZWORK2(IIB:IIE,IJB:IJE,1:IKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
ELSE
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
  PGX_U_M(IIB:IIE,IJB:IJE,1:IKT)= ZWORK1(IIB:IIE,IJB:IJE,1:IKT) / ZWORK2(IIB:IIE,IJB:IJE,1:IKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GX_U_M',1,ZHOOK_HANDLE)
END SUBROUTINE GX_U_M_PHY
!
END MODULE MODE_GRADIENT_U_PHY
