MODULE MODE_GRADIENT_W_PHY
IMPLICIT NONE
CONTAINS
      SUBROUTINE GZ_W_M_PHY(D,PA,PDZZ,PGZ_W_M)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
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
USE SHUMAN_PHY,    ONLY: MZF_PHY, DZF_PHY
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
INTEGER :: IIB,IJB,IIE,IJE
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
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GZ_W_M',0,ZHOOK_HANDLE)
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
CALL DZF_PHY(D,PA,ZWORK1)
CALL MZF_PHY(D,PDZZ,ZWORK2)
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
PGZ_W_M(IIB:IIE,IJB:IJE,1:D%NKT)= ZWORK1(IIB:IIE,IJB:IJE,1:D%NKT)/ZWORK2(IIB:IIE,IJB:IJE,1:D%NKT)    
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GZ_W_M',1,ZHOOK_HANDLE)
END SUBROUTINE GZ_W_M_PHY
!
END MODULE MODE_GRADIENT_W_PHY
