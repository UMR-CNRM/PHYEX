MODULE MODE_GRADIENT_V_PHY
IMPLICIT NONE
CONTAINS
     !     #######################################################
      SUBROUTINE GZ_V_VW_PHY(D,PA,PDZZ,PGZ_V_VW)
!     #######################################################
!
!!****  *GZ_V_VW - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Z
!!                          direction for a variable placed at the 
!!                          V point and the result is placed at
!!                          the VW vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Z cartesian direction for a field PA placed at the 
!     V point. The result is placed at the VW vorticity point.
!
!
!                   dzm(PA) 
!      PGZ_V_VW =   ------  
!                    ____y
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
!!      MYM     : Shuman functions (mean operators)
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
USE SHUMAN_PHY, ONLY : DZM_PHY, MYM_PHY
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PGZ_V_VW ! result UW point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: PA_WORK, PDZZ_WORK
!
INTEGER :: JI,JJ,JK
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_V_VW
!              ---------------------
!
CALL DZM_PHY(D,PA,PA_WORK)
CALL MYM_PHY(D,PDZZ,PDZZ_WORK)
!
!$mnh__expand_array(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
PGZ_V_VW(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT)= PA_WORK(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) &
                                               / PDZZ_WORK(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT)
!$mnh_end_expand_array(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
!----------------------------------------------------------------------------
!
END SUBROUTINE GZ_V_VW_PHY
END MODULE MODE_GRADIENT_V_PHY
