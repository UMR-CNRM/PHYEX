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
USE SHUMAN_PHY, ONLY : DZM_PHY, MXM_PHY
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
INTEGER :: JI,JJ,JK

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
CALL DZM_PHY(D,PA,PA_WORK)
CALL MXM_PHY(D,PDZZ,PDZZ_WORK)
!
!$mnh__expand_array(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
PGZ_U_UW(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT)= PA_WORK(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT) &
                                               / PDZZ_WORK(D%NIBC:D%NIEC,D%NJBC:D%NJEC,1:D%NKT)
!$mnh_end_expand_array(JI=D%NIBC:D%NIEC,JJ=D%NJBC:D%NJEC,JK=1:D%NKT)
!
!----------------------------------------------------------------------------
!
END SUBROUTINE GZ_U_UW_PHY
END MODULE MODE_GRADIENT_U_PHY
