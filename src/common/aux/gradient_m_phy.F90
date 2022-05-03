MODULE MODE_GRADIENT_M_PHY
IMPLICIT NONE
CONTAINS
!     #########################################
      SUBROUTINE GZ_M_W_PHY(D,PY,PDZZ,PGZ_M_W)
!     #########################################
!
!!****  *GZ_M_W * - Compute the gradient along z direction for a 
!!       variable localized at a mass point
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute a gradient along x,y,z 
!     directions for a field PY localized at a mass point. The result PGZ_M_W
!     is localized at a z-flux point (w point)
!
!              
!                    dzm(PY)  
!       PGZ_M_W =    ------- 
!                     d*zz        
!
!!**  METHOD
!!    ------
!!      We employ the Shuman operators to compute the derivatives and the 
!!    averages. The metric coefficients PDZZ are dummy arguments.
!!
!!
!!    EXTERNAL
!!    --------
!!      FUNCTION DZM : compute a finite difference along the z 
!!    direction for a variable at a mass localization
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------  
!!      Module MODI_SHUMAN : interface for the Shuman functions
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (function GZ_M_W)
!!      
!!
!!    AUTHOR
!!    ------
!!	P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification       16/03/95  change the order of the arguments
!!                         19/07/00  inlining(J. Stein)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!              -------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PGZ_M_W  ! result at flux
                                                              ! side
!
INTEGER :: IKT,IKTB,IKTE,IIB,IJB,IIE,IJE
INTEGER :: JI,JJ,JK
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Z
!              -----------------------------
!
IKT=D%NKT
IKTB=D%NKTB              
IKTE=D%NKTE
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
DO JK=IKTB,IKTE 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PGZ_M_W(JI,JJ,JK) =  (PY(JI,JJ,JK)-PY(JI,JJ,JK-D%NKL )) / PDZZ(JI,JJ,JK)
  ENDDO
 ENDDO
ENDDO
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
PGZ_M_W(IIB:IIE,IJB:IJE,D%NKU)=  (PY(IIB:IIE,IJB:IJE,D%NKU)-PY(IIB:IIE,IJB:IJE,D%NKU-D%NKL))  &
                           / PDZZ(IIB:IIE,IJB:IJE,D%NKU)
PGZ_M_W(IIB:IIE,IJB:IJE,D%NKA)= PGZ_M_W(IIB:IIE,IJB:IJE,D%NKU) ! -999.
!$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GZ_M_W_PHY
END MODULE MODE_GRADIENT_M_PHY
