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
!
SUBROUTINE GX_M_M_PHY(D,OFLAT,PA,PDXX,PDZZ,PDZX,PGX_M_M)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #######################################################
!
!!****  *GX_M_M* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian X
!!                          direction for a variable placed at the 
!!                          mass point and the result is placed at
!!                          the mass point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the X cartesian direction for a field PA placed at the 
!     mass point. The result is placed at the mass point.
!
!
!                       (          ______________z )
!                       (          (___x         ) )
!                    1  (    _x    (d*zx dzm(PA) ) ) 
!      PGX_M_M =   ---- (dxf(PA) - (------------)) )
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
!!      MXM,MXF,MZF     : Shuman functions (mean operators)
!!      DXF,DZF         : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_CONF : LFLAT
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
!!      Original    18/07/94
!!                  19/07/00  add the LFLAT switch (J. Stein)
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE SHUMAN_PHY, ONLY: DXF_PHY, MZF_PHY, DZM_PHY, MXF_PHY, MXM_PHY
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),       INTENT(IN)   :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
LOGICAL, INTENT(IN) :: OFLAT
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PGX_M_M ! result mass point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5, ZWORK6, ZMXF_PDXX
!
INTEGER :: IIB,IJB,IIE,IJE
INTEGER :: JI,JJ,JK
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GX_M_M
!              --------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GX_M_M',0,ZHOOK_HANDLE)
!
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
CALL MXF_PHY(D,PDXX,ZMXF_PDXX)
CALL MXM_PHY(D,PA,ZWORK1)
CALL DXF_PHY(D,ZWORK1,ZWORK2)
!
IF (.NOT. OFLAT) THEN
  CALL DZM_PHY(D,PA,ZWORK3)
  CALL MXF_PHY(D,PDZX,ZWORK4)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
  ZWORK5(IIB:IIE,IJB:IJE,1:D%NKT) = ZWORK3(IIB:IIE,IJB:IJE,1:D%NKT) * ZWORK4(IIB:IIE,IJB:IJE,1:D%NKT) &
                                    / PDZZ(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
  CALL MZF_PHY(D,ZWORK5,ZWORK6)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
  PGX_M_M(IIB:IIE,IJB:IJE,1:D%NKT)= (ZWORK2(IIB:IIE,IJB:IJE,1:D%NKT) - ZWORK6(IIB:IIE,IJB:IJE,1:D%NKT)) &
                                    / ZMXF_PDXX(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
ELSE
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
  PGX_M_M(IIB:IIE,IJB:IJE,1:D%NKT)= ZWORK2(IIB:IIE,IJB:IJE,1:D%NKT) / ZMXF_PDXX(IIB:IIE,IJB:IJE,1:D%NKT) 
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GX_M_M',1,ZHOOK_HANDLE)
END SUBROUTINE GX_M_M_PHY
!
      SUBROUTINE GY_M_M_PHY(D,OFLAT,PA,PDYY,PDZZ,PDZY,PGY_M_M)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     #######################################################
!
!!****  *GY_M_M* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Y
!!                          direction for a variable placed at the 
!!                          mass point and the result is placed at
!!                          the mass point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Y cartesian direction for a field PA placed at the 
!     mass point. The result is placed at the mass point.
!
!
!                       (          ______________z )
!                       (          (___y         ) )
!                    1  (    _y    (d*zy dzm(PA) ) ) 
!      PGY_M_M =   ---- (dyf(PA) - (------------)) )
!                  ___y (          (             ) )
!                  d*yy (          (      d*zz   ) )     
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
!!      MYM,MYF,MZF     : Shuman functions (mean operators)
!!      DYF,DZF         : Shuman functions (finite difference operators)
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_CONF : LFLAT
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
!!      Original    18/07/94
!!                  19/07/00  add the LFLAT switch (J. Stein)
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE SHUMAN_PHY, ONLY: DYF_PHY, MZF_PHY, DZM_PHY, MYF_PHY, MYM_PHY
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),        INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
LOGICAL, INTENT(IN) :: OFLAT
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),INTENT(OUT) :: PGY_M_M ! result mass point
REAL, DIMENSION(D%NIT,D%NJT,D%NKT) :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5, ZMYF_PDYY
!
INTEGER :: IIB,IJB,IIE,IJE
INTEGER :: JI,JJ,JK
!
!*       0.2   declaration of local variables
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GY_M_M',0,ZHOOK_HANDLE)
!
IIE=D%NIEC
IIB=D%NIBC
IJE=D%NJEC
IJB=D%NJBC
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_M_M
!              --------------------
!
CALL MYM_PHY(D,PA,ZWORK1)
CALL DYF_PHY(D,ZWORK1,ZWORK2)
CALL MYF_PHY(D,PDYY,ZMYF_PDYY)
!
IF (.NOT. OFLAT) THEN
  !
  CALL DZM_PHY(D,PA,ZWORK3)
  CALL MYF_PHY(D,PDZY,ZWORK4)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
  ZWORK5(IIB:IIE,IJB:IJE,1:D%NKT) = ZWORK4(IIB:IIE,IJB:IJE,1:D%NKT) * ZWORK3(IIB:IIE,IJB:IJE,1:D%NKT) &
                                   / PDZZ(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
  CALL MZF_PHY(D,ZWORK5,ZWORK4)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
  PGY_M_M(IIB:IIE,IJB:IJE,1:D%NKT)= (ZWORK2(IIB:IIE,IJB:IJE,1:D%NKT)-ZWORK4(IIB:IIE,IJB:IJE,1:D%NKT)) &
                                    /ZMYF_PDYY(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
ELSE
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
  PGY_M_M(IIB:IIE,IJB:IJE,1:D%NKT) = ZWORK2(IIB:IIE,IJB:IJE,1:D%NKT)/ZMYF_PDYY(IIB:IIE,IJB:IJE,1:D%NKT)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:D%NKT)    
ENDIF  
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GY_M_M',1,ZHOOK_HANDLE)
END SUBROUTINE GY_M_M_PHY
END MODULE MODE_GRADIENT_M_PHY
