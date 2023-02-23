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
INTEGER :: IKT,IKTB,IKTE,IIB,IJB,IIE,IJE,IKA,IKU,IKL
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
IKT=D%NKT
IKA=D%NKA
IKU=D%NKU
IKL=D%NKL
!
DO JK=IKTB,IKTE 
  DO JJ=IJB,IJE 
    DO JI=IIB,IIE 
      PGZ_M_W(JI,JJ,JK) =  (PY(JI,JJ,JK)-PY(JI,JJ,JK-IKL )) / PDZZ(JI,JJ,JK)
  ENDDO
 ENDDO
ENDDO
!$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE)
PGZ_M_W(IIB:IIE,IJB:IJE,IKU)=  (PY(IIB:IIE,IJB:IJE,IKU)-PY(IIB:IIE,IJB:IJE,IKU-IKL))  &
                           / PDZZ(IIB:IIE,IJB:IJE,IKU)
PGZ_M_W(IIB:IIE,IJB:IJE,IKA)= PGZ_M_W(IIB:IIE,IJB:IJE,IKU) ! -999.
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
USE MODE_SHUMAN_PHY, ONLY: DXF_PHY, MZF_PHY, DZM_PHY, MXF_PHY, MXM_PHY
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
INTEGER :: IIB,IJB,IIE,IJE,IKT
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
IKT=D%NKT
!
CALL MXF_PHY(D,PDXX,ZMXF_PDXX)
CALL MXM_PHY(D,PA,ZWORK1)
CALL DXF_PHY(D,ZWORK1,ZWORK2)
!
IF (.NOT. OFLAT) THEN
  CALL DZM_PHY(D,PA,ZWORK3)
  CALL MXF_PHY(D,PDZX,ZWORK4)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
  ZWORK5(IIB:IIE,IJB:IJE,:) = ZWORK3(IIB:IIE,IJB:IJE,:) * ZWORK4(IIB:IIE,IJB:IJE,:) &
                                    / PDZZ(IIB:IIE,IJB:IJE,:)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
  CALL MZF_PHY(D,ZWORK5,ZWORK6)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
  PGX_M_M(IIB:IIE,IJB:IJE,:)= (ZWORK2(IIB:IIE,IJB:IJE,:) - ZWORK6(IIB:IIE,IJB:IJE,:)) &
                                    / ZMXF_PDXX(IIB:IIE,IJB:IJE,:)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
ELSE
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
  PGX_M_M(IIB:IIE,IJB:IJE,:)= ZWORK2(IIB:IIE,IJB:IJE,:) / ZMXF_PDXX(IIB:IIE,IJB:IJE,:) 
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
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
USE MODE_SHUMAN_PHY, ONLY: DYF_PHY, MZF_PHY, DZM_PHY, MYF_PHY, MYM_PHY
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
INTEGER :: IIB,IJB,IIE,IJE,IKT
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
IKT=D%NKT
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
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
  ZWORK5(IIB:IIE,IJB:IJE,:) = ZWORK4(IIB:IIE,IJB:IJE,:) * ZWORK3(IIB:IIE,IJB:IJE,:) &
                                   / PDZZ(IIB:IIE,IJB:IJE,:)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
  CALL MZF_PHY(D,ZWORK5,ZWORK4)
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
  PGY_M_M(IIB:IIE,IJB:IJE,:)= (ZWORK2(IIB:IIE,IJB:IJE,:)-ZWORK4(IIB:IIE,IJB:IJE,:)) &
                                    /ZMYF_PDYY(IIB:IIE,IJB:IJE,:)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
ELSE
  !$mnh_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
  PGY_M_M(IIB:IIE,IJB:IJE,:) = ZWORK2(IIB:IIE,IJB:IJE,:)/ZMYF_PDYY(IIB:IIE,IJB:IJE,:)
  !$mnh_end_expand_array(JI=IIB:IIE,JJ=IJB:IJE,JK=1:IKT)    
ENDIF  
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GY_M_M',1,ZHOOK_HANDLE)
END SUBROUTINE GY_M_M_PHY
!
!     #######################################################
      SUBROUTINE GX_M_U_PHY(D,OFLAT,PY,PDXX,PDZZ,PDZX,PGX_M_U)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################################################
!
!!****  *GX_M_U * - Compute the gradient along x for a variable localized at
!!                  a mass point
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute a gradient along x
!     direction for a field PY localized at a mass point. The result PGX_M_U
!     is localized at a x-flux point (u point).
!
!                    (           ____________z )
!                    (               ________x )
!                 1  (                dzm(PY)  )
!   PGX_M_U =   ---- (dxm(PY) - d*zx --------  )
!               d*xx (                 d*zz    )
!
!
!
!!**  METHOD
!!    ------
!!      We employ the Shuman operators to compute the derivatives and the
!!    averages. The metric coefficients PDXX,PDZX,PDZZ are dummy arguments.
!!
!!
!!    EXTERNAL
!!    --------
!!      FUNCTION DXM: compute a finite difference along the x direction for
!!      a variable at a mass localization
!!      FUNCTION DZM: compute a finite difference along the y direction for
!!      a variable at a mass localization
!!      FUNCTION MXM: compute an average in the x direction for a variable
!!      at a mass localization
!!      FUNCTION MZF: compute an average in the z direction for a variable
!!      at a flux side
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_CONF : LFLAT
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (function GX_M_U)
!!
!!
!!    AUTHOR
!!    ------
!!      P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original           05/07/94
!!      Modification       16/03/95  change the order of the arguments
!!                         19/07/00  add the LFLAT switch  + inlining(J. Stein)
!!                         20/08/00  optimization (J. Escobar)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT_TURB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),        INTENT(IN)  :: D
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(D%NIT*D%NJT*D%NKT),  INTENT(IN)  :: PY      ! variable at the mass point
REAL, DIMENSION(D%NIT*D%NJT*D%NKT),  INTENT(IN)  :: PDXX    ! metric coefficient dyy
REAL, DIMENSION(D%NIT*D%NJT*D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIT*D%NJT*D%NKT),  INTENT(IN)  :: PDZX    ! metric coefficient dzy
!
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PGX_M_U  ! result at flux
                                                              ! side
REAL, DIMENSION(D%NIT*D%NJT*D%NKT) :: ZGX_M_U
REAL, DIMENSION(D%NIT,D%NJT,D%NKT):: ZY, ZDXX
INTEGER  IIU,IKU,JI,JK,IKL, IKA
!
INTEGER :: JJK,IJU
INTEGER :: JIJK,JIJKOR,JIJKEND
INTEGER :: JI_1JK, JIJK_1, JI_1JK_1, JIJKP1, JI_1JKP1
!
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG X
!              -----------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GX_M_U',0,ZHOOK_HANDLE)
IIU=D%NIT
IJU=D%NJT
IKU=D%NKT
IKL=D%NKL
IKA=D%NKA
!
CALL D1D_TO_3D(D,PDXX,ZDXX)
CALL D1D_TO_3D(D,PY,ZY)
!
IF (.NOT. OFLAT) THEN
  JIJKOR  = 1 + JPHEXT + IIU*IJU*(JPVEXT_TURB+1 - 1)
  JIJKEND = IIU*IJU*(IKU-JPVEXT_TURB)
!CDIR NODEP
!OCL NOVREC
  DO JIJK=JIJKOR , JIJKEND
! indexation
    JI_1JK   = JIJK - 1
    JIJK_1   = JIJK     - IIU*IJU*IKL
    JI_1JK_1 = JIJK - 1 - IIU*IJU*IKL
    JIJKP1   = JIJK     + IIU*IJU*IKL
    JI_1JKP1 = JIJK - 1 + IIU*IJU*IKL
!
    ZGX_M_U(JIJK)=                                              &
       (  PY(JIJK)-PY(JI_1JK)                               &
       -(  (PY(JIJK)-PY(JIJK_1))     / PDZZ(JIJK)       &
       +(PY(JI_1JK)-PY(JI_1JK_1)) / PDZZ(JI_1JK)        &
       ) * PDZX(JIJK)* 0.25                                     &
       -(  (PY(JIJKP1)-PY(JIJK))     / PDZZ(JIJKP1)     &
       +(PY(JI_1JKP1)-PY(JI_1JK)) / PDZZ(JI_1JKP1)      &
       ) * PDZX(JIJKP1)* 0.25                                   &
        )  / PDXX(JIJK)
  END DO
!
CALL D1D_TO_3D(D,ZGX_M_U,PGX_M_U)
!
  DO JI=1+JPHEXT,IIU
    PGX_M_U(JI,:,IKU)=  ( ZY(JI,:,IKU)-ZY(JI-1,:,IKU)  )  / ZDXX(JI,:,IKU)
    PGX_M_U(JI,:,IKA)=  -999.
  END DO
ELSE
!  PGX_M_U = DXM(PY) / PDXX
  PGX_M_U(1+1:IIU,:,:) = ( ZY(1+1:IIU,:,:)-ZY(1:IIU-1,:,:) ) &
                             / ZDXX(1+1:IIU,:,:)
!
ENDIF
DO JI=1,JPHEXT
  PGX_M_U(JI,:,:)=PGX_M_U(IIU-2*JPHEXT+JI,:,:) ! for reprod JPHEXT <> 1
END DO  
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GX_M_U',1,ZHOOK_HANDLE)
END SUBROUTINE GX_M_U_PHY
!
      SUBROUTINE GY_M_V_PHY(D,OFLAT,PY,PDYY,PDZZ,PDZY,PGY_M_V)
      USE PARKIND1, ONLY : JPRB
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK
!     ##################################################
!
!!****  *GY_M_V * - Compute the gradient along y for a variable localized at
!!                  a mass point
!!
!!    PURPOSE
!!    -------
!       The purpose of this routine is to compute a gradient along y
!     direction for a field PY localized at a mass point. The result PGY_M_V
!     is localized at a y-flux point (v point).
!
!                    (           ____________z )
!                    (               ________y )
!                 1  (                dzm(PY)  )
!   PGY_M_V =   ---- (dym(PY) - d*zy --------  )
!               d*yy (                 d*zz    )
!
!
!
!
!!**  METHOD
!!    ------
!!      We employ the Shuman operators to compute the derivatives and the
!!    averages. The metric coefficients PDYY,PDZY,PDZZ are dummy arguments.
!!
!!
!!    EXTERNAL
!!    --------
!!      FUNCTION DYM: compute a finite difference along the y direction for
!!      a variable at a mass localization
!!      FUNCTION DZM: compute a finite difference along the y direction for
!!      a variable at a mass localization
!!      FUNCTION MYM: compute an average in the x direction for a variable
!!      at a mass localization
!!      FUNCTION MZF: compute an average in the z direction for a variable
!!      at a flux side
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      MODD_CONF : LFLAT
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (function GY_M_V)
!!
!!
!!    AUTHOR
!!    ------
!!      P. Hereil and J. Stein       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94
!!      Modification       16/03/95  change the order of the arguments
!!                         19/07/00  add the LFLAT switch + inlining(J. Stein)
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY : JPHEXT, JPVEXT_TURB
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!              -------------------------------------
!
TYPE(DIMPHYEX_t),        INTENT(IN)  :: D
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PDYY                   !d*yy
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PDZY                   !d*zy
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PDZZ                   !d*zz
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization

REAL, DIMENSION(D%NIT,D%NJT,D%NKT),INTENT(OUT) :: PGY_M_V  ! result at flux
                                                              ! side
!REAL, DIMENSION(D%NIT*D%NJT*D%NKT) :: ZGY_M_V
!REAL, DIMENSION(D%NIT,D%NJT,D%NKT):: ZY, ZDYY,ZDZZ,ZDZY
INTEGER  IJU,IKU,JI,JJ,JK,IKL, IKA
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Y
!              ----------------------------
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GY_M_V',0,ZHOOK_HANDLE)
IJU=D%NJT
IKU=D%NKT
IKL=D%NKL
IKA=D%NKA
IKU=D%NKU
IF (.NOT. OFLAT) THEN
!  PGY_M_V = (   DYM(PY)  -  MZF (   MYM(  DZM(PY) /PDZZ  ) * PDZY   )   )/PDYY
  DO JK=1+JPVEXT_TURB,IKU-JPVEXT_TURB
    DO JJ=1+JPHEXT,IJU
        PGY_M_V(:,JJ,JK)=                                                 &
           (  PY(:,JJ,JK)-PY(:,JJ-1,JK)                                   &
             -(  (PY(:,JJ,JK)-PY(:,JJ,JK-IKL))     / PDZZ(:,JJ,JK)          &
                +(PY(:,JJ-1,JK)-PY(:,JJ-IKL,JK-IKL)) / PDZZ(:,JJ-1,JK)        &
              ) * PDZY(:,JJ,JK)* 0.25                                     &
             -(  (PY(:,JJ,JK+IKL)-PY(:,JJ,JK))     / PDZZ(:,JJ,JK+IKL)        &
                +(PY(:,JJ-1,JK+IKL)-PY(:,JJ-1,JK)) / PDZZ(:,JJ-1,JK+IKL)      &
              ) * PDZY(:,JJ,JK+IKL)* 0.25                                   &
            )  / PDYY(:,JJ,JK)
    END DO
  END DO
!
  DO JJ=1+JPHEXT,IJU
    PGY_M_V(:,JJ,IKU)=  ( PY(:,JJ,IKU)-PY(:,JJ-1,IKU)  )  / PDYY(:,JJ,IKU)
    PGY_M_V(:,JJ,IKA)=  -999.
  END DO
!
ELSE
!  PGY_M_V = DYM(PY)/PDYY
  PGY_M_V(:,1+1:IJU,:) = ( PY(:,1+1:IJU,:)-PY(:,1:IJU-1,:) ) &
                               / PDYY(:,1+1:IJU,:)
!
ENDIF
DO JJ=1,JPHEXT
  PGY_M_V(:,JJ,:)=PGY_M_V(:,IJU-2*JPHEXT+JJ,:)
END DO
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GY_M_V',1,ZHOOK_HANDLE)
END SUBROUTINE GY_M_V_PHY
!
SUBROUTINE D1D_TO_3D (D,P1D,P3D)
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
TYPE(DIMPHYEX_t),        INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)  :: P1D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(OUT)  :: P3D

P3D = P1D
END SUBROUTINE D1D_TO_3D
END MODULE MODE_GRADIENT_M_PHY
