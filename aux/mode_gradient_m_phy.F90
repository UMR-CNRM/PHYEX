MODULE MODE_GRADIENT_M_PHY

!$ACDC singlecolumn 

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
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation (function GZ_M_W)
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
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PGZ_M_W  ! result at flux
                                                              ! side
!
INTEGER :: IKT,IKTB,IKTE,IIJB,IIJE,IKA,IKU,IKL
INTEGER :: JIJ,JK
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Z
!              -----------------------------
!
IKT=D%NKT
IKTB=D%NKTB              
IKTE=D%NKTE
IIJE=D%NIJE
IIJB=D%NIJB
IKT=D%NKT
IKA=D%NKA
IKU=D%NKU
IKL=D%NKL
!

DO JK=IKTB, IKTE
  DO JIJ=IIJB, IIJE
    PGZ_M_W(JIJ, JK) =  (PY(JIJ, JK)-PY(JIJ, JK - IKL)) &
                               / PDZZ(JIJ, JK)
  END DO
END DO
DO JIJ=IIJB, IIJE
  PGZ_M_W(JIJ, IKU)=  (PY(JIJ, IKU)-PY(JIJ, IKU-IKL))  &
                             / PDZZ(JIJ, IKU)
  PGZ_M_W(JIJ, IKA)= PGZ_M_W(JIJ, IKU) ! -999.
END DO

!
!-------------------------------------------------------------------------------
!
END SUBROUTINE GZ_M_W_PHY
!
SUBROUTINE GX_M_M_PHY(D,OFLAT,PA,PDXX,PDZZ,PDZX,PGX_M_M)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
LOGICAL, INTENT(IN) :: OFLAT
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PGX_M_M ! result mass point
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5, ZWORK6, ZMXF_PDXX
!
INTEGER :: IIJB,IIJE,IKT
INTEGER :: JIJ,JK
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
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GX_M_M',0,ZHOOK_HANDLE)
!
IIJB=D%NIJB
IIJE=D%NIJE
IKT=D%NKT
!
CALL MXF_PHY(D,PDXX,ZMXF_PDXX)
CALL MXM_PHY(D,PA,ZWORK1)
CALL DXF_PHY(D,ZWORK1,ZWORK2)
!
IF (.NOT. OFLAT) THEN
  CALL DZM_PHY(D,PA,ZWORK3)
  CALL MXF_PHY(D,PDZX,ZWORK4)
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK5(JIJ, JK) = ZWORK3(JIJ, JK) * ZWORK4(JIJ, JK) &
                                        / PDZZ(JIJ, JK)
    END DO
  END DO
  
  CALL MZF_PHY(D,ZWORK5,ZWORK6)
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGX_M_M(JIJ, JK)= (ZWORK2(JIJ, JK) - ZWORK6(JIJ, JK)) &
                                        / ZMXF_PDXX(JIJ, JK)
    END DO
  END DO
  
ELSE
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGX_M_M(JIJ, JK)= ZWORK2(JIJ, JK) / ZMXF_PDXX(JIJ, JK) 
    END DO
  END DO
  
END IF
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GX_M_M',1,ZHOOK_HANDLE)
END SUBROUTINE GX_M_M_PHY
!
      SUBROUTINE GY_M_M_PHY(D,OFLAT,PA,PDYY,PDZZ,PDZY,PGY_M_M)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
LOGICAL, INTENT(IN) :: OFLAT
!
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PGY_M_M ! result mass point
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1, ZWORK2, ZWORK3, ZWORK4, ZWORK5, ZMYF_PDYY
!
INTEGER :: IIJB,IIJE,IKT
INTEGER :: JIJ,JK
!
!*       0.2   declaration of local variables
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GY_M_M',0,ZHOOK_HANDLE)
!
IIJE=D%NIJE
IIJB=D%NIJB
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
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK5(JIJ, JK) = ZWORK4(JIJ, JK) * ZWORK3(JIJ, JK) &
                                       / PDZZ(JIJ, JK)
    END DO
  END DO
  
  CALL MZF_PHY(D,ZWORK5,ZWORK4)
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGY_M_M(JIJ, JK)= (ZWORK2(JIJ, JK)-ZWORK4(JIJ, JK)) &
                                        /ZMYF_PDYY(JIJ, JK)
    END DO
  END DO
  
ELSE
  
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGY_M_M(JIJ, JK) = ZWORK2(JIJ, JK)/ZMYF_PDYY(JIJ, JK)
    END DO
  END DO
  
ENDIF  
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GY_M_M',1,ZHOOK_HANDLE)
END SUBROUTINE GY_M_M_PHY
!
!     #######################################################
      SUBROUTINE GZ_M_M_PHY(D, PA,PDZZ, PGZ_M_M)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
!     #######################################################
!
!!****  *GZ_M_M* - Cartesian Gradient operator:
!!                          computes the gradient in the cartesian Z
!!                          direction for a variable placed at the
!!                          mass point and the result is placed at
!!                          the mass point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient
!     along the Z cartesian direction for a field PA placed at the
!     mass point. The result is placed at the mass point.
!
!                 _________z
!                 (dzm(PA))
!      PGZ_M_M =  (------ )
!                 ( d*zz  )
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
!!      MZF     : Shuman functions (mean operators)
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
!!      Original    18/07/94
!-------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODE_SHUMAN_PHY, ONLY: MZF_PHY, DZM_PHY
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
TYPE(DIMPHYEX_t),        INTENT(IN)  :: D
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(OUT) :: PGZ_M_M ! result mass point
!
!
!*       0.2   declaration of local variables
!
INTEGER :: JIJ, JK, IIJB,IIJE,IKT
REAL, DIMENSION(D%NIJT,D%NKT) :: PA_WORK, ZWORK
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_M_M
!              --------------------
!
IF (LHOOK) CALL DR_HOOK('GZ_M_M_PHY',0,ZHOOK_HANDLE)
IIJB=D%NIJB
IIJE=D%NIJE
IKT=D%NKT
!
CALL DZM_PHY(D,PA,PA_WORK)
!

DO JK=1, IKT
  DO JIJ=IIJB, IIJE
    ZWORK(JIJ, JK)= PA_WORK(JIJ, JK) &
                                                   / PDZZ(JIJ, JK)
  END DO
END DO

!
CALL MZF_PHY(D, ZWORK, PGZ_M_M)
!
!----------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GZ_M_M_PHY',1,ZHOOK_HANDLE)
END SUBROUTINE GZ_M_M_PHY
!
!     #######################################################
      SUBROUTINE GX_M_U_PHY(D,OFLAT,PY,PDXX,PDZZ,PDZX,PGX_M_U)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODE_SHUMAN_PHY, ONLY: DXM_PHY, DZM_PHY, MXM_PHY, MZF_PHY
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),        INTENT(IN)  :: D
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PY      ! variable at the mass point
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(D%NIJT,D%NKT),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(OUT) :: PGX_M_U  ! result at flux
                                                              ! side
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1, ZWORK2, ZWORK3
INTEGER :: IIJB,IIJE,IKT,IKA,IKU
INTEGER :: JIJ,JK
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG X
!              -----------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GX_M_U',0,ZHOOK_HANDLE)
IIJB=D%NIJB
IIJE=D%NIJE
IKT=D%NKT
IKA=D%NKA
IKU=D%NKU
!
CALL DXM_PHY(D,PY,ZWORK1)
IF (.NOT. OFLAT) THEN
  CALL DZM_PHY(D,PY,ZWORK2)
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK2(JIJ, JK) = ZWORK2(JIJ, JK) / PDZZ(JIJ, JK)
    END DO
  END DO
  CALL MXM_PHY(D,ZWORK2,ZWORK3)
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK3(JIJ, JK) = ZWORK3(JIJ, JK) * PDZX(JIJ, JK)
    END DO
  END DO
  CALL MZF_PHY(D,ZWORK3,ZWORK2)
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGX_M_U(JIJ, JK) = ( ZWORK1(JIJ, JK) - ZWORK2(JIJ, JK) ) &
                               / PDXX(JIJ, JK)
    END DO
  END DO
  ! Upper level without the slope correction, lower level equal to the upper one
  DO JIJ=IIJB, IIJE
    PGX_M_U(JIJ, IKU) = ZWORK1(JIJ, IKU) / PDXX(JIJ, IKU)
    PGX_M_U(JIJ, IKA) = PGX_M_U(JIJ, IKU) ! -999.
  END DO
ELSE
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGX_M_U(JIJ, JK) = ZWORK1(JIJ, JK) / PDXX(JIJ, JK)
    END DO
  END DO
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GX_M_U',1,ZHOOK_HANDLE)
END SUBROUTINE GX_M_U_PHY
!
      SUBROUTINE GY_M_V_PHY(D,OFLAT,PY,PDYY,PDZZ,PDZY,PGY_M_V)
      USE YOMHOOK , ONLY : LHOOK, DR_HOOK, JPHOOK
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
USE MODE_SHUMAN_PHY, ONLY: DYM_PHY, DZM_PHY, MYM_PHY, MZF_PHY
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),        INTENT(IN)  :: D
LOGICAL, INTENT(IN) :: OFLAT
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PDYY                   !d*yy
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PDZY                   !d*zy
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(D%NIJT,D%NKT), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(D%NIJT,D%NKT),INTENT(OUT) :: PGY_M_V  ! result at flux
                                                              ! side
REAL, DIMENSION(D%NIJT,D%NKT) :: ZWORK1, ZWORK2, ZWORK3
INTEGER :: IIJB,IIJE,IKT,IKA,IKU
INTEGER :: JIJ,JK
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Y
!              ----------------------------
!
REAL(KIND=JPHOOK) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('GY_M_V',0,ZHOOK_HANDLE)
IIJB=D%NIJB
IIJE=D%NIJE
IKT=D%NKT
IKA=D%NKA
IKU=D%NKU
!
CALL DYM_PHY(D,PY,ZWORK1)
IF (.NOT. OFLAT) THEN
  CALL DZM_PHY(D,PY,ZWORK2)
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK2(JIJ, JK) = ZWORK2(JIJ, JK) / PDZZ(JIJ, JK)
    END DO
  END DO
  CALL MYM_PHY(D,ZWORK2,ZWORK3)
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      ZWORK3(JIJ, JK) = ZWORK3(JIJ, JK) * PDZY(JIJ, JK)
    END DO
  END DO
  CALL MZF_PHY(D,ZWORK3,ZWORK2)
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGY_M_V(JIJ, JK) = ( ZWORK1(JIJ, JK) - ZWORK2(JIJ, JK) ) &
                               / PDYY(JIJ, JK)
    END DO
  END DO
  ! Upper level without the slope correction, lower level undefined
  DO JIJ=IIJB, IIJE
    PGY_M_V(JIJ, IKU) = ZWORK1(JIJ, IKU) / PDYY(JIJ, IKU)
    PGY_M_V(JIJ, IKA) = -999.
  END DO
ELSE
  DO JK=1, IKT
    DO JIJ=IIJB, IIJE
      PGY_M_V(JIJ, JK) = ZWORK1(JIJ, JK) / PDYY(JIJ, JK)
    END DO
  END DO
ENDIF
!
!-------------------------------------------------------------------------------
!
IF (LHOOK) CALL DR_HOOK('GY_M_V',1,ZHOOK_HANDLE)
END SUBROUTINE GY_M_V_PHY
!
END MODULE MODE_GRADIENT_M_PHY
