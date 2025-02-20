!MNH_LIC Copyright 1994-2024 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_M
!     ###################### 
!
IMPLICIT NONE
INTERFACE
!
!
FUNCTION GX_M_M(PA,PDXX,PDZZ,PDZX,KKA,KKU,KL)      RESULT(PGX_M_M)
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
INTEGER, INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes (AROME)
INTEGER, INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise (AROME)
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_M_M ! result mass point
!
END FUNCTION GX_M_M
!
#ifdef MNH_OPENACC
SUBROUTINE GX_M_M_DEVICE(PA,PDXX,PDZZ,PDZX,PGX_M_M)
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)),  INTENT(OUT) :: PGX_M_M ! result mass point
!
END SUBROUTINE GX_M_M_DEVICE
#endif
!
FUNCTION GY_M_M(PA,PDYY,PDZZ,PDZY,KKA,KKU,KL)      RESULT(PGY_M_M)
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
INTEGER, INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes (AROME)
INTEGER, INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise (AROME)
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_M_M ! result mass point
!
END FUNCTION GY_M_M
!
#ifdef MNH_OPENACC
SUBROUTINE GY_M_M_DEVICE(PA,PDYY,PDZZ,PDZY,PGY_M_M)
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)),  INTENT(OUT) :: PGY_M_M ! result mass point
!
END SUBROUTINE GY_M_M_DEVICE
#endif
!
!
FUNCTION GZ_M_M(PA,PDZZ,KKA,KKU,KL)      RESULT(PGZ_M_M)
IMPLICIT NONE
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
INTEGER, INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes (AROME)
INTEGER, INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise (AROME)
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_M_M ! result mass point
!
END FUNCTION GZ_M_M
!
!
#ifdef MNH_OPENACC
SUBROUTINE GZ_M_M_DEVICE(PA,PDZZ,PGZ_M_M)
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)),  INTENT(OUT) :: PGZ_M_M ! result mass point
!
END SUBROUTINE GZ_M_M_DEVICE
#endif
!
!
FUNCTION GX_M_U(KKA,KKU,KL,PY,PDXX,PDZZ,PDZX) RESULT(PGX_M_U)
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX                   ! d*xx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX                   ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   ! d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGX_M_U  ! result at flux
                                                              ! side
END FUNCTION GX_M_U
!
!
#ifdef MNH_OPENACC
SUBROUTINE GX_M_U_DEVICE(KKA,KKU,KL,PY,PDXX,PDZZ,PDZX,PGX_M_U)
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL       ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX     ! d*xx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX     ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ     ! d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PY       ! variable at mass localization
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)),INTENT(OUT) :: PGX_M_U ! result at flux side
!
END SUBROUTINE GX_M_U_DEVICE
#endif
      FUNCTION GY_M_V(KKA,KKU,KL,PY,PDYY,PDZZ,PDZY) RESULT(PGY_M_V)
!
IMPLICIT NONE
!
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY                   !d*yy
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY                   !d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGY_M_V  ! result at flux
                                                              ! side
END FUNCTION GY_M_V
!
#ifdef MNH_OPENACC
SUBROUTINE GY_M_V_DEVICE(KKA,KKU,KL,PY,PDYY,PDZZ,PDZY,PGY_M_V)
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL       ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY     !d*yy
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY     !d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ     !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PY       ! variable at mass localization
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PGY_M_V ! result at flux side
!
END SUBROUTINE GY_M_V_DEVICE
#endif

      FUNCTION GZ_M_W(KKA, KKU, KL,PY,PDZZ) RESULT(PGZ_M_W)
!  
IMPLICIT NONE
!
                                                          ! Metric coefficient:
INTEGER,              INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGZ_M_W  ! result at flux
                                                              ! side
!
END FUNCTION GZ_M_W
!
!
#ifdef MNH_OPENACC
SUBROUTINE GZ_M_W_DEVICE(KKA,KKU,KL,PY,PDZZ,PGZ_M_W)
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL       ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ     !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PY       ! variable at mass localization
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)),INTENT(OUT) :: PGZ_M_W ! result at flux side
!
END SUBROUTINE GZ_M_W_DEVICE
#endif
!
END INTERFACE
!
END MODULE MODI_GRADIENT_M
!
!
!
!     #######################################################
      FUNCTION GX_M_M(PA,PDXX,PDZZ,PDZX,KKA,KKU,KL)      RESULT(PGX_M_M)
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
USE MODI_SHUMAN, ONLY: DXF, MZF, DZM, MXF, MXM
USE MODD_CONF, ONLY:LFLAT
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
INTEGER, INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes (AROME)
INTEGER, INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise (AROME)
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_M_M ! result mass point
!
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
IF (.NOT. LFLAT) THEN 
  PGX_M_M(:,:,:)= (DXF(MXM(PA(:,:,:)))   -                     &
                   MZF(MXF(PDZX)*DZM(PA(:,:,:)) &
                  /PDZZ(:,:,:)) )        /MXF(PDXX(:,:,:))
ELSE
  PGX_M_M(:,:,:)=DXF(MXM(PA(:,:,:))) / MXF(PDXX(:,:,:)) 
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GX_M_M
!
!
#ifdef MNH_OPENACC
!     #######################################################
      SUBROUTINE GX_M_M_DEVICE(PA,PDXX,PDZZ,PDZX,PGX_M_M)
!     #######################################################
!
!*       0.    DECLARATIONS
!
!
USE MODI_SHUMAN_DEVICE
USE MODD_CONF,          ONLY: LFLAT
!
USE MODE_MNH_ZWORK,     ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)),  INTENT(OUT) :: PGX_M_M ! result mass point
!
!
!
!*       0.2   declaration of local variables
!
INTEGER                                      :: IIU, IJU, IKU
REAL, DIMENSION(:,:,:), POINTER , CONTIGUOUS :: ZTMP1, ZTMP2, ZTMP3, ZTMP4
!
!----------------------------------------------------------------------------
IIU = SIZE( PA, 1 )
IJU = SIZE( PA, 2 )
IKU = SIZE( PA, 3 )

!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'GX_M_M' )

CALL MNH_MEM_GET( ZTMP1, IIU, IJU, IKU )
CALL MNH_MEM_GET( ZTMP2, IIU, IJU, IKU )
CALL MNH_MEM_GET( ZTMP3, IIU, IJU, IKU )
CALL MNH_MEM_GET( ZTMP4, IIU, IJU, IKU )

!$acc data present( PA, PDXX, PDZZ, PDZX, PGX_M_M, ZTMP1, ZTMP2, ZTMP3, ZTMP4)
!
!*       1.    DEFINITION of GX_M_M
!              --------------------
!
IF (.NOT. LFLAT) THEN
  CALL MXM_DEVICE( PA(:,:,:), ZTMP1(:,:,:) )
  CALL DXF_DEVICE( ZTMP1(:,:,:), ZTMP2(:,:,:) )
  CALL MXF_DEVICE( PDZX, ZTMP1(:,:,:) )
  CALL DZM_DEVICE( PA(:,:,:), ZTMP3(:,:,:) )
  !$acc kernels present(ZTMP4)
  ZTMP4(:,:,:) = ZTMP1(:,:,:) * ZTMP3(:,:,:) / PDZZ(:,:,:)
  !$acc end kernels
  CALL MZF_DEVICE( ZTMP4(:,:,:), ZTMP1(:,:,:) )
  CALL MXF_DEVICE( PDXX(:,:,:), ZTMP3(:,:,:) )
  !$acc kernels
  PGX_M_M(:,:,:)= (ZTMP2(:,:,:) - ZTMP1(:,:,:) ) / ZTMP3(:,:,:)
  !$acc end kernels
ELSE
  CALL MXM_DEVICE( PA(:,:,:), ZTMP1(:,:,:) )
  CALL DXF_DEVICE( ZTMP1(:,:,:), ZTMP2(:,:,:) )
  CALL MXF_DEVICE( PDXX(:,:,:), ZTMP3(:,:,:) )
  !$acc kernels
  PGX_M_M(:,:,:)= ZTMP2(:,:,:) / ZTMP3(:,:,:)
  !$acc end kernels
END IF
!
!$acc end data

!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'GX_M_M' )

!----------------------------------------------------------------------------
!
END SUBROUTINE GX_M_M_DEVICE
#endif

!     #######################################################
      FUNCTION GY_M_M(PA,PDYY,PDZZ,PDZY,KKA,KKU,KL)      RESULT(PGY_M_M)
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
USE MODD_CONF, ONLY:LFLAT
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
INTEGER, INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes (AROME)
INTEGER, INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise (AROME)
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_M_M ! result mass point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_M_M
!              --------------------
!
!
IF (.NOT. LFLAT) THEN 
  PGY_M_M(:,:,:)= (DYF(MYM(PA))-MZF(MYF(PDZY)*DZM(PA)&
                /PDZZ) ) /MYF(PDYY)
ELSE
  PGY_M_M(:,:,:)= DYF(MYM(PA))/MYF(PDYY)
ENDIF  
!
!----------------------------------------------------------------------------
!
END FUNCTION GY_M_M
!
!
#ifdef MNH_OPENACC
!     #######################################################
      SUBROUTINE GY_M_M_DEVICE(PA,PDYY,PDZZ,PDZY,PGY_M_M)
!     #######################################################
!
!*       0.    DECLARATIONS
!
!
USE MODD_CONF,          ONLY:LFLAT
USE MODI_SHUMAN_DEVICE
!
USE MODE_MNH_ZWORK,     ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)),  INTENT(OUT) :: PGY_M_M ! result mass point
!
!*       0.2   declaration of local variables
!
INTEGER                                     :: IIU, IJU, IKU
REAL, DIMENSION(:,:,:), POINTER, CONTIGUOUS :: ZTMP1, ZTMP2, ZTMP3, ZTMP4
!
!----------------------------------------------------------------------------
IIU = SIZE( PA, 1 )
IJU = SIZE( PA, 2 )
IKU = SIZE( PA, 3 )

!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'GY_M_M' )

CALL MNH_MEM_GET( ZTMP1, IIU, IJU, IKU )
CALL MNH_MEM_GET( ZTMP2, IIU, IJU, IKU )
CALL MNH_MEM_GET( ZTMP3, IIU, IJU, IKU )
CALL MNH_MEM_GET( ZTMP4, IIU, IJU, IKU )

!$acc data present( PA, PDYY, PDZZ, PDZY, PGY_M_M, ZTMP1, ZTMP2, ZTMP3, ZTMP4)
!
!*       1.    DEFINITION of GY_M_M
!              --------------------
!
!
IF (.NOT. LFLAT) THEN
  CALL MYM_DEVICE( PA, ZTMP1(:,:,:) )
  CALL DYF_DEVICE( ZTMP1(:,:,:), ZTMP2(:,:,:) )
  CALL MYF_DEVICE( PDZY, ZTMP1(:,:,:) )
  CALL DZM_DEVICE( PA, ZTMP3(:,:,:) )
  !$acc kernels present(ZTMP4)
  ZTMP4(:,:,:) = ZTMP1(:,:,:) * ZTMP3(:,:,:) / PDZZ(:,:,:)
  !$acc end kernels
  CALL MZF_DEVICE( ZTMP4(:,:,:), ZTMP1(:,:,:) )
  CALL MYF_DEVICE( PDYY, ZTMP3(:,:,:) )
  !$acc kernels
  PGY_M_M(:,:,:)= ( ZTMP2(:,:,:) - ZTMP1(:,:,:) ) / ZTMP3(:,:,:)
  !$acc end kernels
ELSE
  CALL MYM_DEVICE( PA, ZTMP1(:,:,:) )
  CALL DYF_DEVICE( ZTMP1(:,:,:), ZTMP2(:,:,:) )
  CALL MYF_DEVICE( PDYY, ZTMP1(:,:,:) )
  !$acc kernels
  PGY_M_M(:,:,:)= ZTMP2(:,:,:) / ZTMP1(:,:,:)
  !$acc end kernels
ENDIF
!
!$acc end data

!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'GY_M_M' )

!----------------------------------------------------------------------------
!
END SUBROUTINE GY_M_M_DEVICE
#endif
!
!
!     #############################################
      FUNCTION GZ_M_M(PA,PDZZ)      RESULT(PGZ_M_M)
!     #############################################
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
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_M_M ! result mass point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_M_M
!              --------------------
!
PGZ_M_M(:,:,:)= MZF( DZM(PA(:,:,:))/PDZZ(:,:,:) )
!
!----------------------------------------------------------------------------
!
END FUNCTION GZ_M_M
!
!
#ifdef MNH_OPENACC
!     #######################################################
      SUBROUTINE GZ_M_M_DEVICE(PA,PDZZ,PGZ_M_M)
!     #######################################################
!
!*       0.    DECLARATIONS
!
!
USE MODI_SHUMAN_DEVICE
!
USE MODE_MNH_ZWORK,     ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
use mode_msg
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the mass point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)),  INTENT(OUT) :: PGZ_M_M ! result mass point
!
!*       0.2   declaration of local variables
!
INTEGER                                      :: IIU, IJU, IKU
REAL, DIMENSION(:,:,:), POINTER , CONTIGUOUS :: ZTMP1, ZTMP2, ZTMP3, ZTMP4
!
!----------------------------------------------------------------------------
IIU = SIZE( PA, 1 )
IJU = SIZE( PA, 2 )
IKU = SIZE( PA, 3 )

!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'GZ_M_M' )

CALL MNH_MEM_GET( ZTMP1, IIU, IJU, IKU )
CALL MNH_MEM_GET( ZTMP2, IIU, IJU, IKU )

!$acc data present( PA, PDZZ, PGZ_M_M, ZTMP1, ZTMP2 )

call Print_msg( NVERB_WARNING, 'GEN', 'GZ_M_M_DEVICE', 'OpenACC: not yet tested' )
!
!*       1.    DEFINITION of GZ_M_M
!              --------------------
!
CALL DZM_DEVICE( PA(:,:,:), ZTMP1(:,:,:) )
!$acc kernels present_cr(ZTMP2)
ZTMP2(:,:,:) = ZTMP1(:,:,:) / PDZZ(:,:,:)
!$acc end kernels
CALL MZF_DEVICE( ZTMP2(:,:,:), PGZ_M_M(:,:,:) )

!$acc end data

!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'GZ_M_M' )

!----------------------------------------------------------------------------
!
END SUBROUTINE GZ_M_M_DEVICE
#endif
!
!
!     ##################################################
      FUNCTION GX_M_U(KKA,KKU,KL,PY,PDXX,PDZZ,PDZX) RESULT(PGX_M_U)
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
!!    P. Hereil and J. Stein       * Meteo France *
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
USE MODI_SHUMAN
USE MODD_CONF, ONLY:LFLAT
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!              ------------------------------------
!
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX                   ! d*xx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX                   ! d*zx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   ! d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGX_M_U  ! result at flux
                                                              ! side
INTEGER  IIU,IKU,JI,JK
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
IIU=SIZE(PY,1)
IJU=SIZE(PY,2)
IKU=SIZE(PY,3)
IF (.NOT. LFLAT) THEN
! PGX_M_U = (   DXM(PY)  -  MZF (   MXM(  DZM(PY) /PDZZ  ) * PDZX   )   )/PDXX
!!  DO JK=1+JPVEXT_TURB,IKU-JPVEXT_TURB
!!    DO JI=1+JPHEXT,IIU
!!        PGX_M_U(JI,:,JK)=                                                 &
!!           (  PY(JI,:,JK)-PY(JI-1,:,JK)                                   &
!!             -(  (PY(JI,:,JK)-PY(JI,:,JK-1))     / PDZZ(JI,:,JK)          &
!!                +(PY(JI-1,:,JK)-PY(JI-1,:,JK-1)) / PDZZ(JI-1,:,JK)        &
!!              ) * PDZX(JI,:,JK)* 0.25                                     &
!!             -(  (PY(JI,:,JK+1)-PY(JI,:,JK))     / PDZZ(JI,:,JK+1)        &
!!                +(PY(JI-1,:,JK+1)-PY(JI-1,:,JK)) / PDZZ(JI-1,:,JK+1)      &
!!              ) * PDZX(JI,:,JK+1)* 0.25                                   &
!!            )  / PDXX(JI,:,JK)
!!    END DO
!!  END DO
  JIJKOR  = 1 + JPHEXT + IIU*IJU*(JPVEXT_TURB+1 - 1)
  JIJKEND = IIU*IJU*(IKU-JPVEXT_TURB)
!CDIR NODEP
!OCL NOVREC
  DO JIJK=JIJKOR , JIJKEND
! indexation
    JI_1JK   = JIJK - 1
    JIJK_1   = JIJK     - IIU*IJU*KL
    JI_1JK_1 = JIJK - 1 - IIU*IJU*KL
    JIJKP1   = JIJK     + IIU*IJU*KL
    JI_1JKP1 = JIJK - 1 + IIU*IJU*KL
!
    PGX_M_U(JIJK,1,1)=                                              &
       (  PY(JIJK,1,1)-PY(JI_1JK,1,1)                               &
       -(  (PY(JIJK,1,1)-PY(JIJK_1,1,1))     / PDZZ(JIJK,1,1)       &
       +(PY(JI_1JK,1,1)-PY(JI_1JK_1,1,1)) / PDZZ(JI_1JK,1,1)        &
       ) * PDZX(JIJK,1,1)* 0.25                                     &
       -(  (PY(JIJKP1,1,1)-PY(JIJK,1,1))     / PDZZ(JIJKP1,1,1)     &
       +(PY(JI_1JKP1,1,1)-PY(JI_1JK,1,1)) / PDZZ(JI_1JKP1,1,1)      &
       ) * PDZX(JIJKP1,1,1)* 0.25                                   &
        )  / PDXX(JIJK,1,1)
  END DO

!
  DO JI=1+JPHEXT,IIU
    PGX_M_U(JI,:,KKU)=  ( PY(JI,:,KKU)-PY(JI-1,:,KKU)  )  / PDXX(JI,:,KKU)
    PGX_M_U(JI,:,KKA)=   PGX_M_U(JI,:,KKU) ! -999.
  END DO
ELSE
!  PGX_M_U = DXM(PY) / PDXX
  PGX_M_U(1+1:IIU,:,:) = ( PY(1+1:IIU,:,:)-PY(1:IIU-1,:,:) ) & ! +JPHEXT
                             / PDXX(1+1:IIU,:,:)
ENDIF
DO JI=1,JPHEXT
  PGX_M_U(JI,:,:)=PGX_M_U(IIU-2*JPHEXT+JI,:,:) ! for reprod JPHEXT <> 1
END DO  
!
!-------------------------------------------------------------------------------
!
END FUNCTION GX_M_U
!
!
#ifdef MNH_OPENACC
!     ##################################################
      SUBROUTINE GX_M_U_DEVICE(KKA,KKU,KL,PY,PDXX,PDZZ,PDZX,PGX_M_U)
!     ##################################################
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF, ONLY:LFLAT
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and result
!              ------------------------------------
!
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL       ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDXX     ! d*xx
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZX     ! d*zx 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ     ! d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PY       ! variable at mass localization
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)),INTENT(OUT) :: PGX_M_U ! result at flux side
!
INTEGER  IIU,IJU,IKU,JI
!
!-------------------------------------------------------------------------------

!$acc data present( PY, PDXX, PDZZ, PDZX, PGX_M_U )
!
!*       1.    COMPUTE THE GRADIENT ALONG X
!              -----------------------------
!

IIU=SIZE(PY,1)
IJU=SIZE(PY,2)
IKU=SIZE(PY,3)
!$acc kernels !!!!present_cr(PGX_M_U)
IF (.NOT. LFLAT) THEN 
PGX_M_U(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB) = &
  (  PY(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)-PY(JPHEXT:IIU-1,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)      &
    -(  (PY(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)-PY(1+JPHEXT:IIU,1:IJU,JPVEXT_TURB:IKU-JPVEXT_TURB-1)) &
    / PDZZ(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)          &
      +(PY(JPHEXT:IIU-1,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)-PY(JPHEXT:IIU-1,1:IJU,JPVEXT_TURB:IKU-JPVEXT_TURB-1))  &
        / PDZZ(JPHEXT:IIU-1,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)        &
    ) * PDZX(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)* 0.25                                     &
    -(  (PY(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB+1:IKU-JPVEXT_TURB+1)-PY(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)) &
    / PDZZ(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB+1:IKU-JPVEXT_TURB+1)        &
      +(PY(JPHEXT:IIU-1,1:IJU,1+JPVEXT_TURB+1:IKU-JPVEXT_TURB+1)-PY(JPHEXT:IIU-1,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)) &
       / PDZZ(JPHEXT:IIU-1,1:IJU,1+JPVEXT_TURB+1:IKU-JPVEXT_TURB+1)      &
    ) * PDZX(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB+1:IKU-JPVEXT_TURB+1)* 0.25                                   &
  )  / PDXX(1+JPHEXT:IIU,1:IJU,1+JPVEXT_TURB:IKU-JPVEXT_TURB)
!
  DO JI=1+JPHEXT,IIU
    PGX_M_U(JI,:,KKU)=  ( PY(JI,:,KKU)-PY(JI-1,:,KKU)  )  / PDXX(JI,:,KKU) 
    PGX_M_U(JI,:,KKA)=   PGX_M_U(JI,:,KKU) ! -999.
  END DO
ELSE
  PGX_M_U(1+1:IIU,:,:) = ( PY(1+1:IIU,:,:)-PY(1:IIU-1,:,:) ) & ! +JPHEXT
                             / PDXX(1+1:IIU,:,:)
ENDIF  
DO JI=1,JPHEXT
  PGX_M_U(JI,:,:)=PGX_M_U(IIU-2*JPHEXT+JI,:,:) ! for reprod JPHEXT <> 1
END DO  
!$acc end kernels

!$acc end data

!-------------------------------------------------------------------------------
!
END SUBROUTINE GX_M_U_DEVICE
#endif
!
!
!     ##################################################
      FUNCTION GY_M_V(KKA,KKU,KL,PY,PDYY,PDZZ,PDZY) RESULT(PGY_M_V)
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
!!    P. Hereil and J. Stein       * Meteo France *
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
USE MODI_SHUMAN
USE MODD_CONF, ONLY:LFLAT
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!              -------------------------------------
!
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY                   !d*yy
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY                   !d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGY_M_V  ! result at flux
                                                              ! side
INTEGER  IJU,IKU,JJ,JK
!
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Y
!              ----------------------------
!
IJU=SIZE(PY,2)
IKU=SIZE(PY,3)
IF (.NOT. LFLAT) THEN 
!  PGY_M_V = (   DYM(PY)  -  MZF (   MYM(  DZM(PY) /PDZZ  ) * PDZY   )   )/PDYY
  DO JK=1+JPVEXT_TURB,IKU-JPVEXT_TURB
    DO JJ=1+JPHEXT,IJU
        PGY_M_V(:,JJ,JK)=                                                 &
           (  PY(:,JJ,JK)-PY(:,JJ-1,JK)                                   &
             -(  (PY(:,JJ,JK)-PY(:,JJ,JK-KL))     / PDZZ(:,JJ,JK)          &
                +(PY(:,JJ-1,JK)-PY(:,JJ-KL,JK-KL)) / PDZZ(:,JJ-1,JK)        &
              ) * PDZY(:,JJ,JK)* 0.25                                     &
             -(  (PY(:,JJ,JK+KL)-PY(:,JJ,JK))     / PDZZ(:,JJ,JK+KL)        &
                +(PY(:,JJ-1,JK+KL)-PY(:,JJ-1,JK)) / PDZZ(:,JJ-1,JK+KL)      &
              ) * PDZY(:,JJ,JK+KL)* 0.25                                   &
            )  / PDYY(:,JJ,JK)
    END DO
  END DO
!
  DO JJ=1+JPHEXT,IJU
    PGY_M_V(:,JJ,KKU)=  ( PY(:,JJ,KKU)-PY(:,JJ-1,KKU)  )  / PDYY(:,JJ,KKU) 
    PGY_M_V(:,JJ,KKA)=  PGY_M_V(:,JJ,KKU) ! -999.
  END DO
ELSE
!  PGY_M_V = DYM(PY)/PDYY
  PGY_M_V(:,1+1:IJU,:) = ( PY(:,1+1:IJU,:)-PY(:,1:IJU-1,:) ) & ! +JPHEXT
                               / PDYY(:,1+1:IJU,:)
ENDIF  
DO JJ=1,JPHEXT
  PGY_M_V(:,JJ,:)=PGY_M_V(:,IJU-2*JPHEXT+JJ,:)
END DO
!
!-------------------------------------------------------------------------------
!
END FUNCTION GY_M_V
!
!
#ifdef MNH_OPENACC
!     ##################################################
      SUBROUTINE GY_M_V_DEVICE(KKA,KKU,KL,PY,PDYY,PDZZ,PDZY,PGY_M_V)
!     ##################################################
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_CONF, ONLY:LFLAT
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!              -------------------------------------
!
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL       ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDYY     !d*yy
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZY     !d*zy 
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ     !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PY       ! variable at mass localization
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PGY_M_V  ! result at flux side
!
INTEGER  IIU,IJU,IKU,JI,JJ,JK
!
!-------------------------------------------------------------------------------

!$acc data present( PY, PDYY, PDZZ, PDZY, PGY_M_V )
!
!*       1.    COMPUTE THE GRADIENT ALONG Y
!              ----------------------------
!
IIU=SIZE(PY,1)
IJU=SIZE(PY,2)
IKU=SIZE(PY,3)
!$acc kernels present_cr(PGY_M_V)
IF (.NOT. LFLAT) THEN 
  DO JK=1+JPVEXT_TURB,IKU-JPVEXT_TURB
    DO JJ=1+JPHEXT,IJU
        PGY_M_V(:,JJ,JK)=                                                 &
           (  PY(:,JJ,JK)-PY(:,JJ-1,JK)                                   &
             -(  (PY(:,JJ,JK)-PY(:,JJ,JK-KL))     / PDZZ(:,JJ,JK)          &
                +(PY(:,JJ-1,JK)-PY(:,JJ-KL,JK-KL)) / PDZZ(:,JJ-1,JK)        &
              ) * PDZY(:,JJ,JK)* 0.25                                     &
             -(  (PY(:,JJ,JK+KL)-PY(:,JJ,JK))     / PDZZ(:,JJ,JK+KL)        &
                +(PY(:,JJ-1,JK+KL)-PY(:,JJ-1,JK)) / PDZZ(:,JJ-1,JK+KL)      &
              ) * PDZY(:,JJ,JK+KL)* 0.25                                   &
            )  / PDYY(:,JJ,JK)
    END DO
  END DO
!
  DO JJ=1+JPHEXT,IJU
    PGY_M_V(:,JJ,KKU)=  ( PY(:,JJ,KKU)-PY(:,JJ-1,KKU)  )  / PDYY(:,JJ,KKU) 
    PGY_M_V(:,JJ,KKA)=  PGY_M_V(:,JJ,KKU) ! -999.
  END DO
ELSE
  PGY_M_V(:,1+1:IJU,:) = ( PY(:,1+1:IJU,:)-PY(:,1:IJU-1,:) ) & ! +JPHEXT
                               / PDYY(:,1+1:IJU,:)
ENDIF  
DO JJ=1,JPHEXT
  !$mnh_do_concurrent ( JI=1:IIU , JK=1:IKU ) 
     PGY_M_V(JI,JJ,JK)=PGY_M_V(JI,IJU-2*JPHEXT+JJ,JK)
  !$mnh_end_do()
END DO
!$acc end kernels

!$acc end data

!-------------------------------------------------------------------------------
!
END SUBROUTINE GY_M_V_DEVICE
#endif
!
!
!     #########################################
      FUNCTION GZ_M_W(KKA,KKU,KL,PY,PDZZ) RESULT(PGZ_M_W)
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
!!    P. Hereil and J. Stein       * Meteo France *
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
USE MODI_SHUMAN
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!              -------------------------------------
!
INTEGER,           INTENT(IN)     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,           INTENT(IN)     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise

                                                          ! Metric coefficient:
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ                   !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)                :: PY       ! variable at mass
                                                              ! localization
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)) :: PGZ_M_W  ! result at flux
                                                              ! side
!
INTEGER :: IKT,IKTB,IKTE
!-------------------------------------------------------------------------------
!
!*       1.    COMPUTE THE GRADIENT ALONG Z
!              -----------------------------
!
IKT=SIZE(PY,3)
IKTB=1+JPVEXT_TURB              
IKTE=IKT-JPVEXT_TURB

PGZ_M_W(:,:,IKTB:IKTE) =  (PY(:,:,IKTB:IKTE)-PY(:,:,IKTB-KL:IKTE-KL))  &
                           / PDZZ(:,:,IKTB:IKTE)
PGZ_M_W(:,:,KKU)=  (PY(:,:,KKU)-PY(:,:,KKU-KL))  &
                           / PDZZ(:,:,KKU)
PGZ_M_W(:,:,KKA)= PGZ_M_W(:,:,KKU) ! -999.
!
!-------------------------------------------------------------------------------
!
END FUNCTION GZ_M_W
!
!
#ifdef MNH_OPENACC
!     #########################################
      SUBROUTINE GZ_M_W_DEVICE(KKA,KKU,KL,PY,PDZZ,PGZ_M_W)
!     #########################################
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
IMPLICIT NONE
!
!*       0.1   Declarations of arguments and results
!              -------------------------------------
!
INTEGER,                INTENT(IN)  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                INTENT(IN)  :: KL       ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PDZZ     !d*zz
!
REAL, DIMENSION(:,:,:), INTENT(IN)  :: PY       ! variable at mass localization
!
REAL, DIMENSION(SIZE(PY,1),SIZE(PY,2),SIZE(PY,3)),INTENT(OUT) :: PGZ_M_W ! result at flux side
!
INTEGER :: IKT,IKTB,IKTE
!-------------------------------------------------------------------------------

!$acc data present( PY, PDZZ, PGZ_M_W )
!
!*       1.    COMPUTE THE GRADIENT ALONG Z
!              -----------------------------
!
IKT=SIZE(PY,3)
IKTB=1+JPVEXT_TURB              
IKTE=IKT-JPVEXT_TURB

!$acc kernels
PGZ_M_W(:,:,IKTB:IKTE) =  (PY(:,:,IKTB:IKTE)-PY(:,:,IKTB-KL:IKTE-KL))  &
                           / PDZZ(:,:,IKTB:IKTE)
PGZ_M_W(:,:,KKU)=  (PY(:,:,KKU)-PY(:,:,KKU-KL))  &
                           / PDZZ(:,:,KKU)
PGZ_M_W(:,:,KKA)= PGZ_M_W(:,:,KKU) ! -999.
!$acc end kernels

!$acc end data

!-------------------------------------------------------------------------------
!
END SUBROUTINE GZ_M_W_DEVICE
#endif
