!MNH_LIC Copyright 1994-2025 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_V
!     ######################
!
IMPLICIT NONE
INTERFACE
!
!           
FUNCTION GY_V_M(PA,PDYY,PDZZ,PDZY, KKA, KKU, KL)      RESULT(PGY_V_M)
IMPLICIT NONE
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_V_M ! result mass point
!
END FUNCTION GY_V_M
!
!
#ifdef MNH_OPENACC
SUBROUTINE GY_V_M_DEVICE(PA,PDYY,PDZZ,PDZY,PGY_V_M_DEVICE)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA       ! variable at the V point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDYY     ! metric coefficient dyy
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ     ! metric coefficient dzz
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZY     ! metric coefficient dzy
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PGY_V_M_DEVICE ! result mass point
!
END SUBROUTINE GY_V_M_DEVICE
#endif
!
!           
FUNCTION GX_V_UV(PA,PDXX,PDZZ,PDZX, KKA, KKU, KL)      RESULT(PGX_V_UV)
IMPLICIT NONE
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_V_UV ! result UV point
!
END FUNCTION GX_V_UV
!
!
SUBROUTINE GX_V_UV_DEVICE(PA,PDXX,PDZZ,PDZX,PGX_V_UV_DEVICE)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA       ! variable at the V point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDXX     ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ     ! metric coefficient dzz
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZX     ! metric coefficient dzx
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PGX_V_UV_DEVICE ! result UV point
!
END SUBROUTINE GX_V_UV_DEVICE
!
!           
FUNCTION GZ_V_VW(PA,PDZZ, KKA, KKU, KL)      RESULT(PGZ_V_VW)
IMPLICIT NONE
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_V_VW ! result VW point
!
END FUNCTION GZ_V_VW
!
!
#ifdef MNH_OPENACC
SUBROUTINE GZ_V_VW_DEVICE(PA,PDZZ,PGZ_V_VW_DEVICE)
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA       ! variable at the V point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ     ! metric coefficient dzz
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PGZ_V_VW_DEVICE ! result VW point
!
END SUBROUTINE GZ_V_VW_DEVICE
#endif
!
!
END INTERFACE
!
END MODULE MODI_GRADIENT_V
!
!     #######################################################
      FUNCTION GY_V_M(PA,PDYY,PDZZ,PDZY, KKA, KKU, KL)      RESULT(PGY_V_M)
!     #######################################################
!
!!****  *GY_V_M* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Y
!!                          direction for a variable placed at the 
!!                          V point and the result is placed at
!!                          the mass point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Y cartesian direction for a field PA placed at the 
!     V point. The result is placed at the mass point.
!
!
!                       (          ______________z )
!                       (          (___________y ) )
!                    1  (          (d*zy dzm(PA) ) ) 
!      PGY_V_M =   ---- (dyf(PA) - (------------)) )
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
!!      MYF,MZF         : Shuman functions (mean operators)
!!      DYF,DZF         : Shuman functions (finite difference operators)
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
USE MODI_SHUMAN
USE MODD_CONF
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
INTEGER,                 INTENT(IN),OPTIONAL  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                 INTENT(IN),OPTIONAL  :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_V_M ! result mass point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_V_M
!              --------------------
!
IF (.NOT. LFLAT) THEN
  PGY_V_M(:,:,:)= (DYF(PA)        -                      &
                   MZF( MYF(PDZY*DZM(PA))/PDZZ )         &
                  ) / MYF(PDYY)
ELSE
  PGY_V_M(:,:,:)= DYF(PA) / MYF(PDYY)
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GY_V_M
!
!
#ifdef MNH_OPENACC
!     #######################################################
      SUBROUTINE GY_V_M_DEVICE(PA,PDYY,PDZZ,PDZY,PGY_V_M_DEVICE)
!     #######################################################
!
!*       0.    DECLARATIONS
!
!
USE MODI_SHUMAN_DEVICE
USE MODD_CONF
USE MODE_MNH_ZWORK, ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA       ! variable at the V point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDYY     ! metric coefficient dyy
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ     ! metric coefficient dzz
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZY     ! metric coefficient dzy
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PGY_V_M_DEVICE ! result mass point
!
REAL, DIMENSION(:,:,:), pointer , contiguous ::  ZTMP1_DEVICE, ZTMP2_DEVICE, ZTMP3_DEVICE
!
INTEGER  :: JIU,JJU,JKU
INTEGER  :: JI,JJ,JK
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------

!$acc data present_crm( PA, PDYY, PDZZ, PDZY, PGY_V_M_DEVICE )

JIU =  size(pa, 1 )
JJU =  size(pa, 2 )
JKU =  size(pa, 3 )

!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'GY_V_M' )

CALL MNH_MEM_GET( ztmp1_device, JIU, JJU, JKU )
CALL MNH_MEM_GET( ztmp2_device, JIU, JJU, JKU )
CALL MNH_MEM_GET( ztmp3_device, JIU, JJU, JKU )

!$acc data present_crm( ztmp1_device, ztmp2_device, ztmp3_device )

!
!*       1.    DEFINITION of GY_V_M_DEVICE
!              --------------------
!
IF (.NOT. LFLAT) THEN
  CALL DYF_DEVICE(PA,ZTMP1_DEVICE)
  CALL DZM_DEVICE( PA, ZTMP2_DEVICE )
  !$acc kernels
  !$mnh_do_concurrent ( JI=1:JIU,JJ=1:JJU,JK=1:JKU)
     ZTMP3_DEVICE(JI,JJ,JK) = PDZY(JI,JJ,JK)*ZTMP2_DEVICE(JI,JJ,JK)
  !$mnh_end_do() !CONCURRENT
  !$acc end kernels
  CALL MYF_DEVICE(ZTMP3_DEVICE,ZTMP2_DEVICE)
  !$acc kernels
  !$mnh_do_concurrent ( JI=1:JIU,JJ=1:JJU,JK=1:JKU)
     ZTMP3_DEVICE(JI,JJ,JK) = ZTMP2_DEVICE(JI,JJ,JK)/PDZZ(JI,JJ,JK)
  !$mnh_end_do() !CONCURRENT
  !$acc end kernels
  CALL MZF_DEVICE( ZTMP3_DEVICE, ZTMP2_DEVICE )
  CALL MYF_DEVICE(PDYY,ZTMP3_DEVICE)
!$acc kernels
  PGY_V_M_DEVICE(:,:,:)= (ZTMP1_DEVICE(:,:,:) - ZTMP2_DEVICE(:,:,:)) / ZTMP3_DEVICE(:,:,:)
!$acc end kernels
ELSE
  CALL DYF_DEVICE(PA,ZTMP1_DEVICE)
  CALL MYF_DEVICE(PDYY,ZTMP2_DEVICE)
!$acc kernels
  PGY_V_M_DEVICE(:,:,:)= ZTMP1_DEVICE(:,:,:) / ZTMP2_DEVICE(:,:,:)
!$acc end kernels
END IF

!$acc end data

!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'GY_V_M' )

!$acc end data

!----------------------------------------------------------------------------
!
END SUBROUTINE GY_V_M_DEVICE
#endif
!
! 
!     #########################################################
      FUNCTION GX_V_UV(PA,PDXX,PDZZ,PDZX, KKA, KKU, KL)      RESULT(PGX_V_UV)
!     #########################################################
!
!!****  *GX_V_UV* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian X
!!                          direction for a variable placed at the 
!!                          V point and the result is placed at
!!                          the UV vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the X cartesian direction for a field PA placed at the 
!     V point. The result is placed at the UV vorticity point.
!
!
!                       (          _________________z )
!                       (          (___y _________x ) )
!                    1  (          (d*zx (dzm(PA))) ) )
!      PGX_V_UV=   ---- (dxm(PA) - (     (------  ) ) )
!                  ___y (          (     ( ___y   ) ) )
!                  d*xx (          (     ( d*zz   ) ) )    
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
!!      MXM,MZF,MYM     : Shuman functions (mean operators)
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
USE MODI_SHUMAN
USE MODD_CONF
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
INTEGER,                 INTENT(IN),OPTIONAL  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                 INTENT(IN),OPTIONAL  :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_V_UV ! result UV point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GX_V_UV
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  PGX_V_UV(:,:,:)= ( DXM(PA)- MZF( MXM( DZM(PA)/&
                    MYM(PDZZ) ) *MYM(PDZX) )   )   / MYM(PDXX)
ELSE
  PGX_V_UV(:,:,:)= DXM(PA) / MYM(PDXX)
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GX_V_UV
!
!
!     #########################################################
      SUBROUTINE GX_V_UV_DEVICE(PA,PDXX,PDZZ,PDZX,PGX_V_UV_DEVICE)
!     #########################################################
!
!*       0.    DECLARATIONS
!
!
USE MODI_SHUMAN_DEVICE
USE MODD_CONF
USE MODE_MNH_ZWORK, ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA       ! variable at the V point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDXX     ! metric coefficient dxx
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ     ! metric coefficient dzz
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZX     ! metric coefficient dzx
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PGX_V_UV_DEVICE ! result UV point
!
!*       0.2   declaration of local variables
!
REAL, DIMENSION(:,:,:), pointer , contiguous ::  ZTMP1_DEVICE, ZTMP2_DEVICE, ZTMP3_DEVICE, ZTMP4_DEVICE
!
INTEGER  :: JIU,JJU,JKU
INTEGER  :: JI,JJ,JK
!
!----------------------------------------------------------------------------

!$acc data present_crm( PA, PDXX, PDZZ, PDZX, PGX_V_UV_DEVICE )

JIU =  size(pa, 1 )
JJU =  size(pa, 2 )
JKU =  size(pa, 3 )

!Pin positions in the pools of MNH memory
#ifdef MNH_OPENACC
CALL MNH_MEM_POSITION_PIN( 'GX_V_UV' )

CALL MNH_MEM_GET( ztmp1_device, JIU, JJU, JKU )
CALL MNH_MEM_GET( ztmp2_device, JIU, JJU, JKU )
CALL MNH_MEM_GET( ztmp3_device, JIU, JJU, JKU )
CALL MNH_MEM_GET( ztmp4_device, JIU, JJU, JKU )
#else
ALLOCATE(ztmp1_device(JIU, JJU, JKU))
ALLOCATE(ztmp2_device(JIU, JJU, JKU))
ALLOCATE(ztmp3_device(JIU, JJU, JKU))
ALLOCATE(ztmp4_device(JIU, JJU, JKU))
#endif 
!$acc data present_crm( ztmp1_device, ztmp2_device, ztmp3_device, ztmp4_device )

!
!*       1.    DEFINITION of GX_V_UV_DEVICE
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  CALL DXM_DEVICE(PA,ZTMP1_DEVICE)
  CALL MYM_DEVICE(PDZZ,ZTMP2_DEVICE)
  CALL DZM_DEVICE( PA, ZTMP3_DEVICE )
  !$acc kernels
  !$mnh_do_concurrent ( JI=1:JIU,JJ=1:JJU,JK=1:JKU)
     ZTMP4_DEVICE(JI,JJ,JK) = ZTMP3_DEVICE(JI,JJ,JK) / ZTMP2_DEVICE(JI,JJ,JK)
  !$mnh_end_do() !CONCURRENT   
  !$acc end kernels
  CALL MXM_DEVICE(ZTMP4_DEVICE,ZTMP2_DEVICE)
  CALL MYM_DEVICE(PDZX,ZTMP3_DEVICE)
  !$acc kernels
  !$mnh_do_concurrent ( JI=1:JIU,JJ=1:JJU,JK=1:JKU)
     ZTMP4_DEVICE(JI,JJ,JK) = ZTMP2_DEVICE(JI,JJ,JK) *ZTMP3_DEVICE(JI,JJ,JK)
  !$mnh_end_do() !CONCURRENT   
  !$acc end kernels
  CALL MZF_DEVICE( ZTMP4_DEVICE, ZTMP2_DEVICE )
  CALL MYM_DEVICE(PDXX,ZTMP3_DEVICE)
  !$acc kernels
  !$mnh_do_concurrent ( JI=1:JIU,JJ=1:JJU,JK=1:JKU)
     PGX_V_UV_DEVICE(JI,JJ,JK)= ( ZTMP1_DEVICE(JI,JJ,JK) - ZTMP2_DEVICE(JI,JJ,JK) ) / ZTMP3_DEVICE(JI,JJ,JK)
  !$mnh_end_do() !CONCURRENT   
  !$acc end kernels
ELSE
  CALL DXM_DEVICE(PA,ZTMP1_DEVICE)
  CALL MYM_DEVICE(PDXX,ZTMP2_DEVICE)
!$acc kernels
  PGX_V_UV_DEVICE(:,:,:)= ZTMP1_DEVICE(:,:,:) / ZTMP2_DEVICE(:,:,:)
!$acc end kernels
END IF

!$acc end data

!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
#ifdef MNH_OPENACC
CALL MNH_MEM_RELEASE( 'GX_V_UV' )
#else
DEALLOCATE(ZTMP1_DEVICE, ZTMP2_DEVICE, ZTMP3_DEVICE, ZTMP4_DEVICE)
#endif

!$acc end data

!----------------------------------------------------------------------------
!
END SUBROUTINE GX_V_UV_DEVICE
!
!
!     #######################################################
      FUNCTION GZ_V_VW(PA,PDZZ, KKA, KKU, KL)      RESULT(PGZ_V_VW)
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
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_V_VW ! result VW point
!
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
PGZ_V_VW(:,:,:)= DZM(PA) / MYM(PDZZ)
!
!----------------------------------------------------------------------------
!
END FUNCTION GZ_V_VW
!
!
#ifdef MNH_OPENACC
!     #######################################################
      SUBROUTINE GZ_V_VW_DEVICE(PA,PDZZ,PGZ_V_VW_DEVICE)
!     #######################################################
!
!*       0.    DECLARATIONS
!
!
USE MODI_SHUMAN_DEVICE
USE MODI_SHUMAN
USE MODE_MNH_ZWORK, ONLY: MNH_MEM_GET, MNH_MEM_POSITION_PIN, MNH_MEM_RELEASE
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
REAL, DIMENSION(:,:,:), INTENT(IN) :: PA       ! variable at the V point
REAL, DIMENSION(:,:,:), INTENT(IN) :: PDZZ     ! metric coefficient dzz
!
REAL, DIMENSION(:,:,:), INTENT(OUT) :: PGZ_V_VW_DEVICE ! result VW point
!
!*       0.2   declaration of local variables
!
REAL, DIMENSION(:,:,:), pointer , contiguous :: ZTMP1_DEVICE, ZTMP2_DEVICE
!
INTEGER  :: JIU,JJU,JKU
!----------------------------------------------------------------------------

!$acc data present_crm( PA, PDZZ, PGZ_V_VW_DEVICE )

JIU =  size(pa, 1 )
JJU =  size(pa, 2 )
JKU =  size(pa, 3 )

!Pin positions in the pools of MNH memory
CALL MNH_MEM_POSITION_PIN( 'GZ_V_VW' )

CALL MNH_MEM_GET( ztmp1_device, JIU, JJU, JKU )
CALL MNH_MEM_GET( ztmp2_device, JIU, JJU, JKU )

!$acc data present_crm( ztmp1_device, ztmp2_device )

!
!*       1.    DEFINITION of GZ_V_VW_DEVICE
!              ---------------------
!
CALL DZM_DEVICE( PA, ZTMP1_DEVICE )
CALL MYM_DEVICE(PDZZ,ZTMP2_DEVICE)
!$acc kernels
PGZ_V_VW_DEVICE(:,:,:)= ZTMP1_DEVICE(:,:,:) / ZTMP2_DEVICE(:,:,:)
!$acc end kernels

!$acc end data

!Release all memory allocated with MNH_MEM_GET calls since last call to MNH_MEM_POSITION_PIN
CALL MNH_MEM_RELEASE( 'GZ_V_VW' )

!$acc end data

!----------------------------------------------------------------------------
!
END SUBROUTINE GZ_V_VW_DEVICE
#endif
