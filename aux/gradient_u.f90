!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_U
!     ######################
!
IMPLICIT NONE
INTERFACE
!
!     
FUNCTION GX_U_M(PA,PDXX,PDZZ,PDZX, KKA, KKU, KL)      RESULT(PGX_U_M)
IMPLICIT NONE
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_U_M ! result mass point
!
END FUNCTION GX_U_M
!
!     
FUNCTION GY_U_UV(PA,PDYY,PDZZ,PDZY, KKA, KKU, KL)      RESULT(PGY_U_UV)
IMPLICIT NONE
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_U_UV ! result UV point
!
END FUNCTION GY_U_UV
!
!     
FUNCTION GZ_U_UW(PA,PDZZ, KKA, KKU, KL)      RESULT(PGZ_U_UW)
IMPLICIT NONE
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_U_UW ! result UW point
!
END FUNCTION GZ_U_UW
!
END INTERFACE
!
END MODULE MODI_GRADIENT_U
!
!
!
!
!     #######################################################
      FUNCTION GX_U_M(PA,PDXX,PDZZ,PDZX, KKA, KKU, KL)      RESULT(PGX_U_M)
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
USE MODI_SHUMAN
USE MODD_CONF
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
INTEGER,                 INTENT(IN),OPTIONAL   :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                 INTENT(IN),OPTIONAL   :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_U_M ! result mass point
!
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
IF (.NOT. LFLAT) THEN
  PGX_U_M(:,:,:)= ( DXF(PA)        -                 &
                    MZF(MXF(PDZX*DZM(PA)) / PDZZ )  &
                  ) / MXF(PDXX)
ELSE
  PGX_U_M(:,:,:)= DXF(PA) /  MXF(PDXX)
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GX_U_M
!
! 
!     #########################################################
      FUNCTION GY_U_UV(PA,PDYY,PDZZ,PDZY, KKA, KKU, KL)      RESULT(PGY_U_UV)
!     #########################################################
!
!!****  *GY_U_UV* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Y
!!                          direction for a variable placed at the 
!!                          U point and the result is placed at
!!                          the UV vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Y cartesian direction for a field PA placed at the 
!     U point. The result is placed at the UV vorticity point.
!
!
!
!                       (          _________________z )
!                       (          (___x _________y ) )
!                    1  (          (d*zy (dzm(PA))) ) )
!      PGY_U_UV=   ---- (dym(PA) - (     (------  ) ) )
!                  ___x (          (     ( ___x   ) ) )
!                  d*yy (          (     ( d*zz   ) ) )    
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
!!      MXM,MYM,MZF     : Shuman functions (mean operators)
!!      DYM,DZM         : Shuman functions (finite difference operators)
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
INTEGER,                 INTENT(IN),OPTIONAL   :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                 INTENT(IN),OPTIONAL   :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_U_UV ! result UV point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_U_UV
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  PGY_U_UV(:,:,:)=  (DYM(PA)- MZF( MYM( DZM(PA)/&
                 MXM(PDZZ) ) *MXM(PDZY) )   ) / MXM(PDYY)
ELSE
  PGY_U_UV(:,:,:)= DYM(PA) / MXM(PDYY)
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GY_U_UV
!
!
!     #######################################################
      FUNCTION GZ_U_UW(PA,PDZZ, KKA, KKU, KL)      RESULT(PGZ_U_UW)
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
USE MODI_SHUMAN
!
IMPLICIT NONE
!
!
!*       0.1   declarations of arguments and result
!
INTEGER,              INTENT(IN),OPTIONAL      :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL      :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the U point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_U_UW ! result UW point
!
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
PGZ_U_UW(:,:,:)= DZM(PA) / MXM(PDZZ)
!
!----------------------------------------------------------------------------
!
END FUNCTION GZ_U_UW
