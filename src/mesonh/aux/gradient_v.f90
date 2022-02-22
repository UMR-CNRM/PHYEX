!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_V
!     ######################
!
INTERFACE
!
!           
FUNCTION GY_V_M(PA,PDYY,PDZZ,PDZY)      RESULT(PGY_V_M)
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_V_M ! result mass point
!
END FUNCTION GY_V_M
!           
FUNCTION GX_V_UV(PA,PDXX,PDZZ,PDZX)      RESULT(PGX_V_UV)
!
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
FUNCTION GZ_V_VW(PA,PDZZ)      RESULT(PGZ_V_VW)
!
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the V point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_V_VW ! result VW point
!
END FUNCTION GZ_V_VW
!
!
END INTERFACE
!
END MODULE MODI_GRADIENT_V
!
!
!
!     #######################################################
      FUNCTION GY_V_M(PA,PDYY,PDZZ,PDZY)      RESULT(PGY_V_M)
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
!     #########################################################
      FUNCTION GX_V_UV(PA,PDXX,PDZZ,PDZX)      RESULT(PGX_V_UV)
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
!     #######################################################
      FUNCTION GZ_V_VW(PA,PDZZ)      RESULT(PGZ_V_VW)
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
