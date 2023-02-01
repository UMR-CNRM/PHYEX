!MNH_LIC Copyright 1994-2020 CNRS, Meteo-France and Universite Paul Sabatier
!MNH_LIC This is part of the Meso-NH software governed by the CeCILL-C licence
!MNH_LIC version 1. See LICENSE, CeCILL-C_V1-en.txt and CeCILL-C_V1-fr.txt
!MNH_LIC for details. version 1.
!-----------------------------------------------------------------
!     ######################
      MODULE MODI_GRADIENT_W
!     ######################
!
INTERFACE
!
!            
FUNCTION GZ_W_M(PA,PDZZ, KKA, KKU, KL)      RESULT(PGZ_W_M)
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_W_M ! result mass point
!
END FUNCTION GZ_W_M
!            
FUNCTION GX_W_UW(PA,PDXX,PDZZ,PDZX, KKA, KKU, KL)      RESULT(PGX_W_UW)
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_W_UW ! result UW point
!
END FUNCTION GX_W_UW
!
!            
FUNCTION GY_W_VW(PA,PDYY,PDZZ,PDZY, KKA, KKU, KL)      RESULT(PGY_W_VW)
!
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dyy
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzy
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_W_VW ! result VW point
!
END FUNCTION GY_W_VW
!
!
END INTERFACE
!
END MODULE MODI_GRADIENT_W
!
!
!
!     #######################################################
      FUNCTION GZ_W_M(PA,PDZZ, KKA, KKU, KL)      RESULT(PGZ_W_M)
!     #######################################################
!
!!****  *GZ_W_M* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Z
!!                          direction for a variable placed at the 
!!                          W point and the result is placed at
!!                          the mass point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Z cartesian direction for a field PA placed at the 
!     W point. The result is placed at the mass point.
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
!!      DZF     : Shuman functions (finite difference operators)
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
INTEGER,                 INTENT(IN),OPTIONAL  :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,                 INTENT(IN),OPTIONAL  :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGZ_W_M ! result mass point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GZ_W_M
!              --------------------
!
PGZ_W_M(:,:,:)= DZF(PA(:,:,:))/(MZF(PDZZ(:,:,:)))
!
!----------------------------------------------------------------------------
!
END FUNCTION GZ_W_M
!
!
!     #########################################################
      FUNCTION GX_W_UW(PA,PDXX,PDZZ,PDZX, KKA, KKU, KL)      RESULT(PGX_W_UW)
!     #########################################################
!
!!****  *GX_W_UW* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian X
!!                          direction for a variable placed at the 
!!                          V point and the result is placed at
!!                          the UW vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the X cartesian direction for a field PA placed at the 
!     W point. The result is placed at the UW vorticity point.
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MXM,MZM,MZF     : Shuman functions (mean operators)
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
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDXX    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZX    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGX_W_UW ! result UW point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GX_W_UW
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  PGX_W_UW(:,:,:)= DXM(PA(:,:,:))/(MZM(PDXX(:,:,:)))    &
                 -DZM(MXM(MZF(PA(:,:,:))))*PDZX(:,:,:)  &
                  /( MZM(PDXX(:,:,:))*MXM(PDZZ(:,:,:)) )
ELSE
  PGX_W_UW(:,:,:)= DXM(PA(:,:,:))/(MZM(PDXX(:,:,:)))
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GX_W_UW
!
!
!     #########################################################
      FUNCTION GY_W_VW(PA,PDYY,PDZZ,PDZY, KKA, KKU, KL)      RESULT(PGY_W_VW)
!     #########################################################
!
!!****  *GY_W_VW* - Cartesian Gradient operator: 
!!                          computes the gradient in the cartesian Y
!!                          direction for a variable placed at the 
!!                          W point and the result is placed at
!!                          the VW vorticity point.
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the discrete gradient 
!     along the Y cartesian direction for a field PA placed at the 
!     W point. The result is placed at the VW vorticity point.
!
!!**  METHOD
!!    ------
!!      The Chain rule of differencing is applied to variables expressed
!!    in the Gal-Chen & Somerville coordinates to obtain the gradient in
!!    the cartesian system
!!        
!!    EXTERNAL
!!    --------
!!      MYM,MZM,MZF     : Shuman functions (mean operators)
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
INTEGER,              INTENT(IN),OPTIONAL     :: KKA, KKU ! near ground and uppest atmosphere array indexes
INTEGER,              INTENT(IN),OPTIONAL     :: KL     ! +1 if grid goes from ground to atmosphere top, -1 otherwise
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PA      ! variable at the W point
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDYY    ! metric coefficient dxx
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZZ    ! metric coefficient dzz
REAL, DIMENSION(:,:,:),  INTENT(IN)  :: PDZY    ! metric coefficient dzx
!
REAL, DIMENSION(SIZE(PA,1),SIZE(PA,2),SIZE(PA,3)) :: PGY_W_VW ! result VW point
!
!
!*       0.2   declaration of local variables
!
!              NONE
!
!----------------------------------------------------------------------------
!
!*       1.    DEFINITION of GY_W_VW
!              ---------------------
!
IF (.NOT. LFLAT) THEN
  PGY_W_VW(:,:,:)= DYM(PA(:,:,:))/(MZM(PDYY(:,:,:)))    &
                 -DZM(MYM(MZF(PA(:,:,:))))*PDZY(:,:,:)  &
                  /( MZM(PDYY(:,:,:))*MYM(PDZZ(:,:,:)) )
ELSE
  PGY_W_VW(:,:,:)= DYM(PA(:,:,:))/(MZM(PDYY(:,:,:)))
END IF
!
!----------------------------------------------------------------------------
!
END FUNCTION GY_W_VW
