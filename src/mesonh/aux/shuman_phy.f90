MODULE MODE_SHUMAN_PHY
IMPLICIT NONE
CONTAINS
!     ###############################
      SUBROUTINE MXM_PHY(D,PA,PMXM)
!     ###############################
!
!!****  *MXM* -  Shuman operator : mean operator in x direction for a 
!!                                 mass variable 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the x direction (I index) for a field PA localized at a mass
!     point. The result is localized at a x-flux point (u point).
!
!!**  METHOD
!!    ------ 
!!        The result PMXM(i,:,:) is defined by 0.5*(PA(i,:,:)+PA(i-1,:,:))
!!    At i=1, PMXM(1,:,:) are replaced by the values of PMXM,
!!    which are the right values in the x-cyclic case. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMXM   ! result at flux localization 

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI,JJ,JK             ! Loop index in x direction
INTEGER :: IIU            ! Size of the array in the x direction
!          
INTEGER :: IJU,IKU
!                     
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MXM
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PA,PMXM)
!$acc loop independent collapse(3)
DO JK = 1, IKU
  DO JJ = 1, IJU
    DO JI = 1 + 1, IIU
      PMXM(JI,JJ,JK) = 0.5*( PA(JI,JJ,JK)+PA(JI-1,JJ,JK) )
    ENDDO
  ENDDO
ENDDO
!
!$acc loop independent collapse(2)
DO JK = 1, IKU
  DO JJ=1,IJU
    PMXM(1,JJ,JK)    = PMXM(IIU-2*JPHEXT+1,JJ,JK)  	!TODO: voir si ce n'est pas plutot JPHEXT+1
  ENDDO
ENDDO
!$acc end kernels
!
END SUBROUTINE MXM_PHY
!-------------------------------------------------------------------------------
!
!     ###############################
      SUBROUTINE MXM2D_PHY(D,PA,PMXM)
!     ###############################
!
!!****  *MXM* -  Shuman operator : mean operator in x direction for a 
!!                                 mass variable 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the x direction (I index) for a field PA localized at a mass
!     point. The result is localized at a x-flux point (u point).
!
!!**  METHOD
!!    ------ 
!!        The result PMXM(i,:,:) is defined by 0.5*(PA(i,:,:)+PA(i-1,:,:))
!!    At i=1, PMXM(1,:,:) are replaced by the values of PMXM,
!!    which are the right values in the x-cyclic case. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PMXM   ! result at flux localization 

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI,JJ            ! Loop index in x direction
INTEGER :: IIU            ! Size of the array in the x direction
!          
INTEGER :: IJU
!                     
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MXM
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
!
!$acc kernels present(PA,PMXM)
!$acc loop independent collapse(2)
  DO JJ = 1, IJU
    DO JI = 1 + 1, IIU
      PMXM(JI,JJ) = 0.5*( PA(JI,JJ)+PA(JI-1,JJ) )
    ENDDO
  ENDDO
!
!$acc loop independent
  DO JJ=1,IJU
    PMXM(1,JJ)    = PMXM(IIU-2*JPHEXT+1,JJ) !TODO: voir si ce n'est pas plutot JPHEXT+1
  ENDDO
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MXM2D_PHY
!
!     ###############################
      SUBROUTINE MZM_PHY(D,PA,PMZM)
!     ###############################
!
!!****  *MZM* -  Shuman operator : mean operator in z direction for a 
!!                                 mass variable 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean
!     along the z direction (K index) for a field PA localized at a mass
!     point. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------ 
!!        The result PMZM(:,:,k) is defined by 0.5*(PA(:,:,k)+PA(:,:,k-1))
!!        At k=1, PMZM(:,:,1) is defined by -999.
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!                   optimisation                 20/08/00 J. Escobar
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMZM   ! result at flux localization 

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK             ! Loop index in z direction
INTEGER :: IKU            ! upper bound in z direction of PA
!           
INTEGER :: IIU,IJU
INTEGER :: JIJ,JI,JJ
!
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZM
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PA,PMZM)
DO JK=2,IKU !TODO: remplacer le 2 par JPHEXT+1 ?
!$mnh_expand_array(JI=1:IIU,JJ=1:IJU)
  PMZM(:,:,JK) = 0.5* ( PA(:,:,JK) + PA(:,:,JK-1) )
!$mnh_end_expand_array(JI=1:IIU,JJ=1:IJU)
END DO
!
!$mnh_expand_array(JI=1:IIU,JJ=1:IJU)
PMZM(:,:,1)    = -999.
!$mnh_end_expand_array(JI=1:IIU,JJ=1:IJU)
!$acc end kernels
!-------------------------------------------------------------------------------
!
END SUBROUTINE MZM_PHY
!     ###############################
      SUBROUTINE MYM2D_PHY(D,PA,PMYM)
!     ###############################
!
!!****  *MYM* -  Shuman operator : mean operator in y direction for a 
!!                                 mass variable 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the y direction (J index) for a field PA localized at a mass
!     point. The result is localized at a y-flux point (v point).
!
!!**  METHOD
!!    ------ 
!!        The result PMYM(:,j,:) is defined by 0.5*(PA(:,j,:)+PA(:,j-1,:))
!!    At j=1, PMYM(:,j,:) are replaced by the values of PMYM,
!!    which are the right values in the y-cyclic case. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1    
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PMYM   ! result at flux localization 

!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ,JI,JK             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!          
INTEGER :: IIU
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MYM
!              ------------------
!
IIU=SIZE(PA,1)
IJU=SIZE(PA,2)
!
!$acc kernels present(PA,PMYM)
!$acc loop independent collapse(2)
  DO JJ = 2, IJU
    DO JI = 1, IIU
      PMYM(JI,JJ) = 0.5*( PA(JI,JJ)+PA(JI,JJ-1) )
    ENDDO
  ENDDO
!
!$acc loop independent collapse(2)
  DO JJ=1,JPHEXT
    DO JI=1,IIU
      PMYM(JI,JJ)  = PMYM(JI,IJU-2*JPHEXT+JJ) ! for reprod JPHEXT <> 1
    ENDDO
  ENDDO
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MYM2D_PHY
!     ###############################
      SUBROUTINE MYM_PHY(D,PA,PMYM)
!     ###############################
!
!!****  *MYM* -  Shuman operator : mean operator in y direction for a 
!!                                 mass variable 
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the y direction (J index) for a field PA localized at a mass
!     point. The result is localized at a y-flux point (v point).
!
!!**  METHOD
!!    ------ 
!!        The result PMYM(:,j,:) is defined by 0.5*(PA(:,j,:)+PA(:,j-1,:))
!!    At j=1, PMYM(:,j,:) are replaced by the values of PMYM,
!!    which are the right values in the y-cyclic case. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1    
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at mass localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMYM   ! result at flux localization 

!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ,JI,JK             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!          
INTEGER :: IIU,IKU
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MYM
!              ------------------
!
IIU=SIZE(PA,1)
IJU=SIZE(PA,2)
IKU=SIZE(PA,3)
!
!$acc kernels present(PA,PMYM)
!$acc loop independent collapse(3)
DO JK = 1, IKU
  DO JJ = 2, IJU
    DO JI = 1, IIU
      PMYM(JI,JJ,JK) = 0.5*( PA(JI,JJ,JK)+PA(JI,JJ-1,JK) )
    ENDDO
  ENDDO
ENDDO
!
!$acc loop independent collapse(3)
DO JK = 1, IKU
  DO JJ=1,JPHEXT
    DO JI=1,IIU
      PMYM(JI,JJ,JK)  = PMYM(JI,IJU-2*JPHEXT+JJ,JK) ! for reprod JPHEXT <> 1
    ENDDO
  ENDDO
ENDDO
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MYM_PHY
!     ###############################
      SUBROUTINE DZM_PHY(D,PA,PDZM)
!     ###############################
!
!!****  *DZM* -  Shuman operator : finite difference operator in z direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the z direction (K index) for a field PA localized at a mass
!     point. The result is localized at a z-flux point (w point).
!
!!**  METHOD
!!    ------ 
!!        The result PDZM(:,j,:) is defined by (PA(:,:,k)-PA(:,:,k-1))
!!        At k=1, PDZM(:,:,k) is defined by -999.
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!                   optimisation                 20/08/00 J. Escobar
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDZM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK,JI,JJ            ! Loop index in z direction
INTEGER :: IKU           ! upper bound in z direction of PA
!
!         
INTEGER :: IIU,IJU
!           
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZM
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PA,PDZM)
!$acc loop independent collapse(3)
DO JK=2,IKU !TODO: remplacer le 1+1 par 1+JPHEXT ?
  DO JJ=1,IJU
    DO JI=1,IIU
      PDZM(JI,JJ,JK) = PA(JI,JJ,JK) - PA(JI,JJ,JK-1) 
    END DO
  END DO
END DO
!
!$mnh_expand_array(JI=1:IIU,JJ=1:IJU)
PDZM(1:IIU,1:IJU,1) = -999.
!$mnh_end_expand_array(JI=1:IIU,JJ=1:IJU)
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DZM_PHY
!     ###############################
      SUBROUTINE MZF_PHY(D,PA,PMZF)
!     ###############################
!
!!****  *MZF* -  Shuman operator : mean operator in z direction for a 
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the z direction (K index) for a field PA localized at a z-flux
!     point (w point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PMZF(:,:,k) is defined by 0.5*(PA(:,:,k)+PA(:,:,k+1))
!!        At k=size(PA,3), PMZF(:,:,k) is defined by -999.
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!                   optimisation                 20/08/00 J. Escobar
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMZF   ! result at mass localization 

!
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK,JI,JJ             ! Loop index in z direction
INTEGER :: IKU          ! upper bound in z direction of PA 
!     
INTEGER :: IIU,IJU
INTEGER :: JIJ
!            
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MZF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PA,PMZF)
PMZF(:,:,1:IKU-1) = 0.5*( PA(:,:,1:IKU-1)+PA(:,:,2:) )
!
!$mnh_expand_array(JI=1:IIU,JJ=1:IJU)
PMZF(1:IIU,1:IJU,IKU) = -999.
!$mnh_end_expand_array(JI=1:IIU,JJ=1:IJU)
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MZF_PHY
!     ###############################
      SUBROUTINE MXF_PHY(D,PA,PMXF)
!     ###############################
!
!!****  *MXF* -  Shuman operator : mean operator in x direction for a 
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the x direction (I index) for a field PA localized at a x-flux
!     point (u point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PMXF(i,:,:) is defined by 0.5*(PA(i,:,:)+PA(i+1,:,:))
!!        At i=size(PA,1), PMXF(i,:,:) are replaced by the values of PMXF,
!!    which are the right values in the x-cyclic case
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1   
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMXF   ! result at mass localization 

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!         
INTEGER :: JJ,JK,IJU,IKU
!          
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MXF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PMXF,PA)
!$acc loop independent collapse(3)
DO JK = 1, IKU
  DO JJ = 1, IJU
    DO JI = 1 + 1, IIU
      PMXF(JI-1,JJ,JK) = 0.5*( PA(JI-1,JJ,JK)+PA(JI,JJ,JK) )
    ENDDO
  ENDDO
ENDDO
!
PMXF(IIU,:,:)    = PMXF(2*JPHEXT,:,:) 
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MXF_PHY
!
!     ###############################
      SUBROUTINE MXF2D_PHY(D,PA,PMXF)
!     ###############################
!
!!****  *MXF* -  Shuman operator : mean operator in x direction for a 
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the x direction (I index) for a field PA localized at a x-flux
!     point (u point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PMXF(i,:,:) is defined by 0.5*(PA(i,:,:)+PA(i+1,:,:))
!!        At i=size(PA,1), PMXF(i,:,:) are replaced by the values of PMXF,
!!    which are the right values in the x-cyclic case
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1   
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PMXF   ! result at mass localization 

!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!         
INTEGER :: JJ,IJU
!          
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MXF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
!
!$acc kernels present(PMXF,PA)
!$acc loop independent collapse(2)
  DO JJ = 1, IJU
    DO JI = 1 + 1, IIU
      PMXF(JI-1,JJ) = 0.5*( PA(JI-1,JJ)+PA(JI,JJ) )
    ENDDO
  ENDDO
!
PMXF(IIU,:)    = PMXF(2*JPHEXT,:) 
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MXF2D_PHY
!     ###############################
      SUBROUTINE MYF_PHY(D,PA,PMYF)
!     ###############################
!
!!****  *MYF* -  Shuman operator : mean operator in y direction for a 
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the y direction (J index) for a field PA localized at a y-flux
!     point (v point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PMYF(i,:,:) is defined by 0.5*(PA(:,j,:)+PA(:,j+1,:))
!!        At j=size(PA,2), PMYF(:,j,:) are replaced by the values of PMYF,
!!    which are the right values in the y-cyclic case
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1   
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PMYF   ! result at mass localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ,JI,JK             ! Loop index in y direction
INTEGER :: IJU            ! upper bound in y direction of PA 
!           
INTEGER :: IIU,IKU
!                
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MYF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PA,PMYF)
!$acc loop collapse(3) independent 
DO JK = 1, IKU
  DO JJ = 1, IJU-1
    DO JI = 1, IIU
      PMYF(JI,JJ,JK) = 0.5*( PA(JI,JJ,JK)+PA(JI,JJ+1,JK) )
    END DO
  END DO
END DO
!
!$acc loop seq
DO JJ=1,JPHEXT
  !$acc loop collapse(2) independent 
  DO JK = 1, IKU
    DO JI = 1, IIU
      PMYF(JI,IJU-JPHEXT+JJ,JK) = PMYF(JI,JPHEXT+JJ,JK) ! for reprod JPHEXT <>
    END DO
  END DO
END DO
!$acc end kernels
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MYF_PHY
!
!     ###############################
      SUBROUTINE MYF2D_PHY(D,PA,PMYF)
!     ###############################
!
!!****  *MYF* -  Shuman operator : mean operator in y direction for a 
!!                                 variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a mean 
!     along the y direction (J index) for a field PA localized at a y-flux
!     point (v point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PMYF(i,:,:) is defined by 0.5*(PA(:,j,:)+PA(:,j+1,:))
!!        At j=size(PA,2), PMYF(:,j,:) are replaced by the values of PMYF,
!!    which are the right values in the y-cyclic case
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    04/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1   
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1  
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PMYF   ! result at mass localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ,JI             ! Loop index in y direction
INTEGER :: IJU            ! upper bound in y direction of PA 
!           
INTEGER :: IIU
!                
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MYF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
!
!$acc kernels present(PA,PMYF)
!$acc loop collapse(2) independent 
DO JJ = 1, IJU-1
  DO JI = 1, IIU
    PMYF(JI,JJ) = 0.5*( PA(JI,JJ)+PA(JI,JJ+1) )
  END DO
END DO
!
!$acc loop seq
DO JJ=1,JPHEXT
  !$acc loop independent 
  DO JI = 1, IIU
    PMYF(JI,IJU-JPHEXT+JJ) = PMYF(JI,JPHEXT+JJ) ! for reprod JPHEXT <>
  END DO
END DO
!$acc end kernels
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE MYF2D_PHY
!     ###############################
      SUBROUTINE DZF_PHY(D,PA,PDZF)
!     ###############################
!
!!****  *DZF* -  Shuman operator : finite difference operator in z direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the z direction (K index) for a field PA localized at a z-flux
!     point (w point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PDZF(:,:,k) is defined by (PA(:,:,k+1)-PA(:,:,k))
!!        At k=size(PA,3), PDZF(:,:,k) is defined by -999.
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      NONE
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!                   optimisation                 20/08/00 J. Escobar
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)  :: PA     ! variable at flux localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDZF   ! result at mass localization 
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JK,JI,JJ           ! Loop index in z direction
INTEGER :: IKU          ! upper bound in z direction of PA 
!
!           
INTEGER :: IIU,IJU
!         
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DZF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PA,PDZF)
!$acc loop independent collapse(3)
DO JK=1,IKU-1 !TODO: remplacer le 1 par JPVEXT ?
  DO JJ=1,IJU
    DO JI=1,IIU
      PDZF(JI,JJ,JK) = PA(JI,JJ,JK+1)-PA(JI,JJ,JK)
    END DO
  END DO
END DO
!
PDZF(:,:,IKU) = -999.
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DZF_PHY
!
!     ###############################
      SUBROUTINE DXF_PHY(D,PA,PDXF)
!     ###############################
!
!!****  *DXF* -  Shuman operator : finite difference operator in x direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the x direction (I index) for a field PA localized at a x-flux
!     point (u point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PDXF(i,:,:) is defined by (PA(i+1,:,:)-PA(i,:,:))
!!        At i=size(PA,1), PDXF(i,:,:) are replaced by the values of PDXF,
!!    which are the right values in the x-cyclic case
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1    
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPHEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDXF   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI,JJ,JK             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!             
INTEGER :: IJU,IKU
!             
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DXF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PA,PDXF)
!$acc loop independent collapse(3)
DO JK=1,IKU
  DO JJ=1,IJU
    DO JI=1+1,IIU
     PDXF(JI-1,JJ,JK) = PA(JI,JJ,JK) - PA(JI-1,JJ,JK) 
    END DO
  END DO
END DO
!
!$acc loop independent collapse(2)
DO JK=1,IKU
  DO JJ=1,IJU
    PDXF(IIU,JJ,JK)    = PDXF(2*JPHEXT,JJ,JK) 
  ENDDO
ENDDO
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DXF_PHY
!
!     ###############################
      SUBROUTINE DXF2D_PHY(D,PA,PDXF)
!     ###############################
!
!!****  *DXF* -  Shuman operator : finite difference operator in x direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the x direction (I index) for a field PA localized at a x-flux
!     point (u point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PDXF(i,:,:) is defined by (PA(i+1,:,:)-PA(i,:,:))
!!        At i=size(PA,1), PDXF(i,:,:) are replaced by the values of PDXF,
!!    which are the right values in the x-cyclic case
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1    
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPHEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PDXF   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!             
INTEGER :: JJ,IJU
!             
!
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DXF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
!
!$acc kernels present(PA,PDXF)
!$acc loop independent collapse(2)
  DO JJ=1,IJU
    DO JI=1+1,IIU
     PDXF(JI-1,JJ) = PA(JI,JJ) - PA(JI-1,JJ) 
    END DO
  END DO
!
!$acc loop independent
  DO JJ=1,IJU
    PDXF(IIU,JJ)    = PDXF(2*JPHEXT,JJ) 
  ENDDO
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DXF2D_PHY
!     ###############################
      SUBROUTINE DXM_PHY(D,PA,PDXM)
!     ###############################
!
!!****  *DXM* -  Shuman operator : finite difference operator in x direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the x direction (I index) for a field PA localized at a mass
!     point. The result is localized at a x-flux point (u point).
!
!!**  METHOD
!!    ------ 
!!        The result PDXM(i,:,:) is defined by (PA(i,:,:)-PA(i-1,:,:))
!!    At i=1, PDXM(1,:,:) are replaced by the values of PDXM,
!!    which are the right values in the x-cyclic case. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDXM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER :: JI,JJ,JK             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!             
INTEGER :: IJU,IKU
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DXM
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present(PA,PDXM)
!$acc loop independent collapse(3)
DO JK=1,IKU
  DO JJ=1,IJU
    DO JI=1+1,IIU !TODO: remplacer le 1 par JPHEXT ?
      PDXM(JI,JJ,JK) = PA(JI,JJ,JK) - PA(JI-1,JJ,JK) 
    END DO
  END DO
END DO
!
!$acc loop independent collapse(2)
DO JK=1,IKU
  DO JJ=1,IJU
    PDXM(1,JJ,JK)    = PDXM(IIU-2*JPHEXT+1,JJ,JK)   !TODO: remplacer -2*JPHEXT+1 par -JPHEXT ?
  ENDDO
ENDDO
!$acc end kernels
!
END SUBROUTINE DXM_PHY
!-------------------------------------------------------------------------------
!
      SUBROUTINE DXM2D_PHY(D,PA,PDXM)
!     ###############################
!
!!****  *DXM* -  Shuman operator : finite difference operator in x direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the x direction (I index) for a field PA localized at a mass
!     point. The result is localized at a x-flux point (u point).
!
!!**  METHOD
!!    ------ 
!!        The result PDXM(i,:,:) is defined by (PA(i,:,:)-PA(i-1,:,:))
!!    At i=1, PDXM(1,:,:) are replaced by the values of PDXM,
!!    which are the right values in the x-cyclic case. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PDXM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!             
INTEGER :: JJ,IJU
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DXM
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
!
!$acc kernels present(PA,PDXM)
!$acc loop independent collapse(2)
  DO JJ=1,IJU
    DO JI=1+1,IIU !TODO: remplacer le 1 par JPHEXT ?
      PDXM(JI,JJ) = PA(JI,JJ) - PA(JI-1,JJ) 
    END DO
  END DO
!
!$acc loop independent
  DO JJ=1,IJU
    PDXM(1,JJ)    = PDXM(IIU-2*JPHEXT+1,JJ)   !TODO: remplacer -2*JPHEXT+1 par -JPHEXT ?
  ENDDO
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DXM2D_PHY
!
!     ###############################
      SUBROUTINE DYM_PHY(D,PA,PDYM)
!     ###############################
!
!!****  *DYM* -  Shuman operator : finite difference operator in y direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the y direction (J index) for a field PA localized at a mass
!     point. The result is localized at a y-flux point (v point).
!
!!**  METHOD
!!    ------ 
!!        The result PDYM(:,j,:) is defined by (PA(:,j,:)-PA(:,j-1,:))
!!    At j=1, PDYM(:,1,:) are replaced by the values of PDYM,
!!    which are the right values in the y-cyclic case. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDYM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ,JI,JK             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!    
INTEGER :: IIU,IKU
!     
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DYM
!              ------------------
!
IIU=SIZE(PA,1)
IJU=SIZE(PA,2)
IKU=SIZE(PA,3)
!
!$acc kernels present(PA,PDYM)
!$acc loop independent collapse(3)
DO JK=1,IKU
  DO JJ=2,IJU !TODO: remplacer le 2 par JPHEXT+1 ?
    DO JI=1,IIU
      PDYM(JI,JJ,JK) = PA(JI,JJ,JK) - PA(JI,JJ-1,JK) 
    END DO
  END DO
END DO
!
!$acc loop seq
DO JJ=1,JPHEXT
  !$acc loop collapse(2) independent 
  DO JK=1,IKU
    DO JI=1,IIU
     PDYM(JI,JJ,JK) = PDYM(JI,IJU-2*JPHEXT+JJ,JK) ! for reprod JPHEXT <> 1
    END DO
  END DO
END DO
!$acc end kernels
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DYM_PHY
!     ###############################
      SUBROUTINE DYM2D_PHY(D,PA,PDYM)
!     ###############################
!
!!****  *DYM* -  Shuman operator : finite difference operator in y direction
!!                                  for a variable at a mass localization
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the y direction (J index) for a field PA localized at a mass
!     point. The result is localized at a y-flux point (v point).
!
!!**  METHOD
!!    ------ 
!!        The result PDYM(:,j,:) is defined by (PA(:,j,:)-PA(:,j-1,:))
!!    At j=1, PDYM(:,1,:) are replaced by the values of PDYM,
!!    which are the right values in the y-cyclic case. 
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_PARAMETERS, ONLY: JPHEXT
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PDYM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ,JI            ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!    
INTEGER :: IIU
!     
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DYM
!              ------------------
!
IIU=SIZE(PA,1)
IJU=SIZE(PA,2)
!
!$acc kernels present(PA,PDYM)
!$acc loop independent collapse(2)
  DO JJ=2,IJU !TODO: remplacer le 2 par JPHEXT+1 ?
    DO JI=1,IIU
      PDYM(JI,JJ) = PA(JI,JJ) - PA(JI,JJ-1) 
    END DO
  END DO
!
!$acc loop seq
DO JJ=1,JPHEXT
  !$acc loop independent 
  DO JI=1,IIU
    PDYM(JI,JJ) = PDYM(JI,IJU-2*JPHEXT+JJ) ! for reprod JPHEXT <> 1
  END DO
END DO
!$acc end kernels
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DYM2D_PHY
!     ###############################
      SUBROUTINE DYF_PHY(D,PA,PDYF)
!     ###############################
!
!!****  *DYF* -  Shuman operator : finite difference operator in y direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the y direction (J index) for a field PA localized at a y-flux
!     point (v point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PDYF(:,j,:) is defined by (PA(:,j+1,:)-PA(:,j,:))
!!        At j=size(PA,2), PDYF(:,j,:) are replaced by the values of PDYM,
!!    which are the right values in the y-cyclic case
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!      V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPHEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDYF   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ,JI,JK            ! Loop index in y direction
INTEGER :: IJU           ! upper bound in y direction of PA 
!
!          
INTEGER :: IIU,IKU
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DYF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
IKU = SIZE(PA,3)
!
!$acc kernels present_crm(PDYF,PA)
!$mnh_expand_array(JI=1:IIU,JJ=1:IJU-1,JK=1:IKU)
  PDYF(1:IIU,1:IJU-1,1:IKU) = PA(1:IIU,2:IJU+1,1:IKU) - PA(1:IIU,1:IJU-1,1:IKU)
!$mnh_end_expand_array(JI=1:IIU,JJ=1:IJU-1,JK=1:IKU)
!
!$acc loop seq
DO JJ=1,JPHEXT
  !$acc loop collapse(2) independent 
  DO JK=1,IKU
    DO JI=1,IIU
      PDYF(JI,IJU-JPHEXT+JJ,JK) = PDYF(JI,JPHEXT+JJ,JK) ! for reprod JPHEXT <>
    END DO
  END DO
END DO
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DYF_PHY
!
!     ###############################
      SUBROUTINE DYF2D_PHY(D,PA,PDYF)
!     ###############################
!
!!****  *DYF* -  Shuman operator : finite difference operator in y direction
!!                                  for a variable at a flux side
!!
!!    PURPOSE
!!    -------
!       The purpose of this function  is to compute a finite difference 
!     along the y direction (J index) for a field PA localized at a y-flux
!     point (v point). The result is localized at a mass point.
!
!!**  METHOD
!!    ------ 
!!        The result PDYF(:,j,:) is defined by (PA(:,j+1,:)-PA(:,j,:))
!!        At j=size(PA,2), PDYF(:,j,:) are replaced by the values of PDYM,
!!    which are the right values in the y-cyclic case
!!    
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      Module MODD_PARAMETERS: declaration of parameter variables
!!        JPHEXT: define the number of marginal points out of the 
!!        physical domain along the horizontal directions.
!!
!!    REFERENCE
!!    ---------
!!      Book2 of documentation of Meso-NH (SHUMAN operators)
!!      Technical specifications Report of The Meso-NH (chapters 3)  
!!
!!
!!    AUTHOR
!!    ------
!!	V. Ducrocq       * Meteo France *
!!
!!    MODIFICATIONS
!!    -------------
!!      Original    05/07/94 
!!      Modification to include the periodic case 13/10/94 J.Stein 
!!                   optimisation                 20/08/00 J. Escobar
!!      correction of in halo/pseudo-cyclic calculation for JPHEXT<> 1 
!!      J.Escobar : 15/09/2015 : WENO5 & JPHEXT <> 1 
!-------------------------------------------------------------------------------
!
!*       0.    DECLARATIONS
!              ------------
!
USE MODD_DIMPHYEX, ONLY: DIMPHYEX_t
USE MODD_PARAMETERS, ONLY: JPHEXT
!
IMPLICIT NONE
!
!*       0.1   Declarations of argument and result
!              ------------------------------------
!
TYPE(DIMPHYEX_t),       INTENT(IN)  :: D
REAL, DIMENSION(D%NIT,D%NJT), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(D%NIT,D%NJT), INTENT(OUT) :: PDYF   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ ,JI           ! Loop index in y direction
INTEGER :: IJU           ! upper bound in y direction of PA 
!
!          
INTEGER :: IIU
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF DYF
!              ------------------
!
IIU = SIZE(PA,1)
IJU = SIZE(PA,2)
!
!$acc kernels present_crm(PDYF,PA)
!$mnh_expand_array(JI=1:IIU,JJ=1:IJU-1)
  PDYF(1:IIU,1:IJU-1) = PA(1:IIU,2:IJU+1) - PA(1:IIU,1:IJU-1)
!$mnh_end_expand_array(JI=1:IIU,JJ=1:IJU-1)
!
!$acc loop seq
DO JJ=1,JPHEXT
!$mnh_expand_array(JI=1:IIU)
   PDYF(1:IIU,IJU-JPHEXT+JJ) = PDYF(1:IIU,JPHEXT+JJ) ! for reprod JPHEXT <> 1
!$mnh_end_expand_array(JI=1:IIU)
END DO
!$acc end kernels
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DYF2D_PHY

END MODULE MODE_SHUMAN_PHY
