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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JJK,IJU,IKU
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + 1 
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PMXM(JIJK,1,1) = 0.5*( PA(JIJK,1,1)+PA(JIJK-1,1,1) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JI=1,JPHEXT
   DO JJK=1,IJU*IKU
      PMXM(JI,JJK,1) = PMXM(IIU-2*JPHEXT+JI,JJK,1) ! for reprod JPHEXT <> 1
   END DO
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JJK,IJU
INTEGER :: JIJ,JIJOR,JIJEND
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
JIJOR  = 1 + 1 
JIJEND = IIU*IJU
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=JIJOR , JIJEND
   PMXM(JIJ,1) = 0.5*( PA(JIJ,1)+PA(JIJ-1,1) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JI=1,JPHEXT
   DO JJ=1,IJU
      PMXM(JI,JJ) = PMXM(IIU-2*JPHEXT+JI,JJ) ! for reprod JPHEXT <> 1
   END DO
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + IIU*IJU
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PMZM(JIJK,1,1) = 0.5*( PA(JIJK,1,1)+PA(JIJK-IIU*IJU,1,1) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=1,IIU*IJU
   PMZM(JIJ,1,1)    = -999.
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!          
INTEGER :: IIU
INTEGER :: JIJ,JIJOR,JIJEND
!            
!-------------------------------------------------------------------------------
!
!*       1.    DEFINITION OF MYM
!              ------------------
!
IIU=SIZE(PA,1)
IJU=SIZE(PA,2)
!
JIJOR  = 1 + IIU
JIJEND = IIU*IJU
!CDIR NODEP
!OCL NOVREC
DO JIJ=JIJOR , JIJEND
   PMYM(JIJ,1) = 0.5*( PA(JIJ,1)+PA(JIJ-IIU,1) )
END DO
!
DO JJ=1,JPHEXT
   PMYM(:,JJ)  = PMYM(:,IJU-2*JPHEXT+JJ) ! for reprod JPHEXT <> 1
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!          
INTEGER :: IIU,IKU
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + IIU
JIJKEND = IIU*IJU*IKU
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PMYM(JIJK,1,1) = 0.5*( PA(JIJK,1,1)+PA(JIJK-IIU,1,1) )
END DO
!
DO JJ=1,JPHEXT
   PMYM(:,JJ,:)  = PMYM(:,IJU-2*JPHEXT+JJ,:) ! for reprod JPHEXT <> 1
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JK            ! Loop index in z direction
INTEGER :: IKU           ! upper bound in z direction of PA
!
!         
INTEGER :: IIU,IJU
INTEGER :: JIJ
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + IIU*IJU
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PDZM(JIJK,1,1) = PA(JIJK,1,1)-PA(JIJK-IIU*IJU,1,1)
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=1,IIU*IJU
   PDZM(JIJ,1,1)    = -999.
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JK             ! Loop index in z direction
INTEGER :: IKU          ! upper bound in z direction of PA 
!     
INTEGER :: IIU,IJU
INTEGER :: JIJ
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + IIU*IJU
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PMZF(JIJK-IIU*IJU,1,1) = 0.5*( PA(JIJK-IIU*IJU,1,1)+PA(JIJK,1,1) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=1,IIU*IJU
   PMZF(JIJ,1,IKU)    = PMZF(JIJ,1,IKU-1) !-999.
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JJK,IJU,IKU
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + 1 
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
  PMXF(JIJK-1,1,1) = 0.5*( PA(JIJK-1,1,1)+PA(JIJK,1,1) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JI=1,JPHEXT
   DO JJK=1,IJU*IKU
      PMXF(IIU-JPHEXT+JI,JJK,1) = PMXF(JPHEXT+JI,JJK,1) ! for reprod JPHEXT <> 1
   END DO
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JIJ,JIJOR,JIJEND
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
JIJOR  = 1 + 1 
JIJEND = IIU*IJU
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=JIJOR , JIJEND
  PMXF(JIJ-1,1) = 0.5*( PA(JIJ-1,1)+PA(JIJ,1) )
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JI=1,JPHEXT
   DO JJ=1,IJU
      PMXF(IIU-JPHEXT+JI,JJ) = PMXF(JPHEXT+JI,JJ) ! for reprod JPHEXT <> 1
   END DO
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! upper bound in y direction of PA 
!           
INTEGER :: IIU,IKU
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + IIU
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PMYF(JIJK-IIU,1,1) = 0.5*( PA(JIJK-IIU,1,1)+PA(JIJK,1,1) )
END DO
!
DO JJ=1,JPHEXT
   PMYF(:,IJU-JPHEXT+JJ,:) = PMYF(:,JPHEXT+JJ,:) ! for reprod JPHEXT <> 1
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! upper bound in y direction of PA 
!           
INTEGER :: IIU
INTEGER :: JIJ,JIJOR,JIJEND
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
JIJOR  = 1 + IIU
JIJEND = IIU*IJU
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=JIJOR , JIJEND
   PMYF(JIJ-IIU,1) = 0.5*( PA(JIJ-IIU,1)+PA(JIJ,1) )
END DO
!
DO JJ=1,JPHEXT
   PMYF(:,IJU-JPHEXT+JJ) = PMYF(:,JPHEXT+JJ) ! for reprod JPHEXT <> 1
END DO
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
!!	V. Ducrocq       * Meteo France *
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
INTEGER :: JK           ! Loop index in z direction
INTEGER :: IKU          ! upper bound in z direction of PA 
!
!           
INTEGER :: IIU,IJU
INTEGER :: JIJ
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + IIU*IJU
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PDZF(JIJK-IIU*IJU,1,1)     = PA(JIJK,1,1)-PA(JIJK-IIU*IJU,1,1)
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JIJ=1,IIU*IJU
   PDZF(JIJ,1,IKU)    = -999.
END DO
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDXF   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!             
INTEGER :: JJK,IJU,IKU
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + 1
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PDXF(JIJK-1,1,1) = PA(JIJK,1,1) - PA(JIJK-1,1,1) 
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JI=1,JPHEXT
   DO JJK=1,IJU*IKU
      PDXF(IIU-JPHEXT+JI,JJK,1) = PDXF(JPHEXT+JI,JJK,1) ! for reprod JPHEXT <> 1
   END DO
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DXF_PHY
!
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDXM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
!
INTEGER :: JI             ! Loop index in x direction
INTEGER :: IIU            ! upper bound in x direction of PA 
!             
INTEGER :: JJK,IJU,IKU
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + 1
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PDXM(JIJK,1,1) = PA(JIJK,1,1) - PA(JIJK-1,1,1) 
END DO
!
!CDIR NODEP
!OCL NOVREC
DO JI=1,JPHEXT
   DO JJK=1,IJU*IKU
      PDXM(JI,JJK,1) = PDXM(IIU-2*JPHEXT+JI,JJK,1) ! for reprod JPHEXT <> 1
   END DO
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DXM_PHY
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT),  INTENT(IN)                :: PA     ! variable at mass
                                                            ! localization
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDYM   ! result at flux
                                                            ! side
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ             ! Loop index in y direction
INTEGER :: IJU            ! Size of the array in the y direction
!
!    
INTEGER :: IIU,IKU
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + IIU
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PDYM(JIJK,1,1)           = PA(JIJK,1,1)  -  PA(JIJK-IIU,1,1) 
END DO
!
DO JJ=1,JPHEXT
   PDYM(:,JJ,:) = PDYM(:,IJU-2*JPHEXT+JJ,:) ! for reprod JPHEXT <> 1
END DO
!
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DYM_PHY
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
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(IN)                :: PA     ! variable at flux
                                                            !  side
REAL, DIMENSION(D%NIT,D%NJT,D%NKT), INTENT(OUT) :: PDYF   ! result at mass
                                                            ! localization 
!
!*       0.2   Declarations of local variables
!              -------------------------------
!
INTEGER :: JJ            ! Loop index in y direction
INTEGER :: IJU           ! upper bound in y direction of PA 
!
!          
INTEGER :: IIU,IKU
INTEGER :: JIJK,JIJKOR,JIJKEND
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
JIJKOR  = 1 + IIU
JIJKEND = IIU*IJU*IKU
!
!CDIR NODEP
!OCL NOVREC
DO JIJK=JIJKOR , JIJKEND
   PDYF(JIJK-IIU,1,1)         = PA(JIJK,1,1)  -  PA(JIJK-IIU,1,1) 
END DO
!
DO JJ=1,JPHEXT
   PDYF(:,IJU-JPHEXT+JJ,:) = PDYF(:,JPHEXT+JJ,:) ! for reprod JPHEXT <> 1
END DO
!
!-------------------------------------------------------------------------------
!
END SUBROUTINE DYF_PHY
!
END MODULE MODE_SHUMAN_PHY
